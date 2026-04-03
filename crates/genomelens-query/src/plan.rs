use genomelens_core::{Error, Result, VariantType, VcfNumber, VcfValueType};
use genomelens_parse::VcfHeader;

use crate::ast::{CmpOp, Column, Expr, FixedColumn, Literal};
use crate::vm::OpCode;

/// Compile an AST filter expression into a flat Vec<OpCode> for the VM.
/// Validates column types and INFO field existence against the VCF header.
pub fn plan(expr: &Expr, header: &VcfHeader) -> Result<Vec<OpCode>> {
    let mut ops = Vec::new();
    compile_expr(expr, header, &mut ops)?;
    Ok(ops)
}

fn compile_expr(expr: &Expr, header: &VcfHeader, ops: &mut Vec<OpCode>) -> Result<()> {
    match expr {
        Expr::Compare { column, op, value } => {
            compile_compare(column, *op, value, header, ops)
        }
        Expr::FlagCheck { field } => {
            compile_flag(field, header, ops)
        }
        Expr::And(lhs, rhs) => {
            compile_expr(lhs, header, ops)?;
            let jump_idx = ops.len();
            ops.push(OpCode::JumpIfFalse(0)); // placeholder
            compile_expr(rhs, header, ops)?;
            ops.push(OpCode::And);
            ops[jump_idx] = OpCode::JumpIfFalse(ops.len()); // backpatch
            Ok(())
        }
        Expr::Or(lhs, rhs) => {
            compile_expr(lhs, header, ops)?;
            let jump_idx = ops.len();
            ops.push(OpCode::JumpIfTrue(0)); // placeholder
            compile_expr(rhs, header, ops)?;
            ops.push(OpCode::Or);
            ops[jump_idx] = OpCode::JumpIfTrue(ops.len()); // backpatch
            Ok(())
        }
        Expr::Not(inner) => {
            compile_expr(inner, header, ops)?;
            ops.push(OpCode::Not);
            Ok(())
        }
    }
}

fn compile_compare(
    column: &Column,
    op: CmpOp,
    value: &Literal,
    header: &VcfHeader,
    ops: &mut Vec<OpCode>,
) -> Result<()> {
    match column {
        Column::Fixed(fixed) => compile_fixed_compare(fixed, op, value, ops),
        Column::Info(field) => compile_info_compare(field, op, value, header, ops),
    }
}

fn compile_fixed_compare(
    col: &FixedColumn,
    op: CmpOp,
    value: &Literal,
    ops: &mut Vec<OpCode>,
) -> Result<()> {
    match col {
        FixedColumn::Chrom => {
            let s = expect_str(value, "CHROM")?;
            ops.push(OpCode::CmpChrom(op, s.into_bytes()));
        }
        FixedColumn::Pos => {
            let v = expect_int(value, "POS")?;
            if v < 0 {
                return Err(Error::InvalidQuery("POS cannot be negative".to_string()));
            }
            ops.push(OpCode::CmpPos(op, v as usize));
        }
        FixedColumn::Qual => {
            let v = coerce_to_f64(value, "QUAL")?;
            ops.push(OpCode::CmpQual(op, v));
        }
        FixedColumn::Filter => {
            let s = expect_str(value, "FILTER")?;
            ops.push(OpCode::CmpFilter(op, s.into_bytes()));
        }
        FixedColumn::Ref => {
            let s = expect_str(value, "REF")?;
            ops.push(OpCode::CmpRef(op, s.into_bytes()));
        }
        FixedColumn::Type => {
            let s = expect_str(value, "TYPE")?;
            let vtype = parse_variant_type(&s)?;
            ops.push(OpCode::CmpType(op, vtype));
        }
    }
    Ok(())
}

fn is_scalar(number: &VcfNumber) -> bool {
    matches!(number, VcfNumber::Count(0) | VcfNumber::Count(1))
}

fn compile_info_compare(
    field: &str,
    op: CmpOp,
    value: &Literal,
    header: &VcfHeader,
    ops: &mut Vec<OpCode>,
) -> Result<()> {
    let key = field.as_bytes().to_vec();

    match header.info_field(field) {
        Some(def) => {
            let scalar = is_scalar(&def.number);
            match def.ty {
                VcfValueType::Integer => {
                    let v = expect_int(value, field)?;
                    if scalar {
                        ops.push(OpCode::CmpInfoInt { key, op, value: v });
                    } else {
                        ops.push(OpCode::AnyCmpInfoInt { key, op, value: v });
                    }
                }
                VcfValueType::Float => {
                    let v = coerce_to_f64(value, field)?;
                    if scalar {
                        ops.push(OpCode::CmpInfoFloat { key, op, value: v });
                    } else {
                        ops.push(OpCode::AnyCmpInfoFloat { key, op, value: v });
                    }
                }
                VcfValueType::String | VcfValueType::Character => {
                    let s = expect_str(value, field)?;
                    if scalar {
                        ops.push(OpCode::CmpInfoStr { key, op, value: s.into_bytes() });
                    } else {
                        ops.push(OpCode::AnyCmpInfoStr { key, op, value: s.into_bytes() });
                    }
                }
                VcfValueType::Flag => {
                    return Err(Error::InvalidQuery(format!(
                        "INFO.{} is a Flag — use it without a comparison operator",
                        field
                    )));
                }
            }
        }
        None => {
            return Err(Error::InvalidQuery(format!(
                "INFO field '{}' not found in header",
                field
            )));
        }
    }
    Ok(())
}

fn compile_flag(field: &str, header: &VcfHeader, ops: &mut Vec<OpCode>) -> Result<()> {
    match header.info_field(field) {
        Some(def) if def.ty == VcfValueType::Flag => {
            ops.push(OpCode::InfoFlag(field.as_bytes().to_vec()));
            Ok(())
        }
        Some(def) => Err(Error::InvalidQuery(format!(
            "INFO.{} is {:?}, not a Flag — use a comparison operator",
            field, def.ty
        ))),
        None => Err(Error::InvalidQuery(format!(
            "INFO field '{}' not found in header",
            field
        ))),
    }
}

fn expect_str(value: &Literal, context: &str) -> Result<String> {
    match value {
        Literal::Str(s) => Ok(s.clone()),
        _ => Err(Error::InvalidQuery(format!(
            "{} requires a string literal",
            context
        ))),
    }
}

fn expect_int(value: &Literal, context: &str) -> Result<i64> {
    match value {
        Literal::Integer(v) => Ok(*v),
        _ => Err(Error::InvalidQuery(format!(
            "{} requires an integer literal",
            context
        ))),
    }
}

fn coerce_to_f64(value: &Literal, context: &str) -> Result<f64> {
    match value {
        Literal::Float(v) => Ok(*v),
        Literal::Integer(v) => Ok(*v as f64),
        Literal::Str(_) => Err(Error::InvalidQuery(format!(
            "{} requires a numeric literal",
            context
        ))),
    }
}

fn parse_variant_type(s: &str) -> Result<VariantType> {
    match s.to_ascii_uppercase().as_str() {
        "SNV" => Ok(VariantType::Snv),
        "INSERTION" => Ok(VariantType::Insertion),
        "DELETION" => Ok(VariantType::Deletion),
        "MNV" => Ok(VariantType::Mnv),
        "COMPLEX" => Ok(VariantType::Complex),
        _ => Err(Error::InvalidQuery(format!(
            "unknown variant type: '{}'. Valid values: SNV, Insertion, Deletion, Mnv, Complex",
            s
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse as parse_filter;

    fn test_header() -> VcfHeader {
        let input = b"##fileformat=VCFv4.3\n\
            ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
            ##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n\
            ##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP\">\n\
            ##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let (header, _) = VcfHeader::parse(input).unwrap();
        header
    }

    #[test]
    fn plan_qual_comparison() {
        let header = test_header();
        let expr = parse_filter("QUAL > 50").unwrap();
        let ops = plan(&expr, &header).unwrap();
        assert_eq!(ops, vec![OpCode::CmpQual(CmpOp::Gt, 50.0)]);
    }

    #[test]
    fn plan_chrom_comparison() {
        let header = test_header();
        let expr = parse_filter("CHROM = 'chr1'").unwrap();
        let ops = plan(&expr, &header).unwrap();
        assert_eq!(ops, vec![OpCode::CmpChrom(CmpOp::Eq, b"chr1".to_vec())]);
    }

    #[test]
    fn plan_info_int() {
        let header = test_header();
        let expr = parse_filter("INFO.DP > 10").unwrap();
        let ops = plan(&expr, &header).unwrap();
        assert_eq!(ops, vec![OpCode::CmpInfoInt { key: b"DP".to_vec(), op: CmpOp::Gt, value: 10 }]);
    }

    #[test]
    fn plan_info_float_multi_value() {
        // AF is Number=A → should emit AnyCmpInfoFloat
        let header = test_header();
        let expr = parse_filter("INFO.AF > 0.05").unwrap();
        let ops = plan(&expr, &header).unwrap();
        assert_eq!(ops, vec![OpCode::AnyCmpInfoFloat { key: b"AF".to_vec(), op: CmpOp::Gt, value: 0.05 }]);
    }

    #[test]
    fn plan_info_flag() {
        let header = test_header();
        let expr = parse_filter("INFO.DB").unwrap();
        let ops = plan(&expr, &header).unwrap();
        assert_eq!(ops, vec![OpCode::InfoFlag(b"DB".to_vec())]);
    }

    #[test]
    fn plan_type_comparison() {
        let header = test_header();
        let expr = parse_filter("TYPE = 'SNV'").unwrap();
        let ops = plan(&expr, &header).unwrap();
        assert_eq!(ops, vec![OpCode::CmpType(CmpOp::Eq, VariantType::Snv)]);
    }

    #[test]
    fn plan_and_postfix() {
        let header = test_header();
        let expr = parse_filter("QUAL > 50 AND CHROM = 'chr1'").unwrap();
        let ops = plan(&expr, &header).unwrap();
        // Postfix with jump: [CmpQual, JumpIfFalse(4), CmpChrom, And]
        assert_eq!(ops.len(), 4);
        assert!(matches!(ops[0], OpCode::CmpQual(..)));
        assert!(matches!(ops[1], OpCode::JumpIfFalse(4)));
        assert!(matches!(ops[2], OpCode::CmpChrom(..)));
        assert!(matches!(ops[3], OpCode::And));
    }

    #[test]
    fn plan_compound_postfix() {
        let header = test_header();
        let expr = parse_filter("CHROM = 'chr1' OR QUAL > 50 AND POS > 100").unwrap();
        let ops = plan(&expr, &header).unwrap();
        // AST: Or(CmpChrom, And(CmpQual, CmpPos))
        // With jumps: [CmpChrom, JumpIfTrue(7), CmpQual, JumpIfFalse(6), CmpPos, And, Or]
        assert_eq!(ops.len(), 7);
        assert!(matches!(ops[0], OpCode::CmpChrom(..)));
        assert!(matches!(ops[1], OpCode::JumpIfTrue(7)));
        assert!(matches!(ops[3], OpCode::JumpIfFalse(6)));
        assert!(matches!(ops[5], OpCode::And));
        assert!(matches!(ops[6], OpCode::Or));
    }

    #[test]
    fn plan_not() {
        let header = test_header();
        let expr = parse_filter("NOT INFO.DB").unwrap();
        let ops = plan(&expr, &header).unwrap();
        assert_eq!(ops, vec![OpCode::InfoFlag(b"DB".to_vec()), OpCode::Not]);
    }

    #[test]
    fn error_info_missing() {
        let header = test_header();
        let expr = parse_filter("INFO.NONEXISTENT > 5").unwrap();
        assert!(plan(&expr, &header).is_err());
    }

    #[test]
    fn error_flag_with_comparison() {
        let header = test_header();
        let expr = parse_filter("INFO.DB = 1").unwrap();
        assert!(plan(&expr, &header).is_err());
    }

    #[test]
    fn error_int_field_with_string() {
        let header = test_header();
        let expr = parse_filter("INFO.DP = 'foo'").unwrap();
        assert!(plan(&expr, &header).is_err());
    }

    #[test]
    fn plan_scalar_vs_multi_value() {
        let header = test_header();
        // DP is Number=1 → scalar CmpInfoInt
        let expr = parse_filter("INFO.DP > 10").unwrap();
        let ops = plan(&expr, &header).unwrap();
        assert!(matches!(ops[0], OpCode::CmpInfoInt { .. }));

        // AF is Number=A → multi-value AnyCmpInfoFloat
        let expr = parse_filter("INFO.AF > 0.5").unwrap();
        let ops = plan(&expr, &header).unwrap();
        assert!(matches!(ops[0], OpCode::AnyCmpInfoFloat { .. }));
    }

    #[test]
    fn error_unknown_variant_type() {
        let header = test_header();
        let expr = parse_filter("TYPE = 'UNKNOWN'").unwrap();
        assert!(plan(&expr, &header).is_err());
    }
}
