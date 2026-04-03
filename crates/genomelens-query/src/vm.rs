use genomelens_core::{Result, VariantType};
use genomelens_parse::VcfRecord;

use crate::ast::CmpOp;

/// Flat bytecode instruction. Each leaf pushes a bool onto the stack.
/// Boolean operators pop operands and push the result.
#[derive(Debug, Clone, PartialEq)]
pub enum OpCode {
    CmpChrom(CmpOp, Vec<u8>),
    CmpPos(CmpOp, usize),
    CmpQual(CmpOp, f64),
    CmpFilter(CmpOp, Vec<u8>),
    CmpRef(CmpOp, Vec<u8>),
    CmpType(CmpOp, VariantType),
    CmpInfoInt { key: Vec<u8>, op: CmpOp, value: i64 },
    CmpInfoFloat { key: Vec<u8>, op: CmpOp, value: f64 },
    CmpInfoStr { key: Vec<u8>, op: CmpOp, value: Vec<u8> },
    /// Push true if the INFO flag key is present on the record.
    /// Missing fields evaluate to false (strict boolean, not SQL ternary NULL).
    InfoFlag(Vec<u8>),
    And,
    Or,
    Not,
    /// Jump to `target` instruction if stack top is false (peek, don't pop).
    JumpIfFalse(usize),
    /// Jump to `target` instruction if stack top is true (peek, don't pop).
    JumpIfTrue(usize),
}

/// Stack-machine evaluator. Instructions are stored in a contiguous `Vec`
/// (cache-friendly), and the evaluation stack is a `Vec<bool>` allocated once
/// at construction and reused across rows via `clear()`.
pub struct Vm {
    ops: Vec<OpCode>,
    stack: Vec<bool>,
}

impl Vm {
    pub fn new(ops: Vec<OpCode>) -> Self {
        Self {
            ops,
            stack: Vec::with_capacity(64),
        }
    }

    /// Evaluate the compiled filter against a single VCF record.
    ///
    /// Missing fields evaluate to `false` (strict boolean semantics).
    /// This means `NOT INFO.DB` returns `true` when DB is absent — this is
    /// intentional closed-world semantics, not SQL ternary NULL logic.
    #[inline]
    pub fn evaluate(&mut self, record: &VcfRecord<'_>) -> Result<bool> {
        self.stack.clear();
        let mut ip: usize = 0;

        while ip < self.ops.len() {
            match &self.ops[ip] {
                OpCode::CmpChrom(cmp, val) => {
                    self.stack.push(cmp_bytes(record.chrom()?, val, *cmp));
                }
                OpCode::CmpPos(cmp, val) => {
                    let pos = record.pos_usize()?;
                    self.stack.push(cmp_ord(&pos, val, *cmp));
                }
                OpCode::CmpQual(cmp, val) => {
                    let result = match record.qual_f64()? {
                        Some(q) => cmp_f64(q, *val, *cmp),
                        None => false,
                    };
                    self.stack.push(result);
                }
                OpCode::CmpFilter(cmp, val) => {
                    self.stack.push(cmp_bytes(record.filter()?, val, *cmp));
                }
                OpCode::CmpRef(cmp, val) => {
                    self.stack.push(cmp_bytes(record.ref_allele()?, val, *cmp));
                }
                OpCode::CmpType(cmp, val) => {
                    self.stack.push(cmp_eq_only(&record.variant_type()?, val, *cmp));
                }
                OpCode::CmpInfoInt { key, op, value } => {
                    let result = match record.info_value(key)? {
                        Some(Some(raw)) => {
                            let v = genomelens_parse::ascii::parse_ascii_i64(raw)?;
                            cmp_ord(&v, value, *op)
                        }
                        _ => false,
                    };
                    self.stack.push(result);
                }
                OpCode::CmpInfoFloat { key, op, value } => {
                    let result = match record.info_value(key)? {
                        Some(Some(raw)) => {
                            let v = genomelens_parse::ascii::parse_ascii_f64(raw)?;
                            cmp_f64(v, *value, *op)
                        }
                        _ => false,
                    };
                    self.stack.push(result);
                }
                OpCode::CmpInfoStr { key, op, value } => {
                    let result = match record.info_value(key)? {
                        Some(Some(raw)) => cmp_bytes(raw, value, *op),
                        _ => false,
                    };
                    self.stack.push(result);
                }
                OpCode::InfoFlag(key) => {
                    self.stack.push(record.info_value(key)?.is_some());
                }
                OpCode::And => {
                    let r = self.stack.pop().unwrap_or(false);
                    let l = self.stack.last_mut().unwrap();
                    *l = *l && r;
                }
                OpCode::Or => {
                    let r = self.stack.pop().unwrap_or(false);
                    let l = self.stack.last_mut().unwrap();
                    *l = *l || r;
                }
                OpCode::Not => {
                    let top = self.stack.last_mut().unwrap();
                    *top = !*top;
                }
                OpCode::JumpIfFalse(target) => {
                    if !*self.stack.last().unwrap() {
                        ip = *target;
                        continue;
                    }
                }
                OpCode::JumpIfTrue(target) => {
                    if *self.stack.last().unwrap() {
                        ip = *target;
                        continue;
                    }
                }
            }
            ip += 1;
        }

        Ok(self.stack.last().copied().unwrap_or(false))
    }
}

/// Tolerance for float equality. Using 1e-6 instead of f64::EPSILON to
/// account for string-to-float parsing artifacts in VCF data.
const FLOAT_EQ_TOLERANCE: f64 = 1e-6;

#[inline]
fn cmp_bytes(a: &[u8], b: &[u8], op: CmpOp) -> bool {
    match op {
        CmpOp::Eq => a == b,
        CmpOp::Neq => a != b,
        CmpOp::Gt => a > b,
        CmpOp::Lt => a < b,
        CmpOp::Gte => a >= b,
        CmpOp::Lte => a <= b,
    }
}

#[inline]
fn cmp_ord<T: Ord>(a: &T, b: &T, op: CmpOp) -> bool {
    match op {
        CmpOp::Eq => a == b,
        CmpOp::Neq => a != b,
        CmpOp::Gt => a > b,
        CmpOp::Lt => a < b,
        CmpOp::Gte => a >= b,
        CmpOp::Lte => a <= b,
    }
}

#[inline]
fn cmp_f64(a: f64, b: f64, op: CmpOp) -> bool {
    match op {
        CmpOp::Eq => (a - b).abs() < FLOAT_EQ_TOLERANCE,
        CmpOp::Neq => (a - b).abs() >= FLOAT_EQ_TOLERANCE,
        CmpOp::Gt => a > b,
        CmpOp::Lt => a < b,
        CmpOp::Gte => a >= b,
        CmpOp::Lte => a <= b,
    }
}

#[inline]
fn cmp_eq_only<T: PartialEq>(a: &T, b: &T, op: CmpOp) -> bool {
    match op {
        CmpOp::Eq => a == b,
        CmpOp::Neq => a != b,
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use genomelens_parse::vcf::record::parse_line;

    fn make_record(line: &[u8]) -> genomelens_parse::VcfRecord<'_> {
        parse_line(line)
    }

    #[test]
    fn vm_chrom_eq() {
        let mut vm = Vm::new(vec![OpCode::CmpChrom(CmpOp::Eq, b"chr1".to_vec())]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\t.");
        assert!(vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_chrom_neq() {
        let mut vm = Vm::new(vec![OpCode::CmpChrom(CmpOp::Eq, b"chr1".to_vec())]);
        let rec = make_record(b"chr2\t100\t.\tA\tG\t50\tPASS\t.");
        assert!(!vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_pos_gt() {
        let mut vm = Vm::new(vec![OpCode::CmpPos(CmpOp::Gt, 100)]);
        let rec = make_record(b"chr1\t200\t.\tA\tG\t50\tPASS\t.");
        assert!(vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_qual_gt() {
        let mut vm = Vm::new(vec![OpCode::CmpQual(CmpOp::Gt, 50.0)]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t99.5\tPASS\t.");
        assert!(vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_qual_missing() {
        let mut vm = Vm::new(vec![OpCode::CmpQual(CmpOp::Gt, 50.0)]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t.\tPASS\t.");
        assert!(!vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_info_int() {
        let mut vm = Vm::new(vec![OpCode::CmpInfoInt {
            key: b"DP".to_vec(),
            op: CmpOp::Gt,
            value: 10,
        }]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=25");
        assert!(vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_info_float() {
        let mut vm = Vm::new(vec![OpCode::CmpInfoFloat {
            key: b"AF".to_vec(),
            op: CmpOp::Gt,
            value: 0.05,
        }]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\tAF=0.3");
        assert!(vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_info_missing() {
        let mut vm = Vm::new(vec![OpCode::CmpInfoInt {
            key: b"DP".to_vec(),
            op: CmpOp::Gt,
            value: 10,
        }]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\t.");
        assert!(!vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_info_flag() {
        let mut vm = Vm::new(vec![OpCode::InfoFlag(b"DB".to_vec())]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=10;DB");
        assert!(vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_info_flag_missing() {
        let mut vm = Vm::new(vec![OpCode::InfoFlag(b"DB".to_vec())]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=10");
        assert!(!vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_type_snv() {
        let mut vm = Vm::new(vec![OpCode::CmpType(CmpOp::Eq, VariantType::Snv)]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\t.");
        assert!(vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_and_postfix() {
        // CHROM = 'chr1' AND QUAL >= 50
        let mut vm = Vm::new(vec![
            OpCode::CmpChrom(CmpOp::Eq, b"chr1".to_vec()),
            OpCode::CmpQual(CmpOp::Gte, 50.0),
            OpCode::And,
        ]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\t.");
        assert!(vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_and_false() {
        let mut vm = Vm::new(vec![
            OpCode::CmpChrom(CmpOp::Eq, b"chr1".to_vec()),
            OpCode::CmpQual(CmpOp::Gt, 90.0),
            OpCode::And,
        ]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\t.");
        assert!(!vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_or_postfix() {
        let mut vm = Vm::new(vec![
            OpCode::CmpChrom(CmpOp::Eq, b"chr1".to_vec()),
            OpCode::CmpChrom(CmpOp::Eq, b"chr2".to_vec()),
            OpCode::Or,
        ]);
        let rec = make_record(b"chr2\t100\t.\tA\tG\t50\tPASS\t.");
        assert!(vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_not_postfix() {
        let mut vm = Vm::new(vec![
            OpCode::CmpChrom(CmpOp::Eq, b"chr2".to_vec()),
            OpCode::Not,
        ]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\t.");
        assert!(vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_complex_expression() {
        // (CHROM = 'chr1' OR CHROM = 'chr2') AND QUAL > 50
        // Postfix: [CmpChrom(chr1), CmpChrom(chr2), Or, CmpQual(>50), And]
        let mut vm = Vm::new(vec![
            OpCode::CmpChrom(CmpOp::Eq, b"chr1".to_vec()),
            OpCode::CmpChrom(CmpOp::Eq, b"chr2".to_vec()),
            OpCode::Or,
            OpCode::CmpQual(CmpOp::Gt, 50.0),
            OpCode::And,
        ]);
        let rec = make_record(b"chr2\t100\t.\tA\tG\t99\tPASS\t.");
        assert!(vm.evaluate(&rec).unwrap());

        let rec2 = make_record(b"chr3\t100\t.\tA\tG\t99\tPASS\t.");
        assert!(!vm.evaluate(&rec2).unwrap());
    }

    #[test]
    fn vm_empty_ops_returns_false() {
        let mut vm = Vm::new(vec![]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\t.");
        assert!(!vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_float_eq_tolerance() {
        // Test that float equality uses 1e-6 tolerance
        let mut vm = Vm::new(vec![OpCode::CmpQual(CmpOp::Eq, 50.0)]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\t.");
        assert!(vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_deep_nesting() {
        // 8 leaves chained with AND — exercises stack depth without overflow
        // (A AND B AND C AND D AND E AND F AND G AND H)
        // Postfix: [A, B, AND, C, AND, D, AND, E, AND, F, AND, G, AND, H, AND]
        let mut ops = vec![OpCode::CmpChrom(CmpOp::Eq, b"chr1".to_vec())];
        for _ in 0..7 {
            ops.push(OpCode::CmpQual(CmpOp::Gt, 10.0));
            ops.push(OpCode::And);
        }
        let mut vm = Vm::new(ops);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\t.");
        assert!(vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_not_flag_missing_is_true() {
        // NOT INFO.DB where DB is absent → true (closed-world assumption)
        let mut vm = Vm::new(vec![
            OpCode::InfoFlag(b"DB".to_vec()),
            OpCode::Not,
        ]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=10");
        assert!(vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_and_short_circuit() {
        // CHROM = 'chr2' AND INFO.DP > 10
        // With jumps: [CmpChrom(chr2), JumpIfFalse(4), CmpInfoInt(DP>10), And]
        // On chr1 record: CHROM fails → jump skips INFO parsing entirely
        let mut vm = Vm::new(vec![
            OpCode::CmpChrom(CmpOp::Eq, b"chr2".to_vec()),
            OpCode::JumpIfFalse(4),
            OpCode::CmpInfoInt { key: b"DP".to_vec(), op: CmpOp::Gt, value: 10 },
            OpCode::And,
        ]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=25");
        assert!(!vm.evaluate(&rec).unwrap());
    }

    #[test]
    fn vm_or_short_circuit() {
        // CHROM = 'chr1' OR INFO.DP > 10
        // With jumps: [CmpChrom(chr1), JumpIfTrue(4), CmpInfoInt(DP>10), Or]
        // On chr1 record: CHROM succeeds → jump skips INFO parsing entirely
        let mut vm = Vm::new(vec![
            OpCode::CmpChrom(CmpOp::Eq, b"chr1".to_vec()),
            OpCode::JumpIfTrue(4),
            OpCode::CmpInfoInt { key: b"DP".to_vec(), op: CmpOp::Gt, value: 10 },
            OpCode::Or,
        ]);
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=5");
        assert!(vm.evaluate(&rec).unwrap());
    }
}
