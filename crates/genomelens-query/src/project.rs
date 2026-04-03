use std::io::Write;

use genomelens_core::Result;
use genomelens_parse::VcfRecord;

use crate::ast::{Column, FixedColumn};

/// Write a column value directly to the output — zero allocation.
pub fn write_column(column: &Column, record: &VcfRecord<'_>, out: &mut impl Write) -> Result<()> {
    match column {
        Column::Fixed(f) => write_fixed(f, record, out),
        Column::Info(field) => write_info(field, record, out),
    }
}

fn write_fixed(col: &FixedColumn, record: &VcfRecord<'_>, out: &mut impl Write) -> Result<()> {
    match col {
        FixedColumn::Chrom => out.write_all(record.chrom()?)?,
        FixedColumn::Pos => out.write_all(record.pos()?)?,
        FixedColumn::Qual => out.write_all(record.qual()?)?,
        FixedColumn::Filter => out.write_all(record.filter()?)?,
        FixedColumn::Ref => out.write_all(record.ref_allele()?)?,
        FixedColumn::Type => write!(out, "{:?}", record.variant_type()?)?,
    }
    Ok(())
}

fn write_info(field: &str, record: &VcfRecord<'_>, out: &mut impl Write) -> Result<()> {
    match record.info_value(field.as_bytes())? {
        Some(Some(val)) => out.write_all(val)?,
        Some(None) => out.write_all(b"true")?,
        None => out.write_all(b".")?,
    }
    Ok(())
}

/// Write a TSV header line directly to the output — zero allocation.
pub fn write_column_header(columns: &[Column], out: &mut impl Write) -> std::io::Result<()> {
    for (i, col) in columns.iter().enumerate() {
        if i > 0 {
            out.write_all(b"\t")?;
        }
        match col {
            Column::Fixed(FixedColumn::Chrom) => out.write_all(b"CHROM")?,
            Column::Fixed(FixedColumn::Pos) => out.write_all(b"POS")?,
            Column::Fixed(FixedColumn::Qual) => out.write_all(b"QUAL")?,
            Column::Fixed(FixedColumn::Filter) => out.write_all(b"FILTER")?,
            Column::Fixed(FixedColumn::Ref) => out.write_all(b"REF")?,
            Column::Fixed(FixedColumn::Type) => out.write_all(b"TYPE")?,
            Column::Info(name) => {
                out.write_all(b"INFO.")?;
                out.write_all(name.as_bytes())?;
            }
        }
    }
    Ok(())
}

/// Extract column value as owned bytes. Thin wrapper over `write_column` for tests.
pub fn extract_column(column: &Column, record: &VcfRecord<'_>) -> Result<Vec<u8>> {
    let mut buf = Vec::new();
    write_column(column, record, &mut buf)?;
    Ok(buf)
}

/// Build TSV header as owned string. Thin wrapper over `write_column_header` for tests.
pub fn column_header(columns: &[Column]) -> String {
    let mut buf = Vec::new();
    write_column_header(columns, &mut buf).expect("writing to Vec cannot fail");
    String::from_utf8(buf).expect("header is always valid UTF-8")
}

#[cfg(test)]
mod tests {
    use super::*;
    use genomelens_parse::vcf::record::parse_line;

    fn make_record(line: &[u8]) -> genomelens_parse::VcfRecord<'_> {
        parse_line(line)
    }

    #[test]
    fn extract_chrom() {
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=10");
        let val = extract_column(&Column::Fixed(FixedColumn::Chrom), &rec).unwrap();
        assert_eq!(val, b"chr1");
    }

    #[test]
    fn extract_pos() {
        let rec = make_record(b"chr1\t12345\t.\tA\tG\t50\tPASS\t.");
        let val = extract_column(&Column::Fixed(FixedColumn::Pos), &rec).unwrap();
        assert_eq!(val, b"12345");
    }

    #[test]
    fn extract_info_value() {
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=25;AF=0.3");
        let val = extract_column(&Column::Info("DP".to_string()), &rec).unwrap();
        assert_eq!(val, b"25");
    }

    #[test]
    fn extract_info_flag() {
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\tDB;DP=10");
        let val = extract_column(&Column::Info("DB".to_string()), &rec).unwrap();
        assert_eq!(val, b"true");
    }

    #[test]
    fn extract_info_missing() {
        let rec = make_record(b"chr1\t100\t.\tA\tG\t50\tPASS\t.");
        let val = extract_column(&Column::Info("DP".to_string()), &rec).unwrap();
        assert_eq!(val, b".");
    }

    #[test]
    fn column_header_formatting() {
        let cols = vec![
            Column::Fixed(FixedColumn::Chrom),
            Column::Fixed(FixedColumn::Pos),
            Column::Info("AF".to_string()),
        ];
        assert_eq!(column_header(&cols), "CHROM\tPOS\tINFO.AF");
    }
}
