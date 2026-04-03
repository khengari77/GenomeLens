use genomelens_core::{
    ContigDef, Error, FilterDef, FormatFieldDef, InfoFieldDef, Result, VcfNumber, VcfValueType,
};
use memchr::memchr;

/// Parsed VCF header. Owns all extracted metadata.
#[derive(Debug, Clone)]
pub struct VcfHeader {
    pub file_format: String,
    pub info_fields: Vec<InfoFieldDef>,
    pub format_fields: Vec<FormatFieldDef>,
    pub filters: Vec<FilterDef>,
    pub contigs: Vec<ContigDef>,
    pub sample_names: Vec<String>,
    pub raw_meta_lines: Vec<String>,
    pub chrom_line: String,
}

impl VcfHeader {
    /// Parse header from buffer. Returns (header, data_start_offset).
    /// `data_start_offset` points to the first byte after the `#CHROM` header line.
    pub fn parse(buf: &[u8]) -> Result<(Self, usize)> {
        let mut header = VcfHeader {
            file_format: String::new(),
            info_fields: Vec::new(),
            format_fields: Vec::new(),
            filters: Vec::new(),
            contigs: Vec::new(),
            sample_names: Vec::new(),
            raw_meta_lines: Vec::new(),
            chrom_line: String::new(),
        };

        let mut pos = 0;

        loop {
            if pos >= buf.len() {
                return Err(Error::InvalidVcf(
                    "missing #CHROM header line".to_string(),
                ));
            }

            // Find end of current line
            let line_end = match memchr(b'\n', &buf[pos..]) {
                Some(offset) => pos + offset,
                None => buf.len(),
            };

            let line = trim_cr(&buf[pos..line_end]);

            if line.starts_with(b"##") {
                header.parse_meta_line(line)?;
                pos = if line_end < buf.len() {
                    line_end + 1
                } else {
                    buf.len()
                };
            } else if line.starts_with(b"#CHROM") {
                header.chrom_line = String::from_utf8_lossy(line).into_owned();
                header.parse_header_line(line)?;
                pos = if line_end < buf.len() {
                    line_end + 1
                } else {
                    buf.len()
                };
                break;
            } else {
                return Err(Error::InvalidVcf(format!(
                    "expected ## or #CHROM line, found: {}",
                    String::from_utf8_lossy(&line[..line.len().min(40)])
                )));
            }
        }

        Ok((header, pos))
    }

    /// Parse header from a streaming BufRead source. Reads lines one at a time.
    /// After this call, the reader is positioned at the first data line.
    pub fn parse_from_reader(reader: &mut impl std::io::BufRead) -> Result<Self> {
        let mut header = VcfHeader {
            file_format: String::new(),
            info_fields: Vec::new(),
            format_fields: Vec::new(),
            filters: Vec::new(),
            contigs: Vec::new(),
            sample_names: Vec::new(),
            raw_meta_lines: Vec::new(),
            chrom_line: String::new(),
        };

        let mut line_buf = Vec::with_capacity(1024);

        loop {
            line_buf.clear();
            let bytes_read = reader.read_until(b'\n', &mut line_buf).map_err(Error::Io)?;
            if bytes_read == 0 {
                return Err(Error::InvalidVcf(
                    "missing #CHROM header line".to_string(),
                ));
            }

            let line = trim_cr(trim_newline(&line_buf));

            if line.starts_with(b"##") {
                header.parse_meta_line(line)?;
            } else if line.starts_with(b"#CHROM") {
                header.chrom_line = String::from_utf8_lossy(line).into_owned();
                header.parse_header_line(line)?;
                break;
            } else if line.is_empty() {
                continue;
            } else {
                return Err(Error::InvalidVcf(format!(
                    "expected ## or #CHROM line, found: {}",
                    String::from_utf8_lossy(&line[..line.len().min(40)])
                )));
            }
        }

        Ok(header)
    }

    /// Look up an INFO field definition by ID.
    pub fn info_field(&self, id: &str) -> Option<&InfoFieldDef> {
        self.info_fields.iter().find(|f| f.id == id)
    }

    /// Look up a FORMAT field definition by ID.
    pub fn format_field(&self, id: &str) -> Option<&FormatFieldDef> {
        self.format_fields.iter().find(|f| f.id == id)
    }

    fn parse_meta_line(&mut self, line: &[u8]) -> Result<()> {
        // Strip leading "##"
        let content = &line[2..];

        if content.starts_with(b"fileformat=") {
            self.file_format =
                String::from_utf8_lossy(&content[b"fileformat=".len()..]).into_owned();
        } else if content.starts_with(b"INFO=<") {
            let inner = extract_angle_bracket(content, b"INFO=<")?;
            let fields = parse_angle_bracket_fields(inner)?;
            self.info_fields.push(InfoFieldDef {
                id: get_field(&fields, b"ID")?,
                number: parse_vcf_number(get_field_bytes(&fields, b"Number")?)?,
                ty: parse_vcf_type(get_field_bytes(&fields, b"Type")?)?,
                description: get_field(&fields, b"Description")?,
            });
        } else if content.starts_with(b"FORMAT=<") {
            let inner = extract_angle_bracket(content, b"FORMAT=<")?;
            let fields = parse_angle_bracket_fields(inner)?;
            self.format_fields.push(FormatFieldDef {
                id: get_field(&fields, b"ID")?,
                number: parse_vcf_number(get_field_bytes(&fields, b"Number")?)?,
                ty: parse_vcf_type(get_field_bytes(&fields, b"Type")?)?,
                description: get_field(&fields, b"Description")?,
            });
        } else if content.starts_with(b"FILTER=<") {
            let inner = extract_angle_bracket(content, b"FILTER=<")?;
            let fields = parse_angle_bracket_fields(inner)?;
            self.filters.push(FilterDef {
                id: get_field(&fields, b"ID")?,
                description: get_field(&fields, b"Description")?,
            });
        } else if content.starts_with(b"contig=<") {
            let inner = extract_angle_bracket(content, b"contig=<")?;
            let fields = parse_angle_bracket_fields(inner)?;
            let id = get_field(&fields, b"ID")?;
            let length = match get_field_bytes_opt(&fields, b"length") {
                Some(bytes) => {
                    let s = std::str::from_utf8(bytes).map_err(|_| {
                        Error::InvalidVcf("invalid contig length encoding".to_string())
                    })?;
                    Some(s.parse::<usize>().map_err(|_| {
                        Error::InvalidVcf(format!("invalid contig length: {s}"))
                    })?)
                }
                None => None,
            };
            self.contigs.push(ContigDef { id, length });
        }

        self.raw_meta_lines
            .push(String::from_utf8_lossy(line).into_owned());
        Ok(())
    }

    fn parse_header_line(&mut self, line: &[u8]) -> Result<()> {
        // Split on tabs
        let cols: Vec<&[u8]> = line.split(|&b| b == b'\t').collect();
        // Columns 0-8 are fixed: #CHROM POS ID REF ALT QUAL FILTER INFO (FORMAT)
        // Columns 9+ are sample names
        if cols.len() < 8 {
            return Err(Error::InvalidVcf(format!(
                "expected at least 8 columns in #CHROM line, found {}",
                cols.len()
            )));
        }

        // Sample names start at column 9 (if FORMAT column exists at 8)
        if cols.len() > 9 {
            for col in &cols[9..] {
                self.sample_names
                    .push(String::from_utf8_lossy(col).into_owned());
            }
        }

        Ok(())
    }
}

/// Trim trailing `\n` from a line.
fn trim_newline(line: &[u8]) -> &[u8] {
    if line.last() == Some(&b'\n') {
        &line[..line.len() - 1]
    } else {
        line
    }
}

/// Trim trailing `\r` from a line.
fn trim_cr(line: &[u8]) -> &[u8] {
    if line.last() == Some(&b'\r') {
        &line[..line.len() - 1]
    } else {
        line
    }
}

/// Extract content between `<` and `>` from a prefix like `INFO=<...>`.
fn extract_angle_bracket<'a>(content: &'a [u8], prefix: &[u8]) -> Result<&'a [u8]> {
    let inner = &content[prefix.len()..];
    // Find closing `>`
    if inner.last() == Some(&b'>') {
        Ok(&inner[..inner.len() - 1])
    } else {
        Err(Error::InvalidVcf(
            "unclosed angle bracket in header line".to_string(),
        ))
    }
}

/// Parse structured VCF meta-line fields: `ID=DP,Number=1,Type=Integer,Description="Total Depth"`.
/// Handles quoted values containing commas.
fn parse_angle_bracket_fields(inner: &[u8]) -> Result<Vec<(&[u8], &[u8])>> {
    let mut fields = Vec::new();
    let mut pos = 0;

    while pos < inner.len() {
        // Find the `=` separator for key
        let eq_pos = match memchr(b'=', &inner[pos..]) {
            Some(offset) => pos + offset,
            None => break,
        };

        let key = &inner[pos..eq_pos];
        let value_start = eq_pos + 1;

        // Check if value is quoted
        let (value, next_pos) = if value_start < inner.len() && inner[value_start] == b'"' {
            // Quoted value: scan for closing quote
            let quote_content_start = value_start + 1;
            let mut end = quote_content_start;
            while end < inner.len() && inner[end] != b'"' {
                end += 1;
            }
            let value = &inner[quote_content_start..end];
            // Skip past closing quote and comma
            let next = if end + 1 < inner.len() && inner[end + 1] == b',' {
                end + 2
            } else {
                end + 1
            };
            (value, next)
        } else {
            // Unquoted value: scan for comma
            let end = match memchr(b',', &inner[value_start..]) {
                Some(offset) => value_start + offset,
                None => inner.len(),
            };
            let value = &inner[value_start..end];
            let next = if end < inner.len() { end + 1 } else { end };
            (value, next)
        };

        fields.push((key, value));
        pos = next_pos;
    }

    Ok(fields)
}

fn get_field_bytes_opt<'a>(fields: &[(&[u8], &'a [u8])], key: &[u8]) -> Option<&'a [u8]> {
    fields.iter().find(|(k, _)| *k == key).map(|(_, v)| *v)
}

fn get_field_bytes<'a>(fields: &[(&[u8], &'a [u8])], key: &[u8]) -> Result<&'a [u8]> {
    get_field_bytes_opt(fields, key).ok_or_else(|| {
        Error::InvalidVcf(format!(
            "missing required field: {}",
            String::from_utf8_lossy(key)
        ))
    })
}

fn get_field(fields: &[(&[u8], &[u8])], key: &[u8]) -> Result<String> {
    get_field_bytes(fields, key).map(|v| String::from_utf8_lossy(v).into_owned())
}

fn parse_vcf_type(s: &[u8]) -> Result<VcfValueType> {
    match s {
        b"Integer" => Ok(VcfValueType::Integer),
        b"Float" => Ok(VcfValueType::Float),
        b"Flag" => Ok(VcfValueType::Flag),
        b"Character" => Ok(VcfValueType::Character),
        b"String" => Ok(VcfValueType::String),
        other => Err(Error::InvalidVcf(format!(
            "unknown VCF type: {}",
            String::from_utf8_lossy(other)
        ))),
    }
}

fn parse_vcf_number(s: &[u8]) -> Result<VcfNumber> {
    match s {
        b"A" => Ok(VcfNumber::PerAltAllele),
        b"R" => Ok(VcfNumber::PerAllele),
        b"G" => Ok(VcfNumber::PerGenotype),
        b"." => Ok(VcfNumber::Unbounded),
        digits => {
            let n = std::str::from_utf8(digits)
                .map_err(|_| Error::InvalidVcf("invalid Number encoding".to_string()))?
                .parse::<usize>()
                .map_err(|_| {
                    Error::InvalidVcf(format!(
                        "invalid Number value: {}",
                        String::from_utf8_lossy(digits)
                    ))
                })?;
            Ok(VcfNumber::Count(n))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINIMAL_HEADER: &[u8] = b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    #[test]
    fn parse_minimal_header() {
        let (header, offset) = VcfHeader::parse(MINIMAL_HEADER).unwrap();
        assert_eq!(header.file_format, "VCFv4.3");
        assert!(header.info_fields.is_empty());
        assert!(header.sample_names.is_empty());
        assert_eq!(offset, MINIMAL_HEADER.len());
    }

    #[test]
    fn parse_full_header() {
        let input = b"##fileformat=VCFv4.3\n\
            ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n\
            ##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n\
            ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
            ##FILTER=<ID=q10,Description=\"Quality below 10\">\n\
            ##contig=<ID=chr1,length=248956422>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n";

        let (header, _offset) = VcfHeader::parse(input).unwrap();

        assert_eq!(header.info_fields.len(), 2);
        assert_eq!(header.info_fields[0].id, "DP");
        assert_eq!(header.info_fields[0].ty, VcfValueType::Integer);
        assert_eq!(header.info_fields[0].number, VcfNumber::Count(1));

        assert_eq!(header.info_fields[1].id, "AF");
        assert_eq!(header.info_fields[1].ty, VcfValueType::Float);
        assert_eq!(header.info_fields[1].number, VcfNumber::PerAltAllele);

        assert_eq!(header.format_fields.len(), 1);
        assert_eq!(header.format_fields[0].id, "GT");

        assert_eq!(header.filters.len(), 1);
        assert_eq!(header.filters[0].id, "q10");

        assert_eq!(header.contigs.len(), 1);
        assert_eq!(header.contigs[0].id, "chr1");
        assert_eq!(header.contigs[0].length, Some(248956422));

        assert_eq!(header.sample_names, vec!["Sample1", "Sample2"]);
    }

    #[test]
    fn parse_quoted_commas_in_description() {
        let input = b"##fileformat=VCFv4.3\n\
            ##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence, as predicted by VEP\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

        let (header, _) = VcfHeader::parse(input).unwrap();
        assert_eq!(header.info_fields[0].id, "CSQ");
        assert_eq!(
            header.info_fields[0].description,
            "Consequence, as predicted by VEP"
        );
    }

    #[test]
    fn missing_chrom_line() {
        let input = b"##fileformat=VCFv4.3\n";
        let result = VcfHeader::parse(input);
        assert!(result.is_err());
    }

    #[test]
    fn no_sample_vcf() {
        let input = b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let (header, _) = VcfHeader::parse(input).unwrap();
        assert!(header.sample_names.is_empty());
    }

    #[test]
    fn data_start_offset_correct() {
        let input = b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t50\tPASS\t.\n";
        let (_, offset) = VcfHeader::parse(input).unwrap();
        assert_eq!(&input[offset..offset + 4], b"chr1");
    }

    #[test]
    fn windows_line_endings() {
        let input = b"##fileformat=VCFv4.3\r\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\r\n";
        let (header, _) = VcfHeader::parse(input).unwrap();
        assert_eq!(header.file_format, "VCFv4.3");
    }

    #[test]
    fn info_field_lookup() {
        let input = b"##fileformat=VCFv4.3\n\
            ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

        let (header, _) = VcfHeader::parse(input).unwrap();
        assert!(header.info_field("DP").is_some());
        assert!(header.info_field("NONEXISTENT").is_none());
    }

    #[test]
    fn contig_without_length() {
        let input = b"##fileformat=VCFv4.3\n\
            ##contig=<ID=chrM>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

        let (header, _) = VcfHeader::parse(input).unwrap();
        assert_eq!(header.contigs[0].id, "chrM");
        assert_eq!(header.contigs[0].length, None);
    }
}
