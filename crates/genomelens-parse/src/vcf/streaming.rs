use genomelens_core::{Error, Result};
use std::io::BufRead;

use super::header::VcfHeader;
use super::record::{parse_line, VcfRecord};

/// Streaming VCF reader using the lending iterator pattern.
///
/// Owns an internal line buffer. Each call to `next_record()` yields a
/// `VcfRecord<'_>` that borrows from this buffer. The previous record is
/// invalidated when `next_record()` is called again.
///
/// This enables processing arbitrarily large VCF files without loading
/// them into memory.
pub struct VcfReader<R> {
    reader: R,
    line_buf: Vec<u8>,
}

impl<R: BufRead> VcfReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_buf: Vec::with_capacity(4096),
        }
    }

    /// Read and parse all header lines (## and #CHROM).
    /// After this call, the reader is positioned at the first data line.
    pub fn read_header(&mut self) -> Result<VcfHeader> {
        VcfHeader::parse_from_reader(&mut self.reader)
    }

    /// Yield the next VCF record, borrowing from the internal line buffer.
    ///
    /// Returns `Ok(None)` at EOF. Each call invalidates the previous record.
    /// The caller must process the record before calling `next_record()` again.
    pub fn next_record(&mut self) -> Result<Option<VcfRecord<'_>>> {
        loop {
            self.line_buf.clear();
            let bytes_read = self
                .reader
                .read_until(b'\n', &mut self.line_buf)
                .map_err(Error::Io)?;

            if bytes_read == 0 {
                return Ok(None);
            }

            // Trim trailing \n and \r
            let mut end = self.line_buf.len();
            if end > 0 && self.line_buf[end - 1] == b'\n' {
                end -= 1;
            }
            if end > 0 && self.line_buf[end - 1] == b'\r' {
                end -= 1;
            }

            // Skip empty lines
            if end == 0 {
                continue;
            }

            let line = &self.line_buf[..end];
            return Ok(Some(parse_line(line)));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn make_vcf_bytes(data: &str) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(
            b"##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        );
        buf.extend_from_slice(data.as_bytes());
        buf
    }

    #[test]
    fn streaming_basic() {
        let data = make_vcf_bytes(
            "chr1\t100\t.\tA\tG\t50\tPASS\t.\n\
             chr1\t200\t.\tC\tT\t60\tPASS\t.\n",
        );
        let mut reader = VcfReader::new(Cursor::new(data));
        let header = reader.read_header().unwrap();
        assert_eq!(header.file_format, "VCFv4.3");

        let rec1 = reader.next_record().unwrap().unwrap();
        assert_eq!(rec1.chrom().unwrap(), b"chr1");
        assert_eq!(rec1.pos().unwrap(), b"100");

        let rec2 = reader.next_record().unwrap().unwrap();
        assert_eq!(rec2.chrom().unwrap(), b"chr1");
        assert_eq!(rec2.pos().unwrap(), b"200");

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn streaming_empty_data() {
        let data = make_vcf_bytes("");
        let mut reader = VcfReader::new(Cursor::new(data));
        reader.read_header().unwrap();
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn streaming_skips_empty_lines() {
        let data = make_vcf_bytes(
            "chr1\t100\t.\tA\tG\t50\tPASS\t.\n\n\nchr1\t200\t.\tC\tT\t60\tPASS\t.\n",
        );
        let mut reader = VcfReader::new(Cursor::new(data));
        reader.read_header().unwrap();

        let rec1 = reader.next_record().unwrap().unwrap();
        assert_eq!(rec1.pos().unwrap(), b"100");

        let rec2 = reader.next_record().unwrap().unwrap();
        assert_eq!(rec2.pos().unwrap(), b"200");

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn streaming_with_samples() {
        // Need a header with FORMAT + samples
        let full = b"##fileformat=VCFv4.3\n\
            ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
            chr1\t100\t.\tA\tG\t50\tPASS\tDP=10\tGT\t0/1\t1/1\n";

        let mut reader = VcfReader::new(Cursor::new(full.as_slice()));
        let header = reader.read_header().unwrap();
        assert_eq!(header.sample_names, vec!["S1", "S2"]);

        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.format().unwrap(), Some(b"GT".as_slice()));
        let samples: Vec<_> = rec.sample_iter().unwrap().collect();
        assert_eq!(samples.len(), 2);
    }

    #[test]
    fn streaming_windows_line_endings() {
        let data = b"##fileformat=VCFv4.3\r\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\r\nchr1\t100\t.\tA\tG\t50\tPASS\t.\r\n";
        let mut reader = VcfReader::new(Cursor::new(data.as_slice()));
        let header = reader.read_header().unwrap();
        assert_eq!(header.file_format, "VCFv4.3");

        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.chrom().unwrap(), b"chr1");
    }

    #[test]
    fn streaming_lazy_info() {
        let data = make_vcf_bytes("chr1\t100\t.\tA\tG\t50\tPASS\tDP=14;AF=0.5;DB\n");
        let mut reader = VcfReader::new(Cursor::new(data));
        reader.read_header().unwrap();

        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.info_value(b"DP").unwrap(), Some(Some(b"14".as_slice())));
        assert_eq!(rec.info_value(b"AF").unwrap(), Some(Some(b"0.5".as_slice())));
        assert_eq!(rec.info_value(b"DB").unwrap(), Some(None)); // flag
        assert_eq!(rec.info_value(b"MISSING").unwrap(), None);
    }

    #[test]
    fn streaming_count() {
        let data = make_vcf_bytes(
            "chr1\t100\t.\tA\tG\t50\tPASS\t.\n\
             chr1\t200\t.\tC\tT\t60\tPASS\t.\n\
             chr2\t300\t.\tG\tA\t70\tPASS\t.\n\
             chr2\t400\t.\tT\tC\t80\tPASS\t.\n\
             chr3\t500\t.\tA\tG\t90\tPASS\t.\n",
        );
        let mut reader = VcfReader::new(Cursor::new(data));
        reader.read_header().unwrap();

        let mut count = 0;
        while reader.next_record().unwrap().is_some() {
            count += 1;
        }
        assert_eq!(count, 5);
    }
}
