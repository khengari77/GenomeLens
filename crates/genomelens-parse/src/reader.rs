use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

/// Detected file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileFormat {
    Fasta,
    Vcf,
}

/// Check if a gzip file is BGZF by inspecting the header.
/// BGZF has FEXTRA flag (byte 3 bit 2) and "BC" subfield at bytes 12-13.
fn is_bgzf(path: &Path) -> io::Result<bool> {
    let mut file = File::open(path)?;
    let mut buf = [0u8; 18];
    let n = file.read(&mut buf)?;
    if n < 14 {
        return Ok(false);
    }
    Ok((buf[3] & 0x04) != 0 && buf[12] == b'B' && buf[13] == b'C')
}

/// Opens a file and returns a buffered reader that transparently decompresses gzip/BGZF.
/// Detection is by magic bytes (0x1f 0x8b), not file extension.
/// Uses flate2 MultiGzDecoder for sequential streaming (handles both gzip and BGZF).
/// For indexed random access, use `open_bgzf` instead.
pub fn open_transparent(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    let n = file.read(&mut magic)?;

    // Reopen to read from the beginning
    let file = File::open(path)?;

    if n >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
        // Gzip/BGZF — flate2 handles both for sequential streaming.
        Ok(Box::new(BufReader::new(
            flate2::read::MultiGzDecoder::new(file),
        )))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Detect file format by peeking at the first meaningful byte of a buffered stream.
/// Does NOT consume any bytes — the stream remains at its original position.
pub fn detect_format(reader: &mut impl BufRead) -> genomelens_core::Result<FileFormat> {
    let buf = reader.fill_buf().map_err(genomelens_core::Error::Io)?;
    let first = buf
        .iter()
        .find(|&&b| !matches!(b, b'\n' | b'\r' | b' ' | b'\t'));
    match first {
        Some(b'>') => Ok(FileFormat::Fasta),
        Some(b'#') => Ok(FileFormat::Vcf),
        Some(b) => Err(genomelens_core::Error::Parse {
            offset: 0,
            message: format!("unrecognized format: first byte is {:?}", *b as char),
        }),
        None => Err(genomelens_core::Error::Parse {
            offset: 0,
            message: "empty file".to_string(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn open_plain_text() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.fasta");
        std::fs::write(&path, b">seq1\nACGT\n").unwrap();

        let mut reader = open_transparent(&path).unwrap();
        let mut data = Vec::new();
        reader.read_to_end(&mut data).unwrap();
        assert_eq!(data, b">seq1\nACGT\n");
    }

    #[test]
    fn open_gzip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.fasta.gz");

        let file = File::create(&path).unwrap();
        let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        gz.write_all(b">seq1\nACGT\n").unwrap();
        gz.finish().unwrap();

        let mut reader = open_transparent(&path).unwrap();
        let mut data = Vec::new();
        reader.read_to_end(&mut data).unwrap();
        assert_eq!(data, b">seq1\nACGT\n");
    }

    #[test]
    fn detect_format_fasta() {
        let data = b">seq1\nACGT\n";
        let mut cursor = io::Cursor::new(data.as_slice());
        assert_eq!(detect_format(&mut cursor).unwrap(), FileFormat::Fasta);
    }

    #[test]
    fn detect_format_vcf() {
        let data = b"##fileformat=VCFv4.3\n";
        let mut cursor = io::Cursor::new(data.as_slice());
        assert_eq!(detect_format(&mut cursor).unwrap(), FileFormat::Vcf);
    }

    #[test]
    fn detect_format_does_not_consume() {
        let data = b"##fileformat=VCFv4.3\n";
        let mut cursor = io::Cursor::new(data.as_slice());
        detect_format(&mut cursor).unwrap();
        // Should still be able to read from position 0
        let mut buf = String::new();
        cursor.read_line(&mut buf).unwrap();
        assert_eq!(buf, "##fileformat=VCFv4.3\n");
    }

    #[test]
    fn detect_format_empty() {
        let data = b"";
        let mut cursor = io::Cursor::new(data.as_slice());
        assert!(detect_format(&mut cursor).is_err());
    }

    #[test]
    fn open_transparent_returns_bufread() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.fasta");
        std::fs::write(&path, b">seq1\nACGT\n").unwrap();

        let mut reader = open_transparent(&path).unwrap();
        let mut line = String::new();
        reader.read_line(&mut line).unwrap();
        assert_eq!(line, ">seq1\n");
    }
}
