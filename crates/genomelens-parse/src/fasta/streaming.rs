use genomelens_core::{Error, Result, SequenceSummary};
use std::io::BufRead;

use super::{split_header, FastaRecord};

/// Lightweight stats for a FASTA record, computed inline without buffering the sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecordStats {
    pub id: Vec<u8>,
    pub length: usize,
    pub gc_count: usize,
    pub n_count: usize,
}

/// Streaming FASTA reader using the lending iterator pattern.
///
/// Handles multi-line sequences by accumulating into an internal buffer.
/// Each call to `next_record()` yields a `FastaRecord<'_>` that borrows
/// from this buffer. The previous record is invalidated on the next call.
pub struct FastaReader<R> {
    reader: R,
    line_buf: Vec<u8>,
    /// Accumulated record data: header line + sequence bytes (newlines stripped).
    record_buf: Vec<u8>,
    /// The header line of the NEXT record, read ahead when we hit a new '>'.
    peeked_header: Option<Vec<u8>>,
    /// Whether we've reached EOF.
    done: bool,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_buf: Vec::with_capacity(4096),
            record_buf: Vec::with_capacity(8192),
            peeked_header: None,
            done: false,
        }
    }

    /// Read the next header line (without the leading '>').
    /// Consumes from `peeked_header` first, otherwise reads lines until '>'.
    fn read_header(&mut self) -> Result<Option<Vec<u8>>> {
        if self.done {
            return Ok(None);
        }
        if let Some(header) = self.peeked_header.take() {
            return Ok(Some(header));
        }
        loop {
            self.line_buf.clear();
            let bytes_read = self
                .reader
                .read_until(b'\n', &mut self.line_buf)
                .map_err(Error::Io)?;
            if bytes_read == 0 {
                self.done = true;
                return Ok(None);
            }
            let line = trim_newline(&self.line_buf);
            if line.is_empty() {
                continue;
            }
            if line[0] == b'>' {
                return Ok(Some(line[1..].to_vec()));
            }
            return Err(Error::InvalidFasta(format!(
                "expected '>' but found {:?}",
                line[0] as char
            )));
        }
    }

    /// Read sequence lines until the next '>' or EOF.
    /// Returns each line to the caller via the `line_buf`. The caller must
    /// process `trim_newline(&self.line_buf)` before the next iteration.
    fn read_next_seq_line(&mut self) -> Result<Option<bool>> {
        self.line_buf.clear();
        let bytes_read = self
            .reader
            .read_until(b'\n', &mut self.line_buf)
            .map_err(Error::Io)?;

        if bytes_read == 0 {
            self.done = true;
            return Ok(None);
        }

        let line = trim_newline(&self.line_buf);
        if line.is_empty() {
            return Ok(Some(false)); // skip empty lines
        }

        if line[0] == b'>' {
            self.peeked_header = Some(line[1..].to_vec());
            return Ok(None);
        }

        Ok(Some(true)) // has sequence data
    }

    /// Yield the next FASTA record, borrowing from the internal record buffer.
    ///
    /// Returns `Ok(None)` at EOF. Each call invalidates the previous record.
    /// **Warning:** Buffers the entire sequence. For genome-scale contigs, prefer
    /// `next_stats()` which uses O(1) memory.
    pub fn next_record(&mut self) -> Result<Option<FastaRecord<'_>>> {
        let header_line = match self.read_header()? {
            Some(h) => h,
            None => return Ok(None),
        };

        self.record_buf.clear();
        let header_len = header_line.len();
        self.record_buf.extend_from_slice(&header_line);
        self.record_buf.push(b'\n');

        loop {
            match self.read_next_seq_line()? {
                None => break,
                Some(false) => continue,
                Some(true) => {
                    let line = trim_newline(&self.line_buf);
                    self.record_buf.extend_from_slice(line);
                }
            }
        }

        let header = &self.record_buf[..header_len];
        let seq_start = header_len + 1;
        let seq_region = if seq_start <= self.record_buf.len() {
            &self.record_buf[seq_start..]
        } else {
            &[]
        };

        let (id, description) = split_header(header);

        Ok(Some(FastaRecord {
            id,
            description,
            seq_region,
        }))
    }

    /// Compute the next record's summary without exposing the full record.
    /// Slightly more efficient path for stats-only workflows.
    pub fn next_summary(&mut self) -> Result<Option<SequenceSummary>> {
        match self.next_record()? {
            Some(rec) => Ok(Some(rec.summarize())),
            None => Ok(None),
        }
    }

    /// Count bases inline without buffering the sequence. O(1) memory regardless
    /// of sequence length. Use this for stats-only workflows on large genomes.
    pub fn next_stats(&mut self) -> Result<Option<FastaRecordStats>> {
        let header_line = match self.read_header()? {
            Some(h) => h,
            None => return Ok(None),
        };

        let (id_bytes, _) = split_header(&header_line);
        let id = id_bytes.to_vec();

        let mut length = 0usize;
        let mut gc_count = 0usize;
        let mut n_count = 0usize;

        loop {
            match self.read_next_seq_line()? {
                None => break,
                Some(false) => continue,
                Some(true) => {
                    let line = trim_newline(&self.line_buf);
                    for &b in line {
                        match b {
                            // Skip whitespace and digits (FASTA formatting characters)
                            b' ' | b'\t' | b'\r' | b'0'..=b'9' => {}
                            b'G' | b'C' | b'g' | b'c' => {
                                gc_count += 1;
                                length += 1;
                            }
                            b'N' | b'n' => {
                                n_count += 1;
                                length += 1;
                            }
                            _ => {
                                length += 1;
                            }
                        }
                    }
                }
            }
        }

        Ok(Some(FastaRecordStats {
            id,
            length,
            gc_count,
            n_count,
        }))
    }
}

fn trim_newline(buf: &[u8]) -> &[u8] {
    let mut end = buf.len();
    if end > 0 && buf[end - 1] == b'\n' {
        end -= 1;
    }
    if end > 0 && buf[end - 1] == b'\r' {
        end -= 1;
    }
    &buf[..end]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn streaming_basic() {
        let data = b">seq1 a description\nACGT\n>seq2\nGGCC\n";
        let mut reader = FastaReader::new(Cursor::new(data.as_slice()));

        let rec1 = reader.next_record().unwrap().unwrap();
        assert_eq!(rec1.id, b"seq1");
        assert_eq!(rec1.description, b"a description");
        assert_eq!(rec1.summarize().length, 4);

        let rec2 = reader.next_record().unwrap().unwrap();
        assert_eq!(rec2.id, b"seq2");
        assert_eq!(rec2.summarize().gc_count, 4);

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn streaming_multiline() {
        let data = b">long\nACGT\nACGT\nACGT\n";
        let mut reader = FastaReader::new(Cursor::new(data.as_slice()));

        let rec = reader.next_record().unwrap().unwrap();
        let summary = rec.summarize();
        assert_eq!(summary.length, 12);
        assert_eq!(summary.gc_count, 6); // 3x CG per line

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn streaming_empty() {
        let data = b"";
        let mut reader = FastaReader::new(Cursor::new(data.as_slice()));
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn streaming_windows_line_endings() {
        let data = b">seq1\r\nACGT\r\n>seq2\r\nGGCC\r\n";
        let mut reader = FastaReader::new(Cursor::new(data.as_slice()));

        let rec1 = reader.next_record().unwrap().unwrap();
        assert_eq!(rec1.summarize().length, 4);

        let rec2 = reader.next_record().unwrap().unwrap();
        assert_eq!(rec2.summarize().length, 4);

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn streaming_whitespace_and_digits_ignored() {
        // FASTA spec: spaces, tabs, and digits are formatting characters, not bases
        let data = b">seq1\nAC GT\n1234\nNN\tCC\n";
        let mut reader = FastaReader::new(Cursor::new(data.as_slice()));

        let rec = reader.next_record().unwrap().unwrap();
        let summary = rec.summarize();
        // Only A, C, G, T, N, N, C, C = 8 bases (spaces, tab, digits skipped)
        assert_eq!(summary.length, 8);
        assert_eq!(summary.gc_count, 4); // C, G, C, C
        assert_eq!(summary.n_count, 2);
    }

    #[test]
    fn next_stats_whitespace_and_digits_ignored() {
        let data = b">seq1\nAC GT\n1234\nNN\tCC\n";
        let mut reader = FastaReader::new(Cursor::new(data.as_slice()));

        let stats = reader.next_stats().unwrap().unwrap();
        assert_eq!(stats.length, 8);
        assert_eq!(stats.gc_count, 4);
        assert_eq!(stats.n_count, 2);
    }

    #[test]
    fn streaming_leading_blank_lines() {
        let data = b"\n\n>seq1\nACGT\n";
        let mut reader = FastaReader::new(Cursor::new(data.as_slice()));

        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.id, b"seq1");
        assert_eq!(rec.summarize().length, 4);

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn streaming_empty_sequence() {
        let data = b">empty\n>seq2\nACGT\n";
        let mut reader = FastaReader::new(Cursor::new(data.as_slice()));

        let rec1 = reader.next_record().unwrap().unwrap();
        assert_eq!(rec1.id, b"empty");
        assert_eq!(rec1.summarize().length, 0);

        let rec2 = reader.next_record().unwrap().unwrap();
        assert_eq!(rec2.id, b"seq2");
        assert_eq!(rec2.summarize().length, 4);

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn streaming_count() {
        let data = b">a\nA\n>b\nC\n>c\nG\n>d\nT\n>e\nN\n";
        let mut reader = FastaReader::new(Cursor::new(data.as_slice()));

        let mut count = 0;
        while reader.next_record().unwrap().is_some() {
            count += 1;
        }
        assert_eq!(count, 5);
    }

    #[test]
    fn next_stats_basic() {
        let data = b">seq1\nACGTNNGC\n>seq2\nAAAA\nCCCC\n";
        let mut reader = FastaReader::new(Cursor::new(data.as_slice()));

        let s1 = reader.next_stats().unwrap().unwrap();
        assert_eq!(s1.length, 8);
        assert_eq!(s1.gc_count, 4); // C, G, G, C
        assert_eq!(s1.n_count, 2);

        let s2 = reader.next_stats().unwrap().unwrap();
        assert_eq!(s2.length, 8);
        assert_eq!(s2.gc_count, 4); // CCCC
        assert_eq!(s2.n_count, 0);

        assert!(reader.next_stats().unwrap().is_none());
    }

    #[test]
    fn next_stats_empty() {
        let data = b"";
        let mut reader = FastaReader::new(Cursor::new(data.as_slice()));
        assert!(reader.next_stats().unwrap().is_none());
    }

    #[test]
    fn next_stats_matches_next_record() {
        let data = b">s1 desc one\nACGTNNGC\n>s2\nAAAA\nCCCC\n>s3\nGGGG\n";

        // next_record + summarize path
        let mut rec_reader = FastaReader::new(Cursor::new(data.as_slice()));
        let mut rec_stats = Vec::new();
        while let Some(rec) = rec_reader.next_record().unwrap() {
            let s = rec.summarize();
            rec_stats.push((s.length, s.gc_count, s.n_count));
        }

        // next_stats path (constant memory)
        let mut stats_reader = FastaReader::new(Cursor::new(data.as_slice()));
        let mut inline_stats = Vec::new();
        while let Some(s) = stats_reader.next_stats().unwrap() {
            inline_stats.push((s.length, s.gc_count, s.n_count));
        }

        assert_eq!(rec_stats, inline_stats);
    }

    #[test]
    fn next_stats_empty_sequence() {
        let data = b">empty\n>seq2\nACGT\n";
        let mut reader = FastaReader::new(Cursor::new(data.as_slice()));

        let s1 = reader.next_stats().unwrap().unwrap();
        assert_eq!(s1.length, 0);

        let s2 = reader.next_stats().unwrap().unwrap();
        assert_eq!(s2.length, 4);

        assert!(reader.next_stats().unwrap().is_none());
    }
}

#[cfg(test)]
mod prop_tests {
    use super::*;
    use proptest::prelude::*;
    use std::io::Cursor;

    /// Generate a random valid sequence ID (alphanumeric + underscore).
    fn arb_id() -> impl Strategy<Value = String> {
        "[A-Za-z][A-Za-z0-9_]{0,20}"
    }

    /// Generate a random DNA sequence of given length.
    fn arb_dna_seq(max_len: usize) -> impl Strategy<Value = Vec<u8>> {
        proptest::collection::vec(
            prop_oneof![
                Just(b'A'),
                Just(b'C'),
                Just(b'G'),
                Just(b'T'),
                Just(b'N'),
            ],
            1..=max_len,
        )
    }

    /// Generate a FASTA record as bytes (with line wrapping).
    fn arb_fasta_record_bytes() -> impl Strategy<Value = (String, Vec<u8>, Vec<u8>)> {
        (arb_id(), arb_dna_seq(200), prop_oneof![Just(60usize), Just(70), Just(80)])
            .prop_map(|(id, seq, wrap)| {
                let mut bytes = Vec::new();
                bytes.push(b'>');
                bytes.extend_from_slice(id.as_bytes());
                bytes.push(b'\n');
                for chunk in seq.chunks(wrap) {
                    bytes.extend_from_slice(chunk);
                    bytes.push(b'\n');
                }
                (id, seq, bytes)
            })
    }

    /// Generate a complete FASTA file with multiple records.
    fn arb_fasta_file() -> impl Strategy<Value = (Vec<(String, Vec<u8>)>, Vec<u8>)> {
        proptest::collection::vec(arb_fasta_record_bytes(), 1..10).prop_map(|records| {
            let mut file_bytes = Vec::new();
            let mut metadata = Vec::new();
            for (id, seq, bytes) in records {
                file_bytes.extend_from_slice(&bytes);
                metadata.push((id, seq));
            }
            (metadata, file_bytes)
        })
    }

    proptest! {
        #[test]
        fn record_count_matches((metadata, file_bytes) in arb_fasta_file()) {
            let mut reader = FastaReader::new(Cursor::new(file_bytes));
            let mut count = 0;
            while reader.next_record().unwrap().is_some() {
                count += 1;
            }
            prop_assert_eq!(count, metadata.len());
        }

        #[test]
        fn ids_match((metadata, file_bytes) in arb_fasta_file()) {
            let mut reader = FastaReader::new(Cursor::new(file_bytes));
            for (expected_id, _) in &metadata {
                let rec = reader.next_record().unwrap().unwrap();
                prop_assert_eq!(
                    std::str::from_utf8(rec.id).unwrap(),
                    expected_id.as_str()
                );
            }
        }

        #[test]
        fn lengths_match((metadata, file_bytes) in arb_fasta_file()) {
            let mut reader = FastaReader::new(Cursor::new(file_bytes));
            for (_, expected_seq) in &metadata {
                let rec = reader.next_record().unwrap().unwrap();
                let summary = rec.summarize();
                prop_assert_eq!(summary.length, expected_seq.len());
            }
        }

        #[test]
        fn gc_count_consistent((metadata, file_bytes) in arb_fasta_file()) {
            let mut reader = FastaReader::new(Cursor::new(file_bytes));
            for (_, expected_seq) in &metadata {
                let rec = reader.next_record().unwrap().unwrap();
                let summary = rec.summarize();
                let expected_gc = expected_seq.iter()
                    .filter(|&&b| matches!(b, b'G' | b'C' | b'g' | b'c'))
                    .count();
                prop_assert_eq!(summary.gc_count, expected_gc);
            }
        }

        #[test]
        fn to_view_matches_seq((_metadata, file_bytes) in arb_fasta_file()) {
            let mut reader = FastaReader::new(Cursor::new(file_bytes));
            while let Some(rec) = reader.next_record().unwrap() {
                let view = rec.to_view();
                let summary = rec.summarize();
                prop_assert_eq!(view.seq.len(), summary.length);
            }
        }
    }
}
