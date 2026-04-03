use genomelens_core::{Error, Result, SequenceSummary};
use std::io::BufRead;

use super::{split_header, FastaRecord};

/// Lightweight stats for a FASTA record, computed inline without buffering the sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecordStats {
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

    /// Yield the next FASTA record, borrowing from the internal record buffer.
    ///
    /// Returns `Ok(None)` at EOF. Each call invalidates the previous record.
    pub fn next_record(&mut self) -> Result<Option<FastaRecord<'_>>> {
        if self.done {
            return Ok(None);
        }

        self.record_buf.clear();

        // Get the header line — either from peeked or by reading
        let header_line = if let Some(header) = self.peeked_header.take() {
            header
        } else {
            // Read until we find a '>' line
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
                    break line[1..].to_vec();
                }
                // Non-empty, non-header line before first record — skip or error
                return Err(Error::InvalidFasta(format!(
                    "expected '>' but found {:?}",
                    line[0] as char
                )));
            }
        };

        // Store header in record_buf: [header_bytes \n sequence_bytes...]
        // We'll track where the header ends so we can split later.
        let header_len = header_line.len();
        self.record_buf.extend_from_slice(&header_line);
        self.record_buf.push(b'\n'); // separator

        // Read sequence lines until next '>' or EOF
        loop {
            self.line_buf.clear();
            let bytes_read = self
                .reader
                .read_until(b'\n', &mut self.line_buf)
                .map_err(Error::Io)?;

            if bytes_read == 0 {
                self.done = true;
                break;
            }

            let line = trim_newline(&self.line_buf);
            if line.is_empty() {
                continue;
            }

            if line[0] == b'>' {
                // This is the next record's header — stash it
                self.peeked_header = Some(line[1..].to_vec());
                break;
            }

            // Append sequence bytes (no newlines)
            self.record_buf.extend_from_slice(line);
        }

        // Build the FastaRecord from record_buf
        // Layout: [header_bytes] [\n] [sequence_bytes]
        let header = &self.record_buf[..header_len];
        let seq_start = header_len + 1; // skip the \n separator
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
        if self.done {
            return Ok(None);
        }

        // Get the header line — same logic as next_record()
        let _header_line = if let Some(header) = self.peeked_header.take() {
            header
        } else {
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
                    break line[1..].to_vec();
                }
                return Err(Error::InvalidFasta(format!(
                    "expected '>' but found {:?}",
                    line[0] as char
                )));
            }
        };

        // Count bases inline — never store the sequence
        let mut length = 0usize;
        let mut gc_count = 0usize;
        let mut n_count = 0usize;

        loop {
            self.line_buf.clear();
            let bytes_read = self
                .reader
                .read_until(b'\n', &mut self.line_buf)
                .map_err(Error::Io)?;

            if bytes_read == 0 {
                self.done = true;
                break;
            }

            let line = trim_newline(&self.line_buf);
            if line.is_empty() {
                continue;
            }

            if line[0] == b'>' {
                self.peeked_header = Some(line[1..].to_vec());
                break;
            }

            // Count bases inline
            for &b in line {
                match b {
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

        Ok(Some(FastaRecordStats {
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
