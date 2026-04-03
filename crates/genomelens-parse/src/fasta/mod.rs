mod streaming;

pub use streaming::{FastaReader, FastaRecordStats};

use genomelens_core::{SequenceSummary, SequenceView};

/// Split a header line into id (up to first space/tab) and description (rest, trimmed).
pub(crate) fn split_header(header: &[u8]) -> (&[u8], &[u8]) {
    for (i, &b) in header.iter().enumerate() {
        if b == b' ' || b == b'\t' {
            let desc = &header[i + 1..];
            return (&header[..i], desc);
        }
    }
    (header, &[])
}

/// A parsed FASTA record holding slices into the original buffer.
pub struct FastaRecord<'a> {
    pub id: &'a [u8],
    pub description: &'a [u8],
    /// Raw sequence bytes (newlines stripped in streaming reader).
    pub(crate) seq_region: &'a [u8],
}

impl<'a> FastaRecord<'a> {
    /// Compute summary statistics without allocating a stripped sequence buffer.
    pub fn summarize(&self) -> SequenceSummary {
        let mut length = 0usize;
        let mut gc = 0usize;
        let mut n = 0usize;

        for &b in self.seq_region {
            match b {
                b'\n' | b'\r' => continue,
                b'G' | b'C' | b'g' | b'c' => {
                    gc += 1;
                    length += 1;
                }
                b'N' | b'n' => {
                    n += 1;
                    length += 1;
                }
                _ => {
                    length += 1;
                }
            }
        }

        SequenceSummary {
            id: String::from_utf8_lossy(self.id).into_owned(),
            description: String::from_utf8_lossy(self.description).into_owned(),
            length,
            gc_count: gc,
            n_count: n,
        }
    }

    /// Build a SequenceView by stripping newlines from the sequence region.
    /// This allocates a Vec<u8> for the contiguous sequence.
    pub fn to_view(&self) -> SequenceView<'a> {
        let seq: Vec<u8> = self
            .seq_region
            .iter()
            .copied()
            .filter(|&b| b != b'\n' && b != b'\r')
            .collect();

        SequenceView {
            id: self.id,
            description: self.description,
            seq,
        }
    }

    /// Raw sequence region (may contain newlines).
    pub fn seq_region(&self) -> &'a [u8] {
        self.seq_region
    }
}
