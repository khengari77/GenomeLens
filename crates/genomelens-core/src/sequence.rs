/// Zero-copy view into a decompressed FASTA buffer.
/// Lifetime `'a` is tied to the buffer that holds the decompressed bytes.
pub struct SequenceView<'a> {
    pub id: &'a [u8],
    pub description: &'a [u8],
    pub seq: Vec<u8>,
}

/// Owned summary produced after scanning a FASTA record.
#[derive(Debug, Clone)]
pub struct SequenceSummary {
    pub id: String,
    pub description: String,
    pub length: usize,
    pub gc_count: usize,
    pub n_count: usize,
}

impl SequenceSummary {
    pub fn gc_content(&self) -> f64 {
        if self.length == 0 {
            return 0.0;
        }
        self.gc_count as f64 / self.length as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gc_content_empty() {
        let s = SequenceSummary {
            id: String::new(),
            description: String::new(),
            length: 0,
            gc_count: 0,
            n_count: 0,
        };
        assert_eq!(s.gc_content(), 0.0);
    }

    #[test]
    fn gc_content_all_gc() {
        let s = SequenceSummary {
            id: "test".into(),
            description: String::new(),
            length: 10,
            gc_count: 10,
            n_count: 0,
        };
        assert_eq!(s.gc_content(), 1.0);
    }

    #[test]
    fn gc_content_half() {
        let s = SequenceSummary {
            id: "test".into(),
            description: String::new(),
            length: 100,
            gc_count: 50,
            n_count: 0,
        };
        assert!((s.gc_content() - 0.5).abs() < f64::EPSILON);
    }
}
