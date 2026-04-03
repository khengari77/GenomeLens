use crate::SequenceSummary;

/// Aggregate statistics for a FASTA file.
#[derive(Debug, Clone)]
pub struct FastaStats {
    pub sequence_count: usize,
    pub total_length: usize,
    pub min_length: usize,
    pub max_length: usize,
    pub mean_length: f64,
    pub n50: usize,
    pub gc_content: f64,
    pub n_content: f64,
}

impl FastaStats {
    pub fn from_summaries(summaries: &[SequenceSummary]) -> Self {
        if summaries.is_empty() {
            return Self {
                sequence_count: 0,
                total_length: 0,
                min_length: 0,
                max_length: 0,
                mean_length: 0.0,
                n50: 0,
                gc_content: 0.0,
                n_content: 0.0,
            };
        }

        let sequence_count = summaries.len();
        let total_length: usize = summaries.iter().map(|s| s.length).sum();
        let total_gc: usize = summaries.iter().map(|s| s.gc_count).sum();
        let total_n: usize = summaries.iter().map(|s| s.n_count).sum();
        let min_length = summaries.iter().map(|s| s.length).min().unwrap();
        let max_length = summaries.iter().map(|s| s.length).max().unwrap();
        let mean_length = total_length as f64 / sequence_count as f64;

        let gc_content = if total_length > 0 {
            total_gc as f64 / total_length as f64
        } else {
            0.0
        };
        let n_content = if total_length > 0 {
            total_n as f64 / total_length as f64
        } else {
            0.0
        };

        // N50: sort lengths descending, walk until cumulative >= half total
        let mut lengths: Vec<usize> = summaries.iter().map(|s| s.length).collect();
        lengths.sort_unstable_by(|a, b| b.cmp(a));
        let half = (total_length + 1) / 2; // ceiling division
        let mut cumulative = 0usize;
        let mut n50 = 0;
        for &len in &lengths {
            cumulative += len;
            if cumulative >= half {
                n50 = len;
                break;
            }
        }

        Self {
            sequence_count,
            total_length,
            min_length,
            max_length,
            mean_length,
            n50,
            gc_content,
            n_content,
        }
    }
}

/// Streaming accumulator for FASTA statistics. Collects only lengths (for N50)
/// plus running counters. No String allocations per record.
pub struct FastaStatsAccumulator {
    lengths: Vec<usize>,
    total_gc: usize,
    total_n: usize,
    total_length: usize,
}

impl FastaStatsAccumulator {
    pub fn new() -> Self {
        Self {
            lengths: Vec::new(),
            total_gc: 0,
            total_n: 0,
            total_length: 0,
        }
    }

    /// Add a record's numeric stats. No Strings needed.
    pub fn add(&mut self, length: usize, gc_count: usize, n_count: usize) {
        self.lengths.push(length);
        self.total_gc += gc_count;
        self.total_n += n_count;
        self.total_length += length;
    }

    pub fn finish(mut self) -> FastaStats {
        if self.lengths.is_empty() {
            return FastaStats {
                sequence_count: 0,
                total_length: 0,
                min_length: 0,
                max_length: 0,
                mean_length: 0.0,
                n50: 0,
                gc_content: 0.0,
                n_content: 0.0,
            };
        }

        let sequence_count = self.lengths.len();
        let min_length = *self.lengths.iter().min().unwrap();
        let max_length = *self.lengths.iter().max().unwrap();
        let mean_length = self.total_length as f64 / sequence_count as f64;

        let gc_content = if self.total_length > 0 {
            self.total_gc as f64 / self.total_length as f64
        } else {
            0.0
        };
        let n_content = if self.total_length > 0 {
            self.total_n as f64 / self.total_length as f64
        } else {
            0.0
        };

        // N50
        self.lengths.sort_unstable_by(|a, b| b.cmp(a));
        let half = (self.total_length + 1) / 2;
        let mut cumulative = 0usize;
        let mut n50 = 0;
        for &len in &self.lengths {
            cumulative += len;
            if cumulative >= half {
                n50 = len;
                break;
            }
        }

        FastaStats {
            sequence_count,
            total_length: self.total_length,
            min_length,
            max_length,
            mean_length,
            n50,
            gc_content,
            n_content,
        }
    }
}

impl Default for FastaStatsAccumulator {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_summary(id: &str, length: usize, gc: usize, n: usize) -> SequenceSummary {
        SequenceSummary {
            id: id.to_string(),
            description: String::new(),
            length,
            gc_count: gc,
            n_count: n,
        }
    }

    #[test]
    fn empty_summaries() {
        let stats = FastaStats::from_summaries(&[]);
        assert_eq!(stats.sequence_count, 0);
        assert_eq!(stats.total_length, 0);
        assert_eq!(stats.n50, 0);
    }

    #[test]
    fn single_sequence() {
        let summaries = vec![make_summary("seq1", 100, 40, 5)];
        let stats = FastaStats::from_summaries(&summaries);
        assert_eq!(stats.sequence_count, 1);
        assert_eq!(stats.total_length, 100);
        assert_eq!(stats.min_length, 100);
        assert_eq!(stats.max_length, 100);
        assert_eq!(stats.n50, 100);
        assert!((stats.gc_content - 0.4).abs() < f64::EPSILON);
    }

    #[test]
    fn n50_known_value() {
        // Lengths: 10, 20, 30, 40 → total=100, half=50
        // Sorted desc: 40, 30, 20, 10
        // cumulative: 40 < 50, 40+30=70 >= 50 → N50=30
        let summaries = vec![
            make_summary("a", 10, 0, 0),
            make_summary("b", 20, 0, 0),
            make_summary("c", 30, 0, 0),
            make_summary("d", 40, 0, 0),
        ];
        let stats = FastaStats::from_summaries(&summaries);
        assert_eq!(stats.n50, 30);
    }

    #[test]
    fn n50_equal_lengths() {
        let summaries = vec![
            make_summary("a", 50, 0, 0),
            make_summary("b", 50, 0, 0),
        ];
        let stats = FastaStats::from_summaries(&summaries);
        assert_eq!(stats.n50, 50);
    }

    #[test]
    fn accumulator_empty() {
        let acc = FastaStatsAccumulator::new();
        let stats = acc.finish();
        assert_eq!(stats.sequence_count, 0);
        assert_eq!(stats.total_length, 0);
        assert_eq!(stats.n50, 0);
    }

    #[test]
    fn accumulator_matches_from_summaries() {
        let summaries = vec![
            make_summary("a", 10, 3, 1),
            make_summary("b", 20, 8, 2),
            make_summary("c", 30, 12, 0),
            make_summary("d", 40, 20, 5),
        ];
        let expected = FastaStats::from_summaries(&summaries);

        let mut acc = FastaStatsAccumulator::new();
        for s in &summaries {
            acc.add(s.length, s.gc_count, s.n_count);
        }
        let actual = acc.finish();

        assert_eq!(actual.sequence_count, expected.sequence_count);
        assert_eq!(actual.total_length, expected.total_length);
        assert_eq!(actual.min_length, expected.min_length);
        assert_eq!(actual.max_length, expected.max_length);
        assert_eq!(actual.n50, expected.n50);
        assert!((actual.gc_content - expected.gc_content).abs() < f64::EPSILON);
        assert!((actual.n_content - expected.n_content).abs() < f64::EPSILON);
        assert!((actual.mean_length - expected.mean_length).abs() < f64::EPSILON);
    }

    mod prop {
        use super::*;
        use proptest::prelude::*;

        fn arb_summary() -> impl Strategy<Value = SequenceSummary> {
            (1..10000usize).prop_flat_map(|length| {
                (Just(length), 0..=length).prop_flat_map(move |(length, gc)| {
                    let max_n = length - gc;
                    (Just(length), Just(gc), 0..=max_n)
                })
            })
            .prop_map(|(length, gc, n)| make_summary("seq", length, gc, n))
        }

        proptest! {
            #[test]
            fn total_length_is_sum(summaries in proptest::collection::vec(arb_summary(), 1..20)) {
                let stats = FastaStats::from_summaries(&summaries);
                let expected: usize = summaries.iter().map(|s| s.length).sum();
                prop_assert_eq!(stats.total_length, expected);
            }

            #[test]
            fn n50_is_one_of_the_lengths(summaries in proptest::collection::vec(arb_summary(), 1..20)) {
                let stats = FastaStats::from_summaries(&summaries);
                let lengths: Vec<usize> = summaries.iter().map(|s| s.length).collect();
                prop_assert!(lengths.contains(&stats.n50));
            }

            #[test]
            fn n50_gte_min_length(summaries in proptest::collection::vec(arb_summary(), 1..20)) {
                let stats = FastaStats::from_summaries(&summaries);
                prop_assert!(stats.n50 >= stats.min_length);
            }

            #[test]
            fn gc_content_in_range(summaries in proptest::collection::vec(arb_summary(), 1..20)) {
                let stats = FastaStats::from_summaries(&summaries);
                prop_assert!(stats.gc_content >= 0.0 && stats.gc_content <= 1.0);
            }

            #[test]
            fn accumulator_matches_from_summaries(summaries in proptest::collection::vec(arb_summary(), 1..20)) {
                let expected = FastaStats::from_summaries(&summaries);
                let mut acc = FastaStatsAccumulator::new();
                for s in &summaries {
                    acc.add(s.length, s.gc_count, s.n_count);
                }
                let actual = acc.finish();
                prop_assert_eq!(actual.sequence_count, expected.sequence_count);
                prop_assert_eq!(actual.total_length, expected.total_length);
                prop_assert_eq!(actual.min_length, expected.min_length);
                prop_assert_eq!(actual.max_length, expected.max_length);
                prop_assert_eq!(actual.n50, expected.n50);
                prop_assert!((actual.gc_content - expected.gc_content).abs() < f64::EPSILON);
                prop_assert!((actual.n_content - expected.n_content).abs() < f64::EPSILON);
            }
        }
    }
}
