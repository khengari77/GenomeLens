use crate::variant::VariantType;
use crate::Chromosome;
use std::collections::BTreeSet;

/// Lightweight owned summary of a single VCF record, for stats aggregation.
#[derive(Debug, Clone)]
pub struct VcfRecordSummary {
    pub chrom: String,
    pub variant_type: VariantType,
    pub is_pass: bool,
    pub is_multiallelic: bool,
    /// None if not a biallelic SNV; Some(true) for transition, Some(false) for transversion.
    pub is_transition: Option<bool>,
}

/// Aggregate statistics for a VCF file.
#[derive(Debug, Clone)]
pub struct VcfStats {
    pub variant_count: usize,
    pub sample_count: usize,
    pub chromosomes: BTreeSet<Chromosome>,
    pub snv_count: usize,
    pub insertion_count: usize,
    pub deletion_count: usize,
    pub mnv_count: usize,
    pub complex_count: usize,
    pub multiallelic_count: usize,
    pub pass_count: usize,
    pub ts_count: usize,
    pub tv_count: usize,
}

impl VcfStats {
    pub fn from_summaries(summaries: &[VcfRecordSummary], sample_count: usize) -> Self {
        let mut stats = Self {
            variant_count: summaries.len(),
            sample_count,
            chromosomes: BTreeSet::new(),
            snv_count: 0,
            insertion_count: 0,
            deletion_count: 0,
            mnv_count: 0,
            complex_count: 0,
            multiallelic_count: 0,
            pass_count: 0,
            ts_count: 0,
            tv_count: 0,
        };

        for s in summaries {
            stats.chromosomes.insert(Chromosome::from(s.chrom.clone()));

            match s.variant_type {
                VariantType::Snv => stats.snv_count += 1,
                VariantType::Insertion => stats.insertion_count += 1,
                VariantType::Deletion => stats.deletion_count += 1,
                VariantType::Mnv => stats.mnv_count += 1,
                VariantType::Complex => stats.complex_count += 1,
            }

            if s.is_pass {
                stats.pass_count += 1;
            }
            if s.is_multiallelic {
                stats.multiallelic_count += 1;
            }
            match s.is_transition {
                Some(true) => stats.ts_count += 1,
                Some(false) => stats.tv_count += 1,
                None => {}
            }
        }

        stats
    }

    /// Transition/transversion ratio. Returns None if no transversions observed.
    pub fn ts_tv_ratio(&self) -> Option<f64> {
        if self.tv_count == 0 {
            return None;
        }
        Some(self.ts_count as f64 / self.tv_count as f64)
    }
}

/// Streaming accumulator for VCF statistics. Processes records one at a time
/// without collecting into a Vec. Zero allocation per record (except BTreeSet inserts
/// for new chromosomes).
pub struct VcfStatsAccumulator {
    sample_count: usize,
    variant_count: usize,
    chromosomes: BTreeSet<Chromosome>,
    /// Cache for sorted-file fast path — skip BTreeSet lookup when chrom repeats.
    last_chrom: Option<Chromosome>,
    snv_count: usize,
    insertion_count: usize,
    deletion_count: usize,
    mnv_count: usize,
    complex_count: usize,
    multiallelic_count: usize,
    pass_count: usize,
    ts_count: usize,
    tv_count: usize,
}

impl VcfStatsAccumulator {
    pub fn new(sample_count: usize) -> Self {
        Self {
            sample_count,
            variant_count: 0,
            chromosomes: BTreeSet::new(),
            last_chrom: None,
            snv_count: 0,
            insertion_count: 0,
            deletion_count: 0,
            mnv_count: 0,
            complex_count: 0,
            multiallelic_count: 0,
            pass_count: 0,
            ts_count: 0,
            tv_count: 0,
        }
    }

    /// Process a record's fields directly from byte slices.
    /// Only allocates when a genuinely new chromosome is encountered.
    pub fn add(
        &mut self,
        chrom: &[u8],
        variant_type: VariantType,
        is_pass: bool,
        is_multiallelic: bool,
        is_transition: Option<bool>,
    ) {
        self.variant_count += 1;

        // Fast path: sorted VCF files repeat the same chrom for many consecutive records.
        // Check the cache first (zero allocation), then the BTreeSet (zero allocation
        // via Borrow<str>), and only allocate if the chrom is truly new.
        let chrom_str = std::str::from_utf8(chrom).unwrap_or("");
        let is_cached = self
            .last_chrom
            .as_ref()
            .map_or(false, |c| c.as_str() == chrom_str);
        if !is_cached && !self.chromosomes.contains(chrom_str) {
            let c = Chromosome(chrom_str.to_owned());
            self.last_chrom = Some(c.clone());
            self.chromosomes.insert(c);
        } else if !is_cached {
            // Known chrom but not cached — update cache for next iteration
            self.last_chrom = Some(Chromosome(chrom_str.to_owned()));
        }

        match variant_type {
            VariantType::Snv => self.snv_count += 1,
            VariantType::Insertion => self.insertion_count += 1,
            VariantType::Deletion => self.deletion_count += 1,
            VariantType::Mnv => self.mnv_count += 1,
            VariantType::Complex => self.complex_count += 1,
        }

        if is_pass {
            self.pass_count += 1;
        }
        if is_multiallelic {
            self.multiallelic_count += 1;
        }
        match is_transition {
            Some(true) => self.ts_count += 1,
            Some(false) => self.tv_count += 1,
            None => {}
        }
    }

    pub fn finish(self) -> VcfStats {
        VcfStats {
            variant_count: self.variant_count,
            sample_count: self.sample_count,
            chromosomes: self.chromosomes,
            snv_count: self.snv_count,
            insertion_count: self.insertion_count,
            deletion_count: self.deletion_count,
            mnv_count: self.mnv_count,
            complex_count: self.complex_count,
            multiallelic_count: self.multiallelic_count,
            pass_count: self.pass_count,
            ts_count: self.ts_count,
            tv_count: self.tv_count,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_summary(
        chrom: &str,
        vtype: VariantType,
        pass: bool,
        multi: bool,
        ts: Option<bool>,
    ) -> VcfRecordSummary {
        VcfRecordSummary {
            chrom: chrom.to_string(),
            variant_type: vtype,
            is_pass: pass,
            is_multiallelic: multi,
            is_transition: ts,
        }
    }

    #[test]
    fn empty_summaries() {
        let stats = VcfStats::from_summaries(&[], 0);
        assert_eq!(stats.variant_count, 0);
        assert_eq!(stats.chromosomes.len(), 0);
        assert!(stats.ts_tv_ratio().is_none());
    }

    #[test]
    fn accumulator_empty() {
        let acc = VcfStatsAccumulator::new(0);
        let stats = acc.finish();
        assert_eq!(stats.variant_count, 0);
        assert_eq!(stats.chromosomes.len(), 0);
        assert!(stats.ts_tv_ratio().is_none());
    }

    #[test]
    fn accumulator_matches_from_summaries() {
        let summaries = vec![
            make_summary("chr1", VariantType::Snv, true, false, Some(true)),
            make_summary("chr1", VariantType::Snv, true, false, Some(false)),
            make_summary("chr2", VariantType::Insertion, false, false, None),
            make_summary("chr2", VariantType::Deletion, true, true, None),
        ];
        let expected = VcfStats::from_summaries(&summaries, 3);

        let mut acc = VcfStatsAccumulator::new(3);
        for s in &summaries {
            acc.add(
                s.chrom.as_bytes(),
                s.variant_type.clone(),
                s.is_pass,
                s.is_multiallelic,
                s.is_transition,
            );
        }
        let actual = acc.finish();

        assert_eq!(actual.variant_count, expected.variant_count);
        assert_eq!(actual.sample_count, expected.sample_count);
        assert_eq!(actual.chromosomes.len(), expected.chromosomes.len());
        assert_eq!(actual.snv_count, expected.snv_count);
        assert_eq!(actual.insertion_count, expected.insertion_count);
        assert_eq!(actual.deletion_count, expected.deletion_count);
        assert_eq!(actual.pass_count, expected.pass_count);
        assert_eq!(actual.multiallelic_count, expected.multiallelic_count);
        assert_eq!(actual.ts_count, expected.ts_count);
        assert_eq!(actual.tv_count, expected.tv_count);
    }

    #[test]
    fn basic_counts() {
        let summaries = vec![
            make_summary("chr1", VariantType::Snv, true, false, Some(true)),
            make_summary("chr1", VariantType::Snv, true, false, Some(false)),
            make_summary("chr2", VariantType::Insertion, false, false, None),
            make_summary("chr2", VariantType::Deletion, true, true, None),
        ];
        let stats = VcfStats::from_summaries(&summaries, 3);

        assert_eq!(stats.variant_count, 4);
        assert_eq!(stats.sample_count, 3);
        assert_eq!(stats.chromosomes.len(), 2);
        assert_eq!(stats.snv_count, 2);
        assert_eq!(stats.insertion_count, 1);
        assert_eq!(stats.deletion_count, 1);
        assert_eq!(stats.pass_count, 3);
        assert_eq!(stats.multiallelic_count, 1);
        assert_eq!(stats.ts_count, 1);
        assert_eq!(stats.tv_count, 1);
        assert!((stats.ts_tv_ratio().unwrap() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn accumulator_sorted_chrom_cache() {
        // Sorted input: same chrom repeated, then switches once
        let mut acc = VcfStatsAccumulator::new(0);
        for _ in 0..100 {
            acc.add(b"chr1", VariantType::Snv, true, false, None);
        }
        for _ in 0..50 {
            acc.add(b"chr2", VariantType::Snv, true, false, None);
        }
        let stats = acc.finish();
        assert_eq!(stats.variant_count, 150);
        assert_eq!(stats.chromosomes.len(), 2);
    }

    #[test]
    fn accumulator_unsorted_chroms() {
        // Interleaved chroms — exercises the BTreeSet contains() path
        let mut acc = VcfStatsAccumulator::new(0);
        acc.add(b"chr1", VariantType::Snv, true, false, None);
        acc.add(b"chr2", VariantType::Snv, true, false, None);
        acc.add(b"chr1", VariantType::Snv, true, false, None);
        acc.add(b"chr3", VariantType::Snv, true, false, None);
        acc.add(b"chr2", VariantType::Snv, true, false, None);
        let stats = acc.finish();
        assert_eq!(stats.variant_count, 5);
        assert_eq!(stats.chromosomes.len(), 3);
    }
}
