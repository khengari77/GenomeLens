use crate::variant::VariantType;
use serde::Serialize;
use std::collections::BTreeMap;

/// A single histogram bin with lower/upper bounds and a count.
#[derive(Debug, Clone, Serialize)]
pub struct HistBin {
    pub lower: f64,
    pub upper: f64,
    pub count: usize,
}

/// Per-chromosome variant statistics.
#[derive(Debug, Clone, Serialize)]
pub struct ChromStats {
    pub variant_count: usize,
    pub snv_count: usize,
    pub indel_count: usize,
    pub ts_count: usize,
    pub tv_count: usize,
    pub pass_count: usize,
}

impl ChromStats {
    fn new() -> Self {
        Self {
            variant_count: 0,
            snv_count: 0,
            indel_count: 0,
            ts_count: 0,
            tv_count: 0,
            pass_count: 0,
        }
    }

    pub fn ts_tv_ratio(&self) -> Option<f64> {
        if self.tv_count == 0 {
            return None;
        }
        Some(self.ts_count as f64 / self.tv_count as f64)
    }
}

/// Rich aggregated statistics for a VCF file, designed for dashboard visualization.
#[derive(Debug, Clone, Serialize)]
pub struct VcfDashboardStats {
    pub variant_count: usize,
    pub sample_count: usize,
    pub snv_count: usize,
    pub insertion_count: usize,
    pub deletion_count: usize,
    pub mnv_count: usize,
    pub complex_count: usize,
    pub multiallelic_count: usize,
    pub pass_count: usize,
    pub ts_count: usize,
    pub tv_count: usize,
    pub ts_tv_ratio: Option<f64>,
    pub pass_rate: f64,

    /// Per-chromosome breakdown.
    pub per_chrom: BTreeMap<String, ChromStats>,

    /// Variant type counts for pie charts.
    pub type_counts: BTreeMap<String, usize>,

    /// QUAL score distribution.
    pub qual_histogram: Vec<HistBin>,

    /// Allele frequency spectrum (None if no AF field in header).
    pub af_histogram: Option<Vec<HistBin>>,

    /// Read depth distribution (None if no DP field in header).
    pub dp_histogram: Option<Vec<HistBin>>,
}

/// Streaming accumulator for dashboard statistics. Single pass over the file.
pub struct VcfDashboardAccumulator {
    sample_count: usize,
    has_af: bool,
    has_dp: bool,

    variant_count: usize,
    snv_count: usize,
    insertion_count: usize,
    deletion_count: usize,
    mnv_count: usize,
    complex_count: usize,
    multiallelic_count: usize,
    pass_count: usize,
    ts_count: usize,
    tv_count: usize,

    per_chrom: BTreeMap<String, ChromStats>,
    last_chrom: Option<String>,

    /// Raw QUAL values for histogram binning at finish time.
    qual_values: Vec<f64>,
    /// Raw AF values (all per-allele values flattened).
    af_values: Vec<f64>,
    /// Raw DP values.
    dp_values: Vec<f64>,
}

impl VcfDashboardAccumulator {
    pub fn new(sample_count: usize, has_af: bool, has_dp: bool) -> Self {
        Self {
            sample_count,
            has_af,
            has_dp,
            variant_count: 0,
            snv_count: 0,
            insertion_count: 0,
            deletion_count: 0,
            mnv_count: 0,
            complex_count: 0,
            multiallelic_count: 0,
            pass_count: 0,
            ts_count: 0,
            tv_count: 0,
            per_chrom: BTreeMap::new(),
            last_chrom: None,
            qual_values: Vec::new(),
            af_values: Vec::new(),
            dp_values: Vec::new(),
        }
    }

    /// Feed a record's extracted fields into the accumulator.
    pub fn add(
        &mut self,
        chrom: &str,
        variant_type: &VariantType,
        is_pass: bool,
        is_multiallelic: bool,
        is_transition: Option<bool>,
        qual: Option<f64>,
        af_raw: Option<&[u8]>,
        dp_raw: Option<&[u8]>,
    ) {
        self.variant_count += 1;

        // Per-chrom stats (with last-chrom cache for sorted files)
        let chrom_stats = if self.last_chrom.as_deref() == Some(chrom) {
            self.per_chrom.get_mut(chrom).unwrap()
        } else {
            self.last_chrom = Some(chrom.to_owned());
            self.per_chrom
                .entry(chrom.to_owned())
                .or_insert_with(ChromStats::new)
        };

        chrom_stats.variant_count += 1;

        match variant_type {
            VariantType::Snv => {
                self.snv_count += 1;
                chrom_stats.snv_count += 1;
            }
            VariantType::Insertion => {
                self.insertion_count += 1;
                chrom_stats.indel_count += 1;
            }
            VariantType::Deletion => {
                self.deletion_count += 1;
                chrom_stats.indel_count += 1;
            }
            VariantType::Mnv => {
                self.mnv_count += 1;
            }
            VariantType::Ref => {}
            VariantType::Complex => {
                self.complex_count += 1;
            }
        }

        if is_pass {
            self.pass_count += 1;
            chrom_stats.pass_count += 1;
        }
        if is_multiallelic {
            self.multiallelic_count += 1;
        }
        match is_transition {
            Some(true) => {
                self.ts_count += 1;
                chrom_stats.ts_count += 1;
            }
            Some(false) => {
                self.tv_count += 1;
                chrom_stats.tv_count += 1;
            }
            None => {}
        }

        // Collect distribution values
        if let Some(q) = qual {
            self.qual_values.push(q);
        }

        if self.has_af {
            if let Some(raw) = af_raw {
                for part in raw.split(|&b| b == b',') {
                    if part != b"." {
                        if let Ok(v) = fast_parse_f64(part) {
                            self.af_values.push(v);
                        }
                    }
                }
            }
        }

        if self.has_dp {
            if let Some(raw) = dp_raw {
                if raw != b"." {
                    if let Ok(v) = fast_parse_f64(raw) {
                        self.dp_values.push(v);
                    }
                }
            }
        }
    }

    pub fn finish(self) -> VcfDashboardStats {
        let ts_tv_ratio = if self.tv_count == 0 {
            None
        } else {
            Some(self.ts_count as f64 / self.tv_count as f64)
        };

        let pass_rate = if self.variant_count == 0 {
            0.0
        } else {
            self.pass_count as f64 / self.variant_count as f64
        };

        let mut type_counts = BTreeMap::new();
        if self.snv_count > 0 {
            type_counts.insert("SNV".to_string(), self.snv_count);
        }
        if self.insertion_count > 0 {
            type_counts.insert("Insertion".to_string(), self.insertion_count);
        }
        if self.deletion_count > 0 {
            type_counts.insert("Deletion".to_string(), self.deletion_count);
        }
        if self.mnv_count > 0 {
            type_counts.insert("MNV".to_string(), self.mnv_count);
        }
        if self.complex_count > 0 {
            type_counts.insert("Complex".to_string(), self.complex_count);
        }

        let qual_histogram = build_histogram(&self.qual_values, 50);
        let af_histogram = if self.has_af && !self.af_values.is_empty() {
            Some(build_fixed_histogram(&self.af_values, 0.0, 1.0, 50))
        } else {
            None
        };
        let dp_histogram = if self.has_dp && !self.dp_values.is_empty() {
            Some(build_histogram(&self.dp_values, 50))
        } else {
            None
        };

        VcfDashboardStats {
            variant_count: self.variant_count,
            sample_count: self.sample_count,
            snv_count: self.snv_count,
            insertion_count: self.insertion_count,
            deletion_count: self.deletion_count,
            mnv_count: self.mnv_count,
            complex_count: self.complex_count,
            multiallelic_count: self.multiallelic_count,
            pass_count: self.pass_count,
            ts_count: self.ts_count,
            tv_count: self.tv_count,
            ts_tv_ratio,
            pass_rate,
            per_chrom: self.per_chrom,
            type_counts,
            qual_histogram,
            af_histogram,
            dp_histogram,
        }
    }
}

/// Build a histogram with automatic range detection.
fn build_histogram(values: &[f64], num_bins: usize) -> Vec<HistBin> {
    if values.is_empty() || num_bins == 0 {
        return Vec::new();
    }

    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for &v in values {
        if v < min {
            min = v;
        }
        if v > max {
            max = v;
        }
    }

    if (max - min).abs() < f64::EPSILON {
        return vec![HistBin {
            lower: min,
            upper: max,
            count: values.len(),
        }];
    }

    build_fixed_histogram(values, min, max, num_bins)
}

/// Build a histogram with a fixed range.
fn build_fixed_histogram(values: &[f64], min: f64, max: f64, num_bins: usize) -> Vec<HistBin> {
    let bin_width = (max - min) / num_bins as f64;
    let mut counts = vec![0usize; num_bins];

    for &v in values {
        let idx = ((v - min) / bin_width) as usize;
        let idx = idx.min(num_bins - 1); // clamp last value into final bin
        counts[idx] += 1;
    }

    counts
        .into_iter()
        .enumerate()
        .map(|(i, count)| HistBin {
            lower: min + i as f64 * bin_width,
            upper: min + (i + 1) as f64 * bin_width,
            count,
        })
        .collect()
}

/// Minimal f64 parser for ASCII byte slices.
fn fast_parse_f64(bytes: &[u8]) -> std::result::Result<f64, ()> {
    std::str::from_utf8(bytes)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .ok_or(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_accumulator() {
        let acc = VcfDashboardAccumulator::new(0, false, false);
        let stats = acc.finish();
        assert_eq!(stats.variant_count, 0);
        assert_eq!(stats.pass_rate, 0.0);
        assert!(stats.ts_tv_ratio.is_none());
        assert!(stats.per_chrom.is_empty());
        assert!(stats.qual_histogram.is_empty());
        assert!(stats.af_histogram.is_none());
        assert!(stats.dp_histogram.is_none());
    }

    #[test]
    fn basic_accumulation() {
        let mut acc = VcfDashboardAccumulator::new(10, true, true);

        acc.add("chr1", &VariantType::Snv, true, false, Some(true), Some(50.0), Some(b"0.3"), Some(b"100"));
        acc.add("chr1", &VariantType::Snv, true, false, Some(false), Some(30.0), Some(b"0.1,0.2"), Some(b"80"));
        acc.add("chr2", &VariantType::Insertion, false, false, None, Some(10.0), None, Some(b"50"));
        acc.add("chr2", &VariantType::Deletion, true, true, None, None, Some(b".,0.5"), Some(b"."));

        let stats = acc.finish();
        assert_eq!(stats.variant_count, 4);
        assert_eq!(stats.sample_count, 10);
        assert_eq!(stats.snv_count, 2);
        assert_eq!(stats.insertion_count, 1);
        assert_eq!(stats.deletion_count, 1);
        assert_eq!(stats.pass_count, 3);
        assert_eq!(stats.multiallelic_count, 1);
        assert_eq!(stats.ts_count, 1);
        assert_eq!(stats.tv_count, 1);
        assert_eq!(stats.ts_tv_ratio.unwrap(), 1.0);
        assert!((stats.pass_rate - 0.75).abs() < f64::EPSILON);

        // Per-chrom
        assert_eq!(stats.per_chrom.len(), 2);
        let chr1 = &stats.per_chrom["chr1"];
        assert_eq!(chr1.variant_count, 2);
        assert_eq!(chr1.snv_count, 2);
        let chr2 = &stats.per_chrom["chr2"];
        assert_eq!(chr2.variant_count, 2);
        assert_eq!(chr2.indel_count, 2);

        // Type counts
        assert_eq!(stats.type_counts["SNV"], 2);
        assert_eq!(stats.type_counts["Insertion"], 1);
        assert_eq!(stats.type_counts["Deletion"], 1);

        // Histograms populated
        assert!(!stats.qual_histogram.is_empty());
        assert!(stats.af_histogram.is_some());
        // AF values: 0.3, 0.1, 0.2, 0.5 (dot skipped)
        let af_total: usize = stats.af_histogram.unwrap().iter().map(|b| b.count).sum();
        assert_eq!(af_total, 4);
        assert!(stats.dp_histogram.is_some());
        // DP values: 100, 80, 50 (dot skipped)
        let dp_total: usize = stats.dp_histogram.unwrap().iter().map(|b| b.count).sum();
        assert_eq!(dp_total, 3);
    }

    #[test]
    fn histogram_single_value() {
        let bins = build_histogram(&[42.0, 42.0, 42.0], 10);
        assert_eq!(bins.len(), 1);
        assert_eq!(bins[0].count, 3);
    }

    #[test]
    fn histogram_empty() {
        let bins = build_histogram(&[], 10);
        assert!(bins.is_empty());
    }

    #[test]
    fn fixed_histogram_all_bins() {
        let values: Vec<f64> = (0..100).map(|i| i as f64 / 100.0).collect();
        let bins = build_fixed_histogram(&values, 0.0, 1.0, 10);
        assert_eq!(bins.len(), 10);
        let total: usize = bins.iter().map(|b| b.count).sum();
        assert_eq!(total, 100);
    }

    #[test]
    fn serializes_to_json() {
        let acc = VcfDashboardAccumulator::new(0, false, false);
        let stats = acc.finish();
        let json = serde_json::to_string(&stats).unwrap();
        assert!(json.contains("\"variant_count\":0"));
    }
}
