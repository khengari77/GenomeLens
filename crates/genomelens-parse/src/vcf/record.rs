use std::cell::{Cell, RefCell};

use genomelens_core::{Error, Result, VariantType, VcfRecordSummary};
use memchr::memchr;

/// A parsed INFO key-value entry, borrowing from the raw VCF line.
struct InfoEntry<'a> {
    key: &'a [u8],
    /// `None` for boolean flags (key present without `=`).
    value: Option<&'a [u8]>,
}

/// A single VCF data line with lazy, on-demand column parsing.
///
/// Only the raw byte slice is stored at construction time. Tab offsets are
/// located lazily via `memchr` when a column accessor is first called, and
/// cached for subsequent accesses. This means `WHERE CHROM = 'chr1'` only
/// scans for the first tab instead of all 8+.
///
/// The INFO column is NOT parsed upfront — use `info_value()` for lazy lookup.
pub struct VcfRecord<'a> {
    raw: &'a [u8],
    /// Cached absolute byte offsets of tabs [0..8].
    tab_offsets: [Cell<usize>; 9],
    /// How many tabs have been found so far.
    tabs_found: Cell<usize>,
    /// Lazily parsed INFO key-value pairs. `None` = not yet parsed.
    info_cache: RefCell<Option<Vec<InfoEntry<'a>>>>,
}

impl<'a> VcfRecord<'a> {
    // ── Lazy tab scanning ────────────────────────────────────────────

    /// Ensure that tab offsets [0..=need] have been located.
    /// For example, `ensure_tabs(0)` finds the first tab (needed for CHROM).
    fn ensure_tabs(&self, need: usize) -> Result<()> {
        let mut found = self.tabs_found.get();
        if found > need {
            return Ok(());
        }

        let mut start = if found > 0 {
            self.tab_offsets[found - 1].get() + 1
        } else {
            0
        };

        while found <= need {
            match memchr(b'\t', &self.raw[start..]) {
                Some(off) => {
                    self.tab_offsets[found].set(start + off);
                    found += 1;
                    start = start + off + 1;
                }
                None => {
                    self.tabs_found.set(found);
                    return Err(Error::InvalidVcf(format!(
                        "expected at least {} tab-separated columns, found {}",
                        need + 2,
                        found + 1
                    )));
                }
            }
        }
        self.tabs_found.set(found);
        Ok(())
    }

    /// Get the cached offset of tab at `index`.
    #[inline]
    fn tab(&self, index: usize) -> usize {
        self.tab_offsets[index].get()
    }

    /// Try to find tab at `index`. Returns `Ok(Some(offset))` if found,
    /// `Ok(None)` if there are not enough columns (used for optional fields).
    fn try_tab(&self, index: usize) -> Result<Option<usize>> {
        let found = self.tabs_found.get();
        if found > index {
            return Ok(Some(self.tab_offsets[index].get()));
        }

        let mut f = found;
        let mut start = if f > 0 {
            self.tab_offsets[f - 1].get() + 1
        } else {
            0
        };

        while f <= index {
            match memchr(b'\t', &self.raw[start..]) {
                Some(off) => {
                    self.tab_offsets[f].set(start + off);
                    f += 1;
                    start = start + off + 1;
                }
                None => {
                    self.tabs_found.set(f);
                    return Ok(None);
                }
            }
        }
        self.tabs_found.set(f);
        Ok(Some(self.tab_offsets[index].get()))
    }

    // ── Column accessors ─────────────────────────────────────────────

    /// The full raw line (no tab scanning needed).
    #[inline]
    pub fn raw_line(&self) -> &'a [u8] {
        self.raw
    }

    /// CHROM column (column 0). Requires tab 0.
    pub fn chrom(&self) -> Result<&'a [u8]> {
        self.ensure_tabs(0)?;
        Ok(&self.raw[..self.tab(0)])
    }

    /// POS column (column 1). Requires tabs 0-1.
    pub fn pos(&self) -> Result<&'a [u8]> {
        self.ensure_tabs(1)?;
        Ok(&self.raw[self.tab(0) + 1..self.tab(1)])
    }

    /// ID column (column 2). Requires tabs 0-2.
    pub fn id(&self) -> Result<&'a [u8]> {
        self.ensure_tabs(2)?;
        Ok(&self.raw[self.tab(1) + 1..self.tab(2)])
    }

    /// REF column (column 3). Requires tabs 0-3.
    pub fn ref_allele(&self) -> Result<&'a [u8]> {
        self.ensure_tabs(3)?;
        Ok(&self.raw[self.tab(2) + 1..self.tab(3)])
    }

    /// ALT column (column 4). Requires tabs 0-4.
    pub fn alt_alleles(&self) -> Result<&'a [u8]> {
        self.ensure_tabs(4)?;
        Ok(&self.raw[self.tab(3) + 1..self.tab(4)])
    }

    /// QUAL column (column 5). Requires tabs 0-5.
    pub fn qual(&self) -> Result<&'a [u8]> {
        self.ensure_tabs(5)?;
        Ok(&self.raw[self.tab(4) + 1..self.tab(5)])
    }

    /// FILTER column (column 6). Requires tabs 0-6.
    pub fn filter(&self) -> Result<&'a [u8]> {
        self.ensure_tabs(6)?;
        Ok(&self.raw[self.tab(5) + 1..self.tab(6)])
    }

    /// INFO column (column 7). Requires tabs 0-6, tries tab 7 for end boundary.
    pub fn info(&self) -> Result<&'a [u8]> {
        self.ensure_tabs(6)?;
        let start = self.tab(6) + 1;
        let end = match self.try_tab(7)? {
            Some(t7) => t7,
            None => self.raw.len(),
        };
        Ok(&self.raw[start..end])
    }

    /// FORMAT column (column 8, optional). Requires tabs 0-7, tries tab 8.
    pub fn format(&self) -> Result<Option<&'a [u8]>> {
        match self.try_tab(7)? {
            Some(t7) => {
                let start = t7 + 1;
                let end = match self.try_tab(8)? {
                    Some(t8) => t8,
                    None => self.raw.len(),
                };
                Ok(Some(&self.raw[start..end]))
            }
            None => Ok(None),
        }
    }

    /// Samples (everything after tab 8, optional).
    pub fn samples(&self) -> Result<Option<&'a [u8]>> {
        match self.try_tab(8)? {
            Some(t8) => Ok(Some(&self.raw[t8 + 1..])),
            None => Ok(None),
        }
    }

    // ── Derived accessors ────────────────────────────────────────────

    /// Parse POS to usize.
    pub fn pos_usize(&self) -> Result<usize> {
        crate::ascii::parse_ascii_usize(self.pos()?)
    }

    /// Parse QUAL to f64. Returns None if ".".
    pub fn qual_f64(&self) -> Result<Option<f64>> {
        let q = self.qual()?;
        if q == b"." {
            return Ok(None);
        }
        crate::ascii::parse_ascii_f64(q).map(Some)
    }

    /// Iterate over ALT alleles (split on comma).
    pub fn alt_iter(&self) -> Result<impl Iterator<Item = &'a [u8]>> {
        Ok(self.alt_alleles()?.split(|&b| b == b','))
    }

    /// Iterate over FILTER values (split on semicolon).
    pub fn filter_iter(&self) -> Result<impl Iterator<Item = &'a [u8]>> {
        Ok(self.filter()?.split(|&b| b == b';'))
    }

    /// Lazy INFO lookup with per-record caching.
    ///
    /// On the first call, parses the INFO string once and caches all key-value
    /// pairs. Subsequent lookups scan the cache (no re-parsing).
    ///
    /// Returns:
    /// - `None` if the key is not found
    /// - `Some(None)` for boolean flags (key present without `=`)
    /// - `Some(Some(value))` for key=value pairs
    pub fn info_value(&self, key: &[u8]) -> Result<Option<Option<&'a [u8]>>> {
        let info = self.info()?;
        if info == b"." {
            return Ok(None);
        }

        // Populate cache on first access
        {
            let mut cache = self.info_cache.borrow_mut();
            if cache.is_none() {
                let mut entries = Vec::new();
                let mut pos = 0;
                loop {
                    if pos >= info.len() {
                        break;
                    }
                    let field_end = match memchr(b';', &info[pos..]) {
                        Some(offset) => pos + offset,
                        None => info.len(),
                    };
                    let field = &info[pos..field_end];
                    if let Some(eq) = memchr(b'=', field) {
                        entries.push(InfoEntry {
                            key: &field[..eq],
                            value: Some(&field[eq + 1..]),
                        });
                    } else {
                        entries.push(InfoEntry {
                            key: field,
                            value: None,
                        });
                    }
                    pos = field_end + 1;
                }
                *cache = Some(entries);
            }
        }

        // Lookup in cache
        let cache = self.info_cache.borrow();
        let entries = cache.as_ref().unwrap();
        for entry in entries {
            if entry.key == key {
                return Ok(Some(entry.value));
            }
        }
        Ok(None)
    }

    /// Classify variant type by comparing REF and first ALT allele.
    pub fn variant_type(&self) -> Result<VariantType> {
        let first_alt = match self.alt_iter()?.next() {
            Some(alt) if alt != b"." && alt != b"*" => alt,
            _ => return Ok(VariantType::Complex),
        };

        if first_alt.first() == Some(&b'<') {
            return Ok(VariantType::Complex);
        }

        let ref_len = self.ref_allele()?.len();
        let alt_len = first_alt.len();

        Ok(match (ref_len, alt_len) {
            (1, 1) => VariantType::Snv,
            (r, a) if r == a => VariantType::Mnv,
            (r, a) if r < a => VariantType::Insertion,
            (r, a) if r > a => VariantType::Deletion,
            _ => VariantType::Complex,
        })
    }

    /// Check if FILTER is PASS or missing (".").
    pub fn is_pass(&self) -> Result<bool> {
        let f = self.filter()?;
        Ok(f == b"PASS" || f == b".")
    }

    /// Iterate over sample columns (split on tab).
    pub fn sample_iter(&self) -> Result<impl Iterator<Item = &'a [u8]>> {
        Ok(self
            .samples()?
            .into_iter()
            .flat_map(|s| s.split(|&b| b == b'\t')))
    }

    /// Check if the record has more than one ALT allele.
    pub fn is_multiallelic(&self) -> Result<bool> {
        Ok(self.alt_iter()?.nth(1).is_some())
    }

    /// For biallelic SNVs, returns Some(true) for transition, Some(false) for transversion.
    /// Returns None for non-SNV or multiallelic variants.
    pub fn is_transition(&self) -> Result<Option<bool>> {
        if self.variant_type()? == VariantType::Snv && !self.is_multiallelic()? {
            let ref_base = self.ref_allele()?.first().copied();
            let alt_base = self.alt_iter()?.next().and_then(|a| a.first().copied());
            match (ref_base, alt_base) {
                (Some(r), Some(a)) => Ok(Some(is_ts(r, a))),
                _ => Ok(None),
            }
        } else {
            Ok(None)
        }
    }

    /// Build an owned VcfRecordSummary for stats aggregation.
    pub fn summarize(&self) -> Result<VcfRecordSummary> {
        Ok(VcfRecordSummary {
            chrom: String::from_utf8_lossy(self.chrom()?).into_owned(),
            variant_type: self.variant_type()?,
            is_pass: self.is_pass()?,
            is_multiallelic: self.is_multiallelic()?,
            is_transition: self.is_transition()?,
        })
    }
}

/// Check if a base substitution is a transition (A<>G, C<>T).
fn is_ts(ref_base: u8, alt_base: u8) -> bool {
    matches!(
        (ref_base, alt_base),
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
    )
}

/// Construct a lazy VcfRecord from a raw data line.
/// No parsing is performed — tab scanning is deferred to field access.
pub fn parse_line(line: &[u8]) -> VcfRecord<'_> {
    VcfRecord {
        raw: line,
        tab_offsets: core::array::from_fn(|_| Cell::new(0)),
        tabs_found: Cell::new(0),
        info_cache: RefCell::new(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_basic_line() {
        let line = b"chr1\t12345\trs123\tA\tG\t50\tPASS\tDP=14;AF=0.5";
        let rec = parse_line(line);
        assert_eq!(rec.chrom().unwrap(), b"chr1");
        assert_eq!(rec.pos().unwrap(), b"12345");
        assert_eq!(rec.id().unwrap(), b"rs123");
        assert_eq!(rec.ref_allele().unwrap(), b"A");
        assert_eq!(rec.alt_alleles().unwrap(), b"G");
        assert_eq!(rec.qual().unwrap(), b"50");
        assert_eq!(rec.filter().unwrap(), b"PASS");
        assert_eq!(rec.info().unwrap(), b"DP=14;AF=0.5");
        assert!(rec.format().unwrap().is_none());
        assert!(rec.samples().unwrap().is_none());
    }

    #[test]
    fn parse_line_with_samples() {
        let line = b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=10\tGT:DP\t0/1:5\t1/1:8";
        let rec = parse_line(line);
        assert_eq!(rec.format().unwrap(), Some(b"GT:DP".as_slice()));
        assert_eq!(rec.samples().unwrap(), Some(b"0/1:5\t1/1:8".as_slice()));

        let samples: Vec<_> = rec.sample_iter().unwrap().collect();
        assert_eq!(samples.len(), 2);
        assert_eq!(samples[0], b"0/1:5");
        assert_eq!(samples[1], b"1/1:8");
    }

    #[test]
    fn too_few_columns() {
        let line = b"chr1\t100\tA";
        let rec = parse_line(line);
        // Accessing a column that requires more tabs than available should error
        assert!(rec.ref_allele().is_err());
    }

    #[test]
    fn lazy_parsing_only_scans_needed_tabs() {
        let line = b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=14";
        let rec = parse_line(line);
        // Before any access, no tabs have been scanned
        assert_eq!(rec.tabs_found.get(), 0);

        // Accessing CHROM scans only 1 tab
        assert_eq!(rec.chrom().unwrap(), b"chr1");
        assert_eq!(rec.tabs_found.get(), 1);

        // Accessing POS scans up to tab 1 (already have tab 0)
        assert_eq!(rec.pos().unwrap(), b"100");
        assert_eq!(rec.tabs_found.get(), 2);
    }

    #[test]
    fn info_value_key_value() {
        let line = b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=14;AF=0.5;MQ=60";
        let rec = parse_line(line);

        assert_eq!(rec.info_value(b"DP").unwrap(), Some(Some(b"14".as_slice())));
        assert_eq!(
            rec.info_value(b"AF").unwrap(),
            Some(Some(b"0.5".as_slice()))
        );
        assert_eq!(
            rec.info_value(b"MQ").unwrap(),
            Some(Some(b"60".as_slice()))
        );
    }

    #[test]
    fn info_value_flag() {
        let line = b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=14;DB;MQ=60";
        let rec = parse_line(line);
        assert_eq!(rec.info_value(b"DB").unwrap(), Some(None));
    }

    #[test]
    fn info_value_missing() {
        let line = b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=14";
        let rec = parse_line(line);
        assert_eq!(rec.info_value(b"NONEXISTENT").unwrap(), None);
    }

    #[test]
    fn info_value_dot() {
        let line = b"chr1\t100\t.\tA\tG\t50\tPASS\t.";
        let rec = parse_line(line);
        assert_eq!(rec.info_value(b"DP").unwrap(), None);
    }

    #[test]
    fn pos_usize_parse() {
        let line = b"chr1\t12345\t.\tA\tG\t50\tPASS\t.";
        let rec = parse_line(line);
        assert_eq!(rec.pos_usize().unwrap(), 12345);
    }

    #[test]
    fn qual_f64_parse() {
        let line = b"chr1\t100\t.\tA\tG\t99.9\tPASS\t.";
        let rec = parse_line(line);
        assert!((rec.qual_f64().unwrap().unwrap() - 99.9).abs() < f64::EPSILON);
    }

    #[test]
    fn qual_f64_missing() {
        let line = b"chr1\t100\t.\tA\tG\t.\tPASS\t.";
        let rec = parse_line(line);
        assert!(rec.qual_f64().unwrap().is_none());
    }

    #[test]
    fn alt_iter_multiple() {
        let line = b"chr1\t100\t.\tA\tG,T,C\t50\tPASS\t.";
        let rec = parse_line(line);
        let alts: Vec<_> = rec.alt_iter().unwrap().collect();
        assert_eq!(alts, vec![b"G".as_slice(), b"T", b"C"]);
    }

    #[test]
    fn variant_type_snv() {
        let line = b"chr1\t100\t.\tA\tG\t50\tPASS\t.";
        let rec = parse_line(line);
        assert_eq!(rec.variant_type().unwrap(), VariantType::Snv);
    }

    #[test]
    fn variant_type_insertion() {
        let line = b"chr1\t100\t.\tA\tATG\t50\tPASS\t.";
        let rec = parse_line(line);
        assert_eq!(rec.variant_type().unwrap(), VariantType::Insertion);
    }

    #[test]
    fn variant_type_deletion() {
        let line = b"chr1\t100\t.\tATG\tA\t50\tPASS\t.";
        let rec = parse_line(line);
        assert_eq!(rec.variant_type().unwrap(), VariantType::Deletion);
    }

    #[test]
    fn variant_type_mnv() {
        let line = b"chr1\t100\t.\tAT\tGC\t50\tPASS\t.";
        let rec = parse_line(line);
        assert_eq!(rec.variant_type().unwrap(), VariantType::Mnv);
    }

    #[test]
    fn is_pass_true() {
        let line = b"chr1\t100\t.\tA\tG\t50\tPASS\t.";
        let rec = parse_line(line);
        assert!(rec.is_pass().unwrap());
    }

    #[test]
    fn is_pass_dot() {
        let line = b"chr1\t100\t.\tA\tG\t50\t.\t.";
        let rec = parse_line(line);
        assert!(rec.is_pass().unwrap());
    }

    #[test]
    fn is_pass_false() {
        let line = b"chr1\t100\t.\tA\tG\t50\tq10\t.";
        let rec = parse_line(line);
        assert!(!rec.is_pass().unwrap());
    }

    #[test]
    fn summarize_snv_transition() {
        let line = b"chr1\t100\t.\tA\tG\t50\tPASS\t.";
        let rec = parse_line(line);
        let summary = rec.summarize().unwrap();
        assert_eq!(summary.variant_type, VariantType::Snv);
        assert!(summary.is_pass);
        assert!(!summary.is_multiallelic);
        assert_eq!(summary.is_transition, Some(true));
    }

    #[test]
    fn summarize_snv_transversion() {
        let line = b"chr1\t100\t.\tA\tC\t50\tPASS\t.";
        let rec = parse_line(line);
        let summary = rec.summarize().unwrap();
        assert_eq!(summary.is_transition, Some(false));
    }

    #[test]
    fn info_value_prefix_collision() {
        let line = b"chr1\t100\t.\tA\tG\t50\tPASS\tDP=14;D";
        let rec = parse_line(line);
        assert_eq!(
            rec.info_value(b"DP").unwrap(),
            Some(Some(b"14".as_slice()))
        );
        assert_eq!(rec.info_value(b"D").unwrap(), Some(None));
    }
}

#[cfg(test)]
mod prop_tests {
    use super::*;
    use genomelens_core::VariantType;
    use proptest::prelude::*;

    fn arb_chrom() -> impl Strategy<Value = String> {
        prop_oneof![
            (1..=22u32).prop_map(|n| format!("chr{}", n)),
            Just("chrX".to_string()),
            Just("chrY".to_string()),
            Just("chrM".to_string()),
        ]
    }

    fn arb_base() -> impl Strategy<Value = u8> {
        prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')]
    }

    fn arb_allele(min: usize, max: usize) -> impl Strategy<Value = Vec<u8>> {
        proptest::collection::vec(arb_base(), min..=max)
    }

    fn arb_info() -> impl Strategy<Value = String> {
        prop_oneof![
            (1..1000u32).prop_map(|v| format!("DP={}", v)),
            Just(".".to_string()),
        ]
    }

    fn arb_qual() -> impl Strategy<Value = String> {
        prop_oneof![
            (1..10000u32).prop_map(|v| format!("{}", v)),
            Just(".".to_string()),
        ]
    }

    fn arb_filter() -> impl Strategy<Value = String> {
        prop_oneof![
            Just("PASS".to_string()),
            Just(".".to_string()),
            Just("q10".to_string()),
        ]
    }

    fn arb_vcf_line() -> impl Strategy<Value = String> {
        (
            arb_chrom(),
            1..=1000000u32,
            arb_allele(1, 5),
            arb_allele(1, 5),
            arb_qual(),
            arb_filter(),
            arb_info(),
        )
            .prop_map(|(chrom, pos, ref_allele, alt_allele, qual, filter, info)| {
                let ref_str = String::from_utf8(ref_allele).unwrap();
                let alt_str = String::from_utf8(alt_allele).unwrap();
                format!(
                    "{}\t{}\t.\t{}\t{}\t{}\t{}\t{}",
                    chrom, pos, ref_str, alt_str, qual, filter, info
                )
            })
    }

    fn arb_vcf_lines() -> impl Strategy<Value = Vec<String>> {
        proptest::collection::vec(arb_vcf_line(), 1..20)
    }

    proptest! {
        #[test]
        fn record_count_matches_lines(lines in arb_vcf_lines()) {
            let mut parsed = 0;
            for line in &lines {
                let rec = parse_line(line.as_bytes());
                // Access chrom to trigger validation of at least first column
                rec.chrom().unwrap();
                parsed += 1;
            }
            prop_assert_eq!(parsed, lines.len());
        }

        #[test]
        fn all_records_have_valid_chrom(line in arb_vcf_line()) {
            let rec = parse_line(line.as_bytes());
            let chrom = rec.chrom().unwrap();
            prop_assert!(!chrom.is_empty());
            prop_assert!(chrom.starts_with(b"chr"));
        }

        #[test]
        fn all_records_have_parseable_pos(line in arb_vcf_line()) {
            let rec = parse_line(line.as_bytes());
            let pos = rec.pos_usize().unwrap();
            prop_assert!(pos >= 1);
        }

        #[test]
        fn variant_type_consistent_with_allele_lengths(line in arb_vcf_line()) {
            let rec = parse_line(line.as_bytes());
            let vtype = rec.variant_type().unwrap();
            let ref_len = rec.ref_allele().unwrap().len();
            let alt_len = rec.alt_iter().unwrap().next().unwrap().len();

            match vtype {
                VariantType::Snv => {
                    prop_assert_eq!(ref_len, 1);
                    prop_assert_eq!(alt_len, 1);
                }
                VariantType::Mnv => {
                    prop_assert_eq!(ref_len, alt_len);
                    prop_assert!(ref_len > 1);
                }
                VariantType::Insertion => {
                    prop_assert!(alt_len > ref_len);
                }
                VariantType::Deletion => {
                    prop_assert!(ref_len > alt_len);
                }
                VariantType::Complex => {}
            }
        }

        #[test]
        fn is_pass_consistent_with_filter(line in arb_vcf_line()) {
            let rec = parse_line(line.as_bytes());
            if rec.is_pass().unwrap() {
                let f = rec.filter().unwrap();
                prop_assert!(f == b"PASS" || f == b".");
            }
        }

        #[test]
        fn info_dot_returns_none(line in arb_vcf_line()) {
            let rec = parse_line(line.as_bytes());
            if rec.info().unwrap() == b"." {
                prop_assert_eq!(rec.info_value(b"DP").unwrap(), None);
                prop_assert_eq!(rec.info_value(b"AF").unwrap(), None);
            }
        }
    }
}
