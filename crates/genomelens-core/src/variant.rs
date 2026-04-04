/// The type of a value in an INFO or FORMAT field, as declared in VCF headers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VcfValueType {
    Integer,
    Float,
    Flag,
    Character,
    String,
}

/// The "Number" specifier from an INFO/FORMAT header line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VcfNumber {
    /// A fixed count of values.
    Count(usize),
    /// One value per ALT allele ("A").
    PerAltAllele,
    /// One value per allele including REF ("R").
    PerAllele,
    /// One value per genotype ("G").
    PerGenotype,
    /// Unbounded/unknown number of values (".").
    Unbounded,
}

/// Definition of a single INFO field from a `##INFO=<...>` header line.
#[derive(Debug, Clone)]
pub struct InfoFieldDef {
    pub id: String,
    pub number: VcfNumber,
    pub ty: VcfValueType,
    pub description: String,
}

/// Definition of a single FORMAT field from a `##FORMAT=<...>` header line.
#[derive(Debug, Clone)]
pub struct FormatFieldDef {
    pub id: String,
    pub number: VcfNumber,
    pub ty: VcfValueType,
    pub description: String,
}

/// Definition of a FILTER from a `##FILTER=<...>` header line.
#[derive(Debug, Clone)]
pub struct FilterDef {
    pub id: String,
    pub description: String,
}

/// Definition of a contig from a `##contig=<...>` header line.
#[derive(Debug, Clone)]
pub struct ContigDef {
    pub id: String,
    pub length: Option<usize>,
}

/// Classification of a variant by comparing REF and ALT alleles.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VariantType {
    /// Single nucleotide variant (REF and ALT both length 1).
    Snv,
    /// Insertion (ALT longer than REF).
    Insertion,
    /// Deletion (ALT shorter than REF).
    Deletion,
    /// Multi-nucleotide variant (same length, > 1).
    Mnv,
    /// Monomorphic reference site (ALT is `.` or `*`).
    Ref,
    /// Everything else (symbolic, mixed, structural).
    Complex,
}
