# GenomeLens

GenomeLens is a high-performance query engine for genomic file formats (FASTA, VCF). It provides streaming statistics, filtering, and columnar projection with a bytecode VM and zero-copy parsing.

## Architecture

The project is organized as a Cargo workspace with four crates:

- `genomelens-core` – Shared types, errors, and aggregate statistics accumulators.
- `genomelens-parse` – Streaming FASTA/VCF readers, lazy column access, transparent gzip/BGZF decompression.
- `genomelens-query` – SQL-like query parser, AST to bytecode compiler, stack VM.
- `genomelens-cli` – CLI frontend using clap.

## Building

Requires Rust 1.70+ (edition 2021).

```bash
cargo build --release
```

The binary is `target/release/genomelens`.

## Usage

```
genomelens <COMMAND>
```

Commands:

- `stats <file>` – Print summary statistics for a FASTA or VCF file.
- `query --query <expr> <file>` – Project selected columns and filter records, output TSV.
- `filter --query <expr> <file>` – Filter VCF records, output original VCF (optionally without header).
- `head -n <N> <file>` – Show first N records.

All commands transparently handle `.gz` compressed files (both gzip and BGZF) by magic byte detection.

### Examples

**FASTA statistics**

```bash
genomelens stats testdata/small.fasta
```

```
┌─────────────┬─────────┐
│ Metric      │ Value   │
├─────────────┼─────────┤
│ Sequences   │ 3       │
│ Total length│ 36      │
│ Min length  │ 12      │
│ Max length  │ 12      │
│ Mean length │ 12.0    │
│ N50         │ 12      │
│ GC content  │ 58.33%  │
│ N content   │ 8.33%   │
└─────────────┴─────────┘
```

**VCF statistics**

```bash
genomelens stats testdata/small.vcf
```

**Query with projection and filter**

```bash
genomelens query --query "SELECT CHROM POS INFO.AF WHERE QUAL > 50 AND INFO.AF > 0.2" testdata/small.vcf
```

Output (TSV):

```
CHROM   POS     INFO.AF
chr1    12345   0.5
chr1    67890   0.3
```

**Filter VCF records** (outputs original VCF lines)

```bash
genomelens filter --query "CHROM = 'chr1' AND QUAL > 80" testdata/small.vcf
```

**Head command**

```bash
genomelens head -n 5 testdata/small.vcf
```

## Query Language

The query syntax is a subset of SQL:

- `SELECT` clause (optional) – lists fixed columns or `INFO.field` names. Use `SELECT *` for all columns.
- `WHERE` clause (optional) – filter expression.
- If `SELECT` is omitted, the query is interpreted as a filter expression and outputs full VCF lines (same as `filter` command).

### Supported columns

| Column  | Type      | Notes                         |
|---------|-----------|-------------------------------|
| CHROM   | string    |                                |
| POS     | integer   |                                |
| QUAL    | float     | Missing (`.`) evaluates to false |
| FILTER  | string    |                                |
| REF     | string    |                                |
| TYPE    | string    | `SNV`, `Insertion`, `Deletion`, `MNV`, `Complex` |
| INFO.`<field>` | varies | Type-checked against VCF header |

### Operators

- Comparison: `=`, `!=`, `>`, `<`, `>=`, `<=`
- Logical: `AND`, `OR`, `NOT`
- Parentheses for grouping
- INFO flags can be used as booleans: `INFO.DB` (true if present)

### Literals

- Integers: `42`
- Floats: `3.14`
- Strings: `'chr1'` (single quotes)

## Implementation Details

### Parsing

- FASTA: Streaming line reader that accumulates sequences without storing full record bodies. `next_stats()` counts GC/N content inline, O(1) memory.
- VCF: `VcfRecord` stores a raw byte slice and caches tab offsets lazily using `memchr`. Only requested columns are parsed. The INFO column is scanned on-demand via linear scan over `;` separators – suitable for moderate numbers of INFO fields.
- Decompression: `open_transparent` detects gzip magic bytes (`0x1f 0x8b`) and wraps the reader in `flate2::read::MultiGzDecoder` (handles both gzip and BGZF).

### Query Compilation

1. Lexer (Logos) → token stream.
2. Pratt parser → AST (`Expr`, `Select`).
3. Planner checks column types against VCF header and emits `OpCode` sequence.
4. VM evaluates using a fixed-size stack (32 bools) and instruction pointer. Short-circuit jumps (`JumpIfFalse`, `JumpIfTrue`) avoid unnecessary INFO parsing.

### Zero-Copy Design

- `VcfRecord` and `FastaRecord` borrow from internal buffers. The streaming readers reuse a single buffer across records; previous records are invalidated on `next_record()`.
- Column accessors return byte slices; numeric conversion happens only when required (e.g., `QUAL > 50` parses QUAL to f64).

### Error Handling

- Custom `Error` enum with `Parse`, `Io`, `InvalidFasta`, `InvalidVcf`, `InvalidQuery`.
- Parsing errors include byte offset where possible.

## Testing

Unit tests and property tests (proptest) are included in each crate.

```bash
cargo test --workspace
```

Test data resides in `testdata/` (small FASTA/VCF files and their gzipped counterparts).

## Performance Characteristics

- FASTA stats: Scans each base once; no per-record heap allocation for sequences.
- VCF filtering: For a query `CHROM = 'chr1' AND QUAL > 50`, only the CHROM and QUAL columns are parsed. If CHROM mismatches, QUAL is never accessed (short-circuit evaluation in VM).
- INFO field access: Linear scan of the INFO string per record per distinct field. Suitable for typical VCFs with <100 INFO fields. For large numbers of fields, a precomputed index could be added.

## Limitations

- VCF INFO fields with `Number=G` (per-genotype) are not handled specially; the VM treats them as raw strings.
- No support for symbolic alleles (`<DEL>`) beyond classifying them as `Complex`.
- Query planner requires a complete VCF header to type-check INFO fields. Headerless VCFs are not supported for `INFO.` references.
- FASTA reader does not validate IUPAC ambiguity codes beyond `N`; all other characters are counted as bases but not categorized.

## Project Status

Proof-of-concept quality. Core features (stats, filter, query projection) are implemented and tested against small files. Not yet optimized for multi-gigabyte whole-genome VCFs (no indexing, no chunked parallelism).

## License

[MIT](LICENSE)
