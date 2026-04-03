use clap::{Parser, Subcommand};
use comfy_table::presets::UTF8_FULL;
use comfy_table::Table;
use genomelens_core::{FastaStatsAccumulator, VariantType, VcfStatsAccumulator};
use genomelens_parse::{
    detect_format, open_transparent, FastaReader, FileFormat, VcfReader,
};
use std::io::{BufRead, Write};
use std::path::PathBuf;
use std::process;

#[derive(Parser)]
#[command(name = "genomelens", version, about = "High-performance query engine for genomic files")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Print summary statistics for a FASTA or VCF file
    Stats {
        /// Path to FASTA or VCF file (supports .gz compressed)
        file: PathBuf,
    },
    /// Query VCF records with SELECT projection and WHERE filter (TSV output)
    Query {
        /// Path to VCF file (supports .gz compressed)
        file: PathBuf,
        /// Query expression (e.g. "SELECT CHROM POS QUAL WHERE QUAL > 50")
        #[arg(short, long)]
        query: String,
    },
    /// Filter VCF records matching a query expression
    Filter {
        /// Path to VCF file (supports .gz compressed)
        file: PathBuf,
        /// Query expression (e.g. "QUAL > 50 AND CHROM = 'chr1'")
        #[arg(short, long)]
        query: String,
        /// Suppress VCF header in output
        #[arg(long)]
        no_header: bool,
    },
    /// Print the first N records
    Head {
        /// Path to FASTA or VCF file (supports .gz compressed)
        file: PathBuf,
        /// Number of records to show
        #[arg(short, long, default_value_t = 10)]
        n: usize,
    },
}

fn main() {
    let cli = Cli::parse();

    if let Err(e) = run(cli) {
        eprintln!("error: {e}");
        process::exit(1);
    }
}

fn run(cli: Cli) -> genomelens_core::Result<()> {
    match cli.command {
        Command::Stats { file } => {
            let mut reader = open_transparent(&file)?;
            match detect_format(&mut reader)? {
                FileFormat::Fasta => cmd_stats_fasta(reader),
                FileFormat::Vcf => cmd_stats_vcf(reader),
            }
        }
        Command::Query { file, query } => {
            let reader = open_transparent(&file)?;
            let q = genomelens_query::parse_query(&query)?;
            execute_query(reader, q, true)
        }
        Command::Filter {
            file,
            query,
            no_header,
        } => {
            let reader = open_transparent(&file)?;
            let expr = genomelens_query::parse(&query)?;
            let q = genomelens_query::ast::Query {
                select: genomelens_query::ast::Select::All,
                filter: Some(expr),
            };
            execute_query(reader, q, !no_header)
        }
        Command::Head { file, n } => {
            let mut reader = open_transparent(&file)?;
            match detect_format(&mut reader)? {
                FileFormat::Fasta => cmd_head_fasta(reader, n),
                FileFormat::Vcf => cmd_head_vcf(reader, n),
            }
        }
    }
}

fn cmd_stats_fasta(reader: impl BufRead) -> genomelens_core::Result<()> {
    let mut fasta = FastaReader::new(reader);
    let mut acc = FastaStatsAccumulator::new();

    while let Some(stats) = fasta.next_stats()? {
        acc.add(stats.length, stats.gc_count, stats.n_count);
    }

    let stats = acc.finish();

    let mut table = Table::new();
    table.load_preset(UTF8_FULL);
    table.set_header(vec!["Metric", "Value"]);
    table.add_row(vec!["Sequences", &stats.sequence_count.to_string()]);
    table.add_row(vec!["Total length", &stats.total_length.to_string()]);
    table.add_row(vec!["Min length", &stats.min_length.to_string()]);
    table.add_row(vec!["Max length", &stats.max_length.to_string()]);
    table.add_row(vec!["Mean length", &format!("{:.1}", stats.mean_length)]);
    table.add_row(vec!["N50", &stats.n50.to_string()]);
    table.add_row(vec![
        "GC content",
        &format!("{:.2}%", stats.gc_content * 100.0),
    ]);
    table.add_row(vec![
        "N content",
        &format!("{:.2}%", stats.n_content * 100.0),
    ]);

    println!("{table}");
    Ok(())
}

fn cmd_stats_vcf(reader: impl BufRead) -> genomelens_core::Result<()> {
    let mut vcf = VcfReader::new(reader);
    let header = vcf.read_header()?;
    let mut acc = VcfStatsAccumulator::new(header.sample_names.len());

    while let Some(rec) = vcf.next_record()? {
        acc.add(
            rec.chrom()?,
            rec.variant_type()?,
            rec.is_pass()?,
            rec.is_multiallelic()?,
            rec.is_transition()?,
        );
    }

    let stats = acc.finish();
    let chroms: Vec<_> = stats.chromosomes.iter().map(|c| c.as_str()).collect();

    let mut table = Table::new();
    table.load_preset(UTF8_FULL);
    table.set_header(vec!["Metric", "Value"]);
    table.add_row(vec!["Variants", &stats.variant_count.to_string()]);
    table.add_row(vec!["Samples", &stats.sample_count.to_string()]);
    table.add_row(vec!["Chromosomes", &chroms.join(", ")]);
    table.add_row(vec!["SNVs", &stats.snv_count.to_string()]);
    table.add_row(vec!["Insertions", &stats.insertion_count.to_string()]);
    table.add_row(vec!["Deletions", &stats.deletion_count.to_string()]);
    table.add_row(vec!["MNVs", &stats.mnv_count.to_string()]);
    table.add_row(vec!["Complex", &stats.complex_count.to_string()]);
    table.add_row(vec!["Multiallelic", &stats.multiallelic_count.to_string()]);
    table.add_row(vec!["PASS", &stats.pass_count.to_string()]);
    table.add_row(vec![
        "Ts/Tv ratio",
        &match stats.ts_tv_ratio() {
            Some(r) => format!("{r:.2}"),
            None => "N/A".to_string(),
        },
    ]);

    println!("{table}");
    Ok(())
}

fn cmd_head_fasta(reader: impl BufRead, n: usize) -> genomelens_core::Result<()> {
    let mut fasta = FastaReader::new(reader);

    let mut table = Table::new();
    table.load_preset(UTF8_FULL);
    table.set_header(vec!["#", "ID", "Length", "GC%"]);

    let mut count = 0;
    while let Some(rec) = fasta.next_record()? {
        if count >= n {
            break;
        }
        let summary = rec.summarize();
        table.add_row(vec![
            (count + 1).to_string(),
            summary.id.clone(),
            summary.length.to_string(),
            format!("{:.1}%", summary.gc_content() * 100.0),
        ]);
        count += 1;
    }

    println!("{table}");
    Ok(())
}

fn execute_query(
    reader: impl BufRead,
    query: genomelens_query::ast::Query,
    print_header: bool,
) -> genomelens_core::Result<()> {
    use genomelens_query::ast::Select;

    let mut vcf = VcfReader::new(reader);
    let header = vcf.read_header()?;

    let mut vm = match &query.filter {
        Some(expr) => {
            let ops = genomelens_query::plan(expr, &header)?;
            Some(genomelens_query::Vm::new(ops))
        }
        None => None,
    };

    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());

    match &query.select {
        Select::All => {
            if print_header {
                for line in &header.raw_meta_lines {
                    out.write_all(line.as_bytes())?;
                    out.write_all(b"\n")?;
                }
                out.write_all(header.chrom_line.as_bytes())?;
                out.write_all(b"\n")?;
            }

            while let Some(rec) = vcf.next_record()? {
                let pass = match &mut vm {
                    Some(vm) => vm.evaluate(&rec)?,
                    None => true,
                };
                if pass {
                    out.write_all(rec.raw_line())?;
                    out.write_all(b"\n")?;
                }
            }
        }
        Select::Columns(columns) => {
            genomelens_query::write_column_header(columns, &mut out)?;
            out.write_all(b"\n")?;

            while let Some(rec) = vcf.next_record()? {
                let pass = match &mut vm {
                    Some(vm) => vm.evaluate(&rec)?,
                    None => true,
                };
                if pass {
                    for (i, col) in columns.iter().enumerate() {
                        if i > 0 {
                            out.write_all(b"\t")?;
                        }
                        genomelens_query::write_column(col, &rec, &mut out)?;
                    }
                    out.write_all(b"\n")?;
                }
            }
        }
    }

    out.flush()?;
    Ok(())
}

fn cmd_head_vcf(reader: impl BufRead, n: usize) -> genomelens_core::Result<()> {
    let mut vcf = VcfReader::new(reader);
    let _ = vcf.read_header()?;

    let mut table = Table::new();
    table.load_preset(UTF8_FULL);
    table.set_header(vec![
        "#", "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "TYPE",
    ]);

    let mut count = 0;
    while let Some(rec) = vcf.next_record()? {
        if count >= n {
            break;
        }
        table.add_row(vec![
            (count + 1).to_string(),
            std::str::from_utf8(rec.chrom().unwrap_or(b"?")).unwrap_or("?").to_string(),
            std::str::from_utf8(rec.pos().unwrap_or(b"?")).unwrap_or("?").to_string(),
            std::str::from_utf8(rec.id().unwrap_or(b"?")).unwrap_or("?").to_string(),
            std::str::from_utf8(rec.ref_allele().unwrap_or(b"?")).unwrap_or("?").to_string(),
            std::str::from_utf8(rec.alt_alleles().unwrap_or(b"?")).unwrap_or("?").to_string(),
            std::str::from_utf8(rec.qual().unwrap_or(b"?")).unwrap_or("?").to_string(),
            std::str::from_utf8(rec.filter().unwrap_or(b"?")).unwrap_or("?").to_string(),
            format!("{:?}", rec.variant_type().unwrap_or(VariantType::Complex)),
        ]);
        count += 1;
    }

    println!("{table}");
    Ok(())
}
