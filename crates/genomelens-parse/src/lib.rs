pub mod ascii;
pub mod fasta;
pub mod reader;
pub mod vcf;

pub use fasta::{FastaReader, FastaRecord, FastaRecordStats};
pub use reader::{detect_format, open_transparent, FileFormat};
pub use vcf::{VcfHeader, VcfReader, VcfRecord};
