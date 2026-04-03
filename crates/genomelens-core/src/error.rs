#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("parse error at byte offset {offset}: {message}")]
    Parse { offset: usize, message: String },

    #[error("invalid FASTA: {0}")]
    InvalidFasta(String),

    #[error("invalid VCF: {0}")]
    InvalidVcf(String),

    #[error("invalid query: {0}")]
    InvalidQuery(String),
}

pub type Result<T> = std::result::Result<T, Error>;
