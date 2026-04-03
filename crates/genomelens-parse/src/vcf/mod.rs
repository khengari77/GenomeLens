mod header;
pub mod record;
mod streaming;

pub use header::VcfHeader;
pub use record::VcfRecord;
pub use streaming::VcfReader;
