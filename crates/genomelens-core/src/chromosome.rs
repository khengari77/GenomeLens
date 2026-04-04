/// A chromosome or contig identifier.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord, serde::Serialize)]
#[serde(transparent)]
pub struct Chromosome(pub String);

impl Chromosome {
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl std::borrow::Borrow<str> for Chromosome {
    fn borrow(&self) -> &str {
        &self.0
    }
}

impl<S: Into<String>> From<S> for Chromosome {
    fn from(s: S) -> Self {
        Self(s.into())
    }
}
