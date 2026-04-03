#[derive(Debug, Clone, PartialEq)]
pub enum FixedColumn {
    Chrom,
    Pos,
    Qual,
    Filter,
    Ref,
    Type,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Column {
    Fixed(FixedColumn),
    Info(String),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CmpOp {
    Eq,
    Neq,
    Gt,
    Lt,
    Gte,
    Lte,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Literal {
    Integer(i64),
    Float(f64),
    Str(String),
}

#[derive(Debug, Clone, PartialEq)]
pub enum Expr {
    Compare {
        column: Column,
        op: CmpOp,
        value: Literal,
    },
    FlagCheck {
        field: String,
    },
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
    Not(Box<Expr>),
}

/// Column projection: which columns to include in output.
#[derive(Debug, Clone, PartialEq)]
pub enum Select {
    /// Output all columns (raw VCF line or all fixed columns).
    All,
    /// Output specific columns in order.
    Columns(Vec<Column>),
}

/// A complete query: optional projection + optional filter.
#[derive(Debug, Clone, PartialEq)]
pub struct Query {
    pub select: Select,
    pub filter: Option<Expr>,
}
