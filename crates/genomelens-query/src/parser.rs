use genomelens_core::{Error, Result};
use logos::Logos;

use crate::ast::{CmpOp, Column, Expr, FixedColumn, Literal, Query, Select};
use crate::token::Token;

struct Parser<'a> {
    tokens: Vec<Token<'a>>,
    pos: usize,
}

impl<'a> Parser<'a> {
    fn new(input: &'a str) -> Result<Self> {
        let tokens: std::result::Result<Vec<_>, _> = Token::lexer(input).collect();
        let tokens = tokens.map_err(|_| Error::InvalidQuery("invalid token".to_string()))?;
        Ok(Self { tokens, pos: 0 })
    }

    fn peek(&self) -> Option<&Token<'a>> {
        self.tokens.get(self.pos)
    }

    fn next(&mut self) -> Option<Token<'a>> {
        let tok = self.tokens.get(self.pos).cloned();
        if tok.is_some() {
            self.pos += 1;
        }
        tok
    }

    fn expect(&mut self, expected: &Token<'_>) -> Result<()> {
        match self.next() {
            Some(ref tok) if tok == expected => Ok(()),
            Some(tok) => Err(Error::InvalidQuery(format!(
                "expected {:?}, found {:?}",
                expected, tok
            ))),
            None => Err(Error::InvalidQuery(format!(
                "expected {:?}, found end of input",
                expected
            ))),
        }
    }

    fn at_end(&self) -> bool {
        self.pos >= self.tokens.len()
    }

    // --- Query parsing (SELECT ... WHERE ...) ---

    fn parse_query(&mut self) -> Result<Query> {
        let select = if self.peek_is_keyword("SELECT") {
            self.next(); // consume SELECT
            self.parse_select()?
        } else {
            Select::All
        };

        let filter = if self.peek_is_keyword("WHERE") {
            self.next(); // consume WHERE
            Some(self.parse_expr(0)?)
        } else if !self.at_end() && select == Select::All {
            // No SELECT, no WHERE — treat entire input as a filter expression
            Some(self.parse_expr(0)?)
        } else {
            None
        };

        Ok(Query { select, filter })
    }

    fn parse_select(&mut self) -> Result<Select> {
        // Check for SELECT *
        if let Some(Token::Star) = self.peek() {
            self.next();
            return Ok(Select::All);
        }

        let mut columns = Vec::new();
        columns.push(self.parse_column()?);

        loop {
            match self.peek() {
                Some(Token::Ident(kw)) if kw.eq_ignore_ascii_case("WHERE") => break,
                Some(Token::Comma) => { self.next(); continue; }
                Some(Token::Ident(_)) => {}
                _ => break,
            }
            columns.push(self.parse_column()?);
        }

        Ok(Select::Columns(columns))
    }

    fn parse_column(&mut self) -> Result<Column> {
        match self.next() {
            Some(Token::Ident(name)) => self.resolve_column(name),
            other => Err(Error::InvalidQuery(format!(
                "expected column name, found {:?}",
                other
            ))),
        }
    }

    fn peek_is_keyword(&self, keyword: &str) -> bool {
        matches!(self.peek(), Some(Token::Ident(s)) if s.eq_ignore_ascii_case(keyword))
    }

    // --- Expression parsing (Pratt) ---

    fn parse_expr(&mut self, min_bp: u8) -> Result<Expr> {
        let mut lhs = self.parse_atom()?;

        loop {
            let Some(tok) = self.peek() else { break };

            match tok {
                Token::Ident(kw) if kw.eq_ignore_ascii_case("AND") => {
                    let bp = 2;
                    if bp <= min_bp {
                        break;
                    }
                    self.next();
                    let rhs = self.parse_expr(bp)?;
                    lhs = Expr::And(Box::new(lhs), Box::new(rhs));
                }
                Token::Ident(kw) if kw.eq_ignore_ascii_case("OR") => {
                    let bp = 1;
                    if bp <= min_bp {
                        break;
                    }
                    self.next();
                    let rhs = self.parse_expr(bp)?;
                    lhs = Expr::Or(Box::new(lhs), Box::new(rhs));
                }
                _ => break,
            }
        }

        Ok(lhs)
    }

    fn parse_atom(&mut self) -> Result<Expr> {
        let tok = self.next().ok_or_else(|| {
            Error::InvalidQuery("unexpected end of input".to_string())
        })?;

        match tok {
            Token::Ident(kw) if kw.eq_ignore_ascii_case("NOT") => {
                let inner = self.parse_atom()?;
                Ok(Expr::Not(Box::new(inner)))
            }
            Token::LParen => {
                let expr = self.parse_expr(0)?;
                self.expect(&Token::RParen)?;
                Ok(expr)
            }
            Token::Ident(name) => {
                let column = self.resolve_column(name)?;
                self.parse_comparison_or_flag(column)
            }
            other => Err(Error::InvalidQuery(format!(
                "unexpected token: {:?}",
                other
            ))),
        }
    }

    fn resolve_column(&mut self, name: &str) -> Result<Column> {
        if name.eq_ignore_ascii_case("INFO") {
            if self.peek() == Some(&Token::Dot) {
                self.next(); // consume dot
                match self.next() {
                    Some(Token::Ident(field)) => Ok(Column::Info(field.to_string())),
                    _ => Err(Error::InvalidQuery(
                        "expected field name after INFO.".to_string(),
                    )),
                }
            } else {
                Err(Error::InvalidQuery(
                    "INFO must be followed by .FIELD".to_string(),
                ))
            }
        } else {
            let fixed = match name.to_ascii_uppercase().as_str() {
                "CHROM" => FixedColumn::Chrom,
                "POS" => FixedColumn::Pos,
                "QUAL" => FixedColumn::Qual,
                "FILTER" => FixedColumn::Filter,
                "REF" => FixedColumn::Ref,
                "TYPE" => FixedColumn::Type,
                _ => {
                    return Err(Error::InvalidQuery(format!(
                        "unknown column: {}",
                        name
                    )));
                }
            };
            Ok(Column::Fixed(fixed))
        }
    }

    fn parse_comparison_or_flag(&mut self, column: Column) -> Result<Expr> {
        let op = match self.peek() {
            Some(Token::Eq) => Some(CmpOp::Eq),
            Some(Token::Neq) => Some(CmpOp::Neq),
            Some(Token::Gt) => Some(CmpOp::Gt),
            Some(Token::Lt) => Some(CmpOp::Lt),
            Some(Token::Gte) => Some(CmpOp::Gte),
            Some(Token::Lte) => Some(CmpOp::Lte),
            _ => None,
        };

        if let Some(op) = op {
            self.next(); // consume operator
            let value = self.parse_literal()?;
            Ok(Expr::Compare { column, op, value })
        } else {
            match column {
                Column::Info(field) => Ok(Expr::FlagCheck { field }),
                Column::Fixed(_) => Err(Error::InvalidQuery(
                    "fixed columns require a comparison operator".to_string(),
                )),
            }
        }
    }

    fn parse_literal(&mut self) -> Result<Literal> {
        match self.next() {
            Some(Token::Integer(s)) => {
                let v = s.parse::<i64>().map_err(|_| {
                    Error::InvalidQuery(format!("invalid integer: {}", s))
                })?;
                Ok(Literal::Integer(v))
            }
            Some(Token::Float(s)) => {
                let v = s.parse::<f64>().map_err(|_| {
                    Error::InvalidQuery(format!("invalid float: {}", s))
                })?;
                Ok(Literal::Float(v))
            }
            Some(Token::StringLit(s)) => Ok(Literal::Str(s.to_string())),
            other => Err(Error::InvalidQuery(format!(
                "expected literal value, found {:?}",
                other
            ))),
        }
    }
}

/// Parse a filter expression (WHERE clause only). Backward-compatible API.
pub fn parse(input: &str) -> Result<Expr> {
    let mut parser = Parser::new(input)?;
    let expr = parser.parse_expr(0)?;
    if !parser.at_end() {
        return Err(Error::InvalidQuery(format!(
            "unexpected token after expression: {:?}",
            parser.tokens[parser.pos]
        )));
    }
    Ok(expr)
}

/// Parse a full query: `SELECT cols WHERE expr`, `WHERE expr`, or bare filter expression.
pub fn parse_query(input: &str) -> Result<Query> {
    let mut parser = Parser::new(input)?;
    let query = parser.parse_query()?;
    if !parser.at_end() {
        return Err(Error::InvalidQuery(format!(
            "unexpected token after query: {:?}",
            parser.tokens[parser.pos]
        )));
    }
    Ok(query)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_simple_comparison() {
        let expr = parse("QUAL > 50").unwrap();
        assert_eq!(
            expr,
            Expr::Compare {
                column: Column::Fixed(FixedColumn::Qual),
                op: CmpOp::Gt,
                value: Literal::Integer(50),
            }
        );
    }

    #[test]
    fn parse_string_comparison() {
        let expr = parse("CHROM = 'chr1'").unwrap();
        assert_eq!(
            expr,
            Expr::Compare {
                column: Column::Fixed(FixedColumn::Chrom),
                op: CmpOp::Eq,
                value: Literal::Str("chr1".to_string()),
            }
        );
    }

    #[test]
    fn parse_info_comparison() {
        let expr = parse("INFO.AF > 0.05").unwrap();
        assert_eq!(
            expr,
            Expr::Compare {
                column: Column::Info("AF".to_string()),
                op: CmpOp::Gt,
                value: Literal::Float(0.05),
            }
        );
    }

    #[test]
    fn parse_info_flag() {
        let expr = parse("INFO.DB").unwrap();
        assert_eq!(
            expr,
            Expr::FlagCheck {
                field: "DB".to_string(),
            }
        );
    }

    #[test]
    fn parse_and() {
        let expr = parse("QUAL > 50 AND CHROM = 'chr1'").unwrap();
        match expr {
            Expr::And(lhs, rhs) => {
                assert!(matches!(*lhs, Expr::Compare { .. }));
                assert!(matches!(*rhs, Expr::Compare { .. }));
            }
            _ => panic!("expected And"),
        }
    }

    #[test]
    fn parse_or() {
        let expr = parse("CHROM = 'chr1' OR CHROM = 'chr2'").unwrap();
        assert!(matches!(expr, Expr::Or(_, _)));
    }

    #[test]
    fn parse_not() {
        let expr = parse("NOT INFO.DB").unwrap();
        match expr {
            Expr::Not(inner) => {
                assert!(matches!(*inner, Expr::FlagCheck { .. }));
            }
            _ => panic!("expected Not"),
        }
    }

    #[test]
    fn parse_precedence_and_over_or() {
        let expr = parse("CHROM = 'chr1' OR QUAL > 50 AND POS > 100").unwrap();
        match expr {
            Expr::Or(lhs, rhs) => {
                assert!(matches!(*lhs, Expr::Compare { .. }));
                assert!(matches!(*rhs, Expr::And(_, _)));
            }
            _ => panic!("expected Or at top level"),
        }
    }

    #[test]
    fn parse_parentheses() {
        let expr = parse("(CHROM = 'chr1' OR CHROM = 'chr2') AND QUAL > 50").unwrap();
        match expr {
            Expr::And(lhs, rhs) => {
                assert!(matches!(*lhs, Expr::Or(_, _)));
                assert!(matches!(*rhs, Expr::Compare { .. }));
            }
            _ => panic!("expected And at top level"),
        }
    }

    #[test]
    fn parse_type_comparison() {
        let expr = parse("TYPE = 'SNV'").unwrap();
        assert_eq!(
            expr,
            Expr::Compare {
                column: Column::Fixed(FixedColumn::Type),
                op: CmpOp::Eq,
                value: Literal::Str("SNV".to_string()),
            }
        );
    }

    #[test]
    fn error_unknown_column() {
        assert!(parse("FOOBAR > 5").is_err());
    }

    #[test]
    fn error_missing_operator_on_fixed() {
        assert!(parse("CHROM").is_err());
    }

    #[test]
    fn error_trailing_tokens() {
        assert!(parse("QUAL > 50 QUAL").is_err());
    }

    #[test]
    fn error_empty_input() {
        assert!(parse("").is_err());
    }

    #[test]
    fn case_insensitive_keywords() {
        let expr = parse("qual > 50 and chrom = 'chr1'").unwrap();
        assert!(matches!(expr, Expr::And(_, _)));
    }

    // --- Query parsing tests ---

    #[test]
    fn query_bare_filter() {
        let q = parse_query("QUAL > 50").unwrap();
        assert_eq!(q.select, Select::All);
        assert!(q.filter.is_some());
    }

    #[test]
    fn query_select_where() {
        let q = parse_query("SELECT CHROM POS QUAL WHERE QUAL > 50").unwrap();
        assert_eq!(
            q.select,
            Select::Columns(vec![
                Column::Fixed(FixedColumn::Chrom),
                Column::Fixed(FixedColumn::Pos),
                Column::Fixed(FixedColumn::Qual),
            ])
        );
        assert!(q.filter.is_some());
    }

    #[test]
    fn query_select_info_column() {
        let q = parse_query("SELECT CHROM INFO.AF WHERE QUAL > 50").unwrap();
        assert_eq!(
            q.select,
            Select::Columns(vec![
                Column::Fixed(FixedColumn::Chrom),
                Column::Info("AF".to_string()),
            ])
        );
    }

    #[test]
    fn query_select_only() {
        let q = parse_query("SELECT CHROM POS").unwrap();
        assert_eq!(
            q.select,
            Select::Columns(vec![
                Column::Fixed(FixedColumn::Chrom),
                Column::Fixed(FixedColumn::Pos),
            ])
        );
        assert!(q.filter.is_none());
    }

    #[test]
    fn query_where_only() {
        let q = parse_query("WHERE QUAL > 50").unwrap();
        assert_eq!(q.select, Select::All);
        assert!(q.filter.is_some());
    }

    #[test]
    fn query_select_star() {
        // * is not in our ident regex, so we handle it differently.
        // For now, SELECT without * just lists columns until WHERE.
        // SELECT * would need a special token. Let's test SELECT with columns.
        let q = parse_query("SELECT CHROM WHERE QUAL > 50").unwrap();
        assert_eq!(
            q.select,
            Select::Columns(vec![Column::Fixed(FixedColumn::Chrom)])
        );
    }

    #[test]
    fn query_select_with_commas() {
        let q = parse_query("SELECT CHROM, POS, QUAL WHERE QUAL > 50").unwrap();
        assert_eq!(
            q.select,
            Select::Columns(vec![
                Column::Fixed(FixedColumn::Chrom),
                Column::Fixed(FixedColumn::Pos),
                Column::Fixed(FixedColumn::Qual),
            ])
        );
        assert!(q.filter.is_some());
    }

    #[test]
    fn query_select_commas_no_where() {
        let q = parse_query("SELECT CHROM, POS").unwrap();
        assert_eq!(
            q.select,
            Select::Columns(vec![
                Column::Fixed(FixedColumn::Chrom),
                Column::Fixed(FixedColumn::Pos),
            ])
        );
        assert!(q.filter.is_none());
    }

    #[test]
    fn query_select_star_where() {
        let q = parse_query("SELECT * WHERE QUAL > 50").unwrap();
        assert_eq!(q.select, Select::All);
        assert!(q.filter.is_some());
    }

    #[test]
    fn query_select_star_no_where() {
        let q = parse_query("SELECT *").unwrap();
        assert_eq!(q.select, Select::All);
        assert!(q.filter.is_none());
    }
}
