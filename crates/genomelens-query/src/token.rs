use logos::Logos;

#[derive(Logos, Debug, Clone, PartialEq)]
#[logos(skip r"[ \t]+")]
pub enum Token<'a> {
    #[token("=")]
    Eq,
    #[token("!=")]
    Neq,
    #[token(">=")]
    Gte,
    #[token("<=")]
    Lte,
    #[token(">")]
    Gt,
    #[token("<")]
    Lt,
    #[token("(")]
    LParen,
    #[token(")")]
    RParen,
    #[token(".")]
    Dot,
    #[token(",")]
    Comma,
    #[token("*")]
    Star,

    #[regex(r"[0-9]+\.[0-9]+", |lex| lex.slice())]
    Float(&'a str),
    #[regex(r"[0-9]+", |lex| lex.slice())]
    Integer(&'a str),
    #[regex(r"'[^']*'", |lex| &lex.slice()[1..lex.slice().len()-1])]
    StringLit(&'a str),

    #[regex(r"[A-Za-z_][A-Za-z0-9_]*", |lex| lex.slice())]
    Ident(&'a str),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokenize_basic_query() {
        let input = "QUAL > 50 AND INFO.AF > 0.05";
        let tokens: Vec<_> = Token::lexer(input)
            .map(|r| r.unwrap())
            .collect();

        assert_eq!(
            tokens,
            vec![
                Token::Ident("QUAL"),
                Token::Gt,
                Token::Integer("50"),
                Token::Ident("AND"),
                Token::Ident("INFO"),
                Token::Dot,
                Token::Ident("AF"),
                Token::Gt,
                Token::Float("0.05"),
            ]
        );
    }

    #[test]
    fn tokenize_string_literal() {
        let input = "CHROM = 'chr1'";
        let tokens: Vec<_> = Token::lexer(input)
            .map(|r| r.unwrap())
            .collect();

        assert_eq!(
            tokens,
            vec![
                Token::Ident("CHROM"),
                Token::Eq,
                Token::StringLit("chr1"),
            ]
        );
    }

    #[test]
    fn tokenize_comparison_operators() {
        let input = "= != >= <= > <";
        let tokens: Vec<_> = Token::lexer(input)
            .map(|r| r.unwrap())
            .collect();

        assert_eq!(
            tokens,
            vec![Token::Eq, Token::Neq, Token::Gte, Token::Lte, Token::Gt, Token::Lt]
        );
    }

    #[test]
    fn tokenize_parentheses() {
        let input = "(CHROM = 'chr1' OR CHROM = 'chr2')";
        let tokens: Vec<_> = Token::lexer(input)
            .map(|r| r.unwrap())
            .collect();

        assert_eq!(tokens[0], Token::LParen);
        assert_eq!(tokens[tokens.len() - 1], Token::RParen);
    }

    #[test]
    fn tokenize_not() {
        let input = "NOT INFO.DB";
        let tokens: Vec<_> = Token::lexer(input)
            .map(|r| r.unwrap())
            .collect();

        assert_eq!(
            tokens,
            vec![
                Token::Ident("NOT"),
                Token::Ident("INFO"),
                Token::Dot,
                Token::Ident("DB"),
            ]
        );
    }
}
