//! Fast ASCII-only parsing for genomic data fields.
//! Bypasses UTF-8 validation and std lib parsing overhead since genomic coordinates
//! and quality scores are guaranteed ASCII digits.

use genomelens_core::{Error, Result};

/// Parse ASCII digit bytes directly to usize. No UTF-8 validation.
/// Genomic positions (POS) are always positive integers in ASCII.
pub fn parse_ascii_usize(bytes: &[u8]) -> Result<usize> {
    if bytes.is_empty() {
        return Err(Error::Parse {
            offset: 0,
            message: "empty integer field".to_string(),
        });
    }
    let mut val = 0usize;
    for &b in bytes {
        let digit = b.wrapping_sub(b'0');
        if digit > 9 {
            return Err(Error::Parse {
                offset: 0,
                message: format!("invalid digit byte: 0x{:02x}", b),
            });
        }
        val = val
            .checked_mul(10)
            .and_then(|v| v.checked_add(digit as usize))
            .ok_or_else(|| Error::Parse {
                offset: 0,
                message: "integer overflow".to_string(),
            })?;
    }
    Ok(val)
}

/// Parse ASCII digit bytes to i64. Handles optional leading `-`.
pub fn parse_ascii_i64(bytes: &[u8]) -> Result<i64> {
    if bytes.is_empty() {
        return Err(Error::Parse {
            offset: 0,
            message: "empty integer field".to_string(),
        });
    }
    let mut pos = 0;
    let negative = if bytes[0] == b'-' {
        pos = 1;
        true
    } else {
        false
    };
    if pos >= bytes.len() {
        return Err(Error::Parse {
            offset: 0,
            message: "empty integer field after sign".to_string(),
        });
    }
    let mut val = 0i64;
    while pos < bytes.len() {
        let digit = bytes[pos].wrapping_sub(b'0');
        if digit > 9 {
            return Err(Error::Parse {
                offset: pos,
                message: format!("invalid digit byte: 0x{:02x}", bytes[pos]),
            });
        }
        val = val
            .checked_mul(10)
            .and_then(|v| v.checked_add(digit as i64))
            .ok_or_else(|| Error::Parse {
                offset: 0,
                message: "integer overflow".to_string(),
            })?;
        pos += 1;
    }
    Ok(if negative { -val } else { val })
}

/// Parse ASCII bytes to f64. Delegates to the standard library's highly optimized
/// string-to-float parser (Eisel-Lemire / Ryu) for mathematical correctness.
pub fn parse_ascii_f64(bytes: &[u8]) -> Result<f64> {
    if bytes.is_empty() {
        return Err(Error::Parse {
            offset: 0,
            message: "empty float field".to_string(),
        });
    }
    // VCF fields are guaranteed ASCII, so from_utf8 is essentially free.
    let s = std::str::from_utf8(bytes).map_err(|_| Error::Parse {
        offset: 0,
        message: "invalid UTF-8 in float field".to_string(),
    })?;
    s.parse::<f64>().map_err(|_| Error::Parse {
        offset: 0,
        message: format!("invalid float: {}", s),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_usize_basic() {
        assert_eq!(parse_ascii_usize(b"0").unwrap(), 0);
        assert_eq!(parse_ascii_usize(b"1").unwrap(), 1);
        assert_eq!(parse_ascii_usize(b"12345").unwrap(), 12345);
        assert_eq!(parse_ascii_usize(b"999999999").unwrap(), 999999999);
    }

    #[test]
    fn parse_usize_empty() {
        assert!(parse_ascii_usize(b"").is_err());
    }

    #[test]
    fn parse_usize_invalid() {
        assert!(parse_ascii_usize(b"12a45").is_err());
        assert!(parse_ascii_usize(b"-1").is_err());
        assert!(parse_ascii_usize(b" 5").is_err());
    }

    #[test]
    fn parse_f64_integer() {
        let v = parse_ascii_f64(b"50").unwrap();
        assert!((v - 50.0).abs() < f64::EPSILON);
    }

    #[test]
    fn parse_f64_decimal() {
        let v = parse_ascii_f64(b"99.9").unwrap();
        assert!((v - 99.9).abs() < 1e-10);
    }

    #[test]
    fn parse_f64_small_decimal() {
        let v = parse_ascii_f64(b"0.001").unwrap();
        assert!((v - 0.001).abs() < 1e-10);
    }

    #[test]
    fn parse_f64_negative() {
        let v = parse_ascii_f64(b"-10.5").unwrap();
        assert!((v + 10.5).abs() < 1e-10);
    }

    #[test]
    fn parse_f64_zero() {
        let v = parse_ascii_f64(b"0").unwrap();
        assert!((v).abs() < f64::EPSILON);
    }

    #[test]
    fn parse_f64_empty() {
        assert!(parse_ascii_f64(b"").is_err());
    }

    #[test]
    fn parse_f64_invalid() {
        assert!(parse_ascii_f64(b"1.2.3").is_err());
        assert!(parse_ascii_f64(b"abc").is_err());
    }
}

#[cfg(test)]
mod prop_tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn usize_roundtrip(val in 0..=u32::MAX as usize) {
            let s = val.to_string();
            let parsed = parse_ascii_usize(s.as_bytes()).unwrap();
            prop_assert_eq!(parsed, val);
        }

        #[test]
        fn f64_integer_roundtrip(val in 0..=999999u32) {
            let s = val.to_string();
            let parsed = parse_ascii_f64(s.as_bytes()).unwrap();
            prop_assert!((parsed - val as f64).abs() < f64::EPSILON);
        }
    }
}
