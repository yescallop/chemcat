use crate::*;
use owo_colors::OwoColorize;
use std::fmt::*;

impl Display for ChemEq {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        for (i, term) in self.terms.iter().enumerate() {
            if i == self.left_len {
                f.write_str(" -> ")?;
            } else if i != 0 {
                f.write_str(" + ")?;
            }
            Display::fmt(term, f)?;
        }
        Ok(())
    }
}

impl Display for Term {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result {
        fmt_term(self, f, true)
    }
}

/// Formats a term within a chemical equation recursively.
///
/// This function generally parenthesizes the lists, formats the coefficients
/// as subscripts and the charges as superscripts.
///
/// The parameter `outermost` specifies that the term is an outermost list, dropping the
/// parentheses around it and formatting its coefficient as preceding ASCII numbers instead.
fn fmt_term(term: &Term, f: &mut Formatter<'_>, outermost: bool) -> Result {
    const SUPERSCRIPT_PLUS: u32 = 0x207A;
    // const SUPERSCRIPT_MINUS: u32 = 0x207B;

    match *term {
        Term::Elem { ref name, n } => {
            f.write_str(name)?;
            fmt_coef(n, f, Some(true))?;
        }
        Term::List {
            ref list,
            n,
            charge,
        } => {
            if outermost {
                fmt_coef(n, f, None)?;
                for term in list {
                    fmt_term(term, f, false)?;
                }
            } else {
                f.write_char('(')?;
                for term in list {
                    fmt_term(term, f, false)?;
                }
                f.write_char(')')?;
                fmt_coef(n, f, Some(true))?;
            }
            if charge != 0 {
                let neg = charge < 0;
                fmt_coef(charge.wrapping_abs(), f, Some(false))?;
                f.write_char(char::from_u32(SUPERSCRIPT_PLUS + neg as u32).unwrap())?;
            }
        }
        Term::Electron => f.write_char('e')?,
    }
    Ok(())
}

/// Formats a coefficient in decimal as ASCII numbers, subscripts, or superscripts.
///
/// The parameter `script` accepts values `None` for ASCII numbers,
/// `Some(true)` for subscripts, and `Some(false)` for superscripts.
///
/// The coefficients of 1 are dropped by default and can only be kept for ASCII
/// numbers by specifying the "#" flag in the formatter.
///
/// Non-positive coefficients are only allowed for ASCII numbers and are wrapped to
/// be unsigned for subscripts and superscripts.
fn fmt_coef(n: i32, f: &mut Formatter<'_>, script: Option<bool>) -> Result {
    const SUBSCRIPT_ZERO: u16 = 0x2080;
    const SUPERSCRIPT_TABLE: [u16; 10] = [
        0x2070, 0xB9, 0xB2, 0xB3, 0x2074, 0x2075, 0x2076, 0x2077, 0x2078, 0x2079,
    ];

    if n == 1 && !(script == None && f.alternate()) {
        return Ok(());
    }
    match script {
        None => {
            if n > 0 {
                write!(f, "{}", n.green().bold())
            } else {
                write!(f, "{}", n.yellow().bold())
            }
        }
        Some(sub) => {
            let mut n = n as u32;

            let mut buf = [0u16; 10];
            let mut i = buf.len() - 1;
            loop {
                buf[i] = if sub {
                    SUBSCRIPT_ZERO + (n % 10) as u16
                } else {
                    SUPERSCRIPT_TABLE[(n % 10) as usize]
                };
                n /= 10;
                if n == 0 {
                    break;
                }
                i -= 1;
            }
            while i < buf.len() {
                f.write_char(char::from_u32(buf[i] as u32).unwrap())?;
                i += 1;
            }
            Ok(())
        }
    }
}
