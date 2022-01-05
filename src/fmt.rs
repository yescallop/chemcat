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
        fmt_term(self, f, false)
    }
}

fn fmt_term(term: &Term, f: &mut Formatter<'_>, delimited: bool) -> Result {
    const SUPER_PLUS: u32 = 0x207A;
    // const SUPER_MINUS: u32 = 0x207B;

    match term {
        Term::Elem { name, n } => {
            f.write_str(name)?;
            fmt_coef(*n, f, Some(true))?;
        }
        Term::List { list, n, charge } => {
            if delimited {
                f.write_char('(')?;
                for term in list {
                    fmt_term(term, f, true)?;
                }
                f.write_char(')')?;
                fmt_coef(*n, f, Some(true))?;
            } else {
                fmt_coef(*n, f, None)?;
                for term in list {
                    fmt_term(term, f, true)?;
                }
            }
            if *charge != 0 {
                let neg = *charge < 0;
                fmt_coef(charge.abs(), f, Some(false))?;
                f.write_char(char::from_u32(SUPER_PLUS + neg as u32).unwrap())?;
            }
        }
        Term::Electron => f.write_char('e')?,
    }
    Ok(())
}

// None for no script, Some(true) for subscript, Some(false) for superscript.
fn fmt_coef(mut n: i32, f: &mut Formatter<'_>, script: Option<bool>) -> Result {
    const SUB_START: u16 = 0x2080;
    const SUPER_TABLE: [u16; 10] = [
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
            let mut buf = [0u16; 10];
            let mut i = buf.len() - 1;
            loop {
                buf[i] = if sub {
                    SUB_START + (n % 10) as u16
                } else {
                    SUPER_TABLE[(n % 10) as usize]
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
