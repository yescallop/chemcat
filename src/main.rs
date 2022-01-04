use chemcat::Solution;
use chumsky::prelude::*;
use owo_colors::OwoColorize;
use std::{
    collections::HashMap,
    fmt::{self, Write},
    io, str,
};

fn main() {
    let mut buf = String::new();
    let eq_parser = eq_parser();
    loop {
        println!("----- CHEMCAT -----\nNya! Feed me an equation:");
        io::stdin()
            .read_line(&mut buf)
            .expect("failed to read line");
        // TODO: Elegant exit.

        let eq = eq_parser.parse(buf.trim());
        match eq {
            Ok(eq) => try_balance(eq),
            // TODO: Better error messages.
            Err(e) => println!("Error: {:?}\n", e),
        }

        buf.clear();
    }
}

fn try_balance(mut eq: ChemEq) {
    println!("\nInput interpretation:\n{}", eq);

    let term_cnt = eq.terms.len();
    let mut coef_per_elem = HashMap::<&str, i32>::new();
    let mut eq_per_elem = HashMap::<&str, Vec<i32>>::new();

    println!("\n----- COLLAPSING TERMS -----");

    for (i, term) in eq.terms.iter().enumerate() {
        term.collapse(&mut coef_per_elem, 1);
        println!("{} -> {:?}", term, coef_per_elem);

        let left = i < eq.left_len;
        for (elem, n) in coef_per_elem.drain() {
            let row = eq_per_elem.entry(elem).or_insert_with(|| vec![0; term_cnt]);
            row[i] = if left { n } else { -n };
        }
    }

    println!("\n----- BUILDING MATRIX -----");

    let mut matrix: Vec<_> = eq_per_elem.into_values().collect();
    for row in &matrix {
        println!("{:?}", row);
    }

    println!("\n----- GAUSSIAN ELIMINATION -----");

    // In Gauss we trust.
    chemcat::row_reduce(&mut matrix);
    for row in &matrix {
        println!("{:?}", row);
    }

    println!("\n----- BALANCING -----");

    let sol = chemcat::solve(matrix);
    match sol {
        Solution::None => println!("Nya? This equation has no solution!"),
        Solution::Unique(sol) => {
            let mut ill_coef = false;
            eq.terms
                .iter_mut()
                .zip(sol.into_iter())
                .for_each(|(term, coef)| match term {
                    Term::List { n, .. } => {
                        *n = coef;
                        if coef <= 0 {
                            ill_coef = true;
                        }
                    }
                    _ => unreachable!(),
                });

            // The `#` flag keeps preceding coefficients of 1.
            // FIXME: Let users choose whether to keep them or not.
            println!(
                "{}! Here's the solution:\n{:#}",
                if ill_coef { "Ugh" } else { "Yay" },
                eq
            );

            if ill_coef {
                println!(
                    "\n{}{}",
                    "warning".yellow().bold(),
                    ": Chemcat hates non-positive coefficients".white().bold()
                );
            }
        }
        Solution::Infinite(n) => println!(
            "(UwU) Way too complex for Chemcat!\n\
            It's a linear combination of {} independent solutions.",
            n.green().bold()
        ),
    }

    println!();
}

#[derive(Debug)]
struct ChemEq {
    terms: Vec<Term>,
    left_len: usize,
}

impl fmt::Display for ChemEq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, term) in self.terms.iter().enumerate() {
            if i == self.left_len {
                f.write_str(" -> ")?;
            } else if i != 0 {
                f.write_str(" + ")?;
            }
            fmt::Display::fmt(term, f)?;
        }
        Ok(())
    }
}

#[derive(Debug)]
enum Term {
    Elem {
        name: String,
        n: i32,
    },
    List {
        list: Vec<Term>,
        n: i32,
        charge: i32,
    },
}

impl Term {
    const CHARGE_NAME: &'static str = "c";

    fn collapse<'a>(&'a self, map: &mut HashMap<&'a str, i32>, m: i32) {
        match self {
            Term::Elem { name, n } => {
                let cnt = map.entry(name).or_insert(0);
                *cnt += m * n;
            }
            Term::List { list, n, charge } => {
                list.iter().for_each(|term| term.collapse(map, m * n));
                if *charge != 0 {
                    let cnt = map.entry(Term::CHARGE_NAME).or_insert(0);
                    *cnt += m * charge;
                }
            }
        }
    }

    fn fmt_impl(&self, f: &mut fmt::Formatter<'_>, delimited: bool) -> fmt::Result {
        const SUPER_PLUS: u32 = 0x207A;
        // const SUPER_MINUS: u32 = 0x207B;

        match self {
            Term::Elem { name, n } => {
                f.write_str(name)?;
                Term::fmt_coef(*n, f, Some(true))?;
            }
            Term::List { list, n, charge } => {
                if delimited {
                    f.write_char('(')?;
                    for term in list {
                        term.fmt_impl(f, true)?;
                    }
                    f.write_char(')')?;
                    Term::fmt_coef(*n, f, Some(true))?;
                } else {
                    Term::fmt_coef(*n, f, None)?;
                    for term in list {
                        term.fmt_impl(f, true)?;
                    }
                }
                if *charge != 0 {
                    let neg = *charge < 0;
                    Term::fmt_coef(charge.abs(), f, Some(false))?;
                    f.write_char(char::from_u32(SUPER_PLUS + neg as u32).unwrap())?;
                }
            }
        }
        Ok(())
    }

    // None for no script, Some(true) for subscript, Some(false) for superscript.
    fn fmt_coef(mut n: i32, f: &mut fmt::Formatter<'_>, script: Option<bool>) -> fmt::Result {
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
}

impl fmt::Display for Term {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.fmt_impl(f, false)
    }
}

fn eq_parser() -> impl Parser<char, ChemEq, Error = Simple<char>> {
    let coef = text::int(10)
        .try_map(|s: String, span| match s.parse::<i32>() {
            Ok(0) => Err(Simple::custom(span, "zero coefficient")),
            Ok(n) => Ok(n),
            Err(e) => Err(Simple::custom(span, e)),
        })
        .or(empty().to(1));

    let elem_name = filter(char::is_ascii_uppercase)
        .map(Some)
        .chain::<char, _, _>(filter(char::is_ascii_lowercase).repeated())
        .collect();

    let elem = elem_name
        .then(coef.clone())
        .map(|(name, n)| Term::Elem { name, n });

    let sub_term = recursive(|list_or_elem| {
        let list_content = list_or_elem.repeated().at_least(1);

        list_content
            .clone()
            .delimited_by('(', ')')
            .or(list_content.delimited_by('[', ']'))
            .then(coef.clone())
            .map(|(list, n)| Term::List { list, n, charge: 0 })
            .or(elem)
    })
    .repeated()
    .at_least(1);

    let sub_term_separator = one_of("·.*");

    let charge = coef
        .clone()
        .then(just('+').to(false).or(just('-').to(true)))
        .delimited_by('(', ')')
        .map(|(n, neg)| if neg { -n } else { n })
        .or(empty().to(0));

    let term_without_coef = sub_term
        .clone()
        .then(
            sub_term_separator
                .ignore_then(coef.clone().then(sub_term))
                .repeated(),
        )
        .map(|(mut head, rest)| {
            head.extend(
                rest.into_iter()
                    .map(|(n, list)| Term::List { list, n, charge: 0 }),
            );
            head
        })
        .then(charge)
        .map(|(list, charge)| Term::List { list, n: 1, charge });

    let term = coef.ignore_then(term_without_coef);

    let side = term.separated_by(just('+').padded()).at_least(1);

    let yields = choice((
        just('=').repeated().at_least(1).ignored(),
        just('-').repeated().at_least(1).then(just('>')).ignored(),
        just('→').ignored(),
    ));

    side.clone()
        .then_ignore(yields.padded())
        .then(side)
        .map(|(mut left, right)| {
            let left_len = left.len();
            left.extend(right);
            ChemEq {
                terms: left,
                left_len,
            }
        })
        .then_ignore(end())
}
