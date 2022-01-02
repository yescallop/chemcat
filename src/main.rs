use chemcat::Solution;
use chumsky::prelude::*;
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

        let eq = eq_parser.parse(buf.trim());
        match eq {
            Ok(eq) => try_balance(eq),
            Err(e) => println!("Error: {:?}\n", e),
        }

        buf.clear();
    }
}

fn try_balance(mut eq: ChemEq) {
    println!("\nInput interpretation:\n{}", eq);

    let term_cnt = eq.terms.len();
    let mut eq_per_elem = HashMap::<&str, Vec<i32>>::new();

    println!("\n----- COLLAPSING TERMS -----");

    for (i, term) in eq.terms.iter().enumerate() {
        let collapsed = term.collapse();
        println!("{} -> {:?}", term, collapsed);

        let left = i < eq.left_len;
        for (elem, n) in collapsed {
            let row = eq_per_elem.entry(elem).or_insert_with(|| vec![0; term_cnt]);
            row[i] = if left { n as i32 } else { -(n as i32) };
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
            eq.terms
                .iter_mut()
                .zip(sol.into_iter())
                .for_each(|(term, coef)| match term {
                    Term::List(_, n) => *n = coef,
                    _ => unreachable!(),
                });
            println!("Yay! Here's the solution:\n{}", eq);
        }
        Solution::Infinite(n) => println!(
            "(UwU) Way too complex for Chemcat!\n\
            It's a linear combination of {} independent solutions.",
            n
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
            write!(f, "{}", term)?;
        }
        Ok(())
    }
}

#[derive(Debug)]
enum Term {
    Elem(String, u32),
    List(Vec<Term>, u32),
}

impl Term {
    fn collapse(&self) -> HashMap<&str, u32> {
        let mut map = HashMap::new();
        self.collapse_impl(&mut map, 1);
        map
    }

    fn collapse_impl<'a>(&'a self, map: &mut HashMap<&'a str, u32>, m: u32) {
        match self {
            Term::Elem(elem, n) => {
                let cnt = map.entry(elem).or_insert(0);
                *cnt += m * n;
            }
            Term::List(list, n) => list.iter().for_each(|term| term.collapse_impl(map, m * n)),
        }
    }

    fn fmt_impl(&self, f: &mut fmt::Formatter<'_>, delimited: bool) -> fmt::Result {
        match self {
            Term::Elem(elem, n) => {
                f.write_str(elem)?;
                Self::fmt_coef(*n, f, true)?;
            }
            Term::List(list, n) => {
                if delimited {
                    f.write_char('(')?;
                    for term in list {
                        term.fmt_impl(f, true)?;
                    }
                    f.write_char(')')?;
                    Self::fmt_coef(*n, f, true)?;
                } else {
                    Self::fmt_coef(*n, f, false)?;
                    for term in list {
                        term.fmt_impl(f, true)?;
                    }
                }
            }
        }
        Ok(())
    }

    fn fmt_coef(mut n: u32, f: &mut fmt::Formatter<'_>, subscript: bool) -> fmt::Result {
        if n == 1 {
            Ok(())
        } else if !subscript {
            write!(f, "{}", n)
        } else {
            let mut buf = [0u8; 30];
            let mut i = 27;
            loop {
                buf[i] = 0xE2;
                buf[i + 1] = 0x82;
                buf[i + 2] = 0x80 + (n % 10) as u8;
                n /= 10;
                if n == 0 {
                    break;
                }
                i -= 3;
            }
            // SAFETY: Cats won't care.
            let s = unsafe { str::from_utf8_unchecked(&buf[i..]) };
            f.write_str(s)
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
        .try_map(|s: String, span| match s.parse::<u32>() {
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
        .map(|(elem, n)| Term::Elem(elem, n));

    let delimited_list = recursive(|list| {
        let content = list.repeated().at_least(1);

        content
            .clone()
            .delimited_by('(', ')')
            .or(content.delimited_by('[', ']'))
            .then(coef.clone())
            .map(|(list, n)| Term::List(list, n))
            .or(elem.clone())
    });

    let term_one = elem.or(delimited_list).repeated().at_least(1);

    let separator = one_of("·.*");

    let term = coef
        .clone()
        .ignore_then(term_one.clone())
        .then(separator.ignore_then(coef.then(term_one)).repeated())
        .map(|(mut head, rest)| {
            head.extend(rest.into_iter().map(|(n, list)| Term::List(list, n)));
            Term::List(head, 1)
        });

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
