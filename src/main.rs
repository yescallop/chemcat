#![warn(rust_2018_idioms)]

use chemcat::*;
use chumsky::prelude::*;
use owo_colors::OwoColorize;
use std::{io, process};

fn main() {
    ctrlc::set_handler(|| process::exit(0)).expect("Error setting Ctrl-C handler");

    let mut buf = String::new();
    loop {
        println!("----- CHEMCAT (UwU) -----\nNya! Feed me an equation:");
        io::stdin()
            .read_line(&mut buf)
            .expect("failed to read line");

        let eq = eq_parser().parse(buf.trim()).into_result();
        match eq {
            Ok(eq) => try_balance(eq),
            // TODO: Better error messages.
            Err(e) => println!("Error: {:?}\n", e),
        }

        buf.clear();
    }
}

fn try_balance(mut eq: ChemEq<'_>) {
    println!("\nInput interpretation:\n{}", eq);

    println!("\n----- BUILDING MATRIX -----");

    let mut mat = eq.build_matrix();
    for row in &mat {
        println!("{:?}", row);
    }

    println!("\n----- GAUSSIAN ELIMINATION -----");

    let rank = lin_alg::row_reduce(&mut mat);
    for row in &mat {
        println!("{:?}", row);
    }

    println!("\n----- SOLVING NULL SPACE -----");

    let basis = lin_alg::solve(mat, rank);
    match basis.len() {
        0 => println!("(;w;) This equation has no solution!"),
        1 => {
            let sol = &basis[0];
            let ill_coef = eq.set_coefs(sol);

            // The `#` flag keeps preceding coefficients of 1.
            // FIXME: Let users choose whether to keep them or not.
            println!(
                "{}! Here's the solution:\n{:#}",
                if ill_coef { "(XwX) Ugh" } else { "(UwU) Yay" },
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
        n => {
            println!(
                "(OwO) Way too complex for Chemcat!\n\
                It's a linear combination of {} independent solutions.\n\
                Here's one possible basis:\n",
                n.green().bold()
            );
            for row in &basis {
                eq.set_coefs(row);
                println!("{:#}", eq);
            }
        }
    }

    println!();
}

fn eq_parser<'a>() -> impl Parser<'a, &'a str, ChemEq<'a>, extra::Err<Rich<'a, char>>> {
    let coef = text::int(10)
        .try_map(|s: &str, span| match s.parse::<i32>() {
            Ok(0) => Err(Rich::custom(span, "zero coefficient")),
            Ok(n) => Ok(n),
            Err(e) => Err(Rich::custom(span, e)),
        })
        .or(empty().to(1));

    let elem_name = any()
        .filter(char::is_ascii_uppercase)
        .then(any().filter(char::is_ascii_lowercase).repeated())
        .to_slice();

    let elem = elem_name.then(coef).map(|(name, n)| Term::Elem { name, n });

    let sub_term = recursive(|list_or_elem| {
        let list_content = list_or_elem.repeated().at_least(1).collect::<Vec<_>>();

        list_content
            .clone()
            .delimited_by(just('('), just(')'))
            .or(list_content.delimited_by(just('['), just(']')))
            .then(coef)
            .map(|(list, n)| Term::List { list, n, charge: 0 })
            .or(elem)
    })
    .repeated()
    .at_least(1)
    .collect::<Vec<_>>();

    let sub_term_separator = one_of("·.*");

    let charge = coef
        .then(just('+').to(false).or(just('-').to(true)))
        .delimited_by(just('('), just(')'))
        .map(|(n, neg)| if neg { -n } else { n })
        .or(empty().to(0));

    let term_without_coef = sub_term
        .clone()
        .then(
            sub_term_separator
                .ignore_then(coef.then(sub_term))
                .repeated()
                .collect::<Vec<_>>(),
        )
        .map(|(mut head, rest)| {
            head.extend(
                rest.into_iter()
                    .map(|(n, list)| Term::List { list, n, charge: 0 }),
            );
            head
        })
        .then(charge)
        .or(just("e(-)").map(|_| (vec![Term::Electron], -1)))
        .map(|(list, charge)| Term::List { list, n: 1, charge });

    let term = coef.ignore_then(term_without_coef);

    let side = term
        .separated_by(just('+').padded())
        .at_least(1)
        .collect::<Vec<_>>();

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
