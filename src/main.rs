use chemcat::*;
use chumsky::prelude::*;
use owo_colors::OwoColorize;
use std::{io, thread};

fn main() {
    let thread = thread::current();

    ctrlc::set_handler(move || {
        thread.unpark();
    })
    .expect("Error setting Ctrl-C handler");

    thread::spawn(run);
    thread::park();
}

fn run() {
    let mut buf = String::new();
    let eq_parser = eq_parser();
    loop {
        println!("----- CHEMCAT (UwU) -----\nNya! Feed me an equation:");
        io::stdin()
            .read_line(&mut buf)
            .expect("failed to read line");

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
            .delimited_by(just('('), just(')'))
            .or(list_content.delimited_by(just('['), just(']')))
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
        .delimited_by(just('('), just(')'))
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
        .or(just("e(-)").map(|_| (vec![Term::Electron], -1)))
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
