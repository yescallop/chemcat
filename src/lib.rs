mod fmt;
pub mod lin_alg;

use std::collections::HashMap;

#[derive(Debug)]
pub struct ChemEq {
    pub terms: Vec<Term>,
    pub left_len: usize,
}

impl ChemEq {
    pub fn set_coefs(&mut self, coefs: &[i32]) -> bool {
        let mut ill_coef = false;
        self.terms
            .iter_mut()
            .zip(coefs.iter())
            .for_each(|(term, &coef)| match term {
                Term::List { n, .. } => {
                    *n = coef;
                    if coef <= 0 {
                        ill_coef = true;
                    }
                }
                _ => unreachable!(),
            });
        ill_coef
    }
}

#[derive(Debug)]
pub enum Term {
    Elem {
        name: String,
        n: i32,
    },
    Electron,
    List {
        list: Vec<Term>,
        n: i32,
        charge: i32,
    },
}

impl Term {
    const CHARGE_NAME: &'static str = "c";

    pub fn collapse<'a>(&'a self, map: &mut HashMap<&'a str, i32>, m: i32) {
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
            Term::Electron => {}
        }
    }
}

pub enum Solution {
    None,
    Unique(Vec<i32>),
    Infinite(Vec<Vec<i32>>),
}
