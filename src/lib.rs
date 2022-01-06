#![warn(missing_docs)]
//! Chemcat's brain (UwU).

/// Formatting chemical equations and terms.
mod fmt;
/// Bare bone linear algebra. In Gauss We Trust.
pub mod lin_alg;

use std::collections::HashMap;

/// A chemical equation.
#[derive(Debug)]
pub struct ChemEq {
    /// The immediate terms within the equation, starting from the left-hand side.
    pub terms: Vec<Term>,
    /// The number of terms on the left-hand side.
    pub left_len: usize,
}

impl ChemEq {
    /// Sets the coefficients of this `ChemEq` and returns a boolean
    /// indicating whether they are all positive.
    ///
    /// # Panics
    ///
    /// Panics if the numbers of terms and coefficients are not equal.
    pub fn set_coefs(&mut self, coefs: &[i32]) -> bool {
        assert!(self.terms.len() == coefs.len());

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

    /// Builds a matrix from this `ChemEq`.
    pub fn build_matrix(&self) -> Vec<Vec<i32>> {
        let term_cnt = self.terms.len();
        let mut coef_per_elem = HashMap::<&str, i32>::new();
        let mut eq_per_elem = HashMap::<&str, Vec<i32>>::new();

        for (i, term) in self.terms.iter().enumerate() {
            term.collapse(&mut coef_per_elem, 1);

            let left = i < self.left_len;
            for (elem, n) in coef_per_elem.drain() {
                let row = eq_per_elem.entry(elem).or_insert_with(|| vec![0; term_cnt]);
                row[i] = if left { n } else { -n };
            }
        }

        eq_per_elem.into_values().collect()
    }
}

/// A term within a chemical equation.
#[derive(Debug)]
pub enum Term {
    /// An element.
    Elem {
        /// The name.
        name: String,
        /// The coefficient.
        n: i32,
    },
    /// An electron.
    Electron,
    /// A list of terms.
    List {
        /// The list.
        list: Vec<Term>,
        /// The coefficient.
        n: i32,
        /// The electric charge.
        charge: i32,
    },
}

impl Term {
    /// The name of electric charge as an "element".
    pub const CHARGE_NAME: &'static str = "c";

    /// Collapses the term recursively into "elements" and coefficients
    /// with the given multiplier `m`, and adds the results into the map.
    ///
    /// Electric charge is regarded as an "element" named [`Term::CHARGE_NAME`]
    /// (currently "c", lowercase so that it does not conflict with actual elements).
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
