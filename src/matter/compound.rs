use strum_macros::Display;
use thiserror::Error;

use crate::{
    matter::atom::Atom,
    misc::element_info::{ELEMENTS, ELEMENT_WEIGHTS},
};
use std::{
    collections::HashMap,
    iter::Peekable,
    str::{Chars, FromStr},
};

use super::matter_trait::Matter;

#[derive(Clone, PartialEq, Debug)]
pub struct Compound<'a> {
    // Covalent Compound
    center: Atom,
    atom_map: HashMap<u8, Atom>,
    compound_map: HashMap<u8, &'a Compound<'a>>,
    order: [u8; 6],
    // at max 6 compounds - even an expanded octet can have at MOST 6 covalent bonds
}

impl Matter for Compound<'_> {
    fn to_string(&self) -> String {
        let mut s = String::new();
        for atom in self.atom_map.values() {
            s.push_str(atom.to_string().as_str());
        }
        s
    }
}

#[derive(Debug, Error)]
pub enum CompoundError {
    #[error("Invalid Element received: `{0}`")]
    InvalidElementError(String),
}

fn is_open_delimiter(c: &char) -> bool {
    *c == '(' || *c == '['
}

fn is_close_delimiter(c: &char) -> bool {
    *c == ')' || *c == ']'
}

fn is_delimiter(c: &char) -> bool {
    *c == '(' || *c == ')' || *c == '[' || *c == ']'
}

fn get_element_str(curr: &char, iter: &mut Peekable<Chars>) -> Option<String> {
    match Some(curr.to_string()).map(|mut c| {
        if let Some(nxt) = iter.peek() {
            c.push(*nxt);
            let _ = iter.next();
        }
        c
    }) {
        Some(c) => {
            if ELEMENTS.get(c.as_str()).is_some() {
                Some(c)
            } else {
                None
            }
        }
        None => None,
    }
}

impl<'a> Compound<'_> {
    fn parse(s: &str) -> Result<Self, CompoundError> {
        let center: Option<Atom>;
        let mut atom_map: HashMap<u8, Atom>;
        let mut compound_map: HashMap<u8, &'a Compound<'a>>;
        let mut order: [u8; 6];
        let mut idx: u8 = 0;
        let mut iter = s.chars().peekable();
        while let Some(curr) = iter.next() {
            if !curr.is_ascii_alphabetic() {
                continue;
            }
            // O(n)
            let Some(element_str) = get_element_str(&curr, &mut iter) else {
                continue;
            };
            dbg!(element_str);
        }
        todo!()
    }
}

impl FromStr for Compound<'_> {
    type Err = CompoundError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Compound::parse(s)
    }
}
