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

struct Token {
    str: String,
    pub count: u16,
    pub parts: Vec<Token>,
}

impl Token {
    pub fn new(s: String) -> Self {
        Token {
            str: s,
            count: 0,
            parts: Vec::new(),
        }
    }
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

fn get_element_count(iter: &mut Peekable<Chars>) -> u16 {
    let mut count_str = String::new();
    while let Some(next) = iter.peek() {
        if !next.is_numeric() {
            break;
        }
        if !next.is_whitespace() {
            count_str.push(*next);
        }
        iter.next();
    }
    count_str.parse().unwrap_or(1)
}

fn parse_into_tokens(iter: &mut Peekable<Chars>) -> Result<Vec<Token>, CompoundError> {
    let mut tokens_vec: Vec<Token> = Vec::new();
    let mut carbon_tk: Option<Token> = None;
    let mut curr_tk: Option<Token> = None;
    // TODO: Get mapping for number of max substituents
    // Ex: C: 4 (unless carbanion or carbocation)
    // S: 6 (expanded octet)
    // What about if the substituent is first?
    // Ex: H3C, etc?
    // I need to rethink the way I do parsing.
    // Perhaps I can just search for carbon compounds first then work my way to doing more complex parsing.
    // H2CH2
    // H2O
    // H2O2
    // H2SO4
    // Design an algorithm which can detect this.
    while let Some(curr) = iter.next() {
        if is_open_delimiter(&curr) && curr_tk.is_some() {
            let parts = parse_into_tokens(iter)?;
            //let count = get_element_count(iter);
            if let Some(token) = curr_tk.as_mut() {
                token.parts = parts;
                //token.count = count;
            }
        }
        if is_close_delimiter(&curr) {
            return Ok(tokens_vec);
        }
        let Some(element_str) = get_element_str(&curr, iter) else {
            continue;
        };
        let count = get_element_count(iter);
        dbg!(element_str);
    }
    return Ok(tokens_vec);
}

impl<'a> Compound<'_> {
    fn parse(s: &str) -> Result<Self, CompoundError> {
        let center: Option<Atom>;
        let mut atom_map: HashMap<u8, Atom>;
        let mut compound_map: HashMap<u8, &'a Compound<'a>>;
        let mut order: [u8; 6];
        let mut idx: u8 = 0;
        let tokens = parse_into_tokens(&mut s.chars().peekable())?;
        todo!()
    }
}

impl FromStr for Compound<'_> {
    type Err = CompoundError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Compound::parse(s)
    }
}
