pub mod builder;
pub mod deserializer;

use itertools::sorted;
use std::{
    collections::{BTreeSet, HashMap, LinkedList},
    fmt,
    iter::Peekable,
    str::{Chars, FromStr},
};

use crate::{
    constants::ELEMENTS,
    matter::{atom::Atom, compound::deserializer::Chain},
};

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
pub struct Location {
    pub x: i16,
    pub y: i16,
}

#[allow(dead_code)]
impl Location {
    pub fn new(x: i16, y: i16) -> Self {
        Self { x, y }
    }

    pub fn shift(&self, dx: i16, dy: i16) -> Self {
        Self {
            x: self.x + dx,
            y: self.y + dy,
        }
    }
}

pub struct SimpleCompound {}

#[derive(Clone, PartialEq, Debug)]
/// Represents a Linear Compound
pub struct Compound {
    /// Covalent Compound
    atoms: Vec<Atom>,
    /// Each index represents the location of an atom
    locations: Vec<Location>,
    /// From location, we can compute the Atom's index
    location_to_idx: HashMap<Location, usize>,
    /// should theoretically be a size
    backbone: Vec<usize>,
    /// retains order
    /// ^ Side chains can have side-chains (unfortunately)
    /// ^ sort of like an undirected graph
    side_chains: HashMap<usize, BTreeSet<usize>>,
    // TODO: Ensure values != key or backbone idx
}
// TODO: Pseudo-Dijkstra's longest chain implementation (using largest distance)

impl Compound {
    pub fn new(
        atoms: Vec<Atom>,
        locations: Vec<Location>,
        location_to_idx: HashMap<Location, usize>,
        backbone: Vec<usize>,
        side_chains: HashMap<usize, BTreeSet<usize>>,
    ) -> Self {
        Self {
            atoms,
            locations,
            location_to_idx,
            backbone,
            side_chains,
        }
    }

    fn get_atom(&self, i: usize) -> Option<&Atom> {
        self.atoms.get(i)
    }

    fn get_atom_unsafe(&self, idx: usize) -> &Atom {
        if let Some(atom) = self.atoms.get(idx) {
            atom
        } else {
            panic!(
                "Invalid atom index: {} (length is {})",
                idx,
                self.atoms.len(),
            );
        }
    }

    fn has_side_chain(&self, i: usize) -> bool {
        if let Some(side_chain) = self.side_chains.get(&i) {
            !side_chain.is_empty()
        } else {
            false
        }
    }

    fn get_sidechain_unsafe(&self, i: usize) -> &BTreeSet<usize> {
        self.side_chains.get(&i).expect("Side chain not found")
    }
}

impl fmt::Display for Compound {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&Chain::from(self).to_string())
    }
}

pub type CompoundResult<T> = Result<T, CompoundError>;

#[allow(dead_code)]
#[derive(thiserror::Error, Debug)]
pub enum CompoundError {
    #[error("Compound Parsing Error: {0}")]
    Parsing(String),
    #[error("Unknown Error: {0}")]
    Unknown(String),
}

#[inline]
fn is_delimiter(c: char) -> bool {
    is_open_delimiter(c) || is_close_delimiter(c)
}

#[inline]
fn is_open_delimiter(c: char) -> bool {
    c == '(' || c == '['
}

#[inline]
fn is_close_delimiter(c: char) -> bool {
    c == ')' || c == ']'
}

#[allow(dead_code)]
struct Token {
    str: String,
    pub count: u16,
    pub parts: Vec<Token>,
}

#[allow(dead_code)]
impl Token {
    pub fn new(s: String) -> Self {
        Token {
            str: s,
            count: 0,
            parts: Vec::new(),
        }
    }
}

fn get_element_str(curr: char, it: &mut Peekable<Chars>) -> Option<String> {
    let sym = match it.next() {
        Some(c) => format!("{}{}", curr, c),
        None => curr.to_string(),
    };
    if ELEMENTS.contains(sym.as_str()) {
        Some(sym)
    } else {
        None
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
    count_str.parse::<u16>().unwrap_or(1)
}

#[allow(dead_code)]
fn parse_into_tokens(
    iter: &mut Peekable<Chars>,
) -> Result<Vec<Token>, CompoundError> {
    let tokens_vec: Vec<Token> = Vec::new();
    // let mut carbon_tk: Option<Token> = None;
    let mut curr_tk: Option<Token> = None;
    while let Some(curr) = iter.next() {
        if is_open_delimiter(curr) && curr_tk.is_some() {
            let parts = parse_into_tokens(iter)?;
            // let count = get_element_count(iter);
            if let Some(token) = curr_tk.as_mut() {
                token.parts = parts;
                // token.count = count;
            }
        }
        if is_close_delimiter(curr) {
            return Ok(tokens_vec);
        }
        let Some(element_str) = get_element_str(curr, iter) else {
            continue;
        };
        let _count = get_element_count(iter);
        dbg!(element_str);
    }
    Ok(tokens_vec)
}

impl Compound {
    pub fn parse(_s: &str) -> CompoundResult<Self> {
        todo!()
    }
}

impl FromStr for Compound {
    type Err = CompoundError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Compound::parse(s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compound_to_string_simple() {
        let comp = Compound {
            atoms: vec![
                Atom::new_unchecked(1),
                Atom::new_unchecked(2),
                Atom::new_unchecked(2),
                Atom::new_unchecked(3),
            ],
            locations: vec![
                Location::new(1, 2),
                Location::new(2, 3),
                Location::new(2, 2),
                Location::new(1, 3),
            ],
            location_to_idx: HashMap::from([
                (Location::new(1, 2), 0),
                (Location::new(2, 3), 1),
                (Location::new(2, 2), 2),
                (Location::new(1, 3), 3),
            ]),
            backbone: vec![0, 2],
            side_chains: HashMap::from([(0, BTreeSet::from([1, 3]))]),
        };
        assert_eq!(comp.to_string(), "HHeLiHe");
    }
}
