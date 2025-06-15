use itertools::sorted;
use thiserror::Error;

use crate::{constants::ELEMENTS, matter::atom::Atom};
use core::panic;
use std::{
    collections::{BTreeSet, HashMap, LinkedList},
    fmt,
    iter::Peekable,
    str::{Chars, FromStr},
};

use super::matter::Matter;

#[derive(Clone, PartialEq, Debug, Hash, Eq)]
struct Location {
    x: u16,
    y: u16,
}

#[derive(Clone, PartialEq, Debug)]
/// Represents a Linear Compound
pub struct Compound {
    /// Covalent Compound
    atoms: Vec<Atom>,
    /// Each index represents the location of an atom
    locations: Vec<Location>,
    /// From location, we can compute the Atom's index
    location_to_idx: HashMap<Location, u8>,
    /// should theoretically be a size
    backbone: Vec<u8>,
    /// retains order
    /// ^ Side chains can have side-chains (unfortunately)
    /// ^ sort of like an undirected graph
    side_chains: HashMap<u8, BTreeSet<u8>>,
    // TODO: Ensure values != key or backbone idx
}
// TODO: Pseudo-Dijkstra's longest chain implementation (using largest distance)

impl Compound {
    fn get_atom(&self, i: u8) -> Option<&Atom> {
        self.atoms.get(i as usize)
    }

    fn has_side_chain(&self, i: u8) -> bool {
        self.side_chains.contains_key(&i)
            && !self
                .side_chains
                .get(&i)
                .expect(
                    "Expected side chain set while evaluating its emptiness",
                )
                .is_empty()
    }

    fn get_sidechain_len(&self, backbone_i: u8) -> u8 {
        assert!(
            self.has_side_chain(backbone_i),
            "Side chain was expected while retrieving backbone's side chain length."
        );
        self.side_chains
            .get(&backbone_i)
            .expect("Expected side chain set while evaluating its emptiness")
            .len()
            .try_into()
            .expect("Couldn't convert u8 to usize; Getting side chain length")
    }

    fn side_chain_as_str(&self, backbone_i: u8) -> String {
        fn extract_cmp_str_and_count(inp: &String) -> (String, Option<u8>) {
            let has_delims = inp
                .chars()
                .any(|c| is_open_delimiter(&c) || is_close_delimiter(&c));
            let has_count = inp.chars().any(|c| c.is_digit(10));
            if !has_delims && !has_count {
                (inp.to_owned(), None)
            } else if !has_delims && has_count {
                let mut fmt_str = String::new();
                let mut count: String = String::new();
                for c in inp.chars() {
                    if !c.is_digit(10) && count.is_empty() {
                        fmt_str.push(c);
                    } else if c.is_digit(10) {
                        count.push(c);
                    } else if !count.is_empty() {
                        panic!(
                            "XS string detected (extracting cmp str and count)"
                        );
                    }
                }
                let count_u8: u8 =
                    count.parse().expect("Couldn't convert to u8");
                (fmt_str, Some(count_u8))
            } else if has_delims && !has_count {
                let size = inp.len();
                if size <= 2 {
                    panic!("Cannot handle misformatted compound str");
                }
                let mut fmt_str = String::new();
                for (i, c) in inp.chars().enumerate() {
                    if i == 0 || i == size - 1 {
                        if !(is_open_delimiter(&c) || is_close_delimiter(&c)) {
                            panic!(
                                "Issue mapping out delimiter locations for compound with no count listing"
                            );
                        }
                        continue;
                    }
                    fmt_str.push(c);
                }
                if fmt_str.is_empty() {
                    panic!("No string detected within delimiters");
                }
                (fmt_str, None)
            } else if has_delims && has_count {
                let size = inp.len();
                if size <= 2 {
                    panic!("Cannot handle misformatted compound str");
                }
                let mut fmt_str = String::new();
                let mut count: String = String::new();
                let mut lhs_delim_count = 0;
                let mut rhs_delim_count = 0;
                for c in inp.chars() {
                    // TODO: Check if open & close delims are the same type (brackets vs parens)
                    if is_open_delimiter(&c) {
                        lhs_delim_count += 1;
                    } else if is_close_delimiter(&c) {
                        rhs_delim_count += 1;
                    } else if !c.is_digit(10) && count.is_empty() {
                        fmt_str.push(c);
                    } else if c.is_digit(10)
                        && lhs_delim_count == rhs_delim_count
                    {
                        // Matching # of opening & closing delims
                        count.push(c);
                    } else if !count.is_empty() {
                        panic!("XS string detected (extracting cmp str with delims and count)");
                    }
                }
                if fmt_str.is_empty() {
                    panic!("No string detected within delimiters");
                }
                if count.is_empty() {
                    panic!("Count not found when searching around delims (nesting issues)")
                }
                let count_u8: u8 =
                    count.parse().expect("Couldn't convert to u8");
                (fmt_str, Some(count_u8))
            } else {
                panic!("Undefined behavior"); // Should never execute
            }
        }

        assert!(self.has_side_chain(backbone_i), "Side chain was expected");
        let mut s: String = String::new();
        let mut q: LinkedList<String> = LinkedList::new();
        // queue
        for (atom_str, i) in sorted(
            self.side_chains
                .get(&backbone_i)
                .expect("Side chain was expected")
                .iter()
                .map(|&i| {
                    (
                        self.get_atom(i)
                            .expect("Atom expected at side chain index")
                            .to_string(),
                        i,
                    )
                }),
        ) {
            if !self.has_side_chain(i) && !q.is_empty() {
                // push everything before to main string
                // We have different atom, so we still need to append side chain for previous atom.
                let count = q.len();
                let q_str = q
                    .pop_back()
                    .expect("Error: Couldn't pop back of queue when assessing its length");

                let one_letter = q_str.len() == 1;
                if count == 1 && one_letter {
                    s.push_str(q_str.as_str());
                } else if count == 1 {
                    s.push_str(format!("({})", q_str).as_str());
                } else if one_letter {
                    s.push_str(format!("{}{}", q_str, count).as_str());
                } else {
                    s.push_str(format!("({}){}", q_str, count).as_str());
                }
                q = LinkedList::new(); // resets queue
                s.push_str(atom_str.as_str());
                continue;
            } else if !self.has_side_chain(i) {
                s.push_str(atom_str.as_str());
                continue;
            }
            let sc_str = self.side_chain_as_str(i);
            if q.is_empty() {
                q.push_back(atom_str);
                continue;
            }

            // q isn't empty, 4th case
            //

            if atom_str
                == *q
                    .front()
                    .expect("Couldn't get front of LinkedList even though it wasn't empty")
            {
                q.push_back(atom_str);
                continue;
            }

            // most important case, we have different atom with a side chain, while we also need to
            // append the remainder of the previous back to the end.

            // Append remainder of queue to string & flush queue
            {
                let count = q.len();
                let q_str = q
                    .pop_back()
                    .expect("Error: Couldn't pop back of queue when assessing its length");

                let one_letter = q_str.len() == 1;
                if count == 1 && one_letter {
                    s.push_str(q_str.as_str());
                } else if count == 1 {
                    s.push_str(format!("({})", q_str).as_str());
                } else if one_letter {
                    s.push_str(format!("{}{}", q_str, count).as_str());
                } else {
                    s.push_str(format!("({}){}", q_str, count).as_str());
                }
                q = LinkedList::new(); // resets queue
            }

            // Current side chain str
            let (cmp_str, count) = extract_cmp_str_and_count(&sc_str);
            // Pushes current side chain to queue
            if let Some(n) = count {
                for _ in 0..n {
                    q.push_back(cmp_str.clone()); // lots of memory
                }
            } else {
                q.push_back(cmp_str);
            }
        }
        if !q.is_empty() {
            let count = q.len();
            let q_str = q.pop_back().expect(
                "Error: Couldn't pop back of queue when assessing its length",
            );

            let one_letter = q_str.len() == 1;
            if count == 1 && one_letter {
                s.push_str(q_str.as_str());
            } else if count == 1 {
                s.push_str(format!("({})", q_str).as_str());
            } else if one_letter {
                s.push_str(format!("{}{}", q_str, count).as_str());
            } else {
                s.push_str(format!("({}){}", q_str, count).as_str());
            }
        }
        s
    }

    pub fn as_str(&self) -> String {
        fn wrap_str(s: &mut String, side_chain_str: &str) {
            s.push('(');
            s.push_str(side_chain_str);
            s.push(')')
        }

        let mut s: String = String::new();
        for i in self.backbone.iter() {
            let atom = self
                .atoms
                .get(*i as usize)
                .expect("Atom deleted from indexed");
            s.push_str(atom.to_string().as_str());
            if self.has_side_chain(*i) {
                let sc_length = self.get_sidechain_len(*i);
                assert!(
                    sc_length != 0,
                    "Side chain length cannot be zero if it has a side chain....."
                );
                let sc_str = self.side_chain_as_str(*i);
                if sc_length != 1 {
                    wrap_str(&mut s, &sc_str);
                } else {
                    s.push_str(&sc_str);
                }
            }
        }
        s
    }
}

impl fmt::Display for Compound {
    // Converts Compound to string
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

impl Matter for Compound {
    fn to_string(&self) -> String {
        self.as_str()
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

fn parse_into_tokens(
    iter: &mut Peekable<Chars>,
) -> Result<Vec<Token>, CompoundError> {
    let tokens_vec: Vec<Token> = Vec::new();
    // let mut carbon_tk: Option<Token> = None;
    let mut curr_tk: Option<Token> = None;
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
    Ok(tokens_vec)
}

impl Compound {
    fn parse(s: &str) -> Result<Self, CompoundError> {
        todo!()
    }
}

impl FromStr for Compound {
    type Err = CompoundError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        todo!()
        // Compound::parse(s)
    }
}
