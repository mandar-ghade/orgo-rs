use itertools::sorted;
use std::{
    collections::{BTreeSet, HashMap, LinkedList},
    fmt,
    iter::Peekable,
    str::{Chars, FromStr},
};

use crate::{constants::ELEMENTS, matter::atom::Atom};

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

fn extract_cmp_str_and_count(
    inp: &String,
) -> CompoundResult<(String, Option<u8>)> {
    let has_delims = inp.chars().any(|c| is_delimiter(c));
    let has_count = inp.chars().any(|c| c.is_digit(10));
    match (has_delims, has_count) {
        (false, false) => Ok((inp.to_owned(), None)),
        (false, true) => {
            // Count but no delimiters
            let mut fmt_str = String::new();
            let mut count = String::new();
            for c in inp.chars() {
                if c.is_digit(10) {
                    count.push(c);
                } else if count.is_empty() {
                    fmt_str.push(c);
                } else {
                    return Err(CompoundError::Parsing(
                        "XS string detected (extracting cmp str and count)"
                            .into(),
                    ));
                }
            }
            let count: u8 = count
                .parse()
                .map_err(|_| CompoundError::Parsing("Invalid count.".into()))?;
            Ok((fmt_str, Some(count)))
        }
        (true, false) => {
            // Delimiters but no count
            let size = inp.len();
            if size <= 2 {
                return Err(CompoundError::Parsing(
                    "Misformatted compound str".into(),
                ));
            }
            let mut fmt_str = String::new();
            for (i, c) in inp.chars().enumerate() {
                if i == 0 || i == size - 1 {
                    assert!(is_delimiter(c), "Issue mapping out delimiter locations for compound with no count listing.");
                    continue;
                }
                fmt_str.push(c);
            }
            if fmt_str.is_empty() {
                return Err(CompoundError::Parsing(
                    "No string detected within delimiters.".into(),
                ));
            }
            Ok((fmt_str, None))
        }
        (true, true) => {
            // Both delimiters and a count
            let size = inp.len();
            if size <= 2 {
                return Err(CompoundError::Parsing(
                    "Misformatted compound str.".into(),
                ));
            }
            let mut fmt_str = String::new();
            let mut count = String::new();
            let mut lhs_delim_count = 0;
            let mut rhs_delim_count = 0;
            for c in inp.chars() {
                // TODO: Check if open & close delims match: () vs []
                if is_open_delimiter(c) {
                    lhs_delim_count += 1;
                } else if is_close_delimiter(c) {
                    rhs_delim_count += 1;
                } else if c.is_digit(10) && lhs_delim_count == rhs_delim_count {
                    // TODO: is this right? Shouldn't it be recursive?
                    // Matching # of opening & closing delims
                    count.push(c);
                } else if count.is_empty() {
                    fmt_str.push(c);
                } else {
                    return Err(CompoundError::Parsing(
                        "XS string detected (extracting cmp str and count)"
                            .into(),
                    ));
                }
            }
            if fmt_str.is_empty() {
                return Err(CompoundError::Parsing(
                    "No string detected within delimiters.".into(),
                ));
            }
            assert!(
                !count.is_empty(),
                "Count not found when searching around delims (nesting issues)"
            );
            let count_u8: u8 = count
                .parse()
                .map_err(|_| CompoundError::Parsing("Invalid count.".into()))?;
            Ok((fmt_str, Some(count_u8)))
        }
    }
}

impl Compound {
    fn get_atom(&self, i: u8) -> Option<&Atom> {
        self.atoms.get(i as usize)
    }

    fn has_side_chain(&self, i: u8) -> bool {
        if let Some(side_chain) = self.side_chains.get(&i) {
            side_chain.is_empty()
        } else {
            false
        }
    }

    fn get_sidechain_len(&self, backbone_i: u8) -> CompoundResult<u8> {
        let side_chain =
            self.side_chains.get(&backbone_i).ok_or_else(|| {
                CompoundError::Unknown("Side chain not found.".into())
            })?;
        side_chain.len().try_into().map_err(|_| {
            CompoundError::Unknown("Invalid side chain length.".into())
        })
    }

    fn side_chain_as_str(&self, backbone_i: u8) -> Result<String, fmt::Error> {
        let mut s = String::new();
        self.write_side_chain(&mut s, backbone_i)?;
        Ok(s)
    }

    fn write_side_chain<W: fmt::Write>(
        &self,
        w: &mut W,
        backbone_i: u8,
    ) -> fmt::Result {
        let side_chain = self
            .side_chains
            .get(&backbone_i)
            .expect("Side chain not found");
        let mut q: LinkedList<String> = LinkedList::new();
        // queue
        for &i in sorted(side_chain) {
            let atom_str = self
                .get_atom(i)
                .expect("Atom expected at side chain index")
                .to_string();
            if !self.has_side_chain(i) && !q.is_empty() {
                // push everything before to main string
                // We have different atom, so we still need to append side chain for previous atom.
                let count = q.len();
                let q_str = q.pop_back().unwrap();

                let one_letter = q_str.len() == 1;
                if count == 1 && one_letter {
                    write!(w, "{}", q_str)?;
                } else if count == 1 {
                    write!(w, "({})", q_str)?;
                } else if one_letter {
                    write!(w, "{}{}", q_str, count)?;
                } else {
                    write!(w, "({}){}", q_str, count)?;
                }
                q = LinkedList::new(); // resets queue
                write!(w, "{}", atom_str)?;
                continue;
            } else if !self.has_side_chain(i) {
                write!(w, "{}", atom_str)?;
                continue;
            }
            let sc_str = self.side_chain_as_str(i)?;

            let Some(q_front) = q.front() else {
                q.push_back(atom_str);
                continue;
            };
            // q isn't empty, 4th case
            if atom_str == *q_front {
                q.push_back(atom_str);
                continue;
            }

            // most important case, we have different atom with a side chain, while we also need to
            // append the remainder of the previous back to the end.

            // Append remainder of queue to string & flush queue
            let count = q.len();
            let q_str = q.pop_back().expect(
                "Couldn't pop back of queue when assessing its length.",
            );
            let one_letter = q_str.len() == 1;
            if count == 1 && one_letter {
                write!(w, "{}", q_str)?;
            } else if count == 1 {
                write!(w, "({})", q_str)?;
            } else if one_letter {
                write!(w, "{}{}", q_str, count)?;
            } else {
                write!(w, "({}){}", q_str, count)?;
            }
            q = LinkedList::new(); // reset queue

            // Current side chain str
            let (cmp_str, count) = extract_cmp_str_and_count(&sc_str)
                .expect("extract_cmp_str_and_count failed.");
            // Pushes current side chain to queue
            if let Some(n) = count {
                for _ in 0..n {
                    q.push_back(cmp_str.clone()); // lots of memory
                }
            } else {
                q.push_back(cmp_str);
            }
        }
        if let Some(q_str) = q.pop_back() {
            let count = q.len() + 1;
            let one_letter = q_str.len() == 1;
            if count == 1 && one_letter {
                write!(w, "{}", q_str)?;
            } else if count == 1 {
                write!(w, "({})", q_str)?;
            } else if one_letter {
                write!(w, "{}{}", q_str, count)?;
            } else {
                write!(w, "({}){}", q_str, count)?;
            }
        }
        Ok(())
    }
}

impl fmt::Display for Compound {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in self.backbone.iter() {
            let atom = self
                .atoms
                .get(*i as usize)
                .expect("Atom deleted from indexed");
            write!(f, "{}", atom)?;
            if self.has_side_chain(*i) {
                let sc_length = self.get_sidechain_len(*i).unwrap();
                assert!(
                    sc_length != 0,
                    "Side chain length cannot be zero if it has a side chain."
                );
                let sc_str = self.side_chain_as_str(*i)?;
                if sc_length != 1 {
                    write!(f, "({})", sc_str)?;
                } else {
                    write!(f, "{}", sc_str)?;
                }
            }
        }
        Ok(())
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
    fn parse(_s: &str) -> Result<Self, CompoundError> {
        todo!()
    }
}

impl FromStr for Compound {
    type Err = CompoundError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Compound::parse(s)
    }
}
