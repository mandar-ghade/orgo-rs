use std::{
    collections::{HashMap, VecDeque},
    fmt,
};

use itertools::Itertools;
use lazy_static::lazy_static;

use crate::matter::atom::{self, Atom};

#[derive(Debug)]
pub struct Configuration {
    pub vec: Vec<Subshell>,
    pub electrons: u32,
}

/// TODO: Work on exceptions

impl Configuration {
    /// Assuming neutral, unionized
    ///
    pub fn from_electrons(electrons: u32) -> Self {
        Self {
            vec: Vec::new(),
            electrons,
        }
    }

    pub fn from_atom(atom: &Atom) -> Self {
        Self {
            vec: Vec::new(),
            electrons: atom.electrons as u32,
        }
    }

    /// Gets subshell
    /// O(n) subshell search, where n is number of subshells
    fn get_subshell(
        &mut self,
        n: u32,
        orbital: Orbital,
    ) -> Option<&mut Subshell> {
        let l = *ORBITAL_TO_ANGULAR.get(&orbital)?;
        self.vec
            .iter_mut()
            .find(|subshell| subshell.n == n && subshell.l == l)
    }

    fn transfer(
        &mut self,
        electrons: u32,
        lhs: (u32, Orbital),
        rhs: (u32, Orbital),
    ) -> Option<()> {
        // TODO: Delete subshells with 0 electrons
        let (lhs_n, lhs_o) = lhs;
        let (rhs_n, rhs_o) = rhs;
        let four_s = self.get_subshell(lhs_n, lhs_o)?;
        four_s.remove_by(electrons); // can't use transfer func or the borrow checker will defeat me
        let three_d = self.get_subshell(rhs_n, rhs_o)?;
        three_d.fill_by(electrons);
        None
    }

    fn handle_exception(&mut self) -> Option<()> {
        match self.electrons {
            24 | 29 => self.transfer(1, (4, Orbital::S), (3, Orbital::D)),
            41 | 42 | 44 | 45 | 47 => {
                self.transfer(1, (5, Orbital::S), (4, Orbital::D))
            }
            46 => self.transfer(2, (5, Orbital::S), (4, Orbital::D)),
            78 | 79 => self.transfer(1, (6, Orbital::S), (5, Orbital::D)),
            _ => None,
        }
    }

    pub fn build(&mut self) -> Self {
        let mut queue = VecDeque::<Subshell>::new();
        let mut vec = Vec::<Subshell>::new();
        let mut amount_left = self.electrons;
        queue.push_back(Subshell::one_s());
        while !queue.is_empty() {
            let mut front = queue
                .pop_front()
                .expect("Front expected because non-empty queue");
            front.fill_by(amount_left);
            amount_left -= front.current;
            if amount_left == 0 {
                vec.push(front);
                break;
            }
            let l = front.l;
            let mut new_n = front.n;
            let last = front.get_next();
            vec.push(front);
            if l != 0 {
                for new_l in (0..l).rev() {
                    if amount_left == 0 {
                        break;
                    }
                    new_n += 1;
                    let mut subshell =
                        Subshell::generate_subshell(new_n, new_l)
                            .expect("Subshell expected.");
                    subshell.fill_by(amount_left);
                    amount_left -= subshell.current;
                    vec.push(subshell);
                }
            }
            if amount_left != 0 {
                queue.push_back(last.expect("Adjacent subshell expected"));
            }
        }
        self.vec = vec;
        self.handle_exception();
        Self {
            vec: self.vec.clone(),
            electrons: self.electrons,
        }
    }
}

impl fmt::Display for Configuration {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let str = self
            .vec
            .iter()
            .map(|s| format!("{}{}{}", s.n, s.l_name, s.current))
            .join(" ");
        write!(f, "{}", str)
    }
}

#[derive(Debug, Clone)]
pub struct Subshell {
    pub n: u32,
    pub l: u32,
    pub l_name: Orbital,
    pub capacity: u32,
    pub current: u32,
}

impl Subshell {
    pub fn one_s() -> Self {
        Self {
            n: 1,
            l: 0,
            l_name: Orbital::S,
            capacity: 2,
            current: 0,
        }
    }

    pub fn generate_subshell(n: u32, l: u32) -> Option<Self> {
        // F subshell is max
        let orbital = ANGULAR_TO_ORBITAL.get(&l)?;
        let capacity = *ELECTRON_MAXES.get(orbital)?;
        Some(Self {
            n,
            l,
            l_name: orbital.clone(),
            capacity,
            current: 0,
        })
    }

    /// Gets start of next diagonal
    pub fn get_next(&self) -> Option<Self> {
        let remaining = self.n - self.l - 1;
        if remaining == 0 {
            // We move to next principle quantum number once we've ran out of orbitals.
            // l is same because we are moving to next diagonal (don't want to recurse multiple
            // times)
            Self::generate_subshell(self.n + 1, self.l)
        } else if self.l == 3 {
            // we don't go past "f"
            None
        } else {
            // Just incrementing l
            Self::generate_subshell(self.n, self.l + 1)
        }
    }

    fn get_space(&self) -> u32 {
        self.capacity - self.current
    }

    fn fill_by(&mut self, amount: u32) {
        let space = self.get_space();
        if space == 0 {
            return;
        } else if space >= amount {
            self.current += amount;
        } else {
            self.current += space;
        }
    }

    fn remove_by(&mut self, amount: u32) {
        if amount > self.current {
            self.current = 0;
        } else {
            self.current -= amount;
        }
    }
}

#[derive(Hash, PartialEq, Eq, Clone, Debug)]
#[allow(dead_code)]
pub enum Orbital {
    S,
    P,
    D,
    F,
}

impl fmt::Display for Orbital {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::S => write!(f, "s"),
            Self::P => write!(f, "p"),
            Self::D => write!(f, "d"),
            Self::F => write!(f, "f"),
        }
    }
}

lazy_static! {
    pub static ref ANGULAR_TO_ORBITAL: HashMap<u32, Orbital> = HashMap::from([
        (0, Orbital::S),
        (1, Orbital::P),
        (2, Orbital::D),
        (3, Orbital::F),
    ]);
    pub static ref ORBITAL_TO_ANGULAR: HashMap<Orbital, u32> = HashMap::from([
        (Orbital::S, 0),
        (Orbital::P, 1),
        (Orbital::D, 2),
        (Orbital::F, 3),
    ]);
    pub static ref ELECTRON_MAXES: HashMap<Orbital, u32> = HashMap::from([
        (Orbital::S, 2),
        (Orbital::P, 6),
        (Orbital::D, 10),
        (Orbital::F, 14),
    ]);
}
