use std::fmt;

#[derive(Clone, Debug)]
pub enum Chain {
    /// COMPOUND
    Vec(Vec<Chain>, usize),
    /// ATOMS
    KV(String, usize),
}

#[allow(dead_code)]
impl Chain {
    fn incr_count_by(&mut self, count: usize) {
        match self {
            Self::Vec(_, c) => *c += count,
            Self::KV(_, c) => *c += count,
        }
    }

    fn get_count(&self) -> usize {
        match self {
            Self::Vec(_, c) => *c,
            Self::KV(_, c) => *c,
        }
    }

    pub fn group(&self) -> Self {
        match self {
            Self::Vec(chains, count) => {
                let mut new_vec: Vec<Chain> = Vec::new();
                for curr_chain in chains.iter() {
                    let curr = curr_chain.group().clone(); // Recursion yay
                    if let Some(matched) =
                        new_vec.iter_mut().find(|c| c.custom_eq(&curr))
                    {
                        matched.incr_count_by(curr.get_count());
                    } else {
                        new_vec.push(curr);
                    }
                }
                Self::Vec(new_vec, *count)
            }
            Self::KV(_, _) => self.clone(),
        }
    }

    pub fn minimize(&self, factor: usize) -> Self {
        // Convert redundant single-sized vectors into KV (ATOM).
        //
        // Factor = 1 should be default
        match self {
            Self::KV(k, v) => Self::KV(k.clone(), *v * factor),
            Self::Vec(v, c) => {
                if v.len() != 1 {
                    self.clone()
                } else {
                    v.first()
                        .expect("First should've been found")
                        .minimize(*c * factor)
                }
            }
        }
    }

    fn str(&self) -> String {
        let minimized = self.group().minimize(1);
        match minimized {
            Self::KV(k, v) => {
                if v == 1 {
                    k
                } else {
                    format!("{}{}", k, v)
                }
            }
            Self::Vec(v, count) => {
                match (count, v.len()) {
                    // Vec len != 1 due to minimization
                    (_, 0) => {
                        panic!("Vec cannot have length 0");
                    }
                    (1, 1) => {
                        panic!("Grouping & minimalization aren't working correctly");
                    }
                    (_, 1) => {
                        panic!("Minimalization isn't working correctly");
                    }
                    (_, _) => {
                        let mut output_str = String::new();
                        for i in v.iter() {
                            output_str.push_str(&i.str())
                        }
                        if count == 1 {
                            output_str
                        } else {
                            format!("({}){}", output_str, count)
                        }
                    }
                }
            }
        }
    }

    fn custom_eq(&self, other: &Self) -> bool {
        // IF vec, check everything is same
        // If KV, check that identities are the same (Atoms)
        //
        // We have two implementations of equality, a custom equality and a full equality
        match (self, other) {
            (Self::KV(k_lhs, _), Self::KV(k_rhs, _)) => {
                // Counts DO NOT matter (non-vec comparison)
                k_lhs == k_rhs
            }
            (Self::Vec(_, _), Self::Vec(_, _)) => self == other, // EXACT checks
            (_, _) => false,
        }
    }
}

impl PartialEq for Chain {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::KV(k_lhs, v_lhs), Self::KV(k_rhs, v_rhs)) => {
                // Counts matter for side chains to be equal (vec comparison).
                k_lhs == k_rhs && v_lhs == v_rhs
            }
            (Self::Vec(lhs, _), Self::Vec(rhs, _)) => {
                // Exact side chain counts won't matter, their components do.
                lhs.len() == rhs.len()
                    && lhs.iter().all(|l| rhs.contains(l))
                    && rhs.iter().all(|r| lhs.contains(r))
            }
            (_, _) => false,
        }
    }
}

impl fmt::Display for Chain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_grouping_and_minimize_test() {
        // ((H5H5)(H5H5)) => H20
        let inner = Chain::Vec(
            Vec::from([Chain::KV("H".into(), 5), Chain::KV("H".into(), 5)]),
            1,
        );
        let outer = Vec::from([inner.clone(), inner]);
        let lhs = Chain::Vec(outer, 1).group().minimize(1);
        let rhs = Chain::KV("H".into(), 20);
        assert_eq!(
            lhs, rhs,
            "Grouping & minimizing aren't functioning correctly."
        )
    }

    #[test]
    fn condensed_formula_grouping_test() {
        // (H20I2I5H33)2 => (H53I7)2
        let chain = Chain::Vec(
            Vec::from([
                Chain::KV("H".into(), 20),
                Chain::KV("I".into(), 2),
                Chain::KV("I".into(), 5),
                Chain::KV("H".into(), 33),
            ]),
            2,
        );
        assert_eq!(
            chain.to_string(),
            "(H53I7)2".to_string(),
            "Formula is not properly condensed."
        );
    }
}
