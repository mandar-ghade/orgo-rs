use std::fmt::{self, Write};

#[derive(Clone, Debug)]
pub enum Chain {
    /// Compound
    Vec(Vec<Chain>, usize),
    /// Atoms
    KV(String, usize),
}

impl Chain {
    fn reversed(&self) -> Chain {
        // Reverses order of chain
        match self.clone().group().minimize(1) {
            Self::KV(k, v) => Self::KV(k, v),
            Self::Vec(v, c) => {
                let mut reversed_vec = v.clone();
                reversed_vec.reverse();
                reversed_vec =
                    reversed_vec.iter().map(|i| i.reversed()).collect();
                Self::Vec(reversed_vec, c).group().minimize(1)
            }
        }
    }

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

    pub fn group(self) -> Self {
        match self {
            Self::Vec(chains, count) => {
                let mut new_chains = Vec::<Chain>::new();
                for chain in chains {
                    let curr = chain.clone().group();
                    // Vec ordering implies connectivity.
                    // Connects side chains together
                    if let Some(matched) = new_chains.last_mut() {
                        if matched.custom_eq(&curr) {
                            matched.incr_count_by(curr.get_count());
                        } else {
                            new_chains.push(curr);
                        }
                    } else {
                        new_chains.push(curr);
                    }
                }
                Self::Vec(new_chains, count)
            }
            Self::KV(_, _) => self,
        }
    }

    pub fn minimize(self, factor: usize) -> Self {
        // Convert redundant single-sized vectors into KV (ATOM).
        //
        // Factor = 1 should be default
        match self {
            Self::KV(s, c) => Self::KV(s, c * factor),
            Self::Vec(ref v, c) => {
                if v.len() != 1 {
                    self
                } else {
                    v.first()
                        .expect("First should've been found")
                        .clone()
                        .minimize(c * factor)
                }
            }
        }
    }

    fn write_to<W: Write>(&self, w: &mut W) -> fmt::Result {
        match self.clone().group().minimize(1) {
            Self::KV(s, c) => {
                if c == 1 {
                    write!(w, "{}", s)?;
                } else {
                    write!(w, "{}{}", s, c)?;
                }
            }
            Self::Vec(v, count) => {
                assert!(!v.is_empty());
                if v.len() == 1 {
                    if count == 1 {
                        panic!("Grouping & minimalization aren't working correctly");
                    } else {
                        panic!("Minimalization isn't working correctly");
                    }
                }
                if count != 1 {
                    write!(w, "(")?;
                }
                for chain in v {
                    chain.write_to(w)?;
                }
                if count != 1 {
                    write!(w, "){}", count)?;
                }
            }
        }
        Ok(())
    }

    fn custom_eq(&self, other: &Self) -> bool {
        // IF vec, check everything is same
        // If KV, check that identities are the same (Atoms)
        // We have two implementations of equality, a custom equality and a full equality
        match (self, other) {
            (Self::KV(s_left, _), Self::KV(s_right, _)) => {
                // Counts DO NOT matter (non-vec comparison)
                s_left == s_right
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
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.write_to(f)
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
        // (H20H33I2I5)2 => (H53I7)2
        let chain = Chain::Vec(
            Vec::from([
                Chain::KV("H".into(), 20),
                Chain::KV("H".into(), 33),
                Chain::KV("I".into(), 2),
                Chain::KV("I".into(), 5),
            ]),
            2,
        );
        assert_eq!(
            chain.to_string(),
            "(H53I7)2",
            "Formula is not properly condensed."
        );
    }

    #[test]
    fn condense_tert_butanol() {
        // HOC((CH3)(CH3)(CH3)) => HOC(CH3)3
        let chain = Chain::Vec(
            Vec::from([
                Chain::Vec(
                    Vec::from([
                        Chain::KV("H".into(), 1),
                        Chain::KV("O".into(), 1),
                    ]),
                    1,
                ),
                Chain::KV("C".into(), 1),
                Chain::Vec(
                    Vec::from([
                        Chain::Vec(
                            Vec::from([
                                Chain::KV("C".into(), 1),
                                Chain::KV("H".into(), 3),
                            ]),
                            1,
                        ),
                        Chain::Vec(
                            Vec::from([
                                Chain::KV("C".into(), 1),
                                Chain::KV("H".into(), 3),
                            ]),
                            1,
                        ),
                        Chain::Vec(
                            Vec::from([
                                Chain::KV("C".into(), 1),
                                Chain::KV("H".into(), 3),
                            ]),
                            1,
                        ),
                    ]),
                    1,
                ),
            ]),
            1,
        );
        assert_eq!(
            chain.to_string(),
            "HOC(CH3)3",
            "Tert-butanol doesn't condense correctly"
        );
    }

    #[test]
    fn condense_butane() {
        // CH3(CH2)(CH2)CH3 => CH3(CH2)2CH3
        let terminal_cs = Chain::Vec(
            Vec::from([Chain::KV("C".into(), 1), Chain::KV("H".into(), 3)]),
            1,
        );
        let middle_cs = Chain::Vec(
            Vec::from([Chain::KV("C".into(), 1), Chain::KV("H".into(), 2)]),
            1,
        );
        let chain = Chain::Vec(
            Vec::from([
                terminal_cs.clone(),
                middle_cs.clone(),
                middle_cs,
                terminal_cs,
            ]),
            1,
        );
        assert_eq!(
            chain.to_string(),
            "CH3(CH2)2CH3",
            "Failed to condense butane"
        )
    }

    #[test]
    fn condense_1_4_butanediol() {
        // TODO: Make a function that creates a linear chain for Chain
        //
        // HOCH2CH2CH2CH2OH => HO(CH2)4OH
        let hydroxyl_1 = Chain::Vec(
            Vec::from([Chain::KV("H".into(), 1), Chain::KV("O".into(), 1)]),
            1,
        );
        let hydroxyl_4 = Chain::Vec(
            Vec::from([Chain::KV("O".into(), 1), Chain::KV("H".into(), 1)]),
            1,
        );
        let methylene = Chain::Vec(
            Vec::from([Chain::KV("C".into(), 1), Chain::KV("H".into(), 2)]),
            1,
        );
        let chain = Chain::Vec(
            Vec::from([
                hydroxyl_1,
                methylene.clone(),
                methylene.clone(),
                methylene.clone(),
                methylene,
                hydroxyl_4,
            ]),
            1,
        );
        assert_eq!(
            chain.to_string(),
            "HO(CH2)4OH",
            "Could not condense 1,4-butanediol into chemical formula"
        )
    }

    #[test]
    fn test_reverse_with_methane() {
        let methane = Chain::Vec(
            Vec::from([Chain::KV("C".into(), 1), Chain::KV("H".into(), 4)]),
            1,
        )
        .reversed();
        assert_eq!(
            methane.to_string(),
            "H4C".to_string(),
            "Compound reversing doesn't work"
        );
    }
}
