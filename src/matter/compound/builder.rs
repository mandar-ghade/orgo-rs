use std::collections::{BTreeSet, HashMap};

use crate::matter::{
    atom::Atom,
    compound::{Compound, Location},
};

#[derive(Debug)]
pub struct CompoundBuilder {
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

pub type CompoundBuilderResult<T> = Result<T, CompoundBuilderError>;

#[derive(thiserror::Error, strum_macros::Display, Clone, Debug)]
pub enum CompoundBuilderError {
    SideChainError(String),
}

impl CompoundBuilder {
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            backbone: Vec::new(),
            locations: Vec::new(),
            location_to_idx: HashMap::new(),
            side_chains: HashMap::new(),
        }
    }

    #[allow(dead_code)]
    fn recompute_backbone(&self) {
        // TODO: Dijkstra's implementation
        todo!()
    }

    #[allow(dead_code)]
    fn has_side_chain(&self, idx: usize) -> bool {
        if let Some(side_chain) = self.side_chains.get(&idx) {
            !side_chain.is_empty()
        } else {
            false
        }
    }

    #[allow(dead_code)]
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

    #[allow(dead_code)]
    fn get_location_unsafe(&self, idx: usize) -> Location {
        if let Some(&loc) = self.locations.get(idx) {
            loc
        } else {
            panic!(
                "Invalid location index: {} (length is {})",
                idx,
                self.locations.len(),
            );
        }
    }

    fn get_remote_side_chains(
        &self,
        idx: usize,
    ) -> CompoundBuilderResult<Vec<(usize, &Atom)>> {
        // Returns directly adjacent atoms to backbone atom
        // (not entire side chains)
        // TODO: Rework for larger chains
        let Some(side_chain) = self.side_chains.get(&idx) else {
            return Ok(Vec::new());
        };
        let mut chains = Vec::new();
        for &i in side_chain {
            let atom = self.atoms.get(i).ok_or(
                CompoundBuilderError::SideChainError(
                    "Couldn't find adjacent atom to side chain.".into(),
                ),
            )?;
            chains.push((i, atom));
        }
        Ok(chains)
    }

    fn gen_locations(&mut self) -> CompoundBuilderResult<()> {
        // Called after octets are completed.
        let mut locations = Vec::new();
        let mut locations_to_idx = HashMap::new();
        let backbone_len = self.backbone.len();
        self.backbone.iter().for_each(|&i| {
            let loc = Location::new(i as i16, 0);
            locations.push(loc);
            locations_to_idx.insert(loc, i);
        });
        for i in 0..backbone_len {
            let base_loc = *locations.get(i).expect("Backbone loc not found.");
            let remote_atoms = self.get_remote_side_chains(i)?;
            if remote_atoms.len() > 4 {
                todo!("Expanded octet prohibited (for now)");
            }
            let choices = Vec::from([(-1, 0), (1, 0), (0, 1), (0, -1)]);
            // LEFT, RIGHT, UP, DOWN
            for (idx, _) in remote_atoms {
                let side_chain_loc = choices
                    .iter()
                    .find_map(|&(dx, dy)| {
                        let loc = base_loc.shift(dx, dy);
                        if locations_to_idx.contains_key(&loc) {
                            None // invalid because taken
                        } else {
                            Some(loc)
                        }
                    })
                    .ok_or_else(|| {
                        CompoundBuilderError::SideChainError(
                            // octet likely exceed for side chain atom
                            "Location assignment failed".into(),
                        )
                    })?;
                assert!(
                    locations.len() == idx,
                    "Location vector length - remote atom's index mismatch",
                );
                locations.push(side_chain_loc.clone());
                locations_to_idx.insert(side_chain_loc, idx);
            }
        }
        self.locations = locations;
        self.location_to_idx = locations_to_idx;
        Ok(())
    }

    fn satisfy_backbone_octets(&mut self) {
        // ONLY for linear compound generation
        // Satisfy backbone octets
        let backbone_len = self.backbone.len();
        let mut working_idx = backbone_len;
        for i in 0..backbone_len {
            let mut primary_chain_atoms = 0;
            if i < backbone_len - 1 {
                // something comes after it
                primary_chain_atoms += 1;
            }
            if i > 0 {
                // something comes before it
                primary_chain_atoms += 1;
            }
            let mut side_chain_len =
                self.side_chains.get(&i).map(|c| c.len()).unwrap_or(0);
            while primary_chain_atoms + side_chain_len < 4 {
                self.atoms.push(Atom::hydrogen());
                self.side_chains.entry(i).or_default().insert(working_idx);
                working_idx += 1;
                side_chain_len =
                    self.side_chains.get(&i).map(|c| c.len()).unwrap_or(0);
                // I could simplify this line
            }
        }
    }

    pub fn build(&mut self) -> Compound {
        Compound::new(
            self.atoms.clone(),
            self.locations.clone(),
            self.location_to_idx.clone(),
            self.backbone.clone(),
            self.side_chains.clone(),
        )
    }

    pub fn linear_chain(
        &mut self,
        count: usize,
    ) -> CompoundBuilderResult<&mut Self> {
        self.atoms.clear();
        self.backbone.clear();
        for i in 0..count {
            self.atoms.push(Atom::carbon());
            self.backbone.push(i);
        }
        // self.side_chains.clear();
        self.satisfy_backbone_octets();
        self.gen_locations()?;
        Ok(self)
    }
}
