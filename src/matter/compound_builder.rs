use std::collections::{BTreeSet, HashMap};

use strum_macros::Display;
use thiserror::Error;

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
    location_to_idx: HashMap<Location, u8>,
    /// should theoretically be a size
    backbone: Vec<u8>,
    /// retains order
    /// ^ Side chains can have side-chains (unfortunately)
    /// ^ sort of like an undirected graph
    side_chains: HashMap<u8, BTreeSet<u8>>,
    // TODO: Ensure values != key or backbone idx
}

#[derive(Error, Clone, Debug, Display)]
pub enum CompoundBuilderError {
    LocationGenerationErr(String),
    SideChainErr(String),
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
    fn has_side_chain(&self, idx: u8) -> bool {
        if let Some(side_chain) = self.side_chains.get(&idx) {
            !side_chain.is_empty()
        } else {
            false
        }
    }

    fn get_remote_side_chain(
        &self,
        atom_idx: u8,
    ) -> Result<(u8, &Atom), CompoundBuilderError> {
        Ok((
            atom_idx,
            self.atoms
                .get(atom_idx as usize)
                .ok_or(CompoundBuilderError::SideChainErr(
                    "Couldn't find adjacent atom to side chain. (CompoundBuilder issue)".into(),
                ))?,
        ))
    }

    fn get_remote_side_chains(
        &self,
        idx: u8,
    ) -> Result<Vec<(u8, &Atom)>, CompoundBuilderError> {
        // Returns directly adjacent atoms to backbone atom
        // (not entire side chains)
        // TODO: Rework for larger chains
        self.side_chains
            .get(&idx)
            .unwrap_or(&BTreeSet::new())
            .iter()
            .map(|&i| self.get_remote_side_chain(i))
            .collect()
    }

    fn gen_locations(&mut self) -> Result<(), CompoundBuilderError> {
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
            let loc_i = *locations.get(i).ok_or_else(|| {
                CompoundBuilderError::LocationGenerationErr(
                    "Backbone loc not found (CompoundBuilder issue)".into(),
                )
            })?;
            let remote_atoms = self.get_remote_side_chains(i as u8)?;
            if remote_atoms.len() > 4 {
                todo!("Expanded octet prohibited (for now)");
            }
            let choices = Vec::from([(-1, 0), (1, 0), (0, 1), (0, -1)]);
            // LEFT, RIGHT, UP, DOWN
            for (idx, _) in remote_atoms {
                let valid_loc = choices.iter().find_map(|&(dx, dy)| {
                    let loc = loc_i.shift(dx, dy);
                    if locations_to_idx.contains_key(&loc) {
                        None // invalid because taken
                    } else {
                        Some(loc)
                    }
                });
                let Some(side_chain_loc) = valid_loc else {
                    return Err(CompoundBuilderError::SideChainErr(
                        "Location assignment failed (octet likely exceed for side chain atom)"
                            .into(),
                    ));
                };
                if locations.len() != idx as usize {
                    return Err(CompoundBuilderError::LocationGenerationErr(
                        "Location vector length - remote atom's index mismatch"
                            .into(),
                    ));
                }
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
        let mut working_idx = backbone_len as u8;
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
            let mut side_chain_len = self
                .side_chains
                .get(&(i as u8))
                .map(|c| c.len())
                .unwrap_or(0);
            while primary_chain_atoms + side_chain_len < 4 {
                let atom = Atom::new(1).unwrap();
                self.atoms.push(atom);
                self.side_chains
                    .entry(i as u8)
                    .or_default()
                    .insert(working_idx);
                working_idx += 1;
                side_chain_len = self
                    .side_chains
                    .get(&(i as u8))
                    .map(|c| c.len())
                    .unwrap_or(0);
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
        count: u8,
    ) -> Result<&mut Self, CompoundBuilderError> {
        self.atoms.clear();
        self.backbone.clear();
        for i in 0..count {
            self.backbone.push(i);
            self.atoms.push(Atom::new(6).unwrap());
        }
        self.side_chains.clear();
        self.satisfy_backbone_octets();
        self.gen_locations()?;
        Ok(self)
    }
}
