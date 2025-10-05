use crate::matter::{atom::Atom, compound::deserializer::Chain};

trait Compound {
    fn get_chain() -> Chain<Atom>;
    fn get_backbone() -> Chain<Atom>;
}

pub trait LinearChain: Compound {
    // Get chain length
    fn chain_len() -> usize;
}

pub enum CompoundTypes {
    LinearChain,
    AminoAcid,
    Polypeptide,
}
//
// pub enum AminoAcid {
//     Ala,
//     Val,
//     Leu,
//     Iso,
//     Gly,
//     Tyr,
//     Trp,
//     Phe,
// }
