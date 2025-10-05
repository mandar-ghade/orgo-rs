use crate::{
    matter::{
        atom::Atom,
        compound::{
            builder::{CompoundBuilder, CompoundBuilderResult},
            deserializer::Chain,
            Compound, CompoundResult,
        },
    },
    other::qm_model::Configuration,
};

mod constants;
mod matter;
mod other;

// CH3(CH2)5CH3
// CH3(CH(CH3)CH2)CH3

fn run_cmp_builder() -> CompoundBuilderResult<Compound> {
    let mut cmp_builder = CompoundBuilder::new();
    Ok(cmp_builder.linear_chain(6)?.brominate(6)?.build())
}

fn main() {
    let atom = Atom::from_str_unchecked("Cr");
    let cfg = Configuration::from_atom(&atom).build();
    dbg!(&cfg);
    dbg!(&cfg.to_string());
    // let cmp = run_cmp_builder().expect("Compound was expected");
    // dbg!(&cmp);
    // dbg!(&cmp.to_string());
}
