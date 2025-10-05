use crate::{
    matter::{
        atom::Atom,
        compound::{
            builder::{CompoundBuilder, CompoundBuilderResult},
            Compound,
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
    let mut atom = Atom::from_str_unchecked("Al");
    dbg!(atom.get_config());
    // let cmp = run_cmp_builder().expect("Compound was expected");
    // dbg!(&cmp);
    // dbg!(&cmp.to_string());
}
