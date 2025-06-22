use crate::matter::compound::{builder::CompoundBuilder, deserializer::Chain};

mod constants;
mod matter;

// CH3(CH2)5CH3
// CH3(CH(CH3)CH2)CH3

fn main() {
    let cmp = CompoundBuilder::new()
        .linear_chain(6)
        .expect("Linear chain didn't work")
        .build();
    let chain = Chain::from(&cmp);
    dbg!(chain.to_string());
}
