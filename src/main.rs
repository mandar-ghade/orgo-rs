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
    dbg!(&cmp);
    let chain = Chain::from(&cmp).group().minimize();
    dbg!(&chain);
    dbg!(chain.to_string());
}
