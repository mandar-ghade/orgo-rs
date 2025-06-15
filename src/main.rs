use matter::compound_builder::CompoundBuilder;
mod constants;
mod matter;

// CH3(CH2)5CH3
// CH3(CH(CH3)CH2)CH3

fn main() {
    let cmp = CompoundBuilder::new().linear_chain(6).build();
    dbg!(&cmp);
    dbg!(cmp.to_string());
    // Compound::parse("(CH3)2CH(CH3)ZnCl").unwrap();
}
