use matter::compound::Compound;
mod constants;
mod matter;

// CH3(CH2)5CH3
// CH3(CH(CH3)CH2)CH3

fn main() {
    Compound::parse("(CH3)2CH(CH3)ZnCl").unwrap();
}
