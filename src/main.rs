use std::str::FromStr;

use matter::compound::{Compound, CompoundError};
mod matter;
mod misc;

fn parse_str(input_str: String) -> Result<Compound, CompoundError> {
    Compound::from_str(input_str.as_str())
}

// CH3(CH2)5CH3
// CH3(CH(CH3)CH2)CH3

fn main() {
    parse_str("(CH3)2CH(CH3)ZnCl".into()).expect("Error");
}
