use std::fmt;

use super::element::Element;

#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    element: Element,
    neutrons: u8,
    electrons: u8, // same as protons
}

#[allow(unused)]
impl Atom {
    pub fn new(element_num: u8) -> Option<Self> {
        let element = Element::new(element_num)?;
        Some(Atom {
            electrons: element.number,
            element,
            neutrons: 0,
        })
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.element.as_str())
    }
}
