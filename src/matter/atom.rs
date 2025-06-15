use crate::matter::element::Element;

use super::matter_trait::Matter;

#[allow(unused)]
#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    element: Element,
    neutrons: u8,
    electrons: u8, // same as protons
}

#[allow(unused)]
impl Atom {
    pub fn new(element_num: u8) -> Self {
        let element = Element::new(element_num);
        Atom {
            electrons: element.number,
            element,
            neutrons: 0,
        }
    }
}

impl Matter for Atom {
    fn to_string(&self) -> String {
        self.element.as_str().to_string()
    }
}
