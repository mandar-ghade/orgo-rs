use std::fmt;

use super::element::Element;

#[derive(Clone, Debug)]
pub struct Atom {
    element: Element,
    pub neutrons: u8,
    pub electrons: u8, // same as protons
}

#[allow(dead_code)]
impl Atom {
    pub fn new(element_num: u8) -> Option<Self> {
        let element = Element::new(element_num)?;
        Some(Atom {
            electrons: element.number,
            element,
            neutrons: 0,
        })
    }

    pub fn new_unchecked(element_num: u8) -> Self {
        let element = Element::new_unchecked(element_num);
        Atom {
            electrons: element.number,
            element,
            neutrons: 0,
        }
    }

    pub fn hydrogen() -> Self {
        Atom::new_unchecked(1)
    }

    pub fn carbon() -> Self {
        Atom::new_unchecked(6)
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.element.as_str())
    }
}

impl PartialEq for Atom {
    fn eq(&self, other: &Self) -> bool {
        self.element == other.element
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atom_element_number_bound_checks() {
        assert!(Atom::new(0).is_none());
        assert!(Atom::new(1).is_some());
        assert!(Atom::new(118).is_some());
        assert!(Atom::new(119).is_none());
        assert!(Atom::new(u8::MAX).is_none());
    }

    #[test]
    fn atom_to_string() {
        assert_eq!(Atom::new_unchecked(1).to_string(), "H");
        assert_eq!(Atom::new_unchecked(2).to_string(), "He");
        assert_eq!(Atom::new_unchecked(3).to_string(), "Li");
        assert_eq!(Atom::new_unchecked(118).to_string(), "Og");
    }
}
