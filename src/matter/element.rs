use crate::constants::{ELEMENTS, ELEMENTS_VEC};

#[derive(Clone, Debug, PartialEq)]
pub struct Element {
    pub number: u8,
}

#[allow(dead_code)]
impl Element {
    pub fn new(number: u8) -> Option<Self> {
        if number != 0 && number as usize <= ELEMENTS_VEC.len() {
            Some(Element { number })
        } else {
            None
        }
    }

    pub fn from_str(input_str: &str) -> Option<Self> {
        if !ELEMENTS.contains(input_str) {
            None
        } else {
            ELEMENTS_VEC.iter().enumerate().find_map(|(i, e_str)| {
                if *e_str == input_str {
                    Some(Element::new((i as u8) + 1))
                } else {
                    None
                }
            })?
        }
    }

    pub fn new_unchecked(number: u8) -> Self {
        Element::new(number).expect(&format!(
            "Invalid element number: {}. `Element::new_unchecked` expects a \
                valid atomic number.",
            number
        ))
    }

    pub fn hydrogen() -> Self {
        Element::new_unchecked(1)
    }

    pub fn carbon() -> Self {
        Element::new_unchecked(6)
    }

    pub fn as_str(&self) -> &str {
        ELEMENTS_VEC.get(self.number as usize - 1).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn element_number_bound_checks() {
        assert!(Element::new(0).is_none());
        assert!(Element::new(1).is_some());
        assert!(Element::new(118).is_some());
        assert!(Element::new(119).is_none());
        assert!(Element::new(u8::MAX).is_none());
    }

    #[test]
    fn element_as_str() {
        assert_eq!(Element::new_unchecked(1).as_str(), "H");
        assert_eq!(Element::new_unchecked(2).as_str(), "He");
        assert_eq!(Element::new_unchecked(3).as_str(), "Li");
        assert_eq!(Element::new_unchecked(118).as_str(), "Og");
    }
}
