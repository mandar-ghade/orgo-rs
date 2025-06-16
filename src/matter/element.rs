use crate::constants::ELEMENTS_VEC;

#[derive(Clone, Debug, PartialEq)]
pub struct Element {
    pub number: u8,
}

impl Element {
    pub fn new(number: u8) -> Option<Self> {
        if number != 0 && number as usize <= ELEMENTS_VEC.len() {
            Some(Element { number })
        } else {
            None
        }
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
        assert_eq!(Element::new(1).unwrap().as_str(), "H");
        assert_eq!(Element::new(2).unwrap().as_str(), "He");
        assert_eq!(Element::new(3).unwrap().as_str(), "Li");
        assert_eq!(Element::new(118).unwrap().as_str(), "Og");
    }
}
