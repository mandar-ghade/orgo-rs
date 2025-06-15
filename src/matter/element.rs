use crate::constants::ELEMENTS_VEC;

#[allow(unused)]
#[derive(Clone, Debug, PartialEq)]
pub struct Element {
    pub number: u8,
}

#[allow(unused)]
impl Element {
    pub fn new(number: u8) -> Self {
        Element { number }
    }

    pub fn as_str(&self) -> &str {
        ELEMENTS_VEC
            .get(self.number as usize)
            .expect("Invalid element number received.")
    }
}
