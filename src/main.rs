//use thiserror::Error;

use std::{
    borrow::BorrowMut,
    cell::{BorrowMutError, RefCell},
    ops::Index,
    rc::Rc,
};

use strum_macros::Display;
use thiserror::Error;

enum Chirality {
    R,
    S,
}

#[allow(unused)]
#[derive(Clone, Debug)]
struct Element {
    number: u8,
}

impl Element {
    fn new(number: u8) -> Self {
        Element { number }
    }
}

#[allow(unused)]
#[derive(Clone, Debug)]
struct Atom {
    element: Element,
    neutrons: u8,
    electrons: u8,   // same as protons
    oxi: Option<u8>, // oxidation state
}

impl Atom {
    fn new(element_num: u8) -> Self {
        let element = Element::new(element_num);
        Atom {
            electrons: element.number,
            element,
            neutrons: 0,
            oxi: None,
        }
    }
}

#[derive(Error, Clone, Debug, Display)]
enum ParsingError {
    MoleculeParsingError(String),
    CompoundParsingError(String),
}

#[derive(Error, Clone, Debug, Display)]
enum CompoundUsageError {
    PushToSubstituents(String),
}

#[allow(unused)]
#[derive(Clone, Debug)]
enum Particle {
    Atom(Atom),
    Molecule(Vec<Atom>),
}

impl Particle {
    fn new_atom(atom: Atom) -> Particle {
        Particle::Atom(atom)
    }

    fn new_molecule(atoms: Vec<Atom>) -> Particle {
        Particle::Molecule(atoms)
    }
}

#[allow(unused)]
#[derive(Clone, Debug)]
struct Compound {
    // this is sort of like an iterator
    center: Atom,
    subst: Vec<Particle>,                    // substituents
    side_chains: Vec<Rc<RefCell<Compound>>>, // previous compounds
}

impl Compound {
    fn new(
        center: Atom,
        subst: Vec<Particle>,
        side_chains: Vec<Rc<RefCell<Compound>>>,
    ) -> Result<Compound, ParsingError> {
        Ok(Compound {
            center,
            subst,
            side_chains,
        })
        .map_err(ParsingError::CompoundParsingError)
    }

    fn as_box(
        center: Atom,
        subst: Vec<Particle>,
        side_chains: Vec<Rc<RefCell<Compound>>>,
    ) -> Result<Box<Compound>, ParsingError> {
        Ok(Box::new(Compound {
            center,
            subst,
            side_chains,
        }))
        .map_err(ParsingError::CompoundParsingError)
    }

    fn as_rc(
        center: Atom,
        subst: Vec<Particle>,
        side_chains: Vec<Rc<RefCell<Compound>>>,
    ) -> Result<Rc<RefCell<Compound>>, ParsingError> {
        Ok(Rc::new(RefCell::new(Compound {
            center,
            subst,
            side_chains,
        })))
        .map_err(ParsingError::CompoundParsingError)
    }

    fn as_rc_with_center(center: Atom) -> Result<Rc<RefCell<Compound>>, ParsingError> {
        Ok(Rc::new(RefCell::new(Compound {
            center,
            subst: vec![],
            side_chains: vec![],
        })))
        .map_err(ParsingError::CompoundParsingError)
    }
}

fn push_subst(cmp: Rc<RefCell<Compound>>, p: Particle) -> Result<(), BorrowMutError> {
    (*cmp).try_borrow_mut()?.subst.push(p);
    Ok(())
}

fn push_side_chain(
    cmp: Rc<RefCell<Compound>>,
    chain: Rc<RefCell<Compound>>,
) -> Result<(), BorrowMutError> {
    (*cmp).try_borrow_mut()?.side_chains.push(chain);
    Ok(())
}

fn main() {
    let cmp = Compound::as_rc_with_center(Atom::new(3)).expect("Compound was expected");
    {
        let _ = push_subst(Rc::clone(&cmp), Particle::new_atom(Atom::new(18)));
        let _ = push_subst(
            Rc::clone(&cmp),
            Particle::new_molecule(vec![Atom::new(7), Atom::new(7)]),
        );
    }
    let cmp2 = Compound::as_rc_with_center(Atom::new(4)).expect("Compound was expected");
    let cmp3 = Compound::as_rc_with_center(Atom::new(5)).expect("Compound was expected");
    let _ = push_side_chain(Rc::clone(&cmp), Rc::clone(&cmp2));
    let _ = push_side_chain(Rc::clone(&cmp), Rc::clone(&cmp3));
    let _ = push_side_chain(Rc::clone(&cmp2), Rc::clone(&cmp));
    {
        let _ = push_subst(Rc::clone(&cmp2), Particle::new_atom(Atom::new(255)));
    }
    dbg!((*cmp).borrow_mut());
}
