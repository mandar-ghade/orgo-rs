use std::collections::{HashMap, HashSet};

use lazy_static::lazy_static;

lazy_static! {
    pub static ref ELEMENTS_VEC: Vec<&'static str> = vec![
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
        "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
        "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
        "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
        "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
        "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
        "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg",
        "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
    ];
    pub static ref ELEMENTS: HashSet<&'static str> = HashSet::from([
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
        "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
        "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
        "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
        "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
        "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
        "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg",
        "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
    ]);
    pub static ref ELEMENT_WEIGHTS: [f64; 118] = [
        1.008, 4.003, 6.941, 9.012, 10.81, 12.01, 14.01, 16.00, 19.00, 20.18, 22.99, 24.31, 26.98,
        28.09, 30.97, 32.07, 35.45, 39.95, 39.10, 40.08, 44.96, 47.88, 50.94, 52.00, 54.94, 55.85,
        58.93, 58.69, 63.55, 65.39, 69.72, 72.59, 74.92, 78.96, 79.90, 83.80, 85.47, 87.62, 88.91,
        91.22, 92.91, 95.94, 99.0, 101.1, 102.9, 106.4, 107.9, 112.4, 114.8, 118.7, 121.8, 127.6,
        126.9, 131.3, 132.9, 137.3, 138.9, 140.1, 140.9, 144.2, 145.0, 150.4, 152.0, 157.2, 158.9,
        162.5, 164.9, 167.3, 168.9, 173.0, 175.0, 178.5, 180.9, 183.9, 186.2, 190.2, 192.2, 195.1,
        197.0, 200.6, 204.4, 207.2, 209.0, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.0, 231.0,
        238.0, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0, 258.0, 259.0, 260.0, 261.0,
        262.0, 266.0, 264.0, 277.0, 268.0, 281.0, 272.0, 285.0, 284.0, 289.0, 288.0, 291.0, 294.0,
        294.0,
    ];
}
