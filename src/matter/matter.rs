use super::atom::Atom;

enum Matter {
    Compound(Vec<Matter>),
    Atom(Atom),
}
