use std::path::Path;

#[derive(Debug, Default)]
#[allow(unused)]
struct Atom {
    element: String,
    atype: String,
    x: f64,
    y: f64,
    z: f64,
    formal_charge: isize,
    real_charge: f64,
    /// explicit H count
    hexp: usize,
    /// total H count
    htot: usize,
    neighbor_count: usize,
    ring_count: usize,
    arom: bool,
    /// potentially aromatic in a query structure
    q_arom: bool,
    stereo_care: bool,
    heavy: bool,
    metal: bool,
    nvalences: usize,
    tag: bool,
    nucleon_number: usize,
    radical_type: usize,
}

#[derive(Debug, Default)]
#[allow(unused)]
struct Bond {
    a1: usize,
    a2: usize,
    btype: char,
    ring_count: usize,
    arom: bool,
    q_arom: bool, //  potentially aromatic in a query structure
    topo: usize,  //  see MDL file description
    stereo: usize,
    mdl_stereo: usize,
}

#[derive(Debug)]
#[allow(unused)]
struct Input {
    mol_name: String,
    mol_comment: String,
    n_c_tot: usize,
    n_o_tot: usize,
    n_n_tot: usize,
    n_heavy: usize,
    heavy_bonds: usize,
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
}

impl Input {
    /// load an SDF file from `filename`
    fn load(filename: impl AsRef<Path>) -> Self {
        let s = std::fs::read_to_string(filename).unwrap();
        let mut lines = s.lines();
        let mol_name = lines.next().unwrap(); // line 1
        let _ = lines.next().unwrap(); // discard line 2
        let mol_comment = lines.next().unwrap(); // line 3
        let mut info = lines.next().unwrap().split_ascii_whitespace();
        let n_atoms: usize = info.next().unwrap().trim().parse().unwrap();
        let n_bonds: usize = info.next().unwrap().trim().parse().unwrap();
        // TODO check for chirality flag in char 15? maybe 14 if pascal is
        // 1-indexed
        // total number of C atoms
        let mut n_c_tot = 0;
        let mut n_o_tot = 0;
        let mut n_n_tot = 0;
        let mut n_heavy = 0;
        let mut atoms = Vec::with_capacity(n_atoms);
        for _ in 0..n_atoms {
            let line = lines.next().unwrap();
            let sp: Vec<_> = line.split_ascii_whitespace().collect();
            let elem_str = &sp[3];
            if elem_str == &"C" {
                n_c_tot += 1;
            }
            if elem_str == &"O" {
                n_o_tot += 1;
            }
            if elem_str == &"N" {
                n_n_tot += 1;
            }
            let new_atom_type = convert_mdl_type(elem_str);
            let x: f64 = sp[0].parse().unwrap();
            let y: f64 = sp[1].parse().unwrap();
            let z: f64 = sp[2].parse().unwrap();

            let chg: f64 = sp[4].parse().unwrap();
            let is_heavy = is_heavy_atom(elem_str);
            // TODO skipping is_metal and is_trueheavyatom
            if is_heavy {
                n_heavy += 1;
            }
            let nvalences = get_valence(elem_str);
            // TODO skipping some deuterium and tritium stuff for now
            atoms.push(Atom {
                element: elem_str.to_string(),
                atype: new_atom_type,
                x,
                y,
                z,
                formal_charge: chg.round() as isize,
                real_charge: chg,
                nvalences,
                ..Default::default()
            });
        }

        let mut bonds = Vec::with_capacity(n_bonds);
        for line in lines.take(n_bonds) {
            let sp: Vec<_> = line.split_ascii_whitespace().collect();
            let a1 = sp[0].parse::<usize>().unwrap() - 1;
            let a2 = sp[1].parse::<usize>().unwrap() - 1;
            bonds.push(Bond {
                a1,
                a2,
                btype: match sp[2] {
                    "1" => 'S', // single
                    "2" => 'D', // double
                    "3" => 'T', // triple
                    "4" => 'A', // aromatic
                    "5" => 'l', // single or double
                    "6" => 's', // single or aromatic
                    "7" => 'd', // double or aromatic
                    "8" => 'a', // any
                    "9" => 'a', // any in JSME;  v0.5b
                    _ => unimplemented!(),
                },
                // TODO skipping aromaticity reading and topology
                ..Default::default()
            });
        }

        let mut heavy_bonds = 0;
        for bond in &bonds {
            if atoms[bond.a1].heavy && atoms[bond.a2].heavy {
                heavy_bonds += 1;
            }
        }
        Self {
            mol_name: mol_name.into(),
            mol_comment: mol_comment.into(),
            n_c_tot,
            n_o_tot,
            n_n_tot,
            n_heavy,
            heavy_bonds,
            atoms,
            bonds,
        }
    }
}

fn get_valence(elem_str: &str) -> usize {
    match elem_str {
        "H" => 1,
        "D" => 1, // v0.3n
        "C" => 4,
        "N" => 3,
        "O" => 2,
        "S" => 2,
        "SE" => 2,
        "TE" => 2,
        "P" => 3,
        "F" => 1,
        "CL" => 1,
        "BR" => 1,
        "I" => 1,
        "B" => 3,
        "LI" => 1,
        "NA" => 1,
        "K" => 1,
        "CA" => 2,
        "SR" => 2,
        "MG" => 2,
        "FE" => 3,
        "MN" => 2,
        "HG" => 2,
        "SI" => 4,
        "SN" => 4,
        "ZN" => 2,
        "CU" => 2,
        "A" => 4,
        "Q" => 4,
        _ => unimplemented!(),
    }
}

fn is_heavy_atom(elem_str: &str) -> bool {
    elem_str == "H"
}

fn convert_mdl_type(elem_str: &str) -> String {
    match elem_str {
        "H" => "H",
        "C" => "C3",
        "O" => "O2",
        "N" => "N3",
        "F" => "F",
        "Cl" => "CL",
        "Br" => "BR",
        "I" => "I",
        "Al" => "AL",
        "ANY" => "A",
        "Ca" => "CA",
        "Du" => "DU",
        "K" => "K",
        "Li" => "LI",
        "LP" => "LP",
        "Na" => "NA",
        "S" => "S3",
        "Si" => "SI",
        "P" => "P4",
        "A" => "A",
        "Q" => "Q",
        _ => "DU",
    }
    .to_owned()
}

fn main() {
    dbg!(Input::load("/tmp/input.sdf"));
}
