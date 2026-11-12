//! Molecular structure data model.
//!
//! [`Structure`] is the neutral in-memory representation consumed by the SES
//! engine. Parser backends (PDB, mmCIF, XYZ, XYZR) all produce this type; no
//! third-party parser types leak into the public API.

use std::fmt;
use std::ops::Range;

use crate::error::{Result, SesError};

/// Chemical element, reduced to what the surface code needs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Element {
    /// Hydrogen.
    H,
    /// Carbon.
    C,
    /// Nitrogen.
    N,
    /// Oxygen.
    O,
    /// Phosphorus.
    P,
    /// Sulfur.
    S,
    /// Selenium.
    Se,
    /// Fluorine.
    F,
    /// Chlorine.
    Cl,
    /// Bromine.
    Br,
    /// Iodine.
    I,
    /// Sodium.
    Na,
    /// Magnesium.
    Mg,
    /// Potassium.
    K,
    /// Calcium.
    Ca,
    /// Iron.
    Fe,
    /// Zinc.
    Zn,
    /// Any other element, identified by atomic number.
    Other(u8),
}

impl Element {
    /// Parse an element symbol (case-insensitive), e.g. `"C"`, `"ZN"`, `"se"`.
    pub fn from_symbol(s: &str) -> Option<Self> {
        let s = s.trim();
        if s.is_empty() {
            return None;
        }
        // Normalize to Title case: first upper, rest lower.
        let mut chars = s.chars();
        let first = chars.next()?.to_ascii_uppercase();
        let rest: String = chars.as_str().to_ascii_lowercase();
        let symbol = format!("{first}{rest}");
        Some(match symbol.as_str() {
            "H" => Self::H,
            "C" => Self::C,
            "N" => Self::N,
            "O" => Self::O,
            "P" => Self::P,
            "S" => Self::S,
            "Se" => Self::Se,
            "F" => Self::F,
            "Cl" => Self::Cl,
            "Br" => Self::Br,
            "I" => Self::I,
            "Na" => Self::Na,
            "Mg" => Self::Mg,
            "K" => Self::K,
            "Ca" => Self::Ca,
            "Fe" => Self::Fe,
            "Zn" => Self::Zn,
            _ => {
                // Accept atomic numbers as digits.
                if let Ok(z) = s.parse::<u8>() {
                    Self::Other(z)
                } else {
                    return None;
                }
            }
        })
    }

    /// Standard element symbol (Title case).
    pub fn symbol(&self) -> String {
        match self {
            Self::H => "H".into(),
            Self::C => "C".into(),
            Self::N => "N".into(),
            Self::O => "O".into(),
            Self::P => "P".into(),
            Self::S => "S".into(),
            Self::Se => "Se".into(),
            Self::F => "F".into(),
            Self::Cl => "Cl".into(),
            Self::Br => "Br".into(),
            Self::I => "I".into(),
            Self::Na => "Na".into(),
            Self::Mg => "Mg".into(),
            Self::K => "K".into(),
            Self::Ca => "Ca".into(),
            Self::Fe => "Fe".into(),
            Self::Zn => "Zn".into(),
            Self::Other(z) => format!("Z{z}"),
        }
    }

    /// Atomic number.
    pub fn atomic_number(&self) -> u8 {
        match self {
            Self::H => 1,
            Self::C => 6,
            Self::N => 7,
            Self::O => 8,
            Self::P => 15,
            Self::S => 16,
            Self::Se => 34,
            Self::F => 9,
            Self::Cl => 17,
            Self::Br => 35,
            Self::I => 53,
            Self::Na => 11,
            Self::Mg => 12,
            Self::K => 19,
            Self::Ca => 20,
            Self::Fe => 26,
            Self::Zn => 30,
            Self::Other(z) => *z,
        }
    }

    /// `true` for hydrogen (including deuterium/tritium treated as H).
    pub fn is_hydrogen(&self) -> bool {
        matches!(self, Self::H | Self::Other(1))
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

/// Residue identifier: chain, residue name, sequence number, insertion code.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ResidueId {
    /// Chain identifier (e.g. `"A"`). Empty when the format has no chains (XYZ).
    pub chain: String,
    /// Residue name (e.g. `"ALA"`, `"HOH"`).
    pub name: String,
    /// Residue sequence number.
    pub seq: i32,
    /// Insertion code (PDB col 27), empty when absent.
    pub icode: char,
}

impl fmt::Display for ResidueId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}{}{}", self.chain, self.name, self.seq, self.icode)
    }
}

/// A single atom with the data the SES engine needs.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Atom {
    /// Cartesian position in Ångström.
    pub position: [f64; 3],
    /// van der Waals radius in Ångström.
    pub radius: f64,
    /// Chemical element, when known.
    pub element: Option<Element>,
    /// Serial number from the source file (1-based as in PDB).
    pub serial: u32,
    /// Atom name (PDB cols 13–16), trimmed.
    pub name: String,
    /// Residue identity, when the format provides one.
    pub residue: Option<ResidueId>,
    /// `true` for HETATM records / non-polymer atoms.
    pub hetero: bool,
    /// Occupancy (PDB), default 1.0.
    pub occupancy: f64,
    /// Alternate location indicator, empty when absent.
    pub altloc: char,
}

impl Atom {
    /// Create an atom with default occupancy/altloc and no residue.
    pub fn new(position: [f64; 3], radius: f64, element: Option<Element>, serial: u32) -> Self {
        Self {
            position,
            radius,
            element,
            serial,
            name: String::new(),
            residue: None,
            hetero: false,
            occupancy: 1.0,
            altloc: '\0',
        }
    }

    /// Squared distance to a point (avoids `sqrt` in inner loops).
    pub fn distance_squared_to(&self, p: [f64; 3]) -> f64 {
        let dx = self.position[0] - p[0];
        let dy = self.position[1] - p[1];
        let dz = self.position[2] - p[2];
        dx * dx + dy * dy + dz * dz
    }
}

/// Provenance of a parsed structure.
#[derive(Debug, Clone, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SourceInfo {
    /// Path or description of the input.
    pub path: String,
    /// Detected format name (`"pdb"`, `"mmcif"`, `"xyz"`, `"xyzr"`).
    pub format: String,
    /// Number of MODEL blocks found (1 when the format has no models).
    pub model_count: usize,
}

/// A parsed molecular structure.
///
/// Atoms are stored in file order. Multi-model inputs (NMR ensembles) partition
/// the atom list into per-model ranges exposed by [`Structure::models`].
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Structure {
    atoms: Vec<Atom>,
    models: Vec<Range<usize>>,
    source: SourceInfo,
}

impl Structure {
    /// Build a structure from a flat atom list (single model).
    pub fn from_atoms(atoms: Vec<Atom>, source: SourceInfo) -> Self {
        let n = atoms.len();
        Self {
            atoms,
            models: std::iter::once(0..n).collect(),
            source,
        }
    }

    /// Build a multi-model structure. `models` are index ranges into `atoms`.
    pub fn from_models(
        atoms: Vec<Atom>,
        models: impl IntoIterator<Item = Range<usize>>,
        source: SourceInfo,
    ) -> Self {
        Self {
            atoms,
            models: models.into_iter().collect(),
            source,
        }
    }

    /// All atoms of all models.
    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    /// Number of atoms across all models.
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// `true` when the structure has no atoms.
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Per-model atom index ranges.
    pub fn models(&self) -> &[Range<usize>] {
        &self.models
    }

    /// Atoms of model `idx` (0-based).
    pub fn model_atoms(&self, idx: usize) -> Result<&[Atom]> {
        self.models
            .get(idx)
            .map(|r| &self.atoms[r.clone()])
            .ok_or_else(|| {
                SesError::invalid_param(
                    "model",
                    idx,
                    format!("structure has {} model(s)", self.models.len()),
                )
            })
    }

    /// Source metadata.
    pub fn source(&self) -> &SourceInfo {
        &self.source
    }

    /// Axis-aligned bounding box of atom centers and radii: `(min, max)` in Å.
    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for a in &self.atoms {
            for k in 0..3 {
                let r = a.radius;
                lo[k] = lo[k].min(a.position[k] - r);
                hi[k] = hi[k].max(a.position[k] + r);
            }
        }
        if self.atoms.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        (lo, hi)
    }

    /// Replace model ranges after parsing (used by parsers).
    pub(crate) fn set_models(&mut self, models: impl IntoIterator<Item = Range<usize>>) {
        self.models = models.into_iter().collect();
    }

    /// Mutable access to atoms for post-processing (radii, dedup).
    pub(crate) fn atoms_mut(&mut self) -> &mut Vec<Atom> {
        &mut self.atoms
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn element_symbols_roundtrip() {
        for sym in ["H", "C", "N", "O", "S", "Se", "ZN", "cl", "fe"] {
            let e = Element::from_symbol(sym).expect(sym);
            assert_eq!(e.symbol().to_lowercase(), sym.to_lowercase());
        }
        assert!(Element::from_symbol("").is_none());
        assert!(Element::from_symbol("Xx").is_none());
    }

    #[test]
    fn hydrogen_detection() {
        assert!(Element::H.is_hydrogen());
        assert!(!Element::C.is_hydrogen());
    }

    #[test]
    fn structure_models() {
        let atoms = vec![
            Atom::new([0.0; 3], 1.0, Some(Element::C), 1),
            Atom::new([1.0; 3], 1.0, Some(Element::C), 2),
            Atom::new([2.0; 3], 1.0, Some(Element::C), 3),
        ];
        let s = Structure::from_models(
            atoms,
            [0..2, 2..3],
            SourceInfo {
                path: "test".into(),
                format: "pdb".into(),
                model_count: 2,
            },
        );
        assert_eq!(s.len(), 3);
        assert_eq!(s.models().len(), 2);
        assert_eq!(s.model_atoms(0).unwrap().len(), 2);
        assert_eq!(s.model_atoms(1).unwrap().len(), 1);
        assert!(s.model_atoms(2).is_err());
    }

    #[test]
    fn bounding_box_includes_radius() {
        let atoms = vec![Atom::new([0.0; 3], 2.0, Some(Element::C), 1)];
        let s = Structure::from_atoms(atoms, SourceInfo::default());
        let (lo, hi) = s.bounding_box();
        assert_eq!(lo, [-2.0; 3]);
        assert_eq!(hi, [2.0; 3]);
    }
}
