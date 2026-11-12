//! PDB and mmCIF parsing via the `pdbtbx` crate.
//!
//! `pdbtbx` models altLocs as one `Conformer` per alternate location within a
//! residue. This adapter flattens the hierarchy into [`Atom`] records and
//! applies the selection policy.

use std::collections::BTreeMap;
use std::path::Path;

use super::{AltlocPolicy, ModelSelection, ParseOptions, ParseReport, apply_selection};
use crate::error::{Result, SesError};
use crate::structure::{Atom, Element, ResidueId, SourceInfo, Structure};

pub(crate) fn parse_pdb_like(path: &Path, opts: &ParseOptions) -> Result<(Structure, ParseReport)> {
    let fmt = super::detect_format(path)?;
    let mut report = ParseReport::default();

    let path_str = path.to_string_lossy().into_owned();
    let (pdb, errors) = pdbtbx::open(&path_str).map_err(|errs| {
        let msg = errs
            .iter()
            .map(|e| format!("{e:?}"))
            .collect::<Vec<_>>()
            .join("; ");
        SesError::parse(path, None, format!("pdbtbx could not open file: {msg}"))
    })?;

    // Surface the first parse error if the file is clearly broken and empty.
    if pdb.model_count() == 0 {
        let msg = errors
            .iter()
            .map(|e| format!("{e:?}"))
            .collect::<Vec<_>>()
            .join("; ");
        return Err(SesError::parse(
            path,
            None,
            if msg.is_empty() {
                "no models found".to_string()
            } else {
                msg
            },
        ));
    }

    let model_count = pdb.model_count();

    // Collect atoms per model, applying altLoc policy within each residue.
    let mut model_atoms: Vec<Vec<Atom>> = Vec::new();
    let mut serial: u32 = 1;

    for (mi, model) in pdb.models().enumerate() {
        // Model selection: keep all models for now; filter later.
        let mut raw: Vec<Atom> = Vec::new();

        for chain in model.chains() {
            let chain_id = chain.id().to_string();
            for residue in chain.residues() {
                let res_id = ResidueId {
                    chain: chain_id.clone(),
                    name: residue.name().unwrap_or("").to_string(),
                    seq: residue.serial_number() as i32,
                    icode: residue
                        .insertion_code()
                        .and_then(|s| s.chars().next())
                        .unwrap_or(' '),
                };

                // Group atoms by atom name across conformers so we can pick
                // the highest-occupancy altLoc per (residue, atom name).
                // Key: atom name → list of (occupancy, altloc, atom).
                let mut by_name: BTreeMap<String, Vec<(f64, char, Atom)>> = BTreeMap::new();

                for conformer in residue.conformers() {
                    let altloc: char = conformer
                        .alternative_location()
                        .and_then(|s| s.chars().next())
                        .unwrap_or('\0');
                    for atom in conformer.atoms() {
                        let element = map_element(atom.element());
                        let a = Atom {
                            position: [atom.x(), atom.y(), atom.z()],
                            radius: 0.0,
                            element,
                            serial,
                            name: atom.name().to_string(),
                            residue: Some(res_id.clone()),
                            hetero: atom.hetero(),
                            occupancy: atom.occupancy(),
                            altloc,
                        };
                        serial += 1;
                        by_name
                            .entry(a.name.clone())
                            .or_default()
                            .push((a.occupancy, altloc, a));
                    }
                }

                for (_name, mut cands) in by_name {
                    match opts.altloc_policy {
                        AltlocPolicy::All => {
                            for (_, _, a) in cands {
                                raw.push(a);
                            }
                        }
                        AltlocPolicy::HighestOccupancy => {
                            cands.sort_by(|a, b| {
                                b.0.partial_cmp(&a.0)
                                    .unwrap_or(std::cmp::Ordering::Equal)
                                    .then_with(|| a.1.cmp(&b.1))
                            });
                            let keep = cands.len();
                            if let Some((_, _, a)) = cands.into_iter().next() {
                                raw.push(a);
                            }
                            report.altloc_dropped += keep.saturating_sub(1);
                        }
                        AltlocPolicy::First => {
                            // Prefer blank altLoc, then 'A', then first seen.
                            cands.sort_by_key(|(_, alt, _)| match alt {
                                '\0' | ' ' => 0u8,
                                'A' => 1,
                                _ => 2,
                            });
                            let keep = cands.len();
                            if let Some((_, _, a)) = cands.into_iter().next() {
                                raw.push(a);
                            }
                            report.altloc_dropped += keep.saturating_sub(1);
                        }
                    }
                }
            }
        }

        model_atoms.push(raw);
        let _ = mi;
    }

    // Model selection.
    let selected: Vec<Atom> = match opts.model {
        ModelSelection::First => model_atoms.into_iter().next().unwrap_or_default(),
        ModelSelection::Index(i) => model_atoms.into_iter().nth(i).unwrap_or_default(),
        ModelSelection::All => model_atoms.into_iter().flatten().collect(),
    };

    if selected.is_empty() {
        return Err(SesError::EmptySelection);
    }

    let kept = apply_selection(selected, opts, &mut report);
    if kept.is_empty() {
        return Err(SesError::EmptySelection);
    }
    report.atoms_kept = kept.len();

    let n = kept.len();
    let mut structure = Structure::from_atoms(
        kept,
        SourceInfo {
            path: path.display().to_string(),
            format: fmt,
            model_count,
        },
    );
    structure.set_models(Some(0..n));
    Ok((structure, report))
}

/// Map a `pdbtbx` element to our [`Element`].
fn map_element(e: Option<&pdbtbx::Element>) -> Option<Element> {
    let e = e?;
    // pdbtbx::Element has a symbol() method returning &str.
    Element::from_symbol(e.symbol())
}

#[cfg(test)]
mod tests {
    use super::*;

    const PDB_TWO_ATOMS: &str = "\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CB  ALA A   1       1.500   0.000   0.000  1.00  0.00           C
END
";

    #[test]
    fn parse_pdb_two_atoms() {
        let dir = std::env::temp_dir().join("molses_pdb_test");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("two.pdb");
        std::fs::write(&p, PDB_TWO_ATOMS).unwrap();
        let (s, r) = parse_pdb_like(&p, &ParseOptions::default()).unwrap();
        assert_eq!(s.len(), 2);
        assert_eq!(r.atoms_read, 2);
        assert_eq!(s.atoms()[0].element, Some(Element::C));
        assert!((s.atoms()[1].position[0] - 1.5).abs() < 1e-9);
    }

    #[test]
    fn parse_pdb_altloc_highest_occupancy() {
        let dir = std::env::temp_dir().join("molses_pdb_alt");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("alt.pdb");
        let src = "\
ATOM      1  CA AALA A   1       0.000   0.000   0.000  0.30  0.00           C
ATOM      2  CA BALA A   1       0.100   0.000   0.000  0.70  0.00           C
END
";
        std::fs::write(&p, src).unwrap();
        let (s, r) = parse_pdb_like(&p, &ParseOptions::default()).unwrap();
        assert_eq!(s.len(), 1);
        assert_eq!(r.altloc_dropped, 1);
        // Should keep the B altLoc (occupancy 0.7) at x=0.1.
        assert!((s.atoms()[0].position[0] - 0.1).abs() < 1e-9);
    }
}
