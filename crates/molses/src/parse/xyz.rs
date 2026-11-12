//! XYZ format parser (in-house).
//!
//! Classic XYZ: count line, comment line, then `Element x y z` per atom.
//! Multi-frame XYZ (concatenated blocks) becomes a multi-model structure.

use std::path::Path;

use super::{ParseOptions, ParseReport, apply_selection, parse_xyz_line};
use crate::error::{Result, SesError};
use crate::structure::{Atom, SourceInfo, Structure};

pub(crate) fn parse_xyz(path: &Path, opts: &ParseOptions) -> Result<(Structure, ParseReport)> {
    let text = std::fs::read_to_string(path)?;
    let mut report = ParseReport::default();
    let mut all_atoms: Vec<Atom> = Vec::new();
    let mut models: Vec<std::ops::Range<usize>> = Vec::new();

    let mut lines = text.lines().enumerate().peekable();
    let mut serial: u32 = 1;
    let mut model_count = 0usize;

    while lines.peek().is_some() {
        // Skip blank lines between frames.
        while let Some((_, l)) = lines.peek() {
            if l.trim().is_empty() {
                lines.next();
            } else {
                break;
            }
        }
        if lines.peek().is_none() {
            break;
        }

        let (count_line_no, count_line) = lines.next().unwrap();
        let count: usize = count_line.trim().parse().map_err(|_| {
            SesError::parse(
                path,
                Some(count_line_no + 1),
                format!("expected atom count, got `{}`", count_line.trim()),
            )
        })?;

        // Comment line (required by the format; tolerate missing).
        if lines.peek().is_some() {
            lines.next();
        }

        let start = all_atoms.len();
        for _ in 0..count {
            let (lineno, line) = lines.next().ok_or_else(|| {
                SesError::parse(
                    path,
                    Some(count_line_no + 1),
                    "unexpected end of file in atom block",
                )
            })?;
            if line.trim().is_empty() {
                return Err(SesError::parse(
                    path,
                    Some(lineno + 1),
                    "blank line in atom block",
                ));
            }
            let (element, position) = parse_xyz_line(line, path, lineno + 1)?;
            let mut atom = Atom::new(position, 0.0, element, serial);
            serial += 1;
            // Element-symbol fallback for name.
            atom.name = element.map(|e| e.symbol()).unwrap_or_default();
            all_atoms.push(atom);
        }
        models.push(start..all_atoms.len());
        model_count += 1;
    }

    if all_atoms.is_empty() {
        return Err(SesError::parse(path, None, "no atoms found"));
    }

    let kept = apply_selection(all_atoms, opts, &mut report);
    report.atoms_kept = kept.len();
    if kept.is_empty() {
        return Err(SesError::EmptySelection);
    }

    // Selection filters may drop atoms from any frame; for v1 we flatten the
    // kept atoms into a single model (multi-model selection is a stretch goal).
    let n = kept.len();
    let mut structure = Structure::from_atoms(
        kept,
        SourceInfo {
            path: path.display().to_string(),
            format: "xyz".into(),
            model_count,
        },
    );
    structure.set_models(Some(0..n));

    Ok((structure, report))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_single_frame() {
        let dir = std::env::temp_dir().join("molses_xyz_test");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("mol.xyz");
        std::fs::write(
            &p,
            "3\nwater\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n",
        )
        .unwrap();
        let (s, r) = parse_xyz(&p, &ParseOptions::default()).unwrap();
        // Default drops hydrogens → 1 atom kept.
        assert_eq!(r.atoms_read, 3);
        assert_eq!(s.len(), 1);
        assert_eq!(s.atoms()[0].element, Some(crate::structure::Element::O));
    }

    #[test]
    fn parse_keeps_hydrogens_when_asked() {
        let dir = std::env::temp_dir().join("molses_xyz_test2");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("mol.xyz");
        std::fs::write(&p, "3\n\nO 0 0 0\nH 1 0 0\nH 0 1 0\n").unwrap();
        let opts = ParseOptions {
            keep_hydrogens: true,
            ..Default::default()
        };
        let (s, _) = parse_xyz(&p, &opts).unwrap();
        assert_eq!(s.len(), 3);
    }

    #[test]
    fn parse_rejects_bad_count() {
        let dir = std::env::temp_dir().join("molses_xyz_test3");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("bad.xyz");
        std::fs::write(&p, "abc\n\n").unwrap();
        assert!(parse_xyz(&p, &ParseOptions::default()).is_err());
    }
}
