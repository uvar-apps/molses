//! XYZR sphere-file parser (MSMS interchange format).
//!
//! One sphere per line: `x y z r` in Ångström, free format. `#` comments.

use std::path::Path;

use super::{ParseOptions, ParseReport};
use crate::error::{Result, SesError};
use crate::structure::{Atom, SourceInfo, Structure};

pub(crate) fn parse_xyzr(path: &Path, _opts: &ParseOptions) -> Result<(Structure, ParseReport)> {
    let text = std::fs::read_to_string(path)?;
    let mut report = ParseReport::default();
    let mut atoms = Vec::new();
    let mut serial: u32 = 1;

    for (i, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(SesError::parse(
                path,
                Some(i + 1),
                format!("expected `x y z r`, got `{line}`"),
            ));
        }
        let mut vals = [0.0f64; 4];
        for (k, tok) in parts.iter().take(4).enumerate() {
            vals[k] = tok
                .parse()
                .map_err(|_| SesError::parse(path, Some(i + 1), format!("bad number `{tok}`")))?;
        }
        if !vals[3].is_finite() || vals[3] <= 0.0 {
            return Err(SesError::parse(
                path,
                Some(i + 1),
                "radius must be positive",
            ));
        }
        // XYZR has no element; leave it unknown so the radius is used as-is.
        let mut atom = Atom::new([vals[0], vals[1], vals[2]], vals[3], None, serial);
        atom.name = format!("S{serial}");
        atoms.push(atom);
        serial += 1;
        report.atoms_read += 1;
    }

    if atoms.is_empty() {
        return Err(SesError::parse(path, None, "no spheres found"));
    }
    report.atoms_kept = atoms.len();
    let n = atoms.len();
    Ok((
        Structure::from_atoms(
            atoms,
            SourceInfo {
                path: path.display().to_string(),
                format: "xyzr".into(),
                model_count: 1,
            },
        ),
        report,
    ))
    .map(|(mut s, r)| {
        // Ensure model range is valid.
        s.set_models(Some(0..n));
        (s, r)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_xyzr_basic() {
        let dir = std::env::temp_dir().join("molses_xyzr_test");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("spheres.xyzr");
        std::fs::write(&p, "# comment\n0 0 0 1.7\n1.5 0 0 1.7\n\n").unwrap();
        let (s, r) = parse_xyzr(&p, &ParseOptions::default()).unwrap();
        assert_eq!(s.len(), 2);
        assert_eq!(r.atoms_read, 2);
        assert_eq!(s.atoms()[0].radius, 1.7);
    }

    #[test]
    fn parse_xyzr_rejects_zero_radius() {
        let dir = std::env::temp_dir().join("molses_xyzr_test2");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("bad.xyzr");
        std::fs::write(&p, "0 0 0 0\n").unwrap();
        assert!(parse_xyzr(&p, &ParseOptions::default()).is_err());
    }
}
