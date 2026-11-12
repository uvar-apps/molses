//! Structure file parsing.
//!
//! Format dispatch is by extension first, then content sniffing. PDB and mmCIF
//! go through `pdbtbx`; XYZ and XYZR are parsed in-house.

mod pdb;
mod xyz;
mod xyzr;

use std::path::Path;

use crate::error::{Result, SesError};
use crate::structure::{Atom, Element, Structure};

/// Selection policy applied while building a [`Structure`].
#[derive(Debug, Clone)]
pub struct ParseOptions {
    /// Include hydrogen atoms (default: `false`).
    pub keep_hydrogens: bool,
    /// Include waters (HOH/WAT residues) (default: `false`).
    pub include_waters: bool,
    /// Include HETATM records (default: `true`).
    pub include_hetatm: bool,
    /// altLoc policy: 0 = highest occupancy per residue, 1 = first/blank+A, 2 = all.
    pub altloc_policy: AltlocPolicy,
    /// Model selection: `None` = first model, `Some(n)` = 0-based model index,
    /// `All` handled by [`crate::structure::Structure::models`].
    pub model: ModelSelection,
    /// Chain IDs to keep (empty = all).
    pub chains: Vec<String>,
    /// Fail on unknown elements instead of falling back.
    pub strict_elements: bool,
}

/// How to treat alternate location indicators.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AltlocPolicy {
    /// Keep the highest-occupancy altLoc per residue (default).
    #[default]
    HighestOccupancy,
    /// Keep only blank altLoc and altLoc `A`.
    First,
    /// Keep every altLoc as a separate atom.
    All,
}

/// Which model(s) to keep.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ModelSelection {
    /// Keep only the first model (default).
    #[default]
    First,
    /// Keep model `n` (0-based).
    Index(usize),
    /// Keep all models.
    All,
}

impl Default for ParseOptions {
    fn default() -> Self {
        Self {
            keep_hydrogens: false,
            include_waters: false,
            include_hetatm: true,
            altloc_policy: AltlocPolicy::HighestOccupancy,
            model: ModelSelection::First,
            chains: Vec::new(),
            strict_elements: false,
        }
    }
}

/// Statistics gathered while parsing.
#[derive(Debug, Clone, Default)]
pub struct ParseReport {
    /// Atoms read from the file before selection.
    pub atoms_read: usize,
    /// Atoms kept after selection.
    pub atoms_kept: usize,
    /// Atoms dropped as hydrogens.
    pub hydrogens_dropped: usize,
    /// Atoms dropped as waters.
    pub waters_dropped: usize,
    /// Atoms dropped as HETATM (when disabled).
    pub hetatm_dropped: usize,
    /// Atoms dropped by altLoc policy.
    pub altloc_dropped: usize,
    /// Atoms dropped by chain filter.
    pub chain_dropped: usize,
    /// Duplicate centers merged (largest radius kept).
    pub duplicates_merged: usize,
    /// Spheres fully contained in another sphere, dropped.
    pub nested_dropped: usize,
    /// Elements resolved to the fallback radius.
    pub fallback_radii: usize,
}

/// Detect the structure format from path extension and file head.
pub fn detect_format(path: &Path) -> Result<String> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    // Handle `.pdb.gz` etc.
    let ext = if ext == "gz" {
        path.file_stem()
            .and_then(|s| Path::new(s).extension())
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_ascii_lowercase()
    } else {
        ext
    };
    match ext.as_str() {
        "pdb" | "ent" | "brk" => return Ok("pdb".into()),
        "cif" | "mmcif" => return Ok("mmcif".into()),
        "xyz" => return Ok("xyz".into()),
        "xyzr" | "sph" => return Ok("xyzr".into()),
        _ => {}
    }
    // Sniff content.
    let head = read_head(path, 4096)?;
    if head
        .lines()
        .any(|l| l.starts_with("ATOM  ") || l.starts_with("HETATM") || l.starts_with("MODEL "))
    {
        return Ok("pdb".into());
    }
    if head.starts_with("data_") || head.contains("_atom_site.") {
        return Ok("mmcif".into());
    }
    // XYZ: first non-empty line is an integer count.
    if let Some(first) = head.lines().find(|l| !l.trim().is_empty()) {
        if first.trim().parse::<usize>().is_ok() {
            return Ok("xyz".into());
        }
    }
    Err(SesError::UnsupportedFormat(format!(
        "cannot detect format of {}",
        path.display()
    )))
}

fn read_head(path: &Path, n: usize) -> Result<String> {
    use std::io::Read;
    let mut file = std::fs::File::open(path)?;
    let mut buf = vec![0u8; n];
    let read = file.read(&mut buf)?;
    buf.truncate(read);
    Ok(String::from_utf8_lossy(&buf).into_owned())
}

/// Parse a structure file with default options.
pub fn parse_file(path: &Path) -> Result<(Structure, ParseReport)> {
    parse_file_with(path, &ParseOptions::default())
}

/// Parse a structure file with explicit selection options.
pub fn parse_file_with(path: &Path, opts: &ParseOptions) -> Result<(Structure, ParseReport)> {
    let fmt = detect_format(path)?;
    match fmt.as_str() {
        "pdb" | "mmcif" => pdb::parse_pdb_like(path, opts),
        "xyz" => xyz::parse_xyz(path, opts),
        "xyzr" => xyzr::parse_xyzr(path, opts),
        other => Err(SesError::UnsupportedFormat(other.into())),
    }
}

/// Assign van der Waals radii to every atom using the given table.
pub fn assign_radii(
    structure: &mut Structure,
    radii: &crate::radii::RadiiSet,
    strict: bool,
    report: &mut ParseReport,
) -> Result<()> {
    for atom in structure.atoms_mut() {
        match radii.resolve(atom.element, strict) {
            Ok(r) => {
                if atom.radius <= 0.0 {
                    atom.radius = r;
                }
            }
            Err(SesError::UnknownElement { .. }) if !strict => {
                atom.radius = crate::radii::FALLBACK_RADIUS;
                report.fallback_radii += 1;
            }
            Err(e) => return Err(e),
        }
    }
    Ok(())
}

/// Post-parse cleanup: merge duplicate centers, drop fully nested spheres.
///
/// MSMS merges identical centers (keeping the largest radius) and crashes on
/// nested spheres; we merge duplicates and drop nested ones instead.
/// Uses a uniform-grid hash so the cost is O(n) for typical structures.
pub fn deduplicate(structure: &mut Structure, report: &mut ParseReport) {
    let atoms = structure.atoms_mut();
    let n = atoms.len();
    if n == 0 {
        return;
    }

    // --- Pass 1: merge identical centers (within 1e-6 Å) keeping largest radius ---
    // Hash by quantized position; collisions in the same bucket are candidates.
    let quant = 1e-6f64;
    let mut buckets: std::collections::HashMap<(i64, i64, i64), Vec<usize>> =
        std::collections::HashMap::new();
    for (i, a) in atoms.iter().enumerate() {
        let key = (
            (a.position[0] / quant).round() as i64,
            (a.position[1] / quant).round() as i64,
            (a.position[2] / quant).round() as i64,
        );
        buckets.entry(key).or_default().push(i);
    }
    let mut keep = vec![true; n];
    for group in buckets.values() {
        if group.len() < 2 {
            continue;
        }
        for (gi, &i) in group.iter().enumerate() {
            if !keep[i] {
                continue;
            }
            for &j in &group[(gi + 1)..] {
                if !keep[j] {
                    continue;
                }
                if atoms[i].distance_squared_to(atoms[j].position) < 1e-12 {
                    if atoms[j].radius > atoms[i].radius {
                        atoms[i].radius = atoms[j].radius;
                    }
                    keep[j] = false;
                    report.duplicates_merged += 1;
                }
            }
        }
    }

    // --- Pass 2: drop spheres fully contained in another kept sphere ---
    // If sphere i is inside sphere j then d(c_i, c_j) <= r_j - r_i <= r_j_max,
    // so we only need to search within max_radius of c_i.
    let max_r = atoms.iter().map(|a| a.radius).fold(0.0f64, f64::max);
    if max_r > 0.0 {
        let cell = (max_r * 2.0).max(1e-3);
        let mut grid: std::collections::HashMap<(i64, i64, i64), Vec<usize>> =
            std::collections::HashMap::new();
        for (i, a) in atoms.iter().enumerate() {
            if !keep[i] {
                continue;
            }
            let key = (
                (a.position[0] / cell).floor() as i64,
                (a.position[1] / cell).floor() as i64,
                (a.position[2] / cell).floor() as i64,
            );
            grid.entry(key).or_default().push(i);
        }
        let ring = (max_r / cell).ceil() as i64;
        for i in 0..n {
            if !keep[i] {
                continue;
            }
            let a = &atoms[i];
            let base = (
                (a.position[0] / cell).floor() as i64,
                (a.position[1] / cell).floor() as i64,
                (a.position[2] / cell).floor() as i64,
            );
            'outer: for dx in -ring..=ring {
                for dy in -ring..=ring {
                    for dz in -ring..=ring {
                        let key = (base.0 + dx, base.1 + dy, base.2 + dz);
                        let Some(group) = grid.get(&key) else {
                            continue;
                        };
                        for &j in group {
                            if i == j || !keep[j] {
                                continue;
                            }
                            let d2 = atoms[i].distance_squared_to(atoms[j].position);
                            if d2.sqrt() + atoms[i].radius <= atoms[j].radius + 1e-9 {
                                keep[i] = false;
                                report.nested_dropped += 1;
                                break 'outer;
                            }
                        }
                    }
                }
            }
        }
    }

    let mut kept = Vec::with_capacity(n);
    for (i, atom) in atoms.drain(..).enumerate() {
        if keep[i] {
            kept.push(atom);
        }
    }
    *atoms = kept;
    report.atoms_kept = atoms.len();
}

/// Apply selection filters (hydrogens, waters, HETATM, chains).
pub(crate) fn apply_selection(
    atoms: Vec<Atom>,
    opts: &ParseOptions,
    report: &mut ParseReport,
) -> Vec<Atom> {
    let mut out = Vec::with_capacity(atoms.len());
    for a in atoms {
        report.atoms_read += 1;
        if !opts.keep_hydrogens && a.element.map(|e| e.is_hydrogen()).unwrap_or(false) {
            report.hydrogens_dropped += 1;
            continue;
        }
        if !opts.include_waters {
            if let Some(res) = &a.residue {
                if matches!(res.name.as_str(), "HOH" | "WAT" | "DOD" | "H2O") {
                    report.waters_dropped += 1;
                    continue;
                }
            }
        }
        if !opts.include_hetatm && a.hetero {
            report.hetatm_dropped += 1;
            continue;
        }
        if !opts.chains.is_empty() {
            let chain = a.residue.as_ref().map(|r| r.chain.as_str()).unwrap_or("");
            if !opts.chains.iter().any(|c| c == chain) {
                report.chain_dropped += 1;
                continue;
            }
        }
        out.push(a);
    }
    out
}

/// Parse a free-format `element x y z` line (XYZ core).
pub(crate) fn parse_xyz_line(
    line: &str,
    path: &Path,
    lineno: usize,
) -> Result<(Option<Element>, [f64; 3])> {
    let mut it = line.split_whitespace();
    let sym = it
        .next()
        .ok_or_else(|| SesError::parse(path, Some(lineno), "expected element symbol"))?;
    let mut pos = [0.0f64; 3];
    for (k, p) in pos.iter_mut().enumerate() {
        let tok = it.next().ok_or_else(|| {
            SesError::parse(path, Some(lineno), format!("expected coordinate {}", k + 1))
        })?;
        *p = tok
            .parse()
            .map_err(|_| SesError::parse(path, Some(lineno), format!("bad coordinate `{tok}`")))?;
    }
    Ok((Element::from_symbol(sym), pos))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detect_xyz_and_pdb() {
        let dir = std::env::temp_dir().join("molses_fmt_test");
        std::fs::create_dir_all(&dir).unwrap();
        let xyz = dir.join("a.xyz");
        std::fs::write(&xyz, "2\ncomment\nC 0 0 0\nO 1 0 0\n").unwrap();
        assert_eq!(detect_format(&xyz).unwrap(), "xyz");
        let pdb = dir.join("a.pdb");
        std::fs::write(
            &pdb,
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n",
        )
        .unwrap();
        assert_eq!(detect_format(&pdb).unwrap(), "pdb");
        let cif = dir.join("a.cif");
        std::fs::write(&cif, "data_test\nloop_\n_atom_site.id\n").unwrap();
        assert_eq!(detect_format(&cif).unwrap(), "mmcif");
    }
}
