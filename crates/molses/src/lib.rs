//! # molses — Solvent Excluded Surface (Connolly surface) generation
//!
//! `molses` computes a resolution-controllable, watertight triangle mesh of the
//! **Solvent Excluded Surface** (SES) of a molecule — the surface traced by a
//! spherical solvent probe rolling over the van der Waals spheres of its atoms —
//! and exports it as a Wavefront OBJ file with per-vertex normals.
//!
//! The crate is a pure-Rust library; the companion `molses-cli` binary is a
//! thin wrapper around this API.
//!
//! ## Quick start
//!
//! ```no_run
//! use molses::{parse, radii::RadiiSet, SesOptions};
//! # fn main() -> molses::Result<()> {
//! // Parse the structure (format auto-detected from extension/content).
//! let path = std::path::Path::new("1crn.pdb");
//! let (mut structure, mut report) = parse::parse_file(path)?;
//!
//! // Assign van der Waals radii and clean up degenerate atoms.
//! parse::assign_radii(&mut structure, &RadiiSet::Bondi, false, &mut report)?;
//! parse::deduplicate(&mut structure, &mut report);
//!
//! // Compute the SES mesh.
//! let opts = SesOptions {
//!     probe_radius: 1.4,
//!     grid_spacing: 0.35,
//!     ..SesOptions::default()
//! };
//! let mesh = molses::ses(&structure, &opts)?;
//! println!("{} triangles, area {:.1} Å²", mesh.triangle_count(), mesh.area());
//! # Ok(()) }
//! ```
//!
//! ## Writing OBJ output
//!
//! ```no_run
//! # fn main() -> molses::Result<()> {
//! # use molses::{Mesh, io::WriteMeta};
//! # let mesh = Mesh::default();
//! let meta = WriteMeta {
//!     tool: format!("molses {}", molses::version()),
//!     source: "1crn.pdb".into(),
//!     probe_radius: 1.4,
//!     grid_spacing: 0.35,
//!     atoms: 327,
//!     deterministic: true,
//! };
//! let mut f = std::io::BufWriter::new(std::fs::File::create("out.obj")?);
//! molses::io::obj::write_obj(&mut f, &mesh, &meta)?;
//! # Ok(()) }
//! ```
//!
//! ## Pipeline overview
//!
//! The main entry point [`ses`] runs this pipeline:
//!
//! 1. **Parse** — PDB / mmCIF / XYZ / XYZR → [`Structure`] (via [`parse`]).
//! 2. **Radii** — assign Bondi (1964) or custom van der Waals radii.
//! 3. **Field** — paint the inflated union `K = ⋃ B(c_i, r_i+r_p)` on a grid,
//!    compute an exact Euclidean distance transform, and assemble the signed
//!    distance `f_ses = f_sas + r_p`.
//! 4. **Mesh** — extract the zero level set with marching tetrahedra (6-tet
//!    decomposition, crack-free, watertight welding).
//! 5. **Normals** — field-gradient normals (exact SDF normals, smooth shading).
//! 6. **Filter** — connected-component policy (outer shell / all / cavities).
//! 7. **Output** — [`Mesh`] returned to the caller; write OBJ via [`io`].
//!
//! ## Key types
//!
//! | Type | Purpose |
//! |---|---|
//! | [`Structure`] | Parsed molecule: list of [`Atom`]s |
//! | [`SesOptions`] | Probe radius, grid spacing, policies |
//! | [`Mesh`] | Positions, normals, triangle indices |
//! | [`SesError`] | Typed errors (never panics on input) |
//!
//! ## Mathematical definition
//!
//! For atoms `(c_i, r_i)` and probe radius `r_p`, let `K = ⋃ B(c_i, r_i + r_p)`
//! be the SAS cavity and `f_sas` its signed distance function. The SES is the
//! zero level set of `f_ses(p) = f_sas(p) + r_p` (Quan & Stamm 2016, Eq. 3.4).
//!
//! ## Feature flags
//!
//! - `parallel` (default) — Rayon data parallelism.
//! - `simd` (default) — portable SIMD via the `wide` crate.
//! - `serde` — serialization of core types.
//! - `mmap` — memory-mapped fields for very large grids.
//! - `gzip` — transparent reading of `.gz` inputs.

#![warn(missing_docs)]

pub mod error;
pub mod field;
pub mod io;
pub mod mesh;
pub mod parse;
pub mod radii;
pub mod structure;

pub use error::{Result, SesError};
pub use mesh::Mesh;
pub use structure::{Atom, Element, Structure};

/// Which surface to compute.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SurfaceKind {
    /// Solvent Excluded Surface (Connolly / molecular surface). Default.
    #[default]
    Ses,
    /// Solvent Accessible Surface (locus of probe centers).
    Sas,
    /// van der Waals surface (probe radius treated as 0).
    Vdw,
}

/// How to handle connected components of the mesh.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ComponentPolicy {
    /// Keep only the outer shell. Default.
    #[default]
    OuterOnly,
    /// Keep every closed component (outer shell + cavity walls).
    All,
    /// Drop enclosed components whose volume is below the threshold (Å³).
    NoCavities {
        /// Minimum enclosed volume in Å³ for a component to be kept.
        min_volume: f64,
    },
}

/// Options controlling SES generation.
#[derive(Debug, Clone)]
pub struct SesOptions {
    /// Probe radius `r_p` in Å (default 1.4, water).
    pub probe_radius: f64,
    /// Grid spacing `h` in Å (level of detail; default 0.35).
    pub grid_spacing: f64,
    /// Which surface to compute.
    pub surface: SurfaceKind,
    /// Connected-component policy.
    pub components: ComponentPolicy,
    /// Laplacian smoothing iterations (0 = none).
    pub smoothing: u32,
    /// van der Waals radii table.
    pub radii: radii::RadiiSet,
    /// Worker threads (`None` = Rayon default).
    pub threads: Option<usize>,
    /// Soft memory budget in bytes; the engine may raise `h` to fit.
    pub max_memory_bytes: Option<u64>,
    /// Compute per-vertex nearest-atom indices.
    pub nearest_atom: bool,
}

impl Default for SesOptions {
    fn default() -> Self {
        Self {
            probe_radius: 1.4,
            grid_spacing: 0.35,
            surface: SurfaceKind::Ses,
            components: ComponentPolicy::OuterOnly,
            smoothing: 0,
            radii: radii::RadiiSet::Bondi,
            threads: None,
            max_memory_bytes: None,
            nearest_atom: true,
        }
    }
}

impl SesOptions {
    /// Validate parameters, returning a structured error on bad values.
    pub fn validate(&self) -> Result<()> {
        if !self.probe_radius.is_finite() || self.probe_radius <= 0.0 {
            return Err(SesError::invalid_param(
                "probe_radius",
                self.probe_radius,
                "must be finite and > 0",
            ));
        }
        if !self.grid_spacing.is_finite() || self.grid_spacing <= 0.0 {
            return Err(SesError::invalid_param(
                "grid_spacing",
                self.grid_spacing,
                "must be finite and > 0",
            ));
        }
        if self.probe_radius > 10.0 {
            return Err(SesError::invalid_param(
                "probe_radius",
                self.probe_radius,
                "must be ≤ 10 Å",
            ));
        }
        Ok(())
    }
}

/// Compute the requested surface mesh for `structure`.
///
/// This is the primary entry point. It runs the full pipeline: field
/// construction, isosurface extraction, normals, and component filtering.
pub fn ses(structure: &Structure, opts: &SesOptions) -> Result<Mesh> {
    opts.validate()?;
    if structure.is_empty() {
        return Err(SesError::EmptySelection);
    }

    // Effective probe radius depends on the surface kind.
    let probe = match opts.surface {
        SurfaceKind::Ses | SurfaceKind::Sas => opts.probe_radius,
        SurfaceKind::Vdw => 0.0,
    };

    // --- Domain sizing (DESIGN §5.2) ---
    let (lo, hi) = structure.bounding_box();
    let r_max = structure
        .atoms()
        .iter()
        .map(|a| a.radius)
        .fold(0.0f64, f64::max);
    let pad = 2.0 * probe + r_max + 3.0 * opts.grid_spacing;
    let origin = [lo[0] - pad, lo[1] - pad, lo[2] - pad];
    let extent = [
        hi[0] + pad - origin[0],
        hi[1] + pad - origin[1],
        hi[2] + pad - origin[2],
    ];
    let h = opts.grid_spacing;
    let dims = [
        (extent[0] / h).ceil() as usize + 1,
        (extent[1] / h).ceil() as usize + 1,
        (extent[2] / h).ceil() as usize + 1,
    ];

    // Memory budget check (F-PAR-4 / DESIGN §5.3).
    let needed = field::estimate_bytes(dims);
    if let Some(limit) = opts.max_memory_bytes {
        if needed > limit {
            return Err(SesError::ResourceLimit {
                needed: format!("{needed} bytes"),
                limit: format!("{limit} bytes"),
            });
        }
    }

    // --- Field construction: EDT-based signed distance with sub-voxel correction ---
    // EDT gives correct topology (watertight, Euler). Exact edge roots below
    // remove quantization jitter from vertex positions.
    let mut fld = field::ScalarField::new(dims, origin, h)?;
    let paint_probe = match opts.surface {
        SurfaceKind::Vdw => 0.0,
        _ => probe,
    };
    let painted = field::paint::paint(structure.atoms(), paint_probe, &fld);
    let edt_sq = field::edt::edt_squared(&painted.inside, dims);
    field::paint::assemble_signed_distance(&painted, &edt_sq, h, &mut fld);

    // Offset by probe radius for SES (f_ses = f_sas + r_p).
    if matches!(opts.surface, SurfaceKind::Ses) {
        let data = fld.data_mut();
        for v in data.iter_mut() {
            *v += probe as f32;
        }
    }

    // --- Mesh extraction ---
    // Exact edge roots need the atoms (contact spheres). Reentrant toric
    // roots stay off by default (MOLSES_REENTRANT=1) — they need more care.
    let root_atoms = match opts.surface {
        SurfaceKind::Ses | SurfaceKind::Vdw => Some(structure.atoms()),
        SurfaceKind::Sas => None,
    };
    let mut mesh = mesh::polygonize::polygonize(&fld, root_atoms, probe);
    // Put every surface vertex exactly on its analytic patch:
    // contact → VdW sphere; reentrant → torus at r_p from an intersection
    // circle. Corrections are bounded by ~h (we only move vertices already
    // near the target), so topology and area stay put.
    let mut patches = vec![mesh::normals::Patch::None; mesh.positions.len()];
    if matches!(opts.surface, SurfaceKind::Ses | SurfaceKind::Vdw) {
        project_exact_patches(&mut mesh, structure.atoms(), probe, h, &fld, &mut patches);
    }
    // Angle-weighted geometric normals, then blend toward exact radial
    // shading on contact patches (removes polyhedral faceting on snapped
    // vertices without analytic/reentrant seam artifacts).
    mesh::normals::compute_normals(&mut mesh, Some(&fld));
    if matches!(opts.surface, SurfaceKind::Ses | SurfaceKind::Vdw) {
        mesh::normals::apply_patch_normals(&mut mesh, structure.atoms(), probe, &patches);
        mesh::normals::repair_flipped_normals(&mut mesh);
    }
    // Smooth ONLY reentrant vertex positions (Taubin, non-shrinking).
    // Contact-patch vertices are locked to their VdW spheres.
    if opts.smoothing > 0 {
        // Identify contact-patch vertices (locked) vs reentrant (free).
        let locked: Vec<bool> = mesh
            .positions
            .iter()
            .map(|p| {
                structure.atoms().iter().any(|a| {
                    let dx = p[0] as f64 - a.position[0];
                    let dy = p[1] as f64 - a.position[1];
                    let dz = p[2] as f64 - a.position[2];
                    let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                    (dist - a.radius).abs() < h * 0.3
                })
            })
            .collect();
        // Build neighbor list.
        use std::collections::HashMap;
        let n_v = mesh.positions.len();
        let mut neighbors: Vec<Vec<u32>> = vec![Vec::new(); n_v];
        for t in &mesh.indices {
            for k in 0..3 {
                let a = t[k] as usize;
                for j in [t[(k + 1) % 3], t[(k + 2) % 3]] {
                    if !neighbors[a].contains(&j) {
                        neighbors[a].push(j);
                    }
                }
            }
        }
        let lambda = 0.5f32;
        let mu = -0.53f32;
        for _ in 0..opts.smoothing {
            for (w, pass) in [(lambda, 0u8), (mu, 1u8)] {
                let current = mesh.positions.clone();
                for i in 0..n_v {
                    if locked[i] {
                        continue;
                    }
                    let nb = &neighbors[i];
                    if nb.is_empty() {
                        continue;
                    }
                    let mut avg = [0.0f32; 3];
                    for &j in nb {
                        let q = current[j as usize];
                        avg[0] += q[0];
                        avg[1] += q[1];
                        avg[2] += q[2];
                    }
                    let inv = 1.0 / nb.len() as f32;
                    avg = [avg[0] * inv, avg[1] * inv, avg[2] * inv];
                    mesh.positions[i] = [
                        current[i][0] + w * (avg[0] - current[i][0]),
                        current[i][1] + w * (avg[1] - current[i][1]),
                        current[i][2] + w * (avg[2] - current[i][2]),
                    ];
                }
                let _ = pass;
            }
        }
        let _ = HashMap::<u32, u32>::new();
    }

    // --- Component policy (F-SUR-5) ---
    match opts.components {
        ComponentPolicy::All => {} // keep everything
        ComponentPolicy::OuterOnly => {
            let comps = mesh::components::connected_components(&mesh);
            if comps.len() > 1 {
                mesh = mesh::components::filter_components(&mesh, &comps[..1]);
            }
        }
        ComponentPolicy::NoCavities { min_volume } => {
            let comps = mesh::components::connected_components(&mesh);
            let keep: Vec<mesh::components::Component> = comps
                .iter()
                .filter(|c| c.volume.abs() >= min_volume)
                .cloned()
                .collect();
            if keep.len() != comps.len() && !keep.is_empty() {
                mesh = mesh::components::filter_components(&mesh, &keep);
            }
        }
    }

    Ok(mesh)
}

/// Library version string.
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

type V3 = [f64; 3];

#[inline]
fn sub3(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
fn add3(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
fn scale3(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
fn dot3(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
fn len3(a: V3) -> f64 {
    dot3(a, a).sqrt()
}

/// Put every surface vertex exactly on the analytic patch it represents.
///
/// Two patch families occur. A vertex on a **contact** patch lies on the VdW
/// sphere of one atom. A vertex on a **concave** patch lies on the canal torus
/// of an atom pair: the set of points at distance `probe` from the circle a
/// probe traces while touching both atoms, taken on the crevice-facing half of
/// that circle. Correcting instead of re-rooting keeps the tessellation (and
/// therefore the topology) exactly as the mesher produced it.
///
/// A projection is accepted only when the matching probe centre sits outside
/// every other atom's probe-inflated ball. That one visibility rule is what
/// keeps 3-body pockets — whose surface belongs to no 2-atom torus — from
/// being flattened onto one, and it is also what makes a contact vertex under
/// a neighbour stay put.
///
/// Neighbourhoods come from a cell list, so cost is linear in vertices;
/// resolving hosts against every atom pair cost more than the field and the
/// mesher put together.
fn project_exact_patches(
    mesh: &mut Mesh,
    atoms: &[Atom],
    probe: f64,
    h: f64,
    fld: &field::ScalarField,
    patch: &mut [mesh::normals::Patch],
) {
    /// Circle traced by a probe touching both atoms of a pair.
    struct Circle {
        o: V3,
        r: f64,
        u: V3,
    }
    let mut circles: Vec<Circle> = Vec::new();
    // atom -> (partner atom, circle index)
    let mut links: Vec<Vec<(u32, u32)>> = vec![Vec::new(); atoms.len()];
    if probe > 0.0 {
        for i in 0..atoms.len() {
            for j in (i + 1)..atoms.len() {
                let ri = atoms[i].radius + probe;
                let rj = atoms[j].radius + probe;
                let e = sub3(atoms[j].position, atoms[i].position);
                let d = len3(e);
                if d < 1e-9 || d >= ri + rj || d <= (ri - rj).abs() {
                    continue;
                }
                let t = (d * d + ri * ri - rj * rj) / (2.0 * d);
                let rc2 = ri * ri - t * t;
                if rc2 <= 1e-12 {
                    continue;
                }
                let id = circles.len() as u32;
                circles.push(Circle {
                    o: add3(atoms[i].position, scale3(e, t / d)),
                    r: rc2.sqrt(),
                    u: scale3(e, 1.0 / d),
                });
                links[i].push((j as u32, id));
                links[j].push((i as u32, id));
            }
        }
    }

    let r_max = atoms.iter().map(|a| a.radius).fold(0.0f64, f64::max);
    // Covers any possible host pair and any atom that could hide a probe
    // centre belonging to this vertex.
    let reach = r_max + 2.0 * probe + 3.0 * h;
    let cells = field::celllist::CellList::new(atoms, reach, fld.origin(), fld.dims(), h);
    let slack = 0.05 * h;
    let mut near: Vec<u32> = Vec::with_capacity(64);
    let mut nd: Vec<f64> = Vec::with_capacity(64);
    let mut stamp = vec![0u32; atoms.len()];

    for (pi, p) in mesh.positions.iter_mut().enumerate() {
        let q = [p[0] as f64, p[1] as f64, p[2] as f64];
        let epoch = pi as u32 + 1;
        cells.query(q, reach, &mut near);
        nd.clear();
        let mut on_sphere = false;
        let mut host = usize::MAX;
        let mut host_err = f64::INFINITY;
        for &k in &near {
            stamp[k as usize] = epoch;
            let a = &atoms[k as usize];
            let d = len3(sub3(q, a.position));
            let err = (d - a.radius).abs();
            if err < 0.05 * h {
                on_sphere = true;
            }
            if err < host_err {
                host_err = err;
                host = k as usize;
            }
            nd.push(d);
        }

        // `true` when a probe centred at `pc` overlaps no local atom other
        // than the two it is meant to touch.
        let clear = |pc: V3, skip: usize, skip2: usize| -> bool {
            near.iter().all(|&k| {
                let ki = k as usize;
                ki == skip
                    || ki == skip2
                    || len3(sub3(pc, atoms[ki].position)) >= atoms[ki].radius + probe - slack
            })
        };

        // --- contact patch: onto the host sphere ---
        if host != usize::MAX && host_err < 1.5 * h {
            let a = &atoms[host];
            let v = sub3(q, a.position);
            let d = len3(v);
            let uv = scale3(v, 1.0 / d.max(1e-12));
            if d > 1e-12
                && clear(
                    add3(a.position, scale3(uv, a.radius + probe)),
                    host,
                    usize::MAX,
                )
            {
                let np = add3(a.position, scale3(uv, a.radius));
                *p = [np[0] as f32, np[1] as f32, np[2] as f32];
                patch[pi] = mesh::normals::Patch::Contact(host as u32);
                continue;
            }
        }
        // Sitting on a sphere but not projectable, or a VdW-only surface.
        if on_sphere || circles.is_empty() {
            continue;
        }

        // --- canal patch: onto the torus of the best visible pair ---
        // (distance from the circle, nearest circle point, hosts)
        type CanalHit = (f64, V3, usize, usize);
        let mut best: Option<CanalHit> = None;
        let mut best_err = 0.5 * h;
        for (t, &ka) in near.iter().enumerate() {
            // A host cannot be further from the surface point than its own
            // inflated radius; this prunes the pair enumeration hard.
            if nd[t] > atoms[ka as usize].radius + probe + h {
                continue;
            }
            for &(kb, ci) in &links[ka as usize] {
                if kb <= ka || stamp[kb as usize] != epoch {
                    continue;
                }
                let c = &circles[ci as usize];
                let w = sub3(q, c.o);
                let perp = sub3(w, scale3(c.u, dot3(w, c.u)));
                let rho2 = dot3(perp, perp);
                if rho2 >= c.r * c.r || rho2 < 1e-24 {
                    continue; // outer half of the torus carries no SES
                }
                let npt = add3(c.o, scale3(perp, c.r / rho2.sqrt()));
                let err = (len3(sub3(q, npt)) - probe).abs();
                if err < best_err && clear(npt, ka as usize, kb as usize) {
                    best_err = err;
                    best = Some((len3(sub3(q, npt)), npt, ka as usize, kb as usize));
                }
            }
        }
        if let Some((dn, npt, ai, aj)) = best {
            let np = add3(npt, scale3(sub3(q, npt), probe / dn));
            *p = [np[0] as f32, np[1] as f32, np[2] as f32];
            patch[pi] = mesh::normals::Patch::Canal(ai as u32, aj as u32);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn options_validate() {
        let ok = SesOptions::default();
        assert!(ok.validate().is_ok());
        let bad = SesOptions {
            probe_radius: -1.0,
            ..Default::default()
        };
        assert!(bad.validate().is_err());
        let bad2 = SesOptions {
            grid_spacing: 0.0,
            ..Default::default()
        };
        assert!(bad2.validate().is_err());
    }

    #[test]
    fn version_is_semverish() {
        assert!(version().contains('.'));
    }

    /// Regression: when `probe_radius / grid_spacing` is an integer, the EDT
    /// field carries exact zeros on a whole shell. Without tie-breaking the
    /// mesh gets holes. This test runs the full `ses()` pipeline at the
    /// previously failing ratio 1.4 / 0.35 = 4.0 and asserts watertightness.
    #[test]
    fn full_pipeline_watertight_at_integer_ratio() {
        // Two atoms close enough to form reentrant surface, at the
        // historically worst r_p/h = 4.0.
        let atoms = vec![
            Atom::new([0.0, 0.0, 0.0], 1.7, Some(Element::C), 1),
            Atom::new([2.5, 0.0, 0.0], 1.7, Some(Element::C), 2),
        ];
        let structure = Structure::from_atoms(
            atoms,
            structure::SourceInfo {
                path: "test".into(),
                format: "synthetic".into(),
                model_count: 1,
            },
        );
        let opts = SesOptions {
            probe_radius: 1.4,
            grid_spacing: 0.35, // 1.4 / 0.35 = 4.0 exactly
            ..SesOptions::default()
        };
        let mesh = ses(&structure, &opts).unwrap();
        assert!(mesh.triangle_count() > 0, "empty mesh");
        let report = mesh::topology::audit(&mesh);
        assert!(
            report.is_watertight(),
            "mesh has holes at integer r_p/h ratio: {report:?}"
        );
    }

    /// Same check at a non-integer ratio (the historically safe case).
    #[test]
    fn full_pipeline_watertight_at_fractional_ratio() {
        let atoms = vec![
            Atom::new([0.0, 0.0, 0.0], 1.7, Some(Element::C), 1),
            Atom::new([2.5, 0.0, 0.0], 1.7, Some(Element::C), 2),
        ];
        let structure = Structure::from_atoms(
            atoms,
            structure::SourceInfo {
                path: "test".into(),
                format: "synthetic".into(),
                model_count: 1,
            },
        );
        let opts = SesOptions {
            probe_radius: 1.4,
            grid_spacing: 0.5, // 1.4 / 0.5 = 2.8
            ..SesOptions::default()
        };
        let mesh = ses(&structure, &opts).unwrap();
        let report = mesh::topology::audit(&mesh);
        assert!(report.is_watertight(), "{report:?}");
    }
}

#[cfg(test)]
mod geometry_tests {
    use super::*;
    use crate::structure::Element;

    /// Two-atom SES analytic reference (Richards 1977, Connolly 1983).
    ///
    /// Atoms: R=1.7, d=3.4 (just touching), probe r_p=1.4.
    /// The SES toric reentrant patch dips to y=1.192 at x=0 (between atoms).
    /// Contact patches extend to theta_max=56.74° from the bond axis.
    /// Contact endpoint at (-0.768, 1.422) relative to atom 1 at (-1.7, 0).
    fn two_atom_mesh() -> Mesh {
        let atoms = vec![
            Atom::new([-1.7, 0.0, 0.0], 1.7, Some(Element::C), 1),
            Atom::new([1.7, 0.0, 0.0], 1.7, Some(Element::C), 2),
        ];
        let structure = Structure::from_atoms(
            atoms,
            structure::SourceInfo {
                path: "two-atom".into(),
                format: "synthetic".into(),
                model_count: 1,
            },
        );
        let opts = SesOptions {
            probe_radius: 1.4,
            grid_spacing: 0.15, // fine enough for analytic comparison
            ..SesOptions::default()
        };
        ses(&structure, &opts).unwrap()
    }

    /// The SES reentrant surface between two atoms dips to y=1.192 at x=0.
    /// If the mesh has no vertices near that point, the reentrant patch is
    /// missing or wrong.
    #[test]
    fn two_atom_reentrant_patch_exists() {
        let mesh = two_atom_mesh();
        // The toric reentrant inner point is at (0, 1.192, 0).
        // Check for mesh vertices within 0.15 Å of that point.
        let target = [0.0f32, 1.192, 0.0];
        let mut min_dist = f32::MAX;
        for p in &mesh.positions {
            let d = ((p[0] - target[0]).powi(2)
                + (p[1] - target[1]).powi(2)
                + (p[2] - target[2]).powi(2))
            .sqrt();
            min_dist = min_dist.min(d);
        }
        assert!(
            min_dist < 0.15,
            "No mesh vertex near toric reentrant point (0, 1.192, 0). \
             Closest vertex is {min_dist:.3} Å away. \
             The reentrant patch is missing or displaced. \
             This is the key SES geometry test."
        );
    }

    /// The SES must dip INTO the crevice between the atoms.
    /// The VdW surface at z=0 has y=0 at x=0 (atoms just touching).
    /// The SES reentrant surface is at y=1.192 (above the VdW contact).
    /// If the mesh at x=0, z=0 has max y < 1.0, the reentrant patch is too low.
    #[test]
    fn two_atom_reentrant_dips_correctly() {
        let mesh = two_atom_mesh();
        // Find vertices near x=0, z=0 (the bond axis cross-section).
        let mut max_y = f32::MIN;
        let mut count = 0;
        for p in &mesh.positions {
            if p[0].abs() < 0.2 && p[2].abs() < 0.2 {
                if p[1] > max_y {
                    max_y = p[1];
                }
                count += 1;
            }
        }
        assert!(count > 0, "No vertices near the bond axis");
        // The toric patch reaches y=1.192. Contact patches go up to y=1.422.
        // The highest point near x=0 should be around y=1.19 (toric) or higher
        // if contact patches extend there.
        assert!(
            max_y > 1.0 && max_y < 1.6,
            "Bond-axis cross-section max y = {max_y:.3}, expected ~1.19–1.42. \
             The reentrant geometry is wrong."
        );
    }

    /// The SES area of two just-touching atoms = contact patches + toric saddle.
    /// Analytic: contacts = 2·2πR²(1+cos θ_max) = 56.2 Å², toric ≈ 11.9 Å².
    /// Total ≈ 68 Å² (LESS than VdW 72.6 because hidden caps are replaced by
    /// the smaller toric saddle).
    #[test]
    fn two_atom_area_analytic() {
        let mesh = two_atom_mesh();
        let ses_area = mesh.area();
        // Analytic total ≈ 68.1 Å². Allow 5% tolerance for discretization.
        assert!(
            (ses_area - 68.1).abs() / 68.1 < 0.05,
            "SES area {ses_area:.1}, expected 68.1 ± 5%. \
             Contact patches (56.2) + toric saddle (11.9)."
        );
    }

    /// The contact patch on each atom extends to theta_max = 56.74° from the
    /// bond axis. The contact patch endpoint (relative to atom 1 at (-1.7,0,0))
    /// is at (-0.768, 1.422, 0). If the contact/patch boundary is wrong, the
    /// probe geometry is wrong.
    #[test]
    fn two_atom_contact_patch_extent() {
        let mesh = two_atom_mesh();
        // Contact patch endpoint: atom1 center + R*(cos(tmax), sin(tmax), 0)
        // = (-1.7 + 1.7*cos(56.74°), 1.7*sin(56.74°), 0) = (-0.768, 1.422, 0)
        let target = [-0.768f32, 1.422, 0.0];
        let mut min_dist = f32::MAX;
        for p in &mesh.positions {
            let d = ((p[0] - target[0]).powi(2)
                + (p[1] - target[1]).powi(2)
                + (p[2] - target[2]).powi(2))
            .sqrt();
            min_dist = min_dist.min(d);
        }
        assert!(
            min_dist < 0.2,
            "No vertex near contact-patch endpoint (-0.768, 1.422, 0). \
             Closest: {min_dist:.3} Å. The contact/reentrant boundary is wrong."
        );
    }

    /// The mesh must be watertight (every edge shared by exactly 2 triangles).
    #[test]
    fn two_atom_watertight() {
        let mesh = two_atom_mesh();
        let report = mesh::topology::audit(&mesh);
        assert!(report.is_watertight(), "{report:?}");
    }
}
