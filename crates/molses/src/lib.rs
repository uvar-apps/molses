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

    // --- Field construction (M2) ---
    let mut fld = field::ScalarField::new(dims, origin, h)?;
    let painted = field::paint::paint(structure.atoms(), probe, &fld);
    let edt_sq = field::edt::edt_squared(&painted.inside, dims);
    field::paint::assemble_signed_distance(&painted, &edt_sq, h, &mut fld);

    // For SAS the zero set is `f_sas` itself; for SES offset by probe.
    if matches!(opts.surface, SurfaceKind::Ses) {
        let data = fld.data_mut();
        for v in data.iter_mut() {
            *v += probe as f32;
        }
    }

    // Break exact zero ties: when `r_p / h` is an integer, whole EDT shells
    // carry f = 0.0 exactly, which collapses every edge crossing onto a grid
    // vertex and punches holes in the mesh. A smooth deterministic nudge far
    // below the grid spacing moves crossings off grid points without affecting
    // geometry.
    perturb_ties(&mut fld, h as f32 * 1e-4);

    // --- Mesh extraction (M3) ---
    let mut mesh = mesh::polygonize::polygonize(&fld);
    // Gradient normals are the exact surface normals of the SDF and shade
    // smoothly even on a faceted mesh (the default and preferred method).
    mesh::normals::compute_gradient_normals(&mut mesh, &fld);
    if opts.smoothing > 0 {
        // Taubin (λ|μ) smoothing: two Laplacian passes with opposite signs so
        // the mesh smooths without shrinking.
        let lambda = 0.5f32;
        let mu = -0.53f32;
        for _ in 0..opts.smoothing {
            mesh::smooth::laplacian(&mut mesh, 1, lambda);
            mesh::smooth::laplacian(&mut mesh, 1, mu);
        }
        // Recompute normals from the geometry after moving vertices; then
        // smooth them so shading stays clean.
        mesh::normals::compute_gradient_normals(&mut mesh, &fld);
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

/// Add a smooth deterministic perturbation of amplitude `eps` to every sample
/// of `fld`, breaking exact zero ties (see `ses`).
fn perturb_ties(fld: &mut field::ScalarField, eps: f32) {
    let [nx, ny, nz] = fld.dims();
    let data = fld.data_mut();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let i = x + y * nx + z * nx * ny;
                // Cheap smooth hash in [-1, 1].
                let s = (x as f32 * 12.9898 + y as f32 * 78.233 + z as f32 * 37.719).sin();
                data[i] += eps * s;
            }
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
