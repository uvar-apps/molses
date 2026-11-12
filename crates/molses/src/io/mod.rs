//! Surface file writers (M4).

pub mod msms;
pub mod obj;

use std::io::Write;

use crate::error::Result;
use crate::mesh::Mesh;

/// Metadata embedded in output headers.
#[derive(Debug, Clone, Default)]
pub struct WriteMeta {
    /// Tool name and version.
    pub tool: String,
    /// Source path.
    pub source: String,
    /// Probe radius in Å.
    pub probe_radius: f64,
    /// Grid spacing in Å.
    pub grid_spacing: f64,
    /// Number of atoms used.
    pub atoms: usize,
    /// Suppress timestamp lines (for deterministic golden tests).
    pub deterministic: bool,
}

/// Write a mesh in the requested format.
pub fn write_mesh(w: &mut dyn Write, mesh: &Mesh, format: &str, meta: &WriteMeta) -> Result<()> {
    match format {
        "obj" => obj::write_obj(w, mesh, meta),
        "msms" => msms::write_msms(w, mesh, meta),
        other => Err(crate::error::SesError::UnsupportedFormat(other.into())),
    }
}
