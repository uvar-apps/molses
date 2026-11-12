//! MSMS `.vert` / `.face` compatible writer.
//!
//! This writes MSMS-*style* files so VMD/Chimera pipelines can consume our
//! meshes. It does **not** redistribute MSMS. Patch-type fields are best-effort
//! on grid meshes (classified as contact when unknown).

use std::io::Write;

use super::WriteMeta;
use crate::error::Result;
use crate::mesh::Mesh;

/// Write vertices then faces into a single stream (`.vert` layout followed by
/// `.face` layout, separated by comments). For separate files use
/// [`write_msms_split`].
///
/// `.vert` line: `x y z nx ny nz face_number closest_sphere_index patch_type`
/// (1-based indices; patch type 3 = contact on grid meshes by default).
pub fn write_msms(w: &mut dyn Write, mesh: &Mesh, meta: &WriteMeta) -> Result<()> {
    write_msms_vert(w, mesh, meta)?;
    writeln!(w)?;
    write_msms_face(w, mesh, meta)
}

/// Write only the `.vert` section.
pub fn write_msms_vert(w: &mut dyn Write, mesh: &Mesh, meta: &WriteMeta) -> Result<()> {
    let n_v = mesh.vertex_count();
    writeln!(w, "#molses {}", meta.tool)?;
    writeln!(w, "# {} {} {}", n_v, meta.atoms, meta.probe_radius)?;
    writeln!(w, "# x y z nx ny nz face closest_sphere patch_type")?;
    for (i, (p, n)) in mesh.positions.iter().zip(mesh.normals.iter()).enumerate() {
        let closest = mesh.nearest_atom.get(i).map(|a| a + 1).unwrap_or(0);
        writeln!(
            w,
            "{:9.3} {:9.3} {:9.3} {:9.3} {:9.3} {:9.3} {:7} {:7} {:2}",
            p[0], p[1], p[2], n[0], n[1], n[2], 0, closest, 3
        )?;
    }
    Ok(())
}

/// Write only the `.face` section.
pub fn write_msms_face(w: &mut dyn Write, mesh: &Mesh, meta: &WriteMeta) -> Result<()> {
    let n_f = mesh.triangle_count();
    writeln!(w, "#molses {}", meta.tool)?;
    writeln!(w, "# {} {} {}", n_f, meta.atoms, meta.probe_radius)?;
    writeln!(w, "# v1 v2 v3 patch_type face_number")?;
    for t in &mesh.indices {
        writeln!(
            w,
            "{:6} {:6} {:6} {:2} {:6}",
            t[0] + 1,
            t[1] + 1,
            t[2] + 1,
            3,
            0
        )?;
    }
    Ok(())
}

/// Write `.vert` content to `vert_w` and `.face` content to `face_w`.
pub fn write_msms_split(
    vert_w: &mut dyn Write,
    face_w: &mut dyn Write,
    mesh: &Mesh,
    meta: &WriteMeta,
) -> Result<()> {
    write_msms_vert(vert_w, mesh, meta)?;
    write_msms_face(face_w, mesh, meta)
}
