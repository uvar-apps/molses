//! Normal computation.
//!
//! Shading follows the geometry: vertices projected onto an exact analytic
//! patch take that patch's own normal, so shading carries no tessellation
//! noise at all.
//!
//! - [`compute_normals`]: angle-weighted average of adjacent face normals,
//!   then a per-component orientation vote. This is the base field and the
//!   fallback for vertices on no recognised patch.
//! - [`apply_patch_normals`]: overwrite with the exact patch normal —
//!   `normalize(p − c)` on a contact sphere, or the direction towards the
//!   probe-centre circle on a canal torus — then diffuse those values into
//!   unpatched vertices so patch boundaries stay continuous.
//! - [`repair_flipped_normals`]: majority vote over the 1-ring to remove
//!   isolated inverted normals (which render as black specks).
use super::Mesh;
use super::polygonize::nearest_circle_point;
use crate::field::ScalarField;
use crate::structure::Atom;

/// The analytic surface patch a vertex was projected onto.
///
/// Produced by the exact-patch projection in [`crate::ses`] and consumed by
/// [`apply_patch_normals`].
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Patch {
    /// Not on a recognised patch — keep the geometric normal.
    None,
    /// Convex contact patch of atom `i` (sphere of radius `r_i`).
    Contact(u32),
    /// Concave canal patch swept by a probe touching atoms `i` and `j`.
    Canal(u32, u32),
}

/// Outward unit normal of the VdW sphere of atom `k` at `p`.
fn sphere_normal(atoms: &[Atom], k: usize, p: [f32; 3]) -> Option<[f32; 3]> {
    let a = atoms.get(k)?;
    let d = [
        p[0] as f64 - a.position[0],
        p[1] as f64 - a.position[1],
        p[2] as f64 - a.position[2],
    ];
    let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    if len < 1e-12 {
        return None;
    }
    Some([
        (d[0] / len) as f32,
        (d[1] / len) as f32,
        (d[2] / len) as f32,
    ])
}

/// Outward unit normal of the canal torus of pair (`i`, `j`) at `p`.
///
/// The torus is the set of points at distance `probe` from the circle traced
/// by the probe centre, and it faces the solvent, so the outward normal runs
/// from `p` towards the nearest point of that circle.
fn canal_normal(atoms: &[Atom], probe: f64, p: [f32; 3], i: usize, j: usize) -> Option<[f32; 3]> {
    let (ai, aj) = (atoms.get(i)?, atoms.get(j)?);
    let ri = ai.radius + probe;
    let rj = aj.radius + probe;
    let e = [
        aj.position[0] - ai.position[0],
        aj.position[1] - ai.position[1],
        aj.position[2] - ai.position[2],
    ];
    let d = (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt();
    if d < 1e-9 || d >= ri + rj || d <= (ri - rj).abs() {
        return None;
    }
    let u = [e[0] / d, e[1] / d, e[2] / d];
    let t = (d * d + ri * ri - rj * rj) / (2.0 * d);
    let rc = (ri * ri - t * t).max(0.0).sqrt();
    if rc < 1e-9 {
        return None;
    }
    let o = [
        ai.position[0] + t * u[0],
        ai.position[1] + t * u[1],
        ai.position[2] + t * u[2],
    ];
    let (npt, _) = nearest_circle_point([p[0] as f64, p[1] as f64, p[2] as f64], o, rc, u);
    let v = [
        npt[0] - p[0] as f64,
        npt[1] - p[1] as f64,
        npt[2] - p[2] as f64,
    ];
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len < 1e-12 {
        return None;
    }
    Some([
        (v[0] / len) as f32,
        (v[1] / len) as f32,
        (v[2] / len) as f32,
    ])
}

/// Vertex-adjacency lists (unordered, deduplicated).
pub(crate) fn vertex_neighbors(mesh: &Mesh) -> Vec<Vec<u32>> {
    let mut nb: Vec<Vec<u32>> = vec![Vec::new(); mesh.positions.len()];
    for t in &mesh.indices {
        for k in 0..3 {
            let a = t[k] as usize;
            for &j in &[t[(k + 1) % 3], t[(k + 2) % 3]] {
                if !nb[a].contains(&j) {
                    nb[a].push(j);
                }
            }
        }
    }
    nb
}

/// Repair isolated inverted normals by majority vote of the 1-ring.
///
/// A normal that disagrees with most of its neighbours cannot survive on a
/// consistently wound surface — it is a leftover of a degenerate triangle or
/// a wrong patch assignment, and it renders as a black speck. Two passes let
/// a repair cascade through a short run of consecutive bad vertices.
pub fn repair_flipped_normals(mesh: &mut Mesh) {
    let n = mesh.normals.len();
    if n == 0 || mesh.normals.len() != mesh.positions.len() {
        return;
    }
    let nb = vertex_neighbors(mesh);
    for _ in 0..2 {
        let cur = mesh.normals.clone();
        let mut fixed = 0usize;
        for (i, m) in mesh.normals.iter_mut().enumerate() {
            if nb[i].is_empty() {
                continue;
            }
            let mut acc = [0.0f64; 3];
            for &j in &nb[i] {
                let o = cur[j as usize];
                acc[0] += o[0] as f64;
                acc[1] += o[1] as f64;
                acc[2] += o[2] as f64;
            }
            let me = cur[i];
            let dot = acc[0] * me[0] as f64 + acc[1] * me[1] as f64 + acc[2] * me[2] as f64;
            if dot >= 0.0 {
                continue;
            }
            let l = (acc[0] * acc[0] + acc[1] * acc[1] + acc[2] * acc[2]).sqrt();
            if l < 1e-12 {
                continue;
            }
            *m = [
                (acc[0] / l) as f32,
                (acc[1] / l) as f32,
                (acc[2] / l) as f32,
            ];
            fixed += 1;
        }
        if fixed == 0 {
            break;
        }
    }
}

/// Shade the mesh with the *analytic* normal of each exact patch.
///
/// Both patch normals are smooth functions of position, so triangulation
/// noise disappears from the shading — which is what an orange-peel artifact
/// actually is. Vertices whose analytic normal disagrees with the geometric
/// one by more than 60° keep the geometric normal (a wrong patch must never
/// flip a face), and the remaining unpatched vertices (3-body pockets and the
/// transition band) get a harmonic fill from their neighbours so the field
/// stays continuous across patch boundaries.
pub fn apply_patch_normals(mesh: &mut Mesh, atoms: &[Atom], probe: f64, patch: &[Patch]) {
    let n = mesh.positions.len();
    if mesh.normals.len() != n || patch.len() != n {
        return;
    }
    let mut exact = vec![false; n];
    for i in 0..n {
        let cand = match patch[i] {
            Patch::None => continue,
            Patch::Contact(k) => sphere_normal(atoms, k as usize, mesh.positions[i]),
            Patch::Canal(a, b) => {
                canal_normal(atoms, probe, mesh.positions[i], a as usize, b as usize)
            }
        };
        let Some(c) = cand else { continue };
        let g = mesh.normals[i];
        if c[0] * g[0] + c[1] * g[1] + c[2] * g[2] > 0.5 {
            mesh.normals[i] = c;
            exact[i] = true;
        }
    }
    // Harmonic fill: diffuse the (smooth) exact normals into unpatched verts.
    let nb = vertex_neighbors(mesh);
    for _ in 0..16 {
        let cur = mesh.normals.clone();
        for i in 0..n {
            if exact[i] || nb[i].is_empty() {
                continue;
            }
            let mut acc = [0.0f64; 3];
            for &j in &nb[i] {
                let m = cur[j as usize];
                acc[0] += m[0] as f64;
                acc[1] += m[1] as f64;
                acc[2] += m[2] as f64;
            }
            let l = (acc[0] * acc[0] + acc[1] * acc[1] + acc[2] * acc[2]).sqrt();
            if l < 1e-12 {
                continue;
            }
            let mut v = [
                (acc[0] / l) as f32,
                (acc[1] / l) as f32,
                (acc[2] / l) as f32,
            ];
            if v[0] * cur[i][0] + v[1] * cur[i][1] + v[2] * cur[i][2] < 0.0 {
                v = [-v[0], -v[1], -v[2]];
            }
            mesh.normals[i] = v;
        }
    }
}

/// Angle-weighted face normals with optional field-sign orientation.
pub fn compute_normals(mesh: &mut Mesh, field: Option<&ScalarField>) {
    let n = mesh.positions.len();
    let mut normals = vec![[0.0f32; 3]; n];
    for t in &mesh.indices {
        let (i0, i1, i2) = (t[0] as usize, t[1] as usize, t[2] as usize);
        let a = mesh.positions[i0];
        let b = mesh.positions[i1];
        let c = mesh.positions[i2];
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let ba = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let bc = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
        let ca = [a[0] - c[0], a[1] - c[1], a[2] - c[2]];
        let cb = [b[0] - c[0], b[1] - c[1], b[2] - c[2]];
        let cross = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        let angle_at = |u: [f32; 3], v: [f32; 3]| -> f32 {
            let lu = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
            let lv = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
            if lu < 1e-20 || lv < 1e-20 {
                return 0.0;
            }
            let dot = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (lu * lv);
            dot.clamp(-1.0, 1.0).acos()
        };
        let w0 = angle_at(ab, ac);
        let w1 = angle_at(ba, bc);
        let w2 = angle_at(ca, cb);
        for (vi, w) in [(i0, w0), (i1, w1), (i2, w2)] {
            normals[vi][0] += cross[0] * w;
            normals[vi][1] += cross[1] * w;
            normals[vi][2] += cross[2] * w;
        }
    }
    for v in normals.iter_mut() {
        let mut nv = *v;
        let len = ((nv[0] * nv[0] + nv[1] * nv[1] + nv[2] * nv[2]) as f64).sqrt();
        if len > 1e-12 {
            let inv = (1.0 / len) as f32;
            nv = [nv[0] * inv, nv[1] * inv, nv[2] * inv];
        } else {
            nv = [0.0, 0.0, 1.0];
        }
        *v = nv;
    }
    mesh.normals = normals;
    if let Some(fld) = field {
        orient_outward(mesh, fld);
    }
}

/// Make every vertex normal point to the low-field (solvent) side.
///
/// Flipping is decided **per connected component** from the aggregate
/// agreement of its vertices, never per vertex: an isolated flip is what
/// shows up as a black speck in the shading, and probing the field at a
/// single point is unreliable inside crevices where the voxel is ambiguous.
fn orient_outward(mesh: &mut Mesh, field: &ScalarField) {
    let n = mesh.normals.len();
    let nb = vertex_neighbors(mesh);
    let mut seen = vec![false; n];
    let eps = field.spacing();
    let mut stack = Vec::with_capacity(64);
    for seed in 0..n {
        if seen[seed] {
            continue;
        }
        let mut comp: Vec<usize> = Vec::new();
        stack.clear();
        stack.push(seed);
        seen[seed] = true;
        while let Some(i) = stack.pop() {
            comp.push(i);
            for &j in &nb[i] {
                let j = j as usize;
                if !seen[j] {
                    seen[j] = true;
                    stack.push(j);
                }
            }
        }
        // Vote: +1 when the normal points at increasing field value.
        let mut vote = 0i64;
        for &i in &comp {
            let p = mesh.positions[i];
            let nv = mesh.normals[i];
            let outside = field.sample([
                p[0] as f64 + nv[0] as f64 * eps,
                p[1] as f64 + nv[1] as f64 * eps,
                p[2] as f64 + nv[2] as f64 * eps,
            ]);
            let inside = field.sample([
                p[0] as f64 - nv[0] as f64 * eps,
                p[1] as f64 - nv[1] as f64 * eps,
                p[2] as f64 - nv[2] as f64 * eps,
            ]);
            if let (Some(a), Some(b)) = (outside, inside) {
                if a > b {
                    vote += 1;
                } else {
                    vote -= 1;
                }
            }
        }
        if vote < 0 {
            for &i in &comp {
                let v = &mut mesh.normals[i];
                *v = [-v[0], -v[1], -v[2]];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::polygonize::polygonize;

    fn sphere_field(r: f64, h: f64) -> ScalarField {
        let extent = r + 0.6;
        let n = (2.0 * extent / h).ceil() as usize + 1;
        let mut f = ScalarField::new([n, n, n], [-extent; 3], h).unwrap();
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let p = f.position(x, y, z);
                    let d = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                    *f.get_mut(x, y, z).unwrap() = (d - r) as f32;
                }
            }
        }
        f
    }

    /// Regression: every normal on a sphere must point away from the center.
    /// The old code flipped normals toward a single bounding-box corner, which
    /// inverted roughly half of them on any sizable molecule.
    #[test]
    fn sphere_normals_point_outward() {
        let r = 1.0f64;
        let field = sphere_field(r, 0.3);
        let mut mesh = polygonize(&field, None, 1.4);
        compute_normals(&mut mesh, Some(&field));
        let mut flipped = 0;
        for (i, p) in mesh.positions.iter().enumerate() {
            let n = mesh.normals[i];
            let dot = p[0] * n[0] + p[1] * n[1] + p[2] * n[2];
            if dot < 0.0 {
                flipped += 1;
            }
        }
        assert_eq!(flipped, 0, "{flipped} normals point inward on a sphere");
    }

    /// Angle-weighted normals should be unit length.
    #[test]
    fn normals_are_unit_length() {
        let r = 1.0f64;
        let field = sphere_field(r, 0.3);
        let mut mesh = polygonize(&field, None, 1.4);
        compute_normals(&mut mesh, Some(&field));
        for n in &mesh.normals {
            let len = ((n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) as f64).sqrt();
            assert!((len - 1.0).abs() < 1e-4, "normal length {len}");
        }
    }
}
