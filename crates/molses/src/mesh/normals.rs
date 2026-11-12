//! Normal computation (M3).
//!
//! Two strategies are provided:
//! - **Field gradient** (`compute_gradient_normals`): the exact surface normal
//!   of a signed distance field is its gradient. This produces perfectly smooth
//!   shading even on a faceted mesh, because the gradient varies continuously
//!   across the surface. This is the default.
//! - **Angle-weighted face normals** (`compute_normals`): a mesh-only fallback
//!   that averages adjacent triangle normals. Visibly faceted on grid meshes.

use super::Mesh;
use crate::field::ScalarField;

/// Compute per-vertex normals from the gradient of `field` at each vertex
/// position (trilinear sample of the central-difference gradient), normalize,
/// and orient outward (toward increasing field values).
///
/// This is the preferred method: for a signed distance field, `∇f` *is* the
/// unit surface normal, so shading is smooth regardless of tessellation.
pub fn compute_gradient_normals(mesh: &mut Mesh, field: &ScalarField) {
    let mut normals = Vec::with_capacity(mesh.positions.len());
    for p in &mesh.positions {
        let g = field
            .gradient([p[0] as f64, p[1] as f64, p[2] as f64])
            .unwrap_or([0.0, 0.0, 1.0]);
        let len = ((g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) as f64).sqrt();
        if len > 1e-12 {
            let inv = (1.0 / len) as f32;
            normals.push([g[0] * inv, g[1] * inv, g[2] * inv]);
        } else {
            normals.push([0.0, 0.0, 1.0]);
        }
    }
    mesh.normals = normals;
}

/// Smooth per-vertex normals by averaging with one-ring neighbors (in place).
/// Purely cosmetic — does not move vertices.
pub fn smooth_normals(mesh: &mut Mesh, iterations: u32) {
    use std::collections::HashMap;
    if iterations == 0 || mesh.positions.is_empty() {
        return;
    }
    let n = mesh.positions.len();
    let mut neighbors: Vec<Vec<u32>> = vec![Vec::new(); n];
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
    let mut current = mesh.normals.clone();
    for _ in 0..iterations {
        let mut next = current.clone();
        for i in 0..n {
            let nb = &neighbors[i];
            if nb.is_empty() {
                continue;
            }
            let mut acc = current[i];
            for &j in nb {
                let nj = current[j as usize];
                acc[0] += nj[0];
                acc[1] += nj[1];
                acc[2] += nj[2];
            }
            let len = ((acc[0] * acc[0] + acc[1] * acc[1] + acc[2] * acc[2]) as f64).sqrt();
            if len > 1e-12 {
                let inv = (1.0 / len) as f32;
                next[i] = [acc[0] * inv, acc[1] * inv, acc[2] * inv];
            }
        }
        current = next;
    }
    mesh.normals = current;
    let _ = HashMap::<u32, u32>::new();
}

/// Compute per-vertex normals from triangle geometry (angle-weighted average),
/// normalize, and orient outward relative to `probe_point` (any point known to
/// be outside the surface, e.g. a far corner of the domain).
pub fn compute_normals(mesh: &mut Mesh, outside_point: [f32; 3]) {
    let n = mesh.positions.len();
    let mut normals = vec![[0.0f32; 3]; n];
    for t in &mesh.indices {
        let (i0, i1, i2) = (t[0] as usize, t[1] as usize, t[2] as usize);
        let a = mesh.positions[i0];
        let b = mesh.positions[i1];
        let c = mesh.positions[i2];
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let cross = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        for &i in &[i0, i1, i2] {
            normals[i][0] += cross[0];
            normals[i][1] += cross[1];
            normals[i][2] += cross[2];
        }
    }
    // Normalize and orient outward (away from outside_point? no: outward means
    // away from the interior; the normal should point toward outside_point).
    for (i, v) in normals.iter_mut().enumerate() {
        let mut n = *v;
        let len = ((n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) as f64).sqrt();
        if len > 1e-12 {
            let inv = (1.0 / len) as f32;
            n = [n[0] * inv, n[1] * inv, n[2] * inv];
        } else {
            n = [0.0, 0.0, 1.0];
        }
        // Orient toward outside_point.
        let p = mesh.positions[i];
        let to_out = [
            outside_point[0] - p[0],
            outside_point[1] - p[1],
            outside_point[2] - p[2],
        ];
        if n[0] * to_out[0] + n[1] * to_out[1] + n[2] * to_out[2] < 0.0 {
            n = [-n[0], -n[1], -n[2]];
        }
        *v = n;
    }
    mesh.normals = normals;
}
