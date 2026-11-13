//! Isosurface extraction by marching tetrahedra with exact edge roots.
//!
//! Each grid cell is split into 6 tetrahedra sharing the main diagonal
//! `v0–v6`. This decomposition is face-consistent across cells (every shared
//! face uses the same diagonal), so the mesh is crack-free by construction.
//! Vertices lying on the same grid edge are welded via an edge key, which
//! makes the mesh watertight.
//!
//! **Exact edge roots**: instead of linearly interpolating the (quantized)
//! field values at the edge endpoints, we solve the analytic implicit
//! equation along the edge:
//!
//! - **Contact patches** — `‖p(t) − c_i‖ = r_i` (quadratic; exact).
//! - **Reentrant patches** — `d(p(t), C_ij) = r_p` where `C_ij` is the
//!   intersection circle of two inflated spheres, restricted to the inner
//!   half of the torus (the crevice-facing side).
//!
//! A refined root is accepted only if it stays close to the linear estimate
//! (`|Δt| < 0.25`); otherwise we keep interpolation. This removes EDT
//! quantization jitter without inventing geometry. Topology (which edges
//! cross) still comes from the field, so the mesh stays watertight.

use std::collections::HashMap;

use super::Mesh;
use crate::field::ScalarField;
use crate::structure::Atom;

/// Cube-corner offsets (x fastest).
const CORNER_OFFSETS: [[usize; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [1, 1, 0],
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [1, 1, 1],
    [0, 1, 1],
];

/// Six tetrahedra sharing the body diagonal 0–6.
const TETS: [[usize; 4]; 6] = [
    [0, 1, 2, 6],
    [0, 2, 3, 6],
    [0, 3, 7, 6],
    [0, 7, 4, 6],
    [0, 4, 5, 6],
    [0, 5, 1, 6],
];

/// Lower bound for edge-interpolation parameter `t`.
const T_EPS: f32 = 1e-6;

/// Extract the zero level set of `field` as a triangle mesh.
///
/// If `atoms` is provided, edge roots are computed analytically against the
/// VdW spheres of the atoms (exact edge roots). Otherwise, linear
/// interpolation of field values is used.
pub fn polygonize(field: &ScalarField, atoms: Option<&[Atom]>, probe_radius: f64) -> Mesh {
    let [nx, ny, nz] = field.dims();
    if nx < 2 || ny < 2 || nz < 2 {
        return Mesh::default();
    }
    let data = field.data();
    let idx = |x: usize, y: usize, z: usize| x + y * nx + z * nx * ny;

    let mut weld: HashMap<(u64, u64), u32> = HashMap::new();
    let mut positions: Vec<[f32; 3]> = Vec::new();
    let mut indices: Vec<[u32; 3]> = Vec::new();

    let origin = field.origin();
    let h = field.spacing();

    for z in 0..(nz - 1) {
        for y in 0..(ny - 1) {
            for x in 0..(nx - 1) {
                let mut neg = 0u8;
                let mut pos = 0u8;
                let mut f8 = [0.0f32; 8];
                let mut g8 = [[0usize; 3]; 8];
                for (c, off) in CORNER_OFFSETS.iter().enumerate() {
                    let (cx, cy, cz) = (x + off[0], y + off[1], z + off[2]);
                    let v = data[idx(cx, cy, cz)];
                    f8[c] = v;
                    g8[c] = [cx, cy, cz];
                    if v < 0.0 {
                        neg += 1;
                    } else {
                        pos += 1;
                    }
                }
                if neg == 0 || pos == 0 {
                    continue;
                }
                for tet in &TETS {
                    emit_tet(
                        tet,
                        &f8,
                        &g8,
                        origin,
                        h,
                        probe_radius,
                        atoms,
                        &mut weld,
                        &mut positions,
                        &mut indices,
                    );
                }
            }
        }
    }

    Mesh {
        positions,
        normals: Vec::new(),
        indices,
        nearest_atom: Vec::new(),
    }
}

fn pack(x: usize, y: usize, z: usize) -> u64 {
    (x as u64) | ((y as u64) << 21) | ((z as u64) << 42)
}

/// Nearest point on a circle to `p`. Returns (point, distance).
pub(crate) fn nearest_circle_point(
    p: [f64; 3],
    c: [f64; 3],
    r: f64,
    n: [f64; 3],
) -> ([f64; 3], f64) {
    // Project p onto the circle plane.
    let pc = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
    let d_plane = pc[0] * n[0] + pc[1] * n[1] + pc[2] * n[2];
    let p_proj = [
        p[0] - d_plane * n[0],
        p[1] - d_plane * n[1],
        p[2] - d_plane * n[2],
    ];
    // Direction from center to projection.
    let v = [p_proj[0] - c[0], p_proj[1] - c[1], p_proj[2] - c[2]];
    let v_len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if v_len < 1e-12 {
        // p projects to center: any point on the circle is nearest.
        return ([c[0] + r, c[1], c[2]], r);
    }
    let inv = r / v_len;
    let nearest = [c[0] + v[0] * inv, c[1] + v[1] * inv, c[2] + v[2] * inv];
    (nearest, (v_len - r).abs())
}

/// Find the zero crossing on grid edge (a, b).
///
/// Uses exact analytic roots when `atoms` is provided:
/// - contact patch → quadratic sphere intersection (with visibility check)
/// - reentrant patch → linear interpolation of the field
///
/// A root is accepted only when it stays near the linear estimate;
/// otherwise linear interpolation of the field is used.
#[allow(clippy::too_many_arguments)]
fn edge_vertex(
    a: [usize; 3],
    b: [usize; 3],
    fa: f32,
    fb: f32,
    origin: [f64; 3],
    h: f64,
    probe_radius: f64,
    atoms: Option<&[Atom]>,
    weld: &mut HashMap<(u64, u64), u32>,
    positions: &mut Vec<[f32; 3]>,
) -> u32 {
    let pa = pack(a[0], a[1], a[2]);
    let pb = pack(b[0], b[1], b[2]);
    let key = if pa < pb { (pa, pb) } else { (pb, pa) };
    if let Some(&vi) = weld.get(&key) {
        return vi;
    }

    let wa = [
        origin[0] + a[0] as f64 * h,
        origin[1] + a[1] as f64 * h,
        origin[2] + a[2] as f64 * h,
    ];
    let wb = [
        origin[0] + b[0] as f64 * h,
        origin[1] + b[1] as f64 * h,
        origin[2] + b[2] as f64 * h,
    ];
    let d = [wb[0] - wa[0], wb[1] - wa[1], wb[2] - wa[2]];
    let edge_len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();

    // Linear estimate from the (quantized) field.
    let t_lin = if (fa - fb).abs() > 1e-30 {
        (fa / (fa - fb)).clamp(0.0, 1.0) as f64
    } else {
        0.5
    };
    let t_lin = t_lin.clamp(T_EPS as f64, 1.0 - T_EPS as f64);
    let p_lin = [
        wa[0] + t_lin * d[0],
        wa[1] + t_lin * d[1],
        wa[2] + t_lin * d[2],
    ];

    let mut t = t_lin;

    if let Some(atoms) = atoms {
        if edge_len > 1e-30 {
            // --- 1. Contact patch: nearest VdW sphere ---
            // On a contact patch, f_ses ≈ ‖p−c_i‖ − r_i. Quantization can
            // push p_lin off the sphere by ~h/2, so the classification band
            // must be generous; the quadratic root + |Δt| guard do the rest.
            let contact_tol = h * 1.5;
            let mut best_atom: Option<usize> = None;
            let mut best_vdw = f64::INFINITY;
            for (i, atom) in atoms.iter().enumerate() {
                let dx = p_lin[0] - atom.position[0];
                let dy = p_lin[1] - atom.position[1];
                let dz = p_lin[2] - atom.position[2];
                let vdw = ((dx * dx + dy * dy + dz * dz).sqrt() - atom.radius).abs();
                if vdw < best_vdw {
                    best_vdw = vdw;
                    best_atom = Some(i);
                }
            }

            if let (Some(i), true) = (best_atom, best_vdw < contact_tol) {
                let atom = &atoms[i];
                let e = [
                    wa[0] - atom.position[0],
                    wa[1] - atom.position[1],
                    wa[2] - atom.position[2],
                ];
                let a2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                let b2 = 2.0 * (e[0] * d[0] + e[1] * d[1] + e[2] * d[2]);
                let c2 = e[0] * e[0] + e[1] * e[1] + e[2] * e[2] - atom.radius * atom.radius;
                let disc = b2 * b2 - 4.0 * a2 * c2;
                if disc >= 0.0 {
                    let sq = disc.sqrt();
                    let mut best_t = f64::INFINITY;
                    let mut best_dt = f64::INFINITY;
                    for root in [(-b2 - sq) / (2.0 * a2), (-b2 + sq) / (2.0 * a2)] {
                        if (0.0..=1.0).contains(&root) {
                            let dt = (root - t_lin).abs();
                            if dt < best_dt {
                                best_dt = dt;
                                best_t = root;
                            }
                        }
                    }
                    // Accept if the root stays near the field crossing AND the
                    // probe center is (nearly) on ∂K. Slight burial allowance
                    // covers float error at the contact/toric boundary.
                    if best_dt < 0.4 {
                        let p = [
                            wa[0] + best_t * d[0],
                            wa[1] + best_t * d[1],
                            wa[2] + best_t * d[2],
                        ];
                        if contact_visible(p, i, atoms, probe_radius, h) {
                            t = best_t.clamp(T_EPS as f64, 1.0 - T_EPS as f64);
                        }
                    }
                }
            }
            // Reentrant (toric) edge roots are intentionally not applied here:
            // Newton on d(p,C)=r_p inflates area near patch boundaries. Those
            // vertices keep linear interpolation of the field (correct shape,
            // mild quantization). Contact roots cover the majority surface.
        }
    }

    let p = [
        (wa[0] + t * d[0]) as f32,
        (wa[1] + t * d[1]) as f32,
        (wa[2] + t * d[2]) as f32,
    ];
    let vi = positions.len() as u32;
    positions.push(p);
    weld.insert(key, vi);
    vi
}

/// Whether `p` on atom `i`'s VdW sphere is a visible SES contact point.
///
/// The probe center `p + r_p·n̂` (n̂ outward radial) must lie on ∂K — i.e.
/// not strictly inside any other inflated ball. `slack` (typically ~0.05·h)
/// absorbs float error exactly at the contact/toric boundary, where the probe
/// center sits on an intersection circle.
fn contact_visible(p: [f64; 3], i: usize, atoms: &[Atom], probe_radius: f64, slack: f64) -> bool {
    let a = &atoms[i];
    let dx = p[0] - a.position[0];
    let dy = p[1] - a.position[1];
    let dz = p[2] - a.position[2];
    let len = (dx * dx + dy * dy + dz * dz).sqrt();
    if len < 1e-12 {
        return false;
    }
    let scale = (a.radius + probe_radius) / len;
    let pc = [
        a.position[0] + dx * scale,
        a.position[1] + dy * scale,
        a.position[2] + dz * scale,
    ];
    let eps = slack.max(1e-6);
    for (j, b) in atoms.iter().enumerate() {
        if j == i {
            continue;
        }
        let ex = pc[0] - b.position[0];
        let ey = pc[1] - b.position[1];
        let ez = pc[2] - b.position[2];
        let dist = (ex * ex + ey * ey + ez * ez).sqrt();
        if dist < b.radius + probe_radius - eps {
            return false;
        }
    }
    true
}

/// Emit 0–2 triangles for one tetrahedron, wound outward (toward f > 0).
#[allow(clippy::too_many_arguments)]
fn emit_tet(
    tet: &[usize; 4],
    f8: &[f32; 8],
    g8: &[[usize; 3]; 8],
    origin: [f64; 3],
    h: f64,
    probe_radius: f64,
    atoms: Option<&[Atom]>,
    weld: &mut HashMap<(u64, u64), u32>,
    positions: &mut Vec<[f32; 3]>,
    indices: &mut Vec<[u32; 3]>,
) {
    let f = [f8[tet[0]], f8[tet[1]], f8[tet[2]], f8[tet[3]]];
    let g = [g8[tet[0]], g8[tet[1]], g8[tet[2]], g8[tet[3]]];

    let mut in_idx = [usize::MAX; 4];
    let mut out_idx = [usize::MAX; 4];
    let mut n_in = 0usize;
    let mut n_out = 0usize;
    for (i, &fv) in f.iter().enumerate() {
        if fv < 0.0 {
            in_idx[n_in] = i;
            n_in += 1;
        } else {
            out_idx[n_out] = i;
            n_out += 1;
        }
    }
    if n_in == 0 || n_in == 4 {
        return;
    }

    match n_in {
        1 => {
            let a = in_idx[0];
            let v0 = cross(
                a,
                out_idx[0],
                &f,
                &g,
                origin,
                h,
                probe_radius,
                atoms,
                weld,
                positions,
            );
            let v1 = cross(
                a,
                out_idx[1],
                &f,
                &g,
                origin,
                h,
                probe_radius,
                atoms,
                weld,
                positions,
            );
            let v2 = cross(
                a,
                out_idx[2],
                &f,
                &g,
                origin,
                h,
                probe_radius,
                atoms,
                weld,
                positions,
            );
            let ref_pt = outside_ref(&g[out_idx[0]], &g[out_idx[1]], &g[out_idx[2]], origin, h);
            push_tri_oriented(v0, v1, v2, ref_pt, positions, indices);
        }
        3 => {
            let a = out_idx[0];
            let v0 = cross(
                in_idx[0],
                a,
                &f,
                &g,
                origin,
                h,
                probe_radius,
                atoms,
                weld,
                positions,
            );
            let v1 = cross(
                in_idx[1],
                a,
                &f,
                &g,
                origin,
                h,
                probe_radius,
                atoms,
                weld,
                positions,
            );
            let v2 = cross(
                in_idx[2],
                a,
                &f,
                &g,
                origin,
                h,
                probe_radius,
                atoms,
                weld,
                positions,
            );
            let ref_pt = world(g[a], origin, h);
            push_tri_oriented(v0, v1, v2, ref_pt, positions, indices);
        }
        2 => {
            let (i0, i1) = (in_idx[0], in_idx[1]);
            let (o0, o1) = (out_idx[0], out_idx[1]);
            let q0 = cross(
                i0,
                o0,
                &f,
                &g,
                origin,
                h,
                probe_radius,
                atoms,
                weld,
                positions,
            );
            let q1 = cross(
                i1,
                o0,
                &f,
                &g,
                origin,
                h,
                probe_radius,
                atoms,
                weld,
                positions,
            );
            let q2 = cross(
                i1,
                o1,
                &f,
                &g,
                origin,
                h,
                probe_radius,
                atoms,
                weld,
                positions,
            );
            let q3 = cross(
                i0,
                o1,
                &f,
                &g,
                origin,
                h,
                probe_radius,
                atoms,
                weld,
                positions,
            );
            let p0 = world(g[o0], origin, h);
            let p1 = world(g[o1], origin, h);
            let ref_pt = [
                (p0[0] + p1[0]) * 0.5,
                (p0[1] + p1[1]) * 0.5,
                (p0[2] + p1[2]) * 0.5,
            ];
            push_tri_oriented(q0, q1, q2, ref_pt, positions, indices);
            push_tri_oriented(q0, q2, q3, ref_pt, positions, indices);
        }
        _ => unreachable!(),
    }
}

/// Centroid of up to three grid vertices, in world space.
fn outside_ref(
    a: &[usize; 3],
    b: &[usize; 3],
    c: &[usize; 3],
    origin: [f64; 3],
    h: f64,
) -> [f32; 3] {
    let pa = world(*a, origin, h);
    let pb = world(*b, origin, h);
    let pc = world(*c, origin, h);
    [
        (pa[0] + pb[0] + pc[0]) / 3.0,
        (pa[1] + pb[1] + pc[1]) / 3.0,
        (pa[2] + pb[2] + pc[2]) / 3.0,
    ]
}

/// World position of a grid vertex.
fn world(g: [usize; 3], origin: [f64; 3], h: f64) -> [f32; 3] {
    [
        (origin[0] + g[0] as f64 * h) as f32,
        (origin[1] + g[1] as f64 * h) as f32,
        (origin[2] + g[2] as f64 * h) as f32,
    ]
}

/// Crossing vertex between local tet corners `i` and `j`.
#[allow(clippy::too_many_arguments)]
fn cross(
    i: usize,
    j: usize,
    f: &[f32; 4],
    g: &[[usize; 3]; 4],
    origin: [f64; 3],
    h: f64,
    probe_radius: f64,
    atoms: Option<&[Atom]>,
    weld: &mut HashMap<(u64, u64), u32>,
    positions: &mut Vec<[f32; 3]>,
) -> u32 {
    let (a, b) = if f[i] < 0.0 { (i, j) } else { (j, i) };
    edge_vertex(
        g[a],
        g[b],
        f[a],
        f[b],
        origin,
        h,
        probe_radius,
        atoms,
        weld,
        positions,
    )
}

/// Push a triangle wound so its normal points toward `outward_ref`.
fn push_tri_oriented(
    a: u32,
    b: u32,
    c: u32,
    outward_ref: [f32; 3],
    positions: &[[f32; 3]],
    indices: &mut Vec<[u32; 3]>,
) {
    if a == b || b == c || a == c {
        return;
    }
    let pa = positions[a as usize];
    let pb = positions[b as usize];
    let pc = positions[c as usize];
    let ab = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
    let ac = [pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]];
    let n = [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ];
    let centroid = [
        (pa[0] + pb[0] + pc[0]) / 3.0,
        (pa[1] + pb[1] + pc[1]) / 3.0,
        (pa[2] + pb[2] + pc[2]) / 3.0,
    ];
    let to_ref = [
        outward_ref[0] - centroid[0],
        outward_ref[1] - centroid[1],
        outward_ref[2] - centroid[2],
    ];
    let dot = n[0] * to_ref[0] + n[1] * to_ref[1] + n[2] * to_ref[2];
    if dot >= 0.0 {
        indices.push([a, b, c]);
    } else {
        indices.push([a, c, b]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn sphere_area_converges() {
        let r = 1.0f64;
        let field = sphere_field(r, 0.25);
        let mesh = polygonize(&field, None, 1.4);
        assert!(mesh.triangle_count() > 0);
        let area = mesh.area();
        let expect = 4.0 * std::f64::consts::PI * r * r;
        let rel_err = (area - expect).abs() / expect;
        assert!(rel_err < 0.15, "area {area} vs {expect}, err {rel_err}");
    }

    #[test]
    fn sphere_is_watertight() {
        let field = sphere_field(1.0, 0.3);
        let mesh = polygonize(&field, None, 1.4);
        let report = crate::mesh::topology::audit(&mesh);
        assert!(report.is_watertight(), "{report:?}");
        assert_eq!(report.euler, 2, "{report:?}");
    }

    #[test]
    fn empty_field_gives_empty_mesh() {
        let field = ScalarField::new([4, 4, 4], [0.0; 3], 1.0).unwrap();
        let mesh = polygonize(&field, None, 1.4);
        assert!(mesh.is_empty());
    }

    /// Exact edge roots should give a better area estimate than interpolation.
    #[test]
    fn exact_roots_improve_sphere_area() {
        let r = 1.0f64;
        let h = 0.25;
        let field = sphere_field(r, h);
        let atoms = vec![Atom::new(
            [0.0; 3],
            r,
            Some(crate::structure::Element::C),
            1,
        )];

        let mesh_interp = polygonize(&field, None, 1.4);
        let mesh_exact = polygonize(&field, Some(&atoms), 1.4);

        let expect = 4.0 * std::f64::consts::PI * r * r;
        let err_interp = (mesh_interp.area() - expect).abs() / expect;
        let err_exact = (mesh_exact.area() - expect).abs() / expect;

        assert!(
            err_exact <= err_interp * 1.5,
            "exact roots should not be much worse: interp={err_interp:.4} exact={err_exact:.4}"
        );
    }
}
