//! Isosurface extraction by marching tetrahedra.
//!
//! Each grid cell is split into 6 tetrahedra sharing the main diagonal
//! `v0–v6`. This decomposition is face-consistent across cells (every shared
//! face uses the same diagonal), so the mesh is crack-free by construction.
//! Vertices lying on the same grid edge are welded via an edge key, which
//! makes the mesh watertight.
//!
//! `f > 0` is the exterior side, `f < 0` is the solid. Triangles are wound so
//! the geometric normal points toward increasing `f` (outward).

use std::collections::HashMap;

use super::Mesh;
use crate::field::ScalarField;

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

/// Lower bound for edge-interpolation parameter `t`, keeping crossings strictly
/// inside their edge so vertices never coincide with grid points.
const T_EPS: f32 = 1e-6;

/// Extract the zero level set of `field` as a triangle mesh.
pub fn polygonize(field: &ScalarField) -> Mesh {
    let [nx, ny, nz] = field.dims();
    if nx < 2 || ny < 2 || nz < 2 {
        return Mesh::default();
    }
    let data = field.data();
    let idx = |x: usize, y: usize, z: usize| x + y * nx + z * nx * ny;

    // Collision-free weld map: (min packed, max packed) -> vertex index.
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

/// Return the welded vertex index for the zero crossing on grid edge (a, b).
#[allow(clippy::too_many_arguments)]
fn edge_vertex(
    a: [usize; 3],
    b: [usize; 3],
    fa: f32,
    fb: f32,
    origin: [f64; 3],
    h: f64,
    weld: &mut HashMap<(u64, u64), u32>,
    positions: &mut Vec<[f32; 3]>,
) -> u32 {
    let pa = pack(a[0], a[1], a[2]);
    let pb = pack(b[0], b[1], b[2]);
    let key = if pa < pb { (pa, pb) } else { (pb, pa) };
    if let Some(&vi) = weld.get(&key) {
        return vi;
    }
    let denom = fa - fb;
    // Clamp t strictly inside (0, 1) so the interpolated vertex never lands
    // exactly on a grid vertex. When `r_p / h` is an integer, whole shells of
    // the EDT field are exactly zero and every crossing would sit on a grid
    // vertex; coincident-but-distinct vertices then yield degenerate triangles
    // that get dropped, punching holes in the mesh. A tiny interior nudge keeps
    // vertices strictly on edge interiors — distinct, weldable, and invisible.
    let t = if denom.abs() < 1e-30 {
        0.5f64
    } else {
        (fa / denom).clamp(T_EPS, 1.0 - T_EPS) as f64
    };
    let p = [
        (origin[0] + (a[0] as f64 + t * (b[0] as f64 - a[0] as f64)) * h) as f32,
        (origin[1] + (a[1] as f64 + t * (b[1] as f64 - a[1] as f64)) * h) as f32,
        (origin[2] + (a[2] as f64 + t * (b[2] as f64 - a[2] as f64)) * h) as f32,
    ];
    let vi = positions.len() as u32;
    positions.push(p);
    weld.insert(key, vi);
    vi
}

/// Emit 0–2 triangles for one tetrahedron, wound outward (toward f > 0).
#[allow(clippy::too_many_arguments)]
fn emit_tet(
    tet: &[usize; 4],
    f8: &[f32; 8],
    g8: &[[usize; 3]; 8],
    origin: [f64; 3],
    h: f64,
    weld: &mut HashMap<(u64, u64), u32>,
    positions: &mut Vec<[f32; 3]>,
    indices: &mut Vec<[u32; 3]>,
) {
    let f = [f8[tet[0]], f8[tet[1]], f8[tet[2]], f8[tet[3]]];
    let g = [g8[tet[0]], g8[tet[1]], g8[tet[2]], g8[tet[3]]];

    // Inside = f < 0 (solid). Outside = f >= 0.
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
            let v0 = cross(a, out_idx[0], &f, &g, origin, h, weld, positions);
            let v1 = cross(a, out_idx[1], &f, &g, origin, h, weld, positions);
            let v2 = cross(a, out_idx[2], &f, &g, origin, h, weld, positions);
            // Outward = away from the inside vertex → use the centroid of the
            // three outside vertices as the reference point.
            let ref_pt = outside_ref(&g[out_idx[0]], &g[out_idx[1]], &g[out_idx[2]], origin, h);
            push_tri_oriented(v0, v1, v2, ref_pt, positions, indices);
        }
        3 => {
            let a = out_idx[0];
            let v0 = cross(in_idx[0], a, &f, &g, origin, h, weld, positions);
            let v1 = cross(in_idx[1], a, &f, &g, origin, h, weld, positions);
            let v2 = cross(in_idx[2], a, &f, &g, origin, h, weld, positions);
            // Outward = toward the single outside vertex.
            let ref_pt = world(g[a], origin, h);
            push_tri_oriented(v0, v1, v2, ref_pt, positions, indices);
        }
        2 => {
            let (i0, i1) = (in_idx[0], in_idx[1]);
            let (o0, o1) = (out_idx[0], out_idx[1]);
            let q0 = cross(i0, o0, &f, &g, origin, h, weld, positions);
            let q1 = cross(i1, o0, &f, &g, origin, h, weld, positions);
            let q2 = cross(i1, o1, &f, &g, origin, h, weld, positions);
            let q3 = cross(i0, o1, &f, &g, origin, h, weld, positions);
            // Outward = toward the outside half-space; midpoint of the two
            // outside vertices is a reliable anchor.
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

/// Crossing vertex between local tet corners `i` and `j`, oriented inside→outside.
#[allow(clippy::too_many_arguments)]
fn cross(
    i: usize,
    j: usize,
    f: &[f32; 4],
    g: &[[usize; 3]; 4],
    origin: [f64; 3],
    h: f64,
    weld: &mut HashMap<(u64, u64), u32>,
    positions: &mut Vec<[f32; 3]>,
) -> u32 {
    let (a, b) = if f[i] < 0.0 { (i, j) } else { (j, i) };
    edge_vertex(g[a], g[b], f[a], f[b], origin, h, weld, positions)
}

/// World position of a grid vertex.
fn world(g: [usize; 3], origin: [f64; 3], h: f64) -> [f32; 3] {
    [
        (origin[0] + g[0] as f64 * h) as f32,
        (origin[1] + g[1] as f64 * h) as f32,
        (origin[2] + g[2] as f64 * h) as f32,
    ]
}

/// Push a triangle wound so its normal points toward `outward_ref`
/// (a point on the exterior side of the surface).
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
    // Triangle centroid; the normal should point toward outward_ref.
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

    /// Build a signed-distance field of a sphere of radius `r` on a grid of
    /// spacing `h` with `pad` of margin, and polygonize it.
    fn sphere_mesh(r: f64, h: f64, pad: f64) -> Mesh {
        let extent = r + pad;
        let n = (2.0 * extent / h).ceil() as usize + 1;
        let mut field = ScalarField::new([n, n, n], [-extent; 3], h).unwrap();
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let p = field.position(x, y, z);
                    let d = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                    *field.get_mut(x, y, z).unwrap() = (d - r) as f32;
                }
            }
        }
        polygonize(&field)
    }

    #[test]
    fn sphere_area_converges() {
        let r = 1.0f64;
        let mesh = sphere_mesh(r, 0.25, 0.5);
        assert!(mesh.triangle_count() > 0, "empty mesh");
        let area = mesh.area();
        let expect = 4.0 * std::f64::consts::PI * r * r;
        let rel_err = (area - expect).abs() / expect;
        assert!(
            rel_err < 0.15,
            "area {area} vs analytic {expect}, rel_err {rel_err}"
        );
    }

    #[test]
    fn sphere_is_watertight() {
        let mesh = sphere_mesh(1.0, 0.3, 0.6);
        let report = crate::mesh::topology::audit(&mesh);
        assert!(report.is_watertight(), "{report:?}");
        assert_eq!(report.euler, 2, "{report:?}");
    }

    #[test]
    fn sphere_volume_converges() {
        let r = 1.0f64;
        let mesh = sphere_mesh(r, 0.2, 0.5);
        let vol = mesh.volume().abs();
        let expect = 4.0 / 3.0 * std::f64::consts::PI * r * r * r;
        let rel_err = (vol - expect).abs() / expect;
        assert!(rel_err < 0.15, "volume {vol} vs {expect}, err {rel_err}");
    }

    #[test]
    fn empty_field_gives_empty_mesh() {
        let field = ScalarField::new([4, 4, 4], [0.0; 3], 1.0).unwrap();
        let mesh = polygonize(&field);
        assert!(mesh.is_empty());
    }
}
