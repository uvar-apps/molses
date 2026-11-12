//! Painting the occupancy of `K = ⋃ B(c_i, r_i + r_p)` onto a grid.
//!
//! The paint pass also computes the exact exterior distance
//! `min_i (‖p − c_i‖ − R_i)`, which is the signed distance to `∂K` outside `K`.

use super::ScalarField;
use super::celllist::CellList;
use crate::structure::Atom;

/// Result of painting: occupancy of the inflated union and exact exterior distance.
#[derive(Debug)]
pub struct PaintResult {
    /// `true` when the voxel lies inside `K = ⋃ B(c_i, R_i)`.
    pub inside: Vec<bool>,
    /// Exact `min_i (‖p − c_i‖ − R_i)` at every voxel (negative inside some ball).
    pub exterior_dist: Vec<f32>,
}

/// Paint the inflated union of atomic balls onto the field grid.
///
/// `probe_radius` is `r_p`; each atom is inflated to `r_i + r_p`.
pub fn paint(atoms: &[Atom], probe_radius: f64, field: &ScalarField) -> PaintResult {
    let [nx, ny, nz] = field.dims();
    let n = nx * ny * nz;
    let mut inside = vec![false; n];
    let mut exterior_dist = vec![f32::INFINITY; n];

    let max_r = atoms
        .iter()
        .map(|a| a.radius + probe_radius)
        .fold(0.0f64, f64::max);

    let celllist = CellList::new(
        atoms,
        max_r + field.spacing(),
        field.origin(),
        field.dims(),
        field.spacing(),
    );

    let spacing = field.spacing();
    let origin = field.origin();
    let slab = nx * ny;

    // Iterate grid points; query the cell list. This is the SIMD-friendly shape
    // (contiguous x-rows) and is what the `simd` feature accelerates.
    //
    // Correctness: a sphere missing from a query of radius Q can only beat the
    // current best `b` if its center lies within `b + max_r` of p. We therefore
    // expand the query radius until it covers `best + max_r`, which makes the
    // result exact. Far-field points (no atoms in range) expand geometrically
    // until the nearest atom is found.
    //
    // Parallelism: z-slabs are disjoint slices of the output, so
    // `par_chunks_mut` fans out safely. Determinism is preserved: each (x,y,z)
    // writes exactly one slot, so results are bitwise-identical for any thread
    // count.

    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        inside
            .par_chunks_mut(slab)
            .zip(exterior_dist.par_chunks_mut(slab))
            .enumerate()
            .for_each(|(z, (inside_z, dist_z))| {
                paint_slab(
                    z,
                    nx,
                    ny,
                    atoms,
                    probe_radius,
                    max_r,
                    spacing,
                    origin,
                    &celllist,
                    inside_z,
                    dist_z,
                );
            });
    }
    #[cfg(not(feature = "parallel"))]
    {
        for z in 0..nz {
            let (inside_z, rest) = inside.split_at_mut((z + 1) * slab);
            let inside_z = &mut inside_z[z * slab..];
            let (dist_z, _) = exterior_dist.split_at_mut((z + 1) * slab);
            let dist_z = &mut dist_z[z * slab..];
            paint_slab(
                z,
                nx,
                ny,
                &atoms,
                probe_radius,
                max_r,
                spacing,
                origin,
                &celllist,
                inside_z,
                dist_z,
            );
        }
    }

    PaintResult {
        inside,
        exterior_dist,
    }
}

/// Paint one z-slab of the grid into `inside_z` / `dist_z` (slices of length
/// `nx * ny` for a single z).
#[allow(clippy::too_many_arguments)]
fn paint_slab(
    z: usize,
    nx: usize,
    ny: usize,
    atoms: &[Atom],
    probe_radius: f64,
    max_r: f64,
    spacing: f64,
    origin: [f64; 3],
    celllist: &CellList,
    inside_z: &mut [bool],
    dist_z: &mut [f32],
) {
    let mut neighbors = Vec::new();
    let mut nb_x: Vec<f32> = Vec::with_capacity(64);
    let mut nb_y: Vec<f32> = Vec::with_capacity(64);
    let mut nb_z: Vec<f32> = Vec::with_capacity(64);
    let mut nb_r: Vec<f32> = Vec::with_capacity(64);
    for y in 0..ny {
        for x in 0..nx {
            let p = [
                origin[0] + x as f64 * spacing,
                origin[1] + y as f64 * spacing,
                origin[2] + z as f64 * spacing,
            ];
            let mut best = f64::INFINITY;
            let mut any_inside = false;
            let mut q = max_r + spacing;
            for _ in 0..64 {
                celllist.query(p, q, &mut neighbors);
                best = f64::INFINITY;
                any_inside = false;
                nb_x.clear();
                nb_y.clear();
                nb_z.clear();
                nb_r.clear();
                for &ai in &neighbors {
                    let a = &atoms[ai as usize];
                    nb_x.push(a.position[0] as f32);
                    nb_y.push(a.position[1] as f32);
                    nb_z.push(a.position[2] as f32);
                    nb_r.push((a.radius + probe_radius) as f32);
                }
                if !nb_x.is_empty() {
                    let (b, inside) = crate::field::simd::min_distance(
                        p[0], p[1], p[2], &nb_x, &nb_y, &nb_z, &nb_r,
                    );
                    best = b;
                    any_inside = inside;
                }
                if neighbors.is_empty() {
                    q *= 2.0;
                    if q > 1.0e9 {
                        break;
                    }
                    continue;
                }
                if q >= best + max_r {
                    break;
                }
                q = (best + max_r).max(q * 2.0);
            }
            let i = y * nx + x;
            inside_z[i] = any_inside;
            dist_z[i] = if best.is_finite() {
                best as f32
            } else {
                f32::MAX / 2.0
            };
        }
    }
}

/// Assemble the signed distance `f_sas` from occupancy + EDT + exact exterior distances.
///
/// - Free (exterior) voxels: exact `min_i (‖p−c_i‖ − R_i)`.
/// - Solid voxels (inside `K` or enclosed cavities): `−sqrt(edt_squared)`.
pub fn assemble_signed_distance(
    paint: &PaintResult,
    edt_sq: &[f64],
    spacing: f64,
    out: &mut ScalarField,
) {
    let data = out.data_mut();
    for i in 0..data.len() {
        data[i] = if paint.inside[i] {
            // Distance to nearest free voxel, in Å (EDT is in voxel units).
            -(edt_sq[i].max(0.0).sqrt() as f32) * spacing as f32
        } else {
            paint.exterior_dist[i]
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{Atom, Element};

    #[test]
    fn paint_single_sphere() {
        let atoms = vec![Atom::new([0.0; 3], 1.0, Some(Element::C), 1)];
        let field = ScalarField::new([15, 15, 15], [-7.0; 3], 1.0).unwrap();
        let probe = 0.5;
        let pr = paint(&atoms, probe, &field);
        // Center of the grid is at origin: inside the inflated sphere (r=1.5).
        let c = 7 + 7 * 15 + 7 * 15 * 15;
        assert!(pr.inside[c]);
        // A far corner is outside.
        let corner = 0;
        assert!(!pr.inside[corner]);
        // Exterior distance at the corner ≈ distance to sphere.
        let d = pr.exterior_dist[corner];
        assert!(d > 0.0);
    }

    #[test]
    fn signed_distance_outside_is_exact() {
        let atoms = vec![Atom::new([0.0; 3], 1.0, Some(Element::C), 1)];
        let field = ScalarField::new([21, 21, 21], [-10.0; 3], 1.0).unwrap();
        let probe = 0.5; // R = 1.5
        let pr = paint(&atoms, probe, &field);
        // EDT of the solid region to the nearest free voxel.
        let edt = crate::field::edt::edt_squared(&pr.inside, field.dims());
        let mut f = ScalarField::new([21, 21, 21], [-10.0; 3], 1.0).unwrap();
        assemble_signed_distance(&pr, &edt, 1.0, &mut f);

        // Check a point outside: grid point at (0, 0, 0) world = grid (10,10,10)
        // is the center (inside). Grid point at world (5,0,0) = grid (15,10,10).
        let i = 15 + 10 * 21 + 10 * 21 * 21;
        assert!(!pr.inside[i]);
        let f_val = f.data()[i];
        // exact: |p| - R = 5 - 1.5 = 3.5
        assert!((f_val - 3.5).abs() < 1e-4, "got {f_val}");
    }
}
