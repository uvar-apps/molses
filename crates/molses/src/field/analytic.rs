//! Analytic signed distance to the SAS (boundary of the inflated union).
//!
//! Replaces the EDT-based field with a smooth computation free of voxel
//! quantization. The SES is the zero level set of `f_ses = f_sas + r_p`
//! where `f_sas` is the signed distance to `∂K`, `K = ⋃ B(c_i, r_i+r_p)`.
//!
//! - **Outside `K`**: `f_sas = min_i (‖p−c_i‖ − R_i)` (exact).
//! - **Inside `K`**: `f_sas = −d(p, ∂K)`. The nearest visible boundary point
//!   is either the radial projection onto some sphere (when not buried in a
//!   neighbor) or a point on an intersection circle of two inflated spheres
//!   (the crevice case — radial projections are buried there).
//!
//! Neighbor queries use a [`CellList`](super::celllist::CellList) so cost is
//! O(local) per grid point.

use super::ScalarField;
use super::celllist::CellList;
use crate::structure::Atom;

/// Compute `f_ses = f_sas + r_p` on every grid point of `field`.
pub fn build_ses_field(atoms: &[Atom], probe_radius: f64, field: &mut ScalarField) {
    let origin = field.origin();
    let h = field.spacing();
    let [nx, ny, nz] = field.dims();

    let radii: Vec<f64> = atoms.iter().map(|a| a.radius + probe_radius).collect();
    let r_max = radii.iter().cloned().fold(0.0f64, f64::max);

    // Query radius covers: inflated sphere + a little slack for pair circles.
    let query_r = 2.0 * r_max + 2.0 * probe_radius + h;
    let cells = CellList::new(atoms, query_r, origin, [nx, ny, nz], h);

    let mut near: Vec<u32> = Vec::with_capacity(64);
    let mut near_q: Vec<u32> = Vec::with_capacity(64);
    let data = field.data_mut();
    let slab = nx * ny;

    #[cfg(feature = "parallel")]
    {
        use rayon::prelude::*;
        data.par_chunks_mut(slab)
            .enumerate()
            .for_each(|(z, slab_data)| {
                let mut near = Vec::with_capacity(64);
                let mut near_q = Vec::with_capacity(64);
                for y in 0..ny {
                    for x in 0..nx {
                        let p = [
                            origin[0] + x as f64 * h,
                            origin[1] + y as f64 * h,
                            origin[2] + z as f64 * h,
                        ];
                        cells.query(p, query_r, &mut near);
                        let val =
                            point_ses(p, atoms, &radii, &near, &cells, &mut near_q, probe_radius);
                        slab_data[x + y * nx] = val as f32;
                    }
                }
            });
        let _ = (&mut near, &mut near_q);
    }
    #[cfg(not(feature = "parallel"))]
    {
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    let p = [
                        origin[0] + x as f64 * h,
                        origin[1] + y as f64 * h,
                        origin[2] + z as f64 * h,
                    ];
                    cells.query(p, query_r, &mut near);
                    let val = point_ses(p, atoms, &radii, &near, &cells, &mut near_q, probe_radius);
                    data[x + y * nx + z * slab] = val as f32;
                }
            }
        }
    }
}

/// `f_ses(p) = f_sas(p) + r_p`.
fn point_ses(
    p: [f64; 3],
    atoms: &[Atom],
    radii: &[f64],
    near: &[u32],
    cells: &CellList,
    near_q: &mut Vec<u32>,
    probe_radius: f64,
) -> f64 {
    let f_sas = point_sas(p, atoms, radii, near, cells, near_q);
    f_sas + probe_radius
}

/// Signed distance from `p` to `∂K` (negative inside).
fn point_sas(
    p: [f64; 3],
    atoms: &[Atom],
    radii: &[f64],
    near: &[u32],
    cells: &CellList,
    near_q: &mut Vec<u32>,
) -> f64 {
    if near.is_empty() {
        return 1e3;
    }

    // --- Outside test: min signed distance to inflated spheres ---
    let mut best_out = f64::INFINITY;
    let mut best_i = 0usize;
    let mut best_dc = 0.0f64;
    for &i in near {
        let i = i as usize;
        let a = &atoms[i];
        let dx = p[0] - a.position[0];
        let dy = p[1] - a.position[1];
        let dz = p[2] - a.position[2];
        let dc = (dx * dx + dy * dy + dz * dz).sqrt();
        let signed = dc - radii[i];
        if signed < best_out {
            best_out = signed;
            best_i = i;
            best_dc = dc;
        }
    }
    if best_out >= 0.0 {
        return best_out;
    }

    // Fast path: if even the *closest* sphere surface is nearer than r_p,
    // then depth ≤ -best_out < r_p, so f_ses = -depth + r_p > 0.
    // (Only valid for the SES offset caller — here we return f_sas, and
    // the caller adds r_p. Depth ≤ -best_out; if -best_out < r_p then
    // depth < r_p ⇒ f_ses > 0. We still need f_sas itself, so compute
    // exact depth when it could matter: depth near r_p.)
    // depth ≤ -best_out always. If -best_out is huge (deep interior),
    // sign of f_ses is negative as long as depth > r_p. We don't know
    // depth ≥ ?, but a safe sufficient condition: -best_out > r_p does NOT
    // prove depth > r_p. So always compute exact depth for correctness.

    // --- Inside K: distance to nearest visible boundary ---
    let mut best_depth = f64::INFINITY;
    let mut found = false;

    // Visibility of a candidate boundary point q: re-query cells at q so we
    // see atoms that neighbor q but not p.
    let mut vis_q = |q: [f64; 3], host: usize| -> bool {
        cells.query(q, 1.0, near_q);
        for &j in near_q.iter() {
            let j = j as usize;
            if j == host {
                continue;
            }
            let b = &atoms[j];
            let ex = q[0] - b.position[0];
            let ey = q[1] - b.position[1];
            let ez = q[2] - b.position[2];
            if (ex * ex + ey * ey + ez * ez).sqrt() < radii[j] - 1e-9 {
                return false;
            }
        }
        true
    };

    if best_dc > 1e-12 {
        let a = &atoms[best_i];
        let r = radii[best_i];
        let inv = r / best_dc;
        let x = a.position[0] + (p[0] - a.position[0]) * inv;
        let y = a.position[1] + (p[1] - a.position[1]) * inv;
        let z = a.position[2] + (p[2] - a.position[2]) * inv;
        if vis_q([x, y, z], best_i) {
            best_depth = r - best_dc;
            found = true;
        }
    } else {
        best_depth = radii[best_i];
        found = true;
    }

    for &i in near {
        let i = i as usize;
        if i == best_i {
            continue;
        }
        let a = &atoms[i];
        let dx = p[0] - a.position[0];
        let dy = p[1] - a.position[1];
        let dz = p[2] - a.position[2];
        let dc = (dx * dx + dy * dy + dz * dz).sqrt();
        let r = radii[i];
        if dc < 1e-12 || dc >= r {
            continue;
        }
        let inv = r / dc;
        let q = [
            a.position[0] + dx * inv,
            a.position[1] + dy * inv,
            a.position[2] + dz * inv,
        ];
        if vis_q(q, i) {
            let depth = r - dc;
            if depth < best_depth {
                best_depth = depth;
                found = true;
            }
        }
    }

    // Intersection circles — only if their nearest point is visible on ∂K.
    for (ai, &ii) in near.iter().enumerate() {
        let i = ii as usize;
        for &jj in near.iter().skip(ai + 1) {
            let j = jj as usize;
            let ri = radii[i];
            let rj = radii[j];
            let ci = atoms[i].position;
            let cj = atoms[j].position;
            let dxc = cj[0] - ci[0];
            let dyc = cj[1] - ci[1];
            let dzc = cj[2] - ci[2];
            let dij = (dxc * dxc + dyc * dyc + dzc * dzc).sqrt();
            if dij < 1e-12 || dij >= ri + rj || dij <= (ri - rj).abs() {
                continue;
            }
            let t = (ri * ri - rj * rj + dij * dij) / (2.0 * dij);
            let rc2 = ri * ri - t * t;
            if rc2 <= 0.0 {
                continue;
            }
            let rc = rc2.sqrt();
            let inv_t = t / dij;
            let c = [
                ci[0] + inv_t * dxc,
                ci[1] + inv_t * dyc,
                ci[2] + inv_t * dzc,
            ];
            let n = [dxc / dij, dyc / dij, dzc / dij];
            let (nearest, dist) = nearest_on_circle(p, c, rc, n);
            // Skip clearly farther candidates before the visibility query.
            if dist >= best_depth {
                continue;
            }
            if !vis_q(nearest, i) {
                continue;
            }
            if dist < best_depth {
                best_depth = dist;
                found = true;
            }
        }
    }

    if found {
        return -best_depth;
    }
    -best_out.abs()
}

/// Nearest point on a circle to `p`, and the distance.
fn nearest_on_circle(p: [f64; 3], c: [f64; 3], r: f64, n: [f64; 3]) -> ([f64; 3], f64) {
    let pc = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
    let d_plane = pc[0] * n[0] + pc[1] * n[1] + pc[2] * n[2];
    let rx = pc[0] - d_plane * n[0];
    let ry = pc[1] - d_plane * n[1];
    let rz = pc[2] - d_plane * n[2];
    let rho = (rx * rx + ry * ry + rz * rz).sqrt();
    if rho < 1e-12 {
        // Project to any circle point in-plane: pick +x orthogonal to n.
        let mut t = [1.0, 0.0, 0.0];
        if n[0].abs() > 0.9 {
            t = [0.0, 1.0, 0.0];
        }
        // t -= (t·n)n, normalize
        let dot = t[0] * n[0] + t[1] * n[1] + t[2] * n[2];
        t = [t[0] - dot * n[0], t[1] - dot * n[1], t[2] - dot * n[2]];
        let tl = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
        let nearest = [
            c[0] + r * t[0] / tl,
            c[1] + r * t[1] / tl,
            c[2] + r * t[2] / tl,
        ];
        let dx = p[0] - nearest[0];
        let dy = p[1] - nearest[1];
        let dz = p[2] - nearest[2];
        return (nearest, (dx * dx + dy * dy + dz * dz).sqrt());
    }
    let inv = r / rho;
    let nearest = [c[0] + rx * inv, c[1] + ry * inv, c[2] + rz * inv];
    let dx = p[0] - nearest[0];
    let dy = p[1] - nearest[1];
    let dz = p[2] - nearest[2];
    (nearest, (dx * dx + dy * dy + dz * dz).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::ScalarField;
    use crate::field::celllist::CellList;
    use crate::structure::Element;

    fn ses_at(p: [f64; 3], atoms: &[Atom], radii: &[f64], cells: &CellList, rp: f64) -> f64 {
        let mut near = Vec::new();
        let mut near_q = Vec::new();
        let qr = 2.0 * radii.iter().cloned().fold(0.0f64, f64::max) + 2.0 * rp + 1.0;
        cells.query(p, qr, &mut near);
        point_ses(p, atoms, radii, &near, cells, &mut near_q, rp)
    }

    #[test]
    fn single_sphere_ses_is_vdw() {
        let atoms = vec![Atom::new([0.0; 3], 1.7, Some(Element::C), 1)];
        let radii = vec![1.7 + 1.4];
        let mut field = ScalarField::new([21, 21, 21], [-3.0; 3], 0.3).unwrap();
        build_ses_field(&atoms, 1.4, &mut field);
        let cells = CellList::new(&atoms, 10.0, field.origin(), field.dims(), field.spacing());
        let v = ses_at([1.7, 0.0, 0.0], &atoms, &radii, &cells, 1.4);
        assert!(v.abs() < 0.05, "f_ses on VdW sphere = {v}");
        assert!(ses_at([5.0, 0.0, 0.0], &atoms, &radii, &cells, 1.4) > 0.0);
        assert!(ses_at([0.0, 0.0, 0.0], &atoms, &radii, &cells, 1.4) < 0.0);
    }

    #[test]
    fn two_atom_toric_point_is_zero() {
        let r = 1.7;
        let rp = 1.4;
        let r_inf = r + rp;
        let y_q = (r_inf * r_inf - 1.7f64 * 1.7).sqrt();
        let p = [0.0, y_q - rp, 0.0];
        let atoms = vec![
            Atom::new([-1.7, 0.0, 0.0], r, Some(Element::C), 1),
            Atom::new([1.7, 0.0, 0.0], r, Some(Element::C), 2),
        ];
        let radii = vec![r_inf, r_inf];
        let mut field = ScalarField::new([41, 41, 41], [-5.0; 3], 0.25).unwrap();
        build_ses_field(&atoms, rp, &mut field);
        let cells = CellList::new(&atoms, 12.0, field.origin(), field.dims(), field.spacing());
        let v = ses_at(p, &atoms, &radii, &cells, rp);
        assert!(
            v.abs() < 0.08,
            "f_ses at toric point ({:.3}, {:.3}, 0) = {v}, expected 0",
            p[0],
            p[1]
        );
    }

    #[test]
    fn protein_atom_center_is_inside() {
        let pdb = format!("{}/../../testdata/1crn.pdb", env!("CARGO_MANIFEST_DIR"));
        let (mut structure, mut report) =
            crate::parse::parse_file(std::path::Path::new(&pdb)).unwrap();
        crate::parse::assign_radii(
            &mut structure,
            &crate::radii::RadiiSet::Bondi,
            false,
            &mut report,
        )
        .unwrap();
        crate::parse::deduplicate(&mut structure, &mut report);
        let probe = 1.4;
        let h = 0.5;
        let (lo, hi) = structure.bounding_box();
        let r_max = structure
            .atoms()
            .iter()
            .map(|a| a.radius)
            .fold(0.0f64, f64::max);
        let pad = 2.0 * probe + r_max + 3.0 * h;
        let origin = [lo[0] - pad, lo[1] - pad, lo[2] - pad];
        let extent = [
            hi[0] + pad - origin[0],
            hi[1] + pad - origin[1],
            hi[2] + pad - origin[2],
        ];
        let dims = [
            (extent[0] / h).ceil() as usize + 1,
            (extent[1] / h).ceil() as usize + 1,
            (extent[2] / h).ceil() as usize + 1,
        ];
        let mut field = ScalarField::new(dims, origin, h).unwrap();
        build_ses_field(structure.atoms(), probe, &mut field);
        let neg = field.data().iter().filter(|&&v| v < 0.0).count();
        assert!(neg > 1000, "too few interior samples: {neg}");

        let radii: Vec<f64> = structure.atoms().iter().map(|a| a.radius + probe).collect();
        let cells = CellList::new(
            structure.atoms(),
            2.0 * r_max + 2.0 * probe + h,
            field.origin(),
            field.dims(),
            field.spacing(),
        );
        let a0 = structure.atoms()[0].position;
        let v = ses_at(a0, structure.atoms(), &radii, &cells, probe);
        assert!(v < 0.0, "atom center should be inside SES solid, f={v}");
    }
}
