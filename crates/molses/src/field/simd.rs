//! SIMD-accelerated distance kernels for the paint pass.
//!
//! The hot loop computes `min_i (‖p − c_i‖ − R_i)` over a set of neighbor
//! atoms. With the `simd` feature this processes 4 atoms per iteration using
//! `wide::f32x4`; the scalar fallback processes one at a time.

/// Compute `min_i (‖p − c_i‖ − R_i)` and `any_i (‖p − c_i‖ ≤ R_i)` over the
/// given neighbor atoms. Positions and radii are read from the provided slices
/// (already gathered by the caller).
///
/// Returns `(best_distance, any_inside)`.
pub fn min_distance(
    px: f64,
    py: f64,
    pz: f64,
    atom_x: &[f32],
    atom_y: &[f32],
    atom_z: &[f32],
    atom_r: &[f32],
) -> (f64, bool) {
    debug_assert_eq!(atom_x.len(), atom_y.len());
    debug_assert_eq!(atom_x.len(), atom_z.len());
    debug_assert_eq!(atom_x.len(), atom_r.len());
    let n = atom_x.len();
    let mut best = f64::INFINITY;
    let mut any_inside = false;

    #[cfg(feature = "simd")]
    {
        use wide::f32x4;
        let px4 = f32x4::splat(px as f32);
        let py4 = f32x4::splat(py as f32);
        let pz4 = f32x4::splat(pz as f32);

        let mut i = 0;
        // Process 4 atoms per iteration.
        while i + 4 <= n {
            let ax = f32x4::new([atom_x[i], atom_x[i + 1], atom_x[i + 2], atom_x[i + 3]]);
            let ay = f32x4::new([atom_y[i], atom_y[i + 1], atom_y[i + 2], atom_y[i + 3]]);
            let az = f32x4::new([atom_z[i], atom_z[i + 1], atom_z[i + 2], atom_z[i + 3]]);
            let ar = f32x4::new([atom_r[i], atom_r[i + 1], atom_r[i + 2], atom_r[i + 3]]);

            let dx = px4 - ax;
            let dy = py4 - ay;
            let dz = pz4 - az;
            let d2 = dx * dx + dy * dy + dz * dz;
            let d = d2.sqrt() - ar; // f32x4

            // Horizontal reduce: min of d, and any d <= 0.
            let arr = d.to_array();
            for &v in &arr {
                let dv = v as f64;
                if dv < best {
                    best = dv;
                }
                if dv <= 0.0 {
                    any_inside = true;
                }
            }
            i += 4;
        }
        // Scalar tail.
        while i < n {
            let dx = px as f32 - atom_x[i];
            let dy = py as f32 - atom_y[i];
            let dz = pz as f32 - atom_z[i];
            let dv = ((dx * dx + dy * dy + dz * dz) as f64).sqrt() - atom_r[i] as f64;
            if dv < best {
                best = dv;
            }
            if dv <= 0.0 {
                any_inside = true;
            }
            i += 1;
        }
    }

    #[cfg(not(feature = "simd"))]
    {
        for i in 0..n {
            let dx = px - atom_x[i] as f64;
            let dy = py - atom_y[i] as f64;
            let dz = pz - atom_z[i] as f64;
            let dv = (dx * dx + dy * dy + dz * dz).sqrt() - atom_r[i] as f64;
            if dv < best {
                best = dv;
            }
            if dv <= 0.0 {
                any_inside = true;
            }
        }
    }

    (best, any_inside)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn min_distance_single() {
        let (d, inside) = min_distance(0.0, 0.0, 0.0, &[3.0], &[0.0], &[0.0], &[1.5]);
        assert!((d - 1.5).abs() < 1e-5, "d = {d}");
        assert!(!inside);
    }

    #[test]
    fn min_distance_inside() {
        let (d, inside) = min_distance(0.5, 0.0, 0.0, &[0.0], &[0.0], &[0.0], &[1.5]);
        assert!(d < 0.0);
        assert!(inside);
    }

    #[test]
    fn min_distance_many() {
        // 5 atoms: check SIMD + scalar tail agree.
        let xs = [1.0f32, 2.0, 3.0, 4.0, 5.0];
        let ys = [0.0f32; 5];
        let zs = [0.0f32; 5];
        let rs = [1.0f32; 5];
        let (d, _) = min_distance(0.0, 0.0, 0.0, &xs, &ys, &zs, &rs);
        // Nearest atom at x=1, r=1 → d = 0.
        assert!(d.abs() < 1e-5, "d = {d}");
    }
}
