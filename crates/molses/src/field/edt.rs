//! Exact Euclidean distance transform (Felzenszwalb & Huttenlocher 2012).
//!
//! Computes, for every voxel, the squared Euclidean distance to the nearest
//! voxel of a specified set. The algorithm is separable: three 1-D passes.

/// 1-D squared distance transform of a sampled function.
///
/// `f[i]` is the cost at sample `i`; returns `d[i] = min_j (f[j] + (i-j)²)`.
/// The lower envelope of parabolas is built in O(n).
pub fn edt_1d(f: &[f64]) -> Vec<f64> {
    let n = f.len();
    if n == 0 {
        return Vec::new();
    }
    let mut d = vec![0.0f64; n];
    let mut v = vec![0usize; n]; // locations of parabolas in the envelope
    let mut z = vec![0.0f64; n + 1]; // boundaries between parabolas
    let mut k = 0usize;
    v[0] = 0;
    z[0] = f64::NEG_INFINITY;
    z[1] = f64::INFINITY;

    for q in 1..n {
        // Intersection of parabola q with the last parabola in the envelope.
        let mut s = {
            let p = v[k];
            ((f[q] + (q * q) as f64) - (f[p] + (p * p) as f64)) / (2.0 * q as f64 - 2.0 * p as f64)
        };
        while s <= z[k] {
            if k == 0 {
                break;
            }
            k -= 1;
            let p = v[k];
            s = ((f[q] + (q * q) as f64) - (f[p] + (p * p) as f64))
                / (2.0 * q as f64 - 2.0 * p as f64);
        }
        k += 1;
        v[k] = q;
        z[k] = s;
        z[k + 1] = f64::INFINITY;
    }

    k = 0;
    for (q, dq) in d.iter_mut().enumerate() {
        while z[k + 1] < q as f64 {
            k += 1;
        }
        let p = v[k];
        let diff = q as f64 - p as f64;
        *dq = diff * diff + f[p];
    }
    d
}

/// 3-D exact squared EDT over a boolean occupancy grid.
///
/// `solid[i] == true` means the voxel is in the set whose complement we measure
/// distance to... more precisely: returns `dist2[i]` = squared distance from
/// voxel `i` to the **nearest non-solid** voxel. Non-solid voxels get 0.
pub fn edt_squared(solid: &[bool], dims: [usize; 3]) -> Vec<f64> {
    let [nx, ny, nz] = dims;
    let n = nx * ny * nz;
    assert_eq!(solid.len(), n, "occupancy length must equal nx*ny*nz");
    const INF: f64 = 1.0e12;

    // Initialize: 0 at free voxels, INF at solid voxels.
    let mut f: Vec<f64> = solid.iter().map(|&s| if s { INF } else { 0.0 }).collect();

    let idx = |x: usize, y: usize, z: usize| x + y * nx + z * nx * ny;

    // Pass 1: along x (rows).
    let mut line = vec![0.0f64; nx];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                line[x] = f[idx(x, y, z)];
            }
            let d = edt_1d(&line);
            for x in 0..nx {
                f[idx(x, y, z)] = d[x];
            }
        }
    }

    // Pass 2: along y (columns).
    let mut line = vec![0.0f64; ny];
    for z in 0..nz {
        for x in 0..nx {
            for y in 0..ny {
                line[y] = f[idx(x, y, z)];
            }
            let d = edt_1d(&line);
            for y in 0..ny {
                f[idx(x, y, z)] = d[y];
            }
        }
    }

    // Pass 3: along z (pillars).
    let mut line = vec![0.0f64; nz];
    for y in 0..ny {
        for x in 0..nx {
            for z in 0..nz {
                line[z] = f[idx(x, y, z)];
            }
            let d = edt_1d(&line);
            for z in 0..nz {
                f[idx(x, y, z)] = d[z];
            }
        }
    }

    f
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn edt_1d_single_seed() {
        // Free at index 3, cost 0; INF elsewhere.
        let f = vec![1.0e12, 1.0e12, 1.0e12, 0.0, 1.0e12, 1.0e12];
        let d = edt_1d(&f);
        assert!((d[3] - 0.0).abs() < 1e-9);
        assert!((d[2] - 1.0).abs() < 1e-9);
        assert!((d[4] - 1.0).abs() < 1e-9);
        assert!((d[0] - 9.0).abs() < 1e-9);
        assert!((d[5] - 4.0).abs() < 1e-9);
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn edt_1d_matches_brute_force() {
        let f = vec![5.0, 0.0, 3.0, 1.0e12, 1.0e12, 0.0, 2.0];
        let d = edt_1d(&f);
        for i in 0..f.len() {
            let mut best = f64::INFINITY;
            for j in 0..f.len() {
                let cost = f[j] + (i as f64 - j as f64).powi(2);
                best = best.min(cost);
            }
            assert!((d[i] - best).abs() < 1e-9, "i={i}: {} vs {}", d[i], best);
        }
    }

    #[test]
    fn edt_3d_matches_brute_force() {
        let dims = [4, 3, 2];
        let n = dims[0] * dims[1] * dims[2];
        // A couple of free voxels scattered around.
        let mut solid = vec![true; n];
        solid[0] = false; // (0,0,0)
        solid[n - 1] = false; // last
        solid[5] = false;

        let d = edt_squared(&solid, dims);
        let idx = |x: usize, y: usize, z: usize| x + y * dims[0] + z * dims[0] * dims[1];

        for z in 0..dims[2] {
            for y in 0..dims[1] {
                for x in 0..dims[0] {
                    let i = idx(x, y, z);
                    let mut best = f64::INFINITY;
                    for zz in 0..dims[2] {
                        for yy in 0..dims[1] {
                            for xx in 0..dims[0] {
                                let j = idx(xx, yy, zz);
                                if solid[j] {
                                    continue;
                                }
                                let dx = x as f64 - xx as f64;
                                let dy = y as f64 - yy as f64;
                                let dz = z as f64 - zz as f64;
                                best = best.min(dx * dx + dy * dy + dz * dz);
                            }
                        }
                    }
                    let expect = if solid[i] { best } else { 0.0 };
                    assert!(
                        (d[i] - expect).abs() < 1e-9,
                        "voxel {i}: got {}, want {}",
                        d[i],
                        expect
                    );
                }
            }
        }
    }
}
