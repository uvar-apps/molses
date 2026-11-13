//! Scalar field construction for the SES engine (M2).
//!
//! Pipeline: cell list → occupancy painting of `K = ⋃ B(c_i, r_i+r_p)` →
//! exact Euclidean distance transform → signed distance `f_sas` → offset
//! `f_ses = f_sas + r_p`.

pub mod analytic;
pub mod celllist;
pub mod edt;
pub mod paint;
pub mod simd;

use crate::error::{Result, SesError};

/// A dense scalar field sampled on a regular grid.
///
/// Values are `f32` in Ångström. Indexing is x-fastest:
/// `idx = x + y * dims[0] + z * dims[0] * dims[1]`.
#[derive(Debug, Clone)]
pub struct ScalarField {
    dims: [usize; 3],
    origin: [f64; 3],
    spacing: f64,
    data: Vec<f32>,
}

impl ScalarField {
    /// Allocate a zeroed field.
    pub fn new(dims: [usize; 3], origin: [f64; 3], spacing: f64) -> Result<Self> {
        if spacing <= 0.0 || !spacing.is_finite() {
            return Err(SesError::invalid_param(
                "grid_spacing",
                spacing,
                "must be finite and > 0",
            ));
        }
        let total = dims
            .iter()
            .try_fold(1usize, |acc, &d| acc.checked_mul(d))
            .ok_or_else(|| SesError::Geometry {
                stage: "field",
                msg: "grid dimensions overflow".into(),
            })?;
        Ok(Self {
            dims,
            origin,
            spacing,
            data: vec![0.0; total],
        })
    }

    /// Grid dimensions `[nx, ny, nz]`.
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    /// World-space origin of grid point `(0, 0, 0)` in Å.
    pub fn origin(&self) -> [f64; 3] {
        self.origin
    }

    /// Grid spacing `h` in Å.
    pub fn spacing(&self) -> f64 {
        self.spacing
    }

    /// Total number of samples.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// `true` when the field has no samples.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Linear index of grid coordinates.
    pub fn index(&self, x: usize, y: usize, z: usize) -> usize {
        x + y * self.dims[0] + z * self.dims[0] * self.dims[1]
    }

    /// Sample at grid coordinates (bounds-checked).
    pub fn get(&self, x: usize, y: usize, z: usize) -> Option<f32> {
        if x < self.dims[0] && y < self.dims[1] && z < self.dims[2] {
            Some(self.data[self.index(x, y, z)])
        } else {
            None
        }
    }

    /// Mutable sample at grid coordinates.
    pub fn get_mut(&mut self, x: usize, y: usize, z: usize) -> Option<&mut f32> {
        if x < self.dims[0] && y < self.dims[1] && z < self.dims[2] {
            let i = self.index(x, y, z);
            Some(&mut self.data[i])
        } else {
            None
        }
    }

    /// World position of grid point `(x, y, z)`.
    pub fn position(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        [
            self.origin[0] + x as f64 * self.spacing,
            self.origin[1] + y as f64 * self.spacing,
            self.origin[2] + z as f64 * self.spacing,
        ]
    }

    /// Raw data slice.
    pub fn data(&self) -> &[f32] {
        &self.data
    }

    /// Raw mutable data slice.
    pub fn data_mut(&mut self) -> &mut [f32] {
        &mut self.data
    }

    /// Trilinear sample at a world position. Returns `None` outside the grid.
    pub fn sample(&self, p: [f64; 3]) -> Option<f32> {
        let gx = (p[0] - self.origin[0]) / self.spacing;
        let gy = (p[1] - self.origin[1]) / self.spacing;
        let gz = (p[2] - self.origin[2]) / self.spacing;
        if gx < 0.0 || gy < 0.0 || gz < 0.0 {
            return None;
        }
        let x0 = gx.floor() as usize;
        let y0 = gy.floor() as usize;
        let z0 = gz.floor() as usize;
        if x0 + 1 >= self.dims[0] || y0 + 1 >= self.dims[1] || z0 + 1 >= self.dims[2] {
            return None;
        }
        let fx = (gx - x0 as f64) as f32;
        let fy = (gy - y0 as f64) as f32;
        let fz = (gz - z0 as f64) as f32;
        let mut acc = 0.0f32;
        for dz in 0..2 {
            for dy in 0..2 {
                for dx in 0..2 {
                    let w = match (dx, dy, dz) {
                        (0, 0, 0) => (1.0 - fx) * (1.0 - fy) * (1.0 - fz),
                        (1, 0, 0) => fx * (1.0 - fy) * (1.0 - fz),
                        (0, 1, 0) => (1.0 - fx) * fy * (1.0 - fz),
                        (1, 1, 0) => fx * fy * (1.0 - fz),
                        (0, 0, 1) => (1.0 - fx) * (1.0 - fy) * fz,
                        (1, 0, 1) => fx * (1.0 - fy) * fz,
                        (0, 1, 1) => (1.0 - fx) * fy * fz,
                        (1, 1, 1) => fx * fy * fz,
                        _ => unreachable!(),
                    };
                    let v = self.get(x0 + dx, y0 + dy, z0 + dz)?;
                    acc += w * v;
                }
            }
        }
        Some(acc)
    }

    /// Central-difference gradient at a world position (Å per unit, ≈ unit for an SDF).
    pub fn gradient(&self, p: [f64; 3]) -> Option<[f32; 3]> {
        let h = self.spacing;
        let mut g = [0.0f32; 3];
        for (k, gp) in g.iter_mut().enumerate() {
            let mut plus = p;
            let mut minus = p;
            plus[k] += h * 0.5;
            minus[k] -= h * 0.5;
            let a = self.sample(plus)?;
            let b = self.sample(minus)?;
            *gp = (a - b) / h as f32;
        }
        Some(g)
    }
}

/// Estimated memory in bytes for a field of the given dimensions (f32 + occupancy bitset).
pub fn estimate_bytes(dims: [usize; 3]) -> u64 {
    let total = dims[0] as u64 * dims[1] as u64 * dims[2] as u64;
    // f32 field + f32 EDT scratch + 1 byte occupancy ≈ 9 bytes/voxel worst case.
    total.saturating_mul(9)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_constant() {
        let mut f = ScalarField::new([4, 4, 4], [0.0; 3], 1.0).unwrap();
        for v in f.data_mut() {
            *v = 2.5;
        }
        assert!((f.sample([1.5, 1.5, 1.5]).unwrap() - 2.5).abs() < 1e-6);
    }

    #[test]
    fn gradient_of_linear_field() {
        let mut f = ScalarField::new([5, 5, 5], [0.0; 3], 1.0).unwrap();
        for z in 0..5 {
            for y in 0..5 {
                for x in 0..5 {
                    *f.get_mut(x, y, z).unwrap() = x as f32; // linear in x
                }
            }
        }
        let g = f.gradient([2.0, 2.0, 2.0]).unwrap();
        assert!((g[0] - 1.0).abs() < 1e-4, "gx = {}", g[0]);
        assert!(g[1].abs() < 1e-4);
        assert!(g[2].abs() < 1e-4);
    }
}
