//! Triangle mesh representation and processing (M3).

pub mod components;
pub mod normals;
pub mod polygonize;
pub mod quality;
pub mod smooth;
pub mod topology;

/// A triangle mesh of the SES.
#[derive(Debug, Clone, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Mesh {
    /// Vertex positions in Å.
    pub positions: Vec<[f32; 3]>,
    /// Unit normals, parallel to `positions`.
    pub normals: Vec<[f32; 3]>,
    /// Triangle indices (CCW outward), 0-based.
    pub indices: Vec<[u32; 3]>,
    /// Nearest atom index per vertex (optional metadata).
    pub nearest_atom: Vec<u32>,
}

impl Mesh {
    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.positions.len()
    }

    /// Number of triangles.
    pub fn triangle_count(&self) -> usize {
        self.indices.len()
    }

    /// `true` when the mesh has no triangles.
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }

    /// Total surface area in Å² (sum of triangle areas).
    pub fn area(&self) -> f64 {
        let mut acc = 0.0f64;
        for t in &self.indices {
            let a = self.positions[t[0] as usize];
            let b = self.positions[t[1] as usize];
            let c = self.positions[t[2] as usize];
            acc += triangle_area(a, b, c);
        }
        acc
    }

    /// Signed enclosed volume in Å³ (divergence theorem; positive if outward).
    pub fn volume(&self) -> f64 {
        let mut acc = 0.0f64;
        for t in &self.indices {
            let a = self.positions[t[0] as usize];
            let b = self.positions[t[1] as usize];
            let c = self.positions[t[2] as usize];
            // Scalar triple product / 6.
            let triple = a[0] as f64 * (b[1] as f64 * c[2] as f64 - b[2] as f64 * c[1] as f64)
                - a[1] as f64 * (b[0] as f64 * c[2] as f64 - b[2] as f64 * c[0] as f64)
                + a[2] as f64 * (b[0] as f64 * c[1] as f64 - b[1] as f64 * c[0] as f64);
            acc += triple / 6.0;
        }
        acc
    }
}

/// Area of a triangle given its vertices.
pub fn triangle_area(a: [f32; 3], b: [f32; 3], c: [f32; 3]) -> f64 {
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let cross = [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ];
    let norm = ((cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]) as f64).sqrt();
    0.5 * norm
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_triangle() {
        let m = Mesh {
            positions: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            normals: vec![[0.0, 0.0, 1.0]; 3],
            indices: vec![[0, 1, 2]],
            nearest_atom: vec![],
        };
        assert!((m.area() - 0.5).abs() < 1e-6);
    }

    #[test]
    fn tetrahedron_volume() {
        // Regular-ish tet with vertices at origin and unit axes: volume = 1/6.
        let m = Mesh {
            positions: vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            normals: vec![[0.0; 3]; 4],
            indices: vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
            nearest_atom: vec![],
        };
        assert!((m.volume().abs() - 1.0 / 6.0).abs() < 1e-6);
    }
}
