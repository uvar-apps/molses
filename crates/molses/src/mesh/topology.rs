//! Mesh topology checks (M3).

use super::Mesh;

/// Result of a watertightness / manifoldness audit.
#[derive(Debug, Clone, Default)]
pub struct TopologyReport {
    /// Number of vertices.
    pub vertices: usize,
    /// Number of triangles.
    pub triangles: usize,
    /// Edges used by exactly two triangles.
    pub manifold_edges: usize,
    /// Edges used by anything other than two triangles.
    pub non_manifold_edges: usize,
    /// Boundary edges (used by exactly one triangle).
    pub boundary_edges: usize,
    /// Euler characteristic `V − E + F`.
    pub euler: i64,
}

impl TopologyReport {
    /// `true` when every edge is used by exactly two triangles.
    pub fn is_watertight(&self) -> bool {
        self.non_manifold_edges == 0 && self.boundary_edges == 0 && self.triangles > 0
    }
}

/// Audit mesh topology.
pub fn audit(mesh: &Mesh) -> TopologyReport {
    use std::collections::HashMap;
    let mut edge_use: HashMap<(u32, u32), u32> = HashMap::new();
    for t in &mesh.indices {
        for k in 0..3 {
            let a = t[k];
            let b = t[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_use.entry(key).or_insert(0) += 1;
        }
    }
    let mut manifold = 0usize;
    let mut non_manifold = 0usize;
    let mut boundary = 0usize;
    for &c in edge_use.values() {
        match c {
            1 => boundary += 1,
            2 => manifold += 1,
            _ => non_manifold += 1,
        }
    }
    let v = mesh.vertex_count() as i64;
    let e = edge_use.len() as i64;
    let f = mesh.triangle_count() as i64;
    TopologyReport {
        vertices: mesh.vertex_count(),
        triangles: mesh.triangle_count(),
        manifold_edges: manifold,
        non_manifold_edges: non_manifold,
        boundary_edges: boundary,
        euler: v - e + f,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tetrahedron_is_watertight() {
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
        let r = audit(&m);
        assert!(r.is_watertight(), "{r:?}");
        assert_eq!(r.euler, 2);
    }

    #[test]
    fn open_triangle_is_not() {
        let m = Mesh {
            positions: vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            normals: vec![[0.0; 3]; 3],
            indices: vec![[0, 1, 2]],
            nearest_atom: vec![],
        };
        let r = audit(&m);
        assert!(!r.is_watertight());
        assert_eq!(r.boundary_edges, 3);
    }
}
