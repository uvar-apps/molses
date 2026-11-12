//! Connected-component analysis of triangle meshes.

use super::Mesh;

/// Information about one connected component of a mesh.
#[derive(Debug, Clone)]
pub struct Component {
    /// Triangle indices belonging to this component (indices into `mesh.indices`).
    pub triangles: Vec<u32>,
    /// Vertex indices used by this component.
    pub vertices: Vec<u32>,
    /// Surface area in Å².
    pub area: f64,
    /// Signed enclosed volume in Å³ (positive if normals point outward).
    pub volume: f64,
}

/// Label connected components (via shared welded vertices) and compute their
/// area and volume. Components are returned in descending area order; index 0
/// is conventionally the outer shell (largest enclosed volume).
pub fn connected_components(mesh: &Mesh) -> Vec<Component> {
    let n_v = mesh.positions.len();
    if n_v == 0 || mesh.indices.is_empty() {
        return Vec::new();
    }

    // Union-find over vertices, unioned by triangle edges.
    let mut parent: Vec<u32> = (0..n_v as u32).collect();
    fn find(parent: &mut [u32], x: u32) -> u32 {
        let mut r = x;
        while parent[r as usize] != r {
            r = parent[r as usize];
        }
        // Path compression.
        let mut c = x;
        while parent[c as usize] != r {
            let next = parent[c as usize];
            parent[c as usize] = r;
            c = next;
        }
        r
    }
    for t in &mesh.indices {
        for k in 0..3 {
            let a = find(&mut parent, t[k]);
            let b = find(&mut parent, t[(k + 1) % 3]);
            if a != b {
                parent[a as usize] = b;
            }
        }
    }

    // Group triangles by component root.
    use std::collections::BTreeMap;
    let mut comp_tris: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for (ti, t) in mesh.indices.iter().enumerate() {
        let root = find(&mut parent, t[0]);
        comp_tris.entry(root).or_default().push(ti as u32);
    }

    let mut comps: Vec<Component> = comp_tris
        .into_values()
        .map(|tris| {
            let mut vset = std::collections::BTreeSet::new();
            let mut area = 0.0f64;
            let mut vol = 0.0f64;
            for &ti in &tris {
                let t = mesh.indices[ti as usize];
                for &vi in &t {
                    vset.insert(vi);
                }
                let a = mesh.positions[t[0] as usize];
                let b = mesh.positions[t[1] as usize];
                let c = mesh.positions[t[2] as usize];
                area += super::triangle_area(a, b, c);
                let triple = a[0] as f64 * (b[1] as f64 * c[2] as f64 - b[2] as f64 * c[1] as f64)
                    - a[1] as f64 * (b[0] as f64 * c[2] as f64 - b[2] as f64 * c[0] as f64)
                    + a[2] as f64 * (b[0] as f64 * c[1] as f64 - b[1] as f64 * c[0] as f64);
                vol += triple / 6.0;
            }
            Component {
                triangles: tris,
                vertices: vset.into_iter().collect(),
                area,
                volume: vol,
            }
        })
        .collect();

    // Sort by |volume| descending so index 0 is the outer shell.
    comps.sort_by(|a, b| {
        b.volume
            .abs()
            .partial_cmp(&a.volume.abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    comps
}

/// Keep only the triangles of the selected components; rebuild a compact mesh
/// with re-indexed vertices and normals.
pub fn filter_components(mesh: &Mesh, keep: &[Component]) -> Mesh {
    use std::collections::HashMap;
    let mut old_to_new: HashMap<u32, u32> = HashMap::new();
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut nearest_atom = Vec::new();
    let mut indices = Vec::new();
    let has_normals = mesh.normals.len() == mesh.positions.len();
    let has_na = mesh.nearest_atom.len() == mesh.positions.len();

    for comp in keep {
        for &vi in &comp.vertices {
            if let std::collections::hash_map::Entry::Vacant(e) = old_to_new.entry(vi) {
                let ni = positions.len() as u32;
                e.insert(ni);
                positions.push(mesh.positions[vi as usize]);
                if has_normals {
                    normals.push(mesh.normals[vi as usize]);
                }
                if has_na {
                    nearest_atom.push(mesh.nearest_atom[vi as usize]);
                }
            }
        }
        for &ti in &comp.triangles {
            let t = mesh.indices[ti as usize];
            indices.push([old_to_new[&t[0]], old_to_new[&t[1]], old_to_new[&t[2]]]);
        }
    }

    Mesh {
        positions,
        normals,
        indices,
        nearest_atom,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_spheres_mesh() -> Mesh {
        // Two disjoint tetrahedra = two components.
        let mut positions = vec![
            // tet A
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            // tet B (offset)
            [5.0, 0.0, 0.0],
            [6.0, 0.0, 0.0],
            [5.0, 1.0, 0.0],
            [5.0, 0.0, 1.0],
        ];
        let mut normals = vec![[0.0; 3]; 8];
        let indices = vec![
            // tet A (outward winding)
            [0, 2, 1],
            [0, 1, 3],
            [0, 3, 2],
            [1, 2, 3],
            // tet B
            [4, 6, 5],
            [4, 5, 7],
            [4, 7, 6],
            [5, 6, 7],
        ];
        Mesh {
            positions: std::mem::take(&mut positions),
            normals: std::mem::take(&mut normals),
            indices,
            nearest_atom: vec![],
        }
    }

    #[test]
    fn two_components_found() {
        let mesh = two_spheres_mesh();
        let comps = connected_components(&mesh);
        assert_eq!(comps.len(), 2);
        assert_eq!(comps[0].triangles.len(), 4);
        assert_eq!(comps[1].triangles.len(), 4);
    }

    #[test]
    fn filter_keeps_one() {
        let mesh = two_spheres_mesh();
        let comps = connected_components(&mesh);
        let filtered = filter_components(&mesh, &comps[..1]);
        assert_eq!(filtered.triangle_count(), 4);
        assert_eq!(filtered.vertex_count(), 4);
    }

    #[test]
    fn component_volume_positive() {
        let mesh = two_spheres_mesh();
        let comps = connected_components(&mesh);
        // Each tet has volume 1/6 with outward winding.
        assert!((comps[0].volume.abs() - 1.0 / 6.0).abs() < 1e-6);
    }
}
