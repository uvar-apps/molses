//! Uniform-grid cell list for atom neighbor queries.

use crate::structure::Atom;

/// Cell size in Å is chosen so a 3×3×3 neighborhood covers any query of radius
/// `query_radius`. Build is O(n) bucket sort; queries are O(1) amortized.
#[derive(Debug, Clone)]
pub struct CellList {
    cell_size: f64,
    dims: [usize; 3],
    origin: [f64; 3],
    /// Head atom index per cell (`u32::MAX` = empty).
    heads: Vec<u32>,
    /// Next-atom indices (linked list per cell).
    next: Vec<u32>,
}

impl CellList {
    /// Build a cell list over `atoms` with the given query radius.
    pub fn new(
        atoms: &[Atom],
        query_radius: f64,
        grid_origin: [f64; 3],
        grid_dims: [usize; 3],
        spacing: f64,
    ) -> Self {
        let cell_size = (query_radius + spacing).max(spacing);
        // Cell grid spans the same domain as the scalar grid.
        let dims = [
            ((grid_dims[0] as f64 * spacing) / cell_size)
                .ceil()
                .max(1.0) as usize,
            ((grid_dims[1] as f64 * spacing) / cell_size)
                .ceil()
                .max(1.0) as usize,
            ((grid_dims[2] as f64 * spacing) / cell_size)
                .ceil()
                .max(1.0) as usize,
        ];
        let nc = dims[0] * dims[1] * dims[2];
        let mut heads = vec![u32::MAX; nc];
        let mut next = vec![u32::MAX; atoms.len()];

        for (i, atom) in atoms.iter().enumerate() {
            let (cx, cy, cz) = Self::cell_of(&grid_origin, cell_size, dims, atom.position);
            let cell = cx + cy * dims[0] + cz * dims[0] * dims[1];
            next[i] = heads[cell];
            heads[cell] = i as u32;
        }

        Self {
            cell_size,
            dims,
            origin: grid_origin,
            heads,
            next,
        }
    }

    fn cell_of(
        origin: &[f64; 3],
        cell_size: f64,
        dims: [usize; 3],
        p: [f64; 3],
    ) -> (usize, usize, usize) {
        let mut c = [0usize; 3];
        for k in 0..3 {
            let t = ((p[k] - origin[k]) / cell_size).floor();
            c[k] = if t < 0.0 {
                0
            } else {
                (t as usize).min(dims[k].saturating_sub(1))
            };
        }
        (c[0], c[1], c[2])
    }

    /// Append atom indices whose cells may contain points within `radius` of `p`.
    pub fn query(&self, p: [f64; 3], radius: f64, out: &mut Vec<u32>) {
        out.clear();
        let r_cells = (radius / self.cell_size).ceil() as i64;
        let (cx, cy, cz) = Self::cell_of(&self.origin, self.cell_size, self.dims, p);
        let cx = cx as i64;
        let cy = cy as i64;
        let cz = cz as i64;
        for dz in -r_cells..=r_cells {
            for dy in -r_cells..=r_cells {
                for dx in -r_cells..=r_cells {
                    let x = cx + dx;
                    let y = cy + dy;
                    let z = cz + dz;
                    if x < 0 || y < 0 || z < 0 {
                        continue;
                    }
                    let (x, y, z) = (x as usize, y as usize, z as usize);
                    if x >= self.dims[0] || y >= self.dims[1] || z >= self.dims[2] {
                        continue;
                    }
                    let cell = x + y * self.dims[0] + z * self.dims[0] * self.dims[1];
                    let mut i = self.heads[cell];
                    while i != u32::MAX {
                        out.push(i);
                        i = self.next[i as usize];
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::Element;

    #[test]
    fn query_finds_neighbors() {
        let atoms: Vec<Atom> = (0..10)
            .map(|i| {
                Atom::new(
                    [i as f64 * 2.0, 0.0, 0.0],
                    1.0,
                    Some(Element::C),
                    i as u32 + 1,
                )
            })
            .collect();
        let cl = CellList::new(&atoms, 3.0, [-5.0; 3], [40, 10, 10], 1.0);
        let mut out = Vec::new();
        cl.query([0.0, 0.0, 0.0], 3.0, &mut out);
        // Should at least include atom 0 and atom 1 (distance 2).
        assert!(out.contains(&0), "missing self: {out:?}");
        assert!(out.contains(&1), "missing neighbor: {out:?}");
    }
}
