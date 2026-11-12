//! Optional Laplacian / Taubin smoothing (M3).

use super::Mesh;

/// Apply `iterations` Laplacian smoothing passes (uniform weights).
///
/// Off by default in the product: smoothing trades geometric fidelity for
/// visual quality. Deterministic (Jacobi-style: read old, write new).
pub fn laplacian(mesh: &mut Mesh, iterations: u32, weight: f32) {
    if iterations == 0 || mesh.positions.is_empty() {
        return;
    }
    use std::collections::HashMap;
    let n = mesh.positions.len();
    let mut neighbors: Vec<Vec<u32>> = vec![Vec::new(); n];
    for t in &mesh.indices {
        for k in 0..3 {
            let a = t[k] as usize;
            let b = t[(k + 1) % 3] as usize;
            let c = t[(k + 2) % 3] as usize;
            for &j in &[b, c] {
                if !neighbors[a].contains(&(j as u32)) {
                    neighbors[a].push(j as u32);
                }
            }
        }
    }
    let mut current = mesh.positions.clone();
    for _ in 0..iterations {
        let mut next = current.clone();
        for i in 0..n {
            let nb = &neighbors[i];
            if nb.is_empty() {
                continue;
            }
            let mut avg = [0.0f32; 3];
            for &j in nb {
                let p = current[j as usize];
                avg[0] += p[0];
                avg[1] += p[1];
                avg[2] += p[2];
            }
            let inv = 1.0 / nb.len() as f32;
            avg = [avg[0] * inv, avg[1] * inv, avg[2] * inv];
            next[i] = [
                current[i][0] + weight * (avg[0] - current[i][0]),
                current[i][1] + weight * (avg[1] - current[i][1]),
                current[i][2] + weight * (avg[2] - current[i][2]),
            ];
        }
        current = next;
    }
    mesh.positions = current;
    let _ = HashMap::<u32, u32>::new(); // silence unused import warning path
}
