//! Mesh smoothness and quality metrics.
//!
//! Provides objective measures of surface quality so that improvements to the
//! geometry engine can be quantified and regression-tested. The key insight:
//! a correct SES shape can still look rough if vertex positions or normals
//! have high-frequency noise from grid quantization.

use super::Mesh;

/// Smoothness analysis of a triangle mesh.
#[derive(Debug, Clone)]
pub struct SmoothnessReport {
    /// Mean angle (radians) between normals of adjacent vertices.
    pub mean_normal_angle: f64,
    /// Max angle (radians) between normals of adjacent vertices.
    pub max_normal_angle: f64,
    /// Mean deviation of each vertex normal from its neighbors' average (0–2).
    pub mean_normal_deviation: f64,
    /// Mean edge length (Å).
    pub mean_edge_length: f64,
    /// Roughness score: mean_normal_angle / mean_edge_length.
    /// For a smooth surface this approximates the curvature (1/R).
    /// For a noisy surface it is much larger. Lower is smoother.
    pub roughness: f64,
    /// Number of vertex pairs compared.
    pub pair_count: usize,
}

/// Compute smoothness metrics from the mesh.
///
/// Adjacent vertex pairs are those sharing a triangle edge. The angle between
/// their normals is the primary noise indicator: on a smooth surface it scales
/// with edge length × curvature; on a noisy surface it has a large
/// high-frequency component independent of edge length.
pub fn analyze(mesh: &Mesh) -> SmoothnessReport {
    if mesh.normals.len() != mesh.positions.len() || mesh.indices.is_empty() {
        return SmoothnessReport {
            mean_normal_angle: 0.0,
            max_normal_angle: 0.0,
            mean_normal_deviation: 0.0,
            mean_edge_length: 0.0,
            roughness: 0.0,
            pair_count: 0,
        };
    }

    use std::collections::HashMap;
    // Collect unique edges (vertex index pairs).
    let mut edges: HashMap<(u32, u32), ()> = HashMap::new();
    for t in &mesh.indices {
        for k in 0..3 {
            let a = t[k];
            let b = t[(k + 1) % 3];
            let key = if a < b { (a, b) } else { (b, a) };
            edges.insert(key, ());
        }
    }

    let mut angle_sum = 0.0f64;
    let mut angle_max: f64 = 0.0;
    let mut len_sum = 0.0f64;
    let mut count = 0usize;

    for &(a, b) in edges.keys() {
        let na = mesh.normals[a as usize];
        let nb = mesh.normals[b as usize];
        let dot = (na[0] * nb[0] + na[1] * nb[1] + na[2] * nb[2]) as f64;
        let angle = dot.clamp(-1.0, 1.0).acos();
        angle_sum += angle;
        angle_max = angle_max.max(angle);

        let pa = mesh.positions[a as usize];
        let pb = mesh.positions[b as usize];
        let dx = (pa[0] - pb[0]) as f64;
        let dy = (pa[1] - pb[1]) as f64;
        let dz = (pa[2] - pb[2]) as f64;
        len_sum += (dx * dx + dy * dy + dz * dz).sqrt();
        count += 1;
    }

    let mean_normal_angle = if count > 0 {
        angle_sum / count as f64
    } else {
        0.0
    };
    let mean_edge_length = if count > 0 {
        len_sum / count as f64
    } else {
        0.0
    };
    let roughness = if mean_edge_length > 1e-12 {
        mean_normal_angle / mean_edge_length
    } else {
        0.0
    };

    // Normal deviation: for each vertex, compare its normal with the average
    // of its one-ring neighbors.
    let n = mesh.positions.len();
    let mut neighbors: Vec<Vec<u32>> = vec![Vec::new(); n];
    for t in &mesh.indices {
        for k in 0..3 {
            let a = t[k] as usize;
            for j in [t[(k + 1) % 3], t[(k + 2) % 3]] {
                if !neighbors[a].contains(&j) {
                    neighbors[a].push(j);
                }
            }
        }
    }
    let mut dev_sum = 0.0f64;
    let mut dev_count = 0usize;
    for (i, nb) in neighbors.iter().enumerate() {
        if nb.is_empty() {
            continue;
        }
        let mut avg = [0.0f64; 3];
        for &j in nb {
            let nj = mesh.normals[j as usize];
            avg[0] += nj[0] as f64;
            avg[1] += nj[1] as f64;
            avg[2] += nj[2] as f64;
        }
        let inv = 1.0 / nb.len() as f64;
        avg = [avg[0] * inv, avg[1] * inv, avg[2] * inv];
        let avg_len = (avg[0] * avg[0] + avg[1] * avg[1] + avg[2] * avg[2]).sqrt();
        if avg_len < 1e-12 {
            continue;
        }
        avg = [avg[0] / avg_len, avg[1] / avg_len, avg[2] / avg_len];
        let ni = mesh.normals[i];
        let dot = (ni[0] as f64 * avg[0] + ni[1] as f64 * avg[1] + ni[2] as f64 * avg[2])
            .clamp(-1.0, 1.0);
        dev_sum += 1.0 - dot; // 0 = identical, 2 = opposite
        dev_count += 1;
    }
    let mean_normal_deviation = if dev_count > 0 {
        dev_sum / dev_count as f64
    } else {
        0.0
    };

    SmoothnessReport {
        mean_normal_angle,
        max_normal_angle: angle_max,
        mean_normal_deviation,
        mean_edge_length,
        roughness,
        pair_count: count,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A smooth sphere should have low roughness (normal angles scale with
    /// edge length / radius). This is the REFERENCE case for smoothness.
    #[test]
    fn sphere_is_smooth() {
        // Build a sphere mesh via polygonize.
        let r = 1.0f64;
        let h = 0.2;
        let extent = r + 0.5;
        let n = (2.0 * extent / h).ceil() as usize + 1;
        let mut field = crate::field::ScalarField::new([n, n, n], [-extent; 3], h).unwrap();
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let p = field.position(x, y, z);
                    let d = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                    *field.get_mut(x, y, z).unwrap() = (d - r) as f32;
                }
            }
        }
        let mut mesh = crate::mesh::polygonize::polygonize(&field, None, 1.4);
        crate::mesh::normals::compute_normals(&mut mesh, Some(&field));

        let report = analyze(&mesh);
        // For a smooth sphere of radius 1.0 with edge length ~h:
        // expected normal angle ≈ h/r ≈ 0.2 rad.
        // Roughness = angle/length ≈ 1/r = 1.0.
        assert!(
            report.roughness < 2.0,
            "Sphere roughness {:.3} should be < 2.0 (≈1/R for R=1). Report: {:?}",
            report.roughness,
            report
        );
        assert!(
            report.mean_normal_deviation < 0.05,
            "Normal deviation {:.4} should be small for a smooth sphere",
            report.mean_normal_deviation
        );
    }

    /// A sphere with random vertex displacement (noise) should have high
    /// roughness. This is the DETECTION case — the metric must catch noise.
    #[test]
    fn detects_noisy_mesh() {
        // Build a sphere mesh.
        let r = 1.0f64;
        let h = 0.2;
        let extent = r + 0.5;
        let n = (2.0 * extent / h).ceil() as usize + 1;
        let mut field = crate::field::ScalarField::new([n, n, n], [-extent; 3], h).unwrap();
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let p = field.position(x, y, z);
                    let d = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                    *field.get_mut(x, y, z).unwrap() = (d - r) as f32;
                }
            }
        }
        let mut mesh = crate::mesh::polygonize::polygonize(&field, None, 1.4);
        crate::mesh::normals::compute_normals(&mut mesh, Some(&field));

        let smooth_report = analyze(&mesh);

        // Add noise to normals (simulating quantization artifacts).
        // Deterministic pseudo-random perturbation with realistic amplitude:
        // EDT quantization can tilt normals by ~0.1–0.3 rad at grid boundaries.
        for (i, n) in mesh.normals.iter_mut().enumerate() {
            let phase = (i as u32).wrapping_mul(2654435761);
            let noise_a = (phase as f32 / u32::MAX as f32 - 0.5) * 0.6;
            let noise_b = ((phase >> 8) as f32 / u32::MAX as f32 - 0.5) * 0.6;
            n[0] += noise_a;
            n[1] += noise_b;
            n[2] += noise_a * 0.3;
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            if len > 1e-12 {
                n[0] /= len;
                n[1] /= len;
                n[2] /= len;
            }
        }

        let noisy_report = analyze(&mesh);

        assert!(
            noisy_report.mean_normal_deviation > smooth_report.mean_normal_deviation * 1.5,
            "Noisy mesh deviation {:.4} should be > 1.5× smooth {:.4}",
            noisy_report.mean_normal_deviation,
            smooth_report.mean_normal_deviation
        );
        assert!(
            noisy_report.roughness > smooth_report.roughness * 1.2,
            "Noisy roughness {:.3} should exceed smooth {:.3} by >20%",
            noisy_report.roughness,
            smooth_report.roughness
        );
    }
}
