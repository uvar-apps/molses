# Lead reconnaissance notes (primary agent)

Working notes from direct searches and PDF extraction. Merged into final docs alongside
`01-classical-theory.md`, `02-modern-implementations.md`, `03-rust-ecosystem.md`.

## Core definitions (verified against primary sources)

### SAS (Solvent Accessible Surface)
- Traced by the **center** of a spherical probe of radius `r_p` rolling over the van der Waals (VdW) spheres.
- Equals the topological boundary of the union of inflated spheres `B(c_i, r_i + r_p)`.
- Source: Lee & Richards 1971, *J. Mol. Biol.* 55, 379–400 (SAS); Richards 1977 *Annu. Rev. Biophys. Bioeng.* 6, 151–176.

### SES (Solvent Excluded Surface / Connolly surface / molecular surface)
- Surface traced by the **outer surface** of the probe as it rolls over the molecule.
- Sanner et al. 1996 definition: *"the topological boundary of the union of all possible probes having no intersection with M"* — equivalently the boundary of the volume the probe body can never occupy (the SES cavity).
- Quan & Stamm 2015 (HAL hal-01205580 / JCP 2016): *"the SES is the boundary of the union of all spherical probes that do not intersect the VdW-balls"*.
- Composed of three patch types (Connolly 1983):
  1. **Convex / contact spherical patches** — on VdW spheres (probe touching one atom).
  2. **Toric / toroidal reentrant patches** — swept by probe rolling over pairs of atoms.
  3. **Concave spherical (spheric reentrant) patches** — on probe surface, probe touching three atoms.
- Greer & Bush 1978, PNAS 75(1), 303–307 coined the "solvent exclusion" surface.

## Key mathematical formulation (for a robust grid engine)

Let `K = ∪_i B(c_i, r_i + r_p)` be the SAS cavity (Minkowski sum of VdW union with probe ball).
Let `f_sas` be the **signed distance function** to `∂K` (negative inside `K`).

**Quan & Stamm, eq. (3.4) / (5.21):**

```
f_ses(p) = f_sas(p) + r_p
```

Zero level set of `f_ses` = SES. Points on SES have signed distance `-r_p` to the SAS
(they sit `r_p` inside the SAS cavity). This is the cleanest basis for a
marching-cubes / surface-nets implementation.

Consequences:
- Outside `K`: `f_sas(p) = min_i (||p-c_i|| - (r_i+r_p))`, so
  `f_ses(p) = min_i (||p-c_i|| - r_i)` = signed distance to the VdW union (exact).
- Inside `K`: true distance-to-boundary is required (not the naive min over inflated
  spheres, which underestimates depth in overlap regions). Options: Euclidean distance
  transform on a grid, or exact union-of-balls closest-point via the Voronoi / medial
  decomposition (Quan & Stamm §4.2 three equivalence statements).
- SES cavity = `{p : f_ses(p) ≤ 0}` = VdW union **plus** crevices/pockets the probe body
  cannot enter (depth from `∂K` greater than `r_p`).

### Singularities
- SES is smooth a.e. but has non-differentiable sets where multiple closest SAS points exist
  (Corollary 5.1, Quan & Stamm): cusps on toroidal patches, singular circles on concave patches.
- MSMS handles radial + non-radial singularities with clipping heuristics; can fail and
  restart by growing problem atom radii by 0.1 Å (see MSMS docs `-no_rest_on_pbr`).
- Grid isosurface methods sidestep analytic singularities: the SDF is still well defined
  (Lipschitz, |∇f|=1 a.e.), the mesh just rounds features below the grid scale.

### Analytic triangulation route (MSMS / Connolly)
- **Reduced surface** `R_r(S)`: faces = fixed probe positions contacting ≥3 atoms;
  free edges/vertices where probe contacts 2/1 atoms. Size `O(n)` for calibrated sphere sets.
- Complexity `O(n log n)` to build reduced surface, `O(n)` for analytic SES description.
- Triangulation density in **vertices per Å²** (MSMS default 1.0; typical 1.0 for >1000 atoms, 3.0 for small).
- MSMS default probe radius **1.5 Å** in the CLI docs (Sanner 1996 used 1.5); water probe commonly **1.4 Å**.

## Performance reference points (Sanner 1996 / MSMS)
- Crambin (1crn, 327 atoms) → 4772 triangles in 0.73 s (HP 9000/735).
- Thermolysin (3tln, 2437 atoms) → 26462 triangles in 4.6 s.
- Glutamine synthetase (2gls, 43632 atoms) → 476665 triangles in 104.53 s.
- Modern CPUs should beat these by orders of magnitude; NanoShaper/JCP-2018 report multicore scalability.

## Grid / implicit route (what we will build first)
- Quan & Stamm: analytic SDF + implicit SES (foundational).
- Cheng et al.? / JCP 2018 "Fast and scalable algorithms for constructing Solvent-Excluded
  Surfaces of Large Biomolecules" — level-set representation on **adaptive Cartesian octree**,
  distributed computing, load-balanced, **first-order accurate**, feature-preserving.
  DOI: 10.1016/j.jcp.2018.07.035. (Abstract via ACM DL / ScienceDirect.)
- EDTSurf (Xu & Zhang 2009, PLOS ONE 4(12):e8140): fast Euclidean Distance Transform to
  classify/build vdW, SAS, and molecular (SES) surfaces on a grid.
- NanoShaper (Decherchi & Rocchia): grid + ray-casting, multicore, robust; probe radius
  and "skin" (smoothness) parameters; also pocket detection. Springer tutorial DOI 10.1007/978-3-319-12211-3_10.
- Marching Cubes: Lorensen & Cline, SIGGRAPH 1987, DOI 10.1145/37402.37422 (case table, linear interpolation).
- Alternatives: Surface Nets / Dual Contouring (manifold, fewer sharper triangles);
  see arXiv:2401.14906 for a high-performance parallel SurfaceNets.

## LOD / resolution controls seen in the wild
| Tool | Parameter | Typical values |
|------|-----------|----------------|
| MSMS | `density` vertices/Å² | 1.0 (large), 3.0 (small) |
| MSMS | `probe_radius` | 1.5 default, 1.4 common for water |
| NanoShaper | probe radius + skin (smoothness) | 1.4–1.6, skin ~0.5 |
| Grid methods | grid spacing h (Å) | 0.2–0.5; error O(h)–O(h²) depending on interpolation |
| JCP 2018 | octree levels | adaptive |

Triangle count for a grid isosurface scales as `~ Area / h²`. A 30 000 Å² protein at h=0.4 Å
is on the order of 10^5–10^6 vertices — comfortable for a Rust rayon+SIMD pipeline.

## Input formats
- **PDB**: fixed-column; ATOM/HETATM; altLoc, insertion codes, occupancy, element column;
  multi-MODEL ensembles; anisou records ignored by most surface tools.
- **mmCIF / PDBx**: category-based, supersedes PDB; large assemblies; struct_conn etc.
- **XYZ**: minimal element + coordinates; no bonds, no residue info; radii must come from element.
- VdW radii: **Bondi 1964**, *J. Phys. Chem.* 68(3), 441–451, DOI 10.1021/j100785a001
  (38 elements; standard fallback table for molecular surfaces).
- Units: Ångström throughout PDB ecosystem; OBJ is unitless — we will document that we write Å.

## OBJ export conventions (for molecular surfaces)
- `v x y z` vertices in Å, `vn nx ny nz` unit normals, `f v//vn v//vn v//vn` triangles.
- Per-vertex normals preferred for smooth shading of SES (area/angle-weighted, or analytic SDF gradient).
- Winding: counter-clockwise when viewed from outside (right-handed, outward normals).
- Reference: Wavefront OBJ; Library of Congress format description; Paul Bourke notes.

## Notable citations collected
1. Lee, B. & Richards, F.M. (1971). *The interpretation of protein structures: Estimation of static accessibility.* J. Mol. Biol. 55, 379–400.
2. Richards, F.M. (1977). *Areas, volumes, packing, and protein structure.* Annu. Rev. Biophys. Bioeng. 6, 151–176.
3. Greer, J. & Bush, B.L. (1978). *Macromolecular shape and surface maps by solvent exclusion.* PNAS 75(1), 303–307. DOI 10.1073/pnas.75.1.303.
4. Connolly, M.L. (1983). *Analytical molecular surface calculation.* J. Appl. Cryst. 16, 548–558. DOI 10.1107/S0021889883010985.
5. Connolly, M.L. (1983). *Solvent-accessible surfaces of proteins and nucleic acids.* Science 221, 709–713.
6. Bondi, A. (1964). *van der Waals volumes and radii.* J. Phys. Chem. 68(3), 441–451. DOI 10.1021/j100785a001.
7. Sanner, M.F., Olson, A.J. & Spehner, J.-C. (1996). *Reduced surface: An efficient way to compute molecular surfaces.* Biopolymers 38, 305–320. (MSMS)
8. Lorensen, W.E. & Cline, H.E. (1987). *Marching cubes: A high resolution 3D surface construction algorithm.* SIGGRAPH '87. DOI 10.1145/37402.37422.
9. Quan, C. & Stamm, B. (2016). *Mathematical analysis and calculation of molecular surfaces.* J. Comput. Phys. 322, 221–242 (preprint HAL hal-01205580). DOI 10.1016/j.jcp.2016.07.007.
10. Cheng, Y. et al.? (2018). *Fast and scalable algorithms for constructing Solvent-Excluded Surfaces of large biomolecules.* J. Comput. Phys. DOI 10.1016/j.jcp.2018.07.035.
11. Xu, D. & Zhang, Y. (2009). *Generating triangulated macromolecular surfaces by Euclidean distance transform.* PLOS ONE 4(12): e8140. (EDTSurf)
12. Decherchi, S. & Rocchia, W. *NanoShaper* tutorial (2015) DOI 10.1007/978-3-319-12211-3_10; NanoShaper–VMD interface, Bioinformatics 35(7), 1241–1243.

## Open questions to resolve in the spec
- Default probe radius: 1.4 Å (water, physicochemical convention) vs 1.5 Å (MSMS default). **Recommend 1.4 Å with `--probe-radius` override.**
- Primary engine: grid SDF + isosurface (robust, LOD-friendly, SIMD/rayon-friendly) vs analytic reduced-surface (exact patches, hard singularities). **Recommend grid-first, analytic as later optional exact mode.**
- LOD parameterization: expose `grid spacing` (Å) as primary, with named presets (`draft` 0.6 / `balanced` 0.35 / `fine` 0.2 / `ultra` 0.15), and optional triangle budget.
- Whether to emit multiple OBJ components (one per connected SES component) or a single mesh; MSMS writes per-component files.
