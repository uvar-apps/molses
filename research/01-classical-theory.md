# Classical Theory and Algorithms for Solvent Excluded Surface (SES)

Research findings for the SES library+CLI project. Every factual claim carries a source.
Uncertainty is marked explicitly where applicable.

---

## 1. Precise Mathematical Definition of SES

### 1.1 The three surfaces: vdW, SAS, and SES

**Van der Waals (vdW) surface.** Given a molecule with M atoms, each atom i has center c_i and van der Waals radius r_i. The vdW surface is the topological boundary of the union of all atomic spheres:

  Ω_vdW = ⋃_{i=1}^{M} B(c_i, r_i)

where B(c, r) = {p ∈ ℝ³ : ‖p − c‖ ≤ r}. The boundary ∂Ω_vdW is the vdW surface. It exhibits cusps where spheres intersect. [Quan & Stamm 2016, J Comput Phys; https://hal.sorbonne-universite.fr/hal-01205580/document]

**Solvent-Accessible Surface (SAS).** Introduced by Lee & Richards (1971). A spherical probe of radius r_p (representing a solvent molecule, typically water at r_p ≈ 1.4 Å) rolls over the vdW surface. The SAS is the locus of the *center* of the probe as it rolls over the vdW surface. Equivalently, the SAS is the boundary of the union of all *inflated* atomic spheres with radii (r_i + r_p):

  Ω_SAS = ⋃_{i=1}^{M} B(c_i, r_i + r_p)

The SAS = ∂Ω_SAS. [Lee & Richards 1971, J Mol Biol 55:379–400; https://doi.org/10.1016/0022-2836(71)90324-X]
[Sanner, Olson & Spehner 1996, Biopolymers 38:305–320; https://biomolmd.org/mw/images/4/49/Sanner96.pdf]
[Quan & Stamm 2016; https://hal.sorbonne-universite.fr/hal-01205580/document]

**Solvent Excluded Surface (SES), also called Connolly surface or molecular surface.** Introduced by Richards (1977). The SES is the boundary of the region from which the *center* of the probe sphere is excluded — equivalently, the boundary of the volume that the probe sphere *cannot* penetrate. More precisely:

  Ω_SES = Ω_vdW ∪ (regions swept by the probe sphere body as it rolls over Ω_vdW)

Or equivalently: Ω_SES is the closure of the complement of the set of all points reachable by the probe sphere center that are *outside* the inflated union. The SES is the boundary of the "cavity" occupied by the solute plus the volume inaccessible to solvent molecules.

A clean formalization (Quan & Stamm 2016): The SES is the set of points p such that the signed distance from p to the SAS equals exactly −r_p. That is:

  SES = {p ∈ ℝ³ : f_sas(p) = −r_p}

where f_sas is the signed distance function to the SAS surface. [Quan & Stamm 2016, Eq. 3.4; https://hal.sorbonne-universite.fr/hal-01205580/document]

Richards (1977) wrote: "Considered together the contact and reentrant surfaces represent a continuous sheet, which might be called the molecular surface." [Richards 1977, Annu Rev Biophys Bioeng 6:151–176; https://biomolmd.org/mw/images/e/e5/Richards77.pdf]

### 1.2 The three parts of the SES

The SES consists of two fundamental types of patches (Richards 1977; Connolly 1983):

1. **Contact surface** (convex spherical patches): Parts of the vdW surface that can actually be touched by the probe sphere. These are spherical patches centered on atoms, with radius r_i (the atomic vdW radius). Each such patch is the portion of the atom's vdW sphere where the probe can make contact.

2. **Reentrant surface** (concave patches): The inward-facing part of the probe sphere when it is simultaneously in contact with more than one atom. The reentrant surface is subdivided into:
   - **Toric reentrant surfaces** (saddle-shaped toroidal patches): The inward-facing part of the probe when it touches exactly *two* atoms simultaneously. These are patches of a torus.
   - **Spheric reentrant surfaces** (concave spherical triangles): The inward-facing part of the probe when it touches *three or more* atoms simultaneously. These are spherical triangles on the probe sphere.

[Richards 1977; https://biomolmd.org/mw/images/e/e5/Richards77.pdf]
[Connolly 1983, J Appl Crystallogr 16:548–558; https://doi.org/10.1107/S0021889883010985]
[Sanner, Olson & Spehner 1996; https://biomolmd.org/mw/images/4/49/Sanner96.pdf]
[Quan & Stamm 2016; https://hal.sorbonne-universite.fr/hal-01205580/document]

Connolly (1983) formally divided the SES into three types of patches: convex spherical patches, saddle-shaped toroidal patches, and concave spherical triangles. The self-intersection among different patches in this division often causes singularities despite that the whole SES is smooth almost everywhere. [Quan & Stamm 2016, §1.1; https://hal.sorbonne-universite.fr/hal-01205580/document]

### 1.3 Formalization via probe exclusion

The SES boundary can be characterized as follows: a point p is on the SES if and only if:
- The probe sphere of radius r_p centered at p is tangent to at least one atomic sphere (contact condition), OR
- The probe sphere centered at p is tangent to two or more atomic spheres simultaneously (reentrant condition), OR
- p lies on the vdW surface and the probe cannot reach p from outside (contact surface).

More precisely, the SES = ∂(Ω_vdW ∪ ⋃_{x ∈ Ω_SAS} B(x, r_p)) — the boundary of the vdW union inflated by the probe ball centered at every SAS point. But since the SAS is itself the boundary of the inflated union, this simplifies.

**Key identity (Quan & Stamm 2016, Eq. 3.4):**

  f_ses(p) = f_sas(p) + r_p

where f_sas(p) is the signed distance function to the SAS. The SES is the zero level set of f_ses. This is the foundation for implicit/grid methods. [https://hal.sorbonne-universite.fr/hal-01205580/document]

### 1.4 Terminology ambiguity

There is significant terminology confusion in the literature:
- "Molecular surface" (Richards 1977) = SES = Connolly surface
- "Solvent-accessible surface" (Lee & Richards 1971) = SAS = Lee-Richards surface
- Some authors use "SAS" to mean the SES (Lange et al. 2019 note this ambiguity explicitly)
- The term "Connolly surface" acknowledges Connolly's role in popularizing and first analytically computing the SES

[Lange, Herbert, Albrecht & You 2019/2020, Mol Phys 118:e1644384; https://par.nsf.gov/servlets/purl/10276515]
[Wikipedia "Accessible surface area"; https://en.wikipedia.org/wiki/Accessible_surface_area]

---

## 2. Original Algorithms with Full Citations

### 2.1 Lee & Richards (1971) — Accessible surface area

- **Authors:** Lee, B.; Richards, F.M.
- **Year:** 1971
- **Title:** "The interpretation of protein structures: estimation of static accessibility"
- **Venue:** Journal of Molecular Biology, 55(3):379–400
- **DOI:** 10.1016/0022-2836(71)90324-X
- **URL:** https://doi.org/10.1016/0022-2836(71)90324-X
- **Contribution:** First definition and computation of the accessible surface area (SAS) using a rolling probe sphere. The SAS is the locus of probe centers rolling over the vdW surface.

### 2.2 Richards (1977) — Molecular surface (SES) definition

- **Authors:** Richards, Frederic M.
- **Year:** 1977
- **Title:** "Areas, volumes, packing and protein structure"
- **Venue:** Annual Review of Biophysics and Bioengineering, 6:151–176
- **DOI:** 10.1146/annurev.bb.06.060177.001055
- **URL:** https://doi.org/10.1146/annurev.bb.06.060177.001055
- **Free PDF:** https://biomolmd.org/mw/images/e/e5/Richards77.pdf
- **Contribution:** Defined the "molecular surface" (= SES) as the combination of contact surface and reentrant surface. Introduced the three-part decomposition (contact, reentrant, and the concept of toric vs spheric reentrant). Noted that the molecular surface area approaches a finite limiting value as probe radius increases, and that the smallest physically reasonable probe is a water molecule of radius 1.4–1.5 Å.

### 2.3 Shrake & Rupley (1973) — Numerical SAS computation

- **Authors:** Shrake, A.; Rupley, J.A.
- **Year:** 1973
- **Title:** "Environment and exposure to solvent of protein atoms. Lysozyme and insulin"
- **Venue:** Journal of Molecular Biology, 79(2):351–371
- **DOI:** 10.1016/0022-2836(73)90011-9
- **URL:** https://doi.org/10.1016/0022-2836(73)90011-9
- **Contribution:** Numerical "rolling ball" algorithm for SAS area computation using point sampling on the vdW surface.

### 2.4 Connolly (1983) — Analytical SES calculation

Three papers by Michael L. Connolly in 1983 are relevant:

**(a) J Appl Crystallogr:**
- **Authors:** Connolly, Michael L.
- **Year:** 1983
- **Title:** "Analytical molecular surface calculation"
- **Venue:** Journal of Applied Crystallography, 16(5):548–558
- **DOI:** 10.1107/S0021889883010985
- **URL:** https://doi.org/10.1107/S0021889883010985
- **Contribution:** First analytical algorithm for computing the SES. The surface is composed of pieces of spheres and tori joined at circular arcs. Spheres, tori, and arcs are defined by analytical expressions in terms of atomic coordinates, vdW radii, and probe radius. The area of each surface piece can be computed analytically.

**(b) Science:**
- **Authors:** Connolly, Michael L.
- **Year:** 1983
- **Title:** "Solvent-Accessible Surfaces of Proteins and Nucleic Acids"
- **Venue:** Science, 221(4612):709–713
- **DOI:** 10.1126/science.6879170
- **URL:** https://doi.org/10.1126/science.6879170
- **Free PDF:** http://www.gersteinlab.org/courses/452/09-spring/pdf/Connolly.pdf
- **Contribution:** Presented both the dot-surface (numerical) method and the new analytical surface method for SES. Popularized the SES visualization for protein chemistry applications.

**(c) Molecular Surface Package (1993):**
- **Authors:** Connolly, Michael L.
- **Year:** 1993
- **Title:** "The molecular surface package"
- **Venue:** Journal of Molecular Graphics, 11(2):139–141
- **DOI:** 10.1016/0263-7855(93)87010-3
- **URL:** https://doi.org/10.1016/0263-7855(93)87010-3
- **Contribution:** The MSMS/PQMS software package (Connolly's implementation). Referenced as "PQMS" in Sanner et al. 1996 for comparison.

### 2.5 Richmond (1984) — Analytical SAS area and excluded volume

- **Authors:** Richmond, Timothy J.
- **Year:** 1984
- **Title:** "Solvent accessible surface area and excluded volume in proteins. Analytical equations for overlapping spheres and implications for the hydrophobic effect"
- **Venue:** Journal of Molecular Biology, 178(1):63–89
- **DOI:** 10.1016/0022-2836(84)90231-6
- **URL:** https://doi.org/10.1016/0022-2836(84)90231-6
- **Contribution:** Derived analytical formulas for SAS area of overlapping spheres and directional derivatives for minimization/docking. Implemented 3D rolling probe. Referenced alongside Connolly for the 3D implementation of the rolling-ball algorithm. [Wikipedia "Accessible surface area"; https://en.wikipedia.org/wiki/Accessible_surface_area]

### 2.6 Edelsbrunner — Weighted alpha shapes and union of balls

**(a) Alpha shapes (original):**
- **Authors:** Edelsbrunner, H.; Mücke, E.P.
- **Year:** 1994
- **Title:** "Three-dimensional alpha shapes"
- **Venue:** ACM Transactions on Graphics, 13(1):43–72
- **DOI:** 10.1145/174462.156635
- **Contribution:** Introduced weighted alpha shapes, providing a dual combinatorial structure for unions of balls.

**(b) Union of balls and its dual shape:**
- **Authors:** Edelsbrunner, Herbert
- **Year:** 1995
- **Title:** "The union of balls and its dual shape"
- **Venue:** Discrete & Computational Geometry, 13(3–4):415–440
- **DOI:** 10.1007/BF02574053
- **URL (ACM):** https://dl.acm.org/doi/abs/10.1145/160985.161139
- **Free PDF:** https://pub.ista.ac.at/~edels/Papers/1995-02-UnionBallsDualShape.pdf
- **Contribution:** Efficient algorithms for computing topological, combinatorial, and metric properties of the union of finitely many spherical balls in ℝ^d. Based on a simplicial complex dual to a Voronoi-cell decomposition of the union. Established the homotopy equivalence of the union of balls and its dual shape (weighted alpha complex). This is the theoretical foundation for topology-correct SES computation.

**(c) Book: Geometry and Topology for Mesh Generation:**
- **Authors:** Edelsbrunner, Herbert
- **Year:** 2001
- **Title:** "Geometry and Topology for Mesh Generation"
- **Publisher:** Cambridge University Press
- **ISBN:** 0521793098
- **URL:** https://www.cambridge.org/core/books/geometry-and-topology-for-mesh-generation/261548A36951B9F3E7E26506C26A447A
- **Contribution:** Comprehensive treatment of alpha shapes, weighted alpha shapes, and their application to molecular surface mesh generation. Covers the relationship between the reduced surface and alpha shapes.

### 2.7 Sanner, Olson & Spehner (1996) — MSMS algorithm

- **Authors:** Sanner, Michel F.; Olson, Arthur J.; Spehner, Jean-Claude
- **Year:** 1996
- **Title:** "Reduced Surface: An Efficient Way to Compute Molecular Surfaces"
- **Venue:** Biopolymers, 38(3):305–320
- **DOI:** 10.1002/(SICI)1097-0282(199603)38:3<305::AID-BIP4>3.0.CO;2-Y
- **URL (Wiley):** https://onlinelibrary.wiley.com/doi/pdf/10.1002/%28SICI%291097-0282%28199603%2938%3A3%3C305%3A%3AAID-BIP4%3E3.0.CO%3B2-Y
- **Free PDF:** https://biomolmd.org/mw/images/4/49/Sanner96.pdf
- **Contribution:** Defined the "reduced surface" (RS) and presented the MSMS program. The RS is equivalent to a weighted alpha shape with alpha = 0 (or alpha = probe radius in some formulations). Four algorithms: (1) compute reduced surface, (2) build analytical SES from RS, (3) remove self-intersecting parts (singularities), (4) triangulate the SES. Complexity O(n log n) for n atoms. Tested on 709 PDB molecules (33 to 43,632 atoms). MSMS computes topologically correct surfaces where Connolly's PQMS failed on large molecules.

**MSMS software page:** https://ccsb.scripps.edu/msms/ and https://sannerlab.github.io/software/msms/

### 2.8 Voronoi-based SES methods

**(a) Eisenhaber & Argos (1993):**
- **Authors:** Eisenhaber, Frank; Argos, Patrick
- **Year:** 1993
- **Title:** "Improved strategy in analytic surface calculation for molecular systems: handling of singularities and computational efficiency"
- **Venue:** Journal of Computational Chemistry, 14(11):1272–1280
- **DOI:** 10.1002/jcc.540141103
- **Contribution:** Addressed singularities (nearly collinear triples of neighboring atoms, multiple vertices) in analytic surface calculations. Used Voronoi-based approaches.

**(b) Boissonnat, Devillers, Duquesne & Yvinec (1994):**
- **Authors:** Boissonnat, J.-D.; Devillers, O.; Duquesne, J.; Yvinec, M.
- **Year:** 1994
- **Title:** "Computing Connolly surfaces"
- **Venue:** Journal of Molecular Graphics, 12(1):61–62
- **DOI:** 10.1016/0263-7855(94)80033-2
- **URL:** https://inserm.hal.science/INRIA/hal-00795091v1
- **Contribution:** Voronoi-based computation of Connolly surfaces.

**(c) Liang, Edelsbrunner, Fu, Sudhakar & Subramaniam (1998):**
- **Authors:** Liang, J.; Edelsbrunner, H.; Fu, P.; Sudhakar, P.; Subramaniam, S.
- **Year:** 1998
- **Title:** "Analytical shape computation of macromolecules: I. molecular area and volume through alpha shape"
- **Venue:** Proteins, 33(1):1–17
- **DOI:** 10.1002/(SICI)1097-0134(19981001)33:1<1::AID-PROT1>3.0.CO;2-O
- **URL (PubMed):** https://pubmed.ncbi.nlm.nih.gov/9741840/
- **Free PDF:** https://pub.ista.ac.at/~edels/Papers/1998-J-04-AnalyticalShapeComputationI.pdf
- **Contribution:** Used alpha shape theory (Edelsbrunner & Mücke 1994) for exact molecular area and volume computation. The reduced surface is equivalent to an alpha shape with alpha value equal to the probe radius. Foundation for CASTp pocket detection.

**(d) Voss & Gerstein (2010) — 3V:**
- **Authors:** Voss, Neil R.; Gerstein, Mark
- **Year:** 2010
- **Title:** "3V: cavity, channel and cleft volume calculator and extractor"
- **Venue:** Nucleic Acids Research, 38(Web Server issue):W555–W562
- **DOI:** 10.1093/nar/gkq395
- **URL:** http://papers.gersteinlab.org/papers/3v/index-all.html
- **Contribution:** Grid-based tool for computing cavity/channel volumes using the SES. Uses the concept of the SES with probe spheres to identify internal voids.

**(e) Lin, Edelsbrunner & Fu (1994) — Sparse critical points:**
- **Authors:** Lin, S.L.; Edelsbrunner, H.; Fu, P.
- **Year:** 1994
- **Title:** "Molecular surface representations by sparse critical points"
- **Venue:** Proteins, 18(1):94–101
- **DOI:** 10.1002/prot.340180111
- **URL:** https://zenodo.org/records/1229355
- **Contribution:** Surface representation using a limited number of critical points at key locations. Critical points characterize surface topology (genus, connected components).

### 2.9 Quan & Stamm (2016) — Complete mathematical characterization

- **Authors:** Quan, Chaoyu; Stamm, Benjamin
- **Year:** 2016
- **Title:** "Mathematical analysis and calculation of molecular surfaces"
- **Venue:** Journal of Computational Physics, 322:221–239
- **DOI:** 10.1016/j.jcp.2016.07.007
- **Free PDF:** https://hal.sorbonne-universite.fr/hal-01205580/document
- **Contribution:** Complete characterization of the SES including singularities. Based on the signed distance function to the SAS. Derived the implicit function f_ses(p) = f_sas(p) + r_p. Characterized all SES singularities explicitly. Provided analytical area/volume formulas. Distinguished exterior SES (eSES) vs complete SES (cSES) including inner holes.

---

## 3. Algorithm Families and Their Tradeoffs

### 3.1 Analytic/parametric triangulation of SES

**How it works (Connolly 1983, Sanner 1996):**

The SES is decomposed into smooth patches:
- **Contact patches:** Spherical patches on atomic vdW spheres, bounded by circular arcs.
- **Toric reentrant patches:** Toroidal patches (saddle-shaped), each associated with a pair of atoms.
- **Spheric reentrant patches:** Spherical triangles on the probe sphere, associated with triples (or more) of atoms.

The reduced surface (RS) is a simplicial complex dual to the SES. Each RS face corresponds to a probe position in contact with ≥3 atoms (a spheric reentrant face). RS edges correspond to probe positions touching exactly 2 atoms (toric reentrant). RS vertices correspond to single-atom contact points.

**Triangulation density control (MSMS):**
- MSMS uses a `-density` parameter (default 1.0 vertices/Å²).
- Typical values: 1.0 for large molecules (>1000 atoms), 3.0 for smaller molecules.
- The triangulation density is the target number of vertices per unit surface area. MSMS generates vertices on the analytical surface patches and connects them into triangles.
- [MSMS documentation; https://ccsb.scripps.edu/msms/documentation/]

**Known issues (MSMS documented bugs):**
- Singular edges forming full circles are not treated.
- The genus of the reduced surface is sometimes wrong.
- Normal vectors of singular vertices point arbitrarily.
- MSMS crashes or produces wrong results if one sphere lies entirely inside another.
- If all vertices of a RS component also belong to another component, MSMS fails to find it.
- Self-intersecting parts of the SES ("singularities") must be removed; poorly handled singular cases cause triangulation failure and holes.
- By default MSMS restarts up to 5 times on failure, increasing problematic atom radii by 0.1 Å each time.
- [MSMS documentation, "Known bugs"; https://ccsb.scripps.edu/msms/documentation/]

**Performance (Sanner et al. 1996):**
- Crambin (327 atoms, 4772 triangles): 0.73 s on HP 9000/735
- Thermolysin (2437 atoms, 26462 triangles): 4.6 s
- Glutamine synthetase (43632 atoms, 476665 triangles): 104.53 s
- Complexity: O(n log n) for n atoms.
- [Sanner et al. 1996; https://biomolmd.org/mw/images/4/49/Sanner96.pdf]

**Strengths:** Exact analytical patches; correct topology when no singularities; good for surface area computation.
**Weaknesses:** Complex singularity handling; can fail on pathological geometries; triangulation density control is somewhat indirect.

### 3.2 Implicit/grid methods

**How it works:**
1. Define an implicit function (usually a signed distance function or indicator function) for the SES.
2. Sample it on a 3D grid.
3. Extract the isosurface using Marching Cubes (or variants: dual contouring, surface nets).

**Key implicit function (Quan & Stamm 2016):**
- f_sas(p) = signed distance from p to the SAS
- f_ses(p) = f_sas(p) + r_p
- SES = {p : f_ses(p) = 0} = {p : f_sas(p) = −r_p}

The SAS is the boundary of the union of inflated spheres B(c_i, r_i + r_p). The signed distance to the SAS can be computed analytically using a Voronoi-type diagram of the inflated spheres. [Quan & Stamm 2016; https://hal.sorbonne-universite.fr/hal-01205580/document]

**EDTSurf (Xu & Zhang 2009):**
- **Authors:** Xu, Dong; Zhang, Yang
- **Year:** 2009
- **Title:** "Generating Triangulated Macromolecular Surfaces by Euclidean Distance Transform"
- **Venue:** PLoS ONE, 4(12):e8140
- **DOI:** 10.1371/journal.pone.0008140
- **URL:** https://zhanggroup.org/EDTSurf/
- Uses fast Euclidean Distance Transform (EDT) on a grid, then Vertex-Connected Marching Cubes.
- Relative error vs analytical: 2–4% depending on grid resolution.
- 1.5–4× lower error than prior grid methods.
- [EDTSurf paper; https://seq2fun.dcmb.med.umich.edu/EDTSurf/EDTSurf.pdf]

**Resolution vs accuracy:**
- Grid spacing h determines the geometric error. For a signed distance function (1-Lipschitz), the isosurface extracted by marching cubes has Hausdorff distance O(h) from the true surface (linear interpolation on edges gives O(h²) for smooth surfaces, but the isosurface topology can be wrong if features are smaller than h).
- Typical grid spacings for molecular surfaces: 0.2–0.5 Å. (See §6.)
- Thin necks and narrow passages can be lost if the grid spacing exceeds the feature width.

**Strengths:** Robust, easy to implement, handles arbitrary geometries, no singularity handling needed.
**Weaknesses:** Grid resolution limits accuracy; memory O(1/h³); topology can be wrong for fine features; smooth approximation at sharp features.

### 3.3 Ray-casting methods

**NanoShaper (Decherchi & Rocchia 2013):**
- **Authors:** Decherchi, Sergio; Rocchia, Walter
- **Year:** 2013
- **Title:** "A general and robust ray-casting-based algorithm for triangulating surfaces at the nanoscale"
- **Venue:** PLoS ONE, 8(4):e59744
- **DOI:** 10.1371/journal.pone.0059744
- **URL:** https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0059744
- **Free (PMC):** https://pmc.ncbi.nlm.nih.gov/articles/PMC3618509/
- Builds molecular surfaces according to several definitions (SES/Connolly-Richards, Skin, Gaussian) within a common framework.
- Uses ray casting + Marching Cubes on analytically sampled vertices.
- Robust to rotation, scale, and atom displacement.
- Handles fullerene cavity detection where MSMS fails.
- NanoShaper tutorial: Decherchi & Rocchia, "Building and Analyzing Molecular Surfaces: A Tutorial on NanoShaper", Methods in Molecular Biology (2015); DOI: 10.1007/978-3-319-12211-3_10

**Ray-casting approach:** Cast rays through the grid; the intersection of a ray with the SES gives a point where the implicit function changes sign. By casting many rays and locating sign changes, one obtains surface samples that can be triangulated. The ray-casting step is used to "color" the grid (inside/outside classification), after which Marching Cubes extracts the isosurface.

**Strengths:** Robust, works for arbitrary surface definitions, handles complex topologies.
**Weaknesses:** Depends on ray density; may miss features not hit by any ray; more complex than pure grid methods.

### 3.4 Gaussian/surface-switching smooth approximations (SWIG)

**SWIG algorithm (Lange & Herbert 2013–2020):**
- **Authors:** Lange, Adrian W.; Herbert, John M.; Albrecht, Benjamin J.; You, Zhi-Qiang
- **Year:** 2019/2020
- **Title:** "Intrinsically smooth discretisation of Connolly's solvent-excluded molecular surface"
- **Venue:** Molecular Physics, 118(6):e1644384
- **DOI:** 10.1080/00268976.2019.1644384
- **URL:** https://www.tandfonline.com/doi/full/10.1080/00268976.2019.1644384
- **Free PDF:** https://par.nsf.gov/servlets/purl/10276515 or https://www.asc.ohio-state.edu/herbert.44/reprints/MP_118_e1644384.pdf

**The SWIG approach:**
1. Start with atom-centered Lebedev quadrature grids on the SAS (or vdW surface).
2. Use a **switching function** to smoothly attenuate grid points as they enter the cavity interior (rather than discarding them abruptly).
3. Replace point charges with **Gaussian functions** to avoid oscillations when points approach each other:
   g_i(r) = q_i (ζ²_i/π)^{3/2} exp(−ζ²_i |r − s_i|²)
4. For SES-SWIG: project SAS discretization points back onto the vdW surface (this automatically eliminates points that should be replaced by reentrant surface points), then add reentrant surface points analytically using Connolly's prescription.

**Properties:**
- The analytic SES is smooth (C¹ continuous) — its first derivatives are continuous.
- SES-SWIG provides a numerically smooth discretization: solvation energy varies smoothly with nuclear coordinates and is analytically differentiable.
- Avoids the discontinuities that plague finite-difference Poisson-Boltzmann solvers on the SES.
- [Lange et al. 2019/2020; https://par.nsf.gov/servlets/purl/10276515]

**Related smooth approximations:**
- **Gaussian surfaces** (Grant & Pickup 1995): Sum of atom-centered Gaussians as an isosurface. Differentiable, free of self-intersections. [J Phys Chem 99:3503–3510; DOI: 10.1021/j100011a016]
- **Skin surfaces** (Edelsbrunner 1999): Smooth blending of spheres. [Edelsbrunner, "Molecular modeling and bioinformatics", 1999]
- **GBMV** (generalized Born molecular volume): Uses switching functions to smooth SES cusps. [Blessing & Friedrich 2006 / Labute 2008 — see Lange et al. references]

---

## 4. The Exact Implicit Function for SES (Key for Marching Cubes)

### 4.1 Formalization

Let the molecule consist of M atoms with centers c_i and vdW radii r_i. Let r_p be the probe radius.

**Step 1: Union of atomic balls**
Ω_vdW = ⋃ B(c_i, r_i)

**Step 2: Union of inflated balls (SAS region)**
Ω_SAS = ⋃ B(c_i, r_i + r_p)

The SAS = ∂Ω_SAS.

**Step 3: Signed distance to SAS**
f_sas(p) = signed distance from p to ∂Ω_SAS
- f_sas(p) < 0 if p is inside Ω_SAS (probe center cannot be here)
- f_sas(p) > 0 if p is outside Ω_SAS (probe center can be here)
- f_sas(p) = 0 on the SAS

**Step 4: SES implicit function**
f_ses(p) = f_sas(p) + r_p

The SES = {p : f_ses(p) = 0} = {p : f_sas(p) = −r_p}.

[Quan & Stamm 2016, Eq. 3.4; https://hal.sorbonne-universite.fr/hal-01205580/document]

### 4.2 Why the simple min formula is NOT the signed distance

A naive approach would use:
f̃_sas(p) = min_i {‖p − c_i‖ − (r_i + r_p)}

This is the *implicit function* of the SAS (zero level set = SAS), but it is NOT a signed distance function. Its gradient magnitude is not 1. Using it directly in f_ses = f̃_sas + r_p does NOT give the correct SES, because the SES is the set at signed distance −r_p from the SAS, not the set where f̃_sas = −r_p.

Quan & Stamm (2016) show that computing the true signed distance to the SAS requires a Voronoi-type diagram of the inflated spheres, with three equivalence statements giving the closest point on the SAS. They provide a constructive algorithm. [Quan & Stamm 2016, §4–5; https://hal.sorbonne-universite.fr/hal-01205580/document]

**Important caveat (Quan & Stamm, Eq. 3.5):** A representation proposed by Pomelli & Tomasi (implicit function as min over triples of atoms) can fail — it may exclude valid SES regions. [Quan & Stamm 2016, Figure 2]

### 4.3 Alternative: SES as boundary of probe-swept region

The SES is the boundary of:
  Ω_SES = Ω_vdW ∪ {p : dist(p, Ω_vdW) ≤ r_p and the probe at p would collide with Ω_vdW}

More precisely, the SES is the boundary of the set of points that the probe sphere *body* occupies as it rolls over the vdW surface. This equals the boundary of the complement of the "excluded region" — the region where the probe center cannot go.

Equivalently: SES = ∂(ℝ³ \ U) where U is the set of probe-center positions that are outside the inflated union AND outside all probe-body regions that touch the vdW surface from inside.

The cleanest formulation remains the signed-distance one: SES = {p : f_sas(p) = −r_p}.

### 4.4 Lipschitz property and grid error bounds

**Signed distance functions are 1-Lipschitz:**
|f_sas(p) − f_sas(q)| ≤ ‖p − q‖ for all p, q.

This follows directly from the definition of distance: |dist(p, S) − dist(q, S)| ≤ ‖p − q‖. The signed version preserves this (the sign is constant on each side, and the absolute value is the unsigned distance). [Wikipedia "Signed distance function"; https://en.wikipedia.org/wiki/Signed_distance_function]

**Error bound for marching cubes:**
For a 1-Lipschitz function f sampled on a grid of spacing h:
- If f is linearly interpolated along grid edges, the zero-crossing point on an edge has error O(h²) for a smooth surface (since f is smooth and the linear approximation has second-order error).
- The Hausdorff distance between the extracted mesh and the true zero level set is O(h) in general (the mesh vertices lie on the true surface if exact edge crossings are used; the mesh faces deviate by O(h²) for smooth regions).
- **Topology guarantee:** If the zero level set has no features (necks, tunnels) narrower than ~2h, marching cubes preserves the topology. Features smaller than h can be missed or merged.

**For SES specifically:**
- The SES is piecewise smooth (spherical, toroidal patches meeting at circular arcs). At the arcs (where patches meet), the surface normal is continuous (C¹) but the curvature is discontinuous.
- Grid spacing h should be small compared to the smallest feature: the narrowest crevice or tunnel that must be resolved.
- For r_p = 1.4 Å and typical protein features, h = 0.2–0.5 Å is standard (see §6).
- The error in surface position is ~O(h) in practice; EDTSurf reports 2–4% relative area error at standard resolutions. [Xu & Zhang 2009; https://seq2fun.dcmb.med.umich.edu/EDTSurf/EDTSurf.pdf]

**Practical rule of thumb:**
- To resolve a feature of width w, need h ≤ w/3 at minimum, preferably h ≤ w/5.
- The SES reentrant tori have minor radius r_p = 1.4 Å; to resolve the toroidal curvature, h ≤ 0.3 Å is desirable.
- The signed distance function is 1-Lipschitz, so the isosurface error scales linearly with h.

---

## 5. Topology Correctness

### 5.1 Genus and connected components

The topology of the SES (number of connected components, genus/number of handles) is determined by the arrangement of atoms and the probe radius. For a calibrated set of n spheres (r_max < 2·r_min and center distances > r_min/2) with an adequate probe radius (r_min < r_p < r_max), the reduced surface has O(n) faces and edges. [Sanner et al. 1996; https://biomolmd.org/mw/images/4/49/Sanner96.pdf]

Edelsbrunner (1995) established the homotopy equivalence between the union of balls and its dual weighted alpha complex, giving algorithms to compute homology groups (and hence genus) of the union. [Edelsbrunner 1995; https://pub.ista.ac.at/~edels/Papers/1995-02-UnionBallsDualShape.pdf]

### 5.2 When grid methods lose thin necks

Grid methods (marching cubes, EDTSurf, NanoShaper) can lose topological features when:
- A tunnel or neck has diameter < ~2h (grid spacing). The grid cannot resolve the passage.
- Two surface sheets come within < h of each other but do not merge — the grid may incorrectly merge them.
- A thin cavity wall has thickness < h — the grid may punch a hole.

This is a fundamental limitation of grid sampling: features smaller than the Nyquist limit (2h) cannot be resolved. To capture all SES topology features, h must be smaller than half the smallest feature width. For proteins, the smallest features are typically the narrowest crevices (where two atoms nearly touch), which can be < 1 Å wide.

### 5.3 When MSMS fails

MSMS documented failure modes (https://ccsb.scripps.edu/msms/documentation/):
1. **Singularities:** Self-intersecting parts of the SES. When the probe touches 4+ atoms simultaneously, the SES patches can self-intersect. MSMS removes these but poorly handled cases cause holes.
2. **Singular edges forming full circles:** Not treated.
3. **Sphere entirely inside another sphere:** Crash or wrong results.
4. **RS component detection failure:** If all vertices of a component belong to another component.
5. **Genus errors:** The genus of the reduced surface is sometimes wrong.
6. **Near-degenerate geometries:** Nearly collinear atom triples cause numerical instabilities (Eisenhaber & Argos 1993).

### 5.4 Topology guarantees

**What IS guaranteed:**
- Edelsbrunner's weighted alpha complex computes the exact topology (homology) of the union of balls. If the SES is computed via the alpha shape dual (alpha = probe radius), the topology of the reduced surface is exact. [Edelsbrunner 1995; Liang et al. 1998]
- MSMS claims topologically correct surfaces for all 709 test molecules (33–43,632 atoms). [Sanner et al. 1996]

**What is NOT guaranteed:**
- Grid methods have no topological guarantee unless the grid resolves all features.
- MSMS has known bugs (see §5.3) and can produce incorrect genus.
- The complete characterization of SES singularities was only provided by Quan & Stamm (2016); earlier methods handled special cases.

**Practical recommendation for the Rust library:**
- For topology-critical applications, use the alpha-shape / reduced-surface approach (Edelsbrunner) which gives exact topology.
- For robustness, use implicit/grid methods with sufficient resolution and verify topology (Euler characteristic, connected components) against the alpha-shape result.
- Always compute Euler characteristic χ = V − E + F and compare with expected value 2(1 − g) per connected component, where g is the genus.

---

## 6. Typical Parameters

### 6.1 Probe radius

- **Default:** r_p = 1.4 Å (radius of a water molecule, approximated as a sphere).
  - Richards (1977): "About the smallest physically reasonable probe is a water molecule considered as a sphere of radius of 1.4 or 1.5 Å."
  - [Richards 1977; https://biomolmd.org/mw/images/e/e5/Richards77.pdf]
- **MSMS default:** 1.5 Å (MSMS uses 1.5 by default; the standard in the field is 1.4 Å for water).
  - [MSMS documentation; https://ccsb.scripps.edu/msms/documentation/]
- **Alternative values:**
  - 1.2 Å: sometimes used in implicit solvation (PCM) contexts. [Quan & Stamm 2016 use r_p = 1.2 Å in examples]
  - 1.5 Å: MSMS default; also used in some older literature.
  - 0 Å: degenerate case = vdW surface.
- **Water probe justification:** 1.4 Å ≈ half the distance to the first peak in the oxygen–oxygen radial distribution function of neat water. [Lange et al. 2019/2020, citing Ref. 6; https://par.nsf.gov/servlets/purl/10276515]

### 6.2 Grid spacing (implicit methods)

- **Typical range:** 0.2–0.5 Å.
- EDTSurf reports 2–4% surface area error depending on grid resolution. [Xu & Zhang 2009]
- For finite-difference Poisson-Boltzmann solvers, grid spacings of 0.2–0.5 Å are standard (finer grids give better accuracy but higher memory cost).
- NanoShaper uses grid-based ray casting; typical grid resolutions are similar.
- **Recommendation:** h = 0.3 Å as a good default; h = 0.2 Å for high accuracy; h = 0.5 Å for fast preview.

### 6.3 Triangle counts (MSMS)

From Sanner et al. 1996 (at density = 1.0 vertex/Å²):
- Crambin (327 atoms, small protein): 4,772 triangles
- Thermolysin (2437 atoms, medium protein): 26,462 triangles
- Glutamine synthetase (43,632 atoms, large complex): 476,665 triangles

**General scaling:** Triangle count ≈ surface_area × density. For a protein of N atoms, the SES surface area is roughly 5–10 N Å² (depending on compactness). At density 1.0:
- N = 100 atoms → ~5,000–10,000 triangles
- N = 1,000 atoms → ~50,000–100,000 triangles
- N = 10,000 atoms → ~500,000–1,000,000 triangles
- N = 40,000 atoms → ~2,000,000–4,000,000 triangles

MSMS output format: `.vert` and `.face` files. The face file classifies each triangle as:
- 1 = toric reentrant face
- 2 = spheric reentrant face
- 3 = contact face

[MSMS documentation; https://ccsb.scripps.edu/msms/documentation/]

### 6.4 MSMS density parameter

- Default: 1.0 vertex/Å²
- Typical: 1.0 for molecules > 1000 atoms; 3.0 for smaller molecules.
- Values from 0.5 to 10 are usable; the value is not validated by MSMS.
- [MSMS documentation; https://ccsb.scripps.edu/msms/documentation/]

### 6.5 van der Waals radii

Standard vdW radii (Bondi 1964 or similar) are used. The choice of radii set affects the surface. Common choices:
- Bondi (1964) radii
- Jorgensen & Severance (1990) TIP3P-compatible radii
- Element-specific radii from the PDB chemical component dictionary

The SES is sensitive to the radius set — changing radii changes the topology and area.

---

## Bibliography

### Primary papers

1. **Lee, B.; Richards, F.M.** (1971). "The interpretation of protein structures: estimation of static accessibility." *Journal of Molecular Biology*, 55(3):379–400. DOI: [10.1016/0022-2836(71)90324-X](https://doi.org/10.1016/0022-2836(71)90324-X)

2. **Shrake, A.; Rupley, J.A.** (1973). "Environment and exposure to solvent of protein atoms. Lysozyme and insulin." *Journal of Molecular Biology*, 79(2):351–371. DOI: [10.1016/0022-2836(73)90011-9](https://doi.org/10.1016/0022-2836(73)90011-9)

3. **Richards, F.M.** (1977). "Areas, volumes, packing and protein structure." *Annual Review of Biophysics and Bioengineering*, 6:151–176. DOI: [10.1146/annurev.bb.06.060177.001055](https://doi.org/10.1146/annurev.bb.06.060177.001055). Free PDF: [biomolmd.org](https://biomolmd.org/mw/images/e/e5/Richards77.pdf)

4. **Connolly, M.L.** (1983). "Analytical molecular surface calculation." *Journal of Applied Crystallography*, 16(5):548–558. DOI: [10.1107/S0021889883010985](https://doi.org/10.1107/S0021889883010985)

5. **Connolly, M.L.** (1983). "Solvent-Accessible Surfaces of Proteins and Nucleic Acids." *Science*, 221(4612):709–713. DOI: [10.1126/science.6879170](https://doi.org/10.1126/science.6879170). Free PDF: [gersteinlab.org](http://www.gersteinlab.org/courses/452/09-spring/pdf/Connolly.pdf)

6. **Richmond, T.J.** (1984). "Solvent accessible surface area and excluded volume in proteins. Analytical equations for overlapping spheres and implications for the hydrophobic effect." *Journal of Molecular Biology*, 178(1):63–89. DOI: [10.1016/0022-2836(84)90231-6](https://doi.org/10.1016/0022-2836(84)90231-6)

7. **Edelsbrunner, H.; Mücke, E.P.** (1994). "Three-dimensional alpha shapes." *ACM Transactions on Graphics*, 13(1):43–72. DOI: [10.1145/174462.156635](https://doi.org/10.1145/174462.156635)

8. **Edelsbrunner, H.** (1995). "The union of balls and its dual shape." *Discrete & Computational Geometry*, 13(3–4):415–440. DOI: [10.1007/BF02574053](https://doi.org/10.1007/BF02574053). Free PDF: [pub.ista.ac.at](https://pub.ista.ac.at/~edels/Papers/1995-02-UnionBallsDualShape.pdf)

9. **Sanner, M.F.; Olson, A.J.; Spehner, J.-C.** (1996). "Reduced Surface: An Efficient Way to Compute Molecular Surfaces." *Biopolymers*, 38(3):305–320. DOI: [10.1002/(SICI)1097-0282(199603)38:3<305::AID-BIP4>3.0.CO;2-Y](https://onlinelibrary.wiley.com/doi/pdf/10.1002/%28SICI%291097-0282%28199603%2938%3A3%3C305%3A%3AAID-BIP4%3E3.0.CO%3B2-Y). Free PDF: [biomolmd.org](https://biomolmd.org/mw/images/4/49/Sanner96.pdf)

10. **Eisenhaber, F.; Argos, P.** (1993). "Improved strategy in analytic surface calculation for molecular systems: handling of singularities and computational efficiency." *Journal of Computational Chemistry*, 14(11):1272–1280. DOI: [10.1002/jcc.540141103](https://doi.org/10.1002/jcc.540141103)

11. **Boissonnat, J.-D.; Devillers, O.; Duquesne, J.; Yvinec, M.** (1994). "Computing Connolly surfaces." *Journal of Molecular Graphics*, 12(1):61–62. DOI: [10.1016/0263-7855(94)80033-2](https://doi.org/10.1016/0263-7855(94)80033-2)

12. **Lin, S.L.; Edelsbrunner, H.; Fu, P.** (1994). "Molecular surface representations by sparse critical points." *Proteins*, 18(1):94–101. DOI: [10.1002/prot.340180111](https://doi.org/10.1002/prot.340180111)

13. **Liang, J.; Edelsbrunner, H.; Fu, P.; Sudhakar, P.; Subramaniam, S.** (1998). "Analytical shape computation of macromolecules: I. molecular area and volume through alpha shape." *Proteins*, 33(1):1–17. DOI: [10.1002/(SICI)1097-0134(19981001)33:1<1::AID-PROT1>3.0.CO;2-O](https://onlinelibrary.wiley.com/doi/abs/10.1002/%28SICI%291097-0134%2819981001%2933%3A1%3C1%3A%3AAID-PROT1%3E3.0.CO%3B2-O). Free PDF: [pub.ista.ac.at](https://pub.ista.ac.at/~edels/Papers/1998-J-04-AnalyticalShapeComputationI.pdf)

14. **Xu, D.; Zhang, Y.** (2009). "Generating Triangulated Macromolecular Surfaces by Euclidean Distance Transform." *PLoS ONE*, 4(12):e8140. DOI: [10.1371/journal.pone.0008140](https://doi.org/10.1371/journal.pone.0008140). Free PDF: [seq2fun.dcmb.med.umich.edu](https://seq2fun.dcmb.med.umich.edu/EDTSurf/EDTSurf.pdf)

15. **Voss, N.R.; Gerstein, M.** (2010). "3V: cavity, channel and cleft volume calculator and extractor." *Nucleic Acids Research*, 38(Web Server issue):W555–W562. DOI: [10.1093/nar/gq395](https://doi.org/10.1093/nar/gq395)

16. **Decherchi, S.; Rocchia, W.** (2013). "A general and robust ray-casting-based algorithm for triangulating surfaces at the nanoscale." *PLoS ONE*, 8(4):e59744. DOI: [10.1371/journal.pone.0059744](https://doi.org/10.1371/journal.pone.0059744). Free (PMC): [pmc.ncbi.nlm.nih.gov](https://pmc.ncbi.nlm.nih.gov/articles/PMC3618509/)

17. **Quan, C.; Stamm, B.** (2016). "Mathematical analysis and calculation of molecular surfaces." *Journal of Computational Physics*, 322:221–239. DOI: [10.1016/j.jcp.2016.07.007](https://doi.org/10.1016/j.jcp.2016.07.007). Free PDF: [hal.sorbonne-universite.fr](https://hal.sorbonne-universite.fr/hal-01205580/document)

18. **Lange, A.W.; Herbert, J.M.; Albrecht, B.J.; You, Z.-Q.** (2019/2020). "Intrinsically smooth discretisation of Connolly's solvent-excluded molecular surface." *Molecular Physics*, 118(6):e1644384. DOI: [10.1080/00268976.2019.1644384](https://doi.org/10.1080/00268976.2019.1644384). Free PDF: [par.nsf.gov](https://par.nsf.gov/servlets/purl/10276515)

19. **Connolly, M.L.** (1993). "The molecular surface package." *Journal of Molecular Graphics*, 11(2):139–141. DOI: [10.1016/0263-7855(93)87010-3](https://doi.org/10.1016/0263-7855(93)87010-3)

### Books and reviews

20. **Edelsbrunner, H.** (2001). *Geometry and Topology for Mesh Generation*. Cambridge University Press. ISBN: 0521793098. URL: [cambridge.org](https://www.cambridge.org/core/books/geometry-and-topology-for-mesh-generation/261548A36951B9F3E7E26506C26A447A)

21. **Gerstein, M.; Richards, F.S.** (2001). "Protein geometry: volumes, areas and distances." In *International Tables for Crystallography*. DOI: [10.1107/97809553602060000885](https://doi.org/10.1107/97809553602060000885)

### Software documentation

22. **MSMS documentation.** Sanner Lab, The Scripps Research Institute. URL: [ccsb.scripps.edu/msms/documentation/](https://ccsb.scripps.edu/msms/documentation/)

23. **NanoShaper tutorial.** Decherchi, S.; Rocchia, W. "Building and Analyzing Molecular Surfaces: A Tutorial on NanoShaper." *Methods in Molecular Biology* (2015). DOI: [10.1007/978-3-319-12211-3_10](https://doi.org/10.1007/978-3-319-12211-3_10)

### Additional references (from Wikipedia and secondary sources)

24. **Wikipedia: Accessible surface area.** URL: [en.wikipedia.org/wiki/Accessible_surface_area](https://en.wikipedia.org/wiki/Accessible_surface_area)

25. **Wikipedia: Signed distance function.** URL: [en.wikipedia.org/wiki/Signed_distance_function](https://en.wikipedia.org/wiki/Signed_distance_function)

26. **Wikipedia: Van der Waals surface.** URL: [en.wikipedia.org/wiki/Van_der_Waals_surface](https://en.wikipedia.org/wiki/Van_der_Waals_surface)

27. **Bondi, A.** (1964). "van der Waals volumes and radii." *Journal of Physical Chemistry*, 68(3):441–451. DOI: [10.1021/j100785a001](https://doi.org/10.1021/j100785a001) *(vdW radii standard — cited here as the standard radii source)*

28. **Petitjean, M.** (1994). "On the Analytical Calculation of van der Waals Surfaces and Volumes: Some Numerical Aspects." *Journal of Computational Chemistry*, 15(5):507–523. DOI: [10.1002/jcc.540150504](https://doi.org/10.1002/jcc.540150504)

29. **Grant, J.A.; Pickup, B.T.** (1995). "A Gaussian description of molecular shape." *Journal of Physical Chemistry*, 99(11):3503–3510. DOI: [10.1021/j100011a016](https://doi.org/10.1021/j100011a016)

---

## Uncertainty and Open Items

1. **Edelsbrunner 1995 PDF text extraction failed** (likely scanned image). The content was summarized from the abstract and secondary sources. The key result (homotopy equivalence of union of balls and dual shape) is well-established.

2. **The exact formulation of the reduced surface as an alpha shape** is stated in Sanner et al. 1996 as "equivalent to an alpha-shape with an alpha value equal to the probe radius" but Liang et al. 1998 and EDTSurf say "alpha value equal to zero." The precise alpha value convention may differ between authors. Marked as uncertain.

3. **MSMS default probe radius is 1.5 Å** (per documentation) while the standard water probe is 1.4 Å. This discrepancy is noted in the documentation itself.

4. **Triangle count estimates for general N** (§6.3) are rough extrapolations from the three data points in Sanner et al. 1996, not from a cited scaling law.

5. **The Pomelli & Tomasi implicit function** mentioned in Quan & Stamm (2016) as Eq. 3.5 was not independently verified; its failure mode is documented in Quan & Stamm Figure 2.

6. **Voss et al. Voronoi-based SES** — the specific Voss/Voronoi algorithm papers (beyond 3V) were not fully retrieved; the 3V tool and vossvolvox GitHub are documented but the underlying algorithm papers should be consulted directly.

7. **Grid spacing 0.2–0.5 Å** is stated as "typical" based on EDTSurf and general PB-solver practice; a rigorous source for this specific range was not found in the retrieved papers.
