# SES / SAS Modern Implementations, Performance, and Mesh Generation

Research stream 02 for the Rust SES library+CLI product spec. Every factual claim carries a source
(title + URL). Uncertainty is marked explicitly. Last updated: 2026 (this research pass).

Scope: modern (and still-used) SES/SAS software, grid-based algorithm practice, mesh quality,
concurrency/vectorization patterns, file I/O conventions, and Rust ecosystem.

---

## 1. Modern SES/SAS software: algorithms, status, limitations

### 1.1 MSMS (Sanner) — analytic Reduced-Surface SES, still the reference

- **What it is.** MSMS computes the *Reduced Surface* (RS, a combinatorial skeleton of the
  rolling-probe configuration) and from it the analytic SES patches (convex spherical,
  concave spherical "reentrant", and toroidal saddle patches), then triangulates the patches
  at a prescribed vertex density. Author: Michel F. Sanner (PhD thesis, Scripps).
  Source: MSMS home page — https://ccsb.scripps.edu/msms/ ; Sanner Lab —
  https://sannerlab.github.io/software/msms/
- **Canonical citation.** Sanner, M. F., Olson, A. J. & Spehner, J.-C. (1996). "Reduced
  Surface: An Efficient Way to Compute Molecular Surfaces." *Biopolymers* 38:305–320.
  doi:10.1002/(SICI)1097-0282(199603)38:3<305::AID-BIP4>3.0.CO;2-Y
  (Also the 1995 SCG abstract "Fast and Robust Computation of Molecular Surfaces".)
  As of Nov 2026 the 1996 paper has 2,721 citations (Google Scholar, per the MSMS page).
- **Algorithm family.** Analytic / combinatorial (not grid-based). The four sub-algorithms of
  MSMS and their complexity are analyzed in Sanner96 (PDF:
  https://biomolmd.org/mw/images/4/49/Sanner96.pdf).
- **CLI parameters (documentation: https://sannerlab.github.io/software/msms/documentation.html).**
  - Input: `-if filename.xyzr` — one sphere per line: `x y z r` (Å). Since v2.6.2 duplicate
    centers are merged (largest kept).
  - `-probe_radius` (default **1.5 Å**; 0.5–10 safe range).
  - `-density` triangulation density (default **1.0 vertex/Å²**; recommended 1.0 for >1000
    atoms, 3.0 for smaller molecules). This is the "resolution/LOD" knob for MSMS.
  - `-of prefix` writes `prefix.vert` / `prefix.face` (per component: `prefix_3.vert` etc.).
  - `-af` per-atom SES and SAS areas; `-surface <tses,ases>` (triangulated vs analytic areas
    only); `-all_components` (cavities too); `-one_cavity` (single component via seed atoms).
  - `-noh` skips r=1.2 spheres; `-free_vertices`; error recovery: on failure MSMS restarts
    up to 5 times, inflating problem-atom radii by 0.1 Å (`-no_rest_on_pbr`, `-no_rest`).
  - Server mode over Unix socket (`-socket`, `-sinetd`) — historically used by VMD/Chimera.
  - Since v2.6.2, `-all_components` runs partly in parallel, controlled by `OMP_NUM_THREADS`.
- **Output format (important for our I/O design).** `.vert` line: `x y z nx ny nz face_number
  closest_sphere_index patch_type` (1=toric reentrant, 2=spheric reentrant, 3=contact).
  `.face` line: `v1 v2 v3 patch_type face_number`. Coordinates and normals are in Å, 1-based
  indices. Normals are **per-vertex** and point outward. Source: MSMS documentation (URL above).
- **Status.** Still maintained in the sense of being hosted and licensed: "MSMS turns 30!
  After 30 years, the MSMS water probe is still rolling." Free for non-commercial use under a
  "Sannerlab MIT-like Source-Available License"; commercial use requires a separate license.
  Used by VMD, Chimera, PMV, and by MaSIF (Gainza et al., Nature Methods 2020) for geometric
  deep learning. Source: https://ccsb.scripps.edu/msms/
- **Limitations (documented + comparative).**
  - Known bugs (MSMS docs): genus of reduced surface sometimes wrong; singular edges forming
    full circles unhandled; **normals of singular vertices point arbitrarily** to one probe
    center; crashes/wrong results if one sphere lies entirely inside another; a component may
    be missed if all its RS vertices also belong to another component (workaround `-one_cavity`).
  - Scalability: "MSMS crashes on most structures larger than 10,000 atoms" (Tom Goddard,
    "Molecular Surface Algorithms", Oct 2013 —
    https://www.cgl.ucsf.edu/chimera/data/surface-oct2013/surface.html).
  - Robustness: in NanoShaper tests, MSMS failed on the highly symmetric fullerene (C60)
    cavity detection and can fail to converge (it restarts from new probe seeds; "if the
    failure repeats itself a number of times, MSMS exits") — Decherchi & Rocchia 2013
    (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0059744).
  - Per-component restart makes the *total* process time variable; NanoShaper paper reports
    MSMS slower than EDTSurf and NanoShaper on their benchmark, with some runs not converging.
  - Positives: exact patches → analytic areas; careful triangulation → ~3× fewer triangles
    than a 0.5 Å distance-grid MC mesh at comparable visual quality (Goddard 2013, see §3).

### 1.2 NanoShaper (Decherchi & Rocchia) — ray-cast + alpha-shape analytic SES, CGAL-backed

- **What it is.** A C++ suite that builds molecular surfaces under several definitions
  (Connolly–Richards SES, VdW, SAS, Skin, Gaussian/"Blobby") in one framework, triangulates
  them by a ray-casting + "Analytical Intersections Marching Cubes" scheme, detects/fills
  cavities, estimates area/volume, colors grids for PDEs (DelPhi integration), and detects
  pockets. Implemented with CGAL (weighted Delaunay tetrahedralization + alpha shapes),
  exact predicates / inexact constructions, double precision.
  Sources: Decherchi & Rocchia, "A general and Robust Ray-Casting-Based Algorithm for
  Triangulating Surfaces at the Nanoscale", PLoS ONE 8(4):e59744, 2013 —
  https://pmc.ncbi.nlm.nih.gov/articles/PMC3618509/ ; NanoShaper 0.7 User Guide (PDF in
  https://github.com/shi-scala/NanoShaper and https://gitlab.iit.it/SDecherchi/nanoshaper).
- **SES construction (not grid-based for the surface shape).** The SES is built *analytically*
  using **alpha-shape theory** (Liang et al. 1998) over the weighted points (atom spheres):
  CGAL weighted Delaunay tetrahedralization + alpha shape gives the topologically correct
  patch decomposition; patch trimming spheres/planes follow Krone et al. 2009 / Chen et al.
  2010. Degenerate configurations are cured by random atom perturbation up to 0.0001 Å (the
  least significant PDB digit). Ray–sphere and ray–torus intersections are analytic (torus via
  Sturm-sequence root finding). VdW = SES with probe radius 0; SAS = VdW of radii inflated by
  the probe radius. (Decherchi & Rocchia 2013, Methods.)
- **Triangulation.** Grid-consistent "Analytical Intersections Marching Cubes": grid rays
  (3 orthogonal directions) give in/out coloring + volume + cavity flood-fill; edge rays give
  *analytic* surface/edge intersections so MC vertices lie on the surface (no interpolation);
  MC tables produce the mesh. A cheaper "Bisecting Marching Cubes" variant (grid rays only,
  1/3 the ray cost, less memory) trades accuracy for speed. Laplacian smoothing optional.
  Checksum criterion (even intersection count per ray) + random recast for numerical failures.
  (Decherchi & Rocchia 2013.)
- **Performance claims (all from Decherchi & Rocchia 2013 unless noted).**
  - Skin surface of a 100k-atom protein built ~3× faster than Lindow et al. single-threaded
    (hardware-adjusted); triangulated ≥10× faster than CGAL's Kruithof algorithm.
  - 1VSZ human adenovirus, 180,000 atoms (with H): SES + Skin triangulated (6.2M / 7.0M
    triangles at **scale = 2 grids/Å**) on a middle-range workstation. Skin surfacing overall
    234 s (17 s build + 88 s ray-cast + 129 s MC), peak RAM 7.3 GB (dev version with octrees:
    ~4 GB).
  - Vs. EDTSurf and MSMS: NanoShaper faster than EDTSurf, which is faster than MSMS on their
    benchmark; EDTSurf time ~constant because it scales with grid size (typically ~300³),
    not atom count. NanoShaper reaches EDTSurf accuracy at roughly half the grid scale.
  - NanoShaper 0.7 guide: triangulating SV40 (1sva) at scale 2.0 Å⁻¹ takes ~2 min including
    normals on a dual 8-core SandyBridge.
  - Grid scale of 2.0 grids/Å is stated as sufficient for good area/volume estimation and
    biophysically relevant features (smallest atomic radii ~1.2 Å).
- **Parallelism.** Ray casting parallelized over shared-memory CPUs (multithreaded). Config
  key `Number Thread`; guide recommends *4 × physical cores* for full CPU utilization.
- **Parameters (0.7 guide, §5).** `Grid_scale` (Å⁻¹ = 1/cube side, e.g. 2.0);
  `Grid_perfil` (molecule max-dimension as % of grid, e.g. 80); `Surface` = `ses`|`skin`|
  `blobby`|`mesh`; `Skin Surface Parameter` s∈(0,1), default 0.45 (1=convex hull, 0=VdW);
  `Blobbyness` default −2.5; `Accurate Triangulation` (analytic vertex sampling on/off);
  `Triangulation`; `Smooth Mesh` (Laplacian); `Compute Vertex Normals` (analytic where
  possible, else averaged face normals — "the mesh will always have all the normals defined");
  `Vertex Atom Info` (nearest atom index per vertex); `Save Mesh MSMS Format`; cavity options
  (`Cavity Detection Filling`, `Conditional Volume Filling Value` in Å³, `Keep Water Shaped
  Cavities`); pocket options (`Pocket Radius Big/Small`, default 3.0/1.4 Å);
  `Max Probes Self Intersections` (default 100) and `Self Intersections Grid Coefficient` for
  the SES self-intersection grid; acceleration-grid settings (`Max ses patches auxiliary grid
  size` etc.).
- **Note on "SWM4 / Del Bino".** The research brief mentioned "NanoShaper (Del Bino/Oliviero?)
  SWM4 / grid-based SES, collision detection". I could **not** verify a "SWM4" SES algorithm
  or authors "Del Bino/Oliviero" in connection with NanoShaper. (SWM4-NDP is a water model,
  unrelated.) What NanoShaper does have for SES is: alpha-shape + analytic patches + a
  *self-intersection grid* for collision detection during Connolly surface computation
  (`Max Probes Self Intersections`, `Self Intersections Grid Coefficient`). Treat "SWM4" as
  **unverified/likely a mix-up**. [UNCERTAIN]
- **Related recent work.** NanoShaperWeb (JCIM 2025, doi:10.1021/acs.jcim.5c00821,
  https://pubs.acs.org/doi/10.1021/acs.jcim.5c00821) wraps NanoShaper (VdW/SES/skin/blobby)
  in a web UI with pocket detection; NanoShaper-VMD interface (Bioinformatics 2019, 35(7):
  1241–1243, https://pmc.ncbi.nlm.nih.gov/articles/PMC6449750/).

### 1.3 EDTSurf (Xu & Zhang) — Euclidean distance transform + MC/VCMC

- **What it is.** Open-source program constructing triangulated VdW, SAS, and molecular
  (SES) surfaces via **fast Euclidean Distance Transform (EDT)** on a voxel grid; also finds
  cavities and computes atom/residue depth. Authors Dong Xu & Yang Zhang (Zhang lab).
  Sources: https://zhanggroup.org/EDTSurf/ ; Xu D. & Zhang Y. (2009) "Generating Triangulated
  Macromolecular Surfaces by Euclidean Distance Transform", PLoS ONE 4(12):e8140 —
  https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0008140
- **Algorithm sketch (from the official page + paper abstract + NanoShaper paper's account).**
  Voxelize; EDT to the SAS solid and to the VdW solid; the molecular solid = voxels within
  probe radius of the SAS solid boundary ("outer green part is voxels which have Euclidean
  distance less than probe radius"); triangulate the resulting scalar field with **Marching
  Cubes (MC)** or **Vertex-Connected Marching Cubes (VCMC, default)**; optional Laplacian
  smoothing (vertex moved to centroid of one-ring). Cavities found from the grid; cavity count
  depends on resolution and probe radius (their figure: box 128 / probe 1.2 → 8 cavities;
  box 256 / probe 1.4 → 8 cavities).
- **CLI (from https://zhanggroup.org/EDTSurf/).**
  `-i input`; `-o prefix`; `-t triangulation type` 1=MC, 2=VCMC (default 2);
  `-s surface type` 1=VWS, 2=SAS, 3=MS(SES), 0=DEPTH (default 3);
  `-p probe radius` in [0, 2.0], default **1.4**;
  `-f scale factor` in (0, 20.0], default **4.0** — molecule is scaled to fit a bounding box;
  "larger is better but increases memory"; auto-resets if the box would overflow.
  `-c color mode` 1=pure, 2=atom, 3=chain; `-h` inner/outer surface selection (1=both,
  2=outer, 3=inner).
- **Output.** `outname.ply` (triangulated, optionally colored), `outname-cav.pdb`,
  `outname.asa` (per-residue ASA), depth files. PLY is the native export.
- **Performance/quality claims.**
  - NanoShaper paper: EDTSurf "never failed" in their tests and is "in general, a very fast
    algorithm"; time scales with **grid size**, not atoms (they observed grids ~300³ making
    runtime almost constant across molecules). EDTSurf's cavity count oscillates at low scale;
    needs high scale (~8 grids/Å in their tests) to agree with NanoShaper cavity counts.
  - NanoShaper reaches EDTSurf's accuracy at ~half the scale; EDTSurf vertices do not
    necessarily lie exactly on the surface (VCMC + post-smoothing), so its error is partly
    off-surface (Decherchi & Rocchia 2013).
  - Hermosilla et al. 2017 compared their GPU method to EDTSurf (CPU, website build): their
    256³ GPU grid is 2–25× faster than EDTSurf's non-power-of-two grids on 2k–83k-atom cases
    (see §2.4 table).
- **Limitations.** Max scale limited by bounding box (NanoShaper: EDTSurf could not exceed
  ~8 grids/Å on calmodulin 3CLN in their tests). Fullerene symmetric case "has problems" per
  NanoShaper paper. Source code mirror: https://github.com/UnixJunkie/EDTSurf (GPL-style
  permissive notice on the Zhang lab page).

### 1.4 FreeSASA vs ESRF `freesas` — disambiguation (important)

- **FreeSASA (mittinatten).** Open-source **C library + CLI + Python** computing **SASA areas
  only** (no mesh). Implements **Lee & Richards'** (slice/arc) method (default) and
  **Shrake–Rupley** spherical sampling; configurable radii sets (default, ProtOr, NACCESS,
  OONS...), accuracy (number of test points / slice thickness), and output granularity.
  Citation: Mitternacht, S. (2016) "FreeSASA: An open source C library for solvent accessible
  surface area calculations", F1000Research 5:189, doi:10.12688/f1000research.7931.1
  (arXiv:1601.06764 — https://arxiv.org/abs/1601.06764). Repo: https://github.com/mittinatten/freesasa ;
  docs http://freesasa.github.io/. **No SES, no triangulation** — only areas (and per-atom/
  per-residue breakdowns).
- **`freesas` (ESRF / silx-kit).** A Python package for **small-angle scattering (SAXS/SANS)
  data analysis** (autorg, bift, supcomb, cormap...) — *not* a surface code despite the name.
  https://github.com/silx-kit/freesas ; https://pypi.org/project/freesas/ ; described in
  J. Synchrotron Rad. (2022) 29, 1318–1328. **Do not confuse** with FreeSASA.

### 1.5 Library implementations in visualization / electrostatics stacks

- **VMD `surf` (Varshney).** VMD's "Surf" representation "uses the molecular surface solver
  written by Amitabh Varshney when he was at the University of North Carolina": VMD writes
  radii+coordinates to a temp file, runs the `surf` executable with the probe radius, reads
  the result back. Source: VMD User's Guide —
  https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.4/ug/node74.html (also node65 of 1.7.1).
  Algorithm lineage: Varshney, Brooks & Wright, "Linearly Scalable Computation of Smooth
  Molecular Surfaces" / "Computing smooth molecular surfaces", IEEE CG&A 1994 —
  https://ieeexplore.ieee.org/document/310720 (analytic, parallelizable, expected linear
  complexity in atom count; Richards/Connolly smooth molecular surface).
  VMD also wraps **MSMS** and, in recent versions, **NanoShaper** (VMD 1.9.4 UG node76:
  surface type SES/Skin/Blobby, probe radius option, representation method,
  https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.4/ug/node76.html) and has a GPU
  Gaussian/"QuickSurf" representation (Krone et al., EuroVis 2012; not SES).
- **Chimera / ChimeraX.**
  - Chimera historically used MSMS. Tom Goddard's 2013 design review ("Molecular Surface
    Algorithms", https://www.cgl.ucsf.edu/chimera/data/surface-oct2013/surface.html) compares
    MSMS vs distance-grid SES vs Gaussian vs LSMS vs PyMOL, and describes the **distance-grid
    SES** algorithm that was chosen as the MSMS replacement (~150 lines of new code on top of
    Chimera's contouring). Performance: ribosome 1jj2 (100k atoms) at 0.5 Å grid →
    463×463×371 grid (~80M points, 320 MB as f32) in **~30 s unoptimized**, where MSMS fails.
    Meshes at 0.5 Å have ~3× more triangles than MSMS density=2 (3fx2: 22,546 vs 66,828
    triangles) with "many skinny triangles".
  - ChimeraX exposes it as an API: `ses_surface_geometry(xyz, radii, probe_radius=1.4,
    grid_spacing=0.5, sas=False)` returning vertex/normal/triangle arrays (float32 xyz,
    int32 triangles, float32 unit normals) —
    https://chimerax.readthedocs.io/en/latest/modules/surface/surface.html. Implementation:
    `src/bundles/surface/src/gridsurf.py` in https://github.com/RBVI/ChimeraX (LGPL-2.1):
    f32 distance grid padded by 2·probe + rmax + spacing; `sphere_surface_distance` of
    radii+probe; contour at level 0 with `calculate_normals=True` → SAS; then probe spheres of
    radius probe_radius centered at SAS vertices; re-contour → SES; delete connected
    components whose sample vertex is >1.5·probe_radius away from all atom spheres.
    `molsurf.py` (atomic bundle) calls it with `self.probe_radius, self.grid_spacing`.
  - ChimeraX `surface` command `sharpBoundaries` option (default true for SES) "adjusts the
    triangulation of the surface to produce abrupt and relatively straight boundaries between
    per-atom patches" —
    https://opentutorials.cgl.ucsf.edu/docs/user/commands/surface.html
  - LSMS (Can, Chen & Wang, J Mol Graph Model 2006, level-set SES) evaluated by Goddard:
    poor lighting normals (staircase artifacts at all grid sizes), 32-bit demo >4 GB on medium
    proteins, only power-of-2 cube grids up to 512³; at comparable quality estimated 5–40×
    slower than MSMS (his re-analysis: ~7–10× slower at matched quality).
- **PyMOL.** `surface` command; settings `surface_quality`, `surface_solvent`,
  `surface_carve`, `surface_circumscribe`, `surface_ramp_max_mode`; "gaussian resolution"
  smoothing knob. Goddard (2013) inspected PyMOL 1.6 source and concluded it looks like a
  Gaussian/sum-of-balls style surface ("distance maps … carving"), but could not fully
  decipher the 500-line uncommented C routines — **exact algorithm not verified**.
  Sources: https://wiki.pymol.org/index.php?title=Surface (mirror currently erroring;
  cached snippets), https://pymol.org/dokuwiki/doku.php?id=setting:surface,
  Goddard 2013 (URL above). [PyMOL algorithm: UNCERTAIN]
- **PDB2PQR / APBS.** PDB2PQR prepares PQR (charge+radius) files for PB electrostatics
  (Dolinsky et al., NAR 2004, doi:10.1093/nar/gkh381) — it is a *preparation* tool, not a
  surface mesher. APBS consumes surfaces/meshes for FEM/FD; BEM-oriented meshing pipelines
  use e.g. MSMS/NanoShaper/TMSmesh output. See the remeshing literature in §3.
- **MolSurfComp (Quan et al.).** MATLAB package computing the SES with a complete *a priori*
  characterization of SES singularities via a Voronoi-type diagram of balls —
  https://github.com/quanchaoyu/MolSurfComp ; theory: Quan et al., "Mathematical analysis and
  calculation of molecular surfaces", J. Comput. Phys. 319 (2016) 1–23,
  doi:10.1016/j.jcp.2016.07.007. This is the theoretical basis cited by the SESDF GPU work (§1.6).
- **TMSmesh / TMSmesh 2.0.** Manifold triangular meshing of **Gaussian** molecular surfaces
  for >1M-atom biomolecules (Liu, Chen & Lu); 2.0 is >30× faster than 1.0; meshes commonly
  need post-remeshing (see §3). Cited via Khan et al. 2018 (§3) and Liu et al.,
  J. Mol. Model. 2015, 21:113, doi:10.1007/s00894-015-2654-9.

### 1.6 GPU implementations of SES/SAS

- **GPU Contour-Buildup (Krone et al.).** Parallel contour-buildup (CB) for the SES patches;
  CUDA implementation in MegaMol (https://github.com/UniStuttgart-VISUS/megamol, plugin
  `protein_cuda`). Krone, Grottel & Ertl, "Parallel contour-buildup algorithm for the
  molecular surface", BioVis 2011, doi:10.1109/BioVis.2011.6094043; also Krone et al.
  TVCG 2009 (interactive surface dynamics). Characteristics: among the fastest analytic
  options for ≤~10k atoms (Hermosilla et al.), but **large GPU memory**: 1S3S (22k atoms)
  needs 1.9 GB (32 ms on GTX 970); 3K19 (~30k atoms) 2.7 GB — so large virus capsids are
  impractical (Hermosilla et al. 2017, Results).
- **Hermosilla et al. 2017/2018 — grid-based GPU SES with progressive refinement.**
  "Interactive GPU-based Generation of Solvent Excluded Surfaces", Computers & Graphics
  (2018), doi:10.1007/s00371-017-1397-2; preprint:
  https://www.uni-ulm.de/fileadmin/website_uni_ulm/iui.inst.100/institut/Papers/viscom/2017/hermosilla17ses.pdf
  - Algorithm: regular 3D grid of distance values (mipmapped 3D float array, base 512³ max).
    Step 1 classifies grid points (outside / inside SES / on boundary) by probe–atom
    intersection tests using a **spatial hash (Green's method)** with cell size =
    probe radius + max vdW radius → 27-cell neighborhood queries. Step 2 refines boundary
    points' distances. Work is grouped in **8³ bricks**; only bricks containing boundary
    points (plus interior and adjacent bricks to avoid missing cavities/coherence errors)
    are refined. Progressive LOD: coarse level computed in milliseconds (interactive),
    refinement in background (512 bricks/frame). Missing-feature detection re-samples suspect
    bricks at 64³. Rendering by **ray-marching the distance field** (no mesh needed);
    normals via **Sobel filter (26 texture lookups)** or central differences.
  - Performance (Table 1, GTX 970 + i7 + 12 GB): e.g. 2,066 atoms: sort 0.16 ms, base 64³
    5.42 ms, full 256³ refinement 775 ms; 22,367 atoms (1S3S): base 128³ 7.58 ms, full 512³
    930 ms; 545k atoms (1K4R): base 15.89 ms, full 512³ 1,649 ms. Memory: **620 MB** for the
    512³ field + 32 B/atom — vs 1.9–2.7 GB for GPU Contour-Buildup.
  - Accuracy vs analytic SES (Chimera/MSMS, 3,967 atoms, 256³): mean surface deviation
    0.23 Å, RMS 0.27 Å — small outward bias of order the sample spacing.
  - Vs EDTSurf (Table 3): consistently faster (e.g. Traj3 11k atoms: EDTSurf 5.19 s vs 0.20 s
    at 256³; 3EXG 83k atoms: 3.64 s vs 0.16 s).
  - Implementation re-published as **QuickSES** (CUDA): Martinez, Krone & Baaden, "QuickSES:
    A Library for Fast Computation of Solvent Excluded Surfaces", Eurographics 2019,
    doi via HAL hal-02370900 (https://hal.archives-ouvertes.fr/hal-02370900/document);
    repo https://github.com/nezix/QuickSES (MIT). CLI: `-i in.pdb -o out.obj -v voxel_size`
    (default **0.5 Å**, "defines the quality of the mesh"), `-l` Laplacian smoothing passes
    (default 1), `-s slice size` (default 300 sub-grid = GPU memory knob). 3D uniform grid for
    atom neighbors; CUDA Marching Cubes; GPU vertex welding. Output OBJ.
- **SESDF (Plateau-Holleville et al. 2024) — analytic SES on GPU, complete surface.**
  "Efficient GPU computation of large protein Solvent-Excluded Surface", IEEE TVCG (2024),
  doi:10.1109/TVCG.2024.3380100 — https://ieeexplore.ieee.org/document/10477549 ;
  code https://github.com/PlathC/SESDF (reference implementation, C++17 + CUDA ≥11.8,
  benchmark targets `BCSBenchmarkCBSESDF` etc.); project page https://plathc.github.io/publications/sesdf/
  - Based on Quan et al.'s implicit SES / Voronoi-of-balls characterization; builds the
    analytic patches (spherical/toroidal) with precise preallocation of patch counts to cut
    memory; supports **complete surface** (including internal/cavity patches) as well as an
    exterior-only simplified mode.
  - Performance (Table 2 of the paper; AMD Ryzen 5 1600 + RTX 2080; probe 1.4 Å; 1000
    iterations mean): vs MegaMol GPU Contour-Buildup — similar or slightly faster times with
    far less memory. Examples: 1AGA 126 atoms: CB 4.50 ms / 11.92 MB vs SESDF 3.22 ms /
    0.39 MB (exterior 2.21 ms / 0.15 MB). 65,008 atoms (7RGD): CB 16.86 ms / 6.12 GB vs
    SESDF 18.06 ms / 214 MB (exterior 14.64 ms / 105 MB). Largest case **3J3Q, 2,440,800
    atoms: 700 ms / 7.6 GB** complete (324 ms / 3.7 GB exterior). CB cannot run past ~65k
    atoms in their setup (memory). Claim on the project page: "On the same GPU, our method is
    able to compute the surface of 30× larger protein without being restricted to the exterior
    surface." Memory grows linearly with atoms; exterior mode ~2× less memory than complete;
    ~30–60× less than CB on the largest comparable complex.
  - Overheads (allocation + CPU–GPU transfers for precise patch counts) fall from ~60% of
    total time (smallest molecules) to ~25% (largest).
- **GPU marching cubes / blobby.** Kim, Kim & Seong, "GPU accelerated molecular surface
  computing", Appl. Math. Inf. Sci. 6 (2012) 185S–194S (grid SES similar to EDT family);
  D'Agostino et al., "CUDA accelerated blobby molecular surface generation" (PPAM 2012);
  Krone, Stone, Ertl & Schulten, "Fast visualization of Gaussian density surfaces…" EuroVis
  2012 (GPU Gaussian). Cited via Decherchi & Rocchia 2013 refs [20], [29], [30].
- **Rendering-side GPU ray casting of analytic patches** (not meshing): Krone et al. TVCG 2009 /
  MegaMol "RC" path; Hermosilla et al. compare ray-marching vs that ray casting (ray marching
  wins on >100k atoms: e.g. 1K4R 545k — 144 fps at 512³ ray-marching vs 12.6 fps ray casting).

### 1.7 Alpha-shape / Voronoi / QuickHull-family surface tools

- **Alpha shapes** (Edelsbrunner): the classical analytic route to SES/SAS topology and to
  area/volume of unions of balls — Liang, Edelsbrunner, Fu, Sudhakar & Subramaniam (1998)
  "Analytical shape computation of macromolecules: Molecular area and volume through alpha
  shape", Proteins 33:1–17. NanoShaper's SES build is the modern production use (CGAL
  weighted Delaunay + alpha shapes, §1.2). Review of geometric measures of large biomolecules
  (surface, volume, pockets) via alpha shapes: https://pmc.ncbi.nlm.nih.gov/articles/PMC3188685/
  (Edelsbrunner & colleagues' methods).
- **Voronoi/power diagrams of atoms.** Voronota (Olechnovič) — Voronoi-diagram-of-atomic-balls
  analysis toolkit, https://github.com/kliment-olechnovic/voronota . Power-SASA (Rust, §6)
  computes SASA via 3D weighted power (Laguerre) diagrams. Quan et al. 2016 (MolSurfComp)
  uses a Voronoi-type diagram of balls to enumerate SES singularities.
- **Skin surface** (Edelsbrunner): mixed complex of weighted balls; NanoShaper implements it
  (shrinking parameter s); CGAL's Kruithof algorithm is the reference (slow: 4000 atoms →
  100+ s and GBs vs NanoShaper seconds and MBs per Decherchi & Rocchia 2013).
- **QuickHull** is not the basis of any mainstream SES tool found in this survey; convex hull
  appears as the s→1 limit of the Skin surface (NanoShaper guide) and in probe-position hull
  constructions. [The brief's "QuickHull-based surface tools" maps best onto the
  alpha-shape/weighted-Delaunay family above.] [UNCERTAIN / probably conflated]

---

## 2. Grid-based SES algorithm details used in practice

### 2.1 Building the SES scalar field

Three recurring constructions:

1. **Distance-to-union-of-spheres field, thresholded, twice (Chimera/ChimeraX "distance grid").**
   Goddard's description (2013, URL in §1.5) — the cleanest public spec:
   - Build a signed-distance grid to the union of atom spheres **inflated by the probe
     radius**: per sphere, within a slightly-larger-than-sphere subregion, set grid values to
     distance from that sphere's surface (positive outside), taking the **min** over spheres.
   - Contour at level 0 → smooth **SAS**.
   - Place probe spheres (radius = probe) at every **vertex of the SAS contour** (not at grid
     points); rebuild a distance grid to those probe spheres; contour at 0 → SES (plus a
     spurious outer component one probe-diameter out; delete components ≥1.5·probe_radius
     from atom surfaces, tested on one vertex per connected component).
   - Key quality insight: use **floating-point distances, not 0/1 masks** (masks give
     stair-steps); centering probes at contour vertices rather than grid points is what makes
     the result smooth.
   - ChimeraX code implements exactly this (`gridsurf.py`, §1.5) with `contour_surface(...,
     calculate_normals=True)`.

2. **Classification + boundary distance refinement (Hermosilla/QuickSES).** Per grid point:
   outside SES / inside SES / on boundary, by testing probe-at-grid-point vs atoms (collision
   test) with spatial hashing; boundary points then get distances refined from nearest outside
   neighbor. Field is initialized with the probe radius; mipmapped for LOD (§2.3).

3. **EDT-based solids (EDTSurf).** Voxel indicator of the inflated-sphere union (SAS solid);
   Euclidean distance transform of the complement identifies the shell within probe radius;
   the "molecular solid" is carved; MC/VCMC extracts the surface. The grid acts both as the
   distance field and as the cavity detector.

**NanoShaper** is *not* grid-based for SES shape: patches are analytic; the grid only hosts
ray-cast in/out coloring, MC topology, and PDE maps. Its "Analytical Intersections Marching
Cubes" places vertices at analytic ray/surface intersections rather than by linear
interpolation — a third way to get grid topology with exact vertex placement.

### 2.2 Mesh extraction and normals

- **Classic Marching Cubes** (Lorensen & Cline, SIGGRAPH 1987) is the default extractor for
  grid methods (ChimeraX contour_surface, EDTSurf MC mode, QuickSES CUDA MC, NanoShaper MC
  tables). Known issue: ambiguous cases / topological cracks; table variants (Nielson–Hamann,
  Chernyaev) used in practice. Docs of the Rust `isosurface` crate summarize MC tradeoffs:
  "Pretty fast. Extraction has no dependencies on neighbouring chunks. Cons: produces a lot of
  small triangle slivers; can't accurately reproduce sharp corners; cracks between chunks of
  differing levels of detail." https://docs.rs/isosurface/latest/isosurface/marching_cubes/
- **VCMC (Vertex-Connected Marching Cubes)** — EDTSurf default (`-t 2`). Vertices shared
  across cubes ("vertex connected"), giving better connectivity than naive MC; Decherchi &
  Rocchia note VCMC+smoothing moves vertices off the true surface.
- **Analytical-Intersections MC** (NanoShaper) — MC topology, ray-cast exact vertices.
- **No mesh at all**: Hermosilla et al. ray-march the distance field directly for rendering
  (normals from Sobel filter over the field). Relevant if the product ever renders without
  exporting.

**Normal computation options observed:**
- Gradient of the scalar field at the vertex (central differences; Sobel 3D for smoother
  normals — Hermosilla). Works whenever the field is a distance field.
- Analytic normals from the surface definition: MSMS per-vertex normals from the analytic
  patch geometry (with arbitrary direction at singular vertices — documented bug);
  NanoShaper `Compute Vertex Normals` analytic where possible, falling back to averaged
  triangle normals.
- Angle-weighted vertex normals (pseudonormals): the standard robust choice for triangle
  meshes — Bærentzen & Aanæs (2005), "Signed Distance Computation Using the Angle Weighted
  Pseudonormal", IEEE TVCG 11(3):243–253, doi:10.1109/TVCG.2005.49. Also used implicitly by
  "average the normals of surrounding triangles" fallbacks.
- Face normals only are insufficient for shading (faceted look); Goddard documents the
  subtlety that **true SES is C¹ but not C²** at sphere/toroid joins (curvature jumps from
  1/R_atom to −1/R_probe), so even exact normals give visible lighting "seams"; smoother
  shading than analytic SES actually requires C² surfaces (Gaussian/Skin) or smoothing.

### 2.3 Resolution / LOD controls (parameter names by tool)

| Tool | Resolution knob | Default / typical | Notes |
|---|---|---|---|
| MSMS | `-density` (vertices/Å²) | 1.0 (large), 3.0 (small) | triangulation density, analytic surface |
| NanoShaper | `Grid_scale` (Å⁻¹) + `Grid_perfil` (%) | 2.0 / 80 | scale = 1/cube side; 2 grids/Å = 0.5 Å |
| EDTSurf | `-f scale factor` | 4.0 (0–20] | box enlargement; bounded by max box |
| ChimeraX | `grid_spacing` (Å) | 0.5 | distance-grid contouring |
| QuickSES | `-v voxel size` (Å) | 0.5 | "defines the quality of the mesh" |
| Hermosilla | mipmapped grid levels, base ≤512³ | power-of-two levels | progressive LOD |
| Goddard guidance | 0.5 Å ≈ MSMS quality; 0.2 Å more detail | — | smoother appearance than MSMS dens=2 |
| NanoShaper pocket mode | `Pocket Radius Big/Small` | 3.0 / 1.4 Å | probe-radius pair, not mesh LOD |

LOD techniques in practice: (a) grid spacing; (b) power-of-two mipmaps with progressive
refinement and smooth LOD transitions (Hermosilla); (c) decimation/smoothing post-passes
(Laplacian: NanoShaper `Smooth Mesh`, EDTSurf default post-Laplacian, QuickSES `-l`); (d) MSMS
density as a pure triangulation knob on an analytic surface. Adaptive octrees appear as a
NanoShaper development note (memory 7.3 GB → 4 GB on 1VSZ) and in the Rust `isosurface` crate
(`linear_hashed_marching_cubes`), but I found **no production SES tool shipping adaptive
octree meshing** as the default. [octrees: limited evidence]

### 2.4 Reported performance numbers (collected)

| System | Hardware | Case | Time | Memory | Source |
|---|---|---|---|---|---|
| NanoShaper SES+Skin | mid workstation | 1VSZ 180k atoms, 2 grids/Å | 234 s (Skin: 17+88+129 s) | 7.3 GB (4 GB octree dev) | Decherchi & Rocchia 2013 |
| NanoShaper | dual 8-core SandyBridge | SV40 (1sva) @2.0 Å⁻¹ | ~2 min incl. normals | — | NanoShaper 0.7 guide |
| Chimera distance grid | single core, unoptimized | 1jj2 100k atoms, 0.5 Å | ~30 s | 320 MB f32 grid (463×463×371) | Goddard 2013 |
| QuickSES/Hermosilla | GTX 970, i7 | 545k atoms (1K4R) | 1.65 s to 512³ (~0.98 Å cell) | 620 MB + 32 B/atom | Hermosilla et al. 2017 Table 1 |
| QuickSES/Hermosilla | GTX 970 | 22k (1S3S) | 0.93 s to 512³ (0.40 Å) | same | same |
| QuickSES base LOD | GTX 970 | any (≤545k) | 5–16 ms | — | same (base level) |
| SESDF complete | RTX 2080 | 3J3Q 2.44M atoms | 700 ms | 7.6 GB | Plateau-Holleville et al. 2024 Table 2 |
| SESDF complete | RTX 2080 | 7RGD 65k atoms | 18 ms | 214 MB | same |
| GPU CB (MegaMol) | GTX 970 / RTX 2080 | 1S3S 22k | 32 ms | 1.9 GB | Hermosilla 2017 / SESDF 2024 |
| EDTSurf CPU | — | 2k–83k atoms | 1.6–5.2 s (grid ~200–300³) | grid-bound | Hermosilla Table 3 |
| MSMS | — | >10k atoms | often fails | — | Goddard 2013 |

Derived rules of thumb from these data:
- Grid memory ≈ 4 bytes × Nx·Ny·Nz; N ≈ (molecule extent + 2·pad)/spacing. At 0.5 Å a
  100k-atom ribosome is ~80M voxels = 320 MB; at 0.2 Å it is ~15× more.
- Analytic GPU (SESDF) scales to millions of atoms in <1 s with memory linear in atoms
  (~3 KB/atom complete, ~1.5 KB/atom exterior at their largest case).
- CPU analytic (MSMS) is fast per-component but fragile; CPU grid (EDTSurf) is robust but
  grid-size-bound; CPU ray-cast analytic (NanoShaper) is robust and reasonably fast,
  multithreaded.

---

## 3. Mesh quality requirements for molecular surfaces

- **Consumers drive the requirements.** Visualization needs smooth shading (good vertex
  normals), reasonable triangle count, and per-atom/per-residue coloring seams. BEM/FEM
  implicit-solvent solvers (AFMPB, mFES, DelPhi-BEM) need **manifold, watertight, preferably
  high-quality triangles**: "Molecular surface mesh generation is a prerequisite for using the
  boundary element method (BEM) and finite element method (FEM)… Molecular surface meshes
  typically have small angles, redundant vertices, and low-quality elements." — Khan, Yan,
  Gui, Lu & Zhang (2018), "Molecular Surface Remeshing with Local Region Refinement",
  IJMS 19(5):1383, doi:10.3390/ijms19051383, https://pmc.ncbi.nlm.nih.gov/articles/PMC5983798/
- **Defects observed on real molecular meshes (Khan et al. 2018 and refs therein):**
  zero/near-zero degree angles, self-intersecting triangle pairs, non-manifold edges/vertices,
  holes, redundant vertices, extreme aspect ratios. Their quality targets: minimal angle ≥30°,
  maximal angle ≤120°, good aspect ratio, regular-vertex ratio, and **preservation of area and
  volume** (they measure volume/area drift and solvation energy vs SMOPT/ISO2mesh).
  Laplacian and Taubin smoothing are the baseline fixes; SMOPT removes redundant vertices and
  self-intersections; MeshFix/PolyMender repair holes and non-manifoldness; cut-and-fill +
  RAR remeshing is their pipeline.
- **Watertightness / topology.** NanoShaper's whole ray-casting framework assumes a closed
  manifold surface and enforces the even-intersection checksum precisely so that grid
  inside/outside is consistent (needed for PDE maps). MC over a well-defined field yields
  watertight meshes if ambiguous cases are resolved consistently; raw MSMS meshes can have
  holes when triangulation fails ("the resulting triangulated surface will not be correct and
  may have holes" — MSMS docs `-no_rest_on_pbr`).
- **Consistent orientation.** Outward orientation is assumed by area/volume integrators;
  MSMS `.vert` normals are outward; ChimeraX returns unit normals from field contouring;
  any exporter should guarantee CCW-outward (or document the convention — see §5).
- **Smooth normals for shading.** See §2.2: gradient/Sobel of distance field, analytic patch
  normals, or angle-weighted pseudonormals (Bærentzen & Aanæs 2005). For true SES, C¹-only
  joins produce visible shading seams (Goddard 2013); optional Laplacian/Taubin smoothing
  trades geometric fidelity for visual smoothness (EDTSurf and NanoShaper both default-ish
  to offering smoothing; NanoShaper recommends it when `Accurate Triangulation` is off or
  BEM quality is needed).
- **Triangle aspect ratio / slivers.** MC "produces a lot of small triangle slivers"
  (Rust `isosurface` docs); Goddard measured 0.5 Å grid meshes at ~3× the triangle count of
  MSMS density=2 with many skinny triangles. Remeshing (RAR, CVT, edge split/collapse/flip)
  is the standard cure when the mesh feeds FEM/BEM.
- **Multiple chains / components.** MSMS writes one `.vert/.face` pair per surface component
  (`prefix_3.*` for component 3) and `-all_components` enumerates cavities; default is the
  external component only. ChimeraX/Chimera distance-grid deletes extra components from the
  double-contouring artifact (and with them, unwanted internal surfaces) via the
  1.5·probe distance test — note this also drops large internal pockets unless handled
  per connected component (Goddard tests "one vertex of each connected surface component").
  NanoShaper detects cavities by flood fill and can fill them below a volume threshold
  (`Conditional Volume Filling Value`, Å³), keeping water-sized or water-shaped cavities
  optionally.
- **Per-atom vs whole-molecule surfaces.** MSMS records `closest_sphere_index` per vertex
  (and per-face analytic-face ids) — this is what VMD/Chimera use for per-atom/per-residue
  coloring. NanoShaper offers `Vertex Atom Info` (nearest atom index per vertex) in OFF+N+A
  or MSMS format. Grid methods give jagged per-atom patch boundaries ("sawtooth edges" —
  Goddard), fixable only by triangle splitting at atom boundaries; Hermosilla colors by
  nearest atom via a coarse 3D texture for smooth property mapping.

---

## 4. Concurrency / vectorization patterns in this domain

- **Neighbor search: spatial hashing / cell lists (universal).** Subdivide space into a
  uniform grid of cells sized ≈ (max interaction radius); sort atoms by cell; query 27
  neighboring cells (3×3×3) for probe/atom collision tests. Used by Hermosilla et al. (cell
  size = probe radius + max vdW radius; method of Green — "Particle simulation using CUDA",
  NVIDIA, cited as [6] in Hermosilla et al. 2017), QuickSES, SESDF (grid of atom cells),
  RustSASA ("highly optimized spatial grid"), and standard MD cell lists generally.
  Complexity O(N) for fixed radius.
- **Domain decomposition / bricks for grids.** Hermosilla et al.: 8³ bricks, workgroup size
  = brick size for GPU coherence; only active bricks (boundary + interior + adjacent) are
  refined at higher LODs; brick lists are downloaded to CPU and re-uploaded per level —
  a classic GPU-friendly work-queue pattern that avoids CPU-side compaction of boundary
  points. QuickSES `-s slice size` (default 300) slices the domain to cap GPU memory.
- **Ray-level parallelism.** NanoShaper: rays are independent → multithreaded over CPU cores
  (guide: 4 threads/physical core recommended); 2D acceleration grids over patch bounding
  boxes avoid per-ray full scans (`Max ses patches auxiliary grid 2d size` etc.).
- **Component-level parallelism.** MSMS `-all_components` parallel with OpenMP since 2.6.2.
- **SIMD.** RustSASA uses **Rayon** (data parallelism) + **Pulp** (runtime-detected SIMD) for
  Shrake–Rupley point tests, giving 5× vs FreeSASA multi-threaded over a proteome and ~43%
  single-threaded (Campbell 2026, §6). Distance/point-in-sphere inner loops are the obvious
  SIMD target in any grid field builder (SoA xyz+r, FMA distance² compare). GPU methods
  instead get throughput from thread-per-voxel / thread-per-patch kernels (Hermosilla's two
  compute shaders; SESDF's staged kernels).
- **GPU vs CPU tradeoffs (as reported).**
  - GPU analytic (SESDF/QuickSES-class): sub-second for 10⁵–10⁶ atoms; memory linear in
    atoms (analytic) or ~620 MB fixed (512³ field); consumer GPUs sufficient (RTX 2080 / GTX 970).
  - GPU CB: fast for ≤~10⁴–10⁵ atoms but memory-hostile (1.9 GB @22k atoms).
  - CPU analytic (MSMS): fast per patch-set but serial fragility (restarts) and >10k-atom
    crashes.
  - CPU grid (Chimera distance grid, EDTSurf): robust, simple, scales with grid volume;
    30 s for 100k atoms unoptimized single-thread is the reference point to beat.
  - CPU multithreaded analytic+raycast (NanoShaper): robustness leader (exact predicates,
    cavity handling), minutes for 180k atoms at 2 grids/Å.
- **Vertex welding / mesh post.** QuickSES implements GPU vertex welding (MC emits per-cube
  vertices; welding merges duplicates). Any parallel MC needs a deterministic or
  parallel-friendly welding step (spatial hash of vertex keys on edge IDs is the usual trick).
- **Precision.** NanoShaper: all double; CGAL exact predicates/inexact constructions.
  GPU paths: f32 typical (ChimeraX grids f32; Hermosilla f32 distance field). For BEM
  downstream, f64 vertices may matter; for rendering, f32 suffices.

---

## 5. File I/O conventions for molecular surfaces

- **Units.** Everything is in **Ångström**. PDB coordinates are Å; MSMS xyzr, EDTSurf, and
  QuickSES all operate directly in Å (QuickSES `-v` "voxel size in Angstrom", output `v %.3f`
  of Å coordinates). OBJ has no unit; the de facto convention for molecular meshes is Å,
  matching the input structure. Visualization tools (VMD, ChimeraX, PyMOL, MeshLab, Blender)
  consume them as unitless numbers consistent with the PDB.
- **OBJ (the format our product targets).**
  - QuickSES writes plain OBJ: `v x y z` per vertex and `f i j k` (1-based) per triangle —
    **no `vn` lines** in the current source (`CudaSurf.cu::writeToObj`). Face winding is
    written as `(triangle.y+1, triangle.x+1, triangle.z+1)`, i.e. a deliberate permutation
    of the stored order (likely to flip/normalize winding; a consumer should still verify
    outward orientation). Source: https://github.com/nezix/QuickSES (CudaSurf.cu).
  - Wavefront OBJ supports per-vertex normals via `vn` + face syntax `f v//vn ...` (or
    `v/vt/vn`). For shading, exporting `vn` per vertex (MSMS-style) is what viewers prefer;
    some pipelines (MathWorks `writeSurfaceMesh` docs) note OBJ writers that drop normals
    entirely — so being explicit matters.
  - Recommended conventions for our OBJ export (inferred from tool practice, not a formal
    standard): coordinates in Å, 1-based faces, CCW winding when viewed from outside (outward
    normals), optional `vn` per vertex, optional `g`/`o` groups per connected component or
    per chain.
- **Companion formats in the domain:**
  - **MSMS `.vert` / `.face`** — the lingua franca for per-vertex normals + closest-atom +
    patch-type metadata (§1.1). VMD/Chimera/PyMOL ecosystems all read it. NanoShaper can
    *write* MSMS format (`Save Mesh MSMS Format`, requires `Compute Vertex Normals`).
  - **OFF (Geomview)** — NanoShaper's native export (`triangulatedSurf.off`), plus custom
    **OFF+N** and **OFF+N+A** extensions storing vertex normals and nearest-atom index
    (NanoShaper 0.7 guide). Also importable (`msms` files auto-converted to `msms.off`).
  - **PLY** — EDTSurf native output (`outname.ply`, with per-vertex color modes pure/atom/
    chain); widely viewable in MeshLab. PLY supports vertex normals and colors natively
    (binary or ASCII).
  - **STL** — common for 3D printing; no normals per vertex, no color, no units; redundant
    vertices; rarely used in molecular workflows but requested by some makers.
  - **X3D / VRML / glTF / FBX** — appear in visualization pipelines (X3D in older web viewers;
    glTF in modern engines). MathWorks `writeSurfaceMesh` supports OBJ/GLTF/GLB/OFF/DAE/…
    with the caveat that OBJ/GLTF writers there only carry vertices+faces (normals computed
    by the viewer).
  - **PDB/PQR/XYZR** are the input side: XYZR (x y z r per line) is MSMS/NanoShaper's native
    input; PQR adds charge (PDB2PQR); mmCIF is the modern PDB replacement and should be a
    first-class input for our CLI.
- **What visualization tools expect.** VMD reads MSMS vert/face and can call surf/msms/
  NanoShaper; ChimeraX computes in-memory (vertices/normals/triangles numpy arrays) and its
  `save` supports PLY/OBJ-style scene export; MeshLab reads PLY/OFF/STL/OBJ with normals.
  For per-residue coloring without extra maps, carrying a per-vertex atom index (MSMS field,
  OFF+N+A) or per-vertex color (EDTSurf PLY) is the established practice.

---

## 6. Rust crates and recent (2020+) open-source implementations

### 6.1 Directly relevant Rust crates

- **RustSASA** (`rust-sasa` / github.com/maxall41/rustsasa). SASA via Shrake–Rupley in Rust
  with **Rayon + Pulp SIMD + optimized spatial grid**; CLI, Python bindings, MDAnalysis
  integration. JOSS paper: Campbell, M. J. (2026), "RustSASA: A Rust Crate for Accelerated
  Solvent Accessible Surface Area Calculations", JOSS 11(117):9537,
  doi:10.21105/joss.09537 — https://joss.theoj.org/papers/10.21105/joss.09537 ;
  repo https://github.com/maxall41/rustsasa . Claims: 5× vs FreeSASA and 63× vs Biopython
  on E. coli proteome multi-threaded (~5 s vs ~28 s vs ~328 s, 8 cores); single-protein
  4.0 ms (parity with FreeSASA); single-threaded proteome 26.0 s vs FreeSASA 46.2 s (~43%
  faster); ~20× vs mdakit-sasa on an MD trajectory. Accuracy: r ≈ 1.0 vs FreeSASA, chain
  RMSE ~44 Å². **Areas only, no meshes.** Reviews by @mittinatten (FreeSASA author) and
  @kliment-olechnovic (Voronota author) — good sign for correctness pedigree.
- **powersasa** (lib.rs/crates/powersasa, LGPL-3.0). SASA + molecular volume via **3D
  weighted power diagram (Laguerre/Voronoi)**; per-atom exposed spherical contours from
  diagram topology; contour integrals for area/volume. Provenance note: LLM-assisted
  translation of legacy C++ from the SIMONA package (original algorithm authors credited).
  https://lib.rs/crates/powersasa (crate page). Small adoption (used by `volar`); treat as
  reference-grade, check license (LGPL) before depending.
- **isosurface** (github.com/swiftcoder/isosurface, Apache-2.0, docs.rs/isosurface).
  Isosurface extraction algorithms with docs linking the academic papers: `marching_cubes`,
  `linear_hashed_marching_cubes` (MC over a linear hashed octree — adaptive LOD!), `point_cloud`,
  `source` traits. No dependencies in the core. MC caveats documented (slivers, no sharp
  corners, LOD cracks). Candidate to study/replace — its API shape (source trait + extractor)
  is a reasonable design model.
- **mcubes** (lib.rs/crates/mcubes, MIT). Simple Marching Cubes from volume data, Rust or
  Python bindings, ~500 LoC, outputs a Mesh/Vertex. Good minimal reference.
- **therealnv6/marching-cubes** (GitHub). Simple MC implementation, reference-grade.
- Supporting crates one would compose: `rayon`, `pulp`/`wide`/`simba` (SIMD), `nalgebra`/
  `glam`, `rstar` (R-tree) or a hand-rolled cell list, `pdbtbx`/`mmtf` (structure parsing),
  `obj`/`ply-rs`/`tobj` (export), `bytemuck`.

### 6.2 Other recent open-source SES/SAS worth knowing (any language)

- **QuickSES** (2019, MIT, CUDA) — §1.6. The most convenient open GPU SES mesh emitter (OBJ).
- **SESDF** (2024, reference impl, CUDA, C++17) — §1.6. The performance/memory frontier for
  complete analytic SES; benchmark harness included.
- **ChimeraX `surface` bundle** (LGPL-2.1, Python+C++ `gridsurf.py`) — §1.5. The best short,
  readable reference implementation of the two-pass distance-grid SES (≈100 lines of logic).
- **NanoShaper** (GPL, C++, CGAL) — §1.2. Source at gitlab.iit.it/SDecherchi/nanoshaper and
  mirror github.com/shi-scala/NanoShaper (v0.7).
- **EDTSurf** — C/C++ sources released (Zhang lab page; mirror github.com/UnixJunkie/EDTSurf).
- **FreeSASA** (LGPL-ish? — check COPYING; C) — area-only but the reference SASA C library.
- **MolSurfComp** (MATLAB) — SES singularities theory + implementation.
- **Voronota** (C++) — Voronoi-of-balls analysis (SASA, contacts), complementary.
- **MaSIF ecosystem** (Python/TensorFlow) consumes MSMS meshes — shows the downstream
  expectation of MSMS-quality triangulation with per-vertex normals.

### 6.3 Competitive positioning notes for our Rust product

- No Rust crate currently computes **SES meshes** (only SASA areas and generic isosurface
  extraction). The niche is open.
- Performance bar to be credible: CPU grid SES ≤30 s for 100k atoms at 0.5 Å (Goddard's
  unoptimized number) and multi-threaded + SIMD should beat that substantially; analytic
  quality bar is MSMS (exact patches) but with MSMS's robustness fixed.
- Differentiators available: mmCIF in, OBJ/PLY/OFF out with `vn` and per-atom index;
  deterministic LOD (grid spacing like ChimeraX's `grid_spacing`, 0.5 Å default); watertight
  manifold guarantee; cavity policy flags like NanoShaper; rayon+SIMD internals (RustSASA
  proves the speedup story is real in Rust).

---

## Open questions / uncertainties

1. "SWM4" and "Del Bino/Oliviero" in the brief could not be verified for NanoShaper. Treated
   as a probable mix-up with SWM4-NDP water model or another tool. [UNCERTAIN]
2. PyMOL's exact surface algorithm (Gaussian vs SES variants) not verified from source;
   Goddard's 2013 reading is inconclusive. [UNCERTAIN]
3. QuickHull-specific molecular surface tools not found; alpha-shape/weighted-Delaunay family
   covers that niche. [probably conflated in the brief]
4. Exact MSMS binary licensing terms ("Sannerlab MIT-like Source-Available License") should be
   read from https://ccsb.scripps.edu/msms/ before any redistribution talk.
5. Adaptive octree meshing for SES exists in research/development notes (NanoShaper dev,
   Rust `isosurface` hashed octree) but no production default found. [limited evidence]
6. Some performance numbers (NanoShaper Table 1 adjusted timings, CGAL Kruithof comparison)
   rely on figures embedded as images in the PLOS paper; text-level numbers cited here are
   from the extracted text (1VSZ totals, 234 s breakdown, etc.).

---

## Bibliography

1. Sanner, M. F., Olson, A. J., Spehner, J.-C. (1996). Reduced Surface: An Efficient Way to
   Compute Molecular Surfaces. *Biopolymers* 38:305–320.
   doi:10.1002/(SICI)1097-0282(199603)38:3<305::AID-BIP4>3.0.CO;2-Y
   — MSMS paper. https://ccsb.scripps.edu/msms/ ; docs:
   https://sannerlab.github.io/software/msms/documentation.html
2. Decherchi, S., Rocchia, W. (2013). A general and Robust Ray-Casting-Based Algorithm for
   Triangulating Surfaces at the Nanoscale. *PLoS ONE* 8(4):e59744.
   doi:10.1371/journal.pone.0059744 — NanoShaper.
   https://pmc.ncbi.nlm.nih.gov/articles/PMC3618509/
3. Decherchi, S., Rocchia, W. NanoShaper 0.7 User Guide (and tutorial chapter "Building and
   Analyzing Molecular Surfaces: A Tutorial on NanoShaper", Springer, 2015,
   doi:10.1007/978-3-319-12211-3_10). Code: https://gitlab.iit.it/SDecherchi/nanoshaper ;
   mirror https://github.com/shi-scala/NanoShaper
4. Xu, D., Zhang, Y. (2009). Generating Triangulated Macromolecular Surfaces by Euclidean
   Distance Transform. *PLoS ONE* 4(12):e8140. doi:10.1371/journal.pone.0008140
   — EDTSurf. https://zhanggroup.org/EDTSurf/
5. Mitternacht, S. (2016). FreeSASA: An open source C library for solvent accessible surface
   area calculations. *F1000Research* 5:189. doi:10.12688/f1000research.7931.1
   — https://arxiv.org/abs/1601.06764 ; https://github.com/mittinatten/freesasa
6. (ESRF) freesas — Small angle scattering tools, Python. https://github.com/silx-kit/freesas ;
   J. Synchrotron Rad. (2022) 29, 1318–1328. (Not a surface code.)
7. Hermosilla, P., Krone, M., Guallar, V., Vázquez, P.-P., Vinacua, À., Ropinski, T. (2018).
   Interactive GPU-based Generation of Solvent Excluded Surfaces. *Computers & Graphics*.
   doi:10.1007/s00371-017-1397-2 — preprint:
   https://www.uni-ulm.de/fileadmin/website_uni_ulm/iui.inst.100/institut/Papers/viscom/2017/hermosilla17ses.pdf
8. Martinez, X., Krone, M., Baaden, M. (2019). QuickSES: A Library for Fast Computation of
   Solvent Excluded Surfaces. Eurographics. https://hal.archives-ouvertes.fr/hal-02370900/document
   — code https://github.com/nezix/QuickSES
9. Plateau-Holleville, C., Maria, M., Mérillou, S., Montes, M. (2024). Efficient GPU
   Computation of Large Protein Solvent-Excluded Surface. *IEEE TVCG*.
   doi:10.1109/TVCG.2024.3380100 — https://ieeexplore.ieee.org/document/10477549 ;
   code https://github.com/PlathC/SESDF ; page https://plathc.github.io/publications/sesdf/
10. Krone, M., Grottel, S., Ertl, T. (2011). Parallel contour-buildup algorithm for the
    molecular surface. *IEEE BioVis 2011*, 17–22. doi:10.1109/BioVis.2011.6094043
    — GPU Contour-Buildup (MegaMol `protein_cuda`).
11. Quan, C., et al. (2016). Mathematical analysis and calculation of molecular surfaces.
    *J. Comput. Phys.* 319:1–23. doi:10.1016/j.jcp.2016.07.007 — MolSurfComp theory;
    https://github.com/quanchaoyu/MolSurfComp
12. Lorensen, W. E., Cline, H. E. (1987). Marching Cubes: A High Resolution 3D Surface
    Construction Algorithm. *SIGGRAPH '87*.
13. Bærentzen, J. A., Aanæs, H. (2005). Signed Distance Computation Using the Angle Weighted
    Pseudonormal. *IEEE TVCG* 11(3):243–253. doi:10.1109/TVCG.2005.49
14. Goddard, T. (2013). Molecular Surface Algorithms (design review, UCSF Chimera).
    https://www.cgl.ucsf.edu/chimera/data/surface-oct2013/surface.html
15. ChimeraX surface module docs + source (`gridsurf.py`).
    https://chimerax.readthedocs.io/en/latest/modules/surface/surface.html ;
    https://github.com/RBVI/ChimeraX
16. Khan, D., Yan, D.-M., Gui, S., Lu, B., Zhang, X. (2018). Molecular Surface Remeshing with
    Local Region Refinement. *Int. J. Mol. Sci.* 19(5):1383. doi:10.3390/ijms19051383
    — https://pmc.ncbi.nlm.nih.gov/articles/PMC5983798/
17. Varshney, A., Brooks, F. P., Wright, W. V. (1994). Linearly Scalable Computation of Smooth
    Molecular Surfaces / Computing Smooth Molecular Surfaces. *IEEE CG&A*.
    https://ieeexplore.ieee.org/document/310720 — VMD `surf` backend.
18. VMD User's Guide, Surf / NanoShaper representations.
    https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.4/ug/node74.html ;
    https://www.ks.uiuc.edu/Research/vmd/vmd-1.9.4/ug/node76.html
19. Liang, J., Edelsbrunner, H., Fu, P., Sudhakar, P., Subramaniam, S. (1998). Analytical
    shape computation of macromolecules: Molecular area and volume through alpha shape.
    *Proteins* 33:1–17.
20. Campbell, M. J. (2026). RustSASA: A Rust Crate for Accelerated Solvent Accessible Surface
    Area Calculations. *JOSS* 11(117):9537. doi:10.21105/joss.09537
    — https://joss.theoj.org/papers/10.21105/joss.09537 ;
    https://github.com/maxall41/rustsasa
21. powersasa (Rust crate). https://lib.rs/crates/powersasa
22. isosurface (Rust crate). https://github.com/swiftcoder/isosurface ;
    https://docs.rs/isosurface/latest/isosurface/
23. mcubes (Rust crate). https://lib.rs/crates/mcubes
24. Can, T., Chen, C.-I., Wang, Y. F. (2006). Efficient molecular surface generation using
    level-set methods. *J. Mol. Graph. Model.* 25(4):442–454. (LSMS; evaluated in [14].)
25. Kim, B., Kim, K. J., Seong, J. (2012). GPU accelerated molecular surface computing.
    *Appl. Math. Inf. Sci.* 6:185S–194S. (cited in [2])
26. Gainza, P., et al. (2020). MaSIF. *Nature Methods*. (MSMS surfaces as ML input; per
    https://ccsb.scripps.edu/msms/)
27. NanoShaperWeb (2025). *J. Chem. Inf. Model.* doi:10.1021/acs.jcim.5c00821
28. NanoShaper-VMD interface (2019). *Bioinformatics* 35(7):1241–1243.
    https://pmc.ncbi.nlm.nih.gov/articles/PMC6449750/
