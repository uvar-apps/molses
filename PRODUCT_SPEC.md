# Product Specification — `molses`

**Solvent Excluded Surface (SES / Connolly surface) generation library and CLI in Rust**

| Field | Value |
|---|---|
| Status | Draft for review (pre-implementation) |
| Date | 2026-11-12 |
| Working name | `molses` (library crate `molses`, CLI binary `molses`) |
| Deliverable | Cargo workspace: reusable library + command-line tool |
| Inputs | PDB, mmCIF/PDBx, XYZ molecular structure files |
| Output | Wavefront OBJ triangle mesh with per-vertex normals |
| Out of scope (this release) | Source code — see companion documents `DESIGN.md` and `IMPLEMENTATION_PLAN.md` |

---

## 1. Purpose

Molecular surfaces are the primary geometric abstraction of a molecule's shape. The
**Solvent Excluded Surface (SES)** — also called the **Connolly surface**, **molecular
surface**, or **Richards surface** — is the surface traced by a spherical solvent probe as
it rolls over the van der Waals (VdW) spheres of the atoms [Richards 1977; Connolly 1983;
Sanner et al. 1996]. It is the boundary of the volume from which solvent molecules are
excluded, and it is the standard geometric input to:

- visualization of molecular shape, pockets, and interaction interfaces;
- implicit-solvent electrostatics (Poisson–Boltzmann / BEM / FEM) [Khan et al. 2018];
- geometric deep learning on proteins (e.g. MaSIF consumes MSMS meshes) [Gainza et al. 2020];
- shape comparison, docking, and cavity/channel analysis [Voss & Gerstein 2010].

Today, high-quality SES generation is dominated by decades-old tools (MSMS, 1996) and by
C/C++/CUDA packages (NanoShaper, EDTSurf, QuickSES, SESDF). There is no production-quality,
well-tested, permissively licensed **Rust** implementation of SES mesh generation that can
be embedded as a library. This product fills that gap.

**Product statement.** *`molses` computes a resolution-controllable, watertight triangle
mesh of the Solvent Excluded Surface of a molecule from PDB, mmCIF, or XYZ input, using
concurrent and SIMD-optimized algorithms, and exports it as an OBJ file with per-vertex
normals — both as a reusable Rust library and as a command-line tool.*

---

## 2. Users and use cases

| User | Use case | Success looks like |
|---|---|---|
| Computational chemist / structural biologist | Produce an SES mesh of a PDB entry for visualization in Blender / VMD / ChimeraX / MeshLab | `molses 1crn.pdb -o 1crn.obj` produces a mesh that opens cleanly and is immediately viewable |
| Method developer (PB/BEM/FEM solvers) | Feed a watertight, oriented, high-quality mesh into a boundary-element solver | Manifold, consistently oriented mesh; stable area/volume under resolution refinement |
| Rust application developer | Embed SES generation in a larger pipeline (pocket detection, docking, ML featurization) | `molses` as a library dependency; documented API; no C toolchain; controllable memory |
| Scientific-computing engineer | Batch-generate surfaces for hundreds of structures with controlled LOD | CLI accepts resolution presets; deterministic output; parallel across models/files |
| Educator / student | Understand rolling-probe SES geometry | Accurate geometry at default parameters; documented references; reproducible |

---

## 3. What is the SES (normative definition)

> This section is normative: it defines exactly which surface the product must compute.
> Implementations that produce a different surface are non-conforming.

### 3.1 Setup

Let a molecule be a finite set of atoms indexed by `i = 1..M`, with van der Waals
centers `c_i ∈ ℝ³` and van der Waals radii `r_i > 0` (in Ångström). Let `r_p > 0` be the
**probe radius** (radius of the spherical solvent molecule; default 1.4 Å for water
[Richards 1977]).

Denote by `B(c, R) = { p ∈ ℝ³ : ‖p − c‖ ≤ R }` the closed ball of radius `R` about `c`.

### 3.2 The three surfaces

1. **Van der Waals surface.** The topological boundary of the union of atomic balls:
   `Ω_vdW = ⋃_i B(c_i, r_i)`, surface = `∂Ω_vdW`.
   [Quan & Stamm 2016, §3]

2. **Solvent Accessible Surface (SAS).** The locus of the *center* of the probe as it rolls
   over the VdW surface. Equivalently, the boundary of the union of inflated balls:
   `Ω_SAS = ⋃_i B(c_i, r_i + r_p)`, surface = `∂Ω_SAS`.
   [Lee & Richards 1971; Sanner et al. 1996, Fig. 1; Quan & Stamm 2016, §3]

3. **Solvent Excluded Surface (SES).** The boundary of the volume the probe *body* can never
   occupy. Equivalently: the surface traced by the outer boundary of the probe sphere as it
   rolls over the VdW surface while remaining tangent to it and non-overlapping all atoms.
   [Richards 1977; Greer & Bush 1978; Sanner et al. 1996]

### 3.3 Composition of the SES

The SES is a closed, orientable, piecewise-smooth surface composed of three kinds of
patches (Connolly 1983; Richards 1977):

| Patch | Geometry | Where the probe is |
|---|---|---|
| **Contact / convex spherical** | Spherical patch of radius `r_i` on atom `i` | Touching exactly one atom |
| **Toric / toroidal reentrant** | Saddle-shaped patch of a torus (minor radius `r_p`) | Touching exactly two atoms |
| **Spheric reentrant / concave** | Spherical triangle of radius `r_p` | Touching three (or more) atoms |

The SES is C¹ (continuous tangent plane) almost everywhere but **not** C²: curvature jumps
at the seams between patch types. This is intrinsic to the SES definition, not a numerical
defect [Goddard 2013; Lange et al. 2019].

### 3.4 Exact implicit characterization (normative for grid engines)

Let `K = Ω_SAS = ⋃_i B(c_i, r_i + r_p)` and let `f_sas : ℝ³ → ℝ` be the **signed distance
function** to `∂K` (negative inside `K`, positive outside). Then:

```
f_ses(p) = f_sas(p) + r_p
```

and

```
SES = { p ∈ ℝ³ : f_ses(p) = 0 } = { p : f_sas(p) = −r_p }.
```

[Quan & Stamm 2016, Eq. 3.4 and Eq. 5.21; free preprint: HAL hal-01205580]

**Critical implementation note.** The naive field `min_i (‖p−c_i‖ − (r_i+r_p))` has the SAS
as its zero level set but is **not** a signed distance function inside `K` (its gradient
magnitude is not 1 in overlap regions). Using it directly as `f_sas` does **not** yield the
correct SES. Quan & Stamm show a counter-example (their Figure 2) to a related naive
formulation. The implementation must use a true signed distance to `K` — see `DESIGN.md`.

### 3.5 Degenerate and singular cases

- SES singularities (cusps on toroidal patches, singular circles on concave patches) occur
  where the probe contacts 4+ atoms or where patches self-intersect. A complete
  characterization is given by Quan & Stamm (2016, §5.2). MSMS handles special cases with
  heuristics and is documented to fail on several of them [MSMS documentation].
- Grid/implicit implementations are **not required** to resolve analytic singularities;
  they must instead produce a mesh converging to the SES as resolution increases, with
  features below the grid scale smoothed at scale `O(h)`.
- Atoms whose balls are strictly contained in another atom's ball must be merged/ignored
  (MSMS merges identical centers and crashes on nested spheres [MSMS docs]; we must not
  crash).

### 3.6 Terminology mapping (for documentation and CLI help)

| This product | Synonyms in the literature |
|---|---|
| SES | Connolly surface; molecular surface; Richards surface; excluded surface; `r`-excluded surface `E_r(S)` |
| SAS | Lee–Richards surface; accessible surface; `r`-accessible surface `A_r(S)` |
| Probe radius `r_p` | `probe_radius`, `probe`, `rp`; MSMS `-probe_radius`; NanoShaper probe; water ≈ 1.4 Å |

---

## 4. Functional requirements

### 4.1 Structure input

| ID | Requirement | Priority |
|---|---|---|
| F-IN-1 | Read **PDB** files (fixed-column, ATOM/HETATM records, multi-MODEL, altLoc, element column). | Must |
| F-IN-2 | Read **mmCIF / PDBx** files (the wwPDB successor format; category-based `_atom_site`). | Must |
| F-IN-3 | Read **XYZ** files (element + x y z per line; optional comment; multi-frame trajectories as a stretch goal). | Must |
| F-IN-4 | Infer element per atom: prefer the element field (`_atom_site.type_symbol` / PDB cols 77–78); fall back to atom-name heuristic (C/H/N/O/S). | Must |
| F-IN-5 | Assign van der Waals radii from a **Bondi (1964)** table keyed by element [Bondi 1964]; support a user-supplied override table (CSV/JSON). | Must |
| F-IN-6 | Atom-selection controls: include/exclude hydrogens, HETATM/ligands, waters (HOH/WAT), specific chains, altLoc policy (highest occupancy / altLoc A / all), model selection (`--model N`, `--model all`, default first). | Must |
| F-IN-7 | Handle insertion codes without corrupting atom identity. | Must |
| F-IN-8 | Tolerate duplicate atom centers (merge, keep largest radius — MSMS behavior) and fully nested spheres (ignore the inner one). | Must |
| F-IN-9 | Report, not crash, on unknown elements (skip with a warning; require explicit radius or skip). | Must |
| F-IN-10 | Accept gzip-compressed inputs (`.pdb.gz`, `.cif.gz`). | Should |
| F-IN-11 | Read MSMS-style `.xyzr` sphere files (x y z r per line) as a low-level interchange input. | Should |

### 4.2 Geometry parameters

| ID | Requirement | Priority |
|---|---|---|
| F-PAR-1 | `--probe-radius <R>` in Å. Default **1.4** (water) [Richards 1977]. Valid range (0, 10]; reject non-finite / non-positive. | Must |
| F-PAR-2 | `--resolution` / `--grid-spacing <h>` in Å, controlling level of detail. Named presets (see §6). Default `balanced` (0.35 Å). | Must |
| F-PAR-3 | `--radii bondi|custom:<file>` — radii table selection. | Must |
| F-PAR-4 | Deterministic results for the same input and parameters (independent of thread count). | Must |
| F-PAR-5 | `--surface ses|vdw|sas` — compute SES (default), or the VdW surface (`probe = 0`), or the SAS (for comparison/area). | Should |

### 4.3 Surface generation

| ID | Requirement | Priority |
|---|---|---|
| F-SUR-1 | Produce a triangle mesh approximating the SES as defined in §3, converging to it as `h → 0`. | Must |
| F-SUR-2 | Mesh is **watertight** (every edge shared by exactly two triangles) and **consistently oriented** (outward normals: CCW when viewed from outside the SES cavity). | Must |
| F-SUR-3 | Per-vertex unit normals, computed so as to be usable for smooth shading (see §8). | Must |
| F-SUR-4 | Level of detail controlled by grid spacing `h`; report estimated geometric error `O(h)` and resulting triangle count. | Must |
| F-SUR-5 | Emit the **outer component** by default; `--all-components` to include internal cavities; `--no-cavities` to drop enclosed voids below a volume threshold. | Must |
| F-SUR-6 | Per-vertex nearest-atom index available in the library API (and in extended output) for per-atom/per-residue coloring. | Should |
| F-SUR-7 | Per-vertex patch-type classification (contact / toric / spheric) as metadata where the engine can provide it; best-effort on grid engines. | Could |
| F-SUR-8 | Optional Laplacian/Taubin smoothing pass (`--smooth <n>`), off by default (trades fidelity for visual quality). | Should |
| F-SUR-9 | Compute and report SES area (mesh-based) and enclosed volume; compare to analytic values on validation cases. | Should |

### 4.4 Output

| ID | Requirement | Priority |
|---|---|---|
| F-OUT-1 | **Wavefront OBJ** export: `v x y z` (Å), `vn nx ny nz` (unit), `f a//a b//b c//c` (1-based, CCW-outward). | Must |
| F-OUT-2 | OBJ header comments recording: tool version, input path, probe radius, grid spacing, atom count, triangle/vertex count, generation timestamp. | Must |
| F-OUT-3 | Optional `g` groups per connected component. | Should |
| F-OUT-4 | Optional **MSMS-compatible `.vert`/`.face`** export (per-vertex xyz + normal + closest-atom + patch type) for interop with VMD/Chimera pipelines. | Should |
| F-OUT-5 | Optional PLY export (vertex normals, optional per-atom color). | Could |
| F-OUT-6 | Writes to stdout or to a file (`-o`); refuse to overwrite unless `--force`. | Must |
| F-OUT-7 | Library API returns an in-memory `Mesh` structure (positions, normals, indices, metadata) without forcing file I/O. | Must |

### 4.5 Concurrency and performance

| ID | Requirement | Priority |
|---|---|---|
| F-PER-1 | Multi-threaded execution on shared-memory CPUs (Rayon), with `--threads <n>` / `RAYON_NUM_THREADS` control. | Must |
| F-PER-2 | SIMD acceleration of the hot distance-field kernels on stable Rust (no nightly features required). | Must |
| F-PER-3 | Deterministic output regardless of thread count (bitwise-identical mesh). | Must |
| F-PER-4 | Memory use predictable and documented; `--max-memory <MB>` soft budget that degrades resolution rather than aborting (warn + raise `h`) when exceeded. | Should |
| F-PER-5 | Serial fallback when built without the `parallel` feature. | Must |
| F-PER-6 | Scalar fallback when built without the `simd` feature. | Must |

### 4.6 Library API (summary; full design in `DESIGN.md`)

| ID | Requirement | Priority |
|---|---|---|
| F-API-1 | Idiomatic Rust API: `Structure` → `SesOptions` → `Mesh`. Pure library, no process-global state, no C dependencies. | Must |
| F-API-2 | Feature flags: `parallel` (rayon), `simd` (portable SIMD), `serde`, `mmap`. Defaults = `parallel` + `simd`. | Must |
| F-API-3 | Typed error enum (`thiserror`), no panics on malformed input. | Must |
| F-API-4 | Works as a plain `cargo add` dependency on stable Rust (MSRV declared, tested). | Must |
| F-API-5 | The CLI must use only public library APIs (proves the library is sufficient). | Must |

---

## 5. Non-functional requirements

| ID | Category | Requirement |
|---|---|---|
| N-1 | Correctness | Mesh Hausdorff distance to analytic SES ≤ `2h` on validation molecules (see §7). |
| N-2 | Correctness | Mesh area within **5%** of MSMS/analytic SES area at `h = 0.3 Å` on validation set; within 2% at `h = 0.15 Å`. (EDTSurf reports 2–4% at comparable resolutions [Xu & Zhang 2009].) |
| N-3 | Robustness | Never panic, hang, or silently emit a non-watertight mesh on any wwPDB entry in the test corpus; always either succeed or return a structured error. |
| N-4 | Performance | 100 000-atom structure at `h = 0.5 Å` completes in **≤ 5 s** on a contemporary 8-core desktop CPU (target: beat Chimera's unoptimized ~30 s single-core reference [Goddard 2013] by ≥ 6×). |
| N-5 | Performance | 10 000-atom structure at `h = 0.5 Å` completes in **≤ 0.5 s**. |
| N-6 | Scalability | 500 000-atom virus capsid at `h = 1.0 Å` completes within available RAM using slab chunking (documented memory model). |
| N-7 | Determinism | Same input + parameters ⇒ byte-identical OBJ (modulo header timestamp, which must be suppressible via `SOURCE_DATE_EPOCH`-style flag). |
| N-8 | Portability | Linux, macOS, Windows; x86-64 and aarch64; no nightly Rust; no C toolchain required. |
| N-9 | Maintainability | ≥ 80% line coverage on the geometry core; golden-file mesh tests; criterion benchmarks in CI. |
| N-10 | Licensing | MIT OR Apache-2.0 dual license; all dependencies permissively licensed; test data CC0 (wwPDB archive). |

---

## 6. Resolution / level-of-detail model

The primary LOD control is the **grid spacing** `h` (Å) — the sampling distance of the
implicit field from which the mesh is extracted. This is the same knob used by ChimeraX
(`grid_spacing`, default 0.5 Å) and QuickSES (`-v`, default 0.5 Å).

**Named presets** (overridable with an explicit number):

| Preset | `h` (Å) | Intended use | Expected triangle count for a 25 000 Å² protein |
|---|---|---|---|
| `draft` | 0.75 | Interactive preview, thumbnails | ~5 × 10⁴ |
| `coarse` | 0.50 | Default visualization (MSMS-comparable density) | ~1 × 10⁵ |
| `balanced` | 0.35 | **Default** — good visual quality, moderate cost | ~2 × 10⁵ |
| `fine` | 0.25 | High-quality rendering, publication figures | ~4 × 10⁵ |
| `ultra` | 0.15 | Near-analytic geometry, BEM/FEM preprocessing | ~1 × 10⁶ |

Triangle count scales as `≈ Area / (h² · √3 / 2)` for a quasi-uniform isosurface mesh; the
table above is indicative. Geometric error is `O(h)` in Hausdorff distance (the field is
1-Lipschitz); measured error is typically well below that bound in smooth regions.

**Secondary controls:**

- `--max-triangles <N>` — auto-select `h` so the estimated triangle count stays under `N`
  (requires a cheap area estimate: e.g. Shrake–Rupley or a coarse first pass).
- `--smooth <n>` — `n` Laplacian/Taubin iterations after extraction (visual LOD, not
  geometric accuracy).

**Why grid spacing and not vertices/Å²?** MSMS's `-density` (vertices/Å²) applies to an
analytic triangulation and cannot be reproduced exactly on a grid mesh. We expose `h`
honestly (it is what the algorithm controls) and document the approximate correspondence
(`h = 0.5 Å` ≈ MSMS `density = 1–2` in visual density; Goddard 2013).

---

## 7. Validation and accuracy criteria

Validation set (all wwPDB / CC0):

| Entry | Atoms (approx.) | Why |
|---|---|---|
| 1CRN (crambin) | 327 | Classic MSMS benchmark [Sanner et al. 1996] |
| 4HHB (hemoglobin) | ~4 700 | Multi-chain, altLocs |
| 1TLP / thermolysin | ~2 400 | Second MSMS benchmark point |
| 1BBL or an NMR ensemble | ~1 000 | Multi-model handling |
| 1CBS | ~800 | Waters + ligand |
| 1K4R or 3J3Q (subset) | large | Scalability smoke test |

Analytic ground truths to test against:

1. **Single atom:** SES = sphere of radius `r + r_p`. Area = `4π(r + r_p)²`. Exact.
2. **Two identical atoms, distance `d`:** analytic area available from Richmond (1984) /
   Connolly (1983) formulas for overlapping spheres.
3. **Comparison to MSMS** on 1CRN and 1TLP at `probe = 1.5`, `density = 1.0` — area and
   Hausdorff distance of our `h = 0.3` mesh vs MSMS `.vert/.face`.
4. **Monotone convergence:** area(h) and volume(h) converge as `h` decreases; |Δ| < 2%
   between `h = 0.25` and `h = 0.15`.
5. **Topological sanity:** Euler characteristic of the outer component = `2 − 2g` for a
   closed orientable surface; watertightness; genus stable under refinement.

---

## 8. Mesh quality requirements

Derived from downstream consumers (visualization and BEM/FEM) [Khan et al. 2018]:

| Property | Requirement |
|---|---|
| Watertight | Every edge belongs to exactly two triangles. |
| Manifold | No non-manifold edges or vertices on the outer component. |
| Orientation | Consistently outward (CCW from outside); right-hand normals. |
| Normals | Per-vertex, unit length; smooth-shading usable. |
| Degenerate triangles | Zero-area triangles removed; sliver count reported. |
| Aspect ratio | Not enforced in v1 (grid meshes have inherent slivers [Goddard 2013]); documented, with a remeshing hook as future work. |
| Vertex density | Controlled by `h` (§6); no unnecessary oversampling. |

**Normals computation (normative):** primary method is the **gradient of the implicit
field** at the vertex (central differences / Sobel on the signed-distance grid), which for
a true SDF is the exact surface normal up to discretization. Fallback: angle-weighted
pseudonormals [Bærentzen & Aanæs 2005]. Note that even exact SES normals show shading
seams at C¹-but-not-C² patch joins [Goddard 2013] — this is a property of the SES, not a
bug.

---

## 9. CLI specification (summary)

```
molses [OPTIONS] <INPUT>

Arguments:
  <INPUT>                 Structure file: .pdb / .ent / .cif / .mmcif / .xyz / .xyzr [.gz]

Options:
  -o, --output <FILE>     Output OBJ path [default: stdout]
      --probe-radius <R>  Probe radius in Å [default: 1.4]
      --resolution <H|P>  Grid spacing in Å or preset: draft|coarse|balanced|fine|ultra
                          [default: balanced]
      --radii <SET>       Van der Waals radii: bondi | custom:<file> [default: bondi]
      --model <SEL>       Model selection: <N> | all | first [default: first]
      --chains <LIST>     Comma-separated chain IDs to include
      --no-hydrogens      Exclude hydrogen atoms [default: true]
      --keep-hydrogens    Include hydrogen atoms
      --no-hetatm         Exclude HETATM records
      --include-waters    Include waters (HOH/WAT) [default: excluded]
      --altloc <POLICY>   highest-occupancy | first | all [default: highest-occupancy]
      --surface <TYPE>    ses | vdw | sas [default: ses]
      --all-components    Include internal cavity surfaces
      --smooth <N>        Laplacian smoothing iterations [default: 0]
      --threads <N>       Worker threads [default: all]
      --max-memory <MB>   Soft memory budget
      --format <FMT>      obj | msms | ply [default: obj]
      --force             Overwrite existing output
      --quiet             Suppress progress output
      --verbose           Detailed progress
  -h, --help              Help
  -V, --version           Version
```

Exit codes: `0` success; `1` invalid input/arguments; `2` I/O error; `3` geometry failure;
`4` resource limit exceeded.

---

## 10. Competitive landscape

| Tool | Language | Algorithm | Strengths | Weaknesses for our users |
|---|---|---|---|---|
| **MSMS** (Sanner 1996) | C++ | Analytic reduced surface → patches → triangulation | Exact patches; de facto interchange (`.vert/.face`); 30 years of use | Non-commercial license; documented bugs & failures >10k atoms [Goddard 2013]; no OBJ with normals |
| **NanoShaper** (Decherchi & Rocchia 2013) | C++/CGAL | Alpha-shape analytic + ray-cast MC | Robust; cavity/pocket detection; multi-definition | GPL; CGAL dependency; heavy; C++ only |
| **EDTSurf** (Xu & Zhang 2009) | C/C++ | EDT + MC/VCMC | Fast; open source | Grid-size bound; PLY native; aging code |
| **QuickSES / Hermosilla 2018** | CUDA | GPU grid + MC | Very fast on GPU | CUDA-only; OBJ without normals; GPU required |
| **SESDF** (Plateau-Holleville 2024) | C++/CUDA | Analytic SES on GPU | State-of-the-art performance & memory | Reference code; CUDA; not a library |
| **ChimeraX `gridsurf`** | Python/C++ | Two-pass distance grid | Clean readable reference | Embedded in ChimeraX (LGPL); not a standalone library |
| **RustSASA** (2026) | Rust | Shrake–Rupley SASA | Excellent Rust precedent (rayon+pulp) | **Areas only — no mesh** |
| **`molframe-surface` 0.7** | Rust | Flood + EDT + 6-tet polygonisation | Closest prior art: SES mesh + OBJ in Rust | Pre-1.0, single maintainer, grid-only, part of a large engine; no analytic mode; little adoption |

**Positioning.** `molses` is the first *dedicated*, permissively licensed (MIT OR
Apache-2.0) Rust **library + CLI** whose sole product is SES mesh generation. Relative to
`molframe-surface` (the only other Rust SES mesher, verified on docs.rs at time of
writing), `molses` differentiates on: focused scope, documented LOD model and accuracy
targets, deterministic output, tested robustness on wwPDB corpus, MSMS-format interop,
and a designed public API. We treat `molframe-surface` as a validation reference, not a
dependency (pre-1.0, single-maintainer risk) [docs.rs/molframe-surface].

---

## 11. Out of scope (v1)

Explicitly **not** in the first release (candidates for later milestones):

- Analytic/parametric Connolly triangulation (exact quadric patches, MSMS-style) — planned
  as an optional `exact` engine in a later release (see `IMPLEMENTATION_PLAN.md` M6).
- GPU / CUDA / Metal backends.
- SASA-only computation as a product feature (though the area computation shares code;
  RustSASA remains the specialist).
- Gaussian / Skin / Blobby surface definitions (NanoShaper-style alternatives).
- Pocket detection, cavity volume listing, channel analysis (NanoShaper / 3V territory).
- Mesh remeshing / decimation for FEM quality (Khan et al. 2018) — a `--smooth` pass is
  the only concession in v1.
- Python bindings (via PyO3) — planned post-1.0.
- Reading trajectories (XTC/DCD); only multi-model PDB/mmCIF ensembles in v1.
- Color/property mapping in the OBJ (per-atom coloring is exposed via API metadata).

---

## 12. Assumptions and open questions

| # | Item | Current stance |
|---|---|---|
| 1 | Default probe radius 1.4 Å (water) vs MSMS's 1.5 Å | **1.4 Å** is the field standard for water [Richards 1977]; MSMS's 1.5 default is documented as historical. CLI override is trivial. |
| 2 | Primary algorithm family: grid SDF vs analytic patches | **Grid SDF first** (robust, LOD-native, SIMD/rayon-friendly, matches user's performance requirement); analytic later. Justified in `DESIGN.md`. |
| 3 | Crate/binary name `molses` | Working name; rename before crates.io publish if the name is taken or unwanted. |
| 4 | Whether `--no-hydrogens` should default to true | Yes: Bondi united-atom practice; X-ray entries often lack H; users with H can opt in. |
| 5 | Exact MSMS binary license | MSMS is source-available, non-commercial by default [ccsb.scripps.edu/msms]; we do **not** redistribute MSMS; we only implement a compatible `.vert/.face` *writer*. |
| 6 | `std::simd` stabilization | Unstable as of 2026-11; use `wide` (or `pulp`) behind a `simd` feature; keep kernels portable. |

---

## 13. References

### Primary theory

1. Lee, B. & Richards, F.M. (1971). The interpretation of protein structures: estimation of static accessibility. *J. Mol. Biol.* 55, 379–400. doi:[10.1016/0022-2836(71)90324-X](https://doi.org/10.1016/0022-2836(71)90324-X)
2. Shrake, A. & Rupley, J.A. (1973). Environment and exposure to solvent of protein atoms. *J. Mol. Biol.* 79, 351–371. doi:[10.1016/0022-2836(73)90011-9](https://doi.org/10.1016/0022-2836(73)90011-9)
3. Richards, F.M. (1977). Areas, volumes, packing, and protein structure. *Annu. Rev. Biophys. Bioeng.* 6, 151–176. doi:[10.1146/annurev.bb.06.060177.001055](https://doi.org/10.1146/annurev.bb.06.060177.001055). Free PDF: <https://biomolmd.org/mw/images/e/e5/Richards77.pdf>
4. Greer, J. & Bush, B.L. (1978). Macromolecular shape and surface maps by solvent exclusion. *PNAS* 75(1), 303–307. doi:[10.1073/pnas.75.1.303](https://doi.org/10.1073/pnas.75.1.303)
5. Connolly, M.L. (1983). Analytical molecular surface calculation. *J. Appl. Cryst.* 16, 548–558. doi:[10.1107/S0021889883010985](https://doi.org/10.1107/S0021889883010985)
6. Connolly, M.L. (1983). Solvent-accessible surfaces of proteins and nucleic acids. *Science* 221, 709–713. doi:[10.1126/science.6879170](https://doi.org/10.1126/science.6879170)
7. Richmond, T.J. (1984). Solvent accessible surface area and excluded volume in proteins. *J. Mol. Biol.* 178, 63–89. doi:[10.1016/0022-2836(84)90231-6](https://doi.org/10.1016/0022-2836(84)90231-6)
8. Bondi, A. (1964). van der Waals volumes and radii. *J. Phys. Chem.* 68(3), 441–451. doi:[10.1021/j100785a001](https://doi.org/10.1021/j100785a001)
9. Edelsbrunner, H. (1995). The union of balls and its dual shape. *Discrete & Comput. Geom.* 13, 415–440. doi:[10.1007/BF02574053](https://doi.org/10.1007/BF02574053)
10. Sanner, M.F., Olson, A.J. & Spehner, J.-C. (1996). Reduced surface: an efficient way to compute molecular surfaces. *Biopolymers* 38, 305–320. doi:[10.1002/(SICI)1097-0282(199603)38:3<305::AID-BIP4>3.0.CO;2-Y](https://doi.org/10.1002/(SICI)1097-0282(199603)38:3%3C305::AID-BIP4%3E3.0.CO;2-Y). Free PDF: <https://biomolmd.org/mw/images/4/49/Sanner96.pdf> (MSMS)
11. Liang, J. et al. (1998). Analytical shape computation of macromolecules: I. molecular area and volume through alpha shape. *Proteins* 33, 1–17. doi:[10.1002/(SICI)1097-0134(19981001)33:1<1::AID-PROT1>3.0.CO;2-O](https://doi.org/10.1002/(SICI)1097-0134(19981001)33:1%3C1::AID-PROT1%3E3.0.CO;2-O)
12. Quan, C. & Stamm, B. (2016). Mathematical analysis and calculation of molecular surfaces. *J. Comput. Phys.* 322, 221–239. doi:[10.1016/j.jcp.2016.07.007](https://doi.org/10.1016/j.jcp.2016.07.007). Free preprint: <https://hal.sorbonne-universite.fr/hal-01205580/document>
13. Edelsbrunner, H. (2001). *Geometry and Topology for Mesh Generation*. Cambridge University Press. ISBN 0521793098.

### Algorithms and implementations

14. Lorensen, W.E. & Cline, H.E. (1987). Marching cubes: a high resolution 3D surface construction algorithm. *SIGGRAPH '87*, 163–169. doi:[10.1145/37402.37422](https://doi.org/10.1145/37402.37422)
15. Xu, D. & Zhang, Y. (2009). Generating triangulated macromolecular surfaces by Euclidean distance transform. *PLoS ONE* 4(12), e8140. doi:[10.1371/journal.pone.0008140](https://doi.org/10.1371/journal.pone.0008140) (EDTSurf)
16. Decherchi, S. & Rocchia, W. (2013). A general and robust ray-casting-based algorithm for triangulating surfaces at the nanoscale. *PLoS ONE* 8(4), e59744. doi:[10.1371/journal.pone.0059744](https://doi.org/10.1371/journal.pone.0059744) (NanoShaper)
17. Goddard, T. (2013). Molecular Surface Algorithms (design review). UCSF Chimera. <https://www.cgl.ucsf.edu/chimera/data/surface-oct2013/surface.html>
18. Bærentzen, J.A. & Aanæs, H. (2005). Signed distance computation using the angle weighted pseudonormal. *IEEE TVCG* 11(3), 243–253. doi:[10.1109/TVCG.2005.49](https://doi.org/10.1109/TVCG.2005.49)
19. Hermosilla, P. et al. (2018). Interactive GPU-based generation of solvent excluded surfaces. *Computers & Graphics*. doi:[10.1007/s00371-017-1397-2](https://doi.org/10.1007/s00371-017-1397-2) (QuickSES algorithm)
20. Martinez, X., Krone, M. & Baaden, M. (2019). QuickSES: a library for fast computation of solvent excluded surfaces. *Eurographics*. <https://hal.archives-ouvertes.fr/hal-02370900/document>. Code: <https://github.com/nezix/QuickSES>
21. Plateau-Holleville, C. et al. (2024). Efficient GPU computation of large protein solvent-excluded surface. *IEEE TVCG*. doi:[10.1109/TVCG.2024.3380100](https://doi.org/10.1109/TVCG.2024.3380100) (SESDF). Code: <https://github.com/PlathC/SESDF>
22. Khan, D. et al. (2018). Molecular surface remeshing with local region refinement. *Int. J. Mol. Sci.* 19(5), 1383. doi:[10.3390/ijms19051383](https://doi.org/10.3390/ijms19051383)
23. Voss, N.R. & Gerstein, M. (2010). 3V: cavity, channel and cleft volume calculator and extractor. *Nucleic Acids Res.* 38, W555–W562. doi:[10.1093/nar/gkq395](https://doi.org/10.1093/nar/gkq395)
24. Lange, A.W., Herbert, J.M., Albrecht, B.J. & You, Z.-Q. (2019). Intrinsically smooth discretisation of Connolly's solvent-excluded molecular surface. *Mol. Phys.* 118, e1644384. doi:[10.1080/00268976.2019.1644384](https://doi.org/10.1080/00268976.2019.1644384)

### Software documentation

25. MSMS documentation, Sanner Lab / Scripps. <https://sannerlab.github.io/software/msms/documentation.html> (and <https://ccsb.scripps.edu/msms/>)
26. ChimeraX `ses_surface_geometry` API. <https://chimerax.readthedocs.io/en/latest/modules/surface/surface.html> (implementation: `gridsurf.py` in <https://github.com/RBVI/ChimeraX>)
27. NanoShaper 0.7 User Guide; tutorial: Decherchi & Rocchia (2015), *Methods Mol. Biol.*, doi:[10.1007/978-3-319-12211-3_10](https://doi.org/10.1007/978-3-319-12211-3_10)
28. EDTSurf: <https://zhanggroup.org/EDTSurf/>
29. Wavefront OBJ format. <https://en.wikipedia.org/wiki/Wavefront_.obj_file>; Library of Congress format description <https://www.loc.gov/preservation/digital/formats/fdd/fdd000507.shtml>

### Rust ecosystem (surveyed 2026-11)

30. `pdbtbx` 0.12 — PDB/mmCIF parser. <https://crates.io/crates/pdbtbx> · <https://github.com/douweschulte/pdbtbx>
31. RustSASA (`rust-sasa`) — SASA, rayon+pulp. Campbell, M.J. (2026), *JOSS* 11(117):9537, doi:[10.21105/joss.09537](https://doi.org/10.21105/joss.09537)
32. `molframe-surface` 0.7 — existing Rust grid SES + OBJ (prior art / validation reference). <https://docs.rs/molframe-surface/latest/molframe_surface/>
33. `glam`, `rayon`, `wide`, `pulp`, `criterion`, `thiserror`, `clap` — supporting crates. (See `DESIGN.md` for versions and rationale.)
34. wwPDB usage policies (PDB archive data is CC0 1.0). <https://www.wwpdb.org/about/usage-policies>

---

*End of product specification. Companion documents: `DESIGN.md` (technical design),
`IMPLEMENTATION_PLAN.md` (work breakdown and milestones), `research/` (full research
findings with source URLs).*
