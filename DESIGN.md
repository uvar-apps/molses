# Design Document — `molses`

**Technical design for the SES/Connolly surface library and CLI**

| Field | Value |
|---|---|
| Status | Draft for review (pre-implementation) |
| Date | 2026-11-12 |
| Companion docs | `PRODUCT_SPEC.md` (requirements), `IMPLEMENTATION_PLAN.md` (work plan), `research/` (evidence) |

---

## 1. Design goals and non-goals

### 1.1 Goals

1. **Correct SES geometry** as normatively defined in `PRODUCT_SPEC.md` §3, with a
   documented accuracy/LOD relationship.
2. **Robustness over cleverness**: never crash on real wwPDB inputs; structured errors;
   degenerate atoms handled explicitly.
3. **Performance through parallelism and SIMD** on stable Rust (Rayon + portable SIMD),
   with deterministic output.
4. **Library-first**: the CLI is a thin shell over a public API another project can depend on.
5. **Small, auditable dependency surface**: pure Rust, permissive licenses, no C toolchain.

### 1.2 Non-goals (v1)

GPU backends; analytic quadric-patch triangulation (designed for, not implemented — §12);
FEM-grade remeshing; Python bindings; trajectory formats.

---

## 2. Algorithm family selection (and why)

Two established families compute the SES:

| Family | Representative | LOD control | Parallelism | Robustness | Effort |
|---|---|---|---|---|---|
| **A. Analytic / reduced surface** | Connolly 1983, MSMS 1996, SESDF 2024 | triangulation density (vertices/Å²) | harder (irregular work) | fragile (singularities) | very high |
| **B. Implicit / grid (signed distance field + isosurface)** | ChimeraX gridsurf, EDTSurf 2009, QuickSES 2018, `molframe-surface` 0.7 | grid spacing `h` (Å) | excellent (regular grids) | robust | moderate |

**Decision: family B (grid SDF) is the v1 engine.** Rationale:

- The user requirements emphasize *concurrency + vector instructions* and *resolution/LOD
  control*. Regular grids are the textbook data-parallel workload; `h` is an honest,
  directly controllable LOD knob.
- The correct implicit function is fully characterized in the literature
  (`f_ses = f_sas + r_p`, Quan & Stamm 2016), so the target surface is well defined.
- Analytic methods must handle an open-ended catalogue of singularities (MSMS documents
  six failure modes and restarts by perturbing radii) [MSMS docs; Quan & Stamm 2016].
  Grid methods convert that problem into a resolution/accuracy tradeoff — a tradeoff the
  user explicitly asked to control.
- Rust prior art confirms feasibility: `molframe-surface` implements flood-fill + EDT +
  6-tet polygonisation in pure Rust with rayon + `wide` (MIT). We design independently
  but validate against it.
- ChimeraX's `gridsurf.py` (~150 lines of logic) is a readable specification of a
  production-quality grid SES [Goddard 2013].

The analytic family is designed *for* (interfaces leave room, §12) and scheduled as a
later optional engine, not in v1.

---

## 3. Mathematical core

### 3.1 Definitions (recap, normative)

Atoms `(c_i, r_i)`, probe radius `r_p`. Inflated radius `R_i = r_i + r_p`.
SAS cavity `K = ⋃ B(c_i, R_i)`. Signed distance to `∂K`:

```
f_sas(p) = +dist(p, ∂K)   if p ∉ K      (probe center can be at p)
           −dist(p, ∂K)   if p ∈ K      (probe center cannot)
```

SES implicit function and zero set:

```
f_ses(p) = f_sas(p) + r_p
SES      = { p : f_ses(p) = 0 }
```

[Quan & Stamm 2016, Eq. 3.4 / 5.21]

### 3.2 Why not the naive min-field

`g(p) = min_i (‖p − c_i‖ − R_i)` has `∂K` as its zero set and *is* the exact signed
distance **outside** `K`, but inside `K` it under-reports depth (it measures distance to
the nearest inflated *sphere*, not to the boundary of the *union*). Feeding `g + r_p` to
an isosurface extractor produces a surface that is too far inward in overlap regions —
i.e. wrong reentrant geometry. Quan & Stamm exhibit analogous failures (their Fig. 2).

**Design consequence:** the grid engine must produce a true signed distance to `K`, not a
min-of-spheres field. We do this with a Euclidean distance transform (EDT) of the
voxelized solid `K` (§5.3), which is exact for the voxel set and 1-Lipschitz-consistent.

### 3.3 Error model

- `f_sas` is 1-Lipschitz (`|∇f| = 1` a.e.).
- Marching-tetrahedra/MC vertices lie on grid edges; linear interpolation of a smooth
  field puts vertices within `O(h²·κ)` of the true surface in smooth regions; mesh faces
  deviate `O(h)` worst-case; **Hausdorff distance `≤ 2h`** is the contract (N-1).
- Topology: features narrower than `≈ 2h` may be lost. Documented; user controls `h`.
- Empirical area error targets (N-2): 5% at `h=0.3`, 2% at `h=0.15` — consistent with
  EDTSurf's reported 2–4% [Xu & Zhang 2009].

### 3.4 Normals

Primary: **gradient of `f_ses`** at the vertex by central differences on the grid
(field is an SDF ⇒ gradient is the unit surface normal). Fallback: angle-weighted
pseudonormals [Bærentzen & Aanæs 2005]. Both are available; gradient is default because
it is independent of mesh tessellation and smoother across triangles.

---

## 4. System architecture

### 4.1 Workspace layout

```
molses/                          # cargo workspace (resolver = "2")
├── Cargo.toml                   # [workspace] members, [workspace.dependencies], lints
├── LICENSE-MIT
├── LICENSE-APACHE
├── README.md
├── docs/
│   ├── PRODUCT_SPEC.md
│   ├── DESIGN.md                # this document
│   └── IMPLEMENTATION_PLAN.md
├── crates/
│   ├── molses/                  # library crate (published)
│   │   ├── Cargo.toml           # features: parallel, simd, serde, mmap
│   │   ├── src/
│   │   │   ├── lib.rs           # public API re-exports
│   │   │   ├── error.rs         # SesError (thiserror)
│   │   │   ├── structure.rs     # Atom, Structure, selection filters
│   │   │   ├── parse/
│   │   │   │   ├── mod.rs       # format detection, dispatch
│   │   │   │   ├── pdb.rs       # via pdbtbx adapter
│   │   │   │   ├── mmcif.rs     # via pdbtbx adapter
│   │   │   │   ├── xyz.rs       # in-house
│   │   │   │   └── xyzr.rs      # MSMS sphere files, in-house
│   │   │   ├── radii.rs         # Bondi table + custom overrides
│   │   │   ├── field/
│   │   │   │   ├── mod.rs       # ScalarField, slab allocation
│   │   │   │   ├── paint.rs     # SIMD min-distance painting
│   │   │   │   ├── edt.rs       # separable Felzenszwalb–Huttenlocher EDT
│   │   │   │   └── celllist.rs  # uniform-grid neighbor search
│   │   │   ├── mesh/
│   │   │   │   ├── mod.rs       # Mesh, MeshBuilder
│   │   │   │   ├── polygonize.rs# marching tetrahedra (6-tet per cell)
│   │   │   │   ├── normals.rs   # gradient / angle-weighted
│   │   │   │   ├── topology.rs  # watertight checks, components, Euler
│   │   │   │   └── smooth.rs    # optional Laplacian/Taubin
│   │   │   └── io/
│   │   │       ├── mod.rs
│   │   │       ├── obj.rs       # writer (v/vn/f)
│   │   │       ├── msms.rs      # .vert/.face writer
│   │   │       └── ply.rs       # optional
│   │   ├── benches/             # criterion
│   │   └── tests/               # golden + property tests
│   └── molses-cli/              # binary crate (published as molses-cli)
│       ├── Cargo.toml           # clap, anyhow, tracing, mimalloc (opt)
│       └── src/
│           ├── main.rs
│           ├── args.rs          # clap derive
│           └── report.rs        # stats/progress output
├── testdata/                    # CC0 wwPDB fixtures + synthetic XYZ
└── .github/workflows/ci.yml
```

Published crate names: library `molses`; binary crate `molses-cli` with `[[bin]] name =
"molses"`. (Names provisional; see PRODUCT_SPEC §12.)

### 4.2 Layering (dependency direction)

```
molses-cli ──▶ molses (public API)
                   │
   ┌───────────────┼───────────────────┐
   ▼               ▼                   ▼
 parse/        field/ ─────────▶ mesh/
 radii.rs      (paint→edt)      (polygonize→normals)
   │               │                   │
   └───────────────┴─────────▶ io/ ◀───┘
```

Rules: `io` never depends on `parse`; `field` never depends on `mesh`; nothing depends on
the CLI. All cross-layer data passes through plain structs (`Structure`, `ScalarField`,
`Mesh`).

### 4.3 Core data types (public API sketch)

```rust
// structure.rs
pub struct Atom {
    pub position: [f64; 3],     // Å
    pub radius: f64,            // Å, VdW
    pub element: Option<Element>,
    pub serial: u32,
    pub chain: Option<String>,
    pub residue: Option<ResidueId>,
    pub hetero: bool,
    pub name: String,
}

pub struct Structure {
    pub atoms: Vec<Atom>,
    pub models: Vec<Range<usize>>, // atom index ranges per model
    pub source: SourceInfo,
}

// radii.rs
pub enum RadiiSet { Bondi, Custom(RadiiTable) }

// lib.rs
pub struct SesOptions {
    pub probe_radius: f64,          // default 1.4
    pub grid_spacing: f64,          // h in Å
    pub radii: RadiiSet,
    pub components: ComponentPolicy, // OuterOnly | All | NoCavities{min_volume}
    pub smoothing: u32,
    pub threads: Option<usize>,
    pub max_memory_bytes: Option<u64>,
}

pub struct Mesh {
    pub positions: Vec<[f32; 3]>,   // Å
    pub normals:   Vec<[f32; 3]>,   // unit, outward
    pub indices:   Vec<[u32; 3]>,   // CCW outward
    pub nearest_atom: Vec<u32>,     // optional, parallel to positions
    pub components: Vec<ComponentInfo>,
    pub stats: MeshStats,           // area, volume, counts
}

pub fn ses(structure: &Structure, opts: &SesOptions) -> Result<Mesh, SesError>;
```

F32 for mesh coordinates (rendering precision, half the memory of f64; PDB coordinates
are parsed as f64 and converted at the grid stage). Area/volume accumulation in f64.

---

## 5. The SES grid pipeline

The v1 engine implements the **signed-distance + isosurface** family, following the
ChimeraX/EDTSurf lineage but with a true-EDT field (not a min-of-spheres field) and
marching tetrahedra (not classic MC) for watertightness.

### 5.1 Pipeline stages

```
 Structure ──▶ [1] select/filter atoms
           ──▶ [2] radii assignment (Bondi / custom)
           ──▶ [3] cell list over inflated spheres
           ──▶ [4] domain sizing: bbox + pad = 2·r_p + r_max + 3h
           ──▶ [5] PAINT     : boolean/occupancy grid of K = ⋃B(c_i, r_i+r_p)
           ──▶ [6] EDT       : exact Euclidean distance transform → f_sas
           ──▶ [7] OFFSET    : f_ses = f_sas + r_p
           ──▶ [8] POLYGONISE: marching tetrahedra at level 0
           ──▶ [9] WELD      : deterministic vertex welding
           ──▶ [10] NORMALS  : field gradient (primary) / angle-weighted
           ──▶ [11] FILTER   : components policy, degenerate triangles
           ──▶ [12] SMOOTH   : optional Laplacian/Taubin
           ──▶ Mesh
```

### 5.2 Stage 3–4: cell list and domain

- Cell size `= max_i(r_i + r_p) + h` so any point within an inflated sphere queries the
  3×3×3 neighborhood of cells (Green's method; used by QuickSES/SESDF/RustSASA).
- Build: bucket sort, `Box<[u32]>` heads + `next[]` — O(N), ~150 lines, in-house
  (`field/celllist.rs`). No tree dependency on the hot path.
- Bounding box of all inflated spheres, expanded by `2·r_p + r_max + 3h` so the SES fits
  with a clean exterior band for the flood/EDT.
- Grid dimensions `n = (⌈extent/h⌉ + 1)` per axis; index `idx = x + y·nx + z·nx·ny`
  (x fastest).

### 5.3 Stage 5–7: field construction (the correctness-critical part)

**Goal:** a discrete signed distance field `F` whose samples approximate `f_sas`.

1. **Paint occupancy of `K`.** For every grid point `p`, `inside = ∃i : ‖p−c_i‖ ≤ R_i`.
   Implementation: iterate atoms; for each, iterate grid points in its bounding sub-box;
   test squared distance (SIMD over 4/8 points at a time, FMA). Parallel over z-slabs
   with disjoint writes (each slab owned by one thread), or per-thread local occupancy
   buffers reduced at slab boundaries.
2. **Flood-fill exterior** from the six grid faces (iterative stack, parallelizable by
   slabs) to distinguish "outside `K`" from "enclosed cavity". Voxels of `K` and enclosed
   cavities are "solid".
3. **Euclidean distance transform.** Compute exact squared EDT of the solid region to
   the nearest free voxel, using the separable Felzenszwalb–Huttenlocher algorithm
   (three 1-D passes, O(n) per pass). Parallelize the passes along the free axes
   (rows independent; use Rayon chunks with deterministic reduction).
4. **Assemble `f_sas` (f32):**
   ```
   if voxel is free (exterior):   f_sas = +min_i(‖p−c_i‖ − R_i)   // exact, from paint pass
   if voxel is solid (in K or cavity): f_sas = −EDT_distance       // true distance to ∂K
   ```
   The exterior branch is the exact analytic distance (min over inflated spheres is exact
   outside the union). The interior branch is the exact discrete distance to the
   boundary of the voxelized solid. The two agree to `O(h)` at the interface.
5. **Offset:** `F = f_sas + r_p`. Zero level of `F` is the SES.

**Precision:** `f32` fields. Coordinates are `f64` at parse time; distance math in the
paint pass in `f64` accumulators where cheap, stored as `f32`. (Matches ChimeraX
`gridsurf` f32 grids and `molframe-surface`; BEM users can request f64 vertices later.)

**Memory:** one `Box<[f32]>` for `F` + one bitset (or `Box<[u8]>)` for occupancy + EDT
scratch (one `f32` plane per pass). `= 5 bytes/voxel` ≈ 200 MB for a 350³ grid. Large
inputs use **slab chunking** (process k z-slabs at a time; polygonisation needs only a
2-slab stencil) — `max_memory_bytes` controls k; if the full grid cannot fit, `h` is
increased with a warning (PRODUCT_SPEC F-PER-4).

### 5.4 Stage 8: polygonisation — marching tetrahedra

- Split each grid cell into **6 tetrahedra** sharing the main diagonal (standard 6-tet
  decomposition). This avoids the axis-bias and ambiguous cases of classic marching
  cubes and is crack-free across cells by construction.
- 16 tet cases → 0, 1, or 2 triangles; linear interpolation of `F` on tet edges places
  vertices.
- Why not classic MC? MC's 256-case table has ambiguous faces producing topological
  cracks unless resolved consistently; sliver triangles are worse [Lorensen & Cline 1987;
  Goddard 2013]. Marching tetrahedra trades slightly more triangles for a simpler
  correctness argument and easier parallelism (independent cells).
- **Determinism:** cells processed in a fixed global order for vertex emission; parallel
  workers write into thread-local buffers keyed by cell index; final concatenation is
  sorted by cell index (stable). Output is bitwise-identical for any thread count
  (N-F-3 / F-PER-3).
- **Vertex welding:** each tet emits raw vertices on cell edges; weld by edge key
  (`(cell_axis, grid_edge_id) → vertex_index`) in a deterministic hash/BTree map during
  the merge phase. Guarantees the mesh is watertight (every interior edge shared exactly
  twice).

### 5.5 Stage 10: normals

1. Primary: `n = ∇F(p)` via central differences of `F` at the welded vertex position
   (trilinear sample of the grid gradient), then normalize. Because `F` is an SDF,
   `|∇F| ≈ 1` and the direction is the true normal up to `O(h)`.
2. Fallback / post-smoothing: angle-weighted pseudonormal of adjacent faces
   [Bærentzen & Aanæs 2005].
3. Orientation fix: ensure normals point away from the solid (outward) by checking the
   sign of `F` at `p + εn` (must be > 0). Flip triangle winding if inconsistent.

### 5.6 Stage 11: component and cavity policy

- Connected-component labeling on the triangle mesh (union-find over welded vertices).
- Default: keep only the component that contains the exterior (the outer shell).
- `--all-components`: keep all closed shells (outer + cavity walls).
- `--no-cavities`: drop components whose enclosed volume < threshold (default 0 = drop
  all inner shells).
- Report per-component area/volume in `MeshStats`.

### 5.7 Optional surfaces (F-PAR-5)

- **VdW surface:** same pipeline with `r_p = 0` (field is the signed distance to `⋃B(c_i,r_i)`).
- **SAS:** contour `f_sas` at 0 directly (skip the offset). Useful for area cross-checks
  against FreeSASA in tests.

---

## 6. Parallelism design

| Stage | Strategy | Determinism mechanism |
|---|---|---|
| Parse | per-model models in parallel (Rayon) | models concatenated in model order |
| Radii/selection | trivially parallel map | input order preserved |
| Cell list | single-threaded bucket sort (O(N), cheap) | natural |
| Paint | parallel over z-slabs; disjoint slab writes | slab ownership fixed |
| EDT | parallel 1-D passes over independent rows/columns | row-major fixed order |
| Polygonise | parallel over cell blocks; thread-local `Vec<Triangle>` | merge sorted by global cell index |
| Weld / normals | parallel over vertices after deterministic index assignment | index order fixed |
| Smoothing | parallel Jacobi-style (read old, write new) | fixed iteration count |

Additional outer parallelism: `--model all` runs one SES job per model concurrently.

**Thread control:** Rayon global pool (`RAYON_NUM_THREADS`, `--threads`). With the
`parallel` feature disabled, identical code paths use `iter()` instead of `par_iter()`
via a small internal `MaybeParallel` adapter (keeps one code shape).

**No atomics on the hot path** — disjoint slab ownership and thread-local buffers only.

---

## 7. SIMD design

`std::simd` is nightly-only as of 2026-11 [doc.rust-lang.org/std/simd]. Decision:

- Use the **`wide`** crate (Zlib/MIT/Apache) behind a `simd` cargo feature: `f32x4`/`f32x8`
  types, auto-selects SSE2/AVX/NEON. Precedent: `molframe-surface` uses `wide` with
  `#![forbid(unsafe_code)]`. Alternative: `pulp` (used by RustSASA) if runtime-dispatch
  quality becomes important — kernels stay portable.
- Scalar fallback is the same code with `f32` lanes (const-generic or macro-generated).
- Hot kernels:
  1. **Paint distance² loop** — 4/8 grid points × 1 atom, or 1 grid point × 4/8 atoms
     (choose per access pattern; painting walks grid-contiguous, so points-major).
  2. **EDT 1-D quadratic-envelope pass** — mostly scalar (branchy); vectorize the
     distance evaluation part; measure before optimizing.
  3. **Gradient/normal sampling** — trivially vectorizable.
- Keep kernels in a shape portable to `core::simd` later (element-wise lane ops only).

---

## 8. Input subsystem

### 8.1 PDB / mmCIF via `pdbtbx`

- Depend on `pdbtbx` 0.12 (MIT) behind a thin adapter module (`parse/pdb.rs`,
  `parse/mmcif.rs`) that converts `pdbtbx::PDB` into our `Structure`. Rationale:
  best-in-class pure-Rust parser (altLoc→Conformer, multi-Model, HETATM flag, element
  inference); used by RustSASA in production.
- **Risk:** `pdbtbx` maintenance is sporadic (author open to handover). Mitigation: the
  adapter is the only touchpoint; golden tests on fixtures would catch a swap; we could
  vendor a fork if abandoned.
- Do **not** expose `pdbtbx` types in the public API.

### 8.2 XYZ / XYZR in-house

- XYZ: `N`, comment line, `Element x y z` (optional extra columns ignored); multi-frame
  (repeated XYZ blocks) supported by producing multiple models. ~100 lines.
- XYZR (MSMS input): `x y z r` free-format lines, `#` comments. Direct sphere input for
  interop and unit tests.
- Format detection: extension first, then content sniffing (count line vs `ATOM`/`data_`).

### 8.3 Selection policy (defaults per PRODUCT_SPEC F-IN-6)

| Filter | Default | Flag |
|---|---|---|
| Hydrogens | excluded | `--keep-hydrogens` |
| Waters (HOH/WAT) | excluded | `--include-waters` |
| HETATM | included | `--no-hetatm` |
| Model | first | `--model N\|all` |
| altLoc | highest occupancy per residue | `--altloc highest-occupancy\|first\|all` |
| Chains | all | `--chains A,B` |

### 8.4 Radii

- Embedded `const` table: Bondi (1964) radii for the elements common in structural
  biology (H, C, N, O, S, P, Se, halogens, alkali/alkaline earth, common metals), with
  the full citation in source. Default fallback for unknown element: 1.7 Å (carbon-like)
  + warning, or hard error with `--strict`.
- `--radii custom:<file>`: CSV `element,radius` / JSON map. Overrides per element.
- Explicit collision handling: duplicate centers merged (largest radius); sphere fully
  contained in another is dropped (log count).

---

## 9. Output subsystem

### 9.1 OBJ writer (primary, F-OUT-1)

```
# molses 0.1.0
# source: 1crn.pdb
# probe_radius: 1.4  grid_spacing: 0.35  atoms: 327
# vertices: 20412  triangles: 40820  area: 4231.7 Å²  volume: 5123.4 Å³
# units: angstrom
v  x y z            # f32 printed with %.6f
vn nx ny nz         # unit
f i//i j//j k//k    # 1-based, CCW outward
```

- Deterministic vertex order = weld order = global cell order.
- Optional `g component_0` groups (F-OUT-3).
- Timestamp line only when `--stamp` or not `SOURCE_DATE_EPOCH`-style suppression
  (for byte-identical golden tests).

### 9.2 MSMS `.vert`/`.face` (F-OUT-4)

Per-vertex: `x y z nx ny nz face_number closest_sphere_index patch_type`.
Per-face: `v1 v2 v3 patch_type face_number`.
Patch type on grid meshes is **best-effort** (classify by curvature sign / nearest-atom
multiplicity); documented as approximate. Purpose: interop with VMD/Chimera pipelines
that read MSMS files.

### 9.3 PLY (optional, F-OUT-5)

ASCII PLY with `x y z nx ny nz` and optional `red green blue` per nearest-atom hash.

---

## 10. Error handling and observability

- Library: `SesError` enum (`thiserror`):
  `Io`, `Parse{path, line, msg}`, `EmptySelection`, `UnknownElement{symbol}`,
  `InvalidParameter{param, value, reason}`, `ResourceLimit{needed, limit}`,
  `Geometry{stage, msg}`. Never panics on input; `debug_assert`s allowed internally.
- CLI: `anyhow` context, maps `SesError` to exit codes (PRODUCT_SPEC §9).
- Logging: `tracing` spans per stage (`parse`, `paint`, `edt`, `polygonise`, `write`);
  CLI installs `tracing-subscriber` with `--quiet/--verbose`.
- Progress: stage timing to stderr (suppressed with `--quiet`).

---

## 11. Testing and validation strategy

| Layer | Tool | What |
|---|---|---|
| Unit | `#[test]` | EDT vs brute force on tiny grids; tet case table completeness; weld watertightness; radii table; XYZ parser edge cases |
| Property | `proptest` | normals unit length; no zero-area triangles; Euler characteristic even; area ≥ 0; convergence under refinement (coarse ⊂ fine Hausdorff bound) |
| Golden | `insta` / file hashes | OBJ hashes for 1CRN at fixed `h` and `SOURCE_DATE_EPOCH`; MSMS `.vert/.face` snapshots |
| Analytic ground truth | `#[test]` | single sphere area `4π(r+r_p)²` (rtol 1% at h=0.25); two-sphere analytic area (Richmond 1984); VdW `r_p=0` exact |
| Cross-tool | ignored-by-default tests | area vs FreeSASA (SAS only); mesh vs MSMS output when the binary is available; vs `molframe-surface` if present in the dev environment |
| Robustness | corpus test | every fixture in `testdata/` + downloaded wwPDB sample set: must succeed or return structured error; watertightness asserted |
| Bench | `criterion` | scaling vs atoms (1k/10k/100k) and vs `h` (0.75→0.15); throughput (voxels/s); parallel scaling (1..N threads) |
| Fuzz | `cargo-fuzz` (optional) | PDB/mmCIF/XYZ parsers |
| CLI | `assert_cmd` | argument matrix, exit codes, stdout OBJ validity |

CI matrix: `{ubuntu, macos, windows} × {stable} × {default features, no-default, simd-off}`.

---

## 12. Extensibility: the future analytic engine (not v1)

The public API is engine-agnostic:

```rust
pub enum Engine { Grid { spacing: f64 }, /* Exact (future) */ }
```

`SesOptions` carries `engine`; v1 implements `Grid` only. The analytic engine (reduced
surface → quadric patches → density-controlled triangulation, MSMS-style) would slot in
behind the same `ses()` entry point, providing resolution-independent areas and exact
patch metadata (PRODUCT_SPEC F-SUR-7). Interfaces (`Mesh`, `SesError`, IO) do not change.

Other planned extensions: GPU field builder (same `ScalarField` abstraction); Python
bindings over the public API; PLY/glTF writers; per-residue area reports (RustSASA-style)
using the same neighbor grid.

---

## 13. Risks and mitigations

| Risk | Impact | Likelihood | Mitigation |
|---|---|---|---|
| `pdbtbx` stagnation | parsing gap | medium | thin adapter; golden tests; vendorable fork |
| Grid memory blowup on huge assemblies at fine `h` | OOM | medium | slab chunking; `--max-memory` degrades `h` with warning; documented memory formula `≈ 5 bytes/voxel` |
| Topology loss on thin necks at coarse `h` | wrong genus | expected (documented) | document `h ≤ w/3` rule; convergence tests; `--resolution fine` |
| Classic MC watertightness traps | non-manifold mesh | low (we use 6-tet) | 6-tet decomposition + welding tests |
| Determinism vs parallelism conflicts | golden-test flakiness | medium | fixed global ordering + sorted merge; CI runs determinism test with 1 and N threads |
| `wide` / SIMD portability (NEON vs AVX) | numeric drift | low | f32 lane ops identical order; CI on aarch64 + x86-64; bitwise golden tests per-arch if needed |
| Name collision `molses` on crates.io | publish block | low | check before publish; rename is cheap pre-1.0 |
| Scope creep toward NanoShaper feature set | delay | high (social) | PRODUCT_SPEC §11 is binding for v1 |

---

## 14. Dependency manifest (v1)

| Crate | Version | License | Role | Feature-gated |
|---|---|---|---|---|
| `pdbtbx` | 0.12 | MIT | PDB/mmCIF parsing | no |
| `glam` | 0.3x | MIT OR Apache-2.0 | Vec3 at API boundary | no |
| `rayon` | 1.x | MIT OR Apache-2.0 | data parallelism | `parallel` |
| `wide` | 1.x | Zlib/MIT/Apache-2.0 | portable SIMD | `simd` |
| `thiserror` | 2.x | MIT OR Apache-2.0 | library errors | no |
| `serde` | 1.x | MIT OR Apache-2.0 | optional serialization | `serde` |
| `memmap2` | 0.9 | MIT/Apache-2.0 | huge-grid escape hatch | `mmap` |
| `clap` | 4.x (derive) | MIT OR Apache-2.0 | CLI (cli crate only) | — |
| `anyhow` | 1.x | MIT OR Apache-2.0 | CLI errors | — |
| `tracing` + `tracing-subscriber` | 0.1/0.3 | MIT | logging (cli) | — |
| `mimalloc` | 0.1 | MIT | CLI allocator (optional) | opt-in |
| `criterion` | 0.8 | MIT/Apache-2.0 | benches (dev) | — |
| `proptest` | 1.x | MIT/Apache-2.0 | property tests (dev) | — |
| `insta` | 1.x | MIT/Apache-2.0 | snapshot tests (dev) | — |
| `assert_cmd` | 2.x | MIT/Apache-2.0 | CLI tests (dev) | — |

Implemented in-house (no dependency): XYZ/XYZR parsers, Bondi radii table, cell list,
EDT, marching tetrahedra, welding, normals, OBJ/MSMS writers. Rationale in
`research/03-rust-ecosystem.md` §6.2.

MSRV: **1.85** (edition 2024), matching RustSASA's precedent. License: **MIT OR Apache-2.0**.

---

## 15. Design traceability

| PRODUCT_SPEC requirement | Design section |
|---|---|
| F-IN-1..11 input formats | §8 |
| F-PAR-1..5 geometry parameters | §5, §8.4 |
| F-SUR-1..9 surface generation | §3, §5 |
| F-OUT-1..7 output | §9 |
| F-PER-1..6 concurrency/performance | §6, §7 |
| F-API-1..5 library API | §4.2, §4.3 |
| N-1..N-3 accuracy/robustness | §3.3, §11 |
| N-4..N-6 performance targets | §11 (benches), §5.3 memory model |
| N-7 determinism | §6, §5.4 |
| LOD model §6 of spec | §5.3–5.4 |
| Mesh quality §8 of spec | §5.4–5.6, §11 |
