# Implementation Plan — `molses`

**Work breakdown, milestones, and acceptance criteria**

| Field | Value |
|---|---|
| Status | Draft for review (pre-implementation) |
| Date | 2026-11-12 |
| Companion docs | `PRODUCT_SPEC.md` (requirements), `DESIGN.md` (technical design) |
| Scope of this document | Plan only — **no code is written yet** |

---

## 0. How to read this plan

- Work is organized into **milestones M0–M7**, each with tasks, deliverables, and
  **acceptance criteria** (testable statements). A milestone is done only when all its
  acceptance criteria pass.
- Task IDs (`T-M1-3`) are stable references for commit messages, PRs, and the task
  tracker.
- Estimates are in **ideal engineering days** for one engineer familiar with Rust and
  computational geometry; they are planning aids, not commitments.
- Ordering is strictly sequential across milestones (each milestone builds on the
  previous), except M6/M7 which are post-1.0 and independent.

Total estimate to v1.0: **~30–40 engineering days**.

```
M0 scaffold (2d) ─▶ M1 input (4d) ─▶ M2 field (6d) ─▶ M3 mesh (5d)
                                                         │
                       M4 IO + CLI (3d) ◀────────────────┘
                              │
                       M5 perf + validation + release (6–10d)
                              │
                       M6 exact engine / M7 extensions (post-1.0)
```

---

## 1. Milestone overview

| ID | Name | Goal | Est. | Depends on |
|---|---|---|---|---|
| M0 | Workspace scaffold | Compiling workspace, CI, lint, empty API surface | 2 d | — |
| M1 | Structure input & radii | PDB/mmCIF/XYZ/XYZR → `Structure` with radii and selection | 4 d | M0 |
| M2 | Distance field engine | Cell list, paint, EDT, `f_ses` grid, memory model | 6 d | M1 |
| M3 | Mesh extraction | Marching tetrahedra, welding, normals, topology | 5 d | M2 |
| M4 | Output + CLI | OBJ/MSMS writers, full CLI, end-to-end runs | 3 d | M3 |
| M5 | Performance, validation, release | SIMD, parallel tuning, golden corpus, benches, 0.1.0→1.0 | 6–10 d | M4 |
| M6 | Exact analytic engine (post-1.0) | Reduced surface + quadric patches | 20+ d | M5 |
| M7 | Extensions (post-1.0) | Python bindings, GPU, PLY/glTF, areas API | TBD | M5 |

---

## 2. M0 — Workspace scaffold (est. 2 days)

### Tasks

| ID | Task | Notes |
|---|---|---|
| T-M0-1 | Create cargo workspace: root `Cargo.toml`, `crates/molses`, `crates/molses-cli`, `LICENSE-MIT`, `LICENSE-APACHE`, `README.md` skeleton | Workspace resolver 2; `rust-version = "1.85"`; edition 2024 |
| T-M0-2 | Define feature flags (`parallel`, `simd`, `serde`, `mmap`) with empty implementations and docs | `cargo build --no-default-features` must pass |
| T-M0-3 | Skeleton modules per DESIGN §4.1 (`error`, `structure`, `parse`, `radii`, `field`, `mesh`, `io`) with `todo!()`-free stubs that compile | Public API types declared (§4.3 of DESIGN) but returning `unimplemented` errors is acceptable at this stage |
| T-M0-4 | CI workflow: build + test + clippy + fmt on {ubuntu, macos, windows} × {default, no-default-features} | GitHub Actions |
| T-M0-5 | `clippy` pedantic-ish lint config (`[workspace.lints]`), `rustfmt.toml`, `deny.toml` (cargo-deny for licenses) | License policy: MIT/Apache/BSD/Zlib only |
| T-M0-6 | Testdata directory + `testdata/README.md` (CC0 wwPDB attribution policy) | Download 1CRN, 4HHB (PDB + mmCIF), synthetic XYZ |
| T-M0-7 | Decision record: crate names confirmed/renamed after crates.io name check | PRODUCT_SPEC §12 item 3 |

### Acceptance criteria

- [ ] `cargo build`, `cargo test`, `cargo clippy -- -D warnings`, `cargo fmt --check` pass in CI on all three OSes.
- [ ] `cargo doc --no-deps` builds with zero warnings.
- [ ] `cargo deny check` passes (no copyleft/unknown licenses).
- [ ] Testdata present: ≥2 PDB, ≥2 mmCIF (same entries), ≥3 XYZ (including a multi-frame one), ≥1 XYZR.
- [ ] Public API compiles as documented in DESIGN §4.3 (types and function signatures stable enough to review).

---

## 3. M1 — Structure input & radii (est. 4 days)

### Tasks

| ID | Task | Notes |
|---|---|---|
| T-M1-1 | `Atom`, `Structure`, `ResidueId`, `SourceInfo` implemented; model ranges; equality/identity semantics | No pdbtbx types leak through |
| T-M1-2 | `parse/pdb.rs` + `parse/mmcif.rs` via `pdbtbx` adapter | Handle MODEL/ENDMDL, altLoc (Conformer), HETATM flag, insertion codes, element column + name fallback |
| T-M1-3 | `parse/xyz.rs` in-house: single + multi-frame XYZ | Element lookup → radii later |
| T-M1-4 | `parse/xyzr.rs` in-house: MSMS sphere input | `#` comments, free format |
| T-M1-5 | Format detection (extension + content sniffing), gzip via `flate2` (optional feature `gzip`) | `.pdb.gz`, `.cif.gz` |
| T-M1-6 | `radii.rs`: Bondi 1964 `const` table with citation comment; `RadiiSet::Custom` CSV/JSON loader; unknown-element policy (warn+default 1.7 Å, or `--strict` error) | |
| T-M1-7 | Selection filters: hydrogens, waters, HETATM, chains, model, altLoc policy (highest-occupancy default) | Exact per PRODUCT_SPEC F-IN-6 |
| T-M1-8 | Degenerate-atom handling: duplicate centers merged (largest radius), fully nested spheres dropped, with counters in a `ParseReport` | F-IN-8 |
| T-M1-9 | Unit tests: fixture-based golden atom counts; edge cases (empty file, no ATOM records, unknown element, insertion codes, altLoc occupancy ties) | |
| T-M1-10 | `SesError::Parse` paths covered by tests (malformed columns, bad XYZ counts) | |

### Acceptance criteria

- [ ] 1CRN PDB and mmCIF produce identical `Structure` (atom count, coordinates within 1e-6, elements).
- [ ] 4HHB altLoc policy: highest-occupancy selection produces the documented atom count.
- [ ] NMR ensemble fixture: `--model all` returns `n` models with correct atom ranges.
- [ ] XYZ and XYZR fixtures parse; multi-frame XYZ yields multiple models.
- [ ] Unknown element with `custom` radii file resolves; without it: warning + fallback (default) or error (`--strict`).
- [ ] Zero panics under `cargo fuzz` smoke run of parsers (≥1 h budget, optional in CI).
- [ ] All M1 unit tests pass on all three OSes.

---

## 4. M2 — Distance field engine (est. 6 days)

*The correctness-critical milestone. Do not proceed to M3 until analytic ground truths pass.*

### Tasks

| ID | Task | Notes |
|---|---|---|
| T-M2-1 | `field/celllist.rs`: uniform grid, bucket-sort build, 27-cell radius queries | ~150 lines; unit tests vs brute force |
| T-M2-2 | Domain sizing: bbox of inflated spheres + pad `2·r_p + r_max + 3h`; grid dims; index math | Tests for index round-trip |
| T-M2-3 | `field/paint.rs`: occupancy of `K = ⋃B(c_i, r_i+r_p)`; also record exact exterior distance `min_i(‖p−c_i‖ − R_i)` during the same pass | f32 storage; SIMD-ready inner loop shape (scalar first, SIMD in M5) |
| T-M2-4 | Exterior flood-fill from grid faces; mark solid = `K` ∪ enclosed cavities | Stack-based; slab-parallel later |
| T-M2-5 | `field/edt.rs`: separable Felzenszwalb–Huttenlocher exact EDT of solid→free distance (squared domain) | Validate vs brute-force EDT on random 16³ grids |
| T-M2-6 | Assemble `f_sas` (exterior = exact min-distance, interior = −EDT), then `f_ses = f_sas + r_p` | The DESIGN §5.3 contract |
| T-M2-7 | Memory model: full-grid allocation + slab-chunked variant; `max_memory_bytes` budget logic (estimate, warn, raise `h`) | Document formula in README |
| T-M2-8 | `ScalarField` accessor API: trilinear sample, gradient sample (for M3 normals) | |
| T-M2-9 | Analytic ground-truth tests (single atom; two atoms; probe 0 = VdW): field values at probe points match closed-form distances (rtol 1e-4 f32) | **Gate for M3** |
| T-M2-10 | Convergence test: field at a fixed set of sample points converges as `h` halves (error ratio ≈ 2 for first order) | |
| T-M2-11 | Parallel paint + EDT (Rayon, slab/row decomposition) with bit-identical results vs serial | Determinism gate |

### Acceptance criteria

- [ ] EDT matches brute force on 100 random grids (max abs error 0 in squared units).
- [ ] Single-atom field: `|f_sas(p) − (‖p−c‖ − R)| < 1e-4` at 1000 random exterior points; interior distances match analytic.
- [ ] Two-atom field matches analytic union-of-two-balls distance at 1000 random points (exterior exact; interior vs brute-force ray sampling within `O(h)`).
- [ ] Flood-fill correctly marks enclosed cavities on a synthetic "hollow shell" fixture.
- [ ] Serial and 8-thread runs produce **bit-identical** fields.
- [ ] Memory estimate function predicts actual allocation within 5%.
- [ ] 10 000 synthetic random atoms at `h = 0.5`: field built in < 1 s (debug-mode warning acceptable; target in release).

---

## 5. M3 — Mesh extraction (est. 5 days)

### Tasks

| ID | Task | Notes |
|---|---|---|
| T-M3-1 | `mesh/polygonize.rs`: 6-tet decomposition per cell; 16-case table; edge interpolation of `F` | No cracks: shared-edge vertex keys consistent |
| T-M3-2 | Deterministic welding: edge-key → vertex index; thread-local emit + sorted merge by global cell index | Bit-identical across thread counts |
| T-M3-3 | `mesh/normals.rs`: gradient-of-field normals (trilinear gradient sample, normalize); orientation fix (outward via `F` sign probe); angle-weighted fallback | |
| T-M3-4 | `mesh/topology.rs`: watertight check (edge-use counts), manifold check, connected components (union-find), Euler characteristic, genus report | Used by tests and `--verbose` stats |
| T-M3-5 | Component policy: outer only / all / volume-threshold cavity drop | F-SUR-5 |
| T-M3-6 | Degenerate-triangle removal (zero area); sliver count reported | |
| T-M3-7 | `mesh/smooth.rs`: Laplacian and Taubin smoothing, fixed iterations, deterministic | F-SUR-8, off by default |
| T-M3-8 | `Mesh`/`MeshStats`: area (f64 sum), signed volume (divergence theorem), counts, per-component stats | F-SUR-9 |
| T-M3-9 | `nearest_atom` per vertex (cell list query) | F-SUR-6 |
| T-M3-10 | Property tests: watertight on all fixtures; unit normals; outward orientation; area monotone under refinement | |
| T-M3-11 | Analytic mesh tests: single sphere area within 1% at `h = 0.25`; two-sphere area vs Richmond/Connolly analytic within 5% at `h = 0.3` | **Gate for M4** |

### Acceptance criteria

- [ ] Single-atom SES mesh: watertight, genus 0, area error < 1% at `h = 0.25`, volume error < 1%.
- [ ] Two-atom mesh: watertight; area within 5% of analytic at `h = 0.3`.
- [ ] 1CRN at `h = 0.35`: watertight outer component; Euler characteristic even; no zero-area triangles.
- [ ] All component policies behave as specified on a hollow-shell synthetic fixture.
- [ ] 1-thread vs 8-thread meshes are byte-identical (after canonical sorting).
- [ ] Normals: `|‖n‖ − 1| < 1e-5`; outward orientation check passes on all fixtures.
- [ ] `nearest_atom` indices within radius+probe of each vertex.

---

## 6. M4 — Output + CLI (est. 3 days)

### Tasks

| ID | Task | Notes |
|---|---|---|
| T-M4-1 | `io/obj.rs`: `v`/`vn`/`f i//i` writer, header comments, `g` groups, stdout/file, `--force` | F-OUT-1..3, F-OUT-6 |
| T-M4-2 | `io/msms.rs`: `.vert`/`.face` writer (best-effort patch type) | F-OUT-4 |
| T-M4-3 | Optional `io/ply.rs` | F-OUT-5 |
| T-M4-4 | `molses-cli`: clap args exactly as PRODUCT_SPEC §9; exit codes; `--quiet/--verbose`; thread control | The CLI must use only public lib APIs (F-API-5) |
| T-M4-5 | Resolution presets (`draft`…`ultra`) + explicit `h`; `--max-triangles` auto-`h` via coarse-area estimate | F-PAR-2 |
| T-M4-6 | Progress/timing report (`report.rs`): atoms, grid dims, memory, stage times, area/volume/counts | |
| T-M4-7 | End-to-end CLI tests (`assert_cmd`): golden OBJ for 1CRN (`SOURCE_DATE_EPOCH` fixed); stdout mode; error exit codes | |
| T-M4-8 | README quickstart + CLI help text polish; example OBJ screenshots (optional) | |

### Acceptance criteria

- [ ] `molses testdata/1crn.pdb -o out.obj` produces a mesh openable in Blender/MeshLab (manual check) and passing `topology.rs` watertightness.
- [ ] `molses --resolution ultra` measurably increases triangle count vs `draft`; `--max-triangles 50000` respects the budget within 10%.
- [ ] Golden OBJ test: byte-identical across runs and across 1 vs 8 threads.
- [ ] MSMS-format output loads in at least one external tool (manual: VMD or Chimera if available; else format-validator test).
- [ ] All documented flags exist and are covered by at least one CLI test.
- [ ] `cargo doc` examples on `ses()` compile (doctests).

---

## 7. M5 — Performance, validation, release (est. 6–10 days)

### Tasks

| ID | Task | Notes |
|---|---|---|
| T-M5-1 | SIMD painting kernel (`wide`, `simd` feature) + scalar fallback parity tests | Bitwise-close (document tolerance) vs scalar |
| T-M5-2 | Parallel scaling pass: profile (flamegraph/perf), fix contention, tune slab sizes | Target N-4/N-5 |
| T-M5-3 | Criterion benches: atoms ×1k/10k/100k; `h` 0.75→0.15; threads 1..N; record baselines | |
| T-M5-4 | Validation suite vs external references: MSMS areas on 1CRN/1TLP (if binary available; else recorded reference numbers from Sanner 1996); FreeSASA SAS cross-check; optional `molframe-surface` comparison | Document results in `docs/VALIDATION.md` |
| T-M5-5 | Robustness corpus: ≥50 wwPDB entries spanning size/classes; must not panic; watertight or structured error | Nightly/weekly CI job |
| T-M5-6 | Memory stress: 500k-atom synthetic capsid at `h = 1.0` with slab chunking | N-6 |
| T-M5-7 | Documentation pass: README, crate docs, `docs/LOD.md` (resolution guide), CHANGELOG, CONTRIBUTING | |
| T-M5-8 | API freeze review: semver-surface audit (public items), rustdoc lint clean | |
| T-M5-9 | Publish 0.1.0 to crates.io (library + CLI); git tag; GitHub release notes | After name check (T-M0-7) |
| T-M5-10 | 1.0 criteria: all PRODUCT_SPEC NFRs met or explicitly waived in CHANGELOG; validation doc reviewed | |

### Acceptance criteria (release gates)

- [ ] **N-1** Hausdorff ≤ 2h: measured on 1CRN/4HHB vs high-res reference mesh.
- [ ] **N-2** area error ≤ 5% @ h=0.3, ≤ 2% @ h=0.15 vs MSMS/analytic.
- [ ] **N-3** zero panics on 50-entry corpus.
- [ ] **N-4** 100k atoms @ h=0.5 ≤ 5 s on reference 8-core machine (document the machine).
- [ ] **N-5** 10k atoms @ h=0.5 ≤ 0.5 s.
- [ ] **N-6** 500k atoms @ h=1.0 completes under slab mode without OOM on 16 GB machine.
- [ ] **N-7** determinism test in CI.
- [ ] **N-8** CI green on linux/macos/windows, x86-64 (aarch64 best-effort).
- [ ] **N-9** coverage ≥ 80% on `field/`, `mesh/`, `parse/`.
- [ ] **N-10** `cargo deny` clean; testdata CC0 documented.
- [ ] `docs/VALIDATION.md` exists with tables (area/volume/timing vs h, vs references).

---

## 8. M6 — Exact analytic engine (post-1.0, est. 20+ days)

*Optional second engine behind `Engine::Exact` — not required for 1.0.*

| ID | Task |
|---|---|
| T-M6-1 | Reduced surface (RS) construction from weighted alpha-shape / probe-contact theory [Sanner 1996; Edelsbrunner 1995; Liang 1998] |
| T-M6-2 | Analytic patch extraction: contact spherical, toric, spheric reentrant; singularity handling per Quan & Stamm 2016 §5.2 |
| T-M6-3 | Density-controlled triangulation (vertices/Å², MSMS-compatible parameter) |
| T-M6-4 | Exact areas/volumes; patch-type metadata (F-SUR-7 full fidelity) |
| T-M6-5 | Cross-validation vs MSMS `.vert/.face` on the corpus; convergence vs grid engine |

Acceptance: patch-complete SES on 1CRN matching MSMS area to < 0.5%; documented
singularity policy; no worse robustness than MSMS on the corpus.

---

## 9. M7 — Extensions (post-1.0, TBD)

| Item | Trigger / notes |
|---|---|
| Python bindings (PyO3/maturin) | User demand for notebook/MDAnalysis interop |
| GPU field builder (wgpu/compute or CUDA feature) | If N-4 targets are insufficient for interactive use |
| PLY/glTF writers; per-residue area API | Consumers request |
| Trajectory support (XTC/DCD) | MD users |
| Mesh remeshing for FEM quality | BEM users (Khan et al. 2018 quality targets) |

---

## 10. Cross-cutting work (every milestone)

| Activity | Standard |
|---|---|
| Git hygiene | Small PRs, one task ID per PR, conventional commits (`feat(m2): ...`) |
| Tests | Every task lands with its tests in the same PR |
| Docs | Public items documented before merge; examples as doctests where practical |
| Determinism | Any parallel code ships with the 1-vs-N-thread identity test |
| Benchmarks | Perf-affecting PRs update criterion baselines in the PR description |
| Review | Geometry/algorithm PRs (M2, M3) get extra scrutiny against `DESIGN.md` §3 |

---

## 11. Risk register (plan-level)

| Risk | When it bites | Contingency |
|---|---|---|
| Two-atom analytic test reveals EDT assembly bug | M2 | Block M3; re-derive `f_sas` assembly; consult Quan & Stamm §4 |
| Marching tetrahedra vertex welding leaks non-manifold edges | M3 | Fall back to per-cell vertex duplication + post-weld with exact edge keys; add fuzz tests on random fields |
| Parallel determinism proves harder than expected (e.g. hash map ordering) | M2/M3 | Use `BTreeMap`/sorted vecs in merge paths; keep `HashMap` only where order is later sorted |
| 100k-atom perf target missed | M5 | Profile first; likely wins: paint SIMD, slab size tuning, `mimalloc` in CLI; last resort raise default `h` |
| `pdbtbx` breaking change / abandonment | M1/M5 | Adapter isolation; vendor fork under `crates/vendor/pdbtbx` if needed |
| crates.io name taken | M0 | Rename before first publish (cheap) |
| Validation vs MSMS impossible (binary not runnable/licensed) | M5 | Use published Sanner-1996 numbers for 1CRN/1TLP + analytic ground truths; cite limits in `VALIDATION.md` |

---

## 12. Definition of Done (v1.0)

v1.0 is done when:

1. All M0–M5 acceptance criteria are checked.
2. `PRODUCT_SPEC.md` §5 NFRs are met or explicitly waived with rationale in `CHANGELOG.md`.
3. `docs/VALIDATION.md` publishes accuracy and performance tables.
4. The library is on crates.io with docs.rs building clean; the CLI is installable via
   `cargo install molses-cli`.
5. A new user can go from `cargo add molses` (or `cargo install molses-cli`) to a viewable
   OBJ in under five minutes, following the README only.

---

## 13. Immediate next actions (when implementation is approved)

1. Confirm/rename crates (`molses` name check on crates.io) — T-M0-7.
2. Create the workspace scaffold and CI (M0, ~2 days).
3. Fetch and license-check the testdata fixtures (T-M0-6).
4. Start M1 (parsers) — the only milestone with no geometry risk.

No code has been written. Implementation begins only on explicit approval of this plan.
