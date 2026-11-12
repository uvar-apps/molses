# 03 — Rust Ecosystem & Technical Building Blocks

Research stream for the SES/Connolly-surface workspace (reusable lib crate + CLI crate; PDB/XYZ/mmCIF in, OBJ triangles+normals out; rayon concurrency + SIMD). Every claim carries a source (title + URL). Uncertainty is marked explicitly. Date of survey: 2026-11 (crate versions as returned by the crates.io API at survey time).

---

## 1. Molecular structure parsing in Rust

### 1.1 pdbtbx — the dominant pure-Rust PDB/mmCIF parser

**Crate:** `pdbtbx` v0.12.0, MIT license, MSRV 1.64 (crates.io metadata) / README badge says 1.67+, ~72k downloads, last crates.io publish 2024-10-14 (docs.rs page dated 29 Oct 2026 — docs rebuild, not a new release).
- Repository: https://github.com/douweschulte/pdbtbx (74 stars, last push 2026-06-20).
- Sources: [pdbtbx on crates.io](https://crates.io/crates/pdbtbx), [pdbtbx docs.rs](https://docs.rs/pdbtbx/latest/pdbtbx/), [pdbtbx README](https://github.com/douweschulte/pdbtbx/blob/master/README.md).

**Capabilities (from README feature matrix):**
- PDB v3.30 + mmCIF v5.338 column sets. ATOM/HETATM, ANISOU, MODEL, SCALE, ORIGX, MTRIX, CRYST1 supported. SEQRES/MASTER partial. Heterogen/secondary-structure/connectivity "envisioned" (>v1.0).
- Optional features (Cargo.toml deps on docs.rs): `rayon` (parallel iterators), `rstar` (spatial lookup over atoms), `serde`, `flate2` (gzip `.pdb.gz`).
- `ReadOptions` unified read logic added in v0.12.0; `unique_conformer_names`, `chains_in_contact` helpers.

**Hierarchy & altlocs.** Structure is `PDB → Model → Chain → Residue → Conformer → Atom` (source: `src/structs/`). `Conformer` carries `alternative_location: Option<String>` (source: [conformer.rs](https://github.com/douweschulte/pdbtbx/blob/master/src/structs/conformer.rs)) — i.e. **altloc is modeled as one Conformer per altloc within a Residue**, analogous to cctbx `atom_group`. NMR multi-model ensembles are handled by multiple `Model`s. `Conformer` also has an optional `modification: Option<(String,String)>` (MODRES-ish).

**Atom fields** (source: [atom.rs](https://github.com/douweschulte/pdbtbx/blob/master/src/structs/atom.rs)):
`counter` (unique), `hetero: bool` (HETATM flag), `serial_number: usize`, `id: String` (mmCIF `_atom_site.id`), `name: String` (uppercased), `x,y,z: f64` (Å), `occupancy: f64`, `b_factor: f64`, `element: Option<Element>`, `charge: isize`, `atf: Option<[[f64;3];3]>` (anisotropic U).

**Element vs atom-name inference** (Atom::new logic, same source): 1) parse element string; 2) fall back to parsing the atom name as an element symbol; 3) fall back to first character of the name if it is one of `C H N O S`. Otherwise `element = None`. This is exactly the classic PDB element-inference heuristic — works for standard residues, will mis-handle some ligand naming (e.g. `CA` ambiguous calcium vs alpha-carbon is resolved by the explicit element column in practice).

**Insertion codes / waters / hydrogens.** Insertion codes live at the `Residue` level (pdbtbx `residue.rs` has `insertion_code`); waters are just HETATM residues with name HOH (no special-case API needed). Hydrogens are ordinary atoms — filtering them is an iterator filter on `element == Some(H)`. *(Uncertainty: I did not open `residue.rs` line-by-line to confirm the exact insertion-code getter name; the field exists per the struct listing and PDB semantics.)*

**Atomic-properties lookup.** README claims "atomic properties lookup" helpers; `src/reference_tables.rs` contains space-group tables, amino-acid lists, backbone-atom name lists — but I did **not** find a van-der-Waals radii table in pdbtbx itself. Radii assignment must come from elsewhere (see §1.3).

**Maintenance caveat (from README):** the author writes *"if anyone finds this project stale and wants to take over moderation and/or main development feel free to reach out… I would be happy to transfer the project"*. Last publish Oct 2024; development is sporadic. Plan for either (a) vendoring/forking, or (b) keeping the parsing layer thin so pdbtbx can be swapped.

**Downstream precedent:** RustSASA uses `pdbtbx` 0.12.0 as its structure parser (source: [RustSASA Cargo.toml](https://raw.githubusercontent.com/maxall41/RustSASA/main/Cargo.toml)) — validating it as the de-facto choice.

### 1.2 Gemmi and other mmCIF/PDB parsers

**Gemmi (C++) — no native Rust bindings found.**
- Gemmi is a C++ macromolecular crystallography library with Python bindings and a WASM build (`project-gemmi/wasm`). Sources: [project-gemmi/gemmi](https://github.com/project-gemmi/gemmi), [Gemmi docs](https://gemmi.readthedocs.io/en/latest/index.html).
- Searched crates.io (`q=gemmi`) and GitHub (`gh search repos "gemmi rust"`): **no `gemmi` / `gemmi-sys` crate exists**. No third-party FFI binding repo surfaced. *(Confidence: high after crates.io + GitHub search; absence-of-evidence, but two independent registries agree.)*
- Gemmi is BSD-3-Clause (upstream LICENSE), so FFI is license-compatible if ever needed, but that would add a C++/CMake build dependency — poor fit for a pure-Rust workspace goal.

**Other parsers on crates.io (survey via [crates.io API](https://crates.io/api/v1/crates?q=pdb+parser)):**

| Crate | Version | License | Notes |
|---|---|---|---|
| `patoz` | 0.2.0 | — | "Protein Data Bank (PDB) file parser and writer with typed records and a structure model" |
| `nom-pdb` | 0.0.9 | — | PDB parser via `nom`; 11k downloads |
| `isr-pdb` | 0.8.1 | — | PDB parser for ISR |
| `knuckles-parse` | 0.2.2 | — | toolkit for parsing PDB records |
| `protein` | 0.0.5 | — | "Working with protein structures" |
| `mmcif-cif` | 0.3.2 | MIT OR Apache-2.0 | Lossless CIF/mmCIF tokenizer, table layer, BinaryCIF decoder |
| `mmcif-model` | 0.3.2 | MIT OR Apache-2.0 | "Correctness-first mmCIF, BinaryCIF and PDB structure IO with typed numbering" |
| `mmcif-ccd` | 0.3.2 | MIT OR Apache-2.0 | **Compiled, versioned index of the Chemical Component Dictionary** |
| `mmcif-geom` | 0.3.2 | MIT OR Apache-2.0 | Neighbour grid queries over atomic coordinates |
| `kira-mmcif` | 0.3.0 | — | "Low-level, streaming mmCIF/BinaryCIF/PDB parser focused on protein coordinates" |
| `ustar-parser` | 0.1.4 | — | STAR format parser (CIF, mmCIF, NMR-STAR) |
| `cif` | 0.3.1 | — | CIF 1.1/2.0 lossless read/write |

The `mmcif-*` family (repo `fuyu-myk/mmcif`, MIT OR Apache-2.0) is notable: it includes a **compiled Chemical Component Dictionary index** (`mmcif-ccd`) and a neighbour-grid module (`mmcif-geom`). It's a younger, smaller-ecosystem alternative to pdbtbx. Downloads are low (85–569) so treat as immature.

**`molc` is NOT a molecular parser.** The crates.io crate `molc` v0.3.0 (105k downloads) is "A mock language client for testing language servers" — same name, different domain. *(This corrects a common confusion; if the product spec referenced "molc" as a molecular crate, that reference is wrong.)*

**`molar` v2.4.2** (Artistic-2.0, `yesint/molar`, 22.8k downloads): "rust library for analysis of MD trajectories and molecular modeling" — source: [crates.io/molar](https://crates.io/crates/molar). Relevant as a structural-analysis peer; Artistic-2.0 is GPL-compatible but not MIT/Apache, so avoid depending on it directly in a MIT/Apache-2.0 workspace without legal review.

**`molcrafts-molrs` family** (`molcrafts-molrs-io`: "File I/O for MolCrafts molecular modeling (PDB, XYZ, LAMMPS, CHGCAR, CUBE, Zarr)"): a newer molecular-modeling toolkit with format readers. Low downloads (~85–1200). Source: [docs.molcrafts.org/molrs](https://docs.molcrafts.org/molrs/reference/rust/).

**rust-bio (`bio` v4.2.2, MIT)** is sequence/alignment focused (FASTA, BWT, HMMs); its "structure" module is minimal/not a 3D molecular structure model. Not suitable as a PDB/mmCIF parser. Source: [crates.io/bio](https://crates.io/crates/bio), [rust-bio GitHub](https://github.com/rust-bio/rust-bio).

### 1.3 XYZ parsers

| Crate | Version | License | Notes |
|---|---|---|---|
| `xyz-parse` | 0.1.3 | MIT OR Apache-2.0 | "Simple parser for the XYZ file format", by aqora-io. Does **not** support extended XYZ yet. Source: [docs.rs/xyz-parse](https://docs.rs/xyz-parse/latest/xyz_parse/), [GitHub](https://github.com/aqora-io/xyz-parse) |
| `xyz_parser` | — | — | "A fast parser for XYZ molecular geometry files" (einseler/xyz_parser). *(Uncertainty: not verified on crates.io; may be git-only.)* Source: [GitHub](https://github.com/einseler/xyz_parser) |
| `xyz-convert` | 0.1.2 | — | XYZ file fixing and MD trajectory conversion utilities. Source: [docs.rs/xyz-convert](https://docs.rs/xyz-convert) |
| `extxyz-ng` | — | — | Low-level extended-XYZ parser (LAMMPS-style), pairs with `ccmat`. Source: [lib.rs/extxyz-ng](https://lib.rs/crates/extxyz-ng) |
| `sdfrust` | — | — | Pure-Rust SDF/MOL2/XYZ parser with Python bindings. Source: [GitHub](https://github.com/HFooladi/sdfrust) |
| `molcrafts-molrs-io` | 0.1.0 | — | includes XYZ reader among PDB/LAMMPS/etc. |

**Recommendation:** XYZ is a trivial format (count line, comment line, then `Element x y z` per line). Given the small surface area and the ecosystem's thin coverage, **implement XYZ parsing in-house** (~100 lines), with optional multi-frame (concatenated XYZ trajectories) and extended-XYZ as later features. Do not take a dependency for this.

### 1.4 Multi-model NMR ensembles, altlocs, hydrogens, HETATM, waters, insertion codes — handling conventions

These are the standard SES pipeline decisions. Established conventions (and how the ecosystem maps to them):

1. **Multi-model NMR ensembles** (PDB `MODEL`/`ENDMDL`): compute the surface per model (independent jobs — embarrassingly parallel across models), or accept `--model N` / `--model all`. pdbtbx exposes multiple `Model`s.
2. **Altlocs:** pick highest-occupancy altloc per residue, or altloc `A`/blank only, or expand all conformers as separate atoms. pdbtbx models altlocs as `Conformer`s so the selection is explicit. MolFrame explicitly preserves altlocs through its data model (source: [MolFrame README](https://github.com/molstation/molframe)).
3. **Hydrogens:** typically excluded from SES (united-atom radii tables assume heavy atoms); must be an explicit flag (`--no-hydrogens` default or `--keep-hydrogens`). Note: many X-ray PDB entries lack H entirely; cryo-EM/MD-derived files often have them.
4. **HETATM / ligands / ions:** included by default for a true molecular surface; `--no-hetatm` to restrict to polymer. pdbtbx `Atom::hetero()` flag makes this a filter.
5. **Waters:** usually excluded (`HOH`/`WAT` residue-name filter) — the probe radius already represents solvent. Flag to include for explicit-solvent surfaces.
6. **Insertion codes:** part of residue identity; must not be dropped when building residue-level outputs (SASA-per-residue etc.). pdbtbx Residue carries insertion codes.
7. **Element vs atom-name inference:** prefer the mmCIF `_atom_site.type_symbol` / PDB columns 77–78 element field; fall back to name heuristic (pdbtbx does this — see §1.1). For radii assignment, element is the only safe key.

### 1.5 Van der Waals radii sources

**Primary literature:** Bondi, A. (1964). *van der Waals Volumes and Radii.* J. Phys. Chem. 68(3), 441–451. DOI: [10.1021/j100785a001](https://pubs.acs.org/doi/abs/10.1021/j100785a001). 38 elements; the standard in structural biology (FreeSASA default, most SASA tools). A compiled table with citations: [KnowledgeDoor — Bondi Van der Waals Radius](https://www.knowledgedoor.com/2/elements_handbook/bondi_van_der_waals_radius.html).

**Alternates used in the wild** (source: [Wikipedia — Van der Waals radius](https://en.wikipedia.org/wiki/Van_der_Waals_radius), [mofography vdw_radii](https://mofography.github.io/mofography/vdw_radii.html)):
- **Alvarez (2013)** — updated, covers more elements (DOI 10.1039/c3dt50548e).
- **CRC / "Mendeleev"** radii.
- **Tsai–Jordan / Eisenberg** sets used by some SASA implementations (FreeSASA supports several parameter sets).

**Rust crates with element→radius data:**

| Crate | Version | License | Content |
|---|---|---|---|
| `elements_rs` | 0.2.8 | MIT | "Source-specific van der Waals radii for chemical elements" + isotopes, orbitals, oxidation states. Source: [docs.rs/elements_rs](https://docs.rs/elements_rs), [GitHub earth-metabolome-initiative/elements-rs](https://github.com/earth-metabolome-initiative/elements-rs) |
| `chemistry_consts` | 0.3.1 | MIT AND BSD-3-Clause | "Chemical reference data, physical constants, and molecular model parameters" (um-univie). Source: [crates.io/chemistry_consts](https://crates.io/crates/chemistry_consts) |
| `molframe-chem` | 0.7.0 | MIT | "Deposited-component chemistry and versioned element data for molframe"; MolFrame wheels compile in **Bondi VdW radii + covalent radii + metal sites (~3 KB)**. Source: [MolFrame README](https://github.com/molstation/molframe) |
| `mendeleev` | 0.8.4 | — | "List of chemical elements, their isotopes, and their properties" (210k downloads) |
| `periodicelements` | 0.1.2 | — | periodic info tool |

**PDB Chemical Component Dictionary (CCD):** the authoritative per-component atom naming/bonding source. `mmcif-ccd` v0.3.2 (MIT OR Apache-2.0) ships a "compiled, versioned index of the Chemical Component Dictionary". Source: [crates.io/mmcif-ccd](https://crates.io/crates/mmcif-ccd). The CCD itself is released by wwPDB under CC0 (see §4.4). Use CCD only if you need per-atom typing beyond element (e.g. distinguishing carboxylate O); for VdW radii, element-level Bondi is the scientific standard.

**Recommendation:** embed a small `const` table of Bondi (1964) radii keyed by element symbol/atomic number, with the full Bondi citation in the source comment, plus an optional override file (JSON/CSV) for custom radii. Do **not** depend on `elements_rs` etc. for 30 numbers — a table is more auditable and avoids license/attribution drift. Support `--radii bondi|alvarez|custom` in the CLI.

---

## 2. Numeric / geometry building blocks

### 2.1 Vec3 math: glam vs nalgebra vs ultraviolet

| Crate | Version | License | Downloads | SIMD story |
|---|---|---|---|---|
| `glam` | 0.34.3 | MIT OR Apache-2.0 | ~186 M | `Vec3A` = 16-byte aligned `__m128`-friendly; optional `scalar-math`, `libm`; by far the fastest in mathbench. Source: [crates.io/glam](https://crates.io/crates/glam), [mathbench-rs](https://github.com/bitshifter/mathbench-rs) |
| `nalgebra` | 0.35.0 | MIT OR Apache-2.0 | ~101 M | generic/dimensional; SIMD via optional `nalgebra-sparse`/wide integration; heavier compile times |
| `ultraviolet` | 0.10.0 | MIT OR Apache-2.0 (b/zlib-ish) | ~1.5 M | explicitly SIMD-oriented ("wide" types on stable Rust); `ultraviolet_wide` feature in mathbench. Source: [mathbench-rs](https://github.com/bitshifter/mathbench-rs) |

**mathbench-rs** (bitshifter, glam author) benchmarks all of these: glam is the performance leader for raw Vec3 ops on stable; ultraviolet is the choice when you want `Vec3x4`-style batched vectors out of the box. Source: [mathbench README](https://github.com/bitshifter/mathbench-rs).

**Recommendation:** **glam** for the public API (`Vec3A`/`Vec3`) — best performance, best ergonomics, MIT/Apache-2.0, 186M downloads means rock-solid. Keep interior grid math in raw `[f32; 3]` / `Box<[f32]>` SoA arrays and only convert at boundaries (this is what molframe-surface does — its `SurfaceTriangle.vertices` is `[[f32;3];3]`, source: [ses.rs](https://raw.githubusercontent.com/molstation/molframe/main/crates/molframe-surface/src/shape/ses.rs)). Use `wide` or `pulp` (below) for the hot loops rather than glam's SIMD.

### 2.2 Spatial acceleration

| Crate | Version | License | Notes |
|---|---|---|---|
| `kiddo` | 6.3.0 | MIT OR Apache-2.0 | k-d tree, cache-aware layouts, optional SIMD query paths; ~8.4 M downloads; 179 GitHub stars, active (pushed 2026-11-11). Sources: [docs.rs/kiddo](https://docs.rs/kiddo/latest/kiddo/), [sdd-org/kiddo](https://github.com/sdd-org/kiddo) |
| `rstar` | 0.14.0 | MIT OR Apache-2.0 | R*-tree, n-dim, pure safe Rust; ~55 M downloads; used by pdbtbx as an optional feature. Source: [lib.rs/rstar](https://lib.rs/crates/rstar) |
| `parry3d` | — | Apache-2.0 | collision-detection focused (from the `dimforge`/`rapier` ecosystem) — overkill for atom neighbor search |
| `spatial_hash` | 0.2.0 | — | 2D spatial hash only |

**For SES, the right structure is a uniform grid / cell list**, not a tree:
- SES algorithms sweep regular grids (marching cubes/tetrahedra) or need radius queries of bounded range (r + probe). A **cell list** with cell width ≈ max(radius + probe) gives O(N) build, O(1) amortized queries, and perfect cache behavior when iterating grid slabs.
- Precedent: RustSASA's JOSS paper explicitly credits "a highly optimized spatial grid" for its speed (source: [RustSASA JOSS paper](https://www.theoj.org/joss-papers/joss.09537/10.21105.joss.09537.pdf)); molframe has a dedicated `molframe-spatial` crate used by `molframe-surface` (source: [molframe-surface deps on docs.rs](https://docs.rs/molframe-surface/latest/molframe_surface/)).
- If you need general kNN (e.g. for geodesic-on-mesh work later), `kiddo` is the best-maintained k-d tree; `rstar` is what pdbtbx already integrates.

**Recommendation:** implement a small in-house cell list (~150 lines: bucket sort into `Box<[u32]>` heads + next-indices) for the atom→grid painting and for neighbor lists. Do not take a tree dependency for the hot path.

### 2.3 Marching cubes / isosurface crates — maturity assessment

| Crate | Version | Downloads | Last update | License | Assessment |
|---|---|---|---|---|---|
| `isosurface` | 0.1.0-alpha.0 | 31k | 2021-01-29 (crates.io) | Apache-2.0 | **Stale/alpha.** Marching cubes + dual contouring + docs linking papers. Nice API (`isosurface::marching_cubes`) but 5 years unmaintained. Source: [docs.rs/isosurface](https://docs.rs/isosurface/latest/isosurface/marching_cubes/), [lib.rs/isosurface](https://lib.rs/crates/isosurface) |
| `marching-cubes` | 0.1.2 | 7.6k | 2023-01-30 | — | Small, fast MC implementation. Low download count. Source: [docs.rs/marching-cubes](https://docs.rs/marching-cubes/latest/marching_cubes/) |
| `mcubes` | 0.1.7 | 11k | 2025-11-10 (repo push) | MIT | Marching cubes with Python bindings; 5 GitHub stars; by David-O'Connor (`lin-alg` Vec3 — a nonstandard math type, friction). Source: [lib.rs/mcubes](https://lib.rs/crates/mcubes), [GitHub](https://github.com/David-OConnor/mcubes) |
| `splashsurf` | 0.14.1 | 24.7k | — | — | SPH surface reconstruction (anisotropic kernel) — not a general MC, but its marching-cubes-on-grid code is a good reference |

**None of these is a compelling dependency for SES.** Reasons:
1. SES needs *signed* distance / indicator fields with exact probe-radius semantics; generic MC crates take a `Fn([f32;3]) -> f32` closure, which forces per-sample evaluation instead of a dense precomputed grid — a big perf loss.
2. Cracked/alpha status (`isosurface` 0.1.0-alpha.0) or tiny user bases.
3. We need control over triangle ordering, normals, degenerate-sliver handling, and chunked evaluation — in-house is smaller than the adapter layer.

**Better algorithmic reference:** molframe-surface polygonises the SES with a **six-tetrahedra decomposition** ("without the axis bias of exposed voxel faces") instead of vanilla MC — source: [ses.rs module doc](https://raw.githubusercontent.com/molstation/molframe/main/crates/molframe-surface/src/shape/ses.rs). Marching tetrahedra / 6-tet decomposition is crack-free, simpler to implement correctly (no 256-case table, only 16 tet cases), and produces watertight meshes. Strong candidate for our implementation.

### 2.4 Mesh types and OBJ writing

| Crate | Version | License | Notes |
|---|---|---|---|
| `tobj` | 4.0.5 | MIT | Tiny OBJ **loader** (also writes? — primarily loader); 3.3 M downloads; `single_index` flag duplicates positions for per-face normals. Source: [docs.rs/tobj](https://docs.rs/tobj), [GitHub Twinklebear/tobj](https://github.com/Twinklebear/tobj) |
| `obj-exporter` | 0.2.0 | — | Wavefront OBJ **exporter**; **last updated 2018-01-12** — abandoned. Source: [docs.rs/obj-exporter](https://docs.rs/obj-exporter/latest/obj_exporter/) |
| molframe-surface `io::write_obj` | 0.7.0 | MIT | ~20-line OBJ writer (v / vn / `f i//i`), refuses to overwrite, deterministic order. Source: [io.rs](https://raw.githubusercontent.com/molstation/molframe/main/crates/molframe-surface/src/io.rs) |

**OBJ conventions for our exporter** (sources: [Wikipedia — Wavefront .obj](https://en.wikipedia.org/wiki/Wavefront_.obj_file), [Blender OBJ manual](https://docs.blender.org/manual/en/latest/files/import_export/obj.html), molframe `write_obj`):
- Vertices `v x y z`; normals `vn nx ny nz` (unit length).
- Faces 1-indexed. Two conventions: `f v//vn` (shared index, molframe style: `f i//i`) or fully separate index sets. Per-vertex normals with `f i//i` is the simplest for isosurface meshes where each vertex gets one averaged normal.
- **Winding:** counter-clockwise when viewed from outside (right-hand rule) so the face normal `(v1-v0)×(v2-v0)` points outward from the molecular surface. Inverted winding is the #1 OBJ interchange bug (source: [StackOverflow — OBJ face vertex order](https://stackoverflow.com/questions/64242380/obj-wavefront-file-issues-due-to-face-vertex-order)).
- Many tools ignore `vn` and recompute normals (source: Blender manual above) — still emit them, since the task requires normals.
- Face normals (flat shading) vs vertex normals (smooth): for a probe-radius SES, smooth vertex normals (area-weighted average of adjacent face normals, then normalized) are standard; optionally emit face normals too.

**Recommendation:** write OBJ in-house (60 lines, deterministic, zero deps) exactly like molframe's `write_obj`, plus a tiny `Mesh { positions: Vec<[f32;3]>, normals: Vec<[f32;3]>, indices: Vec<[u32;3]> }` struct. `obj-exporter` is dead; `tobj` is a loader. Adding `serde`-driven PLY/STL later is trivial.

---

## 3. Performance tooling

### 3.1 Parallelism: rayon / crossbeam / std::thread

| Crate | Version | License | Downloads | Notes |
|---|---|---|---|---|
| `rayon` | 1.12.0 | MIT OR Apache-2.0 | ~662 M | work-stealing; MSRV 1.80. Source: [crates.io/rayon](https://crates.io/crates/rayon) |
| `crossbeam` / `crossbeam-channel` | — | MIT OR Apache-2.0 | — | channels / scoped threads; only needed for producer-consumer pipelines, not data-parallel grids |

**Embarrassingly parallel structure for SES grid work:**
1. **Atom painting** — partition atoms into chunks, paint into disjoint grid boxes or use per-thread buffers + merge. Rayon `par_iter().for_each()` over atoms works if painting is read-only-atomic or tiled by z-slabs.
2. **Distance transform** — separable EDT (Felzenszwalb–Huttenlocher) is parallel per row/column; or the "exterior flood + chamfer/EDT" pipeline molframe uses (`distance_from_exterior`). Parallelize along z: each slab's EDT depends only on local lines + a down-sweep.
3. **Polygonisation** — process independent 2×2×2 cell blocks in `par_iter`; each block emits its triangles into a thread-local `Vec`, then concatenate in deterministic order (sort by cell index) — this is what "deterministic triangle mesh in grid traversal order" means in molframe's `SolventExcludedSurface` docs.
4. **Per-model NMR ensembles** — outer-level parallelism, one model per rayon task.

**Precedents:** RustSASA deps: `rayon = "1.12.0"` (Cargo.toml above). molframe-surface deps: `rayon = "1"` (docs.rs deps list). pdbtbx optional `rayon` feature gives parallel structure iterators.

**Recommendation:** rayon only. Feature-gate it (`parallel` feature → `rayon`; serial fallback via `iter()` vs `par_iter()` is idiomatic — see rayon's own `use_rayon` pattern or the `pariter`/`maybe_rayon` crates if you want a single code path; but two code paths behind a feature is fine and zero-dep).

### 3.2 SIMD on stable Rust in 2025–2026

**`std::simd` / `core::simd` (portable SIMD) is STILL nightly-only.** The stable docs page carries the badge: *"🔬 This is a nightly-only experimental API. (portable_simd#86656)"*. Sources: [std::simd stable docs](https://doc.rust-lang.org/stable/std/simd/index.html) (verified at survey time), [rust-lang/portable-simd](https://github.com/rust-lang/portable-simd). A third-party blog claimed "stabilized in 2025" (rustz2h.com) — **that claim is not supported by the official docs**; trust the std docs. *(Confidence: high. The feature gate `portable_simd` is still listed as experimental.)*

**Stable-Rust SIMD options:**

| Approach | Crate / API | License | Notes |
|---|---|---|---|
| Portable SIMD wrapper | `wide` 1.7.1 | Zlib OR Apache-2.0 OR MIT | `f32x4`/`f32x8` types, auto-selects SSE/AVX/NEON; 95 M downloads; used by **molframe-surface** (`wide = "1.6"` in its deps). Sources: [crates.io/wide](https://crates.io/crates/wide), [pythonspeed — portable SIMD in stable Rust](https://pythonspeed.com/articles/simd-stable-rust/) |
| Safe generic SIMD | `pulp` 0.22.3 | MIT | by sarah-quinones (author of `faer`); runtime feature detection + safe generic closures; 37.8 M downloads; used by **RustSASA** (`pulp = "0.22.2"`). Sources: [docs.rs/pulp](https://docs.rs/pulp/latest/pulp/), [RustSASA Cargo.toml](https://raw.githubusercontent.com/maxall41/RustSASA/main/Cargo.toml) |
| Target intrinsics | `core::arch::{x86_64, aarch64}` | std | full control, `unsafe`, per-arch; use only after profiling |
| `safe_arch` | — | — | thin safe wrappers over `core::arch` intrinsics |
| Auto-vectorization | plain `for` loops over `&[f32]` | — | LLVM often vectorizes simple map/reduce; measure first |

**What the precedents did:**
- **RustSASA uses `pulp`** + rayon + spatial grid (Cargo.toml; JOSS paper: "readily available SIMD instructions via Pulp"). MSRV 1.85, edition 2024.
- **molframe-surface uses `wide`** + rayon (docs.rs dependency list). Also `#![forbid(unsafe_code)]` in its lib.rs — proving a safe SIMD-via-`wide` path is viable.

**Recommendation:** pick **one** of `pulp` or `wide`, behind a `simd` feature flag, with a scalar fallback. `wide` is the lower-friction choice (plain `f32x4` arithmetic, no closures) and is what a directly competing SES crate (molframe-surface) already uses. `pulp` is the choice if you want runtime-detected dispatch + `faer`-style numeric quality. Do **not** gate the whole product on nightly `std::simd` — but keep the hot kernels written so porting to `std::simd` later is mechanical (element-wise ops on fixed-width vectors).

### 3.3 Memory strategies for large grids

An SES grid at 0.5 Å resolution over a 200 Å cube is 400³ = 64 M cells; a `u8` field is 64 MB, `f32` distance is 256 MB. Strategies (sources: general practice; molframe `SurfaceGridOptions` / `max_workspace_bytes` shows the pattern of a caller-declared cell budget, [ses.rs](https://raw.githubusercontent.com/molstation/molframe/main/crates/molframe-surface/src/shape/ses.rs)):

- **`Box<[f32]>` / `Vec<f32>` SoA** — default; one flat allocation per field; index `x + y*nx + z*nx*ny` (x fastest = cache-friendly for x-sweeps).
- **Chunked slabs** — hold only k z-slabs at a time (needed for MC/tet polygonisation which reads a 2-cell-deep stencil). Emit triangles per slab; keep the total memory O(nx*ny*k).
- **`memmap2` 0.9.12** (MIT/Apache-ish, 423 M downloads; [crates.io/memmap2](https://crates.io/crates/memmap2)) — for grids that outgrow RAM; the OS pages in what the slab sweep touches. Good escape hatch for very large assemblies (ribosomes, virus capsids).
- **Cache blocking:** process tiles of ~32³ cells (128 KB for f32 field, fits L2) so the distance transform and polygonisation reuse data. Separable EDT already has good locality per axis sweep.
- **f32 vs f64:** coordinates are parsed as f64 (pdbtbx), but grid fields should be f32 (halves memory, SIMD-friendly, precision adequate for 0.1–1.0 Å resolutions). molframe-surface uses `f32` grid + `f64` area accumulation (`SurfaceTriangle.area: f64`) — a good split.
- **Allocator:** RustSASA optionally uses `mimalloc` for the CLI (`mimalloc = { version = "0.1.52", default-features = false, optional = true }`) — worth mirroring in the CLI crate only, not the lib.

### 3.4 Criterion & benchmarking

- `criterion` 0.8.2 (Apache-2.0 OR MIT, 337 M downloads, MSRV 1.86). Source: [crates.io/criterion](https://crates.io/crates/criterion). RustSASA uses criterion 0.7 (`[[bench]] name = "sasa_benchmarks" harness = false`); molframe-surface dev-deps criterion 0.5.
- Pattern: `benches/` with `harness = false`, `criterion_group!`/`criterion_main!`; benchmark scaling vs atom count and vs grid resolution; use `criterion::black_box`.
- Also `iai-callgrind` / `divan` as lighter alternatives — but criterion is the ecosystem default.

### 3.5 Workspace layout for lib + CLI

Standard, battle-tested pattern (sources: [The Cargo Book — Workspaces](https://doc.rust-lang.org/cargo/reference/workspaces.html), [The Rust Book ch14-03](https://doc.rust-lang.org/book/ch14-03-cargo-workspaces.html), [standard-rust-project-layout](https://github.com/miguelmartens/standard-rust-project-layout)):

```
ses/                          # workspace root
├── Cargo.toml                # [workspace] members, [workspace.dependencies], [workspace.lints]
├── crates/
│   ├── ses-core/             # (or ses / ses-mesh) library crate — the reusable engine
│   │   ├── Cargo.toml        # features: serde, parallel, simd
│   │   ├── src/lib.rs
│   │   ├── src/parse/        # pdbtbx glue, xyz, mmcif
│   │   ├── src/radii/        # Bondi table + overrides
│   │   ├── src/grid/         # cell list, slab allocator
│   │   ├── src/ses/          # flood, EDT, polygonise
│   │   ├── src/mesh/         # Mesh type, normals, OBJ writer
│   │   └── benches/
│   └── ses-cli/              # binary crate, thin wrapper
│       ├── Cargo.toml        # clap derive, anyhow, tracing-subscriber, mimalloc
│       └── src/main.rs
├── tests/                    # workspace-level golden tests (or per-crate)
└── testdata/                 # CC0 PDB/XYZ/mmCIF fixtures
```

RustSASA is a single-crate variant of this (`[[bin]] name = "rust-sasa"` behind `required-features = ["cli"]`) — viable if you want one crate, but the task specifies a workspace with a lib + CLI crate, and the two-crate split is cleaner for downstream reuse.

---

## 4. Project scaffolding conventions

### 4.1 Cargo workspace, features, MSRV, license

- **Workspace:** `[workspace] resolver = "2"`, `members = ["crates/*"]`, shared `[workspace.dependencies]` (pin rayon/clap/etc. once). Source: [Cargo Book — Workspaces](https://doc.rust-lang.org/cargo/reference/workspaces.html).
- **Features** (library crate): `default = []` or `default = ["parallel"]`; `parallel` → `rayon`; `simd` → `wide` (or `pulp`); `serde` → serialize meshes/atom selections; `cli` never in the lib. RustSASA's feature set is a good template (`cli = ["dep:clap", ...]`).
- **MSRV:** declare `rust-version` in every crate. Precedents: rust-sasa 1.85.0 (edition 2024), kiddo 1.89, rstar 1.90, glam 1.85, rayon 1.80, criterion 1.86. **Practical MSRV target: 1.85** (edition 2024) unless you need something newer. pdbtbx's MSRV (1.64/1.67) is lower and won't constrain you.
- **License:** task specifies dual **MIT OR Apache-2.0** — matches glam, rayon, kiddo, rstar, criterion. Note pdbtbx is **MIT-only**, rust-sasa is **MIT-only**, molar is **Artistic-2.0**. Dual-licensing your own crates is compatible with depending on MIT-only crates (MIT is permissive). Put `license = "MIT OR Apache-2.0"` + `LICENSE-MIT` + `LICENSE-APACHE` files; this is the rust-lang norm.

### 4.2 Error handling, logging, CLI parsing

| Concern | Crate | Version | License | Notes |
|---|---|---|---|---|
| Library errors | `thiserror` | 2.0.21 | MIT OR Apache-2.0 | derive `Error`; 1.8 B downloads. molframe-surface uses `thiserror = "2"` |
| App errors | `anyhow` | 1.0.104 | MIT OR Apache-2.0 | CLI-side context + `?`; 1.2 B downloads |
| (alt) library errors | `snafu` | 0.9.x | MIT OR Apache-2.0 | what RustSASA chose — either is fine |
| Logging | `tracing` | 0.1.44 | MIT | 1.0 B downloads; pair with `tracing-subscriber` in the CLI |
| CLI | `clap` | 4.6.7 | MIT OR Apache-2.0 | derive API; 1.3 B downloads; used by RustSASA (`clap = { version = "4.6.1", features = ["derive"] }`) |

Pattern: lib returns `Result<T, SesError>` (thiserror); CLI maps to `anyhow::Error`, logs via `tracing` (`info!` progress, `warn!` skipped altlocs/unknown elements).

### 4.3 Testing strategy

- **Golden files / snapshot tests for meshes:** store small expected OBJ files (or hashes of canonicalized vertex lists) under `tests/golden/`. Assert exact triangle counts + within-epsilon vertex positions. `insta` is the popular snapshot crate if you want reviewed-update workflow. Deterministic output ordering (molframe's "deterministic triangle mesh in grid traversal order") makes exact-match tests possible.
- **Property tests:** `proptest` 1.12.0 (MIT OR Apache-2.0, 242 M downloads; [crates.io/proptest](https://crates.io/crates/proptest)). Good invariants: total area ≥ 0; all normals unit-length; Euler characteristic / watertightness (every edge shared by exactly 2 triangles for a closed SES); area is continuous under small coordinate perturbations; single-atom SES area ≈ 4π(r+probe)² (analytic check).
- **Numerical regression vs FreeSASA:** RustSASA validates against FreeSASA per-atom areas (README "Validation against freesasa") — do the same for SASA cross-checks; for SES there is no FreeSASA equivalent, so validate against analytic single-sphere and two-sphere results, and against MSMS output if available.
- **Criterion benches:** scaling curves (atoms × resolution); gate CI on "no more than X% regression" via `critcmp` or stored baselines.
- **Fuzzing** (optional): `cargo-fuzz` on the parsers (PDB/mmCIF/XYZ inputs are untrusted-ish files).
- Test harness in RustSASA also uses `assert_cmd` + `escargot` for CLI integration tests — good pattern for the CLI crate.

### 4.4 Test data licensing (PDB/XYZ/mmCIF fixtures)

**PDB archive data is CC0 1.0 Universal (Public Domain Dedication).**
- Source: [wwPDB Usage Policies](https://www.wwpdb.org/about/usage-policies) — *"Data files contained in the PDB archive are available under the CC0 1.0 Universal (CC0 1.0) Public Domain Dedication. Users of PDB data are encouraged to attribute the original authors of the PDB structure data where possible."*
- [RCSB PDB Usage Policies](https://www.rcsb.org/pages/usage-policy) confirms: *"All data provided by RCSB PDB programmatic APIs are available under the same license."*
- RCSB news item confirming the CC0 switch: [PDB adopts a standard Creative Commons open source license](https://www.rcsb.org/news/feature/611e8d97ef055f03d1f222c6).

**Practical implications:**
- You may vendor PDB/mmCIF files in `testdata/` and ship them in source tarballs / Docker images freely (CC0). Keep a `testdata/README` attribiting PDB IDs + original authors anyway (scientific courtesy, and wwPDB asks for it).
- The **Chemical Component Dictionary (CCD)** is likewise wwPDB data (CC0).
- **PDB-101 educational materials** are CC BY 4.0 (different from the archive) — source: [PDB-101 How to Cite](https://pdb101.rcsb.org/more/how-to-cite). Don't confuse the two.
- **XYZ files** have no central archive; use synthetic fixtures you generate yourself, or CC0/CC-BY files from supplementary information (check each source). Safest: generate XYZ programmatically in tests.
- Avoid vendoring files whose provenance you can't trace (e.g. random GitHub test files with no license).

Good small test entries (standard choices in the field): 1CRN (crambin, 46 residues, classic), 1BBL, 4HHB (hemoglobin, multi-subunit, altlocs), 1D3Z or any NMR ensemble (multi-model), 1CBS (has waters + ligand). *(Exact IDs to confirm at implementation time; these are widely used fixtures.)*

---

## 5. Existing Rust crates for SAS/SES and molecular surfaces

### 5.1 Direct hits on crates.io (surveyed via [crates.io API](https://crates.io/api/v1/crates?q=solvent+surface) and `q=molecular+surface`, `q=sasa`)

| Crate | Version | License | What it does | SES mesh? |
|---|---|---|---|---|
| **`rust-sasa`** (RustSASA) | 0.9.2 | MIT | SASA per atom via **Shrake–Rupley**; rayon + pulp SIMD + spatial grid; CLI + Python bindings; validated vs FreeSASA. 17.7k downloads. JOSS-published (2026). Sources: [crates.io/rust-sasa](https://crates.io/crates/rust-sasa), [GitHub maxall41/RustSASA](https://github.com/maxall41/rustsasa), [JOSS paper](https://joss.theoj.org/papers/10.21105/joss.09537) | **No** — areas only |
| **`molframe-surface`** | 0.7.0 | MIT | **`solvent_excluded_surface()` → `SolventExcludedSurface { triangles, area }`**, plus Shrake–Rupley SASA, Lee–Richards SASA, cavities, pockets, depth, geodesics, curvatures, **`write_obj`**. Uses rayon + wide SIMD. Part of the MolFrame engine (Rust+Python). Sources: [docs.rs/molframe-surface](https://docs.rs/molframe-surface/latest/molframe_surface/), [ses.rs](https://raw.githubusercontent.com/molstation/molframe/main/crates/molframe-surface/src/shape/ses.rs), [GitHub molstation/molframe](https://github.com/molstation/molframe) | **YES** |
| `freesasa-sys` | 0.1.12 | non-standard (FFI wrapper; FreeSASA itself is MIT) | raw FFI to the FreeSASA C library (SASA, not SES). 20.9k downloads. Sources: [crates.io/freesasa-sys](https://crates.io/crates/freesasa-sys), [FreeSASA GitHub](https://github.com/mittinatten/freesasa) (license: MIT per GitHub API) | No |
| `voronota` / `voronota-ltr` | 0.3.0 / 0.7.1 | MIT | Voronota-LT: tessellation-derived atomic contact areas & volumes (Voronoi). 3.7k / 551 downloads. Source: [GitHub mlund/voronota-ltr](https://github.com/mlund/voronota-ltr) | No (contacts, not SES) |
| `powersasa` | 0.3.1 | — | SASA + molecular volume via power diagrams | No |
| `molframe-chem` | 0.7.0 | MIT | Bondi VdW radii, covalent radii, component chemistry (feeds molframe-surface) | — |
| `molar` | 2.4.2 | Artistic-2.0 | MD trajectory analysis / molecular modeling | No |
| `cosmolkit` | 0.6.2 | MIT | cheminformatics toolkit | No |

### 5.2 The key finding: `molframe-surface` already implements grid-based SES + OBJ

Its SES algorithm (from the module docs + source):
1. Paint probe-grown atoms (radii + probe) onto a regular grid.
2. Flood-fill the exterior of the union of probe-grown spheres from the grid boundary.
3. Compute a distance field from the exterior (separable Euclidean distance transform).
4. Polygonise the level set at distance == probe with a **six-tetrahedra decomposition** (avoids MC axis bias), producing `SurfaceTriangle { vertices: [[f32;3];3], normal: [f32;3], area: f64 }` — contact + re-entrant regions included ("rolling-probe definition").
5. `write_obj` emits `v`/`vn`/`f i//i`.

API surface (from its lib.rs re-exports): `solvent_excluded_surface(positions, radii, probe, resolution)`, `solvent_excluded_surface_with_options(..., SurfaceGridOptions)` with `max_workspace_bytes` budgets, `IndexedSurfaceMesh`, `MeshReport`, `SurfaceFace`, `write_obj`, plus `read_msms_surface` / `read_surf_surface` importers.

**Implications for us:**
- A Rust SES-with-mesh implementation **exists** (since molframe 0.7, MIT, Oct 2026). The claim "no Rust SES crate" is **false as of 2026-11**.
- However: MolFrame is pre-1.0, 0 GitHub stars at survey time (very new, `miguelcsx` / molstation), single-maintainer risk, and it's a large batteries-included engine (molframe-core/chem/spatial surface + Python story). Its SES is **grid/Eulerian** (voxel + tet polygonisation), not analytic/Connolly-triangulation (rolling-probe exact sphere–torus–sphere triangle decomposition). Grid SES has resolution-dependent area error and staircase artifacts unless resolution is fine; analytic SES (Richards/Connolly 1977; Sanner/MSMS) produces exact quadric-surface patches.
- Product positioning: if the product spec calls for **exact analytic SES triangles** (sphere faces, toroidal re-entrant faces, smooth normals from true surface gradients), molframe-surface is a *reference/validation target*, not a drop-in. If a grid SES at configurable resolution is acceptable, molframe-surface is nearly sufficient and the differentiator becomes UX/perf/robustness.

### 5.3 C/C++ SES libraries and Rust bindings

No safe Rust bindings exist for any of these (checked crates.io + GitHub). All are "call from Rust via your own FFI or subprocess" territory:

| Project | Language | Notes | Source |
|---|---|---|---|
| **FreeSASA** | C | SASA only (Shrake–Rupley, Lee–Richards), MIT; `freesasa-sys` raw FFI exists | [GitHub mittinatten/freesasa](https://github.com/mittinatten/freesasa) |
| **MSMS** (Sanner et al.) | C++ | the classic analytic SES triangulator; **non-free license** (Scripps) — cannot redistribute; executable-only distribution. Source: [Sanner, Olson & Spehner, JMS 1996] | — |
| **ESES** (WeilabMSU) | C++/grid | Eulerian SES on a regular grid, parallel; exactly the molframe approach | [GitHub WeilabMSU/ESES](https://github.com/WeilabMSU/ESES) |
| **QuickSES** | CUDA | GPU SES meshes | [GitHub nezix/QuickSES](https://github.com/nezix/QuickSES) |
| **SESDF** | CUDA | GPU large-protein SES (paper reference impl.) | [GitHub PlathC/SESDF](https://github.com/PlathC/SESDF) |
| **sesA** | — | analytic SAS + SES areas | [GitHub wlincong/sesA](https://github.com/wlincong/sesA) |
| **MolSurfComp** | Matlab | SES singularities via Voronoi-of-balls | [GitHub quanchaoyu/MolSurfComp](https://github.com/quanchaoyu/MolSurfComp) |
| **EDTSurf / NanoShaper / SES tooling** | C++ | various; not surveyed in depth | — |

### 5.4 Search summary for "does a Rust SES crate exist?"

- crates.io queries `solvent surface`, `molecular surface`, `sasa`, `connolly`, `ses`: the only mesh-producing SES is **`molframe-surface` 0.7.0**. SASA-only: `rust-sasa`, `freesasa-sys`, `powersasa`, voronota family.
- GitHub code/repo search for "connolly surface rust", "solvent excluded surface", "marching cubes molecular surface rust": produced GPU/Matlab/C++ hits (QuickSES, ESES, MolSurfComp, DeepSES) and no other Rust mesh SES.
- Conclusion: **Rust SES mesh generation is essentially an empty niche except for the very new `molframe-surface`**; analytic (quadric-patch) Connolly triangulation in Rust does **not** exist in any surveyed crate. *(Uncertainty: absence-of-evidence across crates.io + GitHub + DuckDuckGo; a small unpublished repo could be missed.)*

---

## 6. Recommendations

### 6.1 Depend on (external crates)

| Dependency | Why | License compatibility |
|---|---|---|
| `pdbtbx` 0.12 | Best-in-class pure-Rust PDB/mmCIF; altloc→Conformer, models, HETATM, element inference; optional rayon/rstar features; already proven in RustSASA. Plan a thin `Structure` trait wrapper so pdbtbx can be swapped if it goes unmaintained (author is looking for a successor). | MIT — fine under MIT OR Apache-2.0 workspace |
| `glam` | Vec3 at the API boundary; fastest stable math lib (mathbench); 186M downloads. | MIT OR Apache-2.0 |
| `rayon` | Data-parallel grid sweeps, per-model ensembles; ecosystem default. | MIT OR Apache-2.0 |
| `thiserror` (lib) + `anyhow` (cli) | Standard error split. | MIT OR Apache-2.0 |
| `clap` 4 derive (cli) | Standard CLI. | MIT OR Apache-2.0 |
| `tracing` + `tracing-subscriber` (cli) | Logging. | MIT |
| `criterion` (dev) | Benchmarks. | MIT OR Apache-2.0 |
| `proptest` (dev) | Mesh invariants (watertight, unit normals, area positivity). | MIT OR Apache-2.0 |
| `wide` (optional, `simd` feature) | Stable portable SIMD (`f32x4`/`f32x8`); molframe-surface precedent; scalar fallback without it. (`pulp` is the alternative if you prefer runtime dispatch — pick one.) | Zlib OR MIT OR Apache-2.0 |
| `memmap2` (optional) | Huge-grid escape hatch. | MIT OR Apache-2.0 |
| `mimalloc` (cli only, optional) | Allocator throughput, RustSASA precedent. | MIT |

### 6.2 Implement in-house

| Piece | Why in-house |
|---|---|
| **XYZ parser** | 100 lines; ecosystem crates are thin/inconsistent; no dependency value. |
| **Radii table (Bondi 1964 + Alvarez option + custom override)** | 30 constants + citation; more auditable than `elements_rs`; avoids pulling isotope/orbital metadata. |
| **Cell list / uniform grid** | ~150 lines; the hot-path structure; trees (kiddo/rstar) are the wrong shape for radius-bounded grid painting. |
| **SES core: probe-union flood fill + EDT + polygonisation (or analytic quadric patches)** | This *is* the product. `molframe-surface` is MIT so its algorithm is studyable; do not copy code verbatim without attribution — and its grid approach differs from an analytic Connolly triangulator if that is the spec. |
| **Marching cubes / marching tetrahedra** | No maintained crate fits (see §2.3). The 16-case tet table + interpolation is ~300 lines and gives us control over cracks, normals, determinism. |
| **OBJ writer + Mesh type** | 60 lines; `obj-exporter` is dead since 2018, `tobj` is a loader. Deterministic ordering matters for golden tests. |
| **Normal computation** | Area-weighted vertex normals + optional face normals; conventions documented in §2.4. |

### 6.3 Use as reference / validation, not dependency

- **`molframe-surface`** — read its `ses.rs`/`mesh.rs`/`io.rs` as a working algorithm description (flood + EDT + 6-tet). Validate our areas/triangles against it. Do not depend on a pre-1.0 single-maintainer engine. Its `read_msms_surface`/`read_surf_surface` show which external formats matter for interop.
- **`rust-sasa`** — validate SASA numbers (Shrake–Rupley) and mirror its Cargo feature layout / MSRV / criterion bench setup. It is MIT and stable-ish (0.9.2, JOSS-reviewed), so a SASA fast-path dependency is *defensible*, but a Shrake–Rupley implementation is ~200 lines and the product likely wants unified radii/grid code — implement SASA in-house on the same neighbor grid and cross-check against rust-sasa/FreeSASA in tests.
- **FreeSASA (C)** — golden reference for SASA numbers; call as a dev-dependency via `freesasa-sys` or as a subprocess in tests only. Do not ship FFI in the library.
- **ESES (WeilabMSU), QuickSES, MSMS** — algorithmic/paper references; MSMS license forbids redistribution.

### 6.4 Avoid

- `obj-exporter` (abandoned 2018).
- `isosurface` crate (0.1.0-alpha.0, stale 2021).
- `molar` as a direct dependency (Artistic-2.0 vs MIT/Apache-2.0 policy mismatch) unless legal review approves.
- gemmi FFI — does not exist; would add a C++ toolchain to a pure-Rust build.
- nightly `std::simd` as a hard requirement (still experimental).
- `molc` — it is a language-server test harness, not molecular.

### 6.5 Feature-flag sketch

```
ses-core features:
  default = ["parallel"]
  parallel = ["dep:rayon"]
  simd     = ["dep:wide"]
  serde    = ["dep:serde"]
  mmap     = ["dep:memmap2"]

ses-cli features (bin crate):
  default = ["parallel", "simd"]
  + clap, anyhow, tracing-subscriber, mimalloc
```

### 6.6 Open questions / uncertainty register

1. **Analytic vs grid SES** — the product spec must choose. Analytic (Connolly 1977 quadric patches / MSMS-style) gives exact smooth surfaces and resolution-independent areas but is substantially more code (toroidal re-entrant faces, singularities, external numerics). Grid + tet polygonisation (ESES / molframe-style) is simpler, parallelizes trivially, and is resolution-tunable. molframe-surface proves the grid path is tractable in Rust; the analytic path is unclaimed territory in Rust.
2. **pdbtbx longevity** — last publish 2024-10, author open to handover. Mitigate with a parser abstraction layer + golden tests that would catch a swap.
3. **`std::simd` stabilization timeline** — unstable as of 2026-11; revisit before locking SIMD kernels. Keep kernels in a form portable to `core::simd` later.
4. **molframe-surface maturity** — 0 stars / very new; do not bet the product on it, but do track it (MIT, active push 2026-11-01).
5. **CCD-based atom typing** — needed only if element-level radii prove insufficient for exotic ligands; `mmcif-ccd` is an option but low-download/immature.
6. **Insertion-code getter exact API in pdbtbx** — field presence verified at struct level, exact method name not line-verified.

---

## Bibliography

1. **pdbtbx (crates.io).** "pdbtbx — A library to open/edit/save (crystallographic) Protein Data Bank (PDB) and mmCIF files." v0.12.0, MIT. https://crates.io/crates/pdbtbx
2. **pdbtbx docs.rs.** https://docs.rs/pdbtbx/latest/pdbtbx/
3. **pdbtbx README (GitHub).** douweschulte/pdbtbx. https://github.com/douweschulte/pdbtbx/blob/master/README.md
4. **pdbtbx atom.rs source.** https://github.com/douweschulte/pdbtbx/blob/master/src/structs/atom.rs
5. **pdbtbx conformer.rs source.** https://github.com/douweschulte/pdbtbx/blob/master/src/structs/conformer.rs
6. **Gemmi — macromolecular crystallography library (GitHub).** project-gemmi/gemmi. https://github.com/project-gemmi/gemmi
7. **Gemmi documentation.** https://gemmi.readthedocs.io/en/latest/index.html
8. **RustSASA (GitHub).** maxall41/RustSASA — "A Rust library for calculating the SASA for proteins." https://github.com/maxall41/rustsasa
9. **RustSASA Cargo.toml (dependencies).** https://raw.githubusercontent.com/maxall41/RustSASA/main/Cargo.toml
10. **Campbell, M. J. (2026).** "RustSASA: A Rust Crate for Accelerated Solvent Accessible Surface Area Calculations." Journal of Open Source Software, 11(117), 9537. https://doi.org/10.21105/joss.09537 ; PDF: https://www.theoj.org/joss-papers/joss.09537/10.21105.joss.09537.pdf
11. **rust-sasa (crates.io).** v0.9.2, MIT, MSRV 1.85. https://crates.io/crates/rust-sasa
12. **rust_sasa docs.rs.** https://docs.rs/rust-sasa/latest/rust_sasa/
13. **freesasa-sys (crates.io).** raw FFI bindings for FreeSASA. https://crates.io/crates/freesasa-sys
14. **FreeSASA (GitHub).** mittinatten/freesasa — MIT. https://github.com/mittinatten/freesasa
15. **MolFrame (GitHub).** molstation/molframe — "A batteries-included structural bioinformatics engine." MIT. https://github.com/molstation/molframe
16. **molframe-surface docs.rs.** v0.7.0, MIT. https://docs.rs/molframe-surface/latest/molframe_surface/
17. **molframe-surface ses.rs source (SES algorithm).** https://raw.githubusercontent.com/molstation/molframe/main/crates/molframe-surface/src/shape/ses.rs
18. **molframe-surface io.rs source (OBJ writer).** https://raw.githubusercontent.com/molstation/molframe/main/crates/molframe-surface/src/io.rs
19. **mmcif-cif / mmcif-model / mmcif-ccd (crates.io).** fuyu-myk/mmcif family, MIT OR Apache-2.0. https://crates.io/crates/mmcif-ccd ; https://crates.io/crates/mmcif-model
20. **kira-mmcif (crates.io).** streaming mmCIF/BinaryCIF/PDB parser. https://crates.io/crates/kira-mmcif
21. **xyz-parse docs.rs.** https://docs.rs/xyz-parse/latest/xyz_parse/ ; GitHub: https://github.com/aqora-io/xyz-parse
22. **xyz_parser (GitHub).** einseler/xyz_parser. https://github.com/einseler/xyz_parser
23. **sdfrust (GitHub).** HFooladi/sdfrust — SDF/MOL2/XYZ parser. https://github.com/HFooladi/sdfrust
24. **molar (crates.io).** v2.4.2, Artistic-2.0. https://crates.io/crates/molar
25. **molcrafts-molrs docs.** https://docs.molcrafts.org/molrs/reference/rust/
26. **rust-bio (crates.io / GitHub).** https://crates.io/crates/bio ; https://github.com/rust-bio/rust-bio
27. **Bondi, A. (1964).** "van der Waals Volumes and Radii." J. Phys. Chem. 68(3), 441–451. https://pubs.acs.org/doi/abs/10.1021/j100785a001
28. **KnowledgeDoor — Bondi Van der Waals Radius table.** https://www.knowledgedoor.com/2/elements_handbook/bondi_van_der_waals_radius.html
29. **Wikipedia — Van der Waals radius.** https://en.wikipedia.org/wiki/Van_der_Waals_radius
30. **mofography — VdW Radii and Segmentation** (Bondi / CRC / Alvarez comparison). https://mofography.github.io/mofography/vdw_radii.html
31. **elements_rs docs.rs.** https://docs.rs/elements_rs ; GitHub: https://github.com/earth-metabolome-initiative/elements-rs
32. **chemistry_consts (crates.io).** MIT AND BSD-3-Clause. https://crates.io/crates/chemistry_consts
33. **glam (crates.io).** v0.34.3, MIT OR Apache-2.0. https://crates.io/crates/glam
34. **nalgebra (crates.io).** v0.35.0. https://crates.io/crates/nalgebra
35. **ultraviolet (crates.io).** v0.10.0. https://crates.io/crates/ultraviolet
36. **mathbench-rs (GitHub).** bitshifter/mathbench-rs — math library benchmarks. https://github.com/bitshifter/mathbench-rs
37. **Rust Math Libraries: nalgebra vs ndarray vs cgmath vs glam** (blog). https://www.pistack.xyz/posts/2026-07-23-rust-math-libraries-nalgebra-ndarray-cgmath-glam/
38. **kiddo docs.rs.** https://docs.rs/kiddo/latest/kiddo/ ; GitHub: https://github.com/sdd-org/kiddo
39. **rstar (lib.rs).** https://lib.rs/crates/rstar
40. **isosurface docs.rs.** https://docs.rs/isosurface/latest/isosurface/marching_cubes/ ; lib.rs: https://lib.rs/crates/isosurface
41. **marching-cubes docs.rs.** https://docs.rs/marching-cubes/latest/marching_cubes/
42. **mcubes (lib.rs / GitHub).** https://lib.rs/crates/mcubes ; https://github.com/David-OConnor/mcubes
43. **splashsurf (crates.io).** SPH surface reconstruction. https://crates.io/crates/splashsurf
44. **tobj docs.rs.** https://docs.rs/tobj ; GitHub: https://github.com/Twinklebear/tobj
45. **obj-exporter docs.rs.** https://docs.rs/obj-exporter/latest/obj_exporter/
46. **Wikipedia — Wavefront .obj file.** https://en.wikipedia.org/wiki/Wavefront_.obj_file
47. **Blender Manual — Wavefront OBJ.** https://docs.blender.org/manual/en/latest/files/import_export/obj.html
48. **StackOverflow — OBJ face vertex order / winding.** https://stackoverflow.com/questions/64242380/obj-wavefront-file-issues-due-to-face-vertex-order
49. **rayon (crates.io).** v1.12.0, MIT OR Apache-2.0. https://crates.io/crates/rayon
50. **wide (crates.io).** v1.7.1, Zlib OR Apache-2.0 OR MIT. https://crates.io/crates/wide
51. **pulp docs.rs.** v0.22.3, MIT. https://docs.rs/pulp/latest/pulp/
52. **std::simd (stable docs — nightly-only badge).** https://doc.rust-lang.org/stable/std/simd/index.html
53. **rust-lang/portable-simd (GitHub).** https://github.com/rust-lang/portable-simd
54. **pythonspeed — Using portable SIMD in stable Rust.** https://pythonspeed.com/articles/simd-stable-rust/
55. **memmap2 (crates.io).** v0.9.12. https://crates.io/crates/memmap2
56. **criterion (crates.io).** v0.8.2, Apache-2.0 OR MIT. https://crates.io/crates/criterion
57. **thiserror (crates.io).** v2.0.21. https://crates.io/crates/thiserror
58. **anyhow (crates.io).** v1.0.104. https://crates.io/crates/anyhow
59. **tracing (crates.io).** v0.1.44. https://crates.io/crates/tracing
60. **clap (crates.io).** v4.6.7. https://crates.io/crates/clap
61. **proptest (crates.io).** v1.12.0. https://crates.io/crates/proptest
62. **The Cargo Book — Workspaces.** https://doc.rust-lang.org/cargo/reference/workspaces.html
63. **The Rust Book — Cargo Workspaces.** https://doc.rust-lang.org/book/ch14-03-cargo-workspaces.html
64. **standard-rust-project-layout (GitHub).** miguelmartens/standard-rust-project-layout. https://github.com/miguelmartens/standard-rust-project-layout
65. **wwPDB Usage Policies** (CC0 1.0 for PDB archive data). https://www.wwpdb.org/about/usage-policies
66. **RCSB PDB Usage Policies.** https://www.rcsb.org/pages/usage-policy
67. **RCSB news — "PDB adopts a standard Creative Commons open source license."** https://www.rcsb.org/news/feature/611e8d97ef055f03d1f222c6
68. **PDB-101 How to Cite** (CC BY 4.0 for PDB-101 materials). https://pdb101.rcsb.org/more/how-to-cite
69. **ESES (GitHub).** WeilabMSU/ESES — Eulerian SES. https://github.com/WeilabMSU/ESES
70. **QuickSES (GitHub).** nezix/QuickSES — GPU SES meshes (CUDA). https://github.com/nezix/QuickSES
71. **SESDF (GitHub).** PlathC/SESDF — GPU SES reference implementation. https://github.com/PlathC/SESDF
72. **sesA (GitHub).** wlincong/sesA — analytic SAS/SES areas. https://github.com/wlincong/sesA
73. **MolSurfComp (GitHub).** quanchaoyu/MolSurfComp — Matlab SES singularities. https://github.com/quanchaoyu/MolSurfComp
74. **voronota-ltr (crates.io / GitHub).** https://crates.io/crates/voronota-ltr ; https://github.com/mlund/voronota-ltr
75. **cosmolkit (crates.io).** https://crates.io/crates/cosmolkit
76. **crates.io API** (crate metadata queries used throughout): https://crates.io/api/v1/crates
