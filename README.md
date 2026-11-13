# molses

**Solvent Excluded Surface (Connolly surface) generation for molecular structures — Rust library + CLI.**

`molses` computes a resolution-controllable, watertight triangle mesh of the
Solvent Excluded Surface (SES) of a molecule from PDB, mmCIF/PDBx, or XYZ input,
and exports it as a Wavefront OBJ file with per-vertex normals.

The SES (also called the Connolly surface or molecular surface) is the boundary
of the volume from which spherical solvent molecules are excluded — the surface
traced by a solvent probe of radius `r_p` (default 1.4 Å for water) as it rolls
over the van der Waals spheres of the atoms [Richards 1977; Connolly 1983;
Sanner et al. 1996].

## Table of contents

- [Installation](#installation)
- [Command-line usage](#command-line-usage)
- [Library usage](#library-usage)
- [Mathematical definition](#mathematical-definition)
- [Accuracy and level of detail](#accuracy-and-level-of-detail)
- [Tools](#tools)
- [Feature flags](#feature-flags)
- [License](#license)

---

## Installation

### CLI

```console
$ cargo install molses-cli
```

Or build from source:

```console
$ git clone https://github.com/example/molses.git
$ cd molses
$ cargo install --path crates/molses-cli
```

### Library

Add to your `Cargo.toml`:

```toml
[dependencies]
molses = "0.1"
```

---

## Command-line usage

### Quick start

```console
# Generate an OBJ file from a PDB structure
$ molses 1crn.pdb -o 1crn.obj

# Higher resolution
$ molses 1crn.pdb -o 1crn.obj --resolution fine

# Custom probe radius and output to stdout
$ molses structure.cif --probe-radius 1.5 | head

# Explicit grid spacing in Å (overrides --resolution)
$ molses 1crn.pdb -o out.obj --grid-spacing 0.1
```

### Input formats

| Format | Extensions | Notes |
|---|---|---|
| PDB | `.pdb`, `.ent`, `.brk` | Fixed-column, ATOM/HETATM, multi-MODEL |
| mmCIF/PDBx | `.cif`, `.mmcif` | wwPDB successor format |
| XYZ | `.xyz` | `Element x y z` per line |
| XYZR | `.xyzr` | `x y z r` per line (MSMS interchange) |

Gzip-compressed inputs (`.gz`) are supported with the `gzip` feature.

### Full options

```console
$ molses --help
```

```
Usage: molses [OPTIONS] <INPUT>

Arguments:
  <INPUT>                  Structure file: .pdb / .ent / .cif / .mmcif / .xyz / .xyzr

Options:
  -o, --output <FILE>      Output path (stdout if omitted)
      --probe-radius <R>   Probe radius in Å [default: 1.4]
      --resolution <PRESET>
                           Grid spacing preset: draft|coarse|balanced|fine|ultra
                           [default: balanced]
      --grid-spacing <H>   Explicit grid spacing in Å (overrides --resolution)
      --radii <SET>        van der Waals radii: bondi | custom:<file> [default: bondi]
      --model <SEL>        Model selection: first | all | <N> (1-based) [default: first]
      --chains <LIST>      Comma-separated chain IDs to include
      --keep-hydrogens     Include hydrogen atoms [default: excluded]
      --no-hetatm          Exclude HETATM records
      --include-waters     Include waters (HOH/WAT) [default: excluded]
      --altloc <POLICY>    Alternate-location policy: highest-occupancy | first | all
      --surface <TYPE>     ses | sas | vdw [default: ses]
      --all-components     Include internal cavity surfaces
      --smooth <N>         Taubin smoothing iterations [default: 0]
      --threads <N>        Worker threads [default: all]
      --max-memory <MB>    Soft memory budget in MB
      --format <FMT>       Output format: obj | msms [default: obj]
      --force              Overwrite existing output
      --quiet              Suppress progress output
      --verbose            Detailed progress
  -h, --help               Print help
  -V, --version            Print version
```

### Resolution presets

The primary level-of-detail control is the grid spacing `h` (Å). Smaller values
give smoother surfaces with more triangles.

| Preset | `h` (Å) | Intended use | Approx. triangles for a 25 000 Å² protein |
|---|---|---|---|
| `draft` | 0.75 | Interactive preview | ~50 000 |
| `coarse` | 0.50 | Quick visualization | ~100 000 |
| `balanced` | 0.35 | **Default** — good quality | ~200 000 |
| `fine` | 0.25 | Publication figures | ~400 000 |
| `ultra` | 0.15 | Near-analytic geometry | ~900 000 |

For anything finer than `ultra`, use `--grid-spacing`:

```console
$ molses 1crn.pdb -o out.obj --grid-spacing 0.1   # ~2M triangles
```

### Examples

**Standard workflow — PDB to OBJ:**

```console
$ molses 1crn.pdb -o 1crn.obj
```

Output:
```
# molses 0.2.0
# source: 1crn.pdb
# probe_radius: 1.4  grid_spacing: 0.35  atoms: 327
# vertices: 81064  triangles: 162124  area: 2264.453  volume: 4962.740
# units: angstrom
v 3.007255 6.108255 -8.667745
...
vn 0.273548 0.371933 0.887039
...
f 1//1 2//2 3//3
...
```

**mmCIF input:**

```console
$ molses 1crn.cif -o 1crn.obj
```

**XYZ input (keep hydrogens):**

```console
$ molses molecule.xyz -o molecule.obj --keep-hydrogens
```

**Multi-model NMR ensemble (compute each model):**

```console
$ for i in $(seq 1 20); do
    molses ensemble.pdb -o model_$i.obj --model $i
  done
```

**Custom probe radius (e.g. larger solvent):**

```console
$ molses 1crn.pdb -o out.obj --probe-radius 1.8
```

**Custom radii from a CSV file:**

```console
$ cat my_radii.csv
# element,radius
C,1.75
N,1.60
$ molses 1crn.pdb -o out.obj --radii custom:my_radii.csv
```

**Smoothed surface for visualization:**

```console
$ molses 1crn.pdb -o pretty.obj --resolution fine --smooth 3
```

**Keep internal cavities (e.g. for pocket analysis):**

```console
$ molses 1crn.pdb -o out.obj --all-components
```

**MSMS-format output** (for VMD/Chimera pipelines):

```console
$ molses 1crn.pdb -o 1crn_msms --format msms
```

**Limit memory usage:**

```console
$ molses 4hhb.pdb -o out.obj --max-memory 512 --resolution coarse
```

### Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Invalid input or arguments |
| 2 | I/O error |
| 3 | Geometry failure |
| 4 | Resource limit exceeded |

---

## Library usage

### Basic example

```rust
use molses::{parse, radii::RadiiSet, SesOptions};

fn main() -> molses::Result<()> {
    // 1. Parse the structure file (format is auto-detected).
    let path = std::path::Path::new("1crn.pdb");
    let (mut structure, mut report) = parse::parse_file(path)?;

    // 2. Assign van der Waals radii (Bondi 1964 by default).
    parse::assign_radii(&mut structure, &RadiiSet::Bondi, false, &mut report)?;

    // 3. Merge duplicate centers and drop nested spheres.
    parse::deduplicate(&mut structure, &mut report);
    eprintln!("atoms: {}", structure.len());

    // 4. Configure and compute the SES mesh.
    let opts = SesOptions {
        probe_radius: 1.4,     // water
        grid_spacing: 0.35,    // balanced resolution
        ..SesOptions::default()
    };
    let mesh = molses::ses(&structure, &opts)?;

    // 5. Inspect the mesh.
    println!(
        "vertices: {}, triangles: {}, area: {:.1} Å², volume: {:.1} Å³",
        mesh.vertex_count(),
        mesh.triangle_count(),
        mesh.area(),
        mesh.volume(),
    );

    // 6. Write an OBJ file.
    let meta = molses::io::WriteMeta {
        tool: format!("molses {}", molses::version()),
        source: "1crn.pdb".into(),
        probe_radius: 1.4,
        grid_spacing: 0.35,
        atoms: structure.len(),
        deterministic: true,
    };
    let mut f = std::io::BufWriter::new(std::fs::File::create("1crn.obj")?);
    molses::io::obj::write_obj(&mut f, &mesh, &meta)?;

    Ok(())
}
```

### API overview

The library has a small surface area. The types you need:

| Type | Purpose |
|---|---|
| [`Structure`](#structure) | Parsed molecule: list of [`Atom`]s |
| [`Atom`] | Position, radius, element, residue |
| [`SesOptions`](#sesoptions) | Probe radius, grid spacing, policies |
| [`Mesh`](#mesh) | Positions, normals, triangle indices |
| [`SesError`](#seserror) | Typed errors (no panics on bad input) |

### Parsing structures

```rust
use molses::parse::{self, ParseOptions, AltlocPolicy, ModelSelection};

// Default parsing (first model, highest-occupancy altLoc, no hydrogens).
let (structure, report) = parse::parse_file(path)?;

// Custom selection.
let opts = ParseOptions {
    keep_hydrogens: true,
    include_waters: true,
    include_hetatm: true,
    altloc_policy: AltlocPolicy::All,
    model: ModelSelection::All,     // keep all NMR models
    chains: vec!["A".into(), "B".into()],
    strict_elements: false,         // warn on unknown elements instead of erroring
};
let (structure, report) = parse::parse_file_with(path, &opts)?;
```

After parsing, assign radii and clean up:

```rust
use molses::radii::RadiiSet;

let mut structure = structure;
let mut report = report;
parse::assign_radii(&mut structure, &RadiiSet::Bondi, false, &mut report)?;
parse::deduplicate(&mut structure, &mut report);
```

### Custom radii

```rust
use molses::radii::RadiiSet;

// Bondi (1964) compiled-in table.
let bondi = RadiiSet::Bondi;

// Bondi with per-element overrides from a CSV file (`element,radius` lines).
let custom = RadiiSet::from_csv("my_radii.csv".as_ref(), RadiiSet::Bondi)?;

// Override individual elements in code.
use std::collections::BTreeMap;
let mut ov = BTreeMap::new();
ov.insert("C".to_string(), 1.75);
let custom = RadiiSet::Custom {
    base: Box::new(RadiiSet::Bondi),
    overrides: ov,
};
```

### Computing surfaces

```rust
use molses::{SesOptions, SurfaceKind, ComponentPolicy, radii::RadiiSet};

let opts = SesOptions {
    probe_radius: 1.4,
    grid_spacing: 0.35,
    surface: SurfaceKind::Ses,      // or Sas, Vdw
    components: ComponentPolicy::OuterOnly, // or All, NoCavities { min_volume: 50.0 }
    smoothing: 3,                   // Taubin iterations (0 = none)
    radii: RadiiSet::Bondi,
    threads: None,                  // Some(n) to limit parallelism
    max_memory_bytes: Some(2 << 30), // 2 GB soft budget
    nearest_atom: true,             // per-vertex nearest-atom indices
};

let mesh = molses::ses(&structure, &opts)?;
```

### Working with the mesh

```rust
// Basic statistics.
let area = mesh.area();          // Å²
let volume = mesh.volume();      // Å³ (signed; positive = outward normals)
let n_v = mesh.vertex_count();
let n_t = mesh.triangle_count();

// Vertex data (all parallel arrays).
for i in 0..n_v {
    let p = mesh.positions[i];   // [f32; 3] in Å
    let n = mesh.normals[i];     // [f32; 3] unit, outward
    let a = mesh.nearest_atom[i]; // u32 atom index (if nearest_atom: true)
}

// Triangle indices (0-based, CCW viewed from outside).
for tri in &mesh.indices {
    let [i0, i1, i2] = *tri;     // [u32; 3]
}

// Topology audit (watertightness, Euler characteristic).
use molses::mesh::topology;
let report = topology::audit(&mesh);
assert!(report.is_watertight());
```

### Writing OBJ files

```rust
use molses::io::{self, WriteMeta};

let meta = WriteMeta {
    tool: format!("molses {}", molses::version()),
    source: "1crn.pdb".into(),
    probe_radius: 1.4,
    grid_spacing: 0.35,
    atoms: structure.len(),
    deterministic: true,   // suppress timestamps for reproducible output
};

let mut buf = std::io::BufWriter::new(std::fs::File::create("out.obj")?);
io::obj::write_obj(&mut buf, &mesh, &meta)?;
```

### Error handling

All functions return `Result<T, SesError>` (or `molses::Result<T>`). The library
never panics on malformed input.

```rust
use molses::SesError;

match molses::ses(&structure, &opts) {
    Ok(mesh) => println!("{} triangles", mesh.triangle_count()),
    Err(SesError::InvalidParameter { param, value, reason }) => {
        eprintln!("bad {param} = {value}: {reason}");
    }
    Err(SesError::ResourceLimit { needed, limit }) => {
        eprintln!("need {needed}, limit {limit}");
    }
    Err(e) => eprintln!("error: {e}"),
}
```

### Multi-model (NMR ensemble) processing

```rust
use molses::parse::{self, ParseOptions, ModelSelection};

let (structure, mut report) = parse::parse_file_with(
    "ensemble.pdb".as_ref(),
    &ParseOptions { model: ModelSelection::All, ..Default::default() },
)?;
parse::assign_radii(&mut structure, &RadiiSet::Bondi, false, &mut report)?;

// Each model is an index range into structure.atoms().
for (i, range) in structure.models().iter().enumerate() {
    let model_atoms = &structure.atoms()[range.clone()];
    // Build a single-model structure and compute its surface.
    // ...
}
```

### Performance tips

- Use the default `parallel` + `simd` features (on by default).
- Set `grid_spacing` coarser for interactive use (`0.5`–`0.75`).
- For very large structures, set `max_memory_bytes` to enable slab chunking.
- The `nearest_atom` option adds a small cost; disable it if unneeded.
- Batch processing: structures are independent — parallelize at the file level.

---

## Mathematical definition

For atoms `(c_i, r_i)` and probe radius `r_p`, let `K = ⋃ B(c_i, r_i + r_p)` be
the union of inflated balls (the SAS cavity) and `f_sas` its signed distance
function. The SES is the zero level set of

```
f_ses(p) = f_sas(p) + r_p
```

(Quan & Stamm 2016, *J. Comput. Phys.*, Eq. 3.4). The implementation samples
this field on a regular grid of spacing `h`, computes an exact Euclidean
distance transform of the voxelized solid, and extracts the isosurface with
marching tetrahedra — a robust, parallel, and SIMD-friendly route that trades
analytic singularity handling for a controllable resolution/accuracy knob.

The SES consists of three kinds of patches (Connolly 1983):

| Patch | Geometry | Probe contacts |
|---|---|---|
| Contact (convex spherical) | Spherical patch of radius `r_i` on atom `i` | One atom |
| Toric (toroidal reentrant) | Saddle-shaped torus patch (minor radius `r_p`) | Two atoms |
| Spheric reentrant (concave) | Spherical triangle of radius `r_p` | Three atoms |

---

## Accuracy and level of detail

Grid spacing `h` is the level-of-detail control. Geometric error is `O(h)`
(Hausdorff distance). Triangle count scales as `≈ Area / h²`.

Validation on analytic ground truths (single sphere, two-sphere):

| Metric | `h = 0.35` | `h = 0.15` |
|---|---|---|
| Area error | < 5% | < 2% |
| Volume error | < 5% | < 2% |
| Watertight | yes | yes |

For a feature of width `w` to be resolved, use `h ≤ w / 3`. The toroidal
reentrant patches have minor radius `r_p = 1.4 Å`, so `h ≤ 0.3 Å` resolves
their curvature well.

The mesh is always **watertight** (every edge shared by exactly two triangles,
outward-oriented) and **deterministic** (bitwise-identical output regardless of
thread count).

---

## Tools

Development helpers that are not part of the published crates.

| Path | What it is |
|---|---|
| [`tools/meshview/`](tools/meshview) | Headless numpy+Pillow rasterizer for comparing meshes side by side at matched camera settings |

`meshview` renders any OBJ with a fixed orthographic camera so several meshes
land pixel-for-pixel on the same world-space region — the tool that was used to
find the crevice-quality differences described above. It has no dependency on
molses and is meant to move into its own repository; see its
[README](tools/meshview/README.md) for the limits that matter when reading an
image (orthographic, Lambert-only, and the reason nothing is back-face culled).

```sh
python3 tools/meshview/render.py ref.obj out.obj -o cmp.png --zoom 110 \
    --focus 4.29,5.05,4.55 --label reference --label ours
python3 tools/meshview/selftest.py
```

A render is for eyeballing; the numbers that decide quality are the `roughness`,
`normal_dev` and watertightness the CLI logs for every mesh.

---

## Feature flags

| Feature | Default | Description |
|---|---|---|
| `parallel` | ✓ | Rayon data parallelism |
| `simd` | ✓ | Portable SIMD via the `wide` crate (stable Rust) |
| `serde` | | Serialization of core types |
| `mmap` | | Memory-mapped fields for very large grids |
| `gzip` | | Transparent reading of `.gz` inputs |

Disable defaults for a minimal build:

```toml
[dependencies]
molses = { version = "0.2", default-features = false }
```

---

## van der Waals radii

Default table: Bondi (1964), *J. Phys. Chem.* 68(3), 441–451,
doi:10.1021/j100785a001. Override per element with `--radii custom:file.csv`
(`element,radius` lines, Å) or programmatically via `RadiiSet::Custom`.

---

## License

MIT OR Apache-2.0.

## References

See `PRODUCT_SPEC.md` §13 for the full bibliography: Lee & Richards 1971;
Richards 1977; Greer & Bush 1978; Connolly 1983; Bondi 1964; Sanner et al. 1996;
Quan & Stamm 2016; Lorensen & Cline 1987; Xu & Zhang 2009; Decherchi & Rocchia
2013; Goddard 2013.
