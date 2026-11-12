//! `molses` command-line interface.
//!
//! Thin wrapper over the `molses` library: parse arguments, run the pipeline,
//! write the mesh. All geometry lives in the library crate.

use std::io::Write;
use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Parser, ValueEnum};

#[cfg(feature = "mimalloc")]
#[global_allocator]
static ALLOC: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Generate a Solvent Excluded Surface (Connolly surface) mesh from a molecular
/// structure file (PDB, mmCIF, XYZ, or XYZR) and write it as Wavefront OBJ.
#[derive(Parser, Debug)]
#[command(
    name = "molses",
    version,
    about,
    long_about = None
)]
struct Args {
    /// Structure file: .pdb / .ent / .cif / .mmcif / .xyz / .xyzr [.gz]
    input: PathBuf,

    /// Output path (defaults to stdout)
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Probe radius in Å [default: 1.4]
    #[arg(long)]
    probe_radius: Option<f64>,

    /// Grid spacing in Å, or preset name
    #[arg(long)]
    resolution: Option<ResolutionArg>,

    /// Explicit grid spacing in Å (overrides --resolution)
    #[arg(long)]
    grid_spacing: Option<f64>,

    /// van der Waals radii: bondi | custom:FILE
    #[arg(long, default_value = "bondi")]
    radii: String,

    /// Model selection: first | all | N (1-based)
    #[arg(long, default_value = "first")]
    model: String,

    /// Comma-separated chain IDs to include
    #[arg(long)]
    chains: Option<String>,

    /// Include hydrogen atoms
    #[arg(long)]
    keep_hydrogens: bool,

    /// Exclude HETATM records
    #[arg(long)]
    no_hetatm: bool,

    /// Include waters (HOH/WAT) [default: excluded]
    #[arg(long)]
    include_waters: bool,

    /// Alternate-location policy
    #[arg(long, value_enum, default_value = "highest-occupancy")]
    altloc: AltlocArg,

    /// Which surface to compute
    #[arg(long, value_enum, default_value = "ses")]
    surface: SurfaceArg,

    /// Include internal cavity surfaces
    #[arg(long)]
    all_components: bool,

    /// Laplacian smoothing iterations
    #[arg(long, default_value_t = 0)]
    smooth: u32,

    /// Worker threads [default: all]
    #[arg(long)]
    threads: Option<usize>,

    /// Soft memory budget in MB
    #[arg(long)]
    max_memory: Option<u64>,

    /// Output format
    #[arg(long, value_enum, default_value = "obj")]
    format: FormatArg,

    /// Overwrite existing output
    #[arg(long)]
    force: bool,

    /// Suppress progress output
    #[arg(long)]
    quiet: bool,

    /// Detailed progress
    #[arg(long)]
    verbose: bool,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum ResolutionArg {
    /// 0.75 Å
    Draft,
    /// 0.50 Å
    Coarse,
    /// 0.35 Å (default)
    Balanced,
    /// 0.25 Å
    Fine,
    /// 0.15 Å
    Ultra,
    /// Explicit spacing in Å (use --grid-spacing for this; kept for presets)
    #[value(skip)]
    _Number,
}

impl ResolutionArg {
    fn spacing(self) -> f64 {
        match self {
            Self::Draft => 0.75,
            Self::Coarse => 0.50,
            Self::Balanced => 0.35,
            Self::Fine => 0.25,
            Self::Ultra => 0.15,
            Self::_Number => 0.35,
        }
    }
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum AltlocArg {
    #[value(name = "highest-occupancy")]
    HighestOccupancy,
    First,
    All,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum SurfaceArg {
    Ses,
    Sas,
    Vdw,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum FormatArg {
    Obj,
    Msms,
}

fn main() -> ExitCode {
    let args = Args::parse();
    match run(args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::from(exit_code_for(&e))
        }
    }
}

fn exit_code_for(err: &anyhow::Error) -> u8 {
    if let Some(se) = err.downcast_ref::<molses::SesError>() {
        return match se {
            molses::SesError::Io(_) => 2,
            molses::SesError::Parse { .. } | molses::SesError::UnsupportedFormat(_) => 1,
            molses::SesError::EmptySelection
            | molses::SesError::UnknownElement { .. }
            | molses::SesError::InvalidParameter { .. } => 1,
            molses::SesError::ResourceLimit { .. } => 4,
            molses::SesError::Geometry { .. } => 3,
        };
    }
    1
}

fn run(args: Args) -> anyhow::Result<()> {
    if !args.quiet {
        tracing_subscriber::fmt()
            .with_env_filter(
                tracing_subscriber::EnvFilter::try_from_default_env().unwrap_or_else(|_| {
                    if args.verbose {
                        "molses=debug".into()
                    } else {
                        "molses=info".into()
                    }
                }),
            )
            .with_writer(std::io::stderr)
            .init();
    }

    let grid_spacing = match (&args.grid_spacing, &args.resolution) {
        (Some(h), _) => *h,
        (None, Some(r)) => r.spacing(),
        (None, None) => molses::SesOptions::default().grid_spacing,
    };

    let probe_radius = args.probe_radius.unwrap_or(1.4);
    let radii = parse_radii_arg(&args.radii)?;
    let model = parse_model_arg(&args.model)?;

    let parse_opts = molses::parse::ParseOptions {
        keep_hydrogens: args.keep_hydrogens,
        include_waters: args.include_waters,
        include_hetatm: !args.no_hetatm,
        altloc_policy: match args.altloc {
            AltlocArg::HighestOccupancy => molses::parse::AltlocPolicy::HighestOccupancy,
            AltlocArg::First => molses::parse::AltlocPolicy::First,
            AltlocArg::All => molses::parse::AltlocPolicy::All,
        },
        model,
        chains: args
            .chains
            .as_ref()
            .map(|c| c.split(',').map(|s| s.trim().to_string()).collect())
            .unwrap_or_default(),
        strict_elements: false,
    };

    tracing::info!(path = %args.input.display(), "parsing structure");
    let (mut structure, mut report) = molses::parse::parse_file_with(&args.input, &parse_opts)?;
    molses::parse::assign_radii(&mut structure, &radii, false, &mut report)?;
    molses::parse::deduplicate(&mut structure, &mut report);
    tracing::info!(
        read = report.atoms_read,
        kept = report.atoms_kept,
        merged = report.duplicates_merged,
        nested = report.nested_dropped,
        "structure ready"
    );

    let opts = molses::SesOptions {
        probe_radius,
        grid_spacing,
        surface: match args.surface {
            SurfaceArg::Ses => molses::SurfaceKind::Ses,
            SurfaceArg::Sas => molses::SurfaceKind::Sas,
            SurfaceArg::Vdw => molses::SurfaceKind::Vdw,
        },
        components: if args.all_components {
            molses::ComponentPolicy::All
        } else {
            molses::ComponentPolicy::OuterOnly
        },
        smoothing: args.smooth,
        radii: radii.clone(),
        threads: args.threads,
        max_memory_bytes: args.max_memory.map(|mb| mb * 1024 * 1024),
        nearest_atom: true,
    };

    tracing::info!(
        probe_radius,
        grid_spacing,
        atoms = structure.len(),
        "computing surface"
    );
    let mesh = molses::ses(&structure, &opts)?;
    let topo = molses::mesh::topology::audit(&mesh);
    tracing::info!(
        vertices = mesh.vertex_count(),
        triangles = mesh.triangle_count(),
        area = mesh.area(),
        volume = mesh.volume(),
        boundary_edges = topo.boundary_edges,
        non_manifold_edges = topo.non_manifold_edges,
        euler = topo.euler,
        watertight = topo.is_watertight(),
        "mesh ready"
    );

    let meta = molses::io::WriteMeta {
        tool: format!("molses {}", molses::version()),
        source: args.input.display().to_string(),
        probe_radius,
        grid_spacing,
        atoms: structure.len(),
        deterministic: false,
    };

    match &args.output {
        Some(path) => {
            if path.exists() && !args.force {
                anyhow::bail!("output {} exists (use --force)", path.display());
            }
            let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
            molses::io::write_mesh(&mut f, &mesh, format_name(args.format), &meta)?;
            f.flush()?;
            if !args.quiet {
                eprintln!("wrote {}", path.display());
            }
        }
        None => {
            let stdout = std::io::stdout();
            let mut f = std::io::BufWriter::new(stdout.lock());
            molses::io::write_mesh(&mut f, &mesh, format_name(args.format), &meta)?;
            f.flush()?;
        }
    }
    Ok(())
}

fn format_name(f: FormatArg) -> &'static str {
    match f {
        FormatArg::Obj => "obj",
        FormatArg::Msms => "msms",
    }
}

fn parse_radii_arg(s: &str) -> anyhow::Result<molses::radii::RadiiSet> {
    if s.eq_ignore_ascii_case("bondi") {
        return Ok(molses::radii::RadiiSet::Bondi);
    }
    if let Some(path) = s.strip_prefix("custom:") {
        let p = PathBuf::from(path);
        return Ok(molses::radii::RadiiSet::from_csv(
            &p,
            molses::radii::RadiiSet::Bondi,
        )?);
    }
    // Bare path to a CSV.
    let p = PathBuf::from(s);
    if p.exists() {
        return Ok(molses::radii::RadiiSet::from_csv(
            &p,
            molses::radii::RadiiSet::Bondi,
        )?);
    }
    anyhow::bail!("unknown radii set `{s}` (expected `bondi` or `custom:<file>`)")
}

fn parse_model_arg(s: &str) -> anyhow::Result<molses::parse::ModelSelection> {
    Ok(match s.to_ascii_lowercase().as_str() {
        "first" => molses::parse::ModelSelection::First,
        "all" => molses::parse::ModelSelection::All,
        other => {
            let n: usize = other
                .parse()
                .map_err(|_| anyhow::anyhow!("invalid --model `{s}`"))?;
            if n == 0 {
                anyhow::bail!("--model is 1-based; got 0");
            }
            molses::parse::ModelSelection::Index(n - 1)
        }
    })
}
