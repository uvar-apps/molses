//! Error types for the `molses` library.
//!
//! The library never panics on malformed input; all failure modes are reported
//! through [`SesError`].

use std::fmt;
use std::path::PathBuf;

/// Result alias used throughout the crate.
pub type Result<T> = std::result::Result<T, SesError>;

/// All errors produced by `molses`.
#[derive(Debug, thiserror::Error)]
pub enum SesError {
    /// Underlying I/O failure.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// The structure file could not be parsed.
    #[error("parse error in {path}{}: {msg}", line.map(|l| format!(" line {l}")).unwrap_or_default())]
    Parse {
        /// Path (or description) of the offending input.
        path: PathBuf,
        /// 1-based line number, when known.
        line: Option<usize>,
        /// Human-readable explanation.
        msg: String,
    },

    /// Selection filters removed every atom.
    #[error("atom selection is empty (all atoms filtered out)")]
    EmptySelection,

    /// An element symbol has no radius in the active table and no override.
    #[error("unknown element `{symbol}` (no van der Waals radius)")]
    UnknownElement {
        /// The offending element symbol as written in the file.
        symbol: String,
    },

    /// A user-supplied parameter is out of range or not finite.
    #[error("invalid parameter `{param}`: {value} ({reason})")]
    InvalidParameter {
        /// Parameter name as exposed on the API/CLI.
        param: &'static str,
        /// The rejected value.
        value: String,
        /// Why it was rejected.
        reason: String,
    },

    /// A resource budget (memory, triangle count) could not be met.
    #[error("resource limit exceeded: need {needed}, limit {limit}")]
    ResourceLimit {
        /// What was needed (bytes, triangles, …).
        needed: String,
        /// What was allowed.
        limit: String,
    },

    /// A geometry stage failed in a way that is not attributable to input data.
    #[error("geometry failure in stage `{stage}`: {msg}")]
    Geometry {
        /// Pipeline stage name (e.g. `"edt"`, `"polygonise"`).
        stage: &'static str,
        /// Explanation.
        msg: String,
    },

    /// Unsupported file format.
    #[error("unsupported input format: {0}")]
    UnsupportedFormat(String),
}

impl SesError {
    /// Convenience constructor for parse errors.
    pub fn parse(path: impl Into<PathBuf>, line: Option<usize>, msg: impl Into<String>) -> Self {
        Self::Parse {
            path: path.into(),
            line,
            msg: msg.into(),
        }
    }

    /// Convenience constructor for invalid parameters.
    pub fn invalid_param(
        param: &'static str,
        value: impl fmt::Display,
        reason: impl Into<String>,
    ) -> Self {
        Self::InvalidParameter {
            param,
            value: value.to_string(),
            reason: reason.into(),
        }
    }
}
