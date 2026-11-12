//! van der Waals radii assignment.
//!
//! Default table: Bondi (1964), *J. Phys. Chem.* 68(3), 441–451,
//! doi:10.1021/j100785a001 — the standard radii set in structural biology.
//! Values are in Ångström.

use std::collections::BTreeMap;
use std::path::Path;

use crate::error::{Result, SesError};
use crate::structure::Element;

/// Default radius used for elements absent from the table, when not strict.
pub const FALLBACK_RADIUS: f64 = 1.7;

/// Bondi (1964) van der Waals radii in Ångström.
///
/// Covers the elements commonly present in macromolecular structures.
/// For elements without a Bondi value the table falls back to [`FALLBACK_RADIUS`].
pub fn bondi_radius(element: Element) -> Option<f64> {
    let r = match element {
        Element::H => 1.20,
        Element::C => 1.70,
        Element::N => 1.55,
        Element::O => 1.52,
        Element::P => 1.80,
        Element::S => 1.80,
        Element::Se => 1.90,
        Element::F => 1.47,
        Element::Cl => 1.75,
        Element::Br => 1.85,
        Element::I => 1.98,
        Element::Na => 2.27,
        Element::Mg => 1.73,
        Element::K => 2.75,
        Element::Ca => 2.31,
        Element::Fe => 1.66,
        Element::Zn => 1.39,
        Element::Other(_) => return None,
    };
    Some(r)
}

/// Named radii table selection.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[derive(Default)]
pub enum RadiiSet {
    /// Bondi (1964) compiled-in table.
    #[default]
    Bondi,
    /// Bondi with per-element overrides.
    Custom {
        /// Base table to fall back to.
        base: Box<RadiiSet>,
        /// Overrides keyed by element symbol (Title case).
        overrides: BTreeMap<String, f64>,
    },
}

impl RadiiSet {
    /// Radius of an element under this table, if known.
    pub fn radius(&self, element: Element) -> Option<f64> {
        match self {
            Self::Bondi => bondi_radius(element),
            Self::Custom { base, overrides } => {
                let sym = element.symbol();
                if let Some(r) = overrides.get(&sym) {
                    return Some(*r);
                }
                base.radius(element)
            }
        }
    }

    /// Resolve an atom radius: table lookup, else fallback or error.
    pub fn resolve(&self, element: Option<Element>, strict: bool) -> Result<f64> {
        match element {
            Some(e) => {
                if let Some(r) = self.radius(e) {
                    Ok(r)
                } else if strict {
                    Err(SesError::UnknownElement { symbol: e.symbol() })
                } else {
                    Ok(FALLBACK_RADIUS)
                }
            }
            None if strict => Err(SesError::UnknownElement {
                symbol: "<none>".into(),
            }),
            None => Ok(FALLBACK_RADIUS),
        }
    }

    /// Load a CSV file of `element,radius` lines (`#` comments allowed).
    pub fn from_csv(path: &Path, base: RadiiSet) -> Result<Self> {
        let text = std::fs::read_to_string(path)?;
        let mut overrides = BTreeMap::new();
        for (i, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let (elem, rest) = line
                .split_once(',')
                .ok_or_else(|| SesError::parse(path, Some(i + 1), "expected `element,radius`"))?;
            let radius: f64 = rest.trim().parse().map_err(|_| {
                SesError::parse(path, Some(i + 1), format!("bad radius `{}`", rest.trim()))
            })?;
            if !radius.is_finite() || radius <= 0.0 {
                return Err(SesError::parse(
                    path,
                    Some(i + 1),
                    "radius must be finite and positive",
                ));
            }
            let elem = elem.trim();
            // Normalize symbol capitalization.
            let mut cs = elem.chars();
            let normalized = match cs.next() {
                Some(c) => {
                    format!(
                        "{}{}",
                        c.to_ascii_uppercase(),
                        cs.as_str().to_ascii_lowercase()
                    )
                }
                None => {
                    return Err(SesError::parse(path, Some(i + 1), "empty element symbol"));
                }
            };
            overrides.insert(normalized, radius);
        }
        Ok(Self::Custom {
            base: Box::new(base),
            overrides,
        })
    }

    /// Load a JSON object `{"C": 1.7, ...}` of element → radius.
    #[cfg(feature = "serde")]
    pub fn from_json(path: &Path, base: RadiiSet) -> Result<Self> {
        let text = std::fs::read_to_string(path)?;
        let raw: BTreeMap<String, f64> = serde_json::from_str(&text)
            .map_err(|e| SesError::parse(path, None, format!("bad JSON: {e}")))?;
        let mut overrides = BTreeMap::new();
        for (k, v) in raw {
            if !(v > 0.0) || !v.is_finite() {
                return Err(SesError::parse(
                    path,
                    None,
                    format!("radius for {k} must be finite and positive"),
                ));
            }
            let mut cs = k.chars();
            let normalized = match cs.next() {
                Some(c) => format!(
                    "{}{}",
                    c.to_ascii_uppercase(),
                    cs.as_str().to_ascii_lowercase()
                ),
                None => continue,
            };
            overrides.insert(normalized, v);
        }
        Ok(Self::Custom {
            base: Box::new(base),
            overrides,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bondi_common_elements() {
        // Spot-check against Bondi (1964).
        assert_eq!(bondi_radius(Element::H), Some(1.20));
        assert_eq!(bondi_radius(Element::C), Some(1.70));
        assert_eq!(bondi_radius(Element::N), Some(1.55));
        assert_eq!(bondi_radius(Element::O), Some(1.52));
        assert_eq!(bondi_radius(Element::S), Some(1.80));
    }

    #[test]
    fn resolve_fallback_and_strict() {
        let rs = RadiiSet::Bondi;
        assert_eq!(rs.resolve(Some(Element::C), false).unwrap(), 1.70);
        assert_eq!(rs.resolve(None, false).unwrap(), FALLBACK_RADIUS);
        assert!(rs.resolve(None, true).is_err());
        assert!(rs.resolve(Some(Element::Other(120)), true).is_err());
        assert_eq!(
            rs.resolve(Some(Element::Other(120)), false).unwrap(),
            FALLBACK_RADIUS
        );
    }

    #[test]
    fn custom_overrides() {
        let mut ov = BTreeMap::new();
        ov.insert("C".to_string(), 1.80);
        let rs = RadiiSet::Custom {
            base: Box::new(RadiiSet::Bondi),
            overrides: ov,
        };
        assert_eq!(rs.radius(Element::C), Some(1.80));
        assert_eq!(rs.radius(Element::N), Some(1.55));
    }

    #[test]
    fn csv_roundtrip() {
        let dir = std::env::temp_dir().join("molses_radii_test");
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("radii.csv");
        std::fs::write(&p, "# comment\nC,1.85\nn, 1.60\n").unwrap();
        let rs = RadiiSet::from_csv(&p, RadiiSet::Bondi).unwrap();
        assert_eq!(rs.radius(Element::C), Some(1.85));
        assert_eq!(rs.radius(Element::N), Some(1.60));
        assert_eq!(rs.radius(Element::O), Some(1.52));
    }
}
