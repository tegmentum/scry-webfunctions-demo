//! protparam — protein biochemical properties as a WebAssembly component.
//!
//! Analog of the Expasy ProtParam tool (web.expasy.org/protparam). Given an
//! amino-acid sequence, computes:
//!
//! * `length`             — residue count
//! * `molecular_weight`   — average MW in Daltons (residue masses + one H2O)
//! * `theoretical_pi`     — isoelectric point via Bjellqvist pKa + bisection
//! * `gravy`              — grand average of hydropathy (Kyte-Doolittle)
//!
//! This is a self-contained numeric procedure (no external data, no I/O) — the
//! kind of thing SPARQL can't do natively but that is trivial to embed as a
//! sandboxed Wasm component. Complements `blastp` in the demo family: BLAST
//! finds *which* proteins are similar, protparam characterises *what* they
//! are like biochemically.
//!
//! ## WIT surface (post-migration)
//!
//! Under the substrate WIT (`tegmentum:webfunction@0.1.0` base +
//! `stardog:webfunction@0.3.0` overlay), a filter function returns a single
//! `term`. protparam has four independent outputs, so it is a
//! property-function rather than a filter — the base
//! `property-function.evaluate` interface returns
//! `list<binding-row>` where each row carries a `list<term>`.
//!
//! Consumers reach protparam via the tuple form:
//!   `(<url> ?seq) wf:call (?length ?molecularWeight ?theoreticalPi ?gravy)`
//! The host builds `subjects=[?seq]` and expects a one-row result whose
//! `values` positional order matches the object-side variables.

#[allow(warnings)]
mod bindings;

use bindings::exports::stardog::webfunction::doc::Guest as DocGuest;
use bindings::exports::stardog::webfunction::planner::{
    Accuracy, Cardinality, Guest as PlannerGuest,
};
use bindings::exports::tegmentum::webfunction::aggregate::{
    AggregateDescriptor, AggregateState, Guest as AggregateGuest, GuestAggregateState,
};
use bindings::exports::tegmentum::webfunction::extension::{
    FunctionDescriptor, Guest as ExtensionGuest,
};
use bindings::exports::tegmentum::webfunction::property_function::{
    BindingRow, Guest as PropertyFunctionGuest, PropertyDescriptor,
};
use bindings::tegmentum::webfunction::types::{
    Binding as WitBinding, Literal as WitLiteral, Term as WitTerm,
};

const XSD_STRING: &str = "http://www.w3.org/2001/XMLSchema#string";
const XSD_DECIMAL: &str = "http://www.w3.org/2001/XMLSchema#decimal";
const XSD_INTEGER: &str = "http://www.w3.org/2001/XMLSchema#integer";

const FUNCTION_NAME: &str = "protparam";

// Average residue masses (Da) — the mass of the amino acid MINUS one water,
// suitable for polypeptide sum + one water added at the end. Standard values
// used by Expasy ProtParam.
fn residue_mass(aa: u8) -> Option<f64> {
    Some(match aa {
        b'A' => 71.0788,  b'R' => 156.1875, b'N' => 114.1038, b'D' => 115.0886,
        b'C' => 103.1388, b'E' => 129.1155, b'Q' => 128.1307, b'G' => 57.0519,
        b'H' => 137.1411, b'I' => 113.1594, b'L' => 113.1594, b'K' => 128.1741,
        b'M' => 131.1926, b'F' => 147.1766, b'P' => 97.1167,  b'S' => 87.0782,
        b'T' => 101.1051, b'W' => 186.2132, b'Y' => 163.1760, b'V' => 99.1326,
        _ => return None, // Ignore ambiguity codes (B, Z, X, *, ...).
    })
}

// Kyte-Doolittle hydropathy index. GRAVY = mean(hydropathy over residues).
fn hydropathy(aa: u8) -> Option<f64> {
    Some(match aa {
        b'A' => 1.8,  b'R' => -4.5, b'N' => -3.5, b'D' => -3.5,
        b'C' => 2.5,  b'Q' => -3.5, b'E' => -3.5, b'G' => -0.4,
        b'H' => -3.2, b'I' => 4.5,  b'L' => 3.8,  b'K' => -3.9,
        b'M' => 1.9,  b'F' => 2.8,  b'P' => -1.6, b'S' => -0.8,
        b'T' => -0.7, b'W' => -0.9, b'Y' => -1.3, b'V' => 4.2,
        _ => return None,
    })
}

// Bjellqvist pKa values (as used by ProtParam) for titratable groups.
const PKA_NTERM: f64 = 7.5;
const PKA_CTERM: f64 = 3.55;
const PKA_ARG:   f64 = 12.0;
const PKA_HIS:   f64 = 5.98;
const PKA_LYS:   f64 = 10.0;
const PKA_ASP:   f64 = 4.05;
const PKA_GLU:   f64 = 4.45;
const PKA_CYS:   f64 = 9.0;
const PKA_TYR:   f64 = 10.0;

fn count(seq: &[u8], target: u8) -> usize {
    seq.iter().filter(|&&b| b == target).count()
}

/// Net charge of a peptide at the given pH. Positive contributions from
/// N-terminus, R, H, K; negative from C-terminus, D, E, C, Y.
fn net_charge(seq: &[u8], ph: f64) -> f64 {
    let pos = |pka: f64, n: usize| n as f64 / (1.0 + 10f64.powf(ph - pka));
    let neg = |pka: f64, n: usize| -(n as f64) / (1.0 + 10f64.powf(pka - ph));

    pos(PKA_NTERM, 1)
        + pos(PKA_ARG, count(seq, b'R'))
        + pos(PKA_HIS, count(seq, b'H'))
        + pos(PKA_LYS, count(seq, b'K'))
        + neg(PKA_CTERM, 1)
        + neg(PKA_ASP, count(seq, b'D'))
        + neg(PKA_GLU, count(seq, b'E'))
        + neg(PKA_CYS, count(seq, b'C'))
        + neg(PKA_TYR, count(seq, b'Y'))
}

/// Isoelectric point via bisection on pH ∈ [0, 14]. Terminates once the
/// interval is narrower than 0.01 pH units — good to two decimal places,
/// which is standard for ProtParam-style reports.
fn isoelectric_point(seq: &[u8]) -> f64 {
    let (mut low, mut high) = (0.0_f64, 14.0_f64);
    while high - low > 0.01 {
        let mid = (low + high) / 2.0;
        if net_charge(seq, mid) > 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    (low + high) / 2.0
}

fn string_literal(s: &str) -> WitTerm {
    WitTerm::Literal(WitLiteral {
        value: s.into(),
        datatype: Some(XSD_STRING.into()),
        language: None,
    })
}

fn decimal_literal(v: f64) -> WitTerm {
    WitTerm::Literal(WitLiteral {
        value: format!("{:.2}", v),
        datatype: Some(XSD_DECIMAL.into()),
        language: None,
    })
}

fn integer_literal(v: i64) -> WitTerm {
    WitTerm::Literal(WitLiteral {
        value: v.to_string(),
        datatype: Some(XSD_INTEGER.into()),
        language: None,
    })
}

fn sequence_of(arg: &WitTerm) -> Result<&str, String> {
    match arg {
        WitTerm::Literal(l) => Ok(l.value.as_str()),
        WitTerm::NamedNode(_) => Err("protparam: argument must be a literal, got IRI".into()),
        WitTerm::BlankNode(_) => {
            Err("protparam: argument must be a literal, got blank node".into())
        }
        WitTerm::Triple(_) => {
            Err("protparam: argument must be a literal, got quoted triple".into())
        }
    }
}

/// Core computation returning the four properties per sequence.
fn compute(seq: &[u8]) -> (i64, f64, f64, f64) {
    // Only count residues we recognise so ambiguity codes don't contaminate
    // MW / GRAVY; sequence length still reports the raw input length.
    let mut mw = 18.01528_f64; // one water for the terminal H + OH
    let mut hydropathy_sum = 0.0_f64;
    let mut recognized = 0_i64;
    for &b in seq {
        if let Some(m) = residue_mass(b) {
            mw += m;
            recognized += 1;
        }
        if let Some(h) = hydropathy(b) {
            hydropathy_sum += h;
        }
    }
    let gravy = if recognized == 0 {
        0.0
    } else {
        hydropathy_sum / recognized as f64
    };
    let pi = isoelectric_point(seq);
    (seq.len() as i64, mw, pi, gravy)
}

struct Component;

/// Property-function form: return length, molecular_weight, theoretical_pi,
/// gravy — one row of four terms per input sequence.
impl PropertyFunctionGuest for Component {
    fn register_property_functions() -> Vec<PropertyDescriptor> {
        vec![PropertyDescriptor {
            name: FUNCTION_NAME.to_string(),
            // One input on the subject side: the sequence literal.
            subject_arity: 1,
            // Four outputs: length, molecular_weight, theoretical_pi, gravy.
            object_arity: 4,
        }]
    }

    fn evaluate(
        name: String,
        subjects: Vec<WitTerm>,
        _objects: Vec<WitTerm>,
    ) -> Result<Vec<BindingRow>, String> {
        if name != FUNCTION_NAME {
            return Err(format!(
                "protparam: unknown property function '{name}' (only '{FUNCTION_NAME}')"
            ));
        }
        if subjects.len() != 1 {
            return Err(format!(
                "protparam: expected 1 subject arg (sequence), got {}",
                subjects.len()
            ));
        }
        let seq_str = sequence_of(&subjects[0])?;
        let (length, mw, pi, gravy) = compute(seq_str.as_bytes());
        Ok(vec![BindingRow {
            values: vec![
                integer_literal(length),
                decimal_literal(mw),
                decimal_literal(pi),
                decimal_literal(gravy),
            ],
        }])
    }
}

/// Filter stub — protparam's 4-property output does not collapse to a
/// meaningful single scalar. Consumers reach it via the property-function
/// tuple form.
impl ExtensionGuest for Component {
    fn register() -> Vec<FunctionDescriptor> {
        Vec::new()
    }

    fn call(name: String, _args: Vec<WitTerm>) -> Result<WitTerm, String> {
        Err(format!(
            "protparam: no filter function '{name}' — use the property-function form"
        ))
    }
}

/// Aggregate stub.
impl AggregateGuest for Component {
    type AggregateState = UnreachableState;

    fn register_aggregates() -> Vec<AggregateDescriptor> {
        Vec::new()
    }

    fn new_aggregate(name: String) -> Result<AggregateState, String> {
        Err(format!(
            "protparam: unknown aggregate '{name}' (this component provides none)"
        ))
    }
}

pub struct UnreachableState;

impl GuestAggregateState for UnreachableState {
    fn step(&self, _args: Vec<WitTerm>) -> Result<(), String> {
        Err("protparam: aggregate state was never constructed".into())
    }

    fn finish(&self) -> Result<WitTerm, String> {
        Err("protparam: aggregate state was never constructed".into())
    }
}

/// Stardog planner cardinality: protparam is a per-row one-row property
/// function, so output cardinality matches input.
impl PlannerGuest for Component {
    fn cardinality_estimate(input: Cardinality, _args: Vec<WitTerm>) -> Result<Cardinality, String> {
        Ok(Cardinality {
            value: input.value,
            accuracy: Accuracy::Accurate,
        })
    }
}

/// Stardog `doc` self-description.
impl DocGuest for Component {
    fn doc() -> Vec<WitBinding> {
        vec![WitBinding {
            variable: "doc".to_string(),
            value: string_literal(
                "protparam(sequence) -> (length, molecular_weight, theoretical_pi, gravy). \
                 Mirrors Expasy ProtParam: average residue masses for MW (+ one H2O for \
                 termini), Bjellqvist pKa via bisection for the isoelectric point, and \
                 Kyte-Doolittle hydropathy averaged over the sequence for GRAVY. \
                 Ambiguity codes (B, Z, X, *) are ignored in MW/GRAVY but still count \
                 toward the reported sequence length.",
            ),
        }]
    }
}

bindings::export!(Component with_types_in bindings);
