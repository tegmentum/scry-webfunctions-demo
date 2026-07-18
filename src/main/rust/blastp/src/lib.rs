//! blastp — protein local alignment as a WebAssembly component.
//!
//! Wraps rust-bio's Smith-Waterman aligner (`bio::alignment::pairwise::Aligner`)
//! with BLOSUM62 scoring and NCBI's BLASTP-default affine gap penalties
//! (open=-11, extend=-1) to compute the same local-alignment scores NCBI BLAST
//! reports for `blastp`. Adds Karlin-Altschul statistics on top so we report
//! bit scores and E-values in the standard convention.
//!
//! This is the exact algorithm underlying BLAST — BLAST is a *heuristic* for
//! this local-alignment problem, useful for scanning gigabyte databases. For
//! small candidate sets like the demo query (a few dozen proteins), running
//! the exact algorithm directly is both faster and more accurate.
//!
//! Reusable outside this demo: any consumer that has two amino-acid sequences
//! and wants BLASTP-shape output can invoke this component unchanged.
//!
//! ## WIT surface (post-migration)
//!
//! Under the substrate WIT (`tegmentum:webfunction@0.1.0` base +
//! `stardog:webfunction@0.3.0` overlay) `evaluate` no longer returns a
//! `binding-sets`. Filter functions return a single `term`; multi-column
//! output moves to the `property-function` interface.
//!
//! blastp exposes both:
//!
//! * `extension.call("blastp", [query, subject])` -> bit_score as an
//!   `xsd:decimal` term. Matches the filter-form usage in
//!   `coexpression.rq` which binds the single returned value into
//!   `?bitScore`.
//! * `property-function.evaluate("blastp", subjects=[query, subject],
//!   objects=[])` -> one row of five terms:
//!   `[bit_score, raw_score, e_value, identity, align_length]`.
//!   Matches a hypothetical tuple-form consumer that wants all five
//!   outputs in one call.

#[allow(warnings)]
mod bindings;

use bio::alignment::pairwise::Aligner;
use bio::scores::blosum62;

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

// The single registered name — same string for both the filter form
// (`extension.call`) and the tuple/property-function form
// (`property-function.evaluate`). Downstream plugins dispatch to whichever
// entry point the SPARQL syntax indicates.
const FUNCTION_NAME: &str = "blastp";

// NCBI BLASTP defaults: open=-11, extend=-1, BLOSUM62.
const GAP_OPEN: i32 = -11;
const GAP_EXTEND: i32 = -1;

// Karlin-Altschul parameters for BLOSUM62 with 11/1 affine gap penalties.
// Source: NCBI BLAST parameter tables (published statistics for this matrix +
// gap combination). Values are canonical constants, not fitted here.
const K: f64 = 0.041;
const LAMBDA: f64 = 0.267;

struct Component;

fn string_literal(label: &str) -> WitTerm {
    WitTerm::Literal(WitLiteral {
        value: label.into(),
        datatype: Some(XSD_STRING.into()),
        language: None,
    })
}

fn decimal_literal(v: f64) -> WitTerm {
    WitTerm::Literal(WitLiteral {
        value: format!("{:.4}", v),
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

fn scientific_literal(v: f64) -> WitTerm {
    // E-values span many orders of magnitude — scientific notation is what
    // NCBI blastp reports and what humans expect.
    WitTerm::Literal(WitLiteral {
        value: format!("{:.2e}", v),
        datatype: Some(XSD_DECIMAL.into()),
        language: None,
    })
}

fn sequence_of(arg: &WitTerm) -> Result<&str, String> {
    match arg {
        WitTerm::Literal(l) => Ok(l.value.as_str()),
        WitTerm::NamedNode(_) => Err("blastp: argument must be a literal, got IRI".into()),
        WitTerm::BlankNode(_) => {
            Err("blastp: argument must be a literal, got blank node".into())
        }
        WitTerm::Triple(_) => {
            Err("blastp: argument must be a literal, got quoted triple".into())
        }
    }
}

/// Percent identity across the aligned region, from a rust-bio alignment path.
fn percent_identity(alignment: &bio::alignment::Alignment) -> f64 {
    let mut matches = 0usize;
    let mut length = 0usize;
    for op in &alignment.operations {
        match op {
            bio::alignment::AlignmentOperation::Match => {
                matches += 1;
                length += 1;
            }
            bio::alignment::AlignmentOperation::Subst
            | bio::alignment::AlignmentOperation::Del
            | bio::alignment::AlignmentOperation::Ins => {
                length += 1;
            }
            _ => {}
        }
    }
    if length == 0 {
        0.0
    } else {
        (matches as f64 / length as f64) * 100.0
    }
}

/// Core computation: yield the five raw metrics from a pair of sequences.
/// Wrapped by both entry points below so bit-score-only filter callers
/// pay the same cost as full-tuple callers (the alignment is the work).
fn compute(query: &[u8], subject: &[u8]) -> (f64, f64, f64, f64, i64) {
    let mut aligner = Aligner::with_capacity(
        query.len(),
        subject.len(),
        GAP_OPEN,
        GAP_EXTEND,
        blosum62,
    );
    let alignment = aligner.local(query, subject);
    let raw_score = alignment.score as f64;
    // Karlin-Altschul: bit score S' = (λS - lnK) / ln2.
    let bit_score = (LAMBDA * raw_score - K.ln()) / std::f64::consts::LN_2;
    // Expected number of alignments with score ≥ S under the null model:
    // E = K * m * n * exp(-λS). m,n are effective sequence lengths — for
    // pairwise (not database search) we approximate with the raw lengths.
    let m = query.len() as f64;
    let n = subject.len() as f64;
    let e_value = K * m * n * (-LAMBDA * raw_score).exp();
    let identity_pct = percent_identity(&alignment);
    let align_len = alignment.operations.len() as i64;
    (bit_score, raw_score, e_value, identity_pct, align_len)
}

/// Filter form: return the bit score as a single term. Matches how
/// `coexpression.rq` consumes `wf:call(<url>, ?a, ?b)` — the plugin
/// binds the single returned value to `?bitScore`.
fn call_filter(args: &[WitTerm]) -> Result<WitTerm, String> {
    if args.len() != 2 {
        return Err(format!(
            "blastp: expected 2 args (query_seq, subject_seq), got {}",
            args.len()
        ));
    }
    let query = sequence_of(&args[0])?.as_bytes();
    let subject = sequence_of(&args[1])?.as_bytes();
    let (bit_score, _, _, _, _) = compute(query, subject);
    Ok(decimal_literal(bit_score))
}

impl ExtensionGuest for Component {
    fn register() -> Vec<FunctionDescriptor> {
        vec![FunctionDescriptor {
            name: FUNCTION_NAME.to_string(),
            min_arity: 2,
            max_arity: Some(2),
        }]
    }

    fn call(name: String, args: Vec<WitTerm>) -> Result<WitTerm, String> {
        match name.as_str() {
            FUNCTION_NAME => call_filter(&args),
            other => Err(format!("blastp: unknown filter function '{other}'")),
        }
    }
}

/// Property-function form: return all five outputs in one call. The
/// caller passes the two input sequences via `subjects`; `objects` is
/// unused (the SPARQL tuple form takes the output vars from the
/// object-side collection, which the host binds from the returned row).
impl PropertyFunctionGuest for Component {
    fn register_property_functions() -> Vec<PropertyDescriptor> {
        vec![PropertyDescriptor {
            name: FUNCTION_NAME.to_string(),
            // Two inputs on the subject side: query_seq, subject_seq.
            subject_arity: 2,
            // Five outputs: bit_score, raw_score, e_value, identity, align_length.
            object_arity: 5,
        }]
    }

    fn evaluate(
        name: String,
        subjects: Vec<WitTerm>,
        _objects: Vec<WitTerm>,
    ) -> Result<Vec<BindingRow>, String> {
        if name != FUNCTION_NAME {
            return Err(format!(
                "blastp: unknown property function '{name}' (only '{FUNCTION_NAME}')"
            ));
        }
        if subjects.len() != 2 {
            return Err(format!(
                "blastp: expected 2 subject args (query_seq, subject_seq), got {}",
                subjects.len()
            ));
        }
        let query = sequence_of(&subjects[0])?.as_bytes();
        let subject = sequence_of(&subjects[1])?.as_bytes();
        let (bit_score, raw_score, e_value, identity_pct, align_len) = compute(query, subject);
        Ok(vec![BindingRow {
            values: vec![
                decimal_literal(bit_score),
                decimal_literal(raw_score),
                scientific_literal(e_value),
                decimal_literal(identity_pct),
                integer_literal(align_len),
            ],
        }])
    }
}

/// Aggregate stub — blastp has no aggregates. Required by the shared
/// `sparql-extension` world.
impl AggregateGuest for Component {
    type AggregateState = UnreachableState;

    fn register_aggregates() -> Vec<AggregateDescriptor> {
        Vec::new()
    }

    fn new_aggregate(name: String) -> Result<AggregateState, String> {
        Err(format!(
            "blastp: unknown aggregate '{name}' (this component provides none)"
        ))
    }
}

pub struct UnreachableState;

impl GuestAggregateState for UnreachableState {
    fn step(&self, _args: Vec<WitTerm>) -> Result<(), String> {
        Err("blastp: aggregate state was never constructed".into())
    }

    fn finish(&self) -> Result<WitTerm, String> {
        Err("blastp: aggregate state was never constructed".into())
    }
}

/// Stardog planner cardinality: blastp is a per-row scalar (filter
/// form) or a per-row one-row property function (tuple form). Either
/// way the output cardinality matches the input.
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
                "blastp(query_seq, subject_seq): protein local alignment via \
                 rust-bio's Smith-Waterman with BLOSUM62 + NCBI blastp default \
                 gap penalties (open=-11, extend=-1). Filter form returns \
                 bit_score (xsd:decimal). Tuple/property-function form returns \
                 (bit_score, raw_score, e_value, identity, align_length). \
                 Karlin-Altschul statistics use canonical BLOSUM62/11/1 \
                 constants (K=0.041, λ=0.267).",
            ),
        }]
    }
}

bindings::export!(Component with_types_in bindings);
