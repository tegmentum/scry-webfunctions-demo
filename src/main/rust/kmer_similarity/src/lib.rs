//! k-mer Jaccard similarity as a WebAssembly component.
//!
//! Reproduces the "sequence similarity" role of BLAST in the SCRY paper
//! (Stringer et al., ceur-ws.org/Vol-1795/paper30.pdf) with a much simpler,
//! self-contained algorithm suitable for a WASM component:
//!
//! * Extract the set of k-mers (default k=3) from each of two protein
//!   sequences.
//! * Return the Jaccard coefficient |A ∩ B| / |A ∪ B| as an `xsd:decimal`
//!   literal in `[0, 1]`.
//!
//! Real BLAST does gapped local alignment with a substitution matrix —
//! considerably heavier and requires I/O to a sequence database. This
//! component demonstrates the same shape (embed a sequence-scoring
//! procedure inside a SPARQL query at query time) with an algorithm that
//! ports cleanly to a sandboxed pure-Rust component.
//!
//! Input contract: two args, each a `Term::Literal` whose `value` is the
//! protein sequence in single-letter amino-acid code. Any non-literal
//! args are rejected with a descriptive error.
//!
//! Migrated onto the substrate WIT: base tegmentum:webfunction@0.1.0
//! (extension / aggregate / property-function) + Stardog overlay
//! stardog:webfunction@0.3.0 (planner cardinality-estimate + doc). Under
//! the base surface a filter function returns a single `term`, which
//! matches how kmer_similarity was consumed pre-migration
//! (`wf:call(<url>, ?a, ?b)` in filter form binding the single output).

#[allow(warnings)]
mod bindings;

use std::collections::HashSet;

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

const XSD_DECIMAL: &str = "http://www.w3.org/2001/XMLSchema#decimal";
const XSD_STRING: &str = "http://www.w3.org/2001/XMLSchema#string";
const FUNCTION_NAME: &str = "kmer-similarity";

// Peptide k-mer size. 3 keeps the demo intuitive; larger k gives better
// discrimination on longer sequences at higher memory cost.
const K: usize = 3;

struct Component;

fn sequence_of(arg: &WitTerm) -> Result<&str, String> {
    match arg {
        WitTerm::Literal(l) => Ok(l.value.as_str()),
        WitTerm::NamedNode(_) => Err("kmer-similarity: argument must be a literal, got IRI".into()),
        WitTerm::BlankNode(_) => {
            Err("kmer-similarity: argument must be a literal, got blank node".into())
        }
        WitTerm::Triple(_) => {
            Err("kmer-similarity: argument must be a literal, got quoted triple".into())
        }
    }
}

fn kmers(sequence: &str, k: usize) -> HashSet<String> {
    if sequence.len() < k {
        return HashSet::new();
    }
    let bytes = sequence.as_bytes();
    let mut out = HashSet::with_capacity(bytes.len().saturating_sub(k - 1));
    for window in bytes.windows(k) {
        // Sequences are ASCII (single-letter AA codes), so String::from_utf8
        // never fails on a valid window.
        out.insert(String::from_utf8(window.to_vec()).unwrap());
    }
    out
}

fn jaccard(a: &HashSet<String>, b: &HashSet<String>) -> f64 {
    if a.is_empty() && b.is_empty() {
        return 1.0;
    }
    let intersection = a.intersection(b).count() as f64;
    let union = a.union(b).count() as f64;
    if union == 0.0 {
        0.0
    } else {
        intersection / union
    }
}

fn decimal_literal(v: f64) -> WitTerm {
    WitTerm::Literal(WitLiteral {
        value: format!("{:.4}", v),
        datatype: Some(XSD_DECIMAL.to_string()),
        language: None,
    })
}

fn string_literal(s: &str) -> WitTerm {
    WitTerm::Literal(WitLiteral {
        value: s.to_string(),
        datatype: Some(XSD_STRING.to_string()),
        language: None,
    })
}

fn similarity(args: &[WitTerm]) -> Result<WitTerm, String> {
    if args.len() != 2 {
        return Err(format!(
            "kmer-similarity: expected 2 args (query_seq, subject_seq), got {}",
            args.len()
        ));
    }
    let query = sequence_of(&args[0])?;
    let subject = sequence_of(&args[1])?;
    Ok(decimal_literal(jaccard(&kmers(query, K), &kmers(subject, K))))
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
            FUNCTION_NAME => similarity(&args),
            other => Err(format!("kmer-similarity: unknown filter function '{other}'")),
        }
    }
}

/// Aggregate stub — this component has no aggregates. Required by the
/// shared `sparql-extension` world.
impl AggregateGuest for Component {
    type AggregateState = UnreachableState;

    fn register_aggregates() -> Vec<AggregateDescriptor> {
        Vec::new()
    }

    fn new_aggregate(name: String) -> Result<AggregateState, String> {
        Err(format!(
            "kmer-similarity: unknown aggregate '{name}' (this component provides none)"
        ))
    }
}

pub struct UnreachableState;

impl GuestAggregateState for UnreachableState {
    fn step(&self, _args: Vec<WitTerm>) -> Result<(), String> {
        Err("kmer-similarity: aggregate state was never constructed".into())
    }

    fn finish(&self) -> Result<WitTerm, String> {
        Err("kmer-similarity: aggregate state was never constructed".into())
    }
}

/// Property-function stub — this component has no property functions.
/// The single-column output shape is already covered by `extension.call`.
impl PropertyFunctionGuest for Component {
    fn register_property_functions() -> Vec<PropertyDescriptor> {
        Vec::new()
    }

    fn evaluate(
        name: String,
        _subjects: Vec<WitTerm>,
        _objects: Vec<WitTerm>,
    ) -> Result<Vec<BindingRow>, String> {
        Err(format!(
            "kmer-similarity: unknown property function '{name}' (this component provides none)"
        ))
    }
}

/// Stardog planner cardinality hint: kmer-similarity is a scalar
/// per-row function, so its result cardinality equals the input.
impl PlannerGuest for Component {
    fn cardinality_estimate(input: Cardinality, _args: Vec<WitTerm>) -> Result<Cardinality, String> {
        Ok(Cardinality {
            value: input.value,
            accuracy: Accuracy::Accurate,
        })
    }
}

/// Stardog `doc` self-description. Bindings flat (one row of one
/// variable) per the overlay's convention.
impl DocGuest for Component {
    fn doc() -> Vec<WitBinding> {
        vec![WitBinding {
            variable: "doc".to_string(),
            value: string_literal(
                "kmer-similarity(query_seq, subject_seq) -> similarity: \
                 Jaccard coefficient over k=3 amino-acid k-mers. \
                 A simple, WASM-friendly stand-in for the BLAST procedure \
                 in the SCRY paper (Stringer et al., CEUR Vol-1795, paper 30).",
            ),
        }]
    }
}

bindings::export!(Component with_types_in bindings);
