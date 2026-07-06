//! k-mer Jaccard similarity as a WebAssembly component.
//!
//! Reproduces the "sequence similarity" role of BLAST in the SCRY paper
//! (Stringer et al., ceur-ws.org/Vol-1795/paper30.pdf) with a much simpler,
//! self-contained algorithm suitable for a WASM component:
//!
//! * Extract the set of k-mers (default k=3) from each of two protein
//!   sequences.
//! * Return the Jaccard coefficient |A ∩ B| / |A ∪ B| as a single-row
//!   binding-set with one variable `similarity` bound to an `xsd:decimal`
//!   literal in `[0, 1]`.
//!
//! Real BLAST does gapped local alignment with a substitution matrix —
//! considerably heavier and requires I/O to a sequence database. This
//! component demonstrates the same shape (embed a sequence-scoring
//! procedure inside a SPARQL query at query time) with an algorithm that
//! ports cleanly to a sandboxed pure-Rust component.
//!
//! Input contract: two args, each a `Value::Literal` whose `label` is the
//! protein sequence in single-letter amino-acid code. Any non-literal
//! args are rejected with a descriptive error.

wit_bindgen::generate!({
    world: "webfunction",
    path: "wit",
});

use std::collections::HashSet;
use stardog::webfunction::types::{Accuracy, Binding, Literal};

struct Component;

const XSD_STRING: &str = "http://www.w3.org/2001/XMLSchema#string";
const XSD_DECIMAL: &str = "http://www.w3.org/2001/XMLSchema#decimal";

// Peptide k-mer size. 3 keeps the demo intuitive; larger k gives better
// discrimination on longer sequences at higher memory cost.
const K: usize = 3;

fn string_literal(label: &str) -> Value {
    Value::Literal(Literal {
        label: label.into(),
        datatype: XSD_STRING.into(),
        lang: None,
    })
}

fn decimal_literal(value: f64) -> Value {
    Value::Literal(Literal {
        label: format!("{:.4}", value),
        datatype: XSD_DECIMAL.into(),
        lang: None,
    })
}

fn sequence_of(arg: &Value) -> Result<&str, String> {
    match arg {
        Value::Literal(literal) => Ok(literal.label.as_str()),
        _ => Err("kmer_similarity: each argument must be a literal sequence".into()),
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
    if union == 0.0 { 0.0 } else { intersection / union }
}

impl Guest for Component {
    fn evaluate(args: Vec<Value>) -> Result<BindingSets, String> {
        if args.len() != 2 {
            return Err(format!(
                "kmer_similarity: expected 2 args (query_seq, subject_seq), got {}",
                args.len()
            ));
        }
        let query = sequence_of(&args[0])?;
        let subject = sequence_of(&args[1])?;
        let score = jaccard(&kmers(query, K), &kmers(subject, K));

        Ok(BindingSets {
            vars: vec!["similarity".into()],
            rows: vec![vec![Binding {
                name: "similarity".into(),
                value: decimal_literal(score),
            }]],
        })
    }

    fn aggregate_step(_args: Vec<Value>, _mult: u64) -> Result<(), String> {
        Err("kmer_similarity: aggregate-step not implemented".into())
    }

    fn aggregate_finish() -> Result<BindingSets, String> {
        Err("kmer_similarity: aggregate-finish not implemented".into())
    }

    fn cardinality_estimate(_input: Cardinality, _args: Vec<Value>) -> Result<Cardinality, String> {
        Ok(Cardinality {
            value: 1.0,
            accuracy: Accuracy::Accurate,
        })
    }

    fn doc() -> BindingSets {
        BindingSets {
            vars: vec!["doc".into()],
            rows: vec![vec![Binding {
                name: "doc".into(),
                value: string_literal(
                    "kmer_similarity(query_seq, subject_seq) -> similarity: \
                     Jaccard coefficient over k=3 amino-acid k-mers. \
                     A simple, WASM-friendly stand-in for the BLAST procedure \
                     in the SCRY paper (Stringer et al., CEUR Vol-1795, paper 30).",
                ),
            }]],
        }
    }
}

export!(Component);
