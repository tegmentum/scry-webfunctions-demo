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

wit_bindgen::generate!({
    world: "webfunction",
    path: "wit",
});

use bio::alignment::pairwise::Aligner;
use bio::scores::blosum62;

use stardog::webfunction::types::{Accuracy, Binding, Literal};

struct Component;

const XSD_STRING: &str = "http://www.w3.org/2001/XMLSchema#string";
const XSD_DECIMAL: &str = "http://www.w3.org/2001/XMLSchema#decimal";
const XSD_INTEGER: &str = "http://www.w3.org/2001/XMLSchema#integer";

// NCBI BLASTP defaults: open=-11, extend=-1, BLOSUM62.
const GAP_OPEN: i32 = -11;
const GAP_EXTEND: i32 = -1;

// Karlin-Altschul parameters for BLOSUM62 with 11/1 affine gap penalties.
// Source: NCBI BLAST parameter tables (published statistics for this matrix +
// gap combination). Values are canonical constants, not fitted here.
const K: f64 = 0.041;
const LAMBDA: f64 = 0.267;

fn string_literal(label: &str) -> Value {
    Value::Literal(Literal {
        label: label.into(),
        datatype: XSD_STRING.into(),
        lang: None,
    })
}

fn decimal_literal(v: f64) -> Value {
    Value::Literal(Literal {
        label: format!("{:.4}", v),
        datatype: XSD_DECIMAL.into(),
        lang: None,
    })
}

fn integer_literal(v: i64) -> Value {
    Value::Literal(Literal {
        label: v.to_string(),
        datatype: XSD_INTEGER.into(),
        lang: None,
    })
}

fn scientific_literal(v: f64) -> Value {
    // E-values span many orders of magnitude — scientific notation is what
    // NCBI blastp reports and what humans expect.
    Value::Literal(Literal {
        label: format!("{:.2e}", v),
        datatype: XSD_DECIMAL.into(),
        lang: None,
    })
}

fn sequence_of(arg: &Value) -> Result<&str, String> {
    match arg {
        Value::Literal(literal) => Ok(literal.label.as_str()),
        _ => Err("blastp: each argument must be a literal sequence".into()),
    }
}

/// Percent identity across the aligned region, from a rust-bio alignment path.
fn percent_identity(alignment: &bio::alignment::Alignment, x: &[u8], y: &[u8]) -> f64 {
    let mut matches = 0usize;
    let mut length = 0usize;
    let mut xi = alignment.xstart;
    let mut yi = alignment.ystart;
    for op in &alignment.operations {
        match op {
            bio::alignment::AlignmentOperation::Match => {
                matches += 1;
                length += 1;
                xi += 1;
                yi += 1;
            }
            bio::alignment::AlignmentOperation::Subst => {
                length += 1;
                xi += 1;
                yi += 1;
            }
            bio::alignment::AlignmentOperation::Del => {
                length += 1;
                yi += 1;
            }
            bio::alignment::AlignmentOperation::Ins => {
                length += 1;
                xi += 1;
            }
            _ => {}
        }
    }
    // Guard against zero-length alignment (both empty inputs).
    let _ = (xi, yi); // silence unused warnings
    let _ = x;
    let _ = y;
    if length == 0 { 0.0 } else { (matches as f64 / length as f64) * 100.0 }
}

impl Guest for Component {
    fn evaluate(args: Vec<Value>) -> Result<BindingSets, String> {
        if args.len() != 2 {
            return Err(format!(
                "blastp: expected 2 args (query_seq, subject_seq), got {}",
                args.len()
            ));
        }
        let query = sequence_of(&args[0])?.as_bytes();
        let subject = sequence_of(&args[1])?.as_bytes();

        // Smith-Waterman local alignment via rust-bio, with BLOSUM62 substitution
        // scores and NCBI blastp default gap penalties.
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
        let identity_pct = percent_identity(&alignment, query, subject);
        let align_len = alignment.operations.len() as i64;

        Ok(BindingSets {
            vars: vec![
                "bit_score".into(),
                "raw_score".into(),
                "e_value".into(),
                "identity".into(),
                "align_length".into(),
            ],
            rows: vec![vec![
                Binding { name: "bit_score".into(),    value: decimal_literal(bit_score) },
                Binding { name: "raw_score".into(),    value: decimal_literal(raw_score) },
                Binding { name: "e_value".into(),      value: scientific_literal(e_value) },
                Binding { name: "identity".into(),     value: decimal_literal(identity_pct) },
                Binding { name: "align_length".into(), value: integer_literal(align_len) },
            ]],
        })
    }

    fn aggregate_step(_args: Vec<Value>, _mult: u64) -> Result<(), String> {
        Err("blastp: aggregate-step not implemented".into())
    }

    fn aggregate_finish() -> Result<BindingSets, String> {
        Err("blastp: aggregate-finish not implemented".into())
    }

    fn cardinality_estimate(_input: Cardinality, _args: Vec<Value>) -> Result<Cardinality, String> {
        Ok(Cardinality { value: 1.0, accuracy: Accuracy::Accurate })
    }

    fn doc() -> BindingSets {
        BindingSets {
            vars: vec!["doc".into()],
            rows: vec![vec![Binding {
                name: "doc".into(),
                value: string_literal(
                    "blastp(query_seq, subject_seq) -> (bit_score, raw_score, e_value, \
                     identity, align_length). Wraps rust-bio's Smith-Waterman aligner \
                     with BLOSUM62 + NCBI blastp default gap penalties (open=-11, \
                     extend=-1). Karlin-Altschul statistics use canonical BLOSUM62/11/1 \
                     constants (K=0.041, λ=0.267). Bit score is the primary output for \
                     ranking; retrieve as the first column when used via wf:call filter form.",
                ),
            }]],
        }
    }
}

export!(Component);
