# scry-webfunctions-demo

A WebAssembly-component reimplementation of the bioinformatics use case from
Stringer, Meroño-Peñuela, Abeln, van Harmelen & Heringa,
[**SCRY: extending SPARQL with custom data processing methods for the life sciences**](https://ceur-ws.org/Vol-1795/paper30.pdf)
(SWAT4LS 2016, CEUR Vol-1795, paper 30).

The paper's SCRY is a Python SPARQL endpoint you federate into via `SERVICE`.
You register Python procedures under URIs, and at query time SCRY generates
an RDF subgraph on the fly by executing them. Their headline demo runs
[BLAST](https://en.wikipedia.org/wiki/BLAST_(biotechnology))
against hemoglobin β to find homologs, then joins with tissue expression
data from the Human Protein Atlas to count co-expressed tissues.

**This repo reproduces the same query with WebAssembly components in place
of Python.** No separate endpoint, no HTTP hop, one sandboxed `.wasm` that
runs unchanged under [Stardog](https://github.com/tegmentum/stardog-webfunction-plugin),
[Apache Jena](https://github.com/tegmentum/jena-webfunction-plugin), and
[Eclipse RDF4J](https://github.com/tegmentum/rdf4j-webfunction-plugin). The
runner uses RDF4J's in-JVM `MemoryStore` for the shortest path from
`mvn exec:java` to output.

## The bio component family

Two components ship in this repo, both callable from any of the three
webfunction plugins:

| Component | Role | SPARQL surface |
|---|---|---|
| `blastp.wasm` | BLASTP-shape local alignment: `bio::alignment::pairwise::Aligner` + BLOSUM62 + Karlin-Altschul stats. Returns bit score, raw score, E-value, percent identity, alignment length. | `wf:call` filter form for the bit score, tuple form for all five outputs |
| `protparam.wasm` | Expasy-ProtParam-shape biochemical properties: average MW, Bjellqvist pI via bisection, Kyte-Doolittle GRAVY, length. | Tuple form so all four outputs land in one query |

Together they let a single SPARQL query rank homologs by sequence
similarity **and** characterise their biochemistry — the classic
paper-figure workflow. See `src/main/resources/queries/` for both.

## Actual BLASTP, not an approximation

- Gapped local alignment via [rust-bio](https://crates.io/crates/bio)'s
  Smith-Waterman aligner (`bio::alignment::pairwise::Aligner`).
- Standard **BLOSUM62** substitution matrix
  (`bio::scores::blosum62`).
- NCBI BLASTP default affine gap penalties: open=-11, extend=-1.
- Karlin-Altschul statistics with canonical BLOSUM62/11/1 constants
  (K=0.041, λ=0.267) → bit score and E-value on the exact same scale
  NCBI blastp reports.

BLAST is a *heuristic* for exactly this local-alignment problem, useful
when scanning gigabyte databases. For candidate sets small enough to fit
in a SPARQL query (dozens to thousands of sequences), running the exact
algorithm is both faster than BLAST-the-heuristic and *more accurate*.
The component is 130 lines of Rust plus rust-bio's aligner — reusable
outside this demo by anyone who has two amino-acid sequences and wants
BLASTP-shape output.

A second component, `kmer_similarity`, is included as a naive baseline
(k=3 Jaccard) to show the algorithm choice matters — it drops close
homologs like α that BLAST catches clean.

## The query

```sparql
PREFIX wf:   <http://tegmentum.ai/ns/webfunction/>
PREFIX bio:  <http://tegmentum.ai/ns/scry-demo/bio/>
PREFIX up:   <http://purl.uniprot.org/uniprot/>
PREFIX xsd:  <http://www.w3.org/2001/XMLSchema#>

SELECT ?homolog ?homologName ?bitScore
       (COUNT(DISTINCT ?sharedTissue) AS ?coexpressedTissueCount)
WHERE {
    up:P68871 bio:sequence ?querySeq .                # hemoglobin β
    ?homolog a bio:Protein ;
             rdfs:label ?homologName ;
             bio:sequence ?candidateSeq .
    FILTER (?homolog != up:P68871)

    # WebAssembly component — BLASTP bit score in the same JVM.
    BIND (xsd:decimal(wf:call(<file:.../blastp.wasm>,
                              ?querySeq, ?candidateSeq)) AS ?bitScore)
    FILTER (?bitScore >= "20"^^xsd:decimal)

    up:P68871 bio:expressedIn ?sharedTissue .
    ?homolog  bio:expressedIn ?sharedTissue .
}
GROUP BY ?homolog ?homologName ?bitScore
ORDER BY DESC(?bitScore)
```

Sample output on the fetched UniProt sequences + HPA tissue expression:

```
=== BLASTP hemoglobin-β homologs + tissue coexpression ===
homolog                                    name                             bit-score    co-expressed tissues
--------------------------------------------------------------------------------------------------------------
http://purl.uniprot.org/uniprot/P02042     Hemoglobin subunit delta         284.65       1
http://purl.uniprot.org/uniprot/P69891     Hemoglobin subunit gamma-1       229.56       1
http://purl.uniprot.org/uniprot/P69905     Hemoglobin subunit alpha         114.39       34
http://purl.uniprot.org/uniprot/P02144     Myoglobin                        46.98        6

=== Biochemical properties (protparam) ===
protein                            length      MW (Da)       pI    GRAVY
--------------------------------------------------------------------------------
Hemoglobin subunit alpha              142     15257.55     8.73     0.05
Hemoglobin subunit beta               147     15998.41     6.82     0.01
Hemoglobin subunit delta              147     16055.48     7.99    -0.05
Hemoglobin subunit gamma-1            147     16128.41     6.72    -0.12
Insulin (preproinsulin)               110     11980.91     5.22     0.19
Myoglobin                             154     17183.81     7.29    -0.48
```

Interpretation matches the biology:
- **δ (284 bits, 1 tissue)** — near-identical to β at the sequence level
  (~93% identity), expressed above nTPM 100 only in bone marrow.
- **γ-1 (229 bits, 1 tissue)** — the *fetal* analog of β; HPA shows it
  above nTPM 100 only in placenta, correctly the one tissue it shares
  with β's own placental expression.
- **α (114 bits, 34 tissues)** — the tetramer partner of β. α and β have
  only ~45% sequence identity but co-assemble to form adult hemoglobin,
  so RBCs — which circulate everywhere — carry both wherever they go.
  BLASTP catches the sequence relationship cleanly; the naive k-mer
  baseline scores it at 0.06 Jaccard and misses it entirely.
- **Myoglobin (47 bits, 6 tissues)** — weakly significant hit; shares the
  globin fold, muscle-tissue overlap driven by residual blood.

Insulin (pancreas-only) drops out — no shared tissues with β.

## Why webfunctions instead of SCRY?

SCRY set four requirements: generality, reusability, interoperability,
scalability. Compared point-by-point:

| Requirement | SCRY (Python endpoint) | Webfunctions (Wasm components) |
|---|---|---|
| Generality | Any Python + any pip package | Any language that compiles to Wasm Component Model (Rust, C, Go, Zig, JS via ComponentizeJS, …) |
| Reusability | Share Python source or a running endpoint | Ship a `.wasm` binary — reproducible, hash-addressable, no build required |
| Interoperability | Any SPARQL engine that supports `SERVICE` federation | Any of the three engines below via the same plugin, no federation needed |
| Scalability | Depends on Python + separate HTTP endpoint | No process, no HTTP; direct in-JVM invocation with wasmtime cranelift compilation |
| Privacy | Query text and data cross the wire to SCRY | Component runs in the same JVM — nothing leaves the process boundary |
| Sandboxing | Python: none by default (arbitrary code) | Wasm: capability-based, memory-safe, no filesystem/network unless granted |
| Determinism | Python: none guaranteed | Wasm: pure computation is bit-for-bit reproducible |

The one thing SCRY has that this doesn't: SCRY procedures can shell out
to existing native binaries (e.g. the actual `blastp` binary). Webfunctions
components can't do that without WASI process spawning, which isn't
stable yet. So instead of wrapping `blastp` we bring an alignment
implementation into the sandbox — rust-bio's Smith-Waterman, which
implements the exact problem BLAST heuristically approximates.

## Running

Requires: Java 21, Maven, `cargo`, `cargo component`.

```bash
# 1. Build the wasm components.
for c in blastp protparam kmer_similarity; do
    ( cd "src/main/rust/$c" && cargo component build --release )
done

# 2. Run the demo. This uses the tegmentum/rdf4j-webfunction-plugin, so
#    that repo needs to be `mvn install`-ed into your local ~/.m2 first;
#    the demo pom depends on ai.tegmentum.rdf4j:webfunction:0.1.0-SNAPSHOT.
mvn exec:java

# Optional: raise the bit-score threshold (~50 = strong, ~100 = highly).
mvn exec:java -Dwf.threshold=50

# Optional: run the naive baseline component instead.
mvn exec:java -Dwf.blastp.wasm=src/main/rust/kmer_similarity/target/wasm32-wasip1/release/kmer_similarity.wasm -Dwf.threshold=0

# Optional: override the wasm component memory ceiling (default 64 MB, set
# via the wf.memory.max.bytes project property; forwarded to the plugin as
# webfunctions.memory.max.bytes). blastp.wasm's minimum linear memory is
# ~1.1 MB and grows during alignment, so the wasmtime built-in default is
# not enough — the pom default handles this. Raise for larger inputs.
mvn exec:java -Dwf.memory.max.bytes=134217728   # 128 MB
```

## Layout

```
src/main/
├── java/ai/tegmentum/scry/BioDemo.java    # RDF4J runner
├── resources/
│   ├── data/proteins.ttl                  # Synthetic hemoglobin + tissue graph
│   └── queries/coexpression.rq            # The paper's Figure-3 query, adapted
└── rust/
    ├── blastp/                            # rust-bio Smith-Waterman + BLOSUM62
    │   ├── src/lib.rs
    │   └── wit/webfunction.wit
    ├── protparam/                         # MW / pI / GRAVY / length
    │   ├── src/lib.rs
    │   └── wit/webfunction.wit
    └── kmer_similarity/                   # Naive k=3 Jaccard baseline
        ├── src/lib.rs
        └── wit/webfunction.wit
```

## Data

`proteins.ttl` is regenerated by `scripts/build-data.sh`, which:

- Fetches canonical UniProt FASTA sequences from `rest.uniprot.org` (CC BY 4.0).
- Downloads the HPA `rna_tissue_consensus.tsv.zip` from
  proteinatlas.org (CC BY-SA 3.0) and filters each protein's tissue
  expression to nTPM ≥ 100 (roughly the HPA "elevated" cutoff — this
  keeps the query's coexpression counts biologically meaningful rather
  than dominated by circulation-level background).

The generated file is committed so `mvn exec:java` works offline.
Re-run the script to refresh against upstream.

## Related

- Paper: [ceur-ws.org/Vol-1795/paper30.pdf](https://ceur-ws.org/Vol-1795/paper30.pdf)
- The plugin family:
  [tegmentum/stardog-webfunction-plugin](https://github.com/tegmentum/stardog-webfunction-plugin),
  [tegmentum/jena-webfunction-plugin](https://github.com/tegmentum/jena-webfunction-plugin),
  [tegmentum/rdf4j-webfunction-plugin](https://github.com/tegmentum/rdf4j-webfunction-plugin)
- Wasm runtime: [tegmentum/webassembly4j](https://github.com/tegmentum/webassembly4j) (wasmtime provider)
- Alignment library: [rust-bio](https://crates.io/crates/bio)
