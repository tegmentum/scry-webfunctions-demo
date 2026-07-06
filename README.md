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
runner in this repo uses RDF4J's in-JVM `MemoryStore` for the shortest path
from `mvn exec:java` to output.

## The query

```sparql
PREFIX wf:   <http://tegmentum.ai/ns/webfunction/>
PREFIX bio:  <http://tegmentum.ai/ns/scry-demo/bio/>
PREFIX up:   <http://purl.uniprot.org/uniprot/>
PREFIX xsd:  <http://www.w3.org/2001/XMLSchema#>

SELECT ?homolog ?homologName ?similarity
       (COUNT(DISTINCT ?sharedTissue) AS ?coexpressedTissueCount)
WHERE {
    up:P68871 bio:sequence ?querySeq .            # hemoglobin β
    ?homolog a bio:Protein ;
             rdfs:label ?homologName ;
             bio:sequence ?candidateSeq .
    FILTER (?homolog != up:P68871)

    # WebAssembly component — k-mer Jaccard similarity in the same JVM.
    BIND (xsd:decimal(wf:call(<file:.../kmer_similarity.wasm>,
                              ?querySeq, ?candidateSeq)) AS ?similarity)
    FILTER (?similarity >= 0.0)

    up:P68871 bio:expressedIn ?sharedTissue .
    ?homolog  bio:expressedIn ?sharedTissue .
}
GROUP BY ?homolog ?homologName ?similarity
ORDER BY DESC(?similarity)
```

Sample output on the shipped synthetic dataset:

```
homolog                                    name                             similarity   co-expressed tissues
--------------------------------------------------------------------------------------------------------------
http://purl.uniprot.org/uniprot/P02042     Hemoglobin subunit delta         0.7256       1
http://purl.uniprot.org/uniprot/P69891     Hemoglobin subunit gamma-1       0.3028       2
http://purl.uniprot.org/uniprot/P69905     Hemoglobin subunit alpha         0.0606       2
```

Interpretation: delta is β's near-clone (99% conserved in the region we
score), γ-1 is the fetal analog, and α forms adult hemoglobin
heterotetramers with β so they must be co-expressed — but α diverges
enough in sequence that raw 3-mer Jaccard barely picks it up. That last
point is the paper's exact argument for embedding procedures at query
time: use whatever scoring you need, not what your triplestore's built-in
functions provide. Swap `kmer_similarity` for a real BLAST-in-Rust
component (or, in the paper's version, a Python wrapper around
`blastp`) and α climbs to the top.

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

The one thing SCRY has that this doesn't: SCRY procedures can shell out to
existing native binaries (e.g. the actual `blastp` binary). Webfunctions
components can't do that without WASI process spawning, which isn't
stable yet. For the common case — pure computation, or a library that
compiles to Wasm — webfunctions wins on every other axis.

## Running

Requires: Java 21, Maven, `cargo`, `cargo component`.

```bash
# 1. Build the wasm component.
cd src/main/rust/kmer_similarity
cargo component build --release
cd -

# 2. Run the demo. This uses the tegmentum/rdf4j-webfunction-plugin, so
#    that repo needs to be `mvn install`-ed into your local ~/.m2 first;
#    the demo pom depends on ai.tegmentum.rdf4j:webfunction:0.1.0-SNAPSHOT.
mvn exec:java

# Optional: threshold the results.
mvn exec:java -Dwf.threshold=0.15
```

## Layout

```
src/main/
├── java/ai/tegmentum/scry/BioDemo.java   # RDF4J runner
├── resources/
│   ├── data/proteins.ttl                 # Synthetic hemoglobin + tissue graph
│   └── queries/coexpression.rq           # The paper's Figure-3 query, adapted
└── rust/kmer_similarity/                 # k=3 Jaccard component (cargo component)
    ├── src/lib.rs
    └── wit/webfunction.wit               # stardog:webfunction@0.2.0
```

## Dataset caveats

`proteins.ttl` holds real UniProt canonical sequences (P68871, P69905,
P02042, P69891, P02144, P01308) plus hand-curated tissue-expression
triples that approximate the Human Protein Atlas RNA-expression evidence.
The tissue set is deliberately small so the query's result table is
reproducible and human-readable. **Do not use as a substitute for HPA.**

## Related

- Paper: [ceur-ws.org/Vol-1795/paper30.pdf](https://ceur-ws.org/Vol-1795/paper30.pdf)
- The plugin family:
  [tegmentum/stardog-webfunction-plugin](https://github.com/tegmentum/stardog-webfunction-plugin),
  [tegmentum/jena-webfunction-plugin](https://github.com/tegmentum/jena-webfunction-plugin),
  [tegmentum/rdf4j-webfunction-plugin](https://github.com/tegmentum/rdf4j-webfunction-plugin)
- Wasm runtime: [tegmentum/webassembly4j](https://github.com/tegmentum/webassembly4j) (wasmtime provider)
