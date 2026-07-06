package ai.tegmentum.scry;

import org.eclipse.rdf4j.query.BindingSet;
import org.eclipse.rdf4j.query.TupleQueryResult;
import org.eclipse.rdf4j.repository.RepositoryConnection;
import org.eclipse.rdf4j.repository.sail.SailRepository;
import org.eclipse.rdf4j.rio.RDFFormat;
import org.eclipse.rdf4j.sail.memory.MemoryStore;

import java.io.InputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.util.Objects;

/**
 * Reproduces the bioinformatics use case from Stringer et al. (SCRY, CEUR
 * Vol-1795, paper 30) using the tegmentum rdf4j-webfunction-plugin. Given
 * hemoglobin subunit beta as the query protein, computes k-mer Jaccard
 * similarity against every other protein in the graph via the WASM
 * component, then joins with tissue expression data to count co-expressed
 * tissues per homolog.
 *
 * <p>The paper's version federated to a separate Python SCRY endpoint over
 * SPARQL SERVICE; this version runs the same shape entirely in-process, with
 * a sandboxed WebAssembly component invoked via {@code wf:call}. No HTTP,
 * no separate process, one .wasm binary portable across Stardog / Jena /
 * RDF4J unchanged.
 */
public final class BioDemo {

    private static final String QUERY_PROTEIN = "<http://purl.uniprot.org/uniprot/P68871>";
    // BLASTP bit-score threshold. ~20 is weakly significant, ~50 strongly. The
    // demo default of 20 keeps distant hemoglobin homologs visible; override
    // via -Dwf.threshold=... to raise the bar.
    private static final String THRESHOLD_LITERAL =
            "\"" + System.getProperty("wf.threshold", "20")
                    + "\"^^<http://www.w3.org/2001/XMLSchema#decimal>";

    public static void main(final String[] args) throws Exception {
        final Path wasm = locateWasm();
        final String query = renderQuery(wasm);

        final SailRepository repo = new SailRepository(new MemoryStore());
        repo.init();
        try {
            loadData(repo);
            runQuery(repo, query, System.out);
        } finally {
            repo.shutDown();
        }
    }

    private static Path locateWasm() {
        // -Dwf.blastp.wasm=... points at whatever the caller built. Falls
        // back to the in-tree cargo-component output.
        final String override = System.getProperty("wf.blastp.wasm");
        if (override != null && !override.isBlank()) {
            return Path.of(override).toAbsolutePath();
        }
        return Path.of("src/main/rust/blastp/target/wasm32-wasip1/release/blastp.wasm")
                .toAbsolutePath();
    }

    private static String renderQuery(final Path wasm) throws Exception {
        try (InputStream in = Objects.requireNonNull(
                BioDemo.class.getResourceAsStream("/queries/coexpression.rq"),
                "coexpression.rq not on classpath")) {
            final String template = new String(in.readAllBytes(), StandardCharsets.UTF_8);
            return template
                    .replace("${QUERY_PROTEIN}", QUERY_PROTEIN)
                    .replace("${WASM_URL}", wasm.toUri().toString())
                    .replace("${THRESHOLD}", THRESHOLD_LITERAL);
        }
    }

    private static void loadData(final SailRepository repo) throws Exception {
        try (RepositoryConnection conn = repo.getConnection();
             InputStream in = Objects.requireNonNull(
                     BioDemo.class.getResourceAsStream("/data/proteins.ttl"),
                     "proteins.ttl not on classpath")) {
            conn.add(in, "", RDFFormat.TURTLE);
        }
    }

    private static void runQuery(final SailRepository repo,
                                 final String query,
                                 final PrintStream out) {
        try (RepositoryConnection conn = repo.getConnection();
             TupleQueryResult result = conn.prepareTupleQuery(query).evaluate()) {
            out.printf("%-42s %-32s %-12s %s%n",
                    "homolog", "name", "bit-score", "co-expressed tissues");
            out.println("-".repeat(110));
            int rows = 0;
            while (result.hasNext()) {
                final BindingSet row = result.next();
                out.printf("%-42s %-32s %-12s %s%n",
                        row.getValue("homolog").stringValue(),
                        row.getValue("homologName").stringValue(),
                        row.getValue("bitScore").stringValue(),
                        row.getValue("coexpressedTissueCount").stringValue());
                rows++;
            }
            if (rows == 0) {
                out.println("(no homologs above the similarity threshold with any shared tissue)");
            }
        }
    }
}
