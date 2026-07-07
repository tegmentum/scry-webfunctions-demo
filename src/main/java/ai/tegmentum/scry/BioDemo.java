package ai.tegmentum.scry;

import ai.tegmentum.rdf4j.webfunctions.WfEvaluationStrategyFactory;

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
 * Vol-1795, paper 30) using the tegmentum rdf4j-webfunction-plugin. Two
 * demos in one run:
 *
 * <ol>
 *   <li>{@code coexpression.rq} — BLAST hemoglobin β against every other
 *       protein via {@code blastp.wasm} (rust-bio Smith-Waterman + BLOSUM62 +
 *       Karlin-Altschul stats) and count tissues shared with each homolog.
 *       Uses the {@code wf:call} filter form (single-value output).</li>
 *   <li>{@code protparam.rq} — Compute biochemical properties (length, MW,
 *       theoretical pI, GRAVY) for every protein via {@code protparam.wasm}.
 *       Uses the {@code wf:call} tuple-function form (multi-var output)
 *       via {@link WfEvaluationStrategyFactory}.</li>
 * </ol>
 *
 * <p>The paper's version federated to a separate Python SCRY endpoint over
 * SPARQL SERVICE; this version runs the same shape entirely in-process, with
 * sandboxed WebAssembly components invoked via {@code wf:call}. No HTTP,
 * no separate process, one {@code .wasm} binary per procedure portable
 * across Stardog / Jena / RDF4J unchanged.
 */
public final class BioDemo {

    private static final String QUERY_PROTEIN = "<http://purl.uniprot.org/uniprot/P68871>";
    private static final String THRESHOLD_LITERAL =
            "\"" + System.getProperty("wf.threshold", "20")
                    + "\"^^<http://www.w3.org/2001/XMLSchema#decimal>";

    public static void main(final String[] args) throws Exception {
        final SailRepository repo = openRepo();
        try {
            loadData(repo);
            runCoexpression(repo, System.out);
            System.out.println();
            runProtparam(repo, System.out);
        } finally {
            repo.shutDown();
        }
    }

    /**
     * MemoryStore configured with {@link WfEvaluationStrategyFactory} so both
     * the {@code wf:call} filter form (delegates to the standard
     * FunctionRegistry) and the tuple-function form
     * {@code (args) wf:call (?vars)} work in the same session.
     */
    private static SailRepository openRepo() {
        final MemoryStore store = new MemoryStore();
        store.setEvaluationStrategyFactory(new WfEvaluationStrategyFactory(null));
        final SailRepository repo = new SailRepository(store);
        repo.init();
        return repo;
    }

    private static void loadData(final SailRepository repo) throws Exception {
        try (RepositoryConnection conn = repo.getConnection();
             InputStream in = Objects.requireNonNull(
                     BioDemo.class.getResourceAsStream("/data/proteins.ttl"),
                     "proteins.ttl not on classpath")) {
            conn.add(in, "", RDFFormat.TURTLE);
        }
    }

    private static void runCoexpression(final SailRepository repo, final PrintStream out) throws Exception {
        final Path wasm = locateWasm("wf.blastp.wasm",
                "src/main/rust/blastp/target/wasm32-wasip1/release/blastp.wasm");
        final String query = renderQuery("/queries/coexpression.rq", wasm)
                .replace("${QUERY_PROTEIN}", QUERY_PROTEIN)
                .replace("${THRESHOLD}", THRESHOLD_LITERAL);

        out.println("=== BLASTP hemoglobin-β homologs + tissue coexpression ===");
        out.printf("%-42s %-32s %-12s %s%n",
                "homolog", "name", "bit-score", "co-expressed tissues");
        out.println("-".repeat(110));
        try (RepositoryConnection conn = repo.getConnection();
             TupleQueryResult r = conn.prepareTupleQuery(query).evaluate()) {
            int rows = 0;
            while (r.hasNext()) {
                final BindingSet row = r.next();
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

    private static void runProtparam(final SailRepository repo, final PrintStream out) throws Exception {
        final Path wasm = locateWasm("wf.protparam.wasm",
                "src/main/rust/protparam/target/wasm32-wasip1/release/protparam.wasm");
        final String query = renderQuery("/queries/protparam.rq", wasm);

        out.println("=== Biochemical properties (protparam) ===");
        out.printf("%-32s %8s %12s %8s %8s%n",
                "protein", "length", "MW (Da)", "pI", "GRAVY");
        out.println("-".repeat(80));
        try (RepositoryConnection conn = repo.getConnection();
             TupleQueryResult r = conn.prepareTupleQuery(query).evaluate()) {
            while (r.hasNext()) {
                final BindingSet row = r.next();
                out.printf("%-32s %8s %12s %8s %8s%n",
                        row.getValue("name").stringValue(),
                        row.getValue("length").stringValue(),
                        row.getValue("molecularWeight").stringValue(),
                        row.getValue("theoreticalPi").stringValue(),
                        row.getValue("gravy").stringValue());
            }
        }
    }

    private static Path locateWasm(final String propertyKey, final String defaultRelative) {
        final String override = System.getProperty(propertyKey);
        if (override != null && !override.isBlank()) {
            return Path.of(override).toAbsolutePath();
        }
        return Path.of(defaultRelative).toAbsolutePath();
    }

    private static String renderQuery(final String resourcePath, final Path wasm) throws Exception {
        try (InputStream in = Objects.requireNonNull(
                BioDemo.class.getResourceAsStream(resourcePath),
                resourcePath + " not on classpath")) {
            final String template = new String(in.readAllBytes(), StandardCharsets.UTF_8);
            return template.replace("${WASM_URL}", wasm.toUri().toString());
        }
    }
}
