package ro.licenta.genomicsapi.service;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.http.MediaType;
import org.springframework.stereotype.Service;
import org.springframework.web.reactive.function.client.WebClient;

import java.time.Duration;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * ChatbotService — orchestrează comunicarea cu Ollama (LLM local).
 *
 * Triplu strat de protecție pentru a restricționa la domeniul geneticii:
 *  1. Validare cuvinte cheie (în Java, înainte de Ollama)
 *  2. System prompt strict (în prompt-ul către LLM)
 *  3. Limitare lungime răspuns (max 500 tokens)
 *
 * Endpoint Ollama: http://localhost:11434/api/generate
 */
@Service
public class ChatbotService {

    private static final Logger log = LoggerFactory.getLogger(ChatbotService.class);

    private final WebClient ollamaClient;
    private final String model;

    /**
     * Cuvinte cheie din domeniul geneticii medicale.
     * Întrebarea trebuie să conțină măcar unul ca să fie procesată.
     */
    private static final Set<String> GENETICS_KEYWORDS = Set.of(
            // Genetică de bază
            "chromosome", "chromosomes", "cromozom", "cromozomi",
            "gene", "genes", "gena", "genetic", "genetica",
            "dna", "adn", "rna", "arn", "nucleotide", "nucleotida",
            "allele", "alela", "genome", "genom", "genotype", "genotip",
            "phenotype", "fenotip", "exon", "intron", "promoter",
            // Variante și mutații
            "variant", "varianta", "variante", "mutation", "mutatie",
            "snp", "snv", "indel", "insertion", "deletion",
            "frameshift", "missense", "nonsense", "synonymous",
            "splice", "splicing",
            // Clasificare clinică ClinVar
            "pathogenic", "patogen", "patogenic", "patogenica",
            "benign", "benigna", "vus", "uncertain", "incert",
            "likely", "probabil",
            // Baze de date
            "clinvar", "gnomad", "omim", "ensembl", "dbsnp", "cosmic",
            "ucsc", "ncbi",
            // Bioinformatică
            "bam", "vcf", "fasta", "fastq", "depth", "coverage",
            "allele frequency", "af", "read", "reads", "mapping",
            "alignment", "samtools", "pysam",
            // Tool-uri
            "vep", "sift", "polyphen", "cadd", "deepvariant", "gatk",
            "haplotypecaller",
            // Moștenire
            "heterozygous", "homozygous", "heterozigot", "homozigot",
            "carrier", "purtator", "autosomal", "autozomal",
            "recessive", "recesiv", "dominant", "x-linked",
            "inheritance", "ereditar",
            // Gene celebre (exemple)
            "brca1", "brca2", "cftr", "tp53", "apoe", "huntingtin",
            "pnpla3", "sbf1", "dystrophin",
            // Boli genetice comune
            "cancer", "tumor", "tumora", "fibrosis", "fibroza",
            "huntington", "thalassemia", "talasemie",
            "hemophilia", "hemofilie", "sickle",
            // Concepte avansate
            "penetrance", "penetranta", "expressivity",
            "trio", "proband", "pedigree",
            "exome", "exom", "wgs", "wes", "ngs",
            "sequencing", "secventiere",
            // Termeni aplicație
            "report", "raport", "analysis", "analiza", "result", "rezultat"
    );

    /**
     * System prompt strict — restricționează LLM-ul la domeniul nostru.
     */
    private static final String SYSTEM_PROMPT = """
            You are a medical genetics assistant integrated into a clinical variant calling application.
            
            YOU MUST ONLY answer questions about:
            - Genetics: chromosomes, genes, alleles, mutations, DNA, RNA, genome
            - Variant classification: Pathogenic, Likely Pathogenic, VUS (Variant of Uncertain Significance), Benign, Likely Benign
            - Bioinformatics terms: BAM, VCF, FASTA, depth, allele frequency (AF), read coverage, mapping quality
            - Genetic databases: ClinVar, gnomAD, dbSNP, OMIM, Ensembl, COSMIC
            - Inheritance patterns: autosomal dominant/recessive, X-linked, carrier, heterozygous, homozygous
            - Disease-gene associations (general educational information only)
            - Variant Effect Predictor (VEP), SIFT, PolyPhen-2, CADD predictions
            - Tools: DeepVariant, GATK, samtools, pysam
            
            STRICT RULES:
            1. If the user asks ANYTHING outside genetics (weather, recipes, politics, personal questions, etc.),
               respond ONLY with: "I can only help with questions about genetics and variant interpretation. Please ask me about genetic concepts, variant classification, or bioinformatics terminology."
            2. NEVER provide medical diagnosis or treatment advice.
            3. ALWAYS recommend consulting a genetic counselor or physician for personal medical decisions.
            4. Keep responses concise: maximum 3 short paragraphs.
            5. Use clear, accessible language suitable for medical professionals.
            6. If unsure about a specific variant, recommend checking ClinVar directly.
            
            Answer the user's question following these rules strictly.
            """;

    public ChatbotService(
            @Value("${app.ollama.url:http://localhost:11434}") String ollamaUrl,
            @Value("${app.ollama.model:llama3.1:8b}") String model) {
        this.model = model;
        this.ollamaClient = WebClient.builder()
                .baseUrl(ollamaUrl)
                .codecs(c -> c.defaultCodecs().maxInMemorySize(10 * 1024 * 1024))
                .build();
        log.info("ChatbotService inițializat: url={}, model={}", ollamaUrl, model);
    }

    /**
     * Verifică dacă întrebarea conține cuvinte cheie din domeniul geneticii.
     */
    public boolean isInDomain(String question) {
        if (question == null || question.trim().isEmpty()) return false;
        String lower = question.toLowerCase()
                .replace("ă", "a").replace("â", "a").replace("î", "i")
                .replace("ș", "s").replace("ț", "t");
        return GENETICS_KEYWORDS.stream().anyMatch(lower::contains);
    }

    /**
     * Trimite întrebarea către Ollama și returnează răspunsul.
     */
    public String ask(String question) {
        if (!isInDomain(question)) {
            log.info("Întrebare în afara domeniului: {}", question);
            return "I can only help with questions about genetics and variant interpretation. " +
                    "Please ask me about genetic concepts, variant classification (Pathogenic, VUS, Benign), " +
                    "or bioinformatics terminology (BAM, VCF, ClinVar, etc.).";
        }

        log.info("Întrebare validă, trimit la Ollama: {}", question);
        Map<String, Object> request = new HashMap<>();
        request.put("model", model);
        request.put("prompt", SYSTEM_PROMPT + "\n\nUser question: " + question);
        request.put("stream", false);
        Map<String, Object> options = new HashMap<>();
        options.put("temperature", 0.3);
        options.put("num_predict", 500);
        request.put("options", options);

        try {
            long t0 = System.currentTimeMillis();
            Map<String, Object> response = ollamaClient.post()
                    .uri("/api/generate")
                    .contentType(MediaType.APPLICATION_JSON)
                    .bodyValue(request)
                    .retrieve()
                    .bodyToMono(Map.class)
                    .timeout(Duration.ofSeconds(60))
                    .block();

            long elapsed = (System.currentTimeMillis() - t0) / 1000;
            log.info("Ollama a răspuns în {}s", elapsed);

            if (response == null || !response.containsKey("response")) {
                return "Sorry, I couldn't generate a response. Please try again.";
            }

            return ((String) response.get("response")).trim();

        } catch (Exception e) {
            log.error("Eroare apel Ollama", e);
            return "Sorry, the assistant is currently unavailable. Please try again later.";
        }
    }
}
