package ro.licenta.genomicsapi.service;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.stereotype.Service;
import ro.licenta.genomicsapi.model.Variant;

import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.concurrent.TimeUnit;

/**
 * VepAnnotationService — apelează Ensembl VEP via Docker pentru anotarea
 * variantelor cu informații clinice (ClinVar, gnomAD, SIFT, PolyPhen).
 *
 * COMANDA DOCKER ECHIVALENTĂ:
 *   docker run --rm -v C:\vep_data:/data ensemblorg/ensembl-vep \
 *     vep --input_file /data/input/variants.vcf \
 *         --output_file /data/input/variants_annotated.vcf \
 *         --cache --dir_cache /data --offline \
 *         --fasta /data/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
 *         --assembly GRCh38 \
 *         --custom file=/data/plugins/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,fields=CLNSIG%CLNDN%CLNREVSTAT \
 *         --vcf --force_overwrite
 *
 * Output: VCF anotat parsat în List<Variant>.
 */
@Service
public class VepAnnotationService {

    private static final Logger log = LoggerFactory.getLogger(VepAnnotationService.class);

    @Value("${app.vep.data-dir:C:/vep_data}")
    private String vepDataDir;

    @Value("${app.vep.docker-image:ensemblorg/ensembl-vep}")
    private String dockerImage;

    @Value("${app.vep.timeout-minutes:30}")
    private int timeoutMinutes;

    /**
     * Anotează un fișier VCF cu VEP+ClinVar.
     *
     * @param inputVcfPath calea VCF de intrare (poate fi orice locație)
     * @return lista de variante cu anotări complete
     */
    public List<Variant> annotateVcf(String inputVcfPath) throws IOException, InterruptedException {
        log.info("Pornesc anotare VEP pentru: {}", inputVcfPath);

        // Folosim subfolderul /input din vep_data ca punct de intrare/ieșire pentru Docker
        Path vepInputDir = Paths.get(vepDataDir, "input");
        Files.createDirectories(vepInputDir);

        // Copiem VCF-ul de intrare în /input dacă nu e deja acolo
        Path inputVcf = Paths.get(inputVcfPath);
        Path stagedInput = vepInputDir.resolve(inputVcf.getFileName());
        if (!stagedInput.toAbsolutePath().equals(inputVcf.toAbsolutePath())) {
            Files.copy(inputVcf, stagedInput, StandardCopyOption.REPLACE_EXISTING);
            log.info("VCF copiat în zona Docker: {}", stagedInput);
        }

        String inputName = inputVcf.getFileName().toString();
        String outputName = inputName.replace(".vcf", "_annotated.vcf");
        Path outputVcf = vepInputDir.resolve(outputName);

        // Construim comanda docker run
        // Pe Windows, calea volumelor are format C:\vep_data:/data
        // Java ProcessBuilder folosește forward slash care funcționează pe Docker Desktop
        String volumeMount = vepDataDir.replace("\\", "/") + ":/data";

        List<String> command = new ArrayList<>();
        command.add("docker");
        command.add("run");
        command.add("--rm");
        command.add("-v");
        command.add(volumeMount);
        command.add(dockerImage);
        command.add("vep");
        command.add("--input_file");  command.add("/data/input/" + inputName);
        command.add("--output_file"); command.add("/data/input/" + outputName);
        command.add("--cache");
        command.add("--dir_cache");   command.add("/data");
        command.add("--offline");
        command.add("--fasta");
        command.add("/data/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz");
        command.add("--assembly");    command.add("GRCh38");
        command.add("--custom");
        command.add("file=/data/plugins/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,fields=CLNSIG%CLNDN%CLNREVSTAT");
        command.add("--vcf");
        command.add("--force_overwrite");
        command.add("--no_stats");      // skip generare HTML stats (mai rapid)
        command.add("--symbol");        // adaugă simbolul genei
        command.add("--biotype");       // tipul de transcript
        command.add("--hgvs");          // notația HGVS
        command.add("--sift");          command.add("b");   // both prediction + score
        command.add("--polyphen");      command.add("b");

        log.info("Rulez comanda Docker VEP...");
        long t0 = System.currentTimeMillis();

        ProcessBuilder pb = new ProcessBuilder(command);
        pb.redirectErrorStream(true);
        Process process = pb.start();

        // Citim output-ul în timp real (pentru log)
        StringBuilder outputLog = new StringBuilder();
        try (BufferedReader reader = new BufferedReader(
                new InputStreamReader(process.getInputStream()))) {
            String line;
            while ((line = reader.readLine()) != null) {
                outputLog.append(line).append("\n");
                if (line.contains("WARNING") || line.contains("ERROR")) {
                    log.warn("VEP: {}", line);
                }
            }
        }

        boolean finished = process.waitFor(timeoutMinutes, TimeUnit.MINUTES);
        if (!finished) {
            process.destroyForcibly();
            throw new RuntimeException("VEP timeout după " + timeoutMinutes + " minute");
        }

        int exitCode = process.exitValue();
        long durationSec = (System.currentTimeMillis() - t0) / 1000;

        if (exitCode != 0) {
            log.error("VEP a eșuat cu exit code {}. Output: {}", exitCode, outputLog);
            throw new RuntimeException("VEP a eșuat (exit code " + exitCode + ")");
        }

        log.info("VEP completat în {}s", durationSec);

        // Parsăm VCF-ul anotat
        List<Variant> variants = parseAnnotatedVcf(outputVcf);
        log.info("Parsate {} variante din VCF anotat", variants.size());

        return variants;
    }

    /**
     * Parsează un VCF anotat de VEP și extrage informațiile relevante
     * din câmpul CSQ.
     *
     * Format CSQ așteptat (definit în comanda --vcf):
     *   Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|
     *   INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|
     *   Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|
     *   SYMBOL_SOURCE|HGNC_ID|SOURCE|SIFT|PolyPhen|ClinVar|ClinVar_CLNSIG|
     *   ClinVar_CLNDN|ClinVar_CLNREVSTAT
     */
    private List<Variant> parseAnnotatedVcf(Path vcfPath) throws IOException {
        List<Variant> variants = new ArrayList<>();
        List<String> csqFields = null;

        try (BufferedReader reader = Files.newBufferedReader(vcfPath)) {
            String line;
            while ((line = reader.readLine()) != null) {
                // Header: extragem ordinea câmpurilor CSQ
                if (line.startsWith("##INFO=<ID=CSQ")) {
                    csqFields = extractCsqFieldOrder(line);
                    continue;
                }
                if (line.startsWith("#")) continue;

                // Linie de date
                Variant variant = parseVariantLine(line, csqFields);
                if (variant != null) {
                    variants.add(variant);
                }
            }
        }

        return variants;
    }

    /**
     * Extrage ordinea câmpurilor din header-ul ##INFO=<ID=CSQ Format: ...>
     */
    private List<String> extractCsqFieldOrder(String headerLine) {
        int formatIdx = headerLine.indexOf("Format: ");
        if (formatIdx == -1) return Collections.emptyList();
        String formatStr = headerLine.substring(formatIdx + 8);
        // Eliminăm "> de la final
        formatStr = formatStr.replaceAll("[\">]+$", "");
        return Arrays.asList(formatStr.split("\\|"));
    }

    /**
     * Parsează o linie de variantă din VCF anotat.
     * Format: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE
     */
    private Variant parseVariantLine(String line, List<String> csqFields) {
        String[] cols = line.split("\t");
        if (cols.length < 8) return null;

        Variant v = new Variant();
        v.setChrom(cols[0]);
        v.setPos(Integer.parseInt(cols[1]));
        v.setRef(cols[3]);
        v.setAlt(cols[4]);

        // Extragem datele din câmpul INFO
        String info = cols[7];
        Map<String, String> infoMap = parseInfoField(info);

        // Date predict CNN (din comanda noastră)
        if (infoMap.containsKey("AF"))      v.setAf(parseDoubleSafe(infoMap.get("AF")));
        if (infoMap.containsKey("DP"))      v.setDepth(parseIntSafe(infoMap.get("DP")));
        if (infoMap.containsKey("CONF"))    v.setConfidence(parseDoubleSafe(infoMap.get("CONF")));
        if (infoMap.containsKey("GT_PRED")) v.setPredictedClass(infoMap.get("GT_PRED"));

        // Câmpul CSQ — prima anotare (din mai multe, dacă există)
        if (infoMap.containsKey("CSQ") && csqFields != null && !csqFields.isEmpty()) {
            String csqValue = infoMap.get("CSQ").split(",")[0];  // primul transcript
            String[] csqVals = csqValue.split("\\|", -1);

            Map<String, String> csqMap = new HashMap<>();
            for (int i = 0; i < Math.min(csqFields.size(), csqVals.length); i++) {
                csqMap.put(csqFields.get(i), csqVals[i]);
            }

            v.setGeneSymbol(csqMap.getOrDefault("SYMBOL", ""));
            v.setGeneId(csqMap.getOrDefault("Gene", ""));
            v.setConsequence(csqMap.getOrDefault("Consequence", ""));
            v.setImpact(csqMap.getOrDefault("IMPACT", ""));
            v.setBiotype(csqMap.getOrDefault("BIOTYPE", ""));
            v.setHgvsc(csqMap.getOrDefault("HGVSc", ""));
            v.setHgvsp(csqMap.getOrDefault("HGVSp", ""));
            v.setSift(csqMap.getOrDefault("SIFT", ""));
            v.setPolyphen(csqMap.getOrDefault("PolyPhen", ""));
            v.setClinSig(csqMap.getOrDefault("ClinVar_CLNSIG", ""));
            v.setClinDisease(csqMap.getOrDefault("ClinVar_CLNDN", ""));
            v.setClinReviewStatus(csqMap.getOrDefault("ClinVar_CLNREVSTAT", ""));
        }

        // Calculăm clasificarea finală
        v.setFinalClassification(computeFinalClassification(v));

        return v;
    }

    private Map<String, String> parseInfoField(String info) {
        Map<String, String> map = new HashMap<>();
        for (String item : info.split(";")) {
            int eq = item.indexOf('=');
            if (eq > 0) {
                map.put(item.substring(0, eq), item.substring(eq + 1));
            } else {
                map.put(item, "");
            }
        }
        return map;
    }

    /**
     * Logică de clasificare finală pe baza ClinVar + predicții in-silico.
     *
     * Prioritate:
     *  1. ClinVar (dacă există) — sursa cea mai sigură
     *  2. SIFT + PolyPhen pentru variante necunoscute
     *  3. UNKNOWN dacă nu avem nicio informație
     */
    private String computeFinalClassification(Variant v) {
        // 1. ClinVar
        String cs = v.getClinSig();
        if (cs != null && !cs.isEmpty()) {
            String csLower = cs.toLowerCase();
            if (csLower.contains("pathogenic") && !csLower.contains("conflicting")) {
                if (csLower.contains("likely")) return "LIKELY_PATHOGENIC";
                return "PATHOGENIC";
            }
            if (csLower.contains("benign")) {
                if (csLower.contains("likely")) return "LIKELY_BENIGN";
                return "BENIGN";
            }
            if (csLower.contains("uncertain") || csLower.contains("conflicting")) {
                return "VUS";
            }
        }

        // 2. SIFT + PolyPhen (pentru variante necunoscute clinic)
        String sift = v.getSift();
        String poly = v.getPolyphen();
        if (sift != null && !sift.isEmpty() && poly != null && !poly.isEmpty()) {
            boolean siftBad = sift.contains("deleterious");
            boolean polyBad = poly.contains("damaging") || poly.contains("probably_damaging");
            if (siftBad && polyBad) return "LIKELY_PATHOGENIC";
            if (siftBad || polyBad) return "VUS";
            return "LIKELY_BENIGN";
        }

        // 3. Impact înalt (consecințe biologice cunoscute)
        String impact = v.getImpact();
        if ("HIGH".equals(impact)) return "LIKELY_PATHOGENIC";
        if ("MODERATE".equals(impact)) return "VUS";
        if ("LOW".equals(impact)) return "LIKELY_BENIGN";
        if ("MODIFIER".equals(impact)) return "LIKELY_BENIGN";

        return "UNKNOWN";
    }

    // Helper-uri
    private double parseDoubleSafe(String s) {
        try { return Double.parseDouble(s); }
        catch (NumberFormatException e) { return 0.0; }
    }

    private int parseIntSafe(String s) {
        try { return Integer.parseInt(s); }
        catch (NumberFormatException e) { return 0; }
    }
}