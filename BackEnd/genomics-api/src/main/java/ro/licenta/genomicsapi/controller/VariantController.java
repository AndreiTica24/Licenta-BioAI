package ro.licenta.genomicsapi.controller;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.http.MediaType;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;
import ro.licenta.genomicsapi.model.Variant;
import ro.licenta.genomicsapi.service.PythonApiClient;
import ro.licenta.genomicsapi.service.VepAnnotationService;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * VariantController — endpoint-uri pentru variant calling + anotare.
 *
 * POST /api/variants/upload         — upload BAM + pornire analiză AI
 * GET  /api/variants/status/{id}    — status job AI
 * GET  /api/variants/result/{id}    — rezultat JSON din AI
 * GET  /api/variants/vcf/{id}       — download VCF de la AI
 * POST /api/variants/annotate/{id}  — anotare VEP+ClinVar pe VCF (după AI)
 */
@RestController
@RequestMapping("/api/variants")
public class VariantController {

    private static final Logger log = LoggerFactory.getLogger(VariantController.class);

    private final PythonApiClient pythonApiClient;
    private final VepAnnotationService vepService;

    @Value("${app.upload.directory}")
    private String uploadDir;

    @Value("${app.vep.data-dir:C:/vep_data}")
    private String vepDataDir;

    public VariantController(PythonApiClient pythonApiClient,
                             VepAnnotationService vepService) {
        this.pythonApiClient = pythonApiClient;
        this.vepService = vepService;
    }

    /**
     * Upload BAM + pornire analiză.
     * Streaming direct la disk (nu acumulează în RAM).
     */
    @PostMapping(value = "/upload", consumes = MediaType.MULTIPART_FORM_DATA_VALUE)
    public ResponseEntity<Map<String, Object>> uploadBam(
            @RequestParam("file") MultipartFile file,
            @RequestParam(value = "sampleName", required = false, defaultValue = "sample") String sampleName,
            @RequestParam(value = "confidence", required = false, defaultValue = "0.7") double confidence) {

        Map<String, Object> response = new HashMap<>();

        // Validări de bază
        if (file.isEmpty()) {
            response.put("error", "Fișierul e gol");
            return ResponseEntity.badRequest().body(response);
        }

        String originalName = file.getOriginalFilename();
        if (originalName == null || !originalName.toLowerCase().endsWith(".bam")) {
            response.put("error", "Fișierul trebuie să fie .bam");
            return ResponseEntity.badRequest().body(response);
        }

        try {
            // Creăm directorul de upload dacă nu există
            Path uploadPath = Paths.get(uploadDir).toAbsolutePath();
            Files.createDirectories(uploadPath);

            // Numele fișierului salvat (evităm coliziuni cu timestamp)
            String savedName = System.currentTimeMillis() + "_" + originalName;
            Path bamPath = uploadPath.resolve(savedName);

            // STREAMING upload — copiem stream-ul direct la disk
            log.info("Upload BAM: {} ({} MB)", originalName,
                    file.getSize() / (1024 * 1024));

            long t0 = System.currentTimeMillis();
            try (InputStream in = file.getInputStream();
                 OutputStream out = Files.newOutputStream(bamPath)) {
                byte[] buffer = new byte[8 * 1024 * 1024]; // 8 MB chunks
                int bytesRead;
                while ((bytesRead = in.read(buffer)) != -1) {
                    out.write(buffer, 0, bytesRead);
                }
            }
            long uploadTime = (System.currentTimeMillis() - t0) / 1000;
            log.info("BAM salvat în {} ({}s)", bamPath, uploadTime);

            // IMPORTANT: BAM-ul are nevoie de index .bai. Verificăm dacă există.
            // Pentru demo, presupunem că .bai vine separat sau e generat de Python.
            Path baiPath = Paths.get(bamPath.toString() + ".bai");
            boolean hasBai = Files.exists(baiPath);

            // Convertim calea Windows în cale WSL (/mnt/c/...)
            String wslBamPath = toWslPath(bamPath.toString());
            log.info("Cale WSL pentru Python: {}", wslBamPath);

            // Pornim predicția pe Python
            String jobId = pythonApiClient.startPrediction(
                    wslBamPath, sampleName, 4, confidence);

            response.put("status", "uploaded");
            response.put("job_id", jobId);
            response.put("filename", originalName);
            response.put("size_mb", file.getSize() / (1024 * 1024));
            response.put("upload_time_s", uploadTime);
            response.put("has_index", hasBai);
            response.put("message", "BAM încărcat. Analiza a pornit.");

            if (!hasBai) {
                response.put("warning",
                        "Index .bai lipsește — Python îl va genera (poate dura mai mult)");
            }

            return ResponseEntity.ok(response);

        } catch (IOException e) {
            log.error("Eroare upload BAM", e);
            response.put("error", "Eroare la salvarea fișierului: " + e.getMessage());
            return ResponseEntity.internalServerError().body(response);
        } catch (Exception e) {
            log.error("Eroare pornire analiză", e);
            response.put("error", "Eroare la pornirea analizei: " + e.getMessage());
            return ResponseEntity.internalServerError().body(response);
        }
    }

    /**
     * Status job.
     */
    @GetMapping("/status/{jobId}")
    public ResponseEntity<Map<String, Object>> getStatus(@PathVariable String jobId) {
        try {
            Map<String, Object> status = pythonApiClient.getJobStatus(jobId);
            return ResponseEntity.ok(status);
        } catch (Exception e) {
            Map<String, Object> error = new HashMap<>();
            error.put("error", e.getMessage());
            return ResponseEntity.status(404).body(error);
        }
    }

    /**
     * Rezultat JSON (variantele detectate).
     */
    @GetMapping("/result/{jobId}")
    public ResponseEntity<Map<String, Object>> getResult(@PathVariable String jobId) {
        try {
            Map<String, Object> result = pythonApiClient.getJobResultJson(jobId);
            return ResponseEntity.ok(result);
        } catch (Exception e) {
            Map<String, Object> error = new HashMap<>();
            error.put("error", e.getMessage());
            return ResponseEntity.status(400).body(error);
        }
    }

    /**
     * Download VCF.
     */
    @GetMapping(value = "/vcf/{jobId}", produces = MediaType.TEXT_PLAIN_VALUE)
    public ResponseEntity<String> getVcf(@PathVariable String jobId) {
        try {
            String vcf = pythonApiClient.getJobResultVcf(jobId);
            return ResponseEntity.ok()
                    .header("Content-Disposition",
                            "attachment; filename=\"variants_" + jobId + ".vcf\"")
                    .body(vcf);
        } catch (Exception e) {
            return ResponseEntity.status(400).body("Eroare: " + e.getMessage());
        }
    }

    /**
     * Anotează un VCF (de la un job AI completat) cu VEP+ClinVar.
     * Descarcă VCF de la Python, apelează VEP Docker, returnează variante anotate.
     */
    @PostMapping("/annotate/{jobId}")
    public ResponseEntity<Map<String, Object>> annotateVariants(
            @PathVariable String jobId) {

        Map<String, Object> response = new HashMap<>();

        try {
            // 1. Descărcăm VCF-ul de la Python
            log.info("[{}] Descărcăm VCF de la Python...", jobId);
            String vcfContent = pythonApiClient.getJobResultVcf(jobId);

            // 2. Salvăm temporar în vep_data/input/
            Path vepInputDir = Paths.get(vepDataDir, "input");
            Files.createDirectories(vepInputDir);
            Path vcfPath = vepInputDir.resolve("ai_predictions_" + jobId + ".vcf");
            Files.writeString(vcfPath, vcfContent);
            log.info("[{}] VCF salvat: {} ({} variante)", jobId, vcfPath,
                    vcfContent.lines().filter(l -> !l.startsWith("#")).count());

            // 3. Rulăm VEP
            log.info("[{}] Pornesc anotare VEP...", jobId);
            long t0 = System.currentTimeMillis();
            List<Variant> variants = vepService.annotateVcf(vcfPath.toString());
            long vepTime = (System.currentTimeMillis() - t0) / 1000;

            // 4. Statistici
            Map<String, Long> byClassification = variants.stream()
                    .collect(Collectors.groupingBy(
                            v -> v.getFinalClassification() != null
                                    ? v.getFinalClassification() : "UNKNOWN",
                            Collectors.counting()
                    ));

            long withClinvar = variants.stream()
                    .filter(v -> v.getClinSig() != null && !v.getClinSig().isEmpty())
                    .count();

            long withGene = variants.stream()
                    .filter(v -> v.getGeneSymbol() != null && !v.getGeneSymbol().isEmpty())
                    .count();

            // 5. Răspuns
            response.put("job_id", jobId);
            response.put("n_variants", variants.size());
            response.put("with_clinvar", withClinvar);
            response.put("with_gene", withGene);
            response.put("by_classification", byClassification);
            response.put("annotation_time_s", vepTime);
            response.put("variants", variants);

            log.info("[{}] Anotare completă: {} variante, {} cu ClinVar, {} cu genă, {}s",
                    jobId, variants.size(), withClinvar, withGene, vepTime);

            return ResponseEntity.ok(response);

        } catch (Exception e) {
            log.error("Eroare anotare VEP pentru job " + jobId, e);
            response.put("error", e.getMessage());
            return ResponseEntity.internalServerError().body(response);
        }
    }

    /**
     * Convertește o cale Windows (C:\Users\...) în cale WSL (/mnt/c/Users/...).
     */
    private String toWslPath(String windowsPath) {
        String path = windowsPath.replace("\\", "/");
        if (path.length() >= 2 && path.charAt(1) == ':') {
            char driveLetter = Character.toLowerCase(path.charAt(0));
            path = "/mnt/" + driveLetter + path.substring(2);
        }
        return path;
    }
}