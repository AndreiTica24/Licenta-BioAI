package ro.licenta.genomicsapi.controller;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.http.MediaType;
import org.springframework.http.ResponseEntity;
import org.springframework.security.core.annotation.AuthenticationPrincipal;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;
import ro.licenta.genomicsapi.model.AnalysisJob;
import ro.licenta.genomicsapi.model.User;
import ro.licenta.genomicsapi.model.Variant;
import ro.licenta.genomicsapi.repository.JobRepository;
import ro.licenta.genomicsapi.service.PythonApiClient;
import ro.licenta.genomicsapi.service.VepAnnotationService;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.time.LocalDateTime;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * VariantController — endpoint-uri pentru variant calling + anotare.
 *
 * Asociază fiecare job cu un user în DB pentru istoric persistent.
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
    private final JobRepository jobRepository;

    @Value("${app.upload.directory}")
    private String uploadDir;

    @Value("${app.vep.data-dir:C:/vep_data}")
    private String vepDataDir;

    public VariantController(PythonApiClient pythonApiClient,
                             VepAnnotationService vepService,
                             JobRepository jobRepository) {
        this.pythonApiClient = pythonApiClient;
        this.vepService = vepService;
        this.jobRepository = jobRepository;
    }

    /**
     * Upload BAM + pornire analiză.
     * Streaming direct la disk (nu acumulează în RAM).
     * Creează AnalysisJob în DB pentru istoricul user-ului.
     */
    @PostMapping(value = "/upload", consumes = MediaType.MULTIPART_FORM_DATA_VALUE)
    public ResponseEntity<Map<String, Object>> uploadBam(
            @RequestParam("file") MultipartFile file,
            @RequestParam(value = "sampleName", required = false, defaultValue = "sample") String sampleName,
            @RequestParam(value = "confidence", required = false, defaultValue = "0.7") double confidence,
            @AuthenticationPrincipal User currentUser) {

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
            log.info("Upload BAM: {} ({} MB) by user {}", originalName,
                    file.getSize() / (1024 * 1024), currentUser.getEmail());

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

            Path baiPath = Paths.get(bamPath.toString() + ".bai");
            boolean hasBai = Files.exists(baiPath);

            // Convertim calea Windows în cale WSL (/mnt/c/...)
            String wslBamPath = toWslPath(bamPath.toString());
            log.info("Cale WSL pentru Python: {}", wslBamPath);

            // Pornim predicția pe Python
            String jobId = pythonApiClient.startPrediction(
                    wslBamPath, sampleName, 4, confidence);

            // Salvăm AnalysisJob în H2 DB pentru istoricul user-ului
            AnalysisJob job = new AnalysisJob(jobId, currentUser, originalName, sampleName);
            jobRepository.save(job);
            log.info("AnalysisJob salvat în DB: {}", jobId);

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
     * Când job-ul e completat, actualizează AnalysisJob în DB (status, n_variants, vcf_path).
     */
    @GetMapping("/status/{jobId}")
    public ResponseEntity<Map<String, Object>> getStatus(@PathVariable String jobId) {
        try {
            Map<String, Object> status = pythonApiClient.getJobStatus(jobId);

            // Sincronizăm cu DB dacă job-ul e completed/failed
            String statusStr = (String) status.get("status");
            if ("completed".equalsIgnoreCase(statusStr) || "failed".equalsIgnoreCase(statusStr)) {
                jobRepository.findById(jobId).ifPresent(job -> {
                    if (job.getStatus() == AnalysisJob.JobStatus.PROCESSING) {
                        // Update prima dată când vedem că s-a terminat
                        if ("completed".equalsIgnoreCase(statusStr)) {
                            job.setStatus(AnalysisJob.JobStatus.COMPLETED);
                            Object nVariants = status.get("n_variants");
                            Object nCandidates = status.get("n_candidates");
                            if (nVariants != null) job.setNVariants(((Number) nVariants).intValue());
                            if (nCandidates != null) job.setNCandidates(((Number) nCandidates).intValue());
                            // VCF path-ul din Python (din api_jobs/{jobId}/sample.vcf)
                            job.setVcfPath("api_jobs/" + jobId + "/" + job.getSampleName() + ".vcf");
                        } else {
                            job.setStatus(AnalysisJob.JobStatus.FAILED);
                            Object err = status.get("error");
                            if (err != null) job.setErrorMessage(err.toString());
                        }
                        job.setCompletedAt(LocalDateTime.now());
                        jobRepository.save(job);
                        log.info("Job {} actualizat în DB: {}", jobId, job.getStatus());
                    }
                });
            }

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
     * Descarcă VCF de la Python (sau de pe disk dacă există în DB), apelează VEP Docker.
     */
    @PostMapping("/annotate/{jobId}")
    public ResponseEntity<Map<String, Object>> annotateVariants(
            @PathVariable String jobId) {

        Map<String, Object> response = new HashMap<>();

        try {
            // 1. Obținem conținut VCF — încercăm întâi Python (job activ),
            //    apoi disk (job vechi din DB, după restart Python)
            String vcfContent = null;
            try {
                log.info("[{}] Încerc descărcare VCF de la Python...", jobId);
                vcfContent = pythonApiClient.getJobResultVcf(jobId);
                log.info("[{}] VCF obținut de la Python", jobId);
            } catch (Exception pythonErr) {
                log.warn("[{}] Python nu are job-ul, încerc de pe disk", jobId);
                // Caut VCF-ul pe disk prin DB
                AnalysisJob job = jobRepository.findById(jobId).orElse(null);
                if (job == null) {
                    throw new RuntimeException("Job not found in database");
                }
                if (job.getVcfPath() == null) {
                    throw new RuntimeException("Job has no VCF path");
                }
                // VCF e relativ la folderul AI-Module
                Path vcfDiskPath = Paths.get(
                        System.getProperty("user.home"),
                        "Licenta-BioAI", "AI-Module", job.getVcfPath()
                );
                // Fallback: încearcă pe WSL prin /mnt
                if (!Files.exists(vcfDiskPath)) {
                    vcfDiskPath = Paths.get(
                            "\\\\wsl$\\Ubuntu\\home\\andrei\\Licenta-BioAI\\AI-Module",
                            job.getVcfPath().replace("/", "\\")
                    );
                }
                if (!Files.exists(vcfDiskPath)) {
                    throw new RuntimeException("VCF not found on disk: " + job.getVcfPath());
                }
                vcfContent = Files.readString(vcfDiskPath);
                log.info("[{}] VCF citit de pe disk: {}", jobId, vcfDiskPath);
            }

            // 2. Salvăm în vep_data/input/
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