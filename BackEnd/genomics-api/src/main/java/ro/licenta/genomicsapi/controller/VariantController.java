package ro.licenta.genomicsapi.controller;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.http.MediaType;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import org.springframework.web.multipart.MultipartFile;
import ro.licenta.genomicsapi.service.PythonApiClient;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.HashMap;
import java.util.Map;

/**
 * VariantController — endpoint-uri pentru variant calling.
 *
 * POST /api/variants/upload       — upload BAM + pornire analiză
 * GET  /api/variants/status/{id}  — status job
 * GET  /api/variants/result/{id}  — rezultat JSON
 * GET  /api/variants/vcf/{id}     — download VCF
 *
 * NOTĂ ARHITECTURĂ:
 * BAM-ul e salvat pe disk Windows (uploads/), iar Python (WSL) îl citește
 * prin calea /mnt/c/... — conversie făcută de toWslPath().
 */
@RestController
@RequestMapping("/api/variants")
public class VariantController {

    private static final Logger log = LoggerFactory.getLogger(VariantController.class);

    private final PythonApiClient pythonApiClient;

    @Value("${app.upload.directory}")
    private String uploadDir;

    public VariantController(PythonApiClient pythonApiClient) {
        this.pythonApiClient = pythonApiClient;
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
     * Convertește o cale Windows (C:\Users\...) în cale WSL (/mnt/c/Users/...).
     *
     * Exemplu:
     *   C:\Users\GOGUL\...\ploads\file.bam
     *   → /mnt/c/Users/GOGUL/.../uploads/file.bam
     */
    private String toWslPath(String windowsPath) {
        // Înlocuim backslash cu slash
        String path = windowsPath.replace("\\", "/");

        // Detectăm litera de drive (C:, D:, etc.) și o convertim
        if (path.length() >= 2 && path.charAt(1) == ':') {
            char driveLetter = Character.toLowerCase(path.charAt(0));
            path = "/mnt/" + driveLetter + path.substring(2);
        }

        return path;
    }
}