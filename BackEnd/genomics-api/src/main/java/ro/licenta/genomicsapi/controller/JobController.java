package ro.licenta.genomicsapi.controller;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.http.ResponseEntity;
import org.springframework.security.core.annotation.AuthenticationPrincipal;
import org.springframework.web.bind.annotation.*;
import ro.licenta.genomicsapi.model.AnalysisJob;
import ro.licenta.genomicsapi.model.User;
import ro.licenta.genomicsapi.repository.JobRepository;
import ro.licenta.genomicsapi.service.PythonApiClient;

import java.time.LocalDateTime;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

@RestController
@RequestMapping("/api/jobs")
public class JobController {

    private static final Logger log = LoggerFactory.getLogger(JobController.class);

    private final JobRepository jobRepository;
    private final PythonApiClient pythonApiClient;

    public JobController(JobRepository jobRepository, PythonApiClient pythonApiClient) {
        this.jobRepository = jobRepository;
        this.pythonApiClient = pythonApiClient;
    }

    @GetMapping("/my")
    public ResponseEntity<List<Map<String, Object>>> getMyJobs(
            @AuthenticationPrincipal User user) {
        List<AnalysisJob> jobs = jobRepository.findByUserOrderByCreatedAtDesc(user);

        // Sincronizăm joburile PROCESSING cu Python
        for (AnalysisJob job : jobs) {
            if (job.getStatus() == AnalysisJob.JobStatus.PROCESSING) {
                syncWithPython(job);
            }
        }

        return ResponseEntity.ok(jobs.stream().map(this::toDto).toList());
    }

    private void syncWithPython(AnalysisJob job) {
        try {
            Map<String, Object> status = pythonApiClient.getJobStatus(job.getId());
            String statusStr = (String) status.get("status");

            if ("completed".equalsIgnoreCase(statusStr)) {
                job.setStatus(AnalysisJob.JobStatus.COMPLETED);
                Object nVariants = status.get("n_variants");
                Object nCandidates = status.get("n_candidates");
                if (nVariants != null) job.setNVariants(((Number) nVariants).intValue());
                if (nCandidates != null) job.setNCandidates(((Number) nCandidates).intValue());
                job.setVcfPath("api_jobs/" + job.getId() + "/" + job.getSampleName() + ".vcf");
                job.setCompletedAt(LocalDateTime.now());
                jobRepository.save(job);
                log.info("Job {} sincronizat: COMPLETED ({} variante)",
                        job.getId(), job.getNVariants());
            } else if ("failed".equalsIgnoreCase(statusStr)) {
                job.setStatus(AnalysisJob.JobStatus.FAILED);
                Object err = status.get("error");
                if (err != null) job.setErrorMessage(err.toString());
                job.setCompletedAt(LocalDateTime.now());
                jobRepository.save(job);
                log.info("Job {} sincronizat: FAILED", job.getId());
            }
            // else: running/pending — rămâne PROCESSING
        } catch (Exception e) {
            // Python nu mai are job-ul (restart) — îl marcăm ca FAILED dacă e prea vechi
            // Pentru moment, doar log
            log.debug("Nu pot sincroniza job {} cu Python: {}", job.getId(), e.getMessage());
        }
    }

    @GetMapping("/all")
    @org.springframework.security.access.prepost.PreAuthorize("hasRole('ADMIN')")
    public ResponseEntity<List<Map<String, Object>>> getAllJobs() {
        List<AnalysisJob> jobs = jobRepository.findAllByOrderByCreatedAtDesc();
        return ResponseEntity.ok(jobs.stream().map(j -> {
            Map<String, Object> dto = toDto(j);
            dto.put("userEmail", j.getUser().getEmail());
            dto.put("userFullName", j.getUser().getFullName());
            return dto;
        }).toList());
    }

    @DeleteMapping("/{jobId}")
    public ResponseEntity<Map<String, String>> deleteJob(
            @PathVariable String jobId,
            @AuthenticationPrincipal User user) {

        AnalysisJob job = jobRepository.findById(jobId).orElse(null);
        if (job == null) {
            return ResponseEntity.notFound().build();
        }

        boolean isOwner = job.getUser().getId().equals(user.getId());
        boolean isAdmin = user.getRole().name().equals("ADMIN");

        if (!isOwner && !isAdmin) {
            return ResponseEntity.status(403).body(Map.of("error", "Forbidden"));
        }

        jobRepository.delete(job);
        return ResponseEntity.ok(Map.of("message", "Job șters cu succes"));
    }

    private Map<String, Object> toDto(AnalysisJob job) {
        Map<String, Object> dto = new HashMap<>();
        dto.put("id", job.getId());
        dto.put("bamFilename", job.getBamFilename());
        dto.put("sampleName", job.getSampleName());
        dto.put("status", job.getStatus().name());
        dto.put("nCandidates", job.getNCandidates());
        dto.put("nVariants", job.getNVariants());
        dto.put("hasVcf", job.getVcfPath() != null);
        dto.put("createdAt", job.getCreatedAt());
        dto.put("completedAt", job.getCompletedAt());
        dto.put("errorMessage", job.getErrorMessage());
        return dto;
    }
}