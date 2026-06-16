package ro.licenta.genomicsapi.model;

import jakarta.persistence.*;

import java.time.LocalDateTime;

/**
 * AnalysisJob — entitate JPA pentru tabela ANALYSIS_JOBS.
 *
 * Persistă informațiile despre fiecare analiză BAM, asociată unui utilizator.
 * VCF-ul de la AI rămâne salvat pe disk și poate fi re-anotat cu VEP oricând.
 */
@Entity
@Table(name = "analysis_jobs")
public class AnalysisJob {

    @Id
    @Column(length = 36)
    private String id;  // UUID, același ca job_id din Python

    @ManyToOne(fetch = FetchType.LAZY)
    @JoinColumn(name = "user_id", nullable = false)
    private User user;

    @Column(nullable = false)
    private String bamFilename;     // numele original (HG002_chr22.bam)

    @Column(nullable = false)
    private String sampleName;       // numele dat de user

    @Enumerated(EnumType.STRING)
    @Column(nullable = false)
    private JobStatus status = JobStatus.PROCESSING;

    private Integer nCandidates;     // câți candidați înainte de filtrare
    private Integer nVariants;       // câte variante a returnat AI

    @Column(length = 500)
    private String vcfPath;          // calea VCF generat de AI (api_jobs/{id}/sample.vcf)

    @Column(length = 1000)
    private String errorMessage;     // mesaj eroare dacă FAILED

    @Column(nullable = false)
    private LocalDateTime createdAt = LocalDateTime.now();

    private LocalDateTime completedAt;

    public AnalysisJob() {}

    public AnalysisJob(String id, User user, String bamFilename, String sampleName) {
        this.id = id;
        this.user = user;
        this.bamFilename = bamFilename;
        this.sampleName = sampleName;
    }

    public enum JobStatus {
        PROCESSING,
        COMPLETED,
        FAILED
    }

    // ===== Getters / Setters =====

    public String getId() { return id; }
    public void setId(String id) { this.id = id; }

    public User getUser() { return user; }
    public void setUser(User user) { this.user = user; }

    public String getBamFilename() { return bamFilename; }
    public void setBamFilename(String bamFilename) { this.bamFilename = bamFilename; }

    public String getSampleName() { return sampleName; }
    public void setSampleName(String sampleName) { this.sampleName = sampleName; }

    public JobStatus getStatus() { return status; }
    public void setStatus(JobStatus status) { this.status = status; }

    public Integer getNCandidates() { return nCandidates; }
    public void setNCandidates(Integer nCandidates) { this.nCandidates = nCandidates; }

    public Integer getNVariants() { return nVariants; }
    public void setNVariants(Integer nVariants) { this.nVariants = nVariants; }

    public String getVcfPath() { return vcfPath; }
    public void setVcfPath(String vcfPath) { this.vcfPath = vcfPath; }

    public String getErrorMessage() { return errorMessage; }
    public void setErrorMessage(String errorMessage) { this.errorMessage = errorMessage; }

    public LocalDateTime getCreatedAt() { return createdAt; }
    public void setCreatedAt(LocalDateTime createdAt) { this.createdAt = createdAt; }

    public LocalDateTime getCompletedAt() { return completedAt; }
    public void setCompletedAt(LocalDateTime completedAt) { this.completedAt = completedAt; }
}