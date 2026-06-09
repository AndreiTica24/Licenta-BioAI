package ro.licenta.genomicsapi.service;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.http.MediaType;
import org.springframework.stereotype.Service;
import org.springframework.web.reactive.function.client.WebClient;
import reactor.util.retry.Retry;

import java.time.Duration;
import java.util.HashMap;
import java.util.Map;


@Service
public class PythonApiClient {

    private static final Logger log = LoggerFactory.getLogger(PythonApiClient.class);

    private final WebClient webClient;
    private final String apiKey;

    public PythonApiClient(
            @Value("${app.python-api.url}") String pythonApiUrl,
            @Value("${app.python-api.key}") String apiKey,
            @Value("${app.python-api.timeout-seconds}") int timeoutSeconds) {

        this.apiKey = apiKey;

        this.webClient = WebClient.builder()
                .baseUrl(pythonApiUrl)
                .codecs(configurer -> configurer
                        .defaultCodecs()
                        .maxInMemorySize(100 * 1024 * 1024)) // 100 MB
                .build();

        log.info("PythonApiClient inițializat: url={}, timeout={}s",
                pythonApiUrl, timeoutSeconds);
    }

    /**
     * Pornește un job de predicție pe Python.
     * @return job_id (UUID) returnat de Python
     */
    @SuppressWarnings("unchecked")
    public String startPrediction(String bamPath, String sampleName,
                                  int threads, double confidence) {
        // Construim body-ul JSON cu numele snake_case cerute de FastAPI
        Map<String, Object> body = new HashMap<>();
        body.put("bam_path", bamPath);
        body.put("sample_name", sampleName);
        body.put("threads", threads);
        body.put("confidence", confidence);

        log.info("Pornesc predicție Python: bam={}, sample={}", bamPath, sampleName);

        Map<String, Object> response = webClient.post()
                .uri("/predict")
                .header("X-API-Key", apiKey)
                .contentType(MediaType.APPLICATION_JSON)
                .bodyValue(body)
                .retrieve()
                .bodyToMono(Map.class)
                .block();

        if (response == null || !response.containsKey("job_id")) {
            throw new RuntimeException("Python API nu a returnat job_id");
        }

        String jobId = (String) response.get("job_id");
        log.info("Job Python creat: {}", jobId);
        return jobId;
    }

    @SuppressWarnings("unchecked")
    public Map<String, Object> getJobStatus(String jobId) {
        return webClient.get()
                .uri("/jobs/{jobId}", jobId)
                .header("X-API-Key", apiKey)
                .retrieve()
                .bodyToMono(Map.class)
                .block();
    }

    @SuppressWarnings("unchecked")
    public Map<String, Object> getJobResultJson(String jobId) {
        String rawJson = webClient.get()
                .uri("/jobs/{jobId}/result?format=json", jobId)
                .header("X-API-Key", apiKey)
                .retrieve()
                .bodyToMono(String.class)
                .block();

        if (rawJson == null) {
            throw new RuntimeException("Python API nu a returnat rezultatul JSON");
        }

        try {
            com.fasterxml.jackson.databind.ObjectMapper mapper =
                    new com.fasterxml.jackson.databind.ObjectMapper();
            return mapper.readValue(rawJson, Map.class);
        } catch (Exception e) {
            throw new RuntimeException("Eroare parsare JSON: " + e.getMessage(), e);
        }
    }

    public String getJobResultVcf(String jobId) {
        return webClient.get()
                .uri("/jobs/{jobId}/result?format=vcf", jobId)
                .header("X-API-Key", apiKey)
                .retrieve()
                .bodyToMono(String.class)
                .block();
    }

    @SuppressWarnings("unchecked")
    public Map<String, Object> health() {
        return webClient.get()
                .uri("/health")
                .retrieve()
                .bodyToMono(Map.class)
                .timeout(Duration.ofSeconds(5))
                .block();
    }
}