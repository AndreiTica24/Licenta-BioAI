package ro.licenta.genomicsapi.controller;

import org.springframework.beans.factory.annotation.Value;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;
import org.springframework.web.reactive.function.client.WebClient;

import java.util.HashMap;
import java.util.Map;

/**
 * HealthController — endpoint-uri pentru verificarea stării sistemului.
 *
 * GET /api/health         — status backend Java
 * GET /api/health/python  — verifică conexiunea cu Python AI service
 */
@RestController
@RequestMapping("/api/health")
public class HealthController {

    @Value("${app.python-api.url}")
    private String pythonApiUrl;

    /**
     * Test simplu — confirmă că backend-ul Java rulează.
     */
    @GetMapping
    public ResponseEntity<Map<String, Object>> health() {
        Map<String, Object> response = new HashMap<>();
        response.put("status", "healthy");
        response.put("service", "genomics-api");
        response.put("version", "1.0.0");
        response.put("java_version", System.getProperty("java.version"));
        return ResponseEntity.ok(response);
    }

    /**
     * Test conexiune Windows → WSL Python API.
     * Apelează GET /health pe serverul Python.
     */
    @GetMapping("/python")
    public ResponseEntity<Map<String, Object>> pythonHealth() {
        Map<String, Object> response = new HashMap<>();

        try {
            WebClient client = WebClient.create(pythonApiUrl);

            Map<String, Object> pythonResponse = client
                    .get()
                    .uri("/health")
                    .retrieve()
                    .bodyToMono(Map.class)
                    .block();

            response.put("java_status", "healthy");
            response.put("python_url", pythonApiUrl);
            response.put("python_response", pythonResponse);
            response.put("connection", "✅ OK");

            return ResponseEntity.ok(response);

        } catch (Exception e) {
            response.put("java_status", "healthy");
            response.put("python_url", pythonApiUrl);
            response.put("connection", "❌ FAILED");
            response.put("error", e.getMessage());
            response.put("hint", "Verifică dacă uvicorn rulează în WSL pe portul 8000");

            return ResponseEntity.status(503).body(response);
        }
    }
}