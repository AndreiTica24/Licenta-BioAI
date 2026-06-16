package ro.licenta.genomicsapi.controller;

import jakarta.validation.constraints.NotBlank;
import jakarta.validation.Valid;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import ro.licenta.genomicsapi.service.ChatbotService;

import java.util.Map;

/**
 * ChatbotController — endpoint pentru chatbot-ul de genetică.
 *
 * POST /api/chatbot/ask — trimite întrebare, primește răspuns
 *
 * Necesită JWT valid (USER sau ADMIN).
 */
@RestController
@RequestMapping("/api/chatbot")
public class ChatbotController {

    private final ChatbotService chatbotService;

    public ChatbotController(ChatbotService chatbotService) {
        this.chatbotService = chatbotService;
    }

    @PostMapping("/ask")
    public ResponseEntity<Map<String, Object>> ask(@Valid @RequestBody ChatRequest request) {
        String answer = chatbotService.ask(request.question);
        return ResponseEntity.ok(Map.of(
                "question", request.question,
                "answer", answer,
                "inDomain", chatbotService.isInDomain(request.question)
        ));
    }

    public static class ChatRequest {
        @NotBlank
        public String question;
    }
}