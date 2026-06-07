package ro.licenta.genomicsapi.controller;

import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.GetMapping;

/**
 * PageController — servește paginile HTML cu Thymeleaf.
 *
 * Rute:
 *   /            → redirect spre /login (sau /dashboard dacă e logat)
 *   /login       → pagina login + register
 *   /dashboard   → upload BAM + istoric joburi
 *   /result/{id} → raport medical pentru un job
 *   /admin       → panou admin (doar pentru ADMIN)
 *
 * NOTĂ: aceste pagini sunt publice (HTML static).
 * Verificarea JWT se face în JavaScript și pe endpoint-urile /api/.
 * Dacă utilizatorul nu are token, este redirecționat la /login din JS.
 */
@Controller
public class PageController {

    @GetMapping("/")
    public String index() {
        return "redirect:/login";
    }

    @GetMapping("/login")
    public String login() {
        return "login";
    }

    @GetMapping("/dashboard")
    public String dashboard() {
        return "dashboard";
    }

    @GetMapping("/result/{jobId}")
    public String result() {
        return "result";
    }

    @GetMapping("/admin")
    public String admin() {
        return "admin";
    }
}