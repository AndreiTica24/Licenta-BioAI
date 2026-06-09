package ro.licenta.genomicsapi.controller;

import org.springframework.http.ResponseEntity;
import org.springframework.security.access.prepost.PreAuthorize;
import org.springframework.web.bind.annotation.*;
import ro.licenta.genomicsapi.model.User;
import ro.licenta.genomicsapi.repository.UserRepository;

import java.util.List;
import java.util.Map;

@RestController
@RequestMapping("/api/admin")
@PreAuthorize("hasRole('ADMIN')")
public class AdminController {

    private final UserRepository userRepository;

    public AdminController(UserRepository userRepository) {
        this.userRepository = userRepository;
    }

    @GetMapping("/users")
    public ResponseEntity<List<Map<String, Object>>> getAllUsers() {
        List<Map<String, Object>> users = userRepository.findAll().stream()
                .map(u -> {
                    Map<String, Object> m = new java.util.HashMap<>();
                    m.put("id", u.getId());
                    m.put("email", u.getEmail());
                    m.put("fullName", u.getFullName());
                    m.put("role", u.getRole().name());
                    m.put("enabled", u.isEnabled());
                    m.put("createdAt", u.getCreatedAt());
                    return m;
                })
                .toList();
        return ResponseEntity.ok(users);
    }

    @DeleteMapping("/users/{id}")
    public ResponseEntity<Map<String, String>> deleteUser(@PathVariable Long id) {
        if (!userRepository.existsById(id)) {
            return ResponseEntity.notFound().build();
        }
        userRepository.deleteById(id);
        return ResponseEntity.ok(Map.of("message", "Utilizator șters cu succes"));
    }

    @GetMapping("/stats")
    public ResponseEntity<Map<String, Object>> getSystemStats() {
        Map<String, Object> stats = new java.util.HashMap<>();
        stats.put("total_users", userRepository.count());
        stats.put("total_admins", userRepository.findAll().stream()
                .filter(u -> u.getRole().name().equals("ADMIN")).count());
        stats.put("total_normal_users", userRepository.findAll().stream()
                .filter(u -> u.getRole().name().equals("USER")).count());
        return ResponseEntity.ok(stats);
    }
}