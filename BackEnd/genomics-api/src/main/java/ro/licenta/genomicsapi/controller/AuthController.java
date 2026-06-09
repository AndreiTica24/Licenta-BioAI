package ro.licenta.genomicsapi.controller;

import jakarta.validation.Valid;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.http.ResponseEntity;
import org.springframework.security.authentication.AuthenticationManager;
import org.springframework.security.authentication.UsernamePasswordAuthenticationToken;
import org.springframework.security.crypto.password.PasswordEncoder;
import org.springframework.web.bind.annotation.*;
import ro.licenta.genomicsapi.model.AuthDtos.AuthResponse;
import ro.licenta.genomicsapi.model.AuthDtos.LoginRequest;
import ro.licenta.genomicsapi.model.AuthDtos.RegisterRequest;
import ro.licenta.genomicsapi.model.Role;
import ro.licenta.genomicsapi.model.User;
import ro.licenta.genomicsapi.repository.UserRepository;
import ro.licenta.genomicsapi.service.JwtService;

import java.util.HashMap;
import java.util.Map;

@RestController
@RequestMapping("/api/auth")
public class AuthController {

    private static final Logger log = LoggerFactory.getLogger(AuthController.class);

    private final UserRepository userRepository;
    private final PasswordEncoder passwordEncoder;
    private final JwtService jwtService;
    private final AuthenticationManager authenticationManager;

    @Value("${app.jwt.expiration-ms}")
    private long expirationMs;

    public AuthController(UserRepository userRepository,
                          PasswordEncoder passwordEncoder,
                          JwtService jwtService,
                          AuthenticationManager authenticationManager) {
        this.userRepository = userRepository;
        this.passwordEncoder = passwordEncoder;
        this.jwtService = jwtService;
        this.authenticationManager = authenticationManager;
    }

    @PostMapping("/register")
    public ResponseEntity<?> register(@Valid @RequestBody RegisterRequest request) {

        if (userRepository.existsByEmail(request.getEmail())) {
            Map<String, String> error = new HashMap<>();
            error.put("error", "Email-ul deja inregistrat");
            return ResponseEntity.badRequest().body(error);
        }

        User user = new User(
                request.getEmail(),
                passwordEncoder.encode(request.getPassword()),
                request.getFullName(),
                Role.USER
        );

        userRepository.save(user);
        log.info("Utilizator nou înregistrat: {}", user.getEmail());

        String token = jwtService.generateToken(user);

        return ResponseEntity.ok(new AuthResponse(
                token, user.getEmail(), user.getFullName(),
                user.getRole().name(), expirationMs
        ));
    }

    @PostMapping("/login")
    public ResponseEntity<?> login(@Valid @RequestBody LoginRequest request) {

        try {
            authenticationManager.authenticate(
                    new UsernamePasswordAuthenticationToken(
                            request.getEmail(), request.getPassword())
            );

            User user = userRepository.findByEmail(request.getEmail())
                    .orElseThrow();

            String token = jwtService.generateToken(user);
            log.info("Login reușit pentru: {}", user.getEmail());

            return ResponseEntity.ok(new AuthResponse(
                    token, user.getEmail(), user.getFullName(),
                    user.getRole().name(), expirationMs
            ));

        } catch (Exception e) {
            log.warn("Login esuat pentru: {}", request.getEmail());
            Map<String, String> error = new HashMap<>();
            error.put("error", "Email sau parola incorecta");
            return ResponseEntity.status(401).body(error);
        }
    }
}