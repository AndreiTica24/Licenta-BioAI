package ro.licenta.genomicsapi.model;

import jakarta.validation.constraints.Email;
import jakarta.validation.constraints.NotBlank;
import jakarta.validation.constraints.Size;

/**
 * DTO-uri pentru autentificare.
 * Toate într-un singur fișier pentru simplitate.
 */
public class AuthDtos {

    public static class RegisterRequest {
        @Email(message = "Email invalid")
        @NotBlank(message = "Email-ul este obligatoriu")
        private String email;

        @NotBlank(message = "Parola este obligatorie")
        @Size(min = 8, message = "Parola trebuie să aibă minim 8 caractere")
        private String password;

        @NotBlank(message = "Numele este obligatoriu")
        private String fullName;

        public String getEmail() { return email; }
        public void setEmail(String email) { this.email = email; }
        public String getPassword() { return password; }
        public void setPassword(String password) { this.password = password; }
        public String getFullName() { return fullName; }
        public void setFullName(String fullName) { this.fullName = fullName; }
    }

    public static class LoginRequest {
        @Email
        @NotBlank
        private String email;

        @NotBlank
        private String password;

        public String getEmail() { return email; }
        public void setEmail(String email) { this.email = email; }
        public String getPassword() { return password; }
        public void setPassword(String password) { this.password = password; }
    }

    public static class AuthResponse {
        private String token;
        private String email;
        private String fullName;
        private String role;
        private long expiresInMs;

        public AuthResponse() {}

        public AuthResponse(String token, String email, String fullName,
                            String role, long expiresInMs) {
            this.token = token;
            this.email = email;
            this.fullName = fullName;
            this.role = role;
            this.expiresInMs = expiresInMs;
        }

        public String getToken() { return token; }
        public void setToken(String token) { this.token = token; }
        public String getEmail() { return email; }
        public void setEmail(String email) { this.email = email; }
        public String getFullName() { return fullName; }
        public void setFullName(String fullName) { this.fullName = fullName; }
        public String getRole() { return role; }
        public void setRole(String role) { this.role = role; }
        public long getExpiresInMs() { return expiresInMs; }
        public void setExpiresInMs(long expiresInMs) { this.expiresInMs = expiresInMs; }
    }
}