package ro.licenta.genomicsapi.config;

import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.security.authentication.AuthenticationManager;
import org.springframework.security.authentication.AuthenticationProvider;
import org.springframework.security.authentication.dao.DaoAuthenticationProvider;
import org.springframework.security.config.annotation.authentication.configuration.AuthenticationConfiguration;
import org.springframework.security.config.annotation.method.configuration.EnableMethodSecurity;
import org.springframework.security.config.annotation.web.builders.HttpSecurity;
import org.springframework.security.config.annotation.web.configuration.EnableWebSecurity;
import org.springframework.security.config.http.SessionCreationPolicy;
import org.springframework.security.crypto.bcrypt.BCryptPasswordEncoder;
import org.springframework.security.crypto.password.PasswordEncoder;
import org.springframework.security.web.SecurityFilterChain;
import org.springframework.security.web.authentication.UsernamePasswordAuthenticationFilter;
import ro.licenta.genomicsapi.security.JwtAuthenticationFilter;
import ro.licenta.genomicsapi.service.CustomUserDetailsService;

/**
 * SecurityConfig — configurare completă Spring Security.
 *
 * Endpoint-uri publice (fără auth):
 *   /api/auth/**       — register, login
 *   /api/health/**     — verificare stare
 *   /h2-console/**     — DB console (doar pentru dezvoltare)
 *
 * Endpoint-uri protejate (necesită JWT valid):
 *   /api/variants/**   — upload BAM, status, rezultat
 *   /api/admin/**      — doar ADMIN (verificat cu @PreAuthorize)
 */
@Configuration
@EnableWebSecurity
@EnableMethodSecurity  // activează @PreAuthorize pe metode
public class SecurityConfig {

    private final JwtAuthenticationFilter jwtAuthFilter;
    private final CustomUserDetailsService userDetailsService;

    public SecurityConfig(JwtAuthenticationFilter jwtAuthFilter,
                          CustomUserDetailsService userDetailsService) {
        this.jwtAuthFilter = jwtAuthFilter;
        this.userDetailsService = userDetailsService;
    }

    @Bean
    public SecurityFilterChain filterChain(HttpSecurity http) throws Exception {
        http
                .csrf(csrf -> csrf.disable())
                .sessionManagement(s -> s
                        .sessionCreationPolicy(SessionCreationPolicy.STATELESS))

                .authorizeHttpRequests(auth -> auth
                        // Endpoint-uri publice
                        .requestMatchers("/api/auth/**").permitAll()
                        .requestMatchers("/api/health/**").permitAll()
                        .requestMatchers("/h2-console/**").permitAll()
                        .requestMatchers("/error").permitAll()

                        // Pagini HTML și resurse statice
                        .requestMatchers("/", "/login", "/dashboard", "/admin").permitAll()
                        .requestMatchers("/result/**").permitAll()
                        .requestMatchers("/css/**", "/js/**", "/img/**", "/favicon.ico").permitAll()

                        // Admin endpoint-uri
                        .requestMatchers("/api/admin/**").hasRole("ADMIN")

                        // Variante: USER sau ADMIN
                        .requestMatchers("/api/variants/**").hasAnyRole("USER", "ADMIN")

                        // Restul necesită autentificare
                        .anyRequest().authenticated()
                )

                // Pentru H2 console (frame embedding)
                .headers(h -> h.frameOptions(f -> f.disable()))

                .authenticationProvider(authenticationProvider())
                .addFilterBefore(jwtAuthFilter,
                        UsernamePasswordAuthenticationFilter.class);

        return http.build();
    }

    @Bean
    public AuthenticationProvider authenticationProvider() {
        DaoAuthenticationProvider provider = new DaoAuthenticationProvider();
        provider.setUserDetailsService(userDetailsService);
        provider.setPasswordEncoder(passwordEncoder());
        return provider;
    }

    @Bean
    public AuthenticationManager authenticationManager(
            AuthenticationConfiguration config) throws Exception {
        return config.getAuthenticationManager();
    }

    @Bean
    public PasswordEncoder passwordEncoder() {
        return new BCryptPasswordEncoder();
    }
}