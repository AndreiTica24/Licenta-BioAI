package ro.licenta.genomicsapi.config;


import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.security.config.annotation.web.builders.HttpSecurity;
import org.springframework.security.config.annotation.web.configuration.EnableWebSecurity;
import org.springframework.security.config.http.SessionCreationPolicy;
import org.springframework.security.web.SecurityFilterChain;

/**
 * SecurityConfig — configurare Spring Security.
 *
 * VERSIUNE TEMPORARĂ (pas 1): permite acces public la /api/health
 * pentru testarea de bază. JWT și endpoint-uri protejate vor fi
 * adăugate la pasul următor.
 */
@Configuration
@EnableWebSecurity
public class SecurityConfig {

    @Bean
    public SecurityFilterChain filterChain(HttpSecurity http) throws Exception {
        http
                // Dezactivăm CSRF pentru REST API (folosim JWT la pasul următor)
                .csrf(csrf -> csrf.disable())

                // Sesiune stateless (REST API)
                .sessionManagement(session -> session
                        .sessionCreationPolicy(SessionCreationPolicy.STATELESS))

                // Reguli de autorizare
                .authorizeHttpRequests(auth -> auth
                        // /api/health și /api/health/python sunt publice
                        .requestMatchers("/api/health/**").permitAll()
                        // Restul necesită auth (le adăugăm la pasul următor)
                        .anyRequest().permitAll()  // temporar permitem tot
                )

                // Dezactivăm HTTP Basic și form login
                .httpBasic(httpBasic -> httpBasic.disable())
                .formLogin(formLogin -> formLogin.disable());

        return http.build();
    }
}