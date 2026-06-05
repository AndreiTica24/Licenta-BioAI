package ro.licenta.genomicsapi.config;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.springframework.boot.CommandLineRunner;
import org.springframework.security.crypto.password.PasswordEncoder;
import org.springframework.stereotype.Component;
import ro.licenta.genomicsapi.model.Role;
import ro.licenta.genomicsapi.model.User;
import ro.licenta.genomicsapi.repository.UserRepository;

/**
 * DataInitializer — la primul start al aplicației, creează automat un admin
 * pentru demo. Email: admin@genomics.ro, parolă: Admin123!
 *
 * Util pentru ca juriul să poată testa imediat aplicația.
 */
@Component
public class DataInitializer implements CommandLineRunner {

    private static final Logger log = LoggerFactory.getLogger(DataInitializer.class);

    private static final String ADMIN_EMAIL    = "admin@genomics.ro";
    private static final String ADMIN_PASSWORD = "Admin123!";

    private static final String DEMO_USER_EMAIL    = "pacient@test.ro";
    private static final String DEMO_USER_PASSWORD = "Pacient123!";

    private final UserRepository userRepository;
    private final PasswordEncoder passwordEncoder;

    public DataInitializer(UserRepository userRepository,
                           PasswordEncoder passwordEncoder) {
        this.userRepository = userRepository;
        this.passwordEncoder = passwordEncoder;
    }

    @Override
    public void run(String... args) {
        // Admin implicit
        if (!userRepository.existsByEmail(ADMIN_EMAIL)) {
            User admin = new User(
                    ADMIN_EMAIL,
                    passwordEncoder.encode(ADMIN_PASSWORD),
                    "Administrator Sistem",
                    Role.ADMIN
            );
            userRepository.save(admin);
            log.info("");
            log.info("==============================================================");
            log.info("👑 Admin demo creat:");
            log.info("   Email   : {}", ADMIN_EMAIL);
            log.info("   Parolă  : {}", ADMIN_PASSWORD);
            log.info("==============================================================");
        }

        // User demo
        if (!userRepository.existsByEmail(DEMO_USER_EMAIL)) {
            User demo = new User(
                    DEMO_USER_EMAIL,
                    passwordEncoder.encode(DEMO_USER_PASSWORD),
                    "Pacient Demo",
                    Role.USER
            );
            userRepository.save(demo);
            log.info("👤 User demo creat:");
            log.info("   Email   : {}", DEMO_USER_EMAIL);
            log.info("   Parolă  : {}", DEMO_USER_PASSWORD);
            log.info("==============================================================");
            log.info("");
        }
    }
}