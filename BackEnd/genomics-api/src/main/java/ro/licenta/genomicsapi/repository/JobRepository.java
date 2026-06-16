package ro.licenta.genomicsapi.repository;

import org.springframework.data.jpa.repository.JpaRepository;
import org.springframework.stereotype.Repository;
import ro.licenta.genomicsapi.model.AnalysisJob;
import ro.licenta.genomicsapi.model.User;

import java.util.List;

@Repository
public interface JobRepository extends JpaRepository<AnalysisJob, String> {

    List<AnalysisJob> findByUserOrderByCreatedAtDesc(User user);

    List<AnalysisJob> findAllByOrderByCreatedAtDesc();

    long countByUser(User user);
}