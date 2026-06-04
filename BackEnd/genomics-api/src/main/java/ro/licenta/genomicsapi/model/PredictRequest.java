package ro.licenta.genomicsapi.model;

/**
 * PredictRequest — corespunde body-ului așteptat de Python API /predict.
 *
 * Python (FastAPI) așteaptă:
 *   {
 *     "bam_path": "...",
 *     "sample_name": "...",
 *     "threads": 4,
 *     "confidence": 0.7
 *   }
 */
public class PredictRequest {

    private String bamPath;
    private String sampleName;
    private int threads = 4;
    private double confidence = 0.7;

    public PredictRequest() {
    }

    public PredictRequest(String bamPath, String sampleName) {
        this.bamPath = bamPath;
        this.sampleName = sampleName;
    }

    // IMPORTANT: Python așteaptă "bam_path" (snake_case), nu "bamPath".
    // Numele JSON sunt mapate în PythonApiClient prin Map manual,
    // deci aici păstrăm camelCase Java standard.

    public String getBamPath() {
        return bamPath;
    }

    public void setBamPath(String bamPath) {
        this.bamPath = bamPath;
    }

    public String getSampleName() {
        return sampleName;
    }

    public void setSampleName(String sampleName) {
        this.sampleName = sampleName;
    }

    public int getThreads() {
        return threads;
    }

    public void setThreads(int threads) {
        this.threads = threads;
    }

    public double getConfidence() {
        return confidence;
    }

    public void setConfidence(double confidence) {
        this.confidence = confidence;
    }
}