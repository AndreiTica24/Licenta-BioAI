package ro.licenta.genomicsapi.model;

import java.util.ArrayList;
import java.util.List;

/**
 * Variant — reprezintă o variantă genomică completă cu predicție AI
 * și anotare VEP/ClinVar.
 *
 * Câmpurile principale (din modelul AI):
 *   - chrom, pos, ref, alt — coordonatele variantei
 *   - predictedClass — Ref/Het/Hom-Alt din CNN
 *   - confidence — încrederea modelului (0-1)
 *
 * Câmpuri din anotarea VEP:
 *   - geneSymbol — numele genei (ex: BRCA1)
 *   - consequence — efect (ex: missense_variant)
 *   - impact — HIGH/MODERATE/LOW/MODIFIER
 *   - sift, polyphen — predicții impact pe proteină
 *   - hgvsc, hgvsp — notație standard variantă
 *
 * Câmpuri din ClinVar:
 *   - clinSig — clasificare clinică (Pathogenic/Benign/VUS)
 *   - clinDisease — boala asociată
 *   - clinReviewStatus — nivel de încredere ClinVar
 */
public class Variant {

    // ===== Coordonate (din predicția AI) =====
    private String chrom;
    private int pos;
    private String ref;
    private String alt;

    // ===== Predicție model CNN =====
    private String predictedClass;   // Ref / Het / Hom-Alt
    private double confidence;
    private int depth;
    private double af;

    // ===== Anotare VEP =====
    private String geneSymbol;       // ex: BRCA1
    private String geneId;            // ex: ENSG00000012048
    private String consequence;       // ex: missense_variant
    private String impact;            // HIGH | MODERATE | LOW | MODIFIER
    private String biotype;           // protein_coding, etc.
    private String hgvsc;             // ex: NM_007294.4:c.181T>G
    private String hgvsp;             // ex: NP_009225.1:p.Cys61Gly
    private String sift;              // ex: deleterious(0.01)
    private String polyphen;          // ex: probably_damaging(0.987)

    // ===== ClinVar =====
    private String clinSig;           // Pathogenic | Benign | VUS | ...
    private String clinDisease;       // Boala asociată
    private String clinReviewStatus;  // criteria_provided,multiple_submitters,...

    // ===== Clasificare finală (calculată) =====
    private String finalClassification; // PATHOGENIC | LIKELY_PATHOGENIC | VUS | LIKELY_BENIGN | BENIGN | UNKNOWN

    // ===== Constructor =====
    public Variant() {
    }

    // ===== Getters & Setters =====

    public String getChrom() { return chrom; }
    public void setChrom(String chrom) { this.chrom = chrom; }

    public int getPos() { return pos; }
    public void setPos(int pos) { this.pos = pos; }

    public String getRef() { return ref; }
    public void setRef(String ref) { this.ref = ref; }

    public String getAlt() { return alt; }
    public void setAlt(String alt) { this.alt = alt; }

    public String getPredictedClass() { return predictedClass; }
    public void setPredictedClass(String predictedClass) { this.predictedClass = predictedClass; }

    public double getConfidence() { return confidence; }
    public void setConfidence(double confidence) { this.confidence = confidence; }

    public int getDepth() { return depth; }
    public void setDepth(int depth) { this.depth = depth; }

    public double getAf() { return af; }
    public void setAf(double af) { this.af = af; }

    public String getGeneSymbol() { return geneSymbol; }
    public void setGeneSymbol(String geneSymbol) { this.geneSymbol = geneSymbol; }

    public String getGeneId() { return geneId; }
    public void setGeneId(String geneId) { this.geneId = geneId; }

    public String getConsequence() { return consequence; }
    public void setConsequence(String consequence) { this.consequence = consequence; }

    public String getImpact() { return impact; }
    public void setImpact(String impact) { this.impact = impact; }

    public String getBiotype() { return biotype; }
    public void setBiotype(String biotype) { this.biotype = biotype; }

    public String getHgvsc() { return hgvsc; }
    public void setHgvsc(String hgvsc) { this.hgvsc = hgvsc; }

    public String getHgvsp() { return hgvsp; }
    public void setHgvsp(String hgvsp) { this.hgvsp = hgvsp; }

    public String getSift() { return sift; }
    public void setSift(String sift) { this.sift = sift; }

    public String getPolyphen() { return polyphen; }
    public void setPolyphen(String polyphen) { this.polyphen = polyphen; }

    public String getClinSig() { return clinSig; }
    public void setClinSig(String clinSig) { this.clinSig = clinSig; }

    public String getClinDisease() { return clinDisease; }
    public void setClinDisease(String clinDisease) { this.clinDisease = clinDisease; }

    public String getClinReviewStatus() { return clinReviewStatus; }
    public void setClinReviewStatus(String clinReviewStatus) { this.clinReviewStatus = clinReviewStatus; }

    public String getFinalClassification() { return finalClassification; }
    public void setFinalClassification(String finalClassification) { this.finalClassification = finalClassification; }
}