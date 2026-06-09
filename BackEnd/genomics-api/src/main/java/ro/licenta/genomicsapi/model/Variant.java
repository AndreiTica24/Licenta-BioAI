package ro.licenta.genomicsapi.model;

import java.util.ArrayList;
import java.util.List;

public class Variant {

    private String chrom;
    private int pos;
    private String ref;
    private String alt;

    private String predictedClass;
    private double confidence;
    private int depth;
    private double af;

    private String geneSymbol;
    private String geneId;
    private String consequence;
    private String impact;
    private String biotype;
    private String hgvsc;
    private String hgvsp;
    private String sift;
    private String polyphen;

    private String clinSig;
    private String clinDisease;
    private String clinReviewStatus;

    private String finalClassification;

    public Variant() {
    }

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