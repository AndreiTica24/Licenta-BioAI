"""
predict.py — Clasificare candidați cu CNN 1D
============================================================
Încarcă lista de candidați produsă de scan_bam.py, rulează modelul
CNN pe fiecare și generează două fișiere:

  1. predictions.vcf — VCF standard cu variantele detectate
                       (consumabil de VEP, ClinVar etc.)

  2. predictions.json — Format JSON pentru backend-ul Java
                        (lista variantelor cu prob și confidence)

ARGUMENT TEHNIC:
Modelul a fost antrenat pe poziții cu Ref/Het/Hom-Alt în proporție
echilibrată. Pe candidații din scan_bam (pre-filtrați cu AF>=0.15),
distribuția reală e diferită — așteptăm:
  - ~30-40% Ref (zgomot de secvențiere care a trecut filtrul)
  - ~30-35% Het
  - ~25-30% Hom-Alt

Doar predicțiile Het/Hom-Alt cu confidence >= 0.7 sunt scrise în VCF.

Rulare:
    python predict.py \\
        --candidates results/HG004_candidates.tsv \\
        --bam data/HG004_Mother/HG004.hiseq4000.wes-agilent.50x.dedup.grch38.bam \\
        --model checkpoints_unrelated/best_model.pth \\
        --output_vcf results/HG004_predicted.vcf \\
        --output_json results/HG004_predicted.json
"""

import argparse
import json
import logging
import time
from collections import Counter
from pathlib import Path
from typing import List, Dict, Tuple

import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

import pysam

from model import VariantCallerCNN1D
from dataset import _encode_window, WINDOW_SIZE, N_CHANNELS, REFERENCE_FASTA

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


CLASS_NAMES = ["Ref", "Het", "Hom-Alt"]
DEFAULT_CONFIDENCE_THRESHOLD = 0.7


# ============================================================================
# Dataset pentru candidați (citește din TSV)
# ============================================================================

class CandidateDataset(Dataset):
    """
    Dataset pentru candidații salvați de scan_bam.py în TSV.
    Returnează tensori (10, 200) gata de inferență CNN 1D.
    """

    def __init__(self, candidates_tsv: str, bam_path: str, fasta_path: str):
        self.bam_path   = bam_path
        self.fasta_path = fasta_path
        self.candidates: List[Dict] = []

        # Citim TSV-ul
        with open(candidates_tsv) as f:
            header = f.readline().strip().split("\t")
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 6:
                    continue
                self.candidates.append({
                    "chrom":     parts[0],
                    "pos":       int(parts[1]),
                    "ref_base":  parts[2],
                    "depth":     int(parts[3]),
                    "alt_count": int(parts[4]),
                    "AF":        float(parts[5]),
                })

        logger.info(f"Încărcați {len(self.candidates)} candidați din {candidates_tsv}")

    def __len__(self):
        return len(self.candidates)

    def __getitem__(self, idx):
        c = self.candidates[idx]

        bam   = pysam.AlignmentFile(self.bam_path, "rb")
        fasta = pysam.FastaFile(self.fasta_path)
        try:
            x = _encode_window(bam, c["chrom"], c["pos"],
                               window=WINDOW_SIZE, fasta=fasta)
        finally:
            bam.close()
            fasta.close()

        return torch.from_numpy(x), idx


# ============================================================================
# Determinarea bazei alternative din pileup
# ============================================================================

def get_alt_base(bam_path: str, chrom: str, pos: int, ref_base: str) -> str:
    """
    Returnează baza alternativă majoritară la o poziție.
    Folosit pentru a popula coloana ALT în VCF.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        for col in bam.pileup(chrom, pos, pos + 1, truncate=True,
                              min_base_quality=13, min_mapping_quality=20,
                              stepper="all"):
            if col.reference_pos != pos:
                continue
            counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            for r in col.pileups:
                if r.is_del or r.is_refskip:
                    continue
                try:
                    b = r.alignment.query_sequence[r.query_position].upper()
                    if b in counts and b != ref_base:
                        counts[b] += 1
                except (IndexError, TypeError):
                    continue
            if any(counts.values()):
                return max(counts, key=counts.get)
    finally:
        bam.close()
    return "N"


# ============================================================================
# Inferență batch
# ============================================================================

def run_inference(model, dataloader, device, n_total: int) -> List[Dict]:
    """
    Rulează modelul pe toți candidații și returnează lista de predicții.
    """
    model.eval()
    results = []

    t0 = time.time()
    n_done = 0

    with torch.no_grad():
        for x, indices in dataloader:
            x = x.to(device)
            logits = model(x)
            probs  = F.softmax(logits, dim=1).cpu().numpy()
            preds  = np.argmax(probs, axis=1)

            for i, idx in enumerate(indices):
                results.append({
                    "candidate_idx":   int(idx),
                    "predicted_label": int(preds[i]),
                    "predicted_class": CLASS_NAMES[preds[i]],
                    "prob_ref":        float(probs[i, 0]),
                    "prob_het":        float(probs[i, 1]),
                    "prob_hom_alt":    float(probs[i, 2]),
                    "confidence":      float(probs[i].max()),
                })

            n_done += len(indices)
            if n_done % 5000 < len(indices):
                elapsed = time.time() - t0
                eta = elapsed / n_done * (n_total - n_done)
                print(f"  [{n_done:>7,}/{n_total:>7,}] "
                      f"({elapsed:.0f}s, ETA {eta:.0f}s, "
                      f"{n_done/elapsed:.0f} candidați/s)", flush=True)

    # Sortăm rezultatele după candidate_idx pentru a corespunde ordinii inițiale
    results.sort(key=lambda r: r["candidate_idx"])
    return results


# ============================================================================
# Generare VCF
# ============================================================================

def write_vcf(predictions: List[Dict],
              candidates: List[Dict],
              bam_path: str,
              output_vcf: str,
              confidence_threshold: float = DEFAULT_CONFIDENCE_THRESHOLD) -> int:
    """
    Scrie un VCF standard cu variantele clasificate ca Het sau Hom-Alt.
    Returnează numărul de variante scrise.
    """
    Path(output_vcf).parent.mkdir(parents=True, exist_ok=True)

    n_written = 0

    with open(output_vcf, "w") as f:
        # Header VCF
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=VariantCallerCNN1D_v1.0\n")
        f.write("##reference=GRCh38\n")
        f.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency in BAM">\n')
        f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth at position">\n')
        f.write('##INFO=<ID=AO,Number=1,Type=Integer,Description="Alt observation count">\n')
        f.write('##INFO=<ID=CONF,Number=1,Type=Float,Description="Model confidence">\n')
        f.write('##INFO=<ID=GT_PRED,Number=1,Type=String,Description="Predicted genotype">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
        f.write('##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele frequency">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        for pred, cand in zip(predictions, candidates):
            # Skip Ref-uri și confidence sub prag
            if pred["predicted_label"] == 0:
                continue
            if pred["confidence"] < confidence_threshold:
                continue

            # Determinăm baza alternativă din BAM
            alt_base = get_alt_base(bam_path, cand["chrom"], cand["pos"], cand["ref_base"])
            if alt_base == "N":
                continue

            # GT pentru VCF: 0/1 = Het, 1/1 = Hom-Alt
            gt = "0/1" if pred["predicted_label"] == 1 else "1/1"

            # QUAL = -10 * log10(1 - confidence)
            qual = round(min(99.0, -10.0 * np.log10(max(1 - pred["confidence"], 1e-10))), 2)

            chrom_vcf = cand["chrom"]
            pos_vcf   = cand["pos"] + 1  # VCF e 1-indexed

            info = (f"AF={cand['AF']:.4f};"
                    f"DP={cand['depth']};"
                    f"AO={cand['alt_count']};"
                    f"CONF={pred['confidence']:.4f};"
                    f"GT_PRED={pred['predicted_class']}")

            sample = f"{gt}:{cand['depth']}:{cand['AF']:.4f}"

            f.write(f"{chrom_vcf}\t{pos_vcf}\t.\t{cand['ref_base']}\t{alt_base}\t"
                    f"{qual}\tPASS\t{info}\tGT:DP:AF\t{sample}\n")

            n_written += 1

    return n_written


# ============================================================================
# Generare JSON pentru backend Java
# ============================================================================

def write_json(predictions: List[Dict],
               candidates: List[Dict],
               bam_path: str,
               output_json: str,
               confidence_threshold: float = DEFAULT_CONFIDENCE_THRESHOLD) -> Dict:
    """
    Scrie un JSON cu toate variantele detectate.
    Format prietenos pentru backend-ul Java + VEP.
    """
    Path(output_json).parent.mkdir(parents=True, exist_ok=True)

    variants = []

    for pred, cand in zip(predictions, candidates):
        if pred["predicted_label"] == 0:
            continue
        if pred["confidence"] < confidence_threshold:
            continue

        alt_base = get_alt_base(bam_path, cand["chrom"], cand["pos"], cand["ref_base"])
        if alt_base == "N":
            continue

        variants.append({
            "chrom":           cand["chrom"],
            "pos":             cand["pos"] + 1,
            "ref":             cand["ref_base"],
            "alt":             alt_base,
            "predicted_class": pred["predicted_class"],
            "predicted_label": pred["predicted_label"],
            "prob_ref":        round(pred["prob_ref"], 4),
            "prob_het":        round(pred["prob_het"], 4),
            "prob_hom_alt":    round(pred["prob_hom_alt"], 4),
            "confidence":      round(pred["confidence"], 4),
            "depth":           cand["depth"],
            "alt_count":       cand["alt_count"],
            "AF":              cand["AF"],
        })

    n_het = sum(1 for v in variants if v["predicted_class"] == "Het")
    n_hom = sum(1 for v in variants if v["predicted_class"] == "Hom-Alt")

    output = {
        "metadata": {
            "n_variants":           len(variants),
            "n_het":                n_het,
            "n_hom_alt":            n_hom,
            "n_candidates_processed": len(predictions),
            "confidence_threshold": confidence_threshold,
            "model_version":        "VariantCallerCNN1D-v1.0",
        },
        "variants": variants,
    }

    with open(output_json, "w") as f:
        json.dump(output, f, indent=2)

    return output


# ============================================================================
# Pipeline principal
# ============================================================================

def predict(candidates_tsv: str,
            bam_path:       str,
            model_path:     str,
            output_vcf:     str,
            output_json:    str,
            fasta_path:     str = REFERENCE_FASTA,
            batch_size:     int = 256,
            num_workers:    int = 4,
            confidence:     float = DEFAULT_CONFIDENCE_THRESHOLD) -> Dict:
    """Pipeline complet de predicție."""

    # Device
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"\n🚀 Inferență pe: {device}", end="")
    if device.type == "cuda":
        print(f" ({torch.cuda.get_device_name(0)})")
    else:
        print(" (CPU)")

    # Model
    print(f"\n📥 Încărcăm modelul: {model_path}")
    model = VariantCallerCNN1D.load_from_checkpoint(model_path, device=str(device))

    # Dataset + DataLoader
    print(f"\n📦 Încărcăm candidații...")
    dataset = CandidateDataset(candidates_tsv, bam_path, fasta_path)
    loader  = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
        pin_memory=device.type == "cuda",
    )

    # Inferență
    print(f"\n⏳ Clasificare {len(dataset):,} candidați "
          f"(batch={batch_size}, workers={num_workers})...")
    t0 = time.time()
    predictions = run_inference(model, loader, device, len(dataset))
    inference_time = time.time() - t0

    # Distribuție predicții
    pred_dist = Counter(p["predicted_class"] for p in predictions)
    print()
    print("=" * 70)
    print("📊 DISTRIBUȚIE PREDICȚII (toate clasele, fără filtru confidence)")
    print("=" * 70)
    for cls in CLASS_NAMES:
        n = pred_dist[cls]
        pct = n / len(predictions) * 100 if predictions else 0
        bar = "█" * int(pct / 2)
        print(f"   {cls:10s} : {n:>7,} ({pct:5.1f}%)  {bar}")
    print()

    # Generare output VCF
    print(f"💾 Scriem VCF: {output_vcf}")
    n_vcf = write_vcf(predictions, dataset.candidates, bam_path,
                      output_vcf, confidence_threshold=confidence)

    # Generare output JSON
    print(f"💾 Scriem JSON: {output_json}")
    json_out = write_json(predictions, dataset.candidates, bam_path,
                          output_json, confidence_threshold=confidence)

    # Sumar final
    print()
    print("=" * 70)
    print(f"✅ REZULTATE FINALE (confidence >= {confidence})")
    print("=" * 70)
    print(f"   Candidați procesați   : {len(predictions):>7,}")
    print(f"   Variante VCF scrise   : {n_vcf:>7,}")
    print(f"   Het detectate         : {json_out['metadata']['n_het']:>7,}")
    print(f"   Hom-Alt detectate     : {json_out['metadata']['n_hom_alt']:>7,}")
    print(f"   Timp inferență CNN    : {inference_time:>7.1f}s "
          f"({len(predictions)/inference_time:.0f} cand/s)")
    print(f"   VCF                   : {output_vcf}")
    print(f"   JSON                  : {output_json}")
    print()
    print(f"💡 Pasul următor: anotare VEP + ClinVar pe {output_vcf}")

    return {
        "n_candidates":       len(predictions),
        "n_variants_vcf":     n_vcf,
        "n_het":              json_out["metadata"]["n_het"],
        "n_hom_alt":          json_out["metadata"]["n_hom_alt"],
        "inference_time_s":   round(inference_time, 2),
        "confidence_threshold": confidence,
    }


# ============================================================================
# CLI
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(description="Clasificare candidați cu CNN 1D")
    p.add_argument("--candidates", type=str, required=True,
                   help="TSV cu candidați (output scan_bam.py)")
    p.add_argument("--bam",        type=str, required=True,
                   help="Calea către BAM (pentru encoder + alt_base)")
    p.add_argument("--model",      type=str,
                   default="checkpoints_unrelated/best_model.pth",
                   help="Calea către modelul antrenat")
    p.add_argument("--fasta",      type=str, default=REFERENCE_FASTA)
    p.add_argument("--output_vcf",  type=str, required=True)
    p.add_argument("--output_json", type=str, required=True)
    p.add_argument("--batch_size",  type=int, default=256)
    p.add_argument("--num_workers", type=int, default=4)
    p.add_argument("--confidence",  type=float, default=DEFAULT_CONFIDENCE_THRESHOLD,
                   help="Prag minim confidence pentru a scrie în VCF")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    predict(
        candidates_tsv = args.candidates,
        bam_path       = args.bam,
        model_path     = args.model,
        output_vcf     = args.output_vcf,
        output_json    = args.output_json,
        fasta_path     = args.fasta,
        batch_size     = args.batch_size,
        num_workers    = args.num_workers,
        confidence     = args.confidence,
    )