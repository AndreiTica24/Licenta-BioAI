"""
predict.py — Clasificare candidați cu CNN 1D (v2 OPTIMIZAT)
============================================================
Optimizări față de v1:
  1. BAM și FASTA persistente per worker (nu se deschid/închid per candidat)
  2. alt_base calculat inline în Dataset, nu în post-procesare
  3. Scriere VCF/JSON o singură dată la final
  4. Worker init function pentru pysam handles

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
import os
import time
from collections import Counter
from pathlib import Path
from typing import List, Dict

import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

import pysam

from model import VariantCallerCNN1D
from dataset import _encode_window, WINDOW_SIZE, REFERENCE_FASTA

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


CLASS_NAMES = ["Ref", "Het", "Hom-Alt"]
DEFAULT_CONFIDENCE_THRESHOLD = 0.7

# Variabile globale per-worker (deschise o singură dată)
_WORKER_BAM   = None
_WORKER_FASTA = None
_WORKER_PATHS = (None, None)


def _worker_init(bam_path: str, fasta_path: str):
    """Inițializează handles pysam per worker — apelat o singură dată per proces."""
    global _WORKER_BAM, _WORKER_FASTA, _WORKER_PATHS
    _WORKER_BAM   = pysam.AlignmentFile(bam_path, "rb")
    _WORKER_FASTA = pysam.FastaFile(fasta_path)
    _WORKER_PATHS = (bam_path, fasta_path)


# ============================================================================
# Dataset OPTIMIZAT
# ============================================================================

class CandidateDataset(Dataset):
    """
    Returnează (tensor_input, idx, alt_base) pentru fiecare candidat.
    BAM/FASTA persistente per worker via _worker_init().
    """
    def __init__(self, candidates_tsv: str, bam_path: str, fasta_path: str):
        self.bam_path   = bam_path
        self.fasta_path = fasta_path
        self.candidates: List[Dict] = []

        with open(candidates_tsv) as f:
            f.readline()  # header
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

    def _get_alt_base_fast(self, chrom: str, pos: int, ref_base: str) -> str:
        """Determină baza alt majoritară folosind BAM-ul persistent al workerului."""
        global _WORKER_BAM
        if _WORKER_BAM is None:
            return "N"
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        try:
            for col in _WORKER_BAM.pileup(chrom, pos, pos + 1, truncate=True,
                                           min_base_quality=13, min_mapping_quality=20,
                                           stepper="all"):
                if col.reference_pos != pos:
                    continue
                for r in col.pileups:
                    if r.is_del or r.is_refskip:
                        continue
                    try:
                        b = r.alignment.query_sequence[r.query_position].upper()
                        if b in counts and b != ref_base:
                            counts[b] += 1
                    except (IndexError, TypeError):
                        continue
                break
        except (ValueError, KeyError):
            pass
        if any(counts.values()):
            return max(counts, key=counts.get)
        return "N"

    def __getitem__(self, idx):
        global _WORKER_BAM, _WORKER_FASTA

        # Dacă suntem în main process (num_workers=0), deschidem direct
        if _WORKER_BAM is None:
            _worker_init(self.bam_path, self.fasta_path)

        c = self.candidates[idx]
        x = _encode_window(_WORKER_BAM, c["chrom"], c["pos"],
                           window=WINDOW_SIZE, fasta=_WORKER_FASTA)
        alt_base = self._get_alt_base_fast(c["chrom"], c["pos"], c["ref_base"])

        return torch.from_numpy(x), idx, alt_base


# ============================================================================
# Inferență
# ============================================================================

def run_inference(model, dataloader, device, n_total: int):
    """Rulează modelul și colectează predicții + alt_base."""
    model.eval()
    results = [None] * n_total

    t0 = time.time()
    n_done = 0
    LAST_PRINT = [0]

    with torch.no_grad():
        for batch in dataloader:
            x, indices, alt_bases = batch
            x = x.to(device, non_blocking=True)

            logits = model(x)
            probs  = F.softmax(logits, dim=1).cpu().numpy()
            preds  = np.argmax(probs, axis=1)

            for i, idx in enumerate(indices):
                idx_int = int(idx)
                results[idx_int] = {
                    "predicted_label": int(preds[i]),
                    "predicted_class": CLASS_NAMES[preds[i]],
                    "prob_ref":        float(probs[i, 0]),
                    "prob_het":        float(probs[i, 1]),
                    "prob_hom_alt":    float(probs[i, 2]),
                    "confidence":      float(probs[i].max()),
                    "alt_base":        alt_bases[i],
                }

            n_done += len(indices)

            # Progress la fiecare 10k candidați
            if n_done - LAST_PRINT[0] >= 10000 or n_done == n_total:
                elapsed = time.time() - t0
                eta = elapsed / n_done * (n_total - n_done) if n_done > 0 else 0
                rate = n_done / elapsed if elapsed > 0 else 0
                print(f"  [{n_done:>7,}/{n_total:>7,}] "
                      f"({elapsed:.0f}s, ETA {eta:.0f}s, "
                      f"{rate:.0f} candidați/s)", flush=True)
                LAST_PRINT[0] = n_done

    return results


# ============================================================================
# Output VCF + JSON (o singură trecere, fără BAM access)
# ============================================================================

def write_outputs(predictions: List[Dict],
                  candidates: List[Dict],
                  output_vcf: str,
                  output_json: str,
                  confidence_threshold: float) -> Dict:
    Path(output_vcf).parent.mkdir(parents=True, exist_ok=True)

    variants_json = []
    n_vcf = 0

    with open(output_vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=VariantCallerCNN1D_v1.0\n")
        f.write("##reference=GRCh38\n")
        f.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency">\n')
        f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">\n')
        f.write('##INFO=<ID=AO,Number=1,Type=Integer,Description="Alt count">\n')
        f.write('##INFO=<ID=CONF,Number=1,Type=Float,Description="Model confidence">\n')
        f.write('##INFO=<ID=GT_PRED,Number=1,Type=String,Description="Predicted genotype">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
        f.write('##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele frequency">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        for pred, cand in zip(predictions, candidates):
            if pred is None:
                continue
            if pred["predicted_label"] == 0:
                continue
            if pred["confidence"] < confidence_threshold:
                continue

            alt_base = pred["alt_base"]
            if alt_base == "N" or alt_base == cand["ref_base"]:
                continue

            gt   = "0/1" if pred["predicted_label"] == 1 else "1/1"
            qual = round(min(99.0, -10.0 * np.log10(
                max(1 - pred["confidence"], 1e-10))), 2)

            chrom_vcf = cand["chrom"]
            pos_vcf   = cand["pos"] + 1

            info = (f"AF={cand['AF']:.4f};DP={cand['depth']};AO={cand['alt_count']};"
                    f"CONF={pred['confidence']:.4f};GT_PRED={pred['predicted_class']}")
            sample = f"{gt}:{cand['depth']}:{cand['AF']:.4f}"

            f.write(f"{chrom_vcf}\t{pos_vcf}\t.\t{cand['ref_base']}\t{alt_base}\t"
                    f"{qual}\tPASS\t{info}\tGT:DP:AF\t{sample}\n")
            n_vcf += 1

            variants_json.append({
                "chrom":           cand["chrom"],
                "pos":             pos_vcf,
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

    n_het = sum(1 for v in variants_json if v["predicted_class"] == "Het")
    n_hom = sum(1 for v in variants_json if v["predicted_class"] == "Hom-Alt")

    json_output = {
        "metadata": {
            "n_variants":             len(variants_json),
            "n_het":                  n_het,
            "n_hom_alt":              n_hom,
            "n_candidates_processed": len([p for p in predictions if p is not None]),
            "confidence_threshold":   confidence_threshold,
            "model_version":          "VariantCallerCNN1D-v1.0",
        },
        "variants": variants_json,
    }

    Path(output_json).parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, "w") as f:
        json.dump(json_output, f, indent=2)

    return {"n_vcf": n_vcf, "n_het": n_het, "n_hom_alt": n_hom}


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

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"\n🚀 Inferență pe: {device}", end="")
    if device.type == "cuda":
        print(f" ({torch.cuda.get_device_name(0)})")
    else:
        print(" (CPU)")

    print(f"\n📥 Încărcăm modelul: {model_path}")
    model = VariantCallerCNN1D.load_from_checkpoint(model_path, device=str(device))

    print(f"\n📦 Încărcăm candidații...")
    dataset = CandidateDataset(candidates_tsv, bam_path, fasta_path)

    # Custom collate ca să transmitem alt_base ca listă de string-uri
    def collate_fn(batch):
        xs = torch.stack([b[0] for b in batch])
        indices = torch.tensor([b[1] for b in batch])
        alt_bases = [b[2] for b in batch]
        return xs, indices, alt_bases

    loader = DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
        pin_memory=device.type == "cuda",
        persistent_workers=num_workers > 0,
        collate_fn=collate_fn,
        worker_init_fn=lambda worker_id: _worker_init(bam_path, fasta_path),
    )

    print(f"\n⏳ Clasificare {len(dataset):,} candidați "
          f"(batch={batch_size}, workers={num_workers})...")
    t0 = time.time()
    predictions = run_inference(model, loader, device, len(dataset))
    inference_time = time.time() - t0

    pred_dist = Counter(p["predicted_class"] for p in predictions if p is not None)
    print()
    print("=" * 70)
    print("📊 DISTRIBUȚIE PREDICȚII (toate clasele)")
    print("=" * 70)
    total = sum(pred_dist.values())
    for cls in CLASS_NAMES:
        n = pred_dist[cls]
        pct = n / total * 100 if total > 0 else 0
        bar = "█" * int(pct / 2)
        print(f"   {cls:10s} : {n:>7,} ({pct:5.1f}%)  {bar}")
    print()

    print(f"💾 Scriem VCF + JSON (confidence >= {confidence})...")
    out = write_outputs(predictions, dataset.candidates,
                        output_vcf, output_json, confidence)

    print()
    print("=" * 70)
    print(f"✅ REZULTATE FINALE")
    print("=" * 70)
    print(f"   Candidați procesați   : {total:>7,}")
    print(f"   Variante VCF scrise   : {out['n_vcf']:>7,}")
    print(f"   Het detectate         : {out['n_het']:>7,}")
    print(f"   Hom-Alt detectate     : {out['n_hom_alt']:>7,}")
    print(f"   Timp inferență total  : {inference_time:>7.1f}s "
          f"({total/inference_time:.0f} cand/s)")
    print(f"   VCF                   : {output_vcf}")
    print(f"   JSON                  : {output_json}")
    print()
    print(f"💡 Pasul următor: anotare VEP + ClinVar pe {output_vcf}")

    return {
        "n_candidates":         total,
        "n_variants_vcf":       out["n_vcf"],
        "n_het":                out["n_het"],
        "n_hom_alt":            out["n_hom_alt"],
        "inference_time_s":     round(inference_time, 2),
        "confidence_threshold": confidence,
    }


def parse_args():
    p = argparse.ArgumentParser(description="Clasificare candidați cu CNN 1D")
    p.add_argument("--candidates",  type=str, required=True)
    p.add_argument("--bam",         type=str, required=True)
    p.add_argument("--model",       type=str,
                   default="checkpoints_unrelated/best_model.pth")
    p.add_argument("--fasta",       type=str, default=REFERENCE_FASTA)
    p.add_argument("--output_vcf",  type=str, required=True)
    p.add_argument("--output_json", type=str, required=True)
    p.add_argument("--batch_size",  type=int, default=256)
    p.add_argument("--num_workers", type=int, default=4)
    p.add_argument("--confidence",  type=float,
                   default=DEFAULT_CONFIDENCE_THRESHOLD)
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