"""
evaluate.py — Evaluare finală pe HG004-Mother (test set nevăzut)
Generează raport complet + export JSON pentru backend Java.

Rulare:
    python evaluate.py
    python evaluate.py --model checkpoints/best_model.pth --max_samples 5000
"""

import argparse
import json
import logging
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader
from sklearn.metrics import (
    f1_score, precision_score, recall_score,
    classification_report, confusion_matrix,
    roc_auc_score,
)

from model import VariantCallerCNN
from dataset import GenomicDataset

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

CLASS_NAMES = ["Ref", "Het", "Hom-Alt"]

# ---------------------------------------------------------------------------
# Config HG004 (test set)
# ---------------------------------------------------------------------------
DEFAULT_EVAL_CONFIG = {
    "bam_test": "data/HG004_Mother/HG004.hiseq4000.wes-agilent.50x.dedup.grch38.bam",
    "vcf_test": "data/HG004_Mother/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf",
    "bed_test": "data/HG004_Mother/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
    "max_samples": 10_000,
    "batch_size":  64,
    "num_workers": 4,
    "seed":        42,
    "model_path":  "checkpoints/best_model.pth",
    "output_dir":  "results",
}


# ---------------------------------------------------------------------------
# Funcție evaluare cu probabilități
# ---------------------------------------------------------------------------

def evaluate_with_probs(model, loader, device):
    """
    Returnează predicții, labels și probabilități softmax.
    """
    model.eval()
    all_preds  = []
    all_labels = []
    all_probs  = []

    with torch.no_grad():
        for images, labels in loader:
            images = images.to(device)
            logits = model(images)
            probs  = F.softmax(logits, dim=1)
            _, predicted = torch.max(logits, 1)

            all_preds.extend(predicted.cpu().numpy().tolist())
            all_labels.extend(labels.numpy().tolist())
            all_probs.extend(probs.cpu().numpy().tolist())

    return all_labels, all_preds, all_probs


# ---------------------------------------------------------------------------
# Export JSON pentru backend Java
# ---------------------------------------------------------------------------

def export_for_java(dataset, all_labels, all_preds, all_probs, output_path: str):
    """
    Exportă predicțiile în format JSON consumabil de backend-ul Java.
    Structura: Listă de variante cu predicție + probabilități.
    """
    variants = []
    for i, dp in enumerate(dataset.data_points):
        if i >= len(all_preds):
            break

        pred_label = all_preds[i]
        true_label = all_labels[i]
        probs      = all_probs[i]

        # Includem DOAR variantele non-Ref (pred sau adevărat)
        # Backend-ul Java va filtra prin ClinVar
        entry = {
            "chrom":           dp["chrom"],
            "pos":             dp["pos"] + 1,   # 1-indexed (VCF standard)
            "ref":             dp["ref"],
            "alt":             dp["alt"],
            "predicted_class": CLASS_NAMES[pred_label],
            "predicted_label": pred_label,
            "true_class":      CLASS_NAMES[true_label],
            "true_label":      true_label,
            "prob_ref":        round(probs[0], 4),
            "prob_het":        round(probs[1], 4),
            "prob_hom_alt":    round(probs[2], 4),
            "confidence":      round(max(probs), 4),
            "is_variant":      pred_label > 0,   # True dacă e Het sau Hom-Alt
        }
        variants.append(entry)

    # Statistici sumare
    n_variants = sum(1 for v in variants if v["is_variant"])
    n_het      = sum(1 for v in variants if v["predicted_class"] == "Het")
    n_hom      = sum(1 for v in variants if v["predicted_class"] == "Hom-Alt")

    output = {
        "metadata": {
            "total_positions": len(variants),
            "n_variants":      n_variants,
            "n_het":           n_het,
            "n_hom_alt":       n_hom,
            "model_version":   "VariantCallerCNN-v1",
        },
        "variants": variants,
    }

    with open(output_path, "w") as f:
        json.dump(output, f, indent=2)

    logger.info(f"Export JSON: {output_path} ({len(variants)} poziții, {n_variants} variante)")
    return output


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_evaluation(config: dict):
    Path(config["output_dir"]).mkdir(parents=True, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"\n🔬 Evaluare pe: {device}")
    if device.type == "cuda":
        print(f"   GPU: {torch.cuda.get_device_name(0)}")

    # --- Încarcă modelul ---
    print(f"\n📥 Încărcăm modelul: {config['model_path']}")
    model = VariantCallerCNN.load_from_checkpoint(config["model_path"], device=str(device))
    model.eval()

    # --- Dataset test (HG004) ---
    print(f"\n📦 Încărcăm datele de test (HG004-Mother)...")
    test_dataset = GenomicDataset(
        bam_path    = config["bam_test"],
        vcf_path    = config["vcf_test"],
        bed_path    = config.get("bed_test"),
        max_samples = config["max_samples"],
        seed        = config["seed"],
    )

    test_loader = DataLoader(
        test_dataset,
        batch_size  = config["batch_size"],
        shuffle     = False,
        num_workers = config["num_workers"],
        pin_memory  = device.type == "cuda",
    )

    # --- Evaluare ---
    print(f"\n⏳ Evaluăm {len(test_dataset)} exemple...")
    all_labels, all_preds, all_probs = evaluate_with_probs(model, test_loader, device)

    # --- Metrici ---
    f1   = f1_score(all_labels, all_preds, average="macro", zero_division=0)
    prec = precision_score(all_labels, all_preds, average="macro", zero_division=0)
    rec  = recall_score(all_labels, all_preds, average="macro", zero_division=0)
    report = classification_report(
        all_labels, all_preds, target_names=CLASS_NAMES, zero_division=0
    )
    cm = confusion_matrix(all_labels, all_preds, labels=[0, 1, 2])

    # AUC-ROC (One-vs-Rest)
    try:
        probs_np = np.array(all_probs)
        auc = roc_auc_score(
            all_labels, probs_np, multi_class="ovr", average="macro"
        )
    except Exception:
        auc = None

    # --- Afișare ---
    print("\n" + "=" * 70)
    print("📊 REZULTATE FINALE — HG004-Mother (Test Set Nevăzut)")
    print("=" * 70)
    print(f"   F1  macro  : {f1:.4f}")
    print(f"   Precision  : {prec:.4f}")
    print(f"   Recall     : {rec:.4f}")
    if auc is not None:
        print(f"   AUC-ROC    : {auc:.4f}")
    print(f"\n{report}")
    print("Confusion Matrix (rând=adevărat, coloană=prezis):")
    print(f"   {'':10s} " + "  ".join(f"{n:8s}" for n in CLASS_NAMES))
    for i, row in enumerate(cm):
        vals = "  ".join(f"{v:8d}" for v in row)
        print(f"   {CLASS_NAMES[i]:10s} {vals}")
    print("=" * 70)

    # --- Salvare raport text ---
    report_path = Path(config["output_dir"]) / "evaluation_report.txt"
    with open(report_path, "w") as f:
        f.write("RAPORT EVALUARE — VariantCallerCNN\n")
        f.write(f"Model: {config['model_path']}\n")
        f.write(f"Test set: HG004-Mother\n\n")
        f.write(f"F1 macro   : {f1:.4f}\n")
        f.write(f"Precision  : {prec:.4f}\n")
        f.write(f"Recall     : {rec:.4f}\n")
        if auc is not None:
            f.write(f"AUC-ROC    : {auc:.4f}\n")
        f.write(f"\n{report}\n")
        f.write("Confusion Matrix:\n")
        f.write(f"   {'':10s} " + "  ".join(f"{n:8s}" for n in CLASS_NAMES) + "\n")
        for i, row in enumerate(cm):
            vals = "  ".join(f"{v:8d}" for v in row)
            f.write(f"   {CLASS_NAMES[i]:10s} {vals}\n")

    print(f"\n💾 Raport salvat: {report_path}")

    # --- Export JSON pentru backend Java ---
    json_path = Path(config["output_dir"]) / "predictions_for_java.json"
    export_data = export_for_java(
        test_dataset, all_labels, all_preds, all_probs, str(json_path)
    )
    print(f"📤 Export Java: {json_path}")
    print(f"   {export_data['metadata']['n_variants']} variante identificate "
          f"(Het={export_data['metadata']['n_het']}, "
          f"Hom-Alt={export_data['metadata']['n_hom_alt']})")

    # --- Metrici JSON (pentru plots) ---
    metrics_path = Path(config["output_dir"]) / "metrics.json"
    with open(metrics_path, "w") as f:
        json.dump({
            "f1_macro":   round(f1, 6),
            "precision":  round(prec, 6),
            "recall":     round(rec, 6),
            "auc_roc":    round(auc, 6) if auc else None,
            "confusion_matrix": cm.tolist(),
            "class_names": CLASS_NAMES,
        }, f, indent=2)

    print(f"📊 Metrici JSON: {metrics_path}")
    print(f"\n💡 Backend Java poate consuma: {json_path}")
    print(f"   Format: {{chrom, pos, ref, alt, predicted_class, prob_*, confidence, is_variant}}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="VariantCallerCNN — Evaluare")
    parser.add_argument("--model",       type=str, default=None)
    parser.add_argument("--max_samples", type=int, default=None)
    parser.add_argument("--batch_size",  type=int, default=None)
    parser.add_argument("--num_workers", type=int, default=None)
    parser.add_argument("--output_dir",  type=str, default=None)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    config = DEFAULT_EVAL_CONFIG.copy()

    if args.model       is not None: config["model_path"]  = args.model
    if args.max_samples is not None: config["max_samples"] = args.max_samples
    if args.batch_size  is not None: config["batch_size"]  = args.batch_size
    if args.num_workers is not None: config["num_workers"] = args.num_workers
    if args.output_dir  is not None: config["output_dir"]  = args.output_dir

    run_evaluation(config)