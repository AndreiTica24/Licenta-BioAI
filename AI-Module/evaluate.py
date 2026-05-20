"""
evaluate.py — Evaluare finală pe HG004 (test set nevăzut) + export JSON pentru Java.

Rulare:
    python evaluate.py
    python evaluate.py --model checkpoints/best_model.pth --max_samples 5000
"""

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader
from sklearn.metrics import (
    f1_score, precision_score, recall_score, accuracy_score,
    classification_report, confusion_matrix, roc_auc_score,
)

from model import VariantCallerCNN1D
from dataset import GenomicDataset

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

CLASS_NAMES = ["Ref", "Het", "Hom-Alt"]

DEFAULT_EVAL_CONFIG = {
    "bam_test": "data/HG004_Mother/HG004.hiseq4000.wes-agilent.50x.dedup.grch38.bam",
    "vcf_test": "data/HG004_Mother/HG004_exome.vcf",
    "bed_test": None,
    "max_samples": 10_000,
    "batch_size":  128,
    "num_workers": 4,
    "seed":        42,
    "model_path":  "checkpoints/best_model.pth",
    "output_dir":  "results",
}


def evaluate_with_probs(model, loader, device):
    """Returnează predicții + labels + probabilități softmax."""
    model.eval()
    all_preds, all_labels, all_probs = [], [], []

    with torch.no_grad():
        for x, labels in loader:
            x = x.to(device)
            logits = model(x)
            probs  = F.softmax(logits, dim=1)
            _, pred = torch.max(logits, 1)

            all_preds.extend(pred.cpu().numpy().tolist())
            all_labels.extend(labels.numpy().tolist())
            all_probs.extend(probs.cpu().numpy().tolist())

    return all_labels, all_preds, all_probs


def export_for_java(dataset, all_labels, all_preds, all_probs, path: str):
    """Export JSON consumabil de backend-ul Java cu ClinVar."""
    variants = []
    for i, dp in enumerate(dataset.data_points):
        if i >= len(all_preds):
            break

        pred = all_preds[i]
        true = all_labels[i]
        p    = all_probs[i]

        variants.append({
            "chrom":           dp["chrom"],
            "pos":             dp["pos"] + 1,   # 1-indexed VCF standard
            "ref":             dp["ref"],
            "alt":             dp["alt"],
            "predicted_class": CLASS_NAMES[pred],
            "predicted_label": pred,
            "true_class":      CLASS_NAMES[true],
            "true_label":      true,
            "prob_ref":        round(p[0], 4),
            "prob_het":        round(p[1], 4),
            "prob_hom_alt":    round(p[2], 4),
            "confidence":      round(max(p), 4),
            "is_variant":      pred > 0,
        })

    n_variants = sum(1 for v in variants if v["is_variant"])
    n_het      = sum(1 for v in variants if v["predicted_class"] == "Het")
    n_hom      = sum(1 for v in variants if v["predicted_class"] == "Hom-Alt")

    output = {
        "metadata": {
            "total_positions": len(variants),
            "n_variants":      n_variants,
            "n_het":           n_het,
            "n_hom_alt":       n_hom,
            "model_version":   "VariantCallerCNN1D-v1",
        },
        "variants": variants,
    }

    with open(path, "w") as f:
        json.dump(output, f, indent=2)

    logger.info(f"Export JSON: {path} "
                f"({len(variants)} poziții, {n_variants} variante)")
    return output


def run_evaluation(config: dict):
    Path(config["output_dir"]).mkdir(parents=True, exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"\n🔬 Evaluare pe: {device}")
    if device.type == "cuda":
        print(f"   GPU: {torch.cuda.get_device_name(0)}")

    # --- Model ---
    print(f"\n📥 Încărcăm modelul: {config['model_path']}")
    model = VariantCallerCNN1D.load_from_checkpoint(
        config["model_path"], device=str(device)
    )
    model.eval()

    # --- Dataset test ---
    print(f"\n📦 Test (HG004-Mother)...")
    test_ds = GenomicDataset(
        bam_path    = config["bam_test"],
        vcf_path    = config["vcf_test"],
        bed_path    = config.get("bed_test"),
        max_samples = config["max_samples"],
        seed        = config["seed"],
    )

    test_loader = DataLoader(
        test_ds,
        batch_size  = config["batch_size"],
        shuffle     = False,
        num_workers = config["num_workers"],
        pin_memory  = device.type == "cuda",
    )

    # --- Evaluare ---
    print(f"\n⏳ Evaluăm {len(test_ds)} exemple...")
    all_labels, all_preds, all_probs = evaluate_with_probs(model, test_loader, device)

    # --- Metrici ---
    acc  = accuracy_score(all_labels, all_preds)
    f1   = f1_score(all_labels, all_preds, average="macro", zero_division=0)
    prec = precision_score(all_labels, all_preds, average="macro", zero_division=0)
    rec  = recall_score(all_labels, all_preds, average="macro", zero_division=0)
    try:
        probs_np = np.array(all_probs)
        auc = roc_auc_score(all_labels, probs_np,
                            multi_class="ovr", average="macro")
    except Exception:
        auc = None
    report = classification_report(
        all_labels, all_preds,
        labels=[0, 1, 2], target_names=CLASS_NAMES, zero_division=0,
    )
    cm = confusion_matrix(all_labels, all_preds, labels=[0, 1, 2])

    # --- Raport ---
    print("\n" + "=" * 70)
    print("📊 REZULTATE FINALE — HG004-Mother (Test Set Nevăzut)")
    print("=" * 70)
    print(f"   Accuracy   : {acc:.4f}")
    print(f"   F1 macro   : {f1:.4f}")
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

    # --- Raport text ---
    report_path = Path(config["output_dir"]) / "evaluation_report.txt"
    with open(report_path, "w") as f:
        f.write("RAPORT EVALUARE — VariantCallerCNN1D\n")
        f.write(f"Model: {config['model_path']}\n")
        f.write(f"Test set: HG004-Mother\n\n")
        f.write(f"Accuracy   : {acc:.4f}\n")
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
    out = export_for_java(test_ds, all_labels, all_preds, all_probs, str(json_path))
    print(f"📤 Export Java: {json_path}")
    print(f"   {out['metadata']['n_variants']} variante identificate "
          f"(Het={out['metadata']['n_het']}, "
          f"Hom-Alt={out['metadata']['n_hom_alt']})")

    # --- Metrici JSON ---
    metrics_path = Path(config["output_dir"]) / "metrics.json"
    with open(metrics_path, "w") as f:
        json.dump({
            "accuracy":  round(acc, 6),
            "f1_macro":  round(f1, 6),
            "precision": round(prec, 6),
            "recall":    round(rec, 6),
            "auc_roc":   round(auc, 6) if auc else None,
            "confusion_matrix": cm.tolist(),
            "class_names":      CLASS_NAMES,
        }, f, indent=2)

    print(f"📊 Metrici JSON: {metrics_path}")
    print(f"\n💡 Backend Java va consuma: {json_path}")
    print(f"   Format: {{chrom, pos, ref, alt, predicted_class, prob_*, confidence}}")


def parse_args():
    p = argparse.ArgumentParser(description="VariantCallerCNN1D — Evaluare")
    p.add_argument("--model",       type=str, default=None)
    p.add_argument("--max_samples", type=int, default=None)
    p.add_argument("--batch_size",  type=int, default=None)
    p.add_argument("--num_workers", type=int, default=None)
    p.add_argument("--output_dir",  type=str, default=None)
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    config = DEFAULT_EVAL_CONFIG.copy()
    if args.model       is not None: config["model_path"]  = args.model
    if args.max_samples is not None: config["max_samples"] = args.max_samples
    if args.batch_size  is not None: config["batch_size"]  = args.batch_size
    if args.num_workers is not None: config["num_workers"] = args.num_workers
    if args.output_dir  is not None: config["output_dir"]  = args.output_dir
    run_evaluation(config)