"""
train.py — Pipeline antrenament VariantCallerCNN1D
Antrenează pe HG002 (train) și validează pe HG003.

Rulare:
    python train.py
    python train.py --epochs 30 --lr 1e-3 --batch_size 128
"""

import argparse
import json
import logging
import random
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, WeightedRandomSampler
from sklearn.metrics import (
    f1_score, precision_score, recall_score,
    classification_report, confusion_matrix,
)

from model import VariantCallerCNN1D
from dataset import GenomicDataset

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("train.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Config default
# ---------------------------------------------------------------------------
DEFAULT_CONFIG = {
    # === EXPERIMENT VALIDARE: HG003 (tată) → HG004 (mamă) — neînrudiți ===
    "bam_train": "data/HG003_Father/HG003.hiseq4000.wes-agilent.50x.dedup.grch38.bam",
    "vcf_train": "data/HG003_Father/HG003_exome.vcf",
    "bed_train": None,

    "bam_val":   "data/HG004_Mother/HG004.hiseq4000.wes-agilent.50x.dedup.grch38.bam",
    "vcf_val":   "data/HG004_Mother/HG004_exome.vcf",
    "bed_val":   None,

    "max_samples_train": 50_000,
    "max_samples_val":   10_000,

    "batch_size":  128,    # CNN 1D suportă batch-uri mai mari (input mic)
    "epochs":      30,
    "lr":          1e-3,
    "weight_decay": 1e-5,
    "num_workers": 4,
    "seed":        42,

    "output_dir":      "checkpoints_unrelated",
    "best_model_path": "checkpoints_unrelated/best_model.pth",
    "last_model_path": "checkpoints_unrelated/last_model.pth",
    "history_path":    "checkpoints_unrelated/history.json",
}

CLASS_NAMES = ["Ref", "Het", "Hom-Alt"]


def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def format_time(s: float) -> str:
    m, s = divmod(int(s), 60)
    h, m = divmod(m, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


# ---------------------------------------------------------------------------
# Evaluare
# ---------------------------------------------------------------------------

def evaluate(model, loader, criterion, device):
    model.eval()
    total_loss = 0.0
    all_preds, all_labels = [], []

    with torch.no_grad():
        for x, labels in loader:
            x, labels = x.to(device), labels.to(device)
            out = model(x)
            total_loss += criterion(out, labels).item()
            _, pred = torch.max(out, 1)
            all_preds.extend(pred.cpu().numpy())
            all_labels.extend(labels.cpu().numpy())

    avg_loss = total_loss / max(len(loader), 1)
    f1   = f1_score(all_labels, all_preds, average="macro", zero_division=0)
    prec = precision_score(all_labels, all_preds, average="macro", zero_division=0)
    rec  = recall_score(all_labels, all_preds, average="macro", zero_division=0)
    report = classification_report(
        all_labels, all_preds,
        labels=[0, 1, 2],
        target_names=CLASS_NAMES,
        zero_division=0,
    )
    cm = confusion_matrix(all_labels, all_preds, labels=[0, 1, 2])
    return avg_loss, f1, prec, rec, report, cm


# ---------------------------------------------------------------------------
# Loop principal
# ---------------------------------------------------------------------------

def train_model(config: dict):
    set_seed(config["seed"])
    Path(config["output_dir"]).mkdir(parents=True, exist_ok=True)

    # --- Device ---
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info(f"Device: {device}" + (
        f" ({torch.cuda.get_device_name(0)})" if device.type == "cuda" else ""
    ))
    print(f"\n🚀 Antrenament pe: {device}", end="")
    if device.type == "cuda":
        print(f"  ({torch.cuda.get_device_name(0)})")
    else:
        print("  (CPU)")

    # --- Dataset ---
    print("\n📦 Train (HG002-Son)...")
    train_ds = GenomicDataset(
        bam_path    = config["bam_train"],
        vcf_path    = config["vcf_train"],
        bed_path    = config.get("bed_train"),
        max_samples = config["max_samples_train"],
        seed        = config["seed"],
    )

    print("\n📦 Validare (HG003-Father)...")
    val_ds = GenomicDataset(
        bam_path    = config["bam_val"],
        vcf_path    = config["vcf_val"],
        bed_path    = config.get("bed_val"),
        max_samples = config["max_samples_val"],
        seed        = config["seed"],
    )

    # --- WeightedRandomSampler pentru batch-uri echilibrate ---
    label_counts = defaultdict(int)
    for dp in train_ds.data_points:
        label_counts[dp["label"]] += 1

    sample_weights = torch.tensor([
        1.0 / max(label_counts[dp["label"]], 1)
        for dp in train_ds.data_points
    ], dtype=torch.float32)

    sampler = WeightedRandomSampler(
        weights     = sample_weights,
        num_samples = len(sample_weights),
        replacement = True,
    )
    print(f"\n   Distribuție train: Ref={label_counts[0]} | "
          f"Het={label_counts[1]} | Hom={label_counts[2]}")

    train_loader = DataLoader(
        train_ds,
        batch_size  = config["batch_size"],
        sampler     = sampler,
        num_workers = config["num_workers"],
        pin_memory  = device.type == "cuda",
        persistent_workers = config["num_workers"] > 0,
    )
    val_loader = DataLoader(
        val_ds,
        batch_size  = config["batch_size"],
        shuffle     = False,
        num_workers = config["num_workers"],
        pin_memory  = device.type == "cuda",
        persistent_workers = config["num_workers"] > 0,
    )

    # --- Model ---
    model = VariantCallerCNN1D().to(device)
    total_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"\n🧠 VariantCallerCNN1D — {total_params:,} parametri antrenabili")

    # --- Loss + Optimizer + Scheduler ---
    weights = torch.tensor([1.0, 1.0, 1.0], dtype=torch.float32).to(device)
    print(f"⚖️  Ponderi clase: Ref={weights[0]:.2f} | "
          f"Het={weights[1]:.2f} | Hom={weights[2]:.2f}")

    criterion = nn.CrossEntropyLoss(weight=weights, label_smoothing=0.05)
    optimizer = optim.AdamW(
        model.parameters(),
        lr           = config["lr"],
        weight_decay = config["weight_decay"],
    )
    scheduler = optim.lr_scheduler.OneCycleLR(
        optimizer,
        max_lr            = config["lr"],
        epochs            = config["epochs"],
        steps_per_epoch   = len(train_loader),
        pct_start         = 0.1,
        div_factor        = 10,
        final_div_factor  = 100,
    )

    # --- Training loop ---
    best_f1    = 0.0
    no_improve = 0
    PATIENCE   = 10
    history    = {"train_loss": [], "val_loss": [], "f1": [], "lr": []}

    print(f"\n⏳ Antrenament: {config['epochs']} epoci, "
          f"batch={config['batch_size']}, lr={config['lr']}")
    print("=" * 70)

    for epoch in range(config["epochs"]):
        t0 = time.time()

        # --- TRAIN ---
        model.train()
        train_loss  = 0.0
        train_steps = 0

        for batch_idx, (x, labels) in enumerate(train_loader):
            x, labels = x.to(device), labels.to(device)
            optimizer.zero_grad()
            loss = criterion(model(x), labels)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            scheduler.step()
            train_loss  += loss.item()
            train_steps += 1

            if (batch_idx + 1) % 50 == 0:
                print(f"  [{batch_idx+1:4d}/{len(train_loader)}] "
                      f"loss={train_loss/train_steps:.4f}", end="\r")

        avg_train_loss = train_loss / max(len(train_loader), 1)

        # --- VAL ---
        avg_val_loss, f1, prec, rec, report, cm = evaluate(
            model, val_loader, criterion, device
        )

        current_lr = optimizer.param_groups[0]["lr"]
        elapsed    = time.time() - t0

        history["train_loss"].append(round(avg_train_loss, 6))
        history["val_loss"].append(round(avg_val_loss, 6))
        history["f1"].append(round(f1, 6))
        history["lr"].append(current_lr)

        print(f"\n📊 Epoca {epoch+1:02d}/{config['epochs']}  "
              f"[{format_time(elapsed)}]  lr={current_lr:.2e}")
        print(f"   Train Loss : {avg_train_loss:.4f}")
        print(f"   Val   Loss : {avg_val_loss:.4f}")
        print(f"   F1 macro   : {f1:.4f}  |  Precision: {prec:.4f}  |  Recall: {rec:.4f}")
        print(f"   Confusion Matrix:")
        for i, row in enumerate(cm):
            print(f"     {CLASS_NAMES[i]:8s}: {row}")

        if (epoch + 1) % 5 == 0:
            print(f"\n   📋 Raport detaliat (epoca {epoch+1}):\n{report}")

        # Salvare best model după F1
        if f1 > best_f1:
            best_f1    = f1
            no_improve = 0
            torch.save(
                {
                    "epoch":              epoch + 1,
                    "model_state_dict":   model.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                    "val_loss":           avg_val_loss,
                    "f1":                 best_f1,
                    "config":             config,
                    "class_names":        CLASS_NAMES,
                },
                config["best_model_path"],
            )
            print(f"   💾 Best model salvat (F1={best_f1:.4f})")
        else:
            no_improve += 1

        print()

        if no_improve >= PATIENCE:
            print(f"⏹️  Early stopping la epoca {epoch+1} "
                  f"(nicio îmbunătățire în {PATIENCE} epoci)")
            break

    # Salvare model final + istoricul
    torch.save(model.state_dict(), config["last_model_path"])

    with open(config["history_path"], "w") as f:
        json.dump(history, f, indent=2)

    config_path = Path(config["output_dir"]) / "config_used.json"
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)

    print("=" * 70)
    print(f"\n✅ Antrenament finalizat!")
    print(f"   Best model : {config['best_model_path']}  (F1={best_f1:.4f})")
    print(f"   Last model : {config['last_model_path']}")
    print(f"   Istoricul  : {config['history_path']}")
    print(f"\n💡 Pasul următor: python evaluate.py  (test set HG004-Mother)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="VariantCallerCNN1D — Antrenament")
    parser.add_argument("--epochs",      type=int,   default=None)
    parser.add_argument("--lr",          type=float, default=None)
    parser.add_argument("--batch_size",  type=int,   default=None)
    parser.add_argument("--max_train",   type=int,   default=None)
    parser.add_argument("--max_val",     type=int,   default=None)
    parser.add_argument("--num_workers", type=int,   default=None)
    parser.add_argument("--seed",        type=int,   default=None)
    parser.add_argument("--output_dir",  type=str,   default=None)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    config = DEFAULT_CONFIG.copy()

    if args.epochs      is not None: config["epochs"]             = args.epochs
    if args.lr          is not None: config["lr"]                 = args.lr
    if args.batch_size  is not None: config["batch_size"]         = args.batch_size
    if args.max_train   is not None: config["max_samples_train"]  = args.max_train
    if args.max_val     is not None: config["max_samples_val"]    = args.max_val
    if args.num_workers is not None: config["num_workers"]        = args.num_workers
    if args.seed        is not None: config["seed"]               = args.seed
    if args.output_dir  is not None:
        config["output_dir"]      = args.output_dir
        config["best_model_path"] = f"{args.output_dir}/best_model.pth"
        config["last_model_path"] = f"{args.output_dir}/last_model.pth"
        config["history_path"]    = f"{args.output_dir}/history.json"

    train_model(config)