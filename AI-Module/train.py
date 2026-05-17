"""
train.py — Pipeline complet de antrenament
Antrenează VariantCallerCNN pe HG002 (train) și validează pe HG003.

Rulare:
    python train.py
    python train.py --epochs 30 --lr 5e-5 --batch_size 64
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
from torch.utils.data import DataLoader
from sklearn.metrics import (
    f1_score, precision_score, recall_score, classification_report,
    confusion_matrix,
)

from model import VariantCallerCNN
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
    # --- Date ---
    "bam_train": "data/HG002_Son/HG002.hiseq4000.wes-agilent.50x.dedup.grch38.bam",
    "vcf_train": "data/HG002_Son/HG002_exome.vcf",
    "bed_train": None,   # Ref generat din offset față de variante

    "bam_val":   "data/HG003_Father/HG003.hiseq4000.wes-agilent.50x.dedup.grch38.bam",
    "vcf_val":   "data/HG003_Father/HG003_exome.vcf",
    "bed_val":   None,

    # --- Dimensiuni ---
    "max_samples_train": 50_000,
    "max_samples_val":   10_000,

    # --- Hiperparametri ---
    "batch_size":  32,
    "epochs":      30,
    "lr":          1e-3,
    "weight_decay": 1e-5,
    "num_workers": 4,
    "seed":        42,

    # --- Ieșiri ---
    "output_dir":      "checkpoints",
    "best_model_path": "checkpoints/best_model.pth",
    "last_model_path": "checkpoints/last_model.pth",
    "history_path":    "checkpoints/history.json",
}

CLASS_NAMES = ["Ref", "Het", "Hom-Alt"]


# ---------------------------------------------------------------------------
# Utilități
# ---------------------------------------------------------------------------

def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def format_time(seconds: float) -> str:
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


# ---------------------------------------------------------------------------
# Evaluare
# ---------------------------------------------------------------------------

def evaluate(model, loader, criterion, device):
    """
    Returnează:
        avg_loss, f1_macro, precision_macro, recall_macro,
        classification_report_str, confusion_matrix_array
    """
    model.eval()
    total_loss = 0.0
    all_preds  = []
    all_labels = []

    with torch.no_grad():
        for images, labels in loader:
            images, labels = images.to(device), labels.to(device)
            outputs = model(images)
            total_loss += criterion(outputs, labels).item()
            _, predicted = torch.max(outputs, 1)
            all_preds.extend(predicted.cpu().numpy())
            all_labels.extend(labels.cpu().numpy())

    avg_loss = total_loss / len(loader)
    f1   = f1_score(all_labels, all_preds, average="macro", zero_division=0)
    prec = precision_score(all_labels, all_preds, average="macro", zero_division=0)
    rec  = recall_score(all_labels, all_preds, average="macro", zero_division=0)
    report = classification_report(
        all_labels, all_preds, target_names=CLASS_NAMES, zero_division=0
    )
    cm = confusion_matrix(all_labels, all_preds, labels=[0, 1, 2])
    return avg_loss, f1, prec, rec, report, cm


# ---------------------------------------------------------------------------
# Loop principal
# ---------------------------------------------------------------------------

def train_model(config: dict):
    set_seed(config["seed"])

    # Creare director output
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
        print("  (CPU — poate fi lent!)")

    # --------------------------------------------------------------------------
    # 1. Dataset-uri
    # --------------------------------------------------------------------------
    print("\n📦 Încărcăm datele de antrenament (HG002-Son)...")
    train_dataset = GenomicDataset(
        bam_path    = config["bam_train"],
        vcf_path    = config["vcf_train"],
        bed_path    = config.get("bed_train"),
        max_samples = config["max_samples_train"],
        seed        = config["seed"],
    )

    print("\n📦 Încărcăm datele de validare (HG003-Father)...")
    val_dataset = GenomicDataset(
        bam_path    = config["bam_val"],
        vcf_path    = config["vcf_val"],
        bed_path    = config.get("bed_val"),
        max_samples = config["max_samples_val"],
        seed        = config["seed"],
    )

    # --------------------------------------------------------------------------
    # 2. Model
    # --------------------------------------------------------------------------
    model = VariantCallerCNN().to(device)
    total_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"\n🧠 Model: VariantCallerCNN — {total_params:,} parametri antrenabili")

    # --------------------------------------------------------------------------
    # 3. WeightedRandomSampler — fiecare batch are ~1/3 din fiecare clasă
    # --------------------------------------------------------------------------
    label_counts = defaultdict(int)
    for dp in train_dataset.data_points:
        label_counts[dp["label"]] += 1

    # Greutate per exemplu = inversul frecvenței clasei sale
    sample_weights = torch.tensor([
        1.0 / max(label_counts[dp["label"]], 1)
        for dp in train_dataset.data_points
    ], dtype=torch.float32)

    sampler = torch.utils.data.WeightedRandomSampler(
        weights     = sample_weights,
        num_samples = len(sample_weights),
        replacement = True,
    )
    print(f"   Distribuție train: Ref={label_counts[0]} | Het={label_counts[1]} | Hom={label_counts[2]}")
    print(f"   WeightedSampler activ — fiecare batch va fi ~echilibrat per clasă")

    train_loader = DataLoader(
        train_dataset,
        batch_size  = config["batch_size"],
        sampler     = sampler,          # înlocuiește shuffle=True
        num_workers = config["num_workers"],
        pin_memory  = device.type == "cuda",
        persistent_workers = config["num_workers"] > 0,
    )
    val_loader = DataLoader(
        val_dataset,
        batch_size  = config["batch_size"],
        shuffle     = False,
        num_workers = config["num_workers"],
        pin_memory  = device.type == "cuda",
        persistent_workers = config["num_workers"] > 0,
    )

    total = sum(label_counts.values())
    # Ponderi mai agresive pe Het — clasa cea mai greu de învățat
    weights = torch.tensor([1.0, 1.5, 1.0], dtype=torch.float32).to(device)
    print(f"⚖️  Ponderi clase: Ref={weights[0]:.2f} | Het={weights[1]:.2f} | Hom={weights[2]:.2f}")

    criterion = nn.CrossEntropyLoss(weight=weights, label_smoothing=0.05)
    optimizer = optim.AdamW(
        model.parameters(),
        lr           = config["lr"],
        weight_decay = config["weight_decay"],
    )
    # OneCycleLR — warmup rapid + decay, mult mai bun pentru convergență
    scheduler = optim.lr_scheduler.OneCycleLR(
        optimizer,
        max_lr      = config["lr"],
        epochs      = config["epochs"],
        steps_per_epoch = len(train_loader),
        pct_start   = 0.1,   # 10% warmup
        div_factor  = 10,    # lr_start = max_lr/10
        final_div_factor = 100,
    )

    # --------------------------------------------------------------------------
    # 4. Loop de antrenament
    # --------------------------------------------------------------------------
    best_val_loss = float("inf")
    best_f1       = 0.0
    history       = {"train_loss": [], "val_loss": [], "f1": [], "lr": []}
    no_improve    = 0
    PATIENCE      = 10  # mai mult timp să iasă din platou

    print(f"\n⏳ Antrenament: {config['epochs']} epoci, "
          f"batch={config['batch_size']}, lr={config['lr']}\n")
    print("=" * 70)

    for epoch in range(config["epochs"]):
        epoch_start = time.time()

        # --- TRAIN ---
        model.train()
        train_loss  = 0.0
        train_steps = 0

        for batch_idx, (images, labels) in enumerate(train_loader):
            images, labels = images.to(device), labels.to(device)
            optimizer.zero_grad()
            loss = criterion(model(images), labels)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            scheduler.step()   # OneCycleLR: pas per batch
            train_loss  += loss.item()
            train_steps += 1

            # Progress bar simplu
            if (batch_idx + 1) % 50 == 0:
                print(f"  [{batch_idx+1:4d}/{len(train_loader)}] "
                      f"loss={train_loss/train_steps:.4f}", end="\r")

        avg_train_loss = train_loss / len(train_loader)

        # --- VALIDARE ---
        avg_val_loss, f1, prec, rec, report, cm = evaluate(
            model, val_loader, criterion, device
        )
        scheduler_step_done = True  # OneCycleLR se face per batch

        current_lr = optimizer.param_groups[0]["lr"]
        elapsed    = time.time() - epoch_start

        # Istoricul
        history["train_loss"].append(round(avg_train_loss, 6))
        history["val_loss"].append(round(avg_val_loss, 6))
        history["f1"].append(round(f1, 6))
        history["lr"].append(current_lr)

        # --- Afișare epocă ---
        print(f"\n📊 Epoca {epoch+1:02d}/{config['epochs']}  "
              f"[{format_time(elapsed)}]  lr={current_lr:.2e}")
        print(f"   Train Loss : {avg_train_loss:.4f}")
        print(f"   Val   Loss : {avg_val_loss:.4f}")
        print(f"   F1 macro   : {f1:.4f}  |  Precision: {prec:.4f}  |  Recall: {rec:.4f}")
        print(f"   Confusion Matrix:")
        for i, row in enumerate(cm):
            print(f"     {CLASS_NAMES[i]:8s}: {row}")

        # Raport detaliat la fiecare 5 epoci
        if (epoch + 1) % 5 == 0:
            print(f"\n   📋 Raport detaliat (epoca {epoch+1}):\n{report}")

        # --- Salvare cel mai bun model (criteriu: F1 macro) ---
        if f1 > best_f1:
            best_val_loss = avg_val_loss
            best_f1       = f1
            no_improve    = 0
            torch.save(
                {
                    "epoch":              epoch + 1,
                    "model_state_dict":   model.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                    "val_loss":           best_val_loss,
                    "f1":                 best_f1,
                    "config":             config,
                    "class_names":        CLASS_NAMES,
                },
                config["best_model_path"],
            )
            print(f"   💾 Best model salvat (val_loss={best_val_loss:.4f}, F1={best_f1:.4f})")
        else:
            no_improve += 1

        print()

        # Early stopping
        if no_improve >= PATIENCE:
            print(f"⏹️  Early stopping la epoca {epoch+1} "
                  f"(nicio îmbunătățire în {PATIENCE} epoci)")
            break

    # --------------------------------------------------------------------------
    # 5. Salvare model final + istoricul
    # --------------------------------------------------------------------------
    torch.save(model.state_dict(), config["last_model_path"])

    with open(config["history_path"], "w") as f:
        json.dump(history, f, indent=2)

    # Salvăm și config-ul folosit
    config_path = Path(config["output_dir"]) / "config_used.json"
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)

    print("=" * 70)
    print(f"\n✅ Antrenament finalizat!")
    print(f"   Best model : {config['best_model_path']}  "
          f"(val_loss={best_val_loss:.4f}, F1={best_f1:.4f})")
    print(f"   Last model : {config['last_model_path']}")
    print(f"   Istoricul  : {config['history_path']}")
    print(f"\n💡 Pasul următor: python evaluate.py  (test set HG004-Mother)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="VariantCallerCNN — Antrenament")
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

    # Override din CLI
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