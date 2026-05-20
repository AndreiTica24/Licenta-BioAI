"""
dataset.py — GenomicDataset pentru CNN 1D
Variant calling pe secvențe ADN brute (200bp în jurul fiecărei poziții candidate).

Encoding (6 canale × 200 poziții):
  Canal 0: A (one-hot, 0 sau 1)
  Canal 1: C
  Canal 2: G
  Canal 3: T
  Canal 4: depth normalizat (0-1) per coloană
  Canal 5: AF local (alt_count/depth) per coloană

Această reprezentare e standardul din literatura genomic deep learning:
DeepBind, DeepSEA, Basset folosesc CNN 1D pe secvențe one-hot encoded.

Label: 0=Ref, 1=Het, 2=Hom-Alt
"""

import os
import pickle
import hashlib
import random
import logging
from pathlib import Path
from typing import Optional, List, Dict, Tuple
from collections import Counter

import numpy as np
import torch
from torch.utils.data import Dataset

try:
    import pysam
except ImportError:
    raise ImportError("Instalează pysam: pip install pysam")

try:
    import cyvcf2
except ImportError:
    raise ImportError("Instalează cyvcf2: pip install cyvcf2")

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constante encoding
# ---------------------------------------------------------------------------
WINDOW_SIZE   = 200   # fereastră ± 100bp în jurul poziției candidate
N_CHANNELS    = 6     # A, C, G, T, depth, AF
MAX_DEPTH     = 200   # depth de saturație pentru normalizare
MIN_MAPPING_Q = 1     # filtru minim — eliminăm doar read-urile MQ=0
MIN_BASE_Q    = 0     # fără filtru pe BQ

# One-hot indices pentru nucleotide
BASE_IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def _genotype_to_label(gt: Tuple[int, ...]) -> int:
    """Convertește genotype VCF la label 0/1/2."""
    alleles = [a for a in gt if a is not None and a >= 0]
    if not alleles:
        return -1
    if all(a == 0 for a in alleles):
        return 0
    if len(set(alleles)) > 1:
        return 1
    return 2


def _encode_window(bam: pysam.AlignmentFile,
                   chrom: str,
                   pos: int,
                   window: int = WINDOW_SIZE) -> np.ndarray:
    """
    Construiește o reprezentare CNN 1D pentru fereastra centrată pe `pos`.

    Pentru fiecare coloană din fereastră, calculează din pileup:
      - one-hot encoding al bazei majoritare
      - depth normalizat
      - allele frequency local (procent reads cu baza minoritară)

    Returneaza tensor (6, 200).
    """
    half  = window // 2
    start = max(0, pos - half)
    end   = start + window

    # Inițializăm tensor-ul
    img = np.zeros((N_CHANNELS, window), dtype=np.float32)

    # Normalizare contig
    bam_refs = set(bam.references)
    if chrom not in bam_refs:
        alt_c = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        if alt_c in bam_refs:
            chrom = alt_c
        else:
            return img  # contig necunoscut, returnăm zeros

    try:
        for col in bam.pileup(chrom, start, end,
                               truncate=True,
                               min_base_quality=MIN_BASE_Q,
                               min_mapping_quality=MIN_MAPPING_Q,
                               stepper="all",
                               ignore_overlaps=False):

            col_pos = col.reference_pos
            if col_pos < start or col_pos >= end:
                continue
            col_idx = col_pos - start
            if col_idx < 0 or col_idx >= window:
                continue

            # Numărăm bazele la această coloană
            base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            depth = 0

            for r in col.pileups:
                if r.is_del or r.is_refskip:
                    continue
                try:
                    if r.alignment.query_sequence is None:
                        continue
                    base = r.alignment.query_sequence[r.query_position].upper()
                    if base in base_counts:
                        base_counts[base] += 1
                        depth += 1
                except (IndexError, TypeError):
                    continue

            if depth == 0:
                continue

            # Baza majoritară primește one-hot=1
            major_base = max(base_counts, key=base_counts.get)
            img[BASE_IDX[major_base], col_idx] = 1.0

            # Depth normalizat
            img[4, col_idx] = min(depth, MAX_DEPTH) / MAX_DEPTH

            # AF local = (reads non-majoritare) / depth
            non_major = depth - base_counts[major_base]
            img[5, col_idx] = non_major / depth

    except (ValueError, KeyError, AssertionError):
        pass

    return img


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------

class GenomicDataset(Dataset):
    """
    Dataset PyTorch pentru CNN 1D variant calling.

    Parametri
    ---------
    bam_path    : calea la fișierul BAM (indexat .bai)
    vcf_path    : calea la VCF-ul de referință (benchmark)
    bed_path    : (NEFOLOSIT — păstrat pentru compatibilitate API)
    max_samples : numărul maxim de exemple (None = toate)
    seed        : seed pentru reproducibilitate
    """

    def __init__(self,
                 bam_path:    str,
                 vcf_path:    str,
                 bed_path:    Optional[str] = None,
                 max_samples: Optional[int] = None,
                 seed:        int = 42):
        self.bam_path    = bam_path
        self.vcf_path    = vcf_path
        self.bed_path    = bed_path
        self.max_samples = max_samples
        self.seed        = seed
        self.data_points: List[Dict] = []
        self._build_index()

    # ------------------------------------------------------------------
    def _cache_path(self) -> str:
        """Cale unică cache bazată pe parametri."""
        key = f"{self.vcf_path}|{self.max_samples}|{self.seed}|cnn1d_v1"
        key_hash = hashlib.md5(key.encode()).hexdigest()[:10]
        vcf_stem = os.path.splitext(os.path.basename(self.vcf_path))[0]
        cache_dir = os.path.join(os.path.dirname(self.vcf_path), ".cache")
        os.makedirs(cache_dir, exist_ok=True)
        return os.path.join(cache_dir, f"{vcf_stem}_{key_hash}.pkl")

    # ------------------------------------------------------------------
    def _build_index(self):
        """Parsează VCF + generează exemple Ref din offset."""
        cache_path = self._cache_path()

        if os.path.exists(cache_path):
            print(f"   ⚡ Cache găsit! Încărcare rapidă din: {os.path.basename(cache_path)}")
            with open(cache_path, "rb") as f:
                self.data_points = pickle.load(f)
            counts = Counter(p["label"] for p in self.data_points)
            print(f"   ✅ {len(self.data_points)} exemple din cache | "
                  f"Ref={counts[0]} | Het={counts[1]} | Hom={counts[2]}")
            return

        print(f"   🔍 Prima rulare: parsăm VCF...")

        # Detectăm formatul contigurilor din BAM
        _bam = pysam.AlignmentFile(self.bam_path, "rb")
        bam_refs    = set(_bam.references)
        bam_has_chr = any(r.startswith("chr") for r in bam_refs)
        _bam.close()

        vcf    = cyvcf2.VCF(self.vcf_path)
        points = []

        for variant in vcf:
            chrom = variant.CHROM
            if bam_has_chr and not chrom.startswith("chr"):
                chrom = "chr" + chrom
            elif not bam_has_chr and chrom.startswith("chr"):
                chrom = chrom[3:]

            pos = variant.POS - 1  # 0-indexed

            # Doar SNV-uri simple (nu indel, nu multi-alelic)
            if len(variant.ALT) > 1:
                continue
            if len(variant.REF) != 1 or len(variant.ALT[0]) != 1:
                continue

            gt = variant.genotypes[0][:2]
            label = _genotype_to_label(tuple(gt))
            if label == -1:
                continue

            points.append({
                "chrom": chrom,
                "pos":   pos,
                "ref":   variant.REF,
                "alt":   variant.ALT[0],
                "label": label,
            })

        vcf.close()

        # Generăm exemple Ref din offset (poziții vecine fără variantă)
        rng = random.Random(self.seed)
        ref_budget  = len(points) // 2
        variant_pos = {(p["chrom"], p["pos"]) for p in points}
        ref_points  = []

        for dp in rng.sample(points, min(ref_budget * 3, len(points))):
            if len(ref_points) >= ref_budget:
                break
            for offset in [25, -25, 50, -50, 75, -75, 100, -100]:
                new_pos = dp["pos"] + offset
                if new_pos < 0:
                    continue
                if (dp["chrom"], new_pos) not in variant_pos:
                    ref_points.append({
                        "chrom": dp["chrom"],
                        "pos":   new_pos,
                        "ref":   "N",
                        "alt":   ".",
                        "label": 0,
                    })
                    variant_pos.add((dp["chrom"], new_pos))
                    break

        points.extend(ref_points)

        # Shuffle + limitare
        rng = random.Random(self.seed)
        rng.shuffle(points)
        if self.max_samples is not None:
            points = points[:self.max_samples]
        self.data_points = points

        print(f"   💾 Salvăm cache: {os.path.basename(cache_path)}")
        with open(cache_path, "wb") as f:
            pickle.dump(self.data_points, f)

        counts = Counter(p["label"] for p in self.data_points)
        print(f"   ✅ {len(self.data_points)} exemple | "
              f"Ref={counts[0]} | Het={counts[1]} | Hom={counts[2]}")

    # ------------------------------------------------------------------
    def __len__(self) -> int:
        return len(self.data_points)

    # ------------------------------------------------------------------
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        dp = self.data_points[idx]

        bam = pysam.AlignmentFile(self.bam_path, "rb")
        x   = _encode_window(bam, dp["chrom"], dp["pos"], window=WINDOW_SIZE)
        bam.close()

        tensor = torch.from_numpy(x)  # (6, 200)
        label  = torch.tensor(dp["label"], dtype=torch.long)
        return tensor, label