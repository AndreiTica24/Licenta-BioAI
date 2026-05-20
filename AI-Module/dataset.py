"""
dataset.py v3 — GenomicDataset cu encoder corect

Schimbări majore față de v2:
1. Pileup DENS — fără rânduri goale (repetare ciclică dacă depth < height)
2. Read-uri sortate ALT/REF balansat (50/50 când depth e mare)
3. CIGAR parsing direct — mapăm fiecare bază a read-ului la coloana corectă
4. Eliminat AF "trișat" pe coloana centrală — semnalul vine din pattern-ul real
5. Filtru flag_filter=3848 (exclude duplicates, unmapped, secondary, supplementary)

Canale per pixel (6, height=100, width=100):
  0: baza nucleotidică (A=0.25, C=0.5, G=0.75, T=1.0)
  1: este baza alternativă față de referință? (0/1)
  2: calitate bază normalizată (0-1)
  3: calitate mapping normalizată (0-1)
  4: direcție strand (0=reverse, 1=forward)
  5: este deletion? (0/1)
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
PILEUP_WIDTH   = 100
PILEUP_HEIGHT  = 100
N_CHANNELS     = 6
MAX_BASE_QUAL  = 40
MAX_MAP_QUAL   = 60

BASE_ENC = {'A': 0.25, 'C': 0.50, 'G': 0.75, 'T': 1.00,
            'N': 0.00, '-': 0.00, '*': 0.00}


def _genotype_to_label(gt: Tuple[int, ...]) -> int:
    alleles = [a for a in gt if a is not None and a >= 0]
    if not alleles:
        return -1
    if all(a == 0 for a in alleles):
        return 0
    if len(set(alleles)) > 1:
        return 1
    return 2


def _base_at_position(aln: pysam.AlignedSegment, target_pos: int) -> str:
    """
    Returnează baza unui read la o poziție genomică dată, parcurgând CIGAR.
    """
    if aln.cigartuples is None or aln.query_sequence is None:
        return 'N'

    ref_p = aln.reference_start
    q_p   = 0
    seq   = aln.query_sequence

    for op, ln in aln.cigartuples:
        if op == 0 or op == 7 or op == 8:  # M / = / X
            if ref_p <= target_pos < ref_p + ln:
                offset = target_pos - ref_p
                if q_p + offset < len(seq):
                    return seq[q_p + offset].upper()
                return 'N'
            ref_p += ln
            q_p   += ln
        elif op == 1 or op == 4:           # I, S
            q_p += ln
        elif op == 2 or op == 3:           # D, N
            if ref_p <= target_pos < ref_p + ln:
                return '-'
            ref_p += ln

    return 'N'


def _encode_pileup_2d(bam: pysam.AlignmentFile,
                      chrom: str,
                      pos: int,
                      ref_base: str = 'N',
                      height: int = PILEUP_HEIGHT,
                      width:  int = PILEUP_WIDTH) -> np.ndarray:
    half  = width // 2
    start = max(0, pos - half)
    end   = start + width

    img = np.zeros((N_CHANNELS, height, width), dtype=np.float32)

    ref_upper    = ref_base.upper() if ref_base else 'N'
    use_majority = ref_upper in ('N', 'MAJORITY', '')

    # Normalizare contig
    bam_refs = set(bam.references)
    if chrom not in bam_refs:
        alt_c = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        if alt_c in bam_refs:
            chrom = alt_c
        else:
            return img

    try:
        # =================================================================
        # PASUL 1: Colectăm read-urile care trec PRIN POZIȚIA CENTRALĂ
        # =================================================================
        center_reads: List[pysam.AlignedSegment] = []
        for aln in bam.fetch(chrom, pos, pos + 1):
            # Filtrăm: exclude duplicates, unmapped, secondary, supplementary
            if aln.is_duplicate or aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            if aln.mapping_quality < 20:
                continue
            if aln.query_sequence is None:
                continue
            center_reads.append(aln)

        if not center_reads:
            return img

        # =================================================================
        # PASUL 2: Calculăm baza la centru pentru fiecare read
        # =================================================================
        read_center_base = {
            id(aln): _base_at_position(aln, pos)
            for aln in center_reads
        }

        # =================================================================
        # PASUL 3: Determinăm ref_base dacă e MAJORITY
        # =================================================================
        if use_majority:
            counts = Counter(b for b in read_center_base.values() if b in 'ACGT')
            ref_upper = counts.most_common(1)[0][0] if counts else 'N'

        # =================================================================
        # PASUL 4: Sortăm read-urile — ALT primele, REF apoi, balansat
        # =================================================================
        alt_reads, ref_reads, oth_reads = [], [], []
        for aln in center_reads:
            b = read_center_base.get(id(aln), 'N')
            if b in 'ACGT' and b != ref_upper:
                alt_reads.append(aln)
            elif b == ref_upper:
                ref_reads.append(aln)
            else:
                oth_reads.append(aln)

        # Balansare 50/50 când avem prea multe read-uri
        if len(alt_reads) + len(ref_reads) > height:
            target_alt = min(len(alt_reads), height // 2)
            target_ref = min(len(ref_reads), height - target_alt)
            alt_reads = alt_reads[:target_alt]
            ref_reads = ref_reads[:target_ref]

        ordered = (alt_reads + ref_reads + oth_reads)[:height]

        # Dacă depth < height, repetăm ciclic (NU lăsăm rânduri zero)
        if 0 < len(ordered) < height:
            base_list = list(ordered)
            while len(ordered) < height:
                need = height - len(ordered)
                ordered.extend(base_list[:min(need, len(base_list))])
            ordered = ordered[:height]

        # =================================================================
        # PASUL 5: Umplem imaginea, parcurgând CIGAR pentru fiecare read
        # =================================================================
        for row_idx, aln in enumerate(ordered):
            if row_idx >= height:
                break

            seq    = aln.query_sequence
            quals  = aln.query_qualities
            mq     = aln.mapping_quality or 0
            is_fwd = not aln.is_reverse
            mq_n   = min(mq, MAX_MAP_QUAL) / MAX_MAP_QUAL

            if aln.cigartuples is None:
                continue

            ref_p = aln.reference_start
            q_p   = 0

            for op, ln in aln.cigartuples:
                if op == 0 or op == 7 or op == 8:  # M, =, X
                    for k in range(ln):
                        col_idx = ref_p + k - start
                        if 0 <= col_idx < width and q_p + k < len(seq):
                            base = seq[q_p + k].upper()
                            img[0, row_idx, col_idx] = BASE_ENC.get(base, 0.0)
                            img[1, row_idx, col_idx] = float(
                                base != ref_upper and base in 'ACGT'
                            )
                            bq = quals[q_p + k] if quals is not None else 20
                            img[2, row_idx, col_idx] = min(bq, MAX_BASE_QUAL) / MAX_BASE_QUAL
                            img[3, row_idx, col_idx] = mq_n
                            img[4, row_idx, col_idx] = float(is_fwd)
                    ref_p += ln
                    q_p   += ln
                elif op == 1 or op == 4:
                    q_p += ln
                elif op == 2 or op == 3:
                    for k in range(ln):
                        col_idx = ref_p + k - start
                        if 0 <= col_idx < width:
                            img[5, row_idx, col_idx] = 1.0
                            img[3, row_idx, col_idx] = mq_n
                            img[4, row_idx, col_idx] = float(is_fwd)
                    ref_p += ln

    except (ValueError, KeyError, AssertionError):
        pass

    return img


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------

class GenomicDataset(Dataset):
    def __init__(self,
                 bam_path:    str,
                 vcf_path:    str,
                 bed_path:    Optional[str] = None,
                 max_samples: Optional[int] = None,
                 seed:        int = 42,
                 img_height:  int = PILEUP_HEIGHT):
        self.bam_path    = bam_path
        self.vcf_path    = vcf_path
        self.bed_path    = bed_path
        self.max_samples = max_samples
        self.seed        = seed
        self.img_height  = img_height
        self.data_points: List[Dict] = []
        self._build_index()

    def _cache_path(self) -> str:
        # "v3" invalidează cache-ul vechi (encoder schimbat)
        key = f"{self.vcf_path}|{self.bed_path}|{self.max_samples}|{self.seed}|v3"
        key_hash = hashlib.md5(key.encode()).hexdigest()[:10]
        vcf_stem = os.path.splitext(os.path.basename(self.vcf_path))[0]
        cache_dir = os.path.join(os.path.dirname(self.vcf_path), ".cache")
        os.makedirs(cache_dir, exist_ok=True)
        return os.path.join(cache_dir, f"{vcf_stem}_{key_hash}.pkl")

    def _build_index(self):
        cache_path = self._cache_path()

        if os.path.exists(cache_path):
            print(f"   ⚡ Cache găsit! Încărcare rapidă...")
            with open(cache_path, "rb") as f:
                self.data_points = pickle.load(f)
            counts = Counter(p["label"] for p in self.data_points)
            print(f"   ✅ {len(self.data_points)} exemple din cache | "
                  f"Ref={counts[0]} | Het={counts[1]} | Hom={counts[2]}")
            return

        print(f"   🔍 Prima rulare: parsăm VCF...")

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

            pos = variant.POS - 1
            if len(variant.ALT) > 1:
                continue

            gt = variant.genotypes[0][:2]
            label = _genotype_to_label(tuple(gt))
            if label == -1:
                continue

            points.append({
                "chrom": chrom, "pos": pos,
                "ref":   variant.REF,
                "alt":   variant.ALT[0] if variant.ALT else ".",
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
            for offset in [25, -25, 50, -50, 75, -75]:
                new_pos = dp["pos"] + offset
                if new_pos < 0:
                    continue
                if (dp["chrom"], new_pos) not in variant_pos:
                    ref_points.append({
                        "chrom": dp["chrom"], "pos": new_pos,
                        "ref":   "MAJORITY", "alt": ".",
                        "label": 0,
                    })
                    variant_pos.add((dp["chrom"], new_pos))
                    break

        points.extend(ref_points)

        rng = random.Random(self.seed)
        rng.shuffle(points)
        if self.max_samples is not None:
            points = points[:self.max_samples]
        self.data_points = points

        print(f"   💾 Salvăm cache...")
        with open(cache_path, "wb") as f:
            pickle.dump(self.data_points, f)

        counts = Counter(p["label"] for p in self.data_points)
        print(f"   ✅ {len(self.data_points)} exemple | "
              f"Ref={counts[0]} | Het={counts[1]} | Hom={counts[2]}")

    def __len__(self):
        return len(self.data_points)

    def __getitem__(self, idx):
        dp  = self.data_points[idx]
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        img = _encode_pileup_2d(
            bam,
            chrom    = dp["chrom"],
            pos      = dp["pos"],
            ref_base = dp.get("ref", "N"),
            height   = self.img_height,
            width    = PILEUP_WIDTH,
        )
        bam.close()
        return (torch.from_numpy(img),
                torch.tensor(dp["label"], dtype=torch.long))