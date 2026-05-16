"""
dataset.py — GenomicDataset
Transformă un BAM + VCF (+ BED opțional) într-un dataset PyTorch.

Fiecare exemplu = o imagine pileup RGB (6 canale) de dimensiune (6, H, W):
  Canal 0: frecvența alelei alternative (AF)
  Canal 1: acoperire (depth) normalizată
  Canal 2: calitate medie baze
  Canal 3: calitate medie mapping
  Canal 4: strand bias (forward ratio)
  Canal 5: deletions ratio

Label: 0=Ref, 1=Het, 2=Hom-Alt
"""

import os
import pickle
import hashlib
import random
import logging
from pathlib import Path
from typing import Optional, List, Dict, Tuple

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
PILEUP_WIDTH   = 100   # fereastră ± 50 bp în jurul variantei
PILEUP_HEIGHT  = 100   # număr de read-uri afișate (rânduri imagine)
N_CHANNELS     = 6     # canale per pixel
MAX_DEPTH      = 200   # depth de saturație
MAX_BASE_QUAL  = 40
MAX_MAP_QUAL   = 60

# Codificare nucleotide → valoare [0,1]
BASE_ENC = {'A': 0.25, 'C': 0.50, 'G': 0.75, 'T': 1.00,
            'N': 0.00, '-': 0.00, '*': 0.00}


# ---------------------------------------------------------------------------
# Funcții ajutătoare
# ---------------------------------------------------------------------------

def _load_bed_regions(bed_path: str) -> Dict[str, List[Tuple[int, int]]]:
    """Parsează un fișier BED și returnează {chrom: [(start, end), ...]}."""
    regions: Dict[str, List[Tuple[int, int]]] = {}
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            regions.setdefault(chrom, []).append((start, end))
    return regions


def _in_bed(chrom: str, pos: int,
            bed: Optional[Dict[str, List[Tuple[int, int]]]]) -> bool:
    """Verifică dacă (chrom, pos) se află într-o regiune BED."""
    if bed is None:
        return True
    for start, end in bed.get(chrom, []):
        if start <= pos < end:
            return True
    return False


def _genotype_to_label(gt: Tuple[int, ...]) -> int:
    """
    Convertește un genotype pysam/cyvcf2 la label:
      0/0 → 0 (Ref)
      0/1 → 1 (Het)
      1/1 → 2 (Hom-Alt)
    """
    alleles = [a for a in gt if a is not None and a >= 0]
    if not alleles:
        return -1  # date lipsă, va fi filtrat
    if all(a == 0 for a in alleles):
        return 0
    if len(set(alleles)) > 1:
        return 1
    return 2


# ---------------------------------------------------------------------------
# Encoding pileup — reprezentare read-level 2D (corectă)
# ---------------------------------------------------------------------------

def _get_read_base(read: pysam.PileupRead, col_pos: int) -> str:
    """Returnează baza unui read la o poziție de referință dată."""
    if read.is_del or read.is_refskip:
        return '-'
    try:
        seq = read.alignment.query_sequence
        if seq is None:
            return 'N'
        return seq[read.query_position].upper()
    except (IndexError, TypeError):
        return 'N'


def _encode_pileup_2d(bam: pysam.AlignmentFile,
                      chrom: str,
                      pos: int,
                      ref_base: str = 'N',
                      height: int = PILEUP_HEIGHT,
                      width:  int = PILEUP_WIDTH) -> np.ndarray:
    """
    Construiește o imagine pileup read-level de formă (N_CHANNELS, height, width).

    Fiecare rând = un read diferit; fiecare coloană = o poziție genomică.
    Dacă sunt mai puțin de `height` read-uri, rândurile rămase sunt zero.

    Canale per pixel (x=read, y=coloană):
      0: baza nucleotidică encodată  [0,1]
      1: este baza alternativă?      {0,1}  ← SEMNALUL CHEIE pentru AF
      2: calitate bază               [0,1]
      3: calitate mapping            [0,1]
      4: direcție strand             {0=rev, 1=fwd}
      5: este deletion?              {0,1}
    """
    half  = width // 2
    start = max(0, pos - half)
    end   = start + width

    # img[canal, read_idx, col_idx]
    img = np.zeros((N_CHANNELS, height, width), dtype=np.float32)

    ref_base_upper = ref_base.upper() if ref_base else 'N'

    # --- Normalizare contig: BAM poate folosi '1' în loc de 'chr1' ---
    bam_refs = set(bam.references)
    if chrom not in bam_refs:
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        if alt in bam_refs:
            chrom = alt

    try:
        # Colectăm toate pileup-urile o singură dată per coloană
        # Structura: col_reads[col_idx] = list of PileupRead
        col_reads: Dict[int, list] = {}

        for col in bam.pileup(chrom, start, end,
                               truncate=True,
                               min_base_quality=0,
                               min_mapping_quality=0,
                               stepper="all",
                               ignore_overlaps=False):
            col_pos = col.reference_pos
            if col_pos < start or col_pos >= end:
                continue
            idx = col_pos - start
            col_reads[idx] = list(col.pileups)

        # Construim o listă ordonată de read-uri la poziția centrală
        center_idx = pos - start
        center_idx = max(0, min(center_idx, width - 1))

        center_reads = col_reads.get(center_idx, [])

        # Sortăm: mai întâi read-urile cu baza alternativă (maximizăm
        # informația vizibilă în primele rânduri ale imaginii)
        def sort_key(r):
            base = _get_read_base(r, pos)
            is_alt = int(base != ref_base_upper and base not in ('N', '-'))
            return (-is_alt, r.alignment.query_name or "")

        center_reads_sorted = sorted(center_reads, key=sort_key)

        # Construim un index read_name → rând
        read_to_row: Dict[str, int] = {}
        for row_idx, pread in enumerate(center_reads_sorted[:height]):
            name = pread.alignment.query_name or str(row_idx)
            read_to_row[name] = row_idx

        # Umplem imaginea coloană cu coloană
        for col_idx, reads in col_reads.items():
            if col_idx < 0 or col_idx >= width:
                continue
            for pread in reads:
                name = pread.alignment.query_name or ""
                row_idx = read_to_row.get(name)
                if row_idx is None:
                    continue  # read-ul nu apare la poziția centrală

                base = _get_read_base(pread, start + col_idx)

                # Canal 0: nucleotidă
                img[0, row_idx, col_idx] = BASE_ENC.get(base, 0.0)

                # Canal 1: este baza alternativă? (semnalul cheie AF)
                is_alt = int(
                    base != ref_base_upper
                    and base not in ('N', '-', '*')
                    and not pread.is_del
                    and not pread.is_refskip
                )
                img[1, row_idx, col_idx] = float(is_alt)

                # Canal 2: calitate bază
                if not pread.is_del and not pread.is_refskip:
                    try:
                        quals = pread.alignment.query_qualities
                        bq = quals[pread.query_position] if quals is not None else 20
                    except (IndexError, TypeError):
                        bq = 20
                    img[2, row_idx, col_idx] = min(bq, MAX_BASE_QUAL) / MAX_BASE_QUAL

                # Canal 3: calitate mapping
                mq = pread.alignment.mapping_quality or 0
                img[3, row_idx, col_idx] = min(mq, MAX_MAP_QUAL) / MAX_MAP_QUAL

                # Canal 4: strand (1=forward, 0=reverse)
                img[4, row_idx, col_idx] = float(not pread.alignment.is_reverse)

                # Canal 5: deletion
                img[5, row_idx, col_idx] = float(pread.is_del)

    except (ValueError, KeyError, AssertionError):
        pass  # contig lipsă sau interval invalid — returnăm zeros

    return img  # (N_CHANNELS, height, width) = (6, 100, 100)


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------

class GenomicDataset(Dataset):
    """
    Dataset PyTorch pentru variant calling.

    Parametri
    ---------
    bam_path    : calea la fișierul BAM (indexat .bai)
    vcf_path    : calea la VCF-ul de referință (benchmark)
    bed_path    : (opțional) BED cu regiuni de înaltă-confidență
    max_samples : numărul maxim de exemple (None = toate)
    seed        : seed pentru reproducibilitate
    img_height  : înălțimea imaginii pileup 2D
    """

    def __init__(
        self,
        bam_path:    str,
        vcf_path:    str,
        bed_path:    Optional[str] = None,
        max_samples: Optional[int] = None,
        seed:        int = 42,
        img_height:  int = 100,
    ):
        self.bam_path   = bam_path
        self.vcf_path   = vcf_path
        self.bed_path   = bed_path
        self.max_samples = max_samples
        self.seed       = seed
        self.img_height = img_height

        self.data_points: List[Dict] = []
        self._build_index()

    # ------------------------------------------------------------------
    def _cache_path(self) -> str:
        """
        Generează o cale unică pentru fișierul cache bazată pe:
        vcf_path + bed_path + max_samples + seed
        → astfel, dacă schimbi oricare parametru, cache-ul vechi e ignorat automat.
        """
        key = f"{self.vcf_path}|{self.bed_path}|{self.max_samples}|{self.seed}"
        key_hash = hashlib.md5(key.encode()).hexdigest()[:10]
        vcf_stem  = os.path.splitext(os.path.basename(self.vcf_path))[0]
        cache_dir = os.path.join(os.path.dirname(self.vcf_path), ".cache")
        os.makedirs(cache_dir, exist_ok=True)
        return os.path.join(cache_dir, f"{vcf_stem}_{key_hash}.pkl")

    # ------------------------------------------------------------------
    def _build_index(self):
        """
        Parsează VCF-ul și construiește lista de poziții + labels.
        La prima rulare: parsează și salvează cache pe disc.
        La rulările ulterioare: încarcă direct din cache (instant).
        """
        from collections import Counter

        cache_path = self._cache_path()

        # ----------------------------------------------------------------
        # CACHE HIT — încărcare instantă
        # ----------------------------------------------------------------
        if os.path.exists(cache_path):
            logger.info(f"[GenomicDataset] 🚀 Cache găsit: {cache_path}")
            print(f"   ⚡ Cache găsit! Încărcare rapidă din: {os.path.basename(cache_path)}")
            with open(cache_path, "rb") as f:
                self.data_points = pickle.load(f)
            counts = Counter(p["label"] for p in self.data_points)
            print(
                f"   ✅ {len(self.data_points)} exemple din cache | "
                f"Ref={counts[0]} | Het={counts[1]} | Hom={counts[2]}"
            )
            return

        # ----------------------------------------------------------------
        # CACHE MISS — parsare completă (prima rulare)
        # ----------------------------------------------------------------
        logger.info(f"[GenomicDataset] Prima rulare — indexăm VCF: {self.vcf_path}")
        print(f"   🔍 Prima rulare: parsăm VCF (durează câteva minute)...")

        # Detectăm formatul contigurilor din BAM PRIMA dată — înainte de BED și VCF
        _bam_tmp = pysam.AlignmentFile(self.bam_path, "rb")
        bam_refs = set(_bam_tmp.references)
        bam_has_chr = any(r.startswith("chr") for r in bam_refs)
        _bam_tmp.close()

        bed = None
        if self.bed_path and Path(self.bed_path).exists():
            logger.info(f"[GenomicDataset] Încărcăm BED: {self.bed_path}")
            raw_bed = _load_bed_regions(self.bed_path)
            # Normalizăm cheile BED la același format ca BAM-ul
            bed = {}
            for c, intervals in raw_bed.items():
                if bam_has_chr and not c.startswith("chr"):
                    c = "chr" + c
                elif not bam_has_chr and c.startswith("chr"):
                    c = c[3:]
                bed[c] = intervals

        vcf = cyvcf2.VCF(self.vcf_path)
        points = []

        for variant in vcf:
            chrom = variant.CHROM

            # Normalizăm prefixul chr să se potrivească cu BAM-ul
            if bam_has_chr and not chrom.startswith("chr"):
                chrom = "chr" + chrom
            elif not bam_has_chr and chrom.startswith("chr"):
                chrom = chrom[3:]

            pos = variant.POS - 1  # 0-indexed

            # Filtru BED
            if not _in_bed(chrom, pos, bed):
                continue

            # Saltăm variantele multi-alelice complexe
            if len(variant.ALT) > 1:
                continue

            # Genotype din primul sample
            gt = variant.genotypes[0][:2]  # (allele1, allele2)
            label = _genotype_to_label(tuple(gt))
            if label == -1:
                continue

            points.append({
                "chrom": chrom,
                "pos":   pos,
                "ref":   variant.REF,
                "alt":   variant.ALT[0] if variant.ALT else ".",
                "label": label,
            })

        vcf.close()

        # Adăugăm exemple Ref (homozigot referință) din BED
        if bed is not None:
            variant_positions = {(p["chrom"], p["pos"]) for p in points}
            rng = random.Random(self.seed)

            ref_budget = len(points) // 2
            ref_added  = 0

            for bed_chrom, intervals in bed.items():
                # Normalizăm și cromozomul din BED
                chrom = bed_chrom
                if bam_has_chr and not chrom.startswith("chr"):
                    chrom = "chr" + chrom
                elif not bam_has_chr and chrom.startswith("chr"):
                    chrom = chrom[3:]

                if ref_added >= ref_budget:
                    break
                for start, end in intervals:
                    if ref_added >= ref_budget:
                        break
                    sample_n = min(5, end - start)
                    for pos in rng.sample(range(start, end), k=sample_n):
                        if (chrom, pos) not in variant_positions:
                            points.append({
                                "chrom": chrom,
                                "pos":   pos,
                                "ref":   "N",
                                "alt":   ".",
                                "label": 0,
                            })
                            ref_added += 1
                            if ref_added >= ref_budget:
                                break

        # Shuffle + limitare
        rng = random.Random(self.seed)
        rng.shuffle(points)

        if self.max_samples is not None:
            points = points[:self.max_samples]

        self.data_points = points

        # ----------------------------------------------------------------
        # Salvăm cache pe disc pentru rulările viitoare
        # ----------------------------------------------------------------
        logger.info(f"[GenomicDataset] 💾 Salvăm cache: {cache_path}")
        print(f"   💾 Salvăm cache pentru viitor: {os.path.basename(cache_path)}")
        with open(cache_path, "wb") as f:
            pickle.dump(self.data_points, f)

        # Statistici
        counts = Counter(p["label"] for p in self.data_points)
        logger.info(
            f"[GenomicDataset] Total: {len(self.data_points)} exemple | "
            f"Ref={counts[0]} | Het={counts[1]} | Hom-Alt={counts[2]}"
        )
        print(
            f"   ✅ {len(self.data_points)} exemple încărcate | "
            f"Ref={counts[0]} | Het={counts[1]} | Hom={counts[2]}"
        )

    # ------------------------------------------------------------------
    def __len__(self) -> int:
        return len(self.data_points)

    # ------------------------------------------------------------------
    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, torch.Tensor]:
        dp = self.data_points[idx]

        # Deschidem BAM-ul per worker (thread-safe — fiecare worker
        # are propria instanță AlignmentFile)
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        img = _encode_pileup_2d(
            bam,
            chrom    = dp["chrom"],
            pos      = dp["pos"],
            ref_base = dp.get("ref", "N"),   # ← transmitem baza de referință
            height   = self.img_height,
            width    = PILEUP_WIDTH,
        )
        bam.close()

        tensor = torch.from_numpy(img)                       # (6, H, W)
        label  = torch.tensor(dp["label"], dtype=torch.long)
        return tensor, label