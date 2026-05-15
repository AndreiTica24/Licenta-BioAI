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
PILEUP_HEIGHT  = 6    # număr canale / features per poziție
PILEUP_WIDTH   = 100  # fereastră ± 50 bp în jurul variantei
MAX_DEPTH      = 200  # depth de saturație pentru normalizare
MAX_BASE_QUAL  = 40
MAX_MAP_QUAL   = 60


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
# Encoding pileup
# ---------------------------------------------------------------------------

def _encode_pileup(bam: pysam.AlignmentFile,
                   chrom: str, pos: int) -> np.ndarray:
    """
    Extrage o fereastră de PILEUP_WIDTH bp centrată pe `pos`
    și returnează un tensor (PILEUP_HEIGHT, PILEUP_WIDTH).
    """
    half  = PILEUP_WIDTH // 2
    start = max(0, pos - half)
    end   = pos + half

    # Acumulatoare per coloană
    depth     = np.zeros(PILEUP_WIDTH, dtype=np.float32)
    alt_count = np.zeros(PILEUP_WIDTH, dtype=np.float32)
    base_qual = np.zeros(PILEUP_WIDTH, dtype=np.float32)
    map_qual  = np.zeros(PILEUP_WIDTH, dtype=np.float32)
    fwd_count = np.zeros(PILEUP_WIDTH, dtype=np.float32)
    del_count = np.zeros(PILEUP_WIDTH, dtype=np.float32)

    try:
        for col in bam.pileup(chrom, start, end,
                               truncate=True,
                               min_base_quality=0,
                               stepper="all"):
            col_pos = col.reference_pos
            if col_pos < start or col_pos >= end:
                continue
            idx = col_pos - start
            if idx < 0 or idx >= PILEUP_WIDTH:
                continue

            reads = list(col.pileups)
            if not reads:
                continue

            d = len(reads)
            depth[idx] = d

            bq_sum = mq_sum = fwd = dels = alt = 0
            ref_base = None  # nu avem referință, folosim baza majoritară

            for r in reads:
                if r.is_del or r.is_refskip:
                    dels += 1
                    continue
                bq_sum += r.alignment.query_qualities[r.query_position] \
                    if r.alignment.query_qualities is not None else 20
                mq_sum += r.alignment.mapping_quality
                if not r.alignment.is_reverse:
                    fwd += 1

            base_qual[idx] = bq_sum / max(d, 1)
            map_qual[idx]  = mq_sum / max(d, 1)
            fwd_count[idx] = fwd / max(d, 1)
            del_count[idx] = dels / max(d, 1)

    except (ValueError, KeyError):
        pass  # contig lipsă sau interval invalid

    # Calculăm AF la poziția centrală (pos)
    center = pos - start
    if 0 <= center < PILEUP_WIDTH and depth[center] > 0:
        # AF estimat din counts (simplificat: nu avem referința exactă)
        alt_count[center] = 1.0 - fwd_count[center]  # placeholder

    # Normalizare
    img = np.stack([
        alt_count,                                     # canal 0: AF
        np.clip(depth, 0, MAX_DEPTH) / MAX_DEPTH,     # canal 1: depth
        np.clip(base_qual, 0, MAX_BASE_QUAL) / MAX_BASE_QUAL,  # canal 2
        np.clip(map_qual,  0, MAX_MAP_QUAL)  / MAX_MAP_QUAL,   # canal 3
        fwd_count,                                     # canal 4: strand bias
        del_count,                                     # canal 5: del ratio
    ], axis=0).astype(np.float32)

    return img  # shape: (6, PILEUP_WIDTH)


def _encode_pileup_2d(bam: pysam.AlignmentFile,
                      chrom: str, pos: int,
                      height: int = 100) -> np.ndarray:
    """
    Versiune 2D: (6, height, PILEUP_WIDTH) — fiecare canal e replicat
    pe dimensiunea height pentru compatibilitate CNN 2D.
    Poți înlocui cu read-level pileup pentru mai multă informație.
    """
    img_1d = _encode_pileup(bam, chrom, pos)           # (6, W)
    img_2d = np.repeat(img_1d[:, np.newaxis, :], height, axis=1)  # (6, H, W)
    return img_2d


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
    def _build_index(self):
        """Parsează VCF-ul și construiește lista de poziții + labels."""
        logger.info(f"[GenomicDataset] Indexăm VCF: {self.vcf_path}")

        bed = None
        if self.bed_path and Path(self.bed_path).exists():
            logger.info(f"[GenomicDataset] Încărcăm BED: {self.bed_path}")
            bed = _load_bed_regions(self.bed_path)

        vcf = cyvcf2.VCF(self.vcf_path)
        points = []

        for variant in vcf:
            chrom = variant.CHROM
            # Normalizare prefix chr
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom

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

        # Adăugăm exemple Ref (homozigot referință) dacă avem BED
        # (simplu: luăm poziții aleatorii din BED care nu sunt în VCF)
        if bed is not None:
            variant_positions = {(p["chrom"], p["pos"]) for p in points}
            rng = random.Random(self.seed)

            ref_budget = len(points) // 2  # echilibrare parțială
            ref_added  = 0

            for chrom, intervals in bed.items():
                if ref_added >= ref_budget:
                    break
                for start, end in intervals:
                    if ref_added >= ref_budget:
                        break
                    # Câteva poziții aleatorii per interval
                    sample_n = min(5, end - start)
                    for pos in rng.sample(range(start, end), k=sample_n):
                        if (chrom, pos) not in variant_positions:
                            points.append({
                                "chrom": chrom,
                                "pos":   pos,
                                "ref":   "N",
                                "alt":   ".",
                                "label": 0,  # Ref
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

        # Statistici
        from collections import Counter
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

        # Deschidem BAM-ul per worker (thread-safe)
        bam = pysam.AlignmentFile(self.bam_path, "rb")
        img = _encode_pileup_2d(
            bam, dp["chrom"], dp["pos"], height=self.img_height
        )
        bam.close()

        tensor = torch.from_numpy(img)          # (6, H, W)
        label  = torch.tensor(dp["label"], dtype=torch.long)
        return tensor, label