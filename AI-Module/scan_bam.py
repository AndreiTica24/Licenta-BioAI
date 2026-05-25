"""
scan_bam.py — Pre-filtrare BAM pentru variant calling
============================================================
Scanează un fișier BAM exomic și identifică poziții candidate pentru
clasificare ulterioară de către modelul CNN 1D.

ARGUMENT TEHNIC:
Un BAM exomic conține ~60 milioane de poziții cu coverage. Din acestea,
>99.9% sunt poziții Ref (homozigot referință). Procesarea fiecărei poziții
de către modelul AI ar fi computațional ineficientă și ar produce milioane
de false positives.

Acest script aplică pragul STANDARD din literatura de specialitate:
  - depth >= 10 (confidence statistică pentru genotyping)
  - AF >= 0.15 (filtru pentru zgomot de secvențiere ~1-3%)

Referințe:
  - Poplin et al. (2018) "A universal SNP and small-indel variant caller
    using deep neural networks." Nature Biotechnology
  - Kim et al. (2018) "Strelka2: fast and accurate calling of germline
    and somatic variants." Nature Methods

OUTPUT:
  Fișier TSV cu coloanele: chrom, pos, ref_base, depth, alt_count, AF
  Acest fișier va fi consumat de predict.py pentru clasificarea finală.

Rulare:
    python scan_bam.py --bam input.bam --output candidates.tsv
    python scan_bam.py --bam input.bam --output candidates.tsv --threads 8
"""

import argparse
import logging
import multiprocessing as mp
import os
import time
from collections import Counter
from pathlib import Path
from typing import List, Tuple, Dict

import pysam

# ============================================================================
# Configurare logging
# ============================================================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


# ============================================================================
# Constante pentru pre-filtrare (PRAGUL STANDARD)
# ============================================================================
MIN_DEPTH        = 10      # Minimum 10 read-uri (statistic robust)
MIN_AF           = 0.15    # Minimum 15% reads cu baza alternativă (anti-zgomot)
MIN_BASE_QUAL    = 13      # Filtru calitate bază (Q13 = 95% acuratețe)
MIN_MAPPING_QUAL = 20      # Filtru calitate mapping (Q20 = 99% acuratețe)

# Bază reală pentru DNA — alte caractere (N, etc.) sunt ignorate
VALID_BASES = {'A', 'C', 'G', 'T'}

# Cromozomii umani standard (1-22 + X, Y)
HUMAN_CHROMOSOMES = (
    [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"] +
    [str(i) for i in range(1, 23)] + ["X", "Y"]
)


# ============================================================================
# Pre-filtrare la nivel de cromozom (rulează în paralel)
# ============================================================================

def scan_chromosome(args: Tuple[str, str, str]) -> Tuple[str, List[Dict], int]:
    """
    Scanează un singur cromozom și returnează lista de candidați.

    Argumente:
        bam_path : calea către BAM
        chrom    : numele cromozomului (ex: 'chr1' sau '1')

    Returneaza:
        chrom    : numele cromozomului
        candidates: lista de dict cu {chrom, pos, ref_base, depth, alt_count, AF}
        n_positions_scanned : numărul total de poziții scanate
    """
    bam_path, chrom, fasta_path = args

    candidates = []
    n_scanned = 0

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
        fasta = pysam.FastaFile(fasta_path) if fasta_path else None

        # Verificăm că cromozomul există în BAM
        bam_refs = set(bam.references)
        if chrom not in bam_refs:
            bam.close()
            if fasta:
                fasta.close()
            return chrom, [], 0

        # Pentru FASTA, normalizăm numele cromozomului (poate avea sau nu 'chr')
        fasta_chrom = chrom
        if fasta is not None:
            fasta_refs = set(fasta.references)
            if fasta_chrom not in fasta_refs:
                alt_c = (fasta_chrom[3:] if fasta_chrom.startswith("chr")
                         else "chr" + fasta_chrom)
                if alt_c in fasta_refs:
                    fasta_chrom = alt_c
                else:
                    fasta = None  # FASTA nu are acest cromozom

        # ───────────────────────────────────────────────────────────────────
        # Pileup pe întregul cromozom
        # ───────────────────────────────────────────────────────────────────
        for col in bam.pileup(chrom,
                              min_base_quality=MIN_BASE_QUAL,
                              min_mapping_quality=MIN_MAPPING_QUAL,
                              stepper="all",
                              ignore_overlaps=False):

            n_scanned += 1
            pos = col.reference_pos

            # Numărăm bazele observate
            base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            depth = 0

            for read in col.pileups:
                if read.is_del or read.is_refskip:
                    continue
                try:
                    if read.alignment.query_sequence is None:
                        continue
                    base = read.alignment.query_sequence[
                        read.query_position
                    ].upper()
                    if base in VALID_BASES:
                        base_counts[base] += 1
                        depth += 1
                except (IndexError, TypeError):
                    continue

            # FILTRU 1: depth minim
            if depth < MIN_DEPTH:
                continue

            # Determinăm baza de referință
            if fasta is not None:
                try:
                    ref_base = fasta.fetch(fasta_chrom, pos, pos + 1).upper()
                    if ref_base not in VALID_BASES:
                        continue
                except (ValueError, KeyError):
                    # Fallback: baza majoritară din BAM
                    ref_base = max(base_counts, key=base_counts.get)
            else:
                # Fallback: baza majoritară (presupunem că majoritatea reads = ref)
                ref_base = max(base_counts, key=base_counts.get)

            # Calculăm AF (procent reads cu baza alternativă)
            alt_count = depth - base_counts[ref_base]
            af = alt_count / depth

            # FILTRU 2: AF minim
            if af < MIN_AF:
                continue

            # Poziție candidat — o salvăm
            candidates.append({
                "chrom":     chrom,
                "pos":       pos,
                "ref_base":  ref_base,
                "depth":     depth,
                "alt_count": alt_count,
                "AF":        round(af, 4),
            })

        bam.close()
        if fasta:
            fasta.close()

    except Exception as e:
        logger.error(f"Eroare pe cromozomul {chrom}: {e}")

    return chrom, candidates, n_scanned


# ============================================================================
# Pipeline principal
# ============================================================================

def scan_bam(bam_path: str,
             output_path: str,
             fasta_path: str = None,
             threads: int = 4,
             chromosomes: List[str] = None) -> Dict:
    """
    Scanează un BAM complet și salvează candidații într-un TSV.

    Argumente:
        bam_path    : calea către BAM
        output_path : calea de output (TSV)
        fasta_path  : (opțional) FASTA pentru baza de referință
        threads     : numărul de procese paralele
        chromosomes : (opțional) listă explicită de cromozomi de scanat

    Returneaza:
        dict cu statistici (n_candidates, n_scanned, time_seconds)
    """
    bam_path    = str(bam_path)
    output_path = str(output_path)

    # Verificări preliminare
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"BAM nu există: {bam_path}")

    bai_path = bam_path + ".bai"
    if not os.path.exists(bai_path) and not os.path.exists(
        bam_path.replace(".bam", ".bai")
    ):
        raise FileNotFoundError(f"Index .bai lipsă pentru {bam_path}")

    if fasta_path and not os.path.exists(fasta_path):
        logger.warning(f"FASTA nu există: {fasta_path}, baza ref va fi inferată")
        fasta_path = None

    # Detectăm cromozomii din BAM
    bam = pysam.AlignmentFile(bam_path, "rb")
    all_refs = list(bam.references)
    bam.close()

    if chromosomes is None:
        # Selectăm doar cromozomii umani standard (1-22, X, Y)
        chromosomes = [c for c in all_refs if c in HUMAN_CHROMOSOMES]
        # Sortăm: 1-22, apoi X, Y
        def chr_sort_key(c):
            c_clean = c.replace("chr", "")
            if c_clean == "X":  return 23
            if c_clean == "Y":  return 24
            try:                return int(c_clean)
            except:             return 99
        chromosomes.sort(key=chr_sort_key)

    logger.info(f"BAM: {bam_path}")
    logger.info(f"Cromozomi de scanat: {len(chromosomes)} "
                f"({chromosomes[0]} ... {chromosomes[-1]})")
    logger.info(f"Procese paralele: {threads}")
    logger.info(f"Praguri filtre: depth>={MIN_DEPTH}, AF>={MIN_AF}, "
                f"BQ>={MIN_BASE_QUAL}, MQ>={MIN_MAPPING_QUAL}")

    # ───────────────────────────────────────────────────────────────────────
    # Scanare paralelă pe cromozomi
    # ───────────────────────────────────────────────────────────────────────
    print()
    print("=" * 70)
    print(f"🧬 Scanare BAM exomic: {os.path.basename(bam_path)}")
    print("=" * 70)
    print()

    t0 = time.time()
    all_candidates = []
    total_scanned = 0
    per_chrom_stats = []

    # Pregătim argumentele pentru procese
    args_list = [(bam_path, chrom, fasta_path) for chrom in chromosomes]

    if threads > 1:
        with mp.Pool(processes=threads) as pool:
            # imap_unordered = afișăm progres pe măsură ce termină procesele
            for i, (chrom, candidates, n_scanned) in enumerate(
                pool.imap_unordered(scan_chromosome, args_list), 1
            ):
                all_candidates.extend(candidates)
                total_scanned += n_scanned
                per_chrom_stats.append({
                    "chrom": chrom,
                    "scanned": n_scanned,
                    "candidates": len(candidates),
                })
                elapsed = time.time() - t0
                eta = elapsed / i * (len(chromosomes) - i)
                print(f"  [{i:2d}/{len(chromosomes)}] {chrom:8s}  "
                      f"scanat={n_scanned:>10,}  "
                      f"candidați={len(candidates):>7,}  "
                      f"({elapsed:.0f}s, ETA {eta:.0f}s)", flush=True)
    else:
        # Mod single-thread (pentru debug)
        for i, args in enumerate(args_list, 1):
            chrom, candidates, n_scanned = scan_chromosome(args)
            all_candidates.extend(candidates)
            total_scanned += n_scanned
            per_chrom_stats.append({
                "chrom": chrom,
                "scanned": n_scanned,
                "candidates": len(candidates),
            })
            elapsed = time.time() - t0
            print(f"  [{i:2d}/{len(chromosomes)}] {chrom:8s}  "
                  f"scanat={n_scanned:>10,}  "
                  f"candidați={len(candidates):>7,}  "
                  f"({elapsed:.0f}s)")

    elapsed = time.time() - t0

    # ───────────────────────────────────────────────────────────────────────
    # Salvăm rezultatele
    # ───────────────────────────────────────────────────────────────────────
    print()
    print("=" * 70)
    print(f"💾 Salvăm {len(all_candidates):,} candidați în {output_path}")
    print("=" * 70)

    # Sortăm candidații (chrom, pos) pentru consistență
    def sort_key(c):
        chrom = c["chrom"].replace("chr", "")
        if chrom == "X":  chrom_n = 23
        elif chrom == "Y": chrom_n = 24
        else:
            try: chrom_n = int(chrom)
            except: chrom_n = 99
        return (chrom_n, c["pos"])

    all_candidates.sort(key=sort_key)

    # Scriem TSV
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("chrom\tpos\tref_base\tdepth\talt_count\tAF\n")
        for c in all_candidates:
            f.write(f"{c['chrom']}\t{c['pos']}\t{c['ref_base']}\t"
                    f"{c['depth']}\t{c['alt_count']}\t{c['AF']}\n")

    # ───────────────────────────────────────────────────────────────────────
    # Statistici finale
    # ───────────────────────────────────────────────────────────────────────
    print()
    print("=" * 70)
    print("📊 STATISTICI FINALE")
    print("=" * 70)
    print(f"   Poziții scanate     : {total_scanned:>15,}")
    print(f"   Candidați găsiți    : {len(all_candidates):>15,}")
    if total_scanned > 0:
        ratio = len(all_candidates) / total_scanned * 100
        print(f"   Rată reținere       : {ratio:>14.3f}%")
    print(f"   Timp total          : {elapsed:>14.1f}s ({elapsed/60:.1f} min)")
    print(f"   Output              : {output_path}")
    print()

    # Distribuție AF
    if all_candidates:
        afs = [c["AF"] for c in all_candidates]
        bins = [
            (0.15, 0.30, "0.15-0.30"),
            (0.30, 0.45, "0.30-0.45"),
            (0.45, 0.60, "0.45-0.60 (~Het)"),
            (0.60, 0.85, "0.60-0.85"),
            (0.85, 1.01, "0.85-1.00 (~Hom-Alt)"),
        ]
        print("   Distribuție AF candidați:")
        for low, high, label in bins:
            n = sum(1 for af in afs if low <= af < high)
            pct = n / len(afs) * 100
            bar = "█" * int(pct / 2)
            print(f"      {label:24s} {n:>7,} ({pct:5.1f}%)  {bar}")
        print()

    return {
        "n_candidates":      len(all_candidates),
        "n_positions_scanned": total_scanned,
        "time_seconds":      round(elapsed, 2),
        "per_chromosome":    per_chrom_stats,
    }


# ============================================================================
# CLI
# ============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Pre-filtrare BAM pentru variant calling",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--bam",     type=str, required=True,
                   help="Calea către fișierul BAM (cu .bai)")
    p.add_argument("--output",  type=str, required=True,
                   help="Calea de output (TSV)")
    p.add_argument("--fasta",   type=str,
                   default="data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
                   help="Calea către FASTA de referință (default: GRCh38)")
    p.add_argument("--threads", type=int, default=4,
                   help="Numărul de procese paralele (default: 4)")
    p.add_argument("--chroms",  type=str, default=None,
                   help="Cromozomi de scanat separați prin virgulă "
                        "(default: 1-22, X, Y)")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()

    chromosomes = None
    if args.chroms:
        chromosomes = [c.strip() for c in args.chroms.split(",")]

    stats = scan_bam(
        bam_path    = args.bam,
        output_path = args.output,
        fasta_path  = args.fasta if os.path.exists(args.fasta) else None,
        threads     = args.threads,
        chromosomes = chromosomes,
    )

    print(f"✅ Pre-filtrare completă. Următorul pas: python predict.py "
          f"--candidates {args.output} --bam {args.bam}")