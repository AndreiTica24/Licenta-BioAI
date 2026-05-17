"""
check_af2.py — Diagnostic complet: tipuri variante + pileup cu filtre corecte.
"""
import pysam, cyvcf2, numpy as np
from collections import defaultdict, Counter

BAM = "data/HG002_Son/HG002.hiseq4000.wes-agilent.50x.dedup.grch38.bam"
VCF = "data/HG002_Son/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf"

bam      = pysam.AlignmentFile(BAM, "rb")
bam_refs = set(bam.references)
has_chr  = any(r.startswith("chr") for r in bam_refs)
vcf      = cyvcf2.VCF(VCF)

CLASS = {0:"Ref", 1:"Het", 2:"Hom-Alt"}

# ── PARTE 1: distribuția tipurilor de variante din VCF ────────────────────
print("="*60)
print("PARTE 1: Tipuri de variante în primele 500 din VCF")
print("="*60)

type_counts = Counter()
label_counts = Counter()
n = 0
for v in vcf:
    if n >= 500: break
    if len(v.ALT) != 1: continue
    gt = v.genotypes[0][:2]
    al = [a for a in gt if a is not None and a >= 0]
    if not al: continue
    if all(a==0 for a in al):    label=0
    elif len(set(al))>1:         label=1
    else:                        label=2
    label_counts[CLASS[label]] += 1

    ref_len = len(v.REF)
    alt_len = len(v.ALT[0])
    if ref_len == 1 and alt_len == 1:
        vtype = "SNV"
    elif ref_len > alt_len:
        vtype = "DEL"
    elif ref_len < alt_len:
        vtype = "INS"
    else:
        vtype = "MNV"
    type_counts[vtype] += 1
    n += 1

print(f"Tipuri: {dict(type_counts)}")
print(f"Labels: {dict(label_counts)}")
vcf.close()

# ── PARTE 2: AF cu filtre de calitate corecte (MQ≥20, BQ≥13) ────────────
print("\n" + "="*60)
print("PARTE 2: AF recalculat cu filtre corecte (MQ≥20, BQ≥13)")
print("="*60)

vcf = cyvcf2.VCF(VCF)
af_by_label = defaultdict(list)
checked = zero_depth = skip_type = 0

for v in vcf:
    if checked >= 300: break
    if len(v.ALT) != 1: continue
    if len(v.REF) != 1 or len(v.ALT[0]) != 1:
        skip_type += 1
        continue  # doar SNV

    gt = v.genotypes[0][:2]
    al = [a for a in gt if a is not None and a >= 0]
    if not al: continue
    if all(a==0 for a in al):    label=0
    elif len(set(al))>1:         label=1
    else:                        label=2

    chrom = v.CHROM
    if has_chr and not chrom.startswith("chr"):  chrom="chr"+chrom
    elif not has_chr and chrom.startswith("chr"): chrom=chrom[3:]
    pos, ref, alt = v.POS-1, v.REF.upper(), v.ALT[0].upper()

    ref_c=alt_c=other_c=total=0
    try:
        for col in bam.pileup(chrom, pos, pos+1, truncate=True,
                              stepper="samtools",        # mai strict
                              min_base_quality=13,       # filtru BQ
                              min_mapping_quality=20,    # filtru MQ
                              ignore_overlaps=True,      # evită dubla numărare
                              flag_filter=3848):         # exclude dup, unmapped, etc.
            if col.reference_pos != pos: continue
            for r in col.pileups:
                if r.is_del or r.is_refskip: continue
                try:
                    b = r.alignment.query_sequence[r.query_position].upper()
                    total += 1
                    if b == ref:   ref_c += 1
                    elif b == alt: alt_c += 1
                    else:          other_c += 1
                except: pass
    except: pass

    if total == 0:
        zero_depth += 1
        continue

    af = alt_c / total
    af_by_label[label].append(af)
    checked += 1

bam.close(); vcf.close()

print(f"SNV verificate: {checked} | Depth=0: {zero_depth} | Skip non-SNV: {skip_type}\n")
print(f"{'Clasa':<10} {'N':>4} {'AF_mean':>8} {'AF_std':>8} {'AF<0.1':>7} {'0.1-0.7':>8} {'AF>0.7':>8}")
print("-"*55)
for label in [0,1,2]:
    afs = np.array(af_by_label[label])
    if len(afs)==0: print(f"{CLASS[label]:<10}    0"); continue
    low = (afs<0.1).sum()
    mid = ((afs>=0.1)&(afs<=0.7)).sum()
    hi  = (afs>0.7).sum()
    print(f"{CLASS[label]:<10} {len(afs):>4} {afs.mean():>8.3f} {afs.std():>8.3f} {low:>7} {mid:>8} {hi:>8}")

print("\nDistribuție AF Het (toate valorile):")
afs_het = np.array(af_by_label[1])
if len(afs_het) > 0:
    hist, edges = np.histogram(afs_het, bins=[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.01])
    for i, c in enumerate(hist):
        print(f"  {edges[i]:.1f}-{edges[i+1]:.1f}: {'█'*c} {c}")

print("\nDistribuție AF Hom-Alt (toate valorile):")
afs_hom = np.array(af_by_label[2])
if len(afs_hom) > 0:
    hist, edges = np.histogram(afs_hom, bins=[0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.01])
    for i, c in enumerate(hist):
        print(f"  {edges[i]:.1f}-{edges[i+1]:.1f}: {'█'*c} {c}")

# ── PARTE 3: verificare rapidă dacă cromozomii se potrivesc ──────────────
print("\n" + "="*60)
print("PARTE 3: Primele 5 poziții cu depth>0 — verificare manuală")
print("="*60)
bam = pysam.AlignmentFile(BAM, "rb")
vcf = cyvcf2.VCF(VCF)
shown = 0
for v in vcf:
    if shown >= 5: break
    if len(v.ALT)!=1 or len(v.REF)!=1 or len(v.ALT[0])!=1: continue
    chrom = v.CHROM
    if has_chr and not chrom.startswith("chr"):  chrom="chr"+chrom
    elif not has_chr and chrom.startswith("chr"): chrom=chrom[3:]
    pos, ref, alt = v.POS-1, v.REF.upper(), v.ALT[0].upper()
    gt = v.genotypes[0][:2]
    al = [a for a in gt if a is not None and a >= 0]
    if not al: continue
    if all(a==0 for a in al):    label="Ref"
    elif len(set(al))>1:         label="Het"
    else:                        label="Hom"

    bases = Counter()
    total = 0
    try:
        for col in bam.pileup(chrom, pos, pos+1, truncate=True,
                              stepper="samtools", min_base_quality=13,
                              min_mapping_quality=20, flag_filter=3848):
            if col.reference_pos != pos: continue
            for r in col.pileups:
                if r.is_del or r.is_refskip: continue
                try:
                    b = r.alignment.query_sequence[r.query_position].upper()
                    bases[b] += 1
                    total += 1
                except: pass
    except: pass

    if total == 0: continue
    af = bases.get(alt, 0) / total
    print(f"  {chrom}:{v.POS}  REF={ref} ALT={alt}  GT={label}")
    print(f"    Bases în BAM: {dict(bases)}  total={total}")
    print(f"    AF={af:.3f}  {'✅ OK' if (label=='Het' and 0.2<af<0.8) or (label=='Hom' and af>0.8) else '❌ Greșit'}")
    shown += 1

bam.close(); vcf.close()