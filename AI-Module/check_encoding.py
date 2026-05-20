"""
check_encoding.py — Verifică ce vede modelul pentru fiecare clasă.
"""
import numpy as np
import pysam
import cyvcf2
from collections import defaultdict
from dataset import _encode_window, WINDOW_SIZE

BAM = "data/HG002_Son/HG002.hiseq4000.wes-agilent.50x.dedup.grch38.bam"
VCF = "data/HG002_Son/HG002_exome.vcf"

bam = pysam.AlignmentFile(BAM, "rb")
bam_has_chr = any(r.startswith("chr") for r in bam.references)
vcf = cyvcf2.VCF(VCF)

CLASS = {0: "Ref", 1: "Het", 2: "Hom-Alt"}
samples = defaultdict(list)
N = 20

for v in vcf:
    if all(len(s) >= N for s in samples.values()) and len(samples) >= 2:
        break
    if len(v.ALT) != 1 or len(v.REF) != 1 or len(v.ALT[0]) != 1:
        continue
    gt = v.genotypes[0][:2]
    al = [a for a in gt if a is not None and a >= 0]
    if not al: continue
    if all(a == 0 for a in al): label = 0
    elif len(set(al)) > 1:      label = 1
    else:                       label = 2

    if len(samples[label]) >= N:
        continue

    chrom = v.CHROM
    if bam_has_chr and not chrom.startswith("chr"): chrom = "chr" + chrom
    elif not bam_has_chr and chrom.startswith("chr"): chrom = chrom[3:]

    samples[label].append((chrom, v.POS - 1, v.REF, v.ALT[0]))

vcf.close()

print("="*70)
print("DIAGNOSTIC ENCODER CNN 1D")
print("="*70)
print(f"\nNumărul exemplelor per clasă: ", end="")
print({CLASS[k]: len(v) for k, v in samples.items()})

# Pentru fiecare clasă, calculăm statistici la coloana centrală
center = WINDOW_SIZE // 2

for label, examples in sorted(samples.items()):
    if not examples:
        continue
    
    onehot_center = []  # one-hot la centru (suma canalelor 0-3)
    depth_center  = []  # canal 4
    af_center     = []  # canal 5
    nonzero_cols  = []  # câte coloane au depth > 0
    
    for chrom, pos, ref, alt in examples:
        x = _encode_window(bam, chrom, pos)  # (6, 200)
        
        onehot_center.append(x[:4, center].sum())  # 1 dacă există o bază
        depth_center.append(x[4, center])
        af_center.append(x[5, center])
        nonzero_cols.append((x[4, :] > 0).sum())
    
    print(f"\n📊 {CLASS[label]} (N={len(examples)}):")
    print(f"   One-hot center  : mean={np.mean(onehot_center):.3f}  "
          f"(1.0 = baza setată, 0.0 = nimic)")
    print(f"   Depth center    : mean={np.mean(depth_center):.3f}  "
          f"(canal 4, normalizat 0-1)")
    print(f"   AF center       : mean={np.mean(af_center):.3f}  "
          f"std={np.std(af_center):.3f}")
    print(f"   Coloane non-zero: mean={np.mean(nonzero_cols):.1f}/200")
    print(f"   AF distribuție  : ", end="")
    afs = np.array(af_center)
    for low, high, label_b in [(0, 0.05, "0-5%"), (0.05, 0.3, "5-30%"),
                                (0.3, 0.7, "30-70%"), (0.7, 0.95, "70-95%"),
                                (0.95, 1.01, "95-100%")]:
        n = ((afs >= low) & (afs < high)).sum()
        print(f"{label_b}:{n} ", end="")
    print()

# Verificăm dacă imaginile diferă realist între clase
if 1 in samples and 2 in samples:
    chrom_h, pos_h, _, _ = samples[1][0]
    chrom_o, pos_o, _, _ = samples[2][0]
    x_het = _encode_window(bam, chrom_h, pos_h)
    x_hom = _encode_window(bam, chrom_o, pos_o)
    diff = np.abs(x_het - x_hom).mean()
    print(f"\n📐 Diferență medie absolută între Het și Hom: {diff:.4f}")
    print(f"   (trebuie să fie > 0.05 pentru ca modelul să poată distinge)")
    
    # Pe canalul AF specifically
    diff_af = np.abs(x_het[5, :] - x_hom[5, :]).mean()
    print(f"   Diferență pe canalul AF (canal 5): {diff_af:.4f}")

bam.close()
print("\n" + "="*70)