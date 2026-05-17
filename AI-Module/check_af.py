"""
check_af.py — Verifică AF din BAM vs genotip din VCF pentru primele 200 SNV-uri.
"""
import pysam, cyvcf2, numpy as np
from collections import defaultdict

BAM = "data/HG002_Son/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"
VCF = "data/HG002_Son/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf"

bam      = pysam.AlignmentFile(BAM, "rb")
bam_refs = set(bam.references)
has_chr  = any(r.startswith("chr") for r in bam_refs)
vcf      = cyvcf2.VCF(VCF)

af_by_label   = defaultdict(list)
CLASS = {0:"Ref", 1:"Het", 2:"Hom-Alt"}
checked = zero_depth = 0

for v in vcf:
    if checked >= 200: break
    if len(v.ALT) != 1 or len(v.REF) != 1 or len(v.ALT[0]) != 1:
        continue  # doar SNV simplu

    gt = v.genotypes[0][:2]
    al = [a for a in gt if a is not None and a >= 0]
    if not al: continue
    if all(a==0 for a in al):        label=0
    elif len(set(al))>1:             label=1
    else:                            label=2

    chrom = v.CHROM
    if has_chr and not chrom.startswith("chr"):  chrom = "chr"+chrom
    elif not has_chr and chrom.startswith("chr"): chrom = chrom[3:]
    pos, ref, alt = v.POS-1, v.REF.upper(), v.ALT[0].upper()

    ref_c=alt_c=total=0
    try:
        for col in bam.pileup(chrom, pos, pos+1, truncate=True,
                              stepper="all", min_base_quality=0):
            if col.reference_pos != pos: continue
            for r in col.pileups:
                if r.is_del or r.is_refskip: continue
                try:
                    b = r.alignment.query_sequence[r.query_position].upper()
                    total += 1
                    if b==ref: ref_c+=1
                    elif b==alt: alt_c+=1
                except: pass
    except: pass

    if total==0: zero_depth+=1; continue

    af_by_label[label].append(alt_c/total)
    checked+=1

bam.close(); vcf.close()

print(f"Verificate: {checked} | Depth=0: {zero_depth}\n")
print(f"{'Clasa':<10} {'N':>4} {'AF_mean':>8} {'AF<0.1':>7} {'0.1-0.7':>8} {'AF>0.7':>8}")
print("-"*50)
for label in [0,1,2]:
    afs = np.array(af_by_label[label])
    if len(afs)==0: print(f"{CLASS[label]:<10}    0"); continue
    low = (afs<0.1).sum()
    mid = ((afs>=0.1)&(afs<=0.7)).sum()
    hi  = (afs>0.7).sum()
    print(f"{CLASS[label]:<10} {len(afs):>4} {afs.mean():>8.3f} {low:>7} {mid:>8} {hi:>8}")

print("\nPrimele 10 Het cu AF:")
for af in af_by_label[1][:10]:
    print(f"  AF={af:.3f}")
print("\nPrimele 10 Hom-Alt cu AF:")
for af in af_by_label[2][:10]:
    print(f"  AF={af:.3f}")