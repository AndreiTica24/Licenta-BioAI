"""
test_pileup.py — Test direct pysam la primele 10 Hom-Alt din VCF-ul exomic
"""
import pysam
import cyvcf2

BAM = "data/HG002_Son/HG002.hiseq4000.wes-agilent.50x.dedup.grch38.bam"
VCF = "data/HG002_Son/HG002_exome.vcf"

bam = pysam.AlignmentFile(BAM, "rb")
bam_has_chr = any(r.startswith("chr") for r in bam.references)
vcf = cyvcf2.VCF(VCF)

shown = 0
print(f"BAM contigs[:5]: {list(bam.references)[:5]}")
print(f"BAM has 'chr': {bam_has_chr}")
print()

for v in vcf:
    if shown >= 10:
        break
    if len(v.ALT) != 1 or len(v.REF) != 1 or len(v.ALT[0]) != 1:
        continue
    gt = v.genotypes[0][:2]
    al = [a for a in gt if a is not None and a >= 0]
    if not al: continue
    if all(a == 0 for a in al): label = "Ref"
    elif len(set(al)) > 1:      label = "Het"
    else:                       label = "Hom"

    if label != "Hom":
        continue

    chrom = v.CHROM
    if bam_has_chr and not chrom.startswith("chr"): chrom = "chr" + chrom
    elif not bam_has_chr and chrom.startswith("chr"): chrom = chrom[3:]

    # Metoda 1: count direct
    cnt = bam.count(chrom, v.POS - 1, v.POS)

    # Metoda 2: fetch read-uri și verificăm dacă acoperă pos
    reads = list(bam.fetch(chrom, v.POS - 1, v.POS))
    n_reads_cover = sum(1 for r in reads
                       if r.reference_start <= v.POS - 1 < r.reference_end)

    # Metoda 3: pileup cu filtre minime
    n_pileup = 0
    for col in bam.pileup(chrom, v.POS - 1, v.POS,
                          truncate=True, min_base_quality=0,
                          min_mapping_quality=0, stepper="nofilter"):
        if col.reference_pos == v.POS - 1:
            n_pileup = len(col.pileups)
            break

    print(f"{chrom}:{v.POS}  Hom-Alt  REF={v.REF} ALT={v.ALT[0]}  GT={gt}")
    print(f"   bam.count():        {cnt}")
    print(f"   reads in fetch:     {len(reads)} (overlap pos: {n_reads_cover})")
    print(f"   pileup col reads:   {n_pileup}")
    print()

    shown += 1

bam.close()
vcf.close()