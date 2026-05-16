"""
debug_pipeline.py — Diagnosticare rapidă a problemei de training collapse.
Rulează ÎNAINTE de train.py pentru a verifica fiecare componentă.

    python debug_pipeline.py
"""

import sys
import numpy as np
import torch
import pysam
import cyvcf2

# ---------------------------------------------------------------------------
# CONFIG — modifică aceste căi dacă diferă la tine
# ---------------------------------------------------------------------------
BAM_PATH = "data/HG002_Son/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"
VCF_PATH = "data/HG002_Son/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf"
BED_PATH = "data/HG002_Son/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

PASS = "✅"
FAIL = "❌"
WARN = "⚠️ "

errors = []

# ===========================================================================
# TEST 1 — BAM accesibil și indexat
# ===========================================================================
print("\n" + "="*60)
print("TEST 1: BAM accesibil și indexat")
print("="*60)
try:
    bam = pysam.AlignmentFile(BAM_PATH, "rb")
    n_refs = bam.nreferences
    print(f"{PASS} BAM deschis OK — {n_refs} referințe")

    # Primul contig disponibil
    first_ref = bam.references[0] if n_refs > 0 else None
    print(f"   Primul contig: {first_ref}")

    # Verificăm dacă există indexul .bai
    import os
    bai = BAM_PATH + ".bai"
    bai2 = BAM_PATH.replace(".bam", ".bai")
    if os.path.exists(bai) or os.path.exists(bai2):
        print(f"{PASS} Index .bai găsit")
    else:
        print(f"{FAIL} Index .bai LIPSEȘTE → rulează: samtools index {BAM_PATH}")
        errors.append("BAM index lipsă")

    bam.close()
except Exception as e:
    print(f"{FAIL} BAM ERROR: {e}")
    errors.append(f"BAM: {e}")
    sys.exit(1)

# ===========================================================================
# TEST 2 — VCF: primele variante și formatul contig-urilor
# ===========================================================================
print("\n" + "="*60)
print("TEST 2: VCF — format contiguri și genotype")
print("="*60)
try:
    vcf = cyvcf2.VCF(VCF_PATH)
    variants_seen = []
    chrom_formats = set()

    for i, v in enumerate(vcf):
        chrom_formats.add(v.CHROM)
        gt = v.genotypes[0][:2]
        variants_seen.append({
            "chrom": v.CHROM, "pos": v.POS,
            "ref": v.REF, "alt": v.ALT,
            "gt": gt
        })
        if i >= 9:
            break
    vcf.close()

    print(f"{PASS} VCF citit OK — primele 10 variante:")
    for v in variants_seen[:5]:
        print(f"   {v['chrom']}:{v['pos']}  {v['ref']}→{v['alt']}  GT={v['gt']}")

    has_chr = any(c.startswith("chr") for c in chrom_formats)
    no_chr  = any(not c.startswith("chr") for c in chrom_formats)
    print(f"\n   Format contiguri VCF: {'cu chr' if has_chr else 'fără chr'}")

    # Verificăm dacă BAM are același format
    bam = pysam.AlignmentFile(BAM_PATH, "rb")
    bam_refs = set(bam.references[:5])
    bam.close()
    bam_has_chr = any(r.startswith("chr") for r in bam_refs)
    print(f"   Format contiguri BAM: {'cu chr' if bam_has_chr else 'fără chr'}")

    if has_chr != bam_has_chr:
        print(f"{FAIL} MISMATCH contiguri! VCF={'chr...' if has_chr else '1,2,...'} "
              f"BAM={'chr...' if bam_has_chr else '1,2,...'}")
        print(f"   → dataset.py trebuie să normalizeze prefixul 'chr'")
        errors.append("Mismatch chr prefix VCF vs BAM")
    else:
        print(f"{PASS} Format contiguri consistent VCF ↔ BAM")

except Exception as e:
    print(f"{FAIL} VCF ERROR: {e}")
    errors.append(f"VCF: {e}")

# ===========================================================================
# TEST 3 — Pileup real: verificăm că BAM are read-uri la pozițiile din VCF
# ===========================================================================
print("\n" + "="*60)
print("TEST 3: Pileup — BAM are read-uri la pozițiile din VCF?")
print("="*60)
try:
    vcf = cyvcf2.VCF(VCF_PATH)
    test_variants = []
    for v in vcf:
        if len(v.ALT) == 1:
            gt = v.genotypes[0][:2]
            label = -1
            alleles = [a for a in gt if a is not None and a >= 0]
            if alleles:
                if all(a == 0 for a in alleles):   label = 0
                elif len(set(alleles)) > 1:         label = 1
                else:                               label = 2
            if label >= 0:
                chrom = v.CHROM
                if not chrom.startswith("chr"):
                    chrom = "chr" + chrom
                test_variants.append({
                    "chrom": chrom, "pos": v.POS - 1,
                    "ref": v.REF,
                    "alt": v.ALT[0],
                    "label": label
                })
        if len(test_variants) >= 20:
            break
    vcf.close()

    bam = pysam.AlignmentFile(BAM_PATH, "rb")
    bam_chroms = set(bam.references)

    found_reads = 0
    zero_depth  = 0
    chrom_mismatch = 0
    label_counts_test = {0: 0, 1: 0, 2: 0}

    print(f"   Testăm {len(test_variants)} variante:")
    for v in test_variants:
        label_counts_test[v['label']] += 1

        # Verificăm dacă cromozomul există în BAM
        if v['chrom'] not in bam_chroms:
            # Încearcă fără chr
            alt_chrom = v['chrom'].replace("chr", "") if v['chrom'].startswith("chr") else "chr" + v['chrom']
            if alt_chrom in bam_chroms:
                chrom_mismatch += 1
                continue
            else:
                chrom_mismatch += 1
                continue

        depth = bam.count(v['chrom'], v['pos'], v['pos'] + 1)
        if depth > 0:
            found_reads += 1
        else:
            zero_depth += 1

    bam.close()

    print(f"   Label distribution: Ref={label_counts_test[0]} Het={label_counts_test[1]} Hom={label_counts_test[2]}")
    print(f"   Variante cu depth>0  : {found_reads}/{len(test_variants)}")
    print(f"   Variante cu depth=0  : {zero_depth}")
    print(f"   Mismatch contig chr  : {chrom_mismatch}")

    if chrom_mismatch > 0:
        print(f"\n{FAIL} PROBLEMA GĂSITĂ: {chrom_mismatch} variante au contig incompatibil!")
        print(f"   VCF folosește: {test_variants[0]['chrom']}")
        print(f"   BAM are: {list(bam_chroms)[:5]}")
        errors.append(f"Contig mismatch: {chrom_mismatch} variante")
    elif found_reads == 0:
        print(f"\n{FAIL} NICIO variantă nu are read-uri în BAM la pozițiile din VCF!")
        errors.append("Zero depth la toate variantele testate")
    elif found_reads < len(test_variants) * 0.5:
        print(f"\n{WARN} Doar {found_reads}/{len(test_variants)} variante au coverage!")
    else:
        print(f"\n{PASS} {found_reads}/{len(test_variants)} variante au read-uri în BAM")

except Exception as e:
    print(f"{FAIL} Pileup test ERROR: {e}")
    import traceback; traceback.print_exc()
    errors.append(f"Pileup: {e}")

# ===========================================================================
# TEST 4 — Encoding: imaginile sunt diferite pentru clase diferite?
# ===========================================================================
print("\n" + "="*60)
print("TEST 4: Encoding — imaginile diferă între clase?")
print("="*60)
try:
    from dataset import _encode_pileup_2d

    vcf = cyvcf2.VCF(VCF_PATH)
    samples_per_class = {0: None, 1: None, 2: None}

    bam = pysam.AlignmentFile(BAM_PATH, "rb")
    bam_chroms = set(bam.references)

    for v in vcf:
        if len(v.ALT) != 1:
            continue
        gt = v.genotypes[0][:2]
        alleles = [a for a in gt if a is not None and a >= 0]
        if not alleles:
            continue
        if all(a == 0 for a in alleles):   label = 0
        elif len(set(alleles)) > 1:        label = 1
        else:                              label = 2

        if samples_per_class[label] is not None:
            continue

        chrom = v.CHROM
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom

        # Verificăm că cromozomul există
        actual_chrom = chrom
        if chrom not in bam_chroms:
            alt = chrom.replace("chr", "")
            if alt in bam_chroms:
                actual_chrom = alt
            else:
                continue

        depth = bam.count(actual_chrom, v.POS - 1, v.POS)
        if depth == 0:
            continue

        img = _encode_pileup_2d(bam, actual_chrom, v.POS - 1, ref_base=v.REF)
        samples_per_class[label] = img

        if all(s is not None for s in samples_per_class.values()):
            break

    bam.close()
    vcf.close()

    class_names = ["Ref", "Het", "Hom-Alt"]
    print(f"\n   Statistici per imagine (shape={img.shape}):")
    print(f"   {'Clasa':<10} {'Min':>8} {'Max':>8} {'Mean':>8} {'NonZero%':>10} {'Canal1_mean':>12}")

    all_zeros = True
    for label, name in enumerate(class_names):
        img = samples_per_class[label]
        if img is None:
            print(f"   {name:<10} {'N/A — nu s-a găsit un exemplu valid':}")
            continue

        nonzero_pct = np.count_nonzero(img) / img.size * 100
        canal1_mean = img[1].mean()  # canalul "este ALT?"

        print(f"   {name:<10} {img.min():>8.3f} {img.max():>8.3f} "
              f"{img.mean():>8.4f} {nonzero_pct:>9.1f}% {canal1_mean:>12.4f}")

        if img.max() > 0:
            all_zeros = False

    if all_zeros:
        print(f"\n{FAIL} TOATE IMAGINILE SUNT ZEROS!")
        print("   → Problema e în pileup encoding sau mismatch de contig")
        errors.append("Toate imaginile sunt zeros")
    else:
        # Verificăm dacă imaginile diferă între clase
        imgs = [samples_per_class[l] for l in range(3) if samples_per_class[l] is not None]
        if len(imgs) >= 2:
            diff_01 = np.abs(imgs[0] - imgs[1]).mean() if len(imgs) > 1 else 0
            diff_02 = np.abs(imgs[0] - imgs[2]).mean() if len(imgs) > 2 else 0
            diff_12 = np.abs(imgs[1] - imgs[2]).mean() if len(imgs) > 2 else 0
            print(f"\n   Diferențe medii între clase (ar trebui > 0):")
            print(f"   Ref vs Het    : {diff_01:.4f}")
            print(f"   Ref vs Hom    : {diff_02:.4f}")
            print(f"   Het vs Hom    : {diff_12:.4f}")
            if diff_01 < 1e-6 and diff_12 < 1e-6:
                print(f"\n{FAIL} Imaginile sunt IDENTICE între clase! Encodingul nu funcționează.")
                errors.append("Imagini identice între clase")
            else:
                print(f"\n{PASS} Imaginile diferă între clase")

except Exception as e:
    print(f"{FAIL} Encoding test ERROR: {e}")
    import traceback; traceback.print_exc()
    errors.append(f"Encoding: {e}")

# ===========================================================================
# TEST 5 — Model: gradienții circulă?
# ===========================================================================
print("\n" + "="*60)
print("TEST 5: Model — gradienții circulă corect?")
print("="*60)
try:
    from model import VariantCallerCNN
    import torch.nn as nn

    model = VariantCallerCNN()
    model.train()

    # Batch sintetic cu valori diferite per clasă (simulăm date reale)
    x = torch.zeros(6, 6, 100, 100)
    for c in range(6):
        x[c*2, 1, 50, 50]     = 1.0   # canal 1 (ALT) = 1.0 → Hom
        x[c*2+1, 1, 50, 50]   = 0.5   # canal 1 (ALT) = 0.5 → Het
    y = torch.tensor([0, 0, 1, 1, 2, 2])

    optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
    criterion = nn.CrossEntropyLoss()

    losses = []
    for step in range(10):
        optimizer.zero_grad()
        out  = model(x)
        loss = criterion(out, y)
        loss.backward()
        optimizer.step()
        losses.append(loss.item())

    print(f"   Loss la step 0  : {losses[0]:.4f}")
    print(f"   Loss la step 9  : {losses[-1]:.4f}")

    if losses[-1] < losses[0] - 0.01:
        print(f"{PASS} Modelul învață pe date sintetice (loss scade)")
    elif abs(losses[-1] - losses[0]) < 0.001:
        print(f"{FAIL} Modelul NU scade loss-ul nici pe date sintetice!")
        errors.append("Model nu învață nici sintetic")
    else:
        print(f"{WARN} Scădere mică ({losses[0]:.4f}→{losses[-1]:.4f}), posibil normal pe 10 steps")

    # Verificăm gradienți
    total_grad = sum(
        p.grad.abs().sum().item()
        for p in model.parameters()
        if p.grad is not None
    )
    print(f"   Suma gradienților: {total_grad:.4f}")
    if total_grad < 1e-10:
        print(f"{FAIL} Gradienți aproape ZERO!")
        errors.append("Gradienți zero")
    else:
        print(f"{PASS} Gradienți OK")

except Exception as e:
    print(f"{FAIL} Model test ERROR: {e}")
    import traceback; traceback.print_exc()
    errors.append(f"Model: {e}")

# ===========================================================================
# SUMAR FINAL
# ===========================================================================
print("\n" + "="*60)
print("SUMAR DIAGNOSTICARE")
print("="*60)
if not errors:
    print(f"{PASS} Toate testele au trecut! Poți rula train.py")
    print("\nDacă training-ul tot nu merge, încearcă:")
    print("  python train.py --lr 1e-3 --batch_size 16")
else:
    print(f"{FAIL} {len(errors)} problemă(e) găsită(e):\n")
    for i, err in enumerate(errors, 1):
        print(f"  {i}. {err}")

    print("\n--- SOLUȚII SUGERATE ---")
    for err in errors:
        if "mismatch" in err.lower() or "contig" in err.lower():
            print("""
  FIX contig mismatch:
  Verifică în debug output mai sus dacă BAM-ul tău folosește 'chr1' sau '1'.
  În dataset.py, funcția _encode_pileup_2d primește deja chrom normalizat,
  dar BAM-ul poate folosi altă convenție.
  Înlocuiește în _encode_pileup_2d:
    bam.pileup(chrom, ...)
  cu verificare automată:
    actual = chrom if chrom in bam.references else chrom.replace('chr','')
    bam.pileup(actual, ...)
            """)
        if "zeros" in err.lower():
            print("""
  FIX imagini zeros:
  Cel mai probabil cauza e contig mismatch (fix de mai sus).
  Dacă cromozomii se potrivesc dar tot zeros → verifică că BAM-ul
  are coverage (samtools depth -r chr1:100000-200000 file.bam)
            """)
        if "identice" in err.lower():
            print("""
  FIX imagini identice:
  ref_base nu se transmite corect sau encodingul canal 1 e 0 peste tot.
  Verifică că v.REF din VCF ajunge la _encode_pileup_2d ca ref_base.
            """)