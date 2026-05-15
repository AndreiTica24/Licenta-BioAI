import pysam
import os

vcf_path = "data/HG002_Son/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf"
bam_path = "data/HG002_Son/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"

def find_first_covered_mutation():
    print("⏳ Deschidem fișierele...")
    vcf_in = pysam.VariantFile(vcf_path)
    bam_in = pysam.AlignmentFile(bam_path, "rb")

    print("🔎 Căutăm prima mutație care are date (citiri) în fișierul BAM. Poate dura câteva secunde...")
    
    # Parcurgem VCF-ul rând cu rând
    for record in vcf_in:
        chrom = record.chrom
        pos = record.pos
        ref = record.ref
        
        # Ignorăm mutațiile complexe pentru moment, luăm doar prima variantă
        alt = record.alts[0] 
        
        # Traducem cromozomul
        target_chrom = chrom.replace("chr", "") if chrom.startswith("chr") else chrom
        if target_chrom not in bam_in.references:
            target_chrom = "chr" + target_chrom if not chrom.startswith("chr") else chrom

        has_coverage = False
        
        # Verificăm dacă există citiri în BAM la această poziție
        for pileupcolumn in bam_in.pileup(target_chrom, pos-1, pos):
            # Dacă suntem fix pe poziția mutației și avem Depth > 0
            if pileupcolumn.pos == pos - 1 and pileupcolumn.n > 0:
                print(f"\n🎉 SUCCES! Am găsit o mutație acoperită de aparatul de secvențiere:")
                print(f"Cromozom: {chrom}, Poziția: {pos} (Referință: {ref} -> Mutație: {alt})")
                print(f"🔬 Analiză la fața locului:")
                print(f"Număr de citiri (Depth): {pileupcolumn.n}")
                
                bases = []
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                        bases.append(base)
                
                print(f"Primele 20 de litere: {bases[:20]}")
                
                count_ref = bases.count(ref.upper())
                count_alt = bases.count(alt.upper())
                print(f"👉 Dintre acestea: {count_ref} sunt '{ref}' (Normale) și {count_alt} sunt '{alt}' (Mutația).")
                
                has_coverage = True
                break # Ieșim din for-ul de BAM
        
        if has_coverage:
            break # Ieșim din for-ul de VCF odată ce am găsit prima mutație bună

if __name__ == "__main__":
    find_first_covered_mutation()