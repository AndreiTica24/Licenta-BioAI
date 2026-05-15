import pysam
import os

# Definim calea către fișierul VCF al Fiului
# Calea se potrivește cu structura pe care mi-ai arătat-o
vcf_path = "data/HG002_Son/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf"

def explore_vcf(file_path):
    print(f"--- Încercăm să deschidem: {file_path} ---")
    
    # Verificăm dacă fișierul există la adresa dată
    if not os.path.exists(file_path):
        print("❌ Eroare: Nu găsesc fișierul! Verifică dacă numele și calea sunt corecte.")
        return

    try:
        # Deschidem fișierul VCF cu pysam
        vcf_in = pysam.VariantFile(file_path)
        
        print("✅ Fișier deschis cu succes! Iată primele 5 mutații găsite:\n")
        
        # Parcurgem primele 5 înregistrări (mutații)
        for i, record in enumerate(vcf_in):
            if i >= 5:
                break
            
            # record.chrom = Cromozomul
            # record.pos = Poziția exactă
            # record.ref = Litera normală (din referință)
            # record.alts = Litera modificată (mutația)
            print(f"Mutația {i+1}: Cromozom {record.chrom}, Poziția {record.pos}")
            print(f"   Referință: {record.ref}  ->  Mutație: {record.alts[0]}")
            print("-" * 40)
            
    except Exception as e:
        print(f"❌ A apărut o eroare la citire: {e}")

if __name__ == "__main__":
    explore_vcf(vcf_path)