import pysam
import numpy as np
import torch

bam_path = "data/HG002_Son/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"

# Dicționarul de traducere: Biologie -> Matematică
base_to_num = {'A': 0.25, 'C': 0.50, 'G': 0.75, 'T': 1.0, 'N': 0.0}

def create_tensor(chrom, pos, window_size=50, max_depth=50):
    bam_in = pysam.AlignmentFile(bam_path, "rb")
    
    # Lățimea totală a ferestrei (50 la stânga + poziția centrală + 50 la dreapta = 101)
    width = window_size * 2 + 1
    
    # Inițializăm o matrice goală cu zerouri de formă (3 canale, 50 citiri, 101 lățime)
    tensor = np.zeros((3, max_depth, width), dtype=np.float32)
    
    target_chrom = chrom.replace("chr", "") if chrom.startswith("chr") else chrom
    if target_chrom not in bam_in.references:
        target_chrom = "chr" + target_chrom if not chrom.startswith("chr") else chrom

    start_pos = pos - window_size
    end_pos = pos + window_size + 1
    
    print(f"⏳ Extragem matricea pentru {chrom}:{pos} (Zona: {start_pos} - {end_pos})...")

    # Parcurgem fiecare coloană din această fereastră
    for pileupcolumn in bam_in.pileup(target_chrom, start_pos-1, end_pos):
        col_idx = pileupcolumn.pos - (start_pos - 1)
        
        # Dacă ne aflăm în interiorul ferestrei noastre de 101 baze
        if 0 <= col_idx < width:
            read_idx = 0
            for pileupread in pileupcolumn.pileups:
                if read_idx >= max_depth:
                    break # Ne oprim dacă avem mai mult de 50 de citiri
                    
                if not pileupread.is_del and not pileupread.is_refskip:
                    # Extragem valorile brute din BAM
                    base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                    qual = pileupread.alignment.query_qualities[pileupread.query_position]
                    strand = -1.0 if pileupread.alignment.is_reverse else 1.0
                    
                    # Umplem matricea (Tensorul)
                    tensor[0, read_idx, col_idx] = base_to_num.get(base, 0.0)  # Canal 0: Baza
                    tensor[1, read_idx, col_idx] = min(qual / 40.0, 1.0)       # Canal 1: Calitatea
                    tensor[2, read_idx, col_idx] = strand                      # Canal 2: Direcția
                    
                read_idx += 1

    # Convertim matricea de Numpy în Tensor de PyTorch pentru placa video
    pytorch_tensor = torch.from_numpy(tensor)
    return pytorch_tensor

if __name__ == "__main__":
    # Testăm exact pe mutația pe care am găsit-o mai devreme!
    test_chrom = "chr1"
    test_pos = 604358
    
    # Generăm imaginea matematică
    my_tensor = create_tensor(test_chrom, test_pos, window_size=50, max_depth=50)
    
    print("\n✅ Tensor generat cu succes!")
    print(f"Structura Tensorului: {my_tensor.shape} -> (Canale, Citiri, Baze)")
    
    print("\n🔬 Hai să ne uităm la Canalul 0 (Bazele ADN) exact pe coloana centrală (mutația):")
    # Afișăm primele 5 rânduri din coloana de mijloc (indexul 50)
    print(my_tensor[0, :5, 50])