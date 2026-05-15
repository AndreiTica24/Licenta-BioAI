import torch
from torch.utils.data import Dataset, DataLoader
import pysam
import numpy as np
import os

class GenomicDataset(Dataset):
    def __init__(self, bam_path, vcf_path, window_size=50, max_depth=50, max_samples=1000):
        """
        Inițializează clasa citind pozițiile mutațiilor din VCF.
        Folosim max_samples pentru a nu încărca milioane de mutații deodată în faza de test.
        """
        self.bam_path = bam_path
        self.vcf_path = vcf_path
        self.window_size = window_size
        self.max_depth = max_depth
        self.base_to_num = {'A': 0.25, 'C': 0.50, 'G': 0.75, 'T': 1.0, 'N': 0.0}
        
        print(f"⏳ Încărcăm maxim {max_samples} mutații din VCF...")
        self.mutations = []
        
        vcf_in = pysam.VariantFile(self.vcf_path)
        for i, record in enumerate(vcf_in):
            if i >= max_samples:
                break
            
            self.mutations.append({
                'chrom': record.chrom,
                'pos': record.pos,
                # Atribuim eticheta 1 (Mutație prezentă) pentru toate înregistrările din acest VCF de benchmark
                'label': 1 
            })
            
        print(f"✅ Am pregătit o listă de {len(self.mutations)} exemple.")

    def __len__(self):
        # PyTorch are nevoie să știe dimensiunea totală a setului de date
        return len(self.mutations)

    def __getitem__(self, idx):
        # Această funcție este apelată automat de PyTorch în timpul antrenamentului
        # pentru a extrage un singur exemplu la un moment dat.
        
        mut = self.mutations[idx]
        chrom = mut['chrom']
        pos = mut['pos']
        
        bam_in = pysam.AlignmentFile(self.bam_path, "rb")
        width = self.window_size * 2 + 1
        tensor = np.zeros((3, self.max_depth, width), dtype=np.float32)
        
        # Corectăm numele cromozomului (chr1 vs 1)
        target_chrom = chrom.replace("chr", "") if chrom.startswith("chr") else chrom
        if target_chrom not in bam_in.references:
            target_chrom = "chr" + target_chrom if not chrom.startswith("chr") else chrom

        start_pos = pos - self.window_size
        end_pos = pos + self.window_size + 1
        
        # Construim imaginea
        for pileupcolumn in bam_in.pileup(target_chrom, start_pos-1, end_pos):
            col_idx = pileupcolumn.pos - (start_pos - 1)
            
            if 0 <= col_idx < width:
                read_idx = 0
                for pileupread in pileupcolumn.pileups:
                    if read_idx >= self.max_depth: 
                        break
                        
                    if not pileupread.is_del and not pileupread.is_refskip:
                        base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                        qual = pileupread.alignment.query_qualities[pileupread.query_position]
                        strand = -1.0 if pileupread.alignment.is_reverse else 1.0
                        
                        tensor[0, read_idx, col_idx] = self.base_to_num.get(base, 0.0)
                        tensor[1, read_idx, col_idx] = min(qual / 40.0, 1.0)
                        tensor[2, read_idx, col_idx] = strand
                        
                    read_idx += 1
                    
        bam_in.close()
        
        # Returnăm tensorul de date și tensorul etichetei
        x_tensor = torch.from_numpy(tensor)
        y_label = torch.tensor(mut['label'], dtype=torch.long)
        
        return x_tensor, y_label

# --- Blocul de Testare ---
if __name__ == "__main__":
    bam_file = "data/HG002_Son/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"
    vcf_file = "data/HG002_Son/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf"
    
    # 1. Inițializăm Dataset-ul (încărcăm doar 100 pentru viteză)
    my_dataset = GenomicDataset(bam_file, vcf_file, max_samples=100)
    
    # 2. Configurăm DataLoader-ul (Banda Rulantă)
    # batch_size=16 înseamnă că grupăm câte 16 mutații într-un singur pachet pentru eficiența GPU-ului
    my_dataloader = DataLoader(my_dataset, batch_size=16, shuffle=True)
    
    # 3. Tragem de pe bandă primul pachet
    for images, labels in my_dataloader:
        print("\n📦 Am extras primul pachet (Batch):")
        print(f"Forma Tensorului: {images.shape}")
        print(f"Semnificație: (Dimensiune Pachet, Canale, Înălțime/Citiri, Lățime/Baze)")
        print(f"Etichete (Labels): {labels}")
        break