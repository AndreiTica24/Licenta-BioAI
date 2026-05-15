import torch
from torch.utils.data import Dataset, DataLoader
import pysam
import numpy as np
import random

class GenomicDataset(Dataset):
    def __init__(self, bam_path, vcf_path, window_size=50, max_depth=50, max_samples=1000):
        self.bam_path = bam_path
        self.vcf_path = vcf_path
        self.window_size = window_size
        self.max_depth = max_depth
        self.base_to_num = {'A': 0.25, 'C': 0.50, 'G': 0.75, 'T': 1.0, 'N': 0.0}
        
        print(f"⏳ Construim un dataset inteligent ({max_samples} exemple: 50% Mutații, 50% Normale)...")
        self.data_points = []
        
        vcf_in = pysam.VariantFile(self.vcf_path)
        mutations_needed = max_samples // 2
        
        # Ținem minte unde sunt mutațiile ca să nu picăm pe ele accidental
        mutation_positions = set()
        
        # 1. Extragem Mutațiile (Eticheta 1)
        for i, record in enumerate(vcf_in):
            if i >= mutations_needed:
                break
            self.data_points.append({'chrom': record.chrom, 'pos': record.pos, 'label': 1})
            mutation_positions.add(record.pos)
            
        # 2. Generăm Zone Normale (Eticheta 0)
        # Ne deplasăm cu 500 de baze față de mutații pentru a găsi ADN normal
        for mut in self.data_points[:mutations_needed]:
            normal_pos = mut['pos'] + 500
            
            # Ne asigurăm că nu am picat peste o altă mutație din greșeală
            while normal_pos in mutation_positions:
                normal_pos += 50 
                
            self.data_points.append({'chrom': mut['chrom'], 'pos': normal_pos, 'label': 0})
            
        # 3. Amestecăm perfect datele
        random.shuffle(self.data_points)
        print(f"✅ Am pregătit {len(self.data_points)} exemple amestecate! AI-ul nu mai poate trișa.")

    def __len__(self):
        return len(self.data_points)

    def __getitem__(self, idx):
        dp = self.data_points[idx]
        chrom = dp['chrom']
        pos = dp['pos']
        
        bam_in = pysam.AlignmentFile(self.bam_path, "rb")
        width = self.window_size * 2 + 1
        tensor = np.zeros((3, self.max_depth, width), dtype=np.float32)
        
        target_chrom = chrom.replace("chr", "") if chrom.startswith("chr") else chrom
        if target_chrom not in bam_in.references:
            target_chrom = "chr" + target_chrom if not chrom.startswith("chr") else chrom

        start_pos = pos - self.window_size
        end_pos = pos + self.window_size + 1
        
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
        
        x_tensor = torch.from_numpy(tensor)
        y_label = torch.tensor(dp['label'], dtype=torch.long)
        
        return x_tensor, y_label

if __name__ == "__main__":
    bam_file = "data/HG002_Son/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"
    vcf_file = "data/HG002_Son/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf"
    
    my_dataset = GenomicDataset(bam_file, vcf_file, max_samples=20)
    my_dataloader = DataLoader(my_dataset, batch_size=5, shuffle=True)
    
    for images, labels in my_dataloader:
        print(f"\n📦 Pachet de test - Etichete (0=Normal, 1=Mutație): {labels}")
        break