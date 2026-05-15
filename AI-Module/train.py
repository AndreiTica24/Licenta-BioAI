import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader

# Importăm clasele pe care le-am creat în celelalte fișiere
from model import VariantCallerCNN
from dataset import GenomicDataset

def train_model():
    # 1. Configurarea dispozitivului (Să ne asigurăm că folosește RTX-ul)
    device = torch.device("cpu")
    print(f"🚀 Pornim antrenamentul pe: {device}")

    # 2. Căile către date
    bam_file = "data/HG002_Son/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"
    vcf_file = "data/HG002_Son/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf"

    # 3. Pregătim Banda Rulantă (Încărcăm 500 de mutații pentru acest test)
    dataset = GenomicDataset(bam_file, vcf_file, max_samples=500)
    dataloader = DataLoader(dataset, batch_size=32, shuffle=True)

    # 4. Aducem "Creierul" și îl mutăm pe placa video
    model = VariantCallerCNN().to(device)

    # 5. Funcția de "Pedeapsă" și "Profesorul"
    # CrossEntropyLoss calculează cât de mult a greșit AI-ul
    criterion = nn.CrossEntropyLoss()
    # Adam este optimizatorul (cel care ajustează neuronii ca să nu mai greșească)
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    # 6. Bucla de Antrenament (Epocile)
    epochs = 5
    print("\n⏳ Începem învățarea...")
    
    for epoch in range(epochs):
        model.train() # Punem modelul în modul de învățare
        running_loss = 0.0
        
        # Extragem pachete de pe banda rulantă
        for batch_idx, (images, labels) in enumerate(dataloader):
            # Mutăm datele pe placa video
            images = images.to(device)
            labels = labels.to(device)
            
            # Resetăm calculele vechi
            optimizer.zero_grad()
            
            # Pasul Înainte: AI-ul face o predicție
            outputs = model(images)
            
            # Calculăm greșeala (Loss)
            loss = criterion(outputs, labels)
            
            # Pasul Înapoi: AI-ul învață din greșeală (Backpropagation)
            loss.backward()
            
            # Actualizăm neuronii
            optimizer.step()
            
            running_loss += loss.item()
            
        # La finalul fiecărei epoci, afișăm cât de mult a greșit (ideal trebuie să scadă)
        avg_loss = running_loss / len(dataloader)
        print(f"📊 Epoca {epoch + 1}/{epochs} - Eroare medie (Loss): {avg_loss:.4f}")

    print("\n✅ Antrenament complet! Modelul a învățat din datele tale.")

    # Salvăm modelul ca să nu o luăm de la zero data viitoare
    torch.save(model.state_dict(), "variant_caller_model.pth")
    print("💾 Modelul a fost salvat ca 'variant_caller_model.pth'.")

if __name__ == "__main__":
    train_model()