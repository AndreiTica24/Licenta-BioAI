import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from model import VariantCallerCNN
from dataset import GenomicDataset

def train_model():
    device = torch.device("cuda")
    print(f"🚀 Antrenament pe: {torch.cuda.get_device_name(0)}")

    # 1. Căile către date (HG002 pentru Antrenament, HG003 pentru Validare)
    # --- HG002 (Train) ---
    bam_train = "data/HG002_Son/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"
    vcf_train = "data/HG002_Son/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf"
    
    # --- HG003 (Validation) ---
    bam_val = "data/HG003_Father/151002_7001448_0359_AC7F6GANXX_Sample_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup.bam" # Asigură-te că aceste căi sunt corecte
    vcf_val = "data/HG003_Father/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf"

    # 2. Crearea Dataset-urilor
    print("📦 Încărcăm datele de antrenament (HG002)...")
    train_dataset = GenomicDataset(bam_train, vcf_train, max_samples=5000)
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)

    print("\n📦 Încărcăm datele de validare (HG003)...")
    val_dataset = GenomicDataset(bam_val, vcf_val, max_samples=1000)
    val_loader = DataLoader(val_dataset, batch_size=32, shuffle=False)

    # 3. Modelul, Optimizatorul și Criticul
    model = VariantCallerCNN().to(device)
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=0.0001)

    print("\n⏳ Începem procesul de învățare...")
    epochs = 20 # Putem crește numărul de epoci acum că avem GPU

    for epoch in range(epochs):
        # FAZA DE ANTRENAMENT
        model.train()
        train_loss = 0.0
        for images, labels in train_loader:
            images, labels = images.to(device), labels.to(device)
            optimizer.zero_grad()
            outputs = model(images)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            train_loss += loss.item()

        # FAZA DE VALIDARE (Aici vedem dacă modelul funcționează pe HG003)
        model.eval()
        val_loss = 0.0
        correct = 0
        total = 0
        with torch.no_grad():
            for images, labels in val_loader:
                images, labels = images.to(device), labels.to(device)
                outputs = model(images)
                loss = criterion(outputs, labels)
                val_loss += loss.item()
                
                _, predicted = torch.max(outputs.data, 1)
                total += labels.size(0)
                correct += (predicted == labels).sum().item()

        avg_train_loss = train_loss / len(train_loader)
        avg_val_loss = val_loss / len(val_loader)
        accuracy = 100 * correct / total

        print(f"📊 Epoca {epoch+1}/{epochs}")
        print(f"   [Train] Loss: {avg_train_loss:.4f} | [Val HG003] Loss: {avg_val_loss:.4f} | Acuratețe: {accuracy:.2f}%")

    # Salvarea modelului final
    torch.save(model.state_dict(), "variant_caller_model.pth")
    print("\n✅ Gata! Modelul antrenat pe HG002 și validat pe HG003 a fost salvat.")

if __name__ == "__main__":
    train_model()