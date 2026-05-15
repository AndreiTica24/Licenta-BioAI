import torch
import torch.nn as nn
import torch.nn.functional as F

class VariantCallerCNN(nn.Module):
    def __init__(self):
        super(VariantCallerCNN, self).__init__()
        
        # 1. Stratul de intrare (Conv1)
        # Primim 3 canale (Bază, Calitate, Direcție)
        self.conv1 = nn.Conv2d(in_channels=3, out_channels=32, kernel_size=3, padding=1)
        
        # 2. Stratul Conv2
        self.conv2 = nn.Conv2d(32, 64, kernel_size=3, padding=1)
        
        # 3. Stratul Conv3
        self.conv3 = nn.Conv2d(64, 128, kernel_size=3, padding=1)
        
        # 4. Pooling: Reduce dimensiunea imaginii la jumătate
        self.pool = nn.MaxPool2d(2, 2)
        
        # 5. Dropout: Previne overfitting-ul (memorarea datelor)
        self.dropout = nn.Dropout(0.25)
        
        # 6. Straturile de decizie (Fully Connected)
        # Calculul dimensiunii după 3 operații de MaxPool2d pe o imagine de 50x101:
        # 50x101 -> 25x50 -> 12x25 -> 6x12
        # 128 canale * 6 înălțime * 12 lățime = 9216
        self.fc1 = nn.Linear(9216, 512) 
        self.fc2 = nn.Linear(512, 3) # 3 clase: 0=Ref, 1=Het, 2=Hom-Alt

    def forward(self, x):
        x = self.pool(F.relu(self.conv1(x)))
        x = self.pool(F.relu(self.conv2(x)))
        x = self.pool(F.relu(self.conv3(x)))
        
        # Platizăm matricea 3D într-un vector 1D
        x = torch.flatten(x, 1)
        
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        x = self.fc2(x) # Returnăm scorurile pentru cele 3 clase
        
        return x