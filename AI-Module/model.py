"""
model.py — VariantCallerCNN v2 (simplificat)

Față de v1 (cu ResidualBlock + SE):
  - Arhitectură mai simplă: 4 blocuri Conv+BN+ReLU+Pool
  - Mai puțini parametri (~500K vs 2.3M)
  - Convergență mai rapidă pentru semnalul nostru (AF din pattern pileup)

Input:  (B, 6, 100, 100)
Output: logits (B, 3)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


class VariantCallerCNN(nn.Module):
    def __init__(self,
                 in_channels: int = 6,
                 num_classes: int = 3,
                 dropout:     float = 0.3):
        super().__init__()

        # Block 1: 6 → 32, 100×100 → 50×50
        self.block1 = nn.Sequential(
            nn.Conv2d(in_channels, 32, 3, padding=1, bias=False),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True),
            nn.Conv2d(32, 32, 3, padding=1, bias=False),
            nn.BatchNorm2d(32),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(2),
        )

        # Block 2: 32 → 64, 50×50 → 25×25
        self.block2 = nn.Sequential(
            nn.Conv2d(32, 64, 3, padding=1, bias=False),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True),
            nn.Conv2d(64, 64, 3, padding=1, bias=False),
            nn.BatchNorm2d(64),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(2),
        )

        # Block 3: 64 → 128, 25×25 → 12×12
        self.block3 = nn.Sequential(
            nn.Conv2d(64, 128, 3, padding=1, bias=False),
            nn.BatchNorm2d(128),
            nn.ReLU(inplace=True),
            nn.Conv2d(128, 128, 3, padding=1, bias=False),
            nn.BatchNorm2d(128),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(2),
        )

        # Block 4: 128 → 256, 12×12 → 6×6
        self.block4 = nn.Sequential(
            nn.Conv2d(128, 256, 3, padding=1, bias=False),
            nn.BatchNorm2d(256),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(2),
        )

        # Head
        self.global_pool = nn.AdaptiveAvgPool2d(1)
        self.classifier  = nn.Sequential(
            nn.Flatten(),
            nn.Dropout(dropout),
            nn.Linear(256, 128),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(128, num_classes),
        )

        self._init_weights()

    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, mode="fan_out", nonlinearity="relu")
            elif isinstance(m, nn.BatchNorm2d):
                nn.init.ones_(m.weight)
                nn.init.zeros_(m.bias)
            elif isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        assert x.ndim == 4 and x.shape[1] == 6, \
            f"Input așteptat (B,6,H,W), primit {tuple(x.shape)}"
        x = self.block1(x)
        x = self.block2(x)
        x = self.block3(x)
        x = self.block4(x)
        x = self.global_pool(x)
        return self.classifier(x)

    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        with torch.no_grad():
            return F.softmax(self.forward(x), dim=1)

    @staticmethod
    def load_from_checkpoint(path: str, device: str = "cpu") -> "VariantCallerCNN":
        ckpt = torch.load(path, map_location=device)
        model = VariantCallerCNN()
        if "model_state_dict" in ckpt:
            model.load_state_dict(ckpt["model_state_dict"])
            print(f"✅ Model încărcat epoch={ckpt.get('epoch','?')} "
                  f"F1={ckpt.get('f1', 0):.4f}")
        else:
            model.load_state_dict(ckpt)
        model.to(device)
        return model


if __name__ == "__main__":
    model = VariantCallerCNN()
    total = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Parametri: {total:,}")
    x = torch.randn(4, 6, 100, 100)
    print(f"Input:  {x.shape}")
    print(f"Output: {model(x).shape}")