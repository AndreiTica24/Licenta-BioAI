"""
model.py — VariantCallerCNN
Rețea CNN pentru clasificarea variantelor genomice în 3 clase:
  0 = Ref (homozigot referință)
  1 = Het (heterozigot)
  2 = Hom-Alt (homozigot alternativ)

Input:  tensor (B, 6, 100, 100) — 6 canale pileup
Output: logits (B, 3)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


# ---------------------------------------------------------------------------
# Bloc rezidual
# ---------------------------------------------------------------------------

class ResidualBlock(nn.Module):
    """Bloc rezidual standard cu BN + ReLU."""

    def __init__(self, channels: int, dropout: float = 0.1):
        super().__init__()
        self.conv1 = nn.Conv2d(channels, channels, 3, padding=1, bias=False)
        self.bn1   = nn.BatchNorm2d(channels)
        self.conv2 = nn.Conv2d(channels, channels, 3, padding=1, bias=False)
        self.bn2   = nn.BatchNorm2d(channels)
        self.drop  = nn.Dropout2d(dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        residual = x
        out = F.relu(self.bn1(self.conv1(x)), inplace=True)
        out = self.drop(out)
        out = self.bn2(self.conv2(out))
        return F.relu(out + residual, inplace=True)


# ---------------------------------------------------------------------------
# Atenție pe canale (Squeeze-and-Excitation)
# ---------------------------------------------------------------------------

class SEBlock(nn.Module):
    """Squeeze-and-Excitation: re-calibrează importanța canalelor."""

    def __init__(self, channels: int, reduction: int = 4):
        super().__init__()
        self.pool = nn.AdaptiveAvgPool2d(1)
        self.fc   = nn.Sequential(
            nn.Linear(channels, channels // reduction, bias=False),
            nn.ReLU(inplace=True),
            nn.Linear(channels // reduction, channels, bias=False),
            nn.Sigmoid(),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        b, c, _, _ = x.shape
        w = self.pool(x).view(b, c)
        w = self.fc(w).view(b, c, 1, 1)
        return x * w


# ---------------------------------------------------------------------------
# Model principal
# ---------------------------------------------------------------------------

class VariantCallerCNN(nn.Module):
    """
    CNN cu blocuri reziduale și SE-attention pentru variant calling.

    Arhitectură:
      Stem  → 3 × (Conv-BN-ReLU + MaxPool + ResBlock + SE)
      Head  → GlobalAvgPool → FC(256) → Dropout → FC(3)

    Parametri (~1.2M) — potrivit pentru dataset de ~50k exemple.
    """

    def __init__(
        self,
        in_channels:  int   = 6,
        num_classes:  int   = 3,
        base_filters: int   = 32,
        dropout:      float = 0.3,
    ):
        super().__init__()

        f = base_filters  # 32

        # --- Stem ---
        self.stem = nn.Sequential(
            nn.Conv2d(in_channels, f, kernel_size=7, stride=2, padding=3, bias=False),
            nn.BatchNorm2d(f),
            nn.ReLU(inplace=True),
            nn.MaxPool2d(3, stride=2, padding=1),
        )

        # --- Stage 1: 32 → 64 ---
        self.stage1 = nn.Sequential(
            nn.Conv2d(f, f * 2, 3, stride=2, padding=1, bias=False),
            nn.BatchNorm2d(f * 2),
            nn.ReLU(inplace=True),
            ResidualBlock(f * 2, dropout=0.1),
            SEBlock(f * 2),
        )

        # --- Stage 2: 64 → 128 ---
        self.stage2 = nn.Sequential(
            nn.Conv2d(f * 2, f * 4, 3, stride=2, padding=1, bias=False),
            nn.BatchNorm2d(f * 4),
            nn.ReLU(inplace=True),
            ResidualBlock(f * 4, dropout=0.1),
            ResidualBlock(f * 4, dropout=0.1),
            SEBlock(f * 4),
        )

        # --- Stage 3: 128 → 256 ---
        self.stage3 = nn.Sequential(
            nn.Conv2d(f * 4, f * 8, 3, stride=2, padding=1, bias=False),
            nn.BatchNorm2d(f * 8),
            nn.ReLU(inplace=True),
            ResidualBlock(f * 8, dropout=0.1),
            SEBlock(f * 8),
        )

        # --- Head ---
        self.pool    = nn.AdaptiveAvgPool2d(1)
        self.flatten = nn.Flatten()
        self.head    = nn.Sequential(
            nn.Linear(f * 8, 256),
            nn.ReLU(inplace=True),
            nn.Dropout(dropout),
            nn.Linear(256, num_classes),
        )

        self._init_weights()

    # ------------------------------------------------------------------
    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, mode="fan_out",
                                        nonlinearity="relu")
            elif isinstance(m, nn.BatchNorm2d):
                nn.init.ones_(m.weight)
                nn.init.zeros_(m.bias)
            elif isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

    # ------------------------------------------------------------------
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        x: (B, 6, H, W)
        returns: logits (B, 3)
        """
        x = self.stem(x)
        x = self.stage1(x)
        x = self.stage2(x)
        x = self.stage3(x)
        x = self.pool(x)
        x = self.flatten(x)
        return self.head(x)

    # ------------------------------------------------------------------
    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        """Returnează probabilități softmax (B, 3)."""
        with torch.no_grad():
            return F.softmax(self.forward(x), dim=1)

    # ------------------------------------------------------------------
    @staticmethod
    def load_from_checkpoint(path: str,
                             device: str = "cpu") -> "VariantCallerCNN":
        """Încarcă modelul dintr-un checkpoint salvat de train.py."""
        checkpoint = torch.load(path, map_location=device)
        model = VariantCallerCNN()
        if "model_state_dict" in checkpoint:
            model.load_state_dict(checkpoint["model_state_dict"])
            print(f"✅ Model încărcat din epoch {checkpoint.get('epoch', '?')}"
                  f"  (val_loss={checkpoint.get('val_loss', '?'):.4f},"
                  f"  F1={checkpoint.get('f1', '?'):.4f})")
        else:
            model.load_state_dict(checkpoint)
        model.to(device)
        return model


# ---------------------------------------------------------------------------
# Smoke test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    model = VariantCallerCNN()
    total = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Parametri antrenabili: {total:,}")

    x = torch.randn(4, 6, 100, 100)
    out = model(x)
    print(f"Input:  {x.shape}")
    print(f"Output: {out.shape}")   # (4, 3)
    print(f"Probs:  {model.predict_proba(x)[0]}")