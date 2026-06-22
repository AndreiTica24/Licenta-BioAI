"""
model.py — VariantCallerCNN1D
Rețea CNN 1D pentru clasificarea variantelor genomice în 3 clase:
  0 = Ref (homozigot referință)
  1 = Het (heterozigot)
  2 = Hom-Alt (homozigot alternativ)

Arhitectura urmează literatura standardă din genomic deep learning:
DeepBind (Alipanahi et al. 2015), DeepSEA (Zhou & Troyanskaya 2015),
Basset (Kelley et al. 2016).

Input:  tensor (B, 6, 200) — 6 canale per poziție genomică
Output: logits (B, 3)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


class VariantCallerCNN1D(nn.Module):
    def __init__(self,
                 in_channels: int = 10,
                 num_classes: int = 3,
                 dropout:     float = 0.3):
        super().__init__()

        self.block1 = nn.Sequential(
            nn.Conv1d(in_channels, 64, kernel_size=11, padding=5, bias=False),
            nn.BatchNorm1d(64),
            nn.ReLU(inplace=True),
            nn.Conv1d(64, 64, kernel_size=11, padding=5, bias=False),
            nn.BatchNorm1d(64),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(2),
        )

        self.block2 = nn.Sequential(
            nn.Conv1d(64, 128, kernel_size=7, padding=3, bias=False),
            nn.BatchNorm1d(128),
            nn.ReLU(inplace=True),
            nn.Conv1d(128, 128, kernel_size=7, padding=3, bias=False),
            nn.BatchNorm1d(128),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(2),
        )

        self.block3 = nn.Sequential(
            nn.Conv1d(128, 256, kernel_size=5, padding=2, bias=False),
            nn.BatchNorm1d(256),
            nn.ReLU(inplace=True),
            nn.Conv1d(256, 256, kernel_size=5, padding=2, bias=False),
            nn.BatchNorm1d(256),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(2),
        )

        self.global_pool = nn.AdaptiveAvgPool1d(1)
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
            if isinstance(m, nn.Conv1d):
                nn.init.kaiming_normal_(m.weight, mode="fan_out",
                                        nonlinearity="relu")
            elif isinstance(m, nn.BatchNorm1d):
                nn.init.ones_(m.weight)
                nn.init.zeros_(m.bias)
            elif isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        x: (B, 10, 200) — 10 canale (6 din BAM + 4 din referință), 200 poziții
        returns: logits (B, 3)
        """
        assert x.ndim == 3 and x.shape[1] == 10, \
            f"Input așteptat (B, 10, L), primit {tuple(x.shape)}"

        x = self.block1(x)
        x = self.block2(x)
        x = self.block3(x)
        x = self.global_pool(x)
        return self.classifier(x)

    def predict_proba(self, x: torch.Tensor) -> torch.Tensor:
        """Probabilități softmax (B, 3)."""
        with torch.no_grad():
            return F.softmax(self.forward(x), dim=1)

    @staticmethod
    def load_from_checkpoint(path: str,
                             device: str = "cpu") -> "VariantCallerCNN1D":
        """Încarcă model dintr-un checkpoint."""
        ckpt = torch.load(path, map_location=device)
        model = VariantCallerCNN1D()
        if "model_state_dict" in ckpt:
            model.load_state_dict(ckpt["model_state_dict"])
            print(f"✅ Model încărcat epoch={ckpt.get('epoch','?')} "
                  f"F1={ckpt.get('f1', 0):.4f}")
        else:
            model.load_state_dict(ckpt)
        model.to(device)
        return model


# ---------------------------------------------------------------------------
# Smoke test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    model = VariantCallerCNN1D()
    total = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"Parametri antrenabili: {total:,}")

    x = torch.randn(8, 10, 200)
    out = model(x)
    print(f"Input:  {x.shape}")
    print(f"Output: {out.shape}")
    print(f"Probs:  {model.predict_proba(x)[0]}")