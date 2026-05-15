import torch
import pysam
import pandas as pd

print("=== Testare Mediu AI ===")
print(f"Versiune PyTorch: {torch.__version__}")
print(f"Versiune pysam: {pysam.__version__}")

# Verificăm dacă PyTorch vede placa ta video NVIDIA
cuda_available = torch.cuda.is_available()
print(f"Placa video NVIDIA este detectată de AI? : {cuda_available}")

if cuda_available:
    print(f"Numele plăcii video: {torch.cuda.get_device_name(0)}")