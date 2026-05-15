import torch
print(f"Versiune Torch: {torch.__version__}")
print(f"CUDA disponibil: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"Placa video: {torch.cuda.get_device_name(0)}")
    # Testăm un calcul mic
    x = torch.rand(5, 3).cuda()
    print("✅ Calculul pe GPU a reușit!")