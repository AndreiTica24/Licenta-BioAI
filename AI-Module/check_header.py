import pysam

bam_path = "data/HG002_Son/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam"

with pysam.AlignmentFile(bam_path, "rb") as bam:
    print("Numele cromozomilor din fișierul BAM sunt:")
    # Afișăm primii 5 cromozomi din header
    print(bam.references[:5])