# GenomicsAI

**Detectarea și anotarea variantelor exomice folosind rețele neuronale convoluționale**

Aplicație web integrată pentru identificarea și interpretarea clinică a variantelor genetice din date de secvențiere exomică (WES), dezvoltată ca lucrare de licență la Universitatea Politehnica Timișoara, Facultatea de Automatică și Calculatoare.

---

## Autor

**Andrei Tica** — Anul III, Informatică  
Universitatea Politehnica Timișoara, 2026

**Coordonator științific:** Raul Robu

---

## Descrierea aplicației

GenomicsAI este o aplicație web care integrează cinci componente independente pentru identificarea variantelor genetice din fișiere BAM exomice și clasificarea lor clinică automată:

1. **Backend Spring Boot** — orchestrare, autentificare JWT, RBAC, persistență H2
2. **Modul AI Python (FastAPI + PyTorch)** — inferență CNN 1D pentru clasificarea variantelor
3. **Anotare VEP (Docker)** — anotare clinică prin Ensembl VEP + ClinVar
4. **Asistent conversațional (Ollama + Llama 3.1 8B)** — rulat local pentru păstrarea confidențialității
5. **Frontend Thymeleaf** — interfață utilizator pentru medici

**Performanță model:** F1 macro = 0,9657 pe setul de validare HG004 (GIAB).

---

## Cerințe de sistem

### Hardware minim recomandat
- **CPU:** 8 core-uri (Intel i7 / AMD Ryzen 7 sau superior)
- **RAM:** 16 GB minim, 32 GB recomandat
- **GPU:** NVIDIA cu 8 GB VRAM și suport CUDA 12.x (opțional dar recomandat)
- **Storage:** 50 GB spațiu liber (pentru cache VEP, model Llama, date temporare)

### Software necesar
- **Windows 11** cu **WSL2** (Ubuntu 24.04 recomandat) SAU Linux nativ
- **Java 21 LTS**
- **Python 3.10 sau 3.11**
- **Docker Desktop 27.0+**
- **NVIDIA CUDA 12.4** + drivere GPU (pentru inferență GPU)
- **Ollama 0.4.7+**
- **Maven 3.9+**

---

## Instalare și lansare

### 1. Clonare repository

```bash
git clone https://github.com/AndreiTica24/Licenta-BioAI.git
cd Licenta-BioAI
```

### 2. Configurare Modul Python AI

```bash
cd AI-Module

# Creare mediu conda
conda create -n ai-genomics python=3.11
conda activate ai-genomics

# Instalare dependințe
pip install -r requirements.txt

# Verificare GPU (opțional)
python -c "import torch; print(torch.cuda.is_available())"
```

**Descarcă modelul antrenat:**
Modelul `best_model.pth` (~15 MB) trebuie plasat în `AI-Module/checkpoints/`.

### 3. Configurare Backend Spring Boot

```bash
cd BackEnd/genomics-api
mvn clean install
```

**Configurare application.properties** (`src/main/resources/application.properties`):
```properties
# Adresa modulului Python
python.api.url=http://localhost:8000
python.api.key=[cheia API]

# Ollama
ollama.url=http://localhost:11434
ollama.model=llama3.1:8b

# JWT
jwt.secret=[secret HS384]
```

### 4. Configurare Docker VEP

```bash
# Pull imagine VEP oficială
docker pull ensemblorg/ensembl-vep:release_115.0

# Cache VEP + ClinVar (~ 15 GB, descărcare inițială o singură dată)
# Instrucțiuni detaliate în VEP/README.md
```

### 5. Configurare Ollama + Llama 3.1

```bash
# Pornire Ollama server
ollama serve

# Descărcare model (~ 4.7 GB)
ollama pull llama3.1:8b
```

---

## Rulare aplicație

Aplicația necesită pornirea a **trei servicii** în paralel:

### Terminal 1 — Modul Python AI (WSL/Linux)

```bash
cd AI-Module
conda activate ai-genomics
uvicorn predict_api:app --host 0.0.0.0 --port 8000
```

### Terminal 2 — Ollama (Windows/Linux)

```bash
ollama serve
```

### Terminal 3 — Backend Spring Boot (Windows/Linux)

```bash
cd BackEnd/genomics-api
mvn spring-boot:run
```

**Accesare aplicație:** http://localhost:8080

### Conturi demonstrative

- **Admin:** `admin@genomics.ro` / `Admin123!`
- **Utilizator:** `pacient@test.ro` / `Pacient123!`

---

## Testare rapidă

Pentru testare, se poate folosi un fișier BAM de dimensiuni reduse (chr22, ~130 MB) din setul GIAB HG002: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/

Procesarea pe chr22 durează aproximativ 2 minute și generează ~1.900 variante.

---

## Antrenare model (opțional)

Pentru re-antrenarea modelului CNN:

```bash
cd AI-Module
python train.py --config configs/train_config.yaml
```

Antrenarea durează aproximativ 5 ore pe RTX 5070 Laptop (30 epoci, HG003 train / HG004 validation).

---

## Tehnologii utilizate

| Componentă | Tehnologie |
|------------|------------|
| Model AI | PyTorch 2.5, CUDA 12.4 |
| API AI | FastAPI, Uvicorn |
| Backend | Spring Boot 3.3.5, Spring Security, JPA |
| Bază de date | H2 (embedded) |
| Anotare | Ensembl VEP 115, ClinVar 202502 |
| Asistent | Ollama, Llama 3.1 8B |
| Frontend | Thymeleaf, HTML5, JavaScript |
| Containerizare | Docker 27.0 |

---

## Licență

Cod dezvoltat în scop academic pentru lucrarea de licență. Nu este certificat ca dispozitiv medical.

---

## Contact

Andrei Tica — [email institutional]

