"""
predict_api.py — REST API pentru pipeline-ul de variant calling
============================================================
FastAPI service care expune funcționalitatea scan_bam.py + predict.py
ca endpoint HTTP, consumabil de backend-ul Spring Boot.

ARHITECTURA:
  Spring Boot (Java) ──HTTPS+JWT──► predict_api.py (Python) ──► CNN model
                                          │
                                          └──► scan_bam + predict
                                          
SECURITATE:
  - API key validation (X-API-Key header)
  - Job IDs unice pentru fiecare cerere (UUID)
  - Fișiere temporare șterse automat după procesare
  - Logging request-uri pentru audit

ENDPOINTS:
  GET  /health           — verificare stare service
  POST /predict          — pipeline complet BAM → VCF + JSON
  GET  /jobs/{job_id}    — status job în procesare
  GET  /jobs/{job_id}/result — descarcă rezultat finalizat

Rulare:
    uvicorn predict_api:app --host 0.0.0.0 --port 8000
    
    # Cu HTTPS (auto-cert pentru dezvoltare):
    uvicorn predict_api:app --host 0.0.0.0 --port 8443 \\
        --ssl-keyfile=certs/key.pem --ssl-certfile=certs/cert.pem
"""

import asyncio
import logging
import os
import secrets
import shutil
import time
import uuid
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict

from fastapi import FastAPI, HTTPException, Header, BackgroundTasks, UploadFile, File
from fastapi.responses import JSONResponse, FileResponse
from pydantic import BaseModel

import scan_bam
import predict as predict_module


# ============================================================================
# Configurare
# ============================================================================

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("predict_api.log"),
        logging.StreamHandler(),
    ],
)
logger = logging.getLogger(__name__)

# Director pentru fișiere temporare per job
JOBS_DIR = Path("api_jobs")
JOBS_DIR.mkdir(parents=True, exist_ok=True)

# Cale către modelul antrenat
MODEL_PATH = os.environ.get(
    "MODEL_PATH", "checkpoints_unrelated/best_model.pth"
)

# API key — în producție vine din env var, în dev citim/generăm
API_KEY_FILE = Path(".api_key")
if "API_KEY" in os.environ:
    API_KEY = os.environ["API_KEY"]
elif API_KEY_FILE.exists():
    API_KEY = API_KEY_FILE.read_text().strip()
else:
    API_KEY = secrets.token_urlsafe(32)
    API_KEY_FILE.write_text(API_KEY)
    logger.info(f"📔 API key generat și salvat în {API_KEY_FILE}: {API_KEY}")

# Jobs in-memory (production: Redis/database)
jobs_state: Dict[str, Dict] = {}


# ============================================================================
# Modele Pydantic (request/response schema)
# ============================================================================

class PredictRequest(BaseModel):
    bam_path:    str               # cale absolută pe filesystem-ul WSL
    sample_name: Optional[str] = "sample"
    threads:     int   = 4
    confidence:  float = 0.7


class JobStatus(BaseModel):
    job_id:      str
    status:      str               # pending | running | completed | failed
    progress:    Optional[str] = None
    created_at:  str
    completed_at: Optional[str] = None
    n_candidates: Optional[int] = None
    n_variants:   Optional[int] = None
    error:        Optional[str] = None


class HealthResponse(BaseModel):
    status:      str
    model_loaded: bool
    gpu_available: bool
    model_path:   str
    version:      str


# ============================================================================
# Autentificare
# ============================================================================

def verify_api_key(x_api_key: Optional[str] = Header(None)) -> bool:
    """Verifică API key trimisă în header X-API-Key."""
    if not x_api_key or not secrets.compare_digest(x_api_key, API_KEY):
        raise HTTPException(status_code=401, detail="Invalid or missing API key")
    return True


# ============================================================================
# Pipeline procesare (rulează asincron în background)
# ============================================================================

def process_bam(job_id: str, bam_path: str, sample_name: str,
                threads: int, confidence: float):
    """Pipeline complet scan_bam + predict, executat asincron."""
    job = jobs_state[job_id]
    job_dir = JOBS_DIR / job_id

    try:
        job["status"] = "running"
        job["progress"] = "Verificare BAM..."
        logger.info(f"[{job_id}] Pornire procesare: {bam_path}")

        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM nu există: {bam_path}")

        # Verificăm/generăm index .bai dacă lipsește
        bai_path = bam_path + ".bai"
        bai_alt  = bam_path.replace(".bam", ".bai")
        if not os.path.exists(bai_path) and not os.path.exists(bai_alt):
            job["progress"] = "Generare index BAM (.bai)..."
            logger.info(f"[{job_id}] Index lipsă, generez .bai...")
            import pysam
            pysam.index(bam_path)
            logger.info(f"[{job_id}] Index .bai generat")

        # PASUL 1: scan_bam
        job["progress"] = "Pre-filtrare BAM (~15 min)..."
        logger.info(f"[{job_id}] Pre-filtrare...")

        candidates_tsv = str(job_dir / "candidates.tsv")
        t0 = time.time()
        scan_stats = scan_bam.scan_bam(
            bam_path    = bam_path,
            output_path = candidates_tsv,
            fasta_path  = scan_bam.HUMAN_CHROMOSOMES and "data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta",
            threads     = threads,
        )
        scan_time = time.time() - t0
        job["n_candidates"] = scan_stats["n_candidates"]
        logger.info(f"[{job_id}] Scan: {scan_stats['n_candidates']:,} candidați "
                    f"în {scan_time:.0f}s")

        # PASUL 2: predict
        job["progress"] = "Clasificare CNN (~15 min)..."
        logger.info(f"[{job_id}] Inferență CNN...")

        output_vcf  = str(job_dir / f"{sample_name}.vcf")
        output_json = str(job_dir / f"{sample_name}.json")

        t0 = time.time()
        predict_stats = predict_module.predict(
            candidates_tsv = candidates_tsv,
            bam_path       = bam_path,
            model_path     = MODEL_PATH,
            output_vcf     = output_vcf,
            output_json    = output_json,
            batch_size     = 512,
            num_workers    = 8,
            confidence     = confidence,
        )
        predict_time = time.time() - t0

        # PASUL 3: finalizare
        job["status"]       = "completed"
        job["completed_at"] = datetime.now().isoformat()
        job["n_variants"]   = predict_stats["n_variants_vcf"]
        job["n_het"]        = predict_stats["n_het"]
        job["n_hom_alt"]    = predict_stats["n_hom_alt"]
        job["output_vcf"]   = output_vcf
        job["output_json"]  = output_json
        job["scan_time_s"]   = round(scan_time, 1)
        job["predict_time_s"] = round(predict_time, 1)
        job["progress"]     = "Completat"

        logger.info(f"[{job_id}] ✅ Completat: {predict_stats['n_variants_vcf']:,} variante "
                    f"(scan {scan_time:.0f}s + predict {predict_time:.0f}s)")

        # Ștergem candidates.tsv (mare, nu mai e nevoie)
        try:
            os.remove(candidates_tsv)
        except OSError:
            pass

    except Exception as e:
        logger.error(f"[{job_id}] ❌ Eroare: {e}", exc_info=True)
        job["status"]       = "failed"
        job["error"]        = str(e)
        job["completed_at"] = datetime.now().isoformat()


# ============================================================================
# FastAPI app
# ============================================================================

app = FastAPI(
    title="VariantCallerCNN1D API",
    description="REST API pentru variant calling cu CNN 1D",
    version="1.0.0",
)


@app.get("/health", response_model=HealthResponse)
async def health():
    """Endpoint pentru verificare stare. Nu necesită autentificare."""
    import torch
    model_exists = os.path.exists(MODEL_PATH)
    return HealthResponse(
        status         = "healthy" if model_exists else "model_missing",
        model_loaded   = model_exists,
        gpu_available  = torch.cuda.is_available(),
        model_path     = MODEL_PATH,
        version        = "1.0.0",
    )


@app.post("/predict")
async def start_predict(
    request: PredictRequest,
    background_tasks: BackgroundTasks,
    _auth: bool = None,  # validat manual mai jos
    x_api_key: Optional[str] = Header(None),
):
    """
    Pornește pipeline-ul de variant calling pe un BAM.
    Returnează job_id pentru polling ulterior.
    """
    verify_api_key(x_api_key)

    # Verificări input
    if not os.path.exists(request.bam_path):
        raise HTTPException(status_code=404,
                            detail=f"BAM nu există: {request.bam_path}")

    # Creăm job
    job_id = str(uuid.uuid4())
    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    jobs_state[job_id] = {
        "job_id":       job_id,
        "status":       "pending",
        "progress":     "Așteptare procesare...",
        "created_at":   datetime.now().isoformat(),
        "completed_at": None,
        "bam_path":     request.bam_path,
        "sample_name":  request.sample_name,
    }

    # Pornim procesarea în background
    background_tasks.add_task(
        process_bam,
        job_id      = job_id,
        bam_path    = request.bam_path,
        sample_name = request.sample_name or "sample",
        threads     = request.threads,
        confidence  = request.confidence,
    )

    logger.info(f"[{job_id}] Job creat pentru {request.bam_path}")
    return {"job_id": job_id, "status": "pending", "message": "Job pornit"}


@app.get("/jobs/{job_id}", response_model=JobStatus)
async def get_job_status(
    job_id: str,
    x_api_key: Optional[str] = Header(None),
):
    """Returnează statusul unui job."""
    verify_api_key(x_api_key)

    if job_id not in jobs_state:
        raise HTTPException(status_code=404, detail="Job not found")

    job = jobs_state[job_id]
    return JobStatus(
        job_id       = job["job_id"],
        status       = job["status"],
        progress     = job.get("progress"),
        created_at   = job["created_at"],
        completed_at = job.get("completed_at"),
        n_candidates = job.get("n_candidates"),
        n_variants   = job.get("n_variants"),
        error        = job.get("error"),
    )


@app.get("/jobs/{job_id}/result")
async def get_job_result(
    job_id: str,
    format: str = "json",  # "json" | "vcf"
    x_api_key: Optional[str] = Header(None),
):
    """Returnează rezultatul (VCF sau JSON) după job completat."""
    verify_api_key(x_api_key)

    if job_id not in jobs_state:
        raise HTTPException(status_code=404, detail="Job not found")

    job = jobs_state[job_id]
    if job["status"] != "completed":
        raise HTTPException(
            status_code=400,
            detail=f"Job e în status '{job['status']}', nu 'completed'"
        )

    if format == "vcf":
        file_path = job.get("output_vcf")
    else:
        file_path = job.get("output_json")

    if not file_path or not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="Output file not found")

    return FileResponse(
        path     = file_path,
        media_type = "application/octet-stream",
        filename = os.path.basename(file_path),
    )


@app.delete("/jobs/{job_id}")
async def delete_job(
    job_id: str,
    x_api_key: Optional[str] = Header(None),
):
    """Șterge un job și fișierele aferente."""
    verify_api_key(x_api_key)

    if job_id not in jobs_state:
        raise HTTPException(status_code=404, detail="Job not found")

    # Ștergem directorul jobului
    job_dir = JOBS_DIR / job_id
    if job_dir.exists():
        shutil.rmtree(job_dir, ignore_errors=True)

    del jobs_state[job_id]
    logger.info(f"[{job_id}] Job șters")
    return {"message": "Job șters cu succes"}


@app.get("/jobs")
async def list_jobs(x_api_key: Optional[str] = Header(None)):
    """Listează toate joburile (cu detaliile lor)."""
    verify_api_key(x_api_key)
    return {
        "n_jobs": len(jobs_state),
        "jobs":   list(jobs_state.values()),
    }


# ============================================================================
# Startup banner
# ============================================================================

@app.on_event("startup")
async def startup_event():
    import torch
    print("\n" + "=" * 70)
    print("🚀 VariantCallerCNN1D API — Pornit")
    print("=" * 70)
    print(f"   Model              : {MODEL_PATH}")
    print(f"   Model exists       : {os.path.exists(MODEL_PATH)}")
    print(f"   GPU available      : {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"   GPU device         : {torch.cuda.get_device_name(0)}")
    print(f"   Jobs directory     : {JOBS_DIR}")
    print(f"   API Key            : {API_KEY[:8]}...{API_KEY[-4:]} "
          f"(complet în {API_KEY_FILE})")
    print(f"   Docs               : http://localhost:8000/docs")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)