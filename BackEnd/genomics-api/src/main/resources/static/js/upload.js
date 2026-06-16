/* ============================================================================
   upload.js — BAM upload + job polling
   ============================================================================ */

let currentJobId = null;
let pollingInterval = null;

/**
 * Initialize upload area with click + drag-and-drop.
 */
function initUploadArea() {
    const uploadArea = document.getElementById('upload-area');
    const fileInput = document.getElementById('file-input');

    if (!uploadArea || !fileInput) return;

    uploadArea.addEventListener('click', () => fileInput.click());

    uploadArea.addEventListener('dragover', (e) => {
        e.preventDefault();
        uploadArea.classList.add('dragover');
    });

    uploadArea.addEventListener('dragleave', () => {
        uploadArea.classList.remove('dragover');
    });

    uploadArea.addEventListener('drop', (e) => {
        e.preventDefault();
        uploadArea.classList.remove('dragover');
        if (e.dataTransfer.files.length > 0) {
            fileInput.files = e.dataTransfer.files;
            handleFileSelection();
        }
    });

    fileInput.addEventListener('change', handleFileSelection);
}

/**
 * Display selected file info.
 */
function handleFileSelection() {
    const fileInput = document.getElementById('file-input');
    const fileInfo = document.getElementById('file-info');
    const uploadButton = document.getElementById('upload-button');

    if (fileInput.files.length === 0) {
        fileInfo.style.display = 'none';
        uploadButton.disabled = true;
        return;
    }

    const file = fileInput.files[0];

    if (!file.name.toLowerCase().endsWith('.bam')) {
        alert('Please select a BAM file (.bam extension)');
        fileInput.value = '';
        return;
    }

    document.getElementById('file-name').textContent = file.name;
    document.getElementById('file-size').textContent = formatFileSize(file.size);
    fileInfo.style.display = 'block';
    uploadButton.disabled = false;
}

/**
 * Handle upload form submission.
 */
async function handleUploadSubmit(event) {
    event.preventDefault();

    const fileInput = document.getElementById('file-input');
    const sampleName = document.getElementById('sample-name').value || 'sample';
    const confidence = document.getElementById('confidence').value || 0.7;
    const uploadButton = document.getElementById('upload-button');
    const progressSection = document.getElementById('progress-section');
    const progressBar = document.getElementById('progress-bar');
    const progressText = document.getElementById('progress-text');
    const errorBox = document.getElementById('upload-error');

    if (fileInput.files.length === 0) {
        return;
    }

    const file = fileInput.files[0];

    uploadButton.disabled = true;
    progressSection.style.display = 'block';
    errorBox.style.display = 'none';
    progressBar.style.width = '0%';
    progressText.textContent = 'Uploading...';

    try {
        const result = await uploadBam(file, sampleName, confidence, (percent) => {
            progressBar.style.width = percent + '%';
            progressText.textContent = `Uploading... ${percent.toFixed(0)}%`;
        });

        currentJobId = result.job_id;
        progressBar.style.width = '100%';
        progressText.textContent = `Upload complete! Job ${currentJobId.substring(0, 8)}...`;

        // Start polling for status
        startPollingStatus(currentJobId);
        await loadJobs();  // Refresh job list

    } catch (e) {
        errorBox.textContent = e.message;
        errorBox.style.display = 'block';
        uploadButton.disabled = false;
        progressSection.style.display = 'none';
    }
}

/**
 * Poll job status every 5 seconds until completion.
 * Also refreshes the jobs table.
 */
function startPollingStatus(jobId) {
    if (pollingInterval) {
        clearInterval(pollingInterval);
    }

    pollingInterval = setInterval(async () => {
        try {
            const status = await apiGet(`/api/variants/status/${jobId}`);
            updateProgressDisplay(status);

            // Refresh jobs list to sync UI with DB
            await loadJobs();

            if (status.status === 'completed' || status.status === 'failed') {
                clearInterval(pollingInterval);
                pollingInterval = null;

                if (status.status === 'completed') {
                    document.getElementById('progress-text').innerHTML =
                        `✅ Analysis complete! ${status.n_variants} variants detected. ` +
                        `<a href="/result/${jobId}">View Report →</a>`;
                } else {
                    document.getElementById('progress-text').innerHTML =
                        `❌ Analysis failed: ${status.error || 'Unknown error'}`;
                }
                document.getElementById('upload-button').disabled = false;
            }
        } catch (e) {
            console.error('Polling error:', e);
        }
    }, 5000);
}

/**
 * Update progress display based on job status.
 */
function updateProgressDisplay(status) {
    const progressText = document.getElementById('progress-text');
    const progressBar = document.getElementById('progress-bar');

    if (status.status === 'running') {
        progressBar.classList.add('indeterminate');
        progressText.textContent = `🔬 ${status.progress || 'Analyzing...'}`;
    } else if (status.status === 'completed') {
        progressBar.classList.remove('indeterminate');
        progressBar.style.width = '100%';
    } else if (status.status === 'failed') {
        progressBar.classList.remove('indeterminate');
    }
}

/**
 * Load all jobs from backend DB (persistent across restarts).
 */
async function loadJobs() {
    const jobsTable = document.getElementById('jobs-table-body');
    if (!jobsTable) return;

    try {
        const jobs = await apiGet('/api/jobs/my');

        if (jobs.length === 0) {
            jobsTable.innerHTML = '<tr><td colspan="6" style="text-align:center; padding:24px; color:var(--text-muted);">No analyses yet. Upload a BAM file to start.</td></tr>';
            return;
        }

        const rows = jobs.map(job => renderJobRow(job));
        jobsTable.innerHTML = rows.join('');
    } catch (e) {
        console.error('Error loading jobs:', e);
        jobsTable.innerHTML = '<tr><td colspan="6" style="text-align:center; padding:24px; color:var(--danger);">Error loading jobs: ' + e.message + '</td></tr>';
    }
}

/**
 * Render one job row in the table.
 */
function renderJobRow(job) {
    const statusBadge = `<span class="badge badge-status-${job.status.toLowerCase()}">${job.status}</span>`;
    const variants = job.nVariants ? `${job.nVariants.toLocaleString()}` : '-';
    const filename = job.bamFilename || '-';
    let actions = '-';

    if (job.status === 'COMPLETED') {
        actions = `<a href="/result/${job.id}" class="btn btn-primary" style="padding:6px 12px; font-size:13px;">📊 Analyze</a>
                   <button onclick="deleteJob('${job.id}')" class="btn btn-secondary" style="padding:6px 12px; font-size:13px; margin-left:6px;">🗑️</button>`;
    } else if (job.status === 'PROCESSING') {
        actions = `<span style="color:var(--info);"><span class="loader" style="width:14px;height:14px;"></span> Processing...</span>`;
    } else {
        actions = `<button onclick="deleteJob('${job.id}')" class="btn btn-secondary" style="padding:6px 12px; font-size:13px;">🗑️</button>`;
    }

    return `
        <tr>
            <td><code style="font-size:12px;">${job.id.substring(0, 8)}...</code></td>
            <td style="font-size:13px;">${filename}</td>
            <td>${statusBadge}</td>
            <td>${formatDate(job.createdAt)}</td>
            <td>${variants}</td>
            <td class="actions">${actions}</td>
        </tr>
    `;
}

async function deleteJob(jobId) {
    if (!confirm('Delete this analysis from history? This cannot be undone.')) {
        return;
    }
    try {
        await apiDelete('/api/jobs/' + jobId);
        await loadJobs();
    } catch (e) {
        alert('Error: ' + e.message);
    }
}