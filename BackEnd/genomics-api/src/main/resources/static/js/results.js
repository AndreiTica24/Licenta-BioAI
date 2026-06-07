/* ============================================================================
   results.js — Display medical report with VEP annotations
   ============================================================================ */

/**
 * Load and display annotation results for a job.
 */
async function loadResults(jobId) {
    const loadingDiv = document.getElementById('loading');
    const reportDiv = document.getElementById('report');
    const errorBox = document.getElementById('result-error');

    loadingDiv.style.display = 'block';
    reportDiv.style.display = 'none';
    errorBox.style.display = 'none';

    try {
        // First check job status
        const status = await apiGet(`/api/variants/status/${jobId}`);

        if (status.status !== 'completed') {
            throw new Error(`Job is in state "${status.status}". Please wait for completion.`);
        }

        // Trigger VEP annotation
        document.getElementById('loading-text').textContent =
            '🧬 Running clinical annotation (ClinVar + VEP)... This takes ~30 seconds.';

        const result = await apiPost(`/api/variants/annotate/${jobId}`);

        renderReport(result, jobId);
        loadingDiv.style.display = 'none';
        reportDiv.style.display = 'block';

    } catch (e) {
        loadingDiv.style.display = 'none';
        errorBox.textContent = e.message;
        errorBox.style.display = 'block';
    }
}

/**
 * Render the medical report from VEP annotation result.
 */
function renderReport(data, jobId) {
    // Statistics cards
    document.getElementById('stat-total').textContent = data.n_variants.toLocaleString();
    document.getElementById('stat-clinvar').textContent = data.with_clinvar.toLocaleString();
    document.getElementById('stat-genes').textContent = data.with_gene.toLocaleString();
    document.getElementById('stat-time').textContent = data.annotation_time_s + 's';

    // Classification distribution
    const dist = data.by_classification || {};
    renderClassificationChart(dist);

    // Pathogenic / Likely pathogenic variants (most important)
    const importantVariants = data.variants.filter(v =>
        v.finalClassification === 'PATHOGENIC' ||
        v.finalClassification === 'LIKELY_PATHOGENIC'
    );
    renderVariantsTable('important-variants', importantVariants, 50);

    // VUS variants (uncertain)
    const vusVariants = data.variants.filter(v =>
        v.finalClassification === 'VUS'
    );
    renderVariantsTable('vus-variants', vusVariants, 50);

    // All gene-affecting variants summary
    const allWithGene = data.variants.filter(v =>
        v.geneSymbol && v.geneSymbol.length > 0
    );
    renderVariantsTable('all-variants', allWithGene, 100);

    // Update counters
    document.getElementById('count-important').textContent = importantVariants.length;
    document.getElementById('count-vus').textContent = vusVariants.length;
    document.getElementById('count-all').textContent = allWithGene.length;

    // Download links
    document.getElementById('download-vcf').href = `/api/variants/vcf/${jobId}`;
}

/**
 * Render classification distribution as bar chart.
 */
function renderClassificationChart(dist) {
    const container = document.getElementById('classification-chart');
    const order = ['PATHOGENIC', 'LIKELY_PATHOGENIC', 'VUS', 'LIKELY_BENIGN', 'BENIGN', 'UNKNOWN'];
    const colors = {
        'PATHOGENIC': 'var(--pathogenic)',
        'LIKELY_PATHOGENIC': 'var(--likely-pathogenic)',
        'VUS': 'var(--vus)',
        'LIKELY_BENIGN': 'var(--likely-benign)',
        'BENIGN': 'var(--benign)',
        'UNKNOWN': 'var(--text-muted)'
    };

    const total = Object.values(dist).reduce((a, b) => a + b, 0) || 1;

    const html = order.map(category => {
        const count = dist[category] || 0;
        const percent = ((count / total) * 100).toFixed(1);
        const display = category.replace(/_/g, ' ');
        return `
            <div style="margin-bottom: 14px;">
                <div style="display: flex; justify-content: space-between; margin-bottom: 4px;">
                    <span style="font-weight: 500; font-size: 14px;">${display}</span>
                    <span style="color: var(--text-muted); font-size: 13px;">
                        <strong>${count.toLocaleString()}</strong> (${percent}%)
                    </span>
                </div>
                <div style="width: 100%; height: 12px; background: var(--bg); border-radius: 6px; overflow: hidden;">
                    <div style="width: ${percent}%; height: 100%; background: ${colors[category]}; border-radius: 6px;"></div>
                </div>
            </div>
        `;
    }).join('');

    container.innerHTML = html;
}

/**
 * Render variants in a table.
 */
function renderVariantsTable(tableId, variants, maxRows) {
    const tbody = document.getElementById(tableId);
    if (!tbody) return;

    if (variants.length === 0) {
        tbody.innerHTML = '<tr><td colspan="6" style="text-align:center; padding:24px; color:var(--text-muted);">No variants found in this category.</td></tr>';
        return;
    }

    const limit = Math.min(maxRows, variants.length);
    const rows = variants.slice(0, limit).map(v => {
        const classBadge = v.finalClassification
            ? `<span class="badge badge-${v.finalClassification.toLowerCase().replace(/_/g, '-')}">${v.finalClassification.replace(/_/g, ' ')}</span>`
            : '-';
        const clinvarInfo = v.clinSig
            ? `<small style="color:var(--text-muted);">${v.clinSig.substring(0, 50)}</small>`
            : '';
        const disease = v.clinDisease
            ? v.clinDisease.replace(/&/g, ' • ').substring(0, 80)
            : '-';

        return `
            <tr>
                <td><code style="font-size: 12px;">${v.chrom}:${v.pos}</code><br><small>${v.ref}→${v.alt}</small></td>
                <td><strong>${v.geneSymbol || '-'}</strong></td>
                <td><small>${v.consequence ? v.consequence.replace(/_/g, ' ') : '-'}</small></td>
                <td>${classBadge}<br>${clinvarInfo}</td>
                <td><small>${disease}</small></td>
                <td>${(v.confidence * 100).toFixed(1)}%</td>
            </tr>
        `;
    }).join('');

    let html = rows;
    if (variants.length > limit) {
        html += `<tr><td colspan="6" style="text-align:center; padding:12px; color:var(--text-muted);">
            Showing ${limit} of ${variants.length.toLocaleString()} total
        </td></tr>`;
    }

    tbody.innerHTML = html;
}