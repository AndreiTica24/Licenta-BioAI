/* ============================================================================
   api.js — REST API client with JWT auto-injection
   ============================================================================ */

/**
 * Generic API call with JWT in Authorization header.
 * Auto-logout on 401/403.
 */
async function apiCall(path, options = {}) {
    const auth = getAuth();
    if (!auth) {
        window.location.href = '/login';
        throw new Error('Not authenticated');
    }

    const headers = options.headers || {};
    headers['Authorization'] = `Bearer ${auth.token}`;

    if (options.body && !(options.body instanceof FormData)) {
        headers['Content-Type'] = 'application/json';
    }

    const response = await fetch(path, {
        ...options,
        headers
    });

    if (response.status === 401 || response.status === 403) {
        // Token expired or invalid → logout
        logout();
        throw new Error('Session expired');
    }

    return response;
}

/**
 * GET request that returns JSON.
 */
async function apiGet(path) {
    const response = await apiCall(path);
    if (!response.ok) {
        const error = await response.json().catch(() => ({ error: 'Request failed' }));
        throw new Error(error.error || `HTTP ${response.status}`);
    }
    return response.json();
}

/**
 * POST request with JSON body.
 */
async function apiPost(path, body) {
    const response = await apiCall(path, {
        method: 'POST',
        body: body ? JSON.stringify(body) : null
    });
    if (!response.ok) {
        const error = await response.json().catch(() => ({ error: 'Request failed' }));
        throw new Error(error.error || `HTTP ${response.status}`);
    }
    return response.json();
}

/**
 * DELETE request.
 */
async function apiDelete(path) {
    const response = await apiCall(path, { method: 'DELETE' });
    if (!response.ok) {
        const error = await response.json().catch(() => ({ error: 'Request failed' }));
        throw new Error(error.error || `HTTP ${response.status}`);
    }
    return response.json();
}

/**
 * Upload BAM file with progress tracking.
 * Uses XMLHttpRequest to support progress events.
 */
function uploadBam(file, sampleName, confidence, onProgress) {
    return new Promise((resolve, reject) => {
        const auth = getAuth();
        if (!auth) {
            window.location.href = '/login';
            return;
        }

        const formData = new FormData();
        formData.append('file', file);
        formData.append('sampleName', sampleName);
        formData.append('confidence', confidence);

        const xhr = new XMLHttpRequest();

        xhr.upload.addEventListener('progress', (e) => {
            if (e.lengthComputable && onProgress) {
                const percent = (e.loaded / e.total) * 100;
                onProgress(percent, e.loaded, e.total);
            }
        });

        xhr.addEventListener('load', () => {
            if (xhr.status === 200) {
                try {
                    resolve(JSON.parse(xhr.responseText));
                } catch (e) {
                    reject(new Error('Invalid response'));
                }
            } else if (xhr.status === 401 || xhr.status === 403) {
                logout();
                reject(new Error('Session expired'));
            } else {
                try {
                    const err = JSON.parse(xhr.responseText);
                    reject(new Error(err.error || `HTTP ${xhr.status}`));
                } catch (e) {
                    reject(new Error(`HTTP ${xhr.status}`));
                }
            }
        });

        xhr.addEventListener('error', () => reject(new Error('Network error')));
        xhr.addEventListener('abort', () => reject(new Error('Upload cancelled')));

        xhr.open('POST', '/api/variants/upload');
        xhr.setRequestHeader('Authorization', `Bearer ${auth.token}`);
        xhr.send(formData);
    });
}

/**
 * Format file size in human-readable form.
 */
function formatFileSize(bytes) {
    if (bytes < 1024) return bytes + ' B';
    if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(1) + ' KB';
    if (bytes < 1024 * 1024 * 1024) return (bytes / 1024 / 1024).toFixed(1) + ' MB';
    return (bytes / 1024 / 1024 / 1024).toFixed(2) + ' GB';
}

/**
 * Format date in human-readable form.
 */
function formatDate(isoString) {
    if (!isoString) return '-';
    const d = new Date(isoString);
    return d.toLocaleString('en-US', {
        year: 'numeric', month: 'short', day: 'numeric',
        hour: '2-digit', minute: '2-digit'
    });
}