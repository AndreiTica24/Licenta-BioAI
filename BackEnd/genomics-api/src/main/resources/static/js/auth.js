/* ============================================================================
   auth.js — Login, register, JWT management
   ============================================================================ */

const AUTH_KEY = 'genomics_auth';

/**
 * Get current logged-in user info from localStorage.
 * Returns null if not logged in or token expired.
 */
function getAuth() {
    const raw = localStorage.getItem(AUTH_KEY);
    if (!raw) return null;
    try {
        const auth = JSON.parse(raw);
        if (Date.now() > auth.expiresAt) {
            localStorage.removeItem(AUTH_KEY);
            return null;
        }
        return auth;
    } catch (e) {
        localStorage.removeItem(AUTH_KEY);
        return null;
    }
}

/**
 * Store auth response in localStorage.
 */
function setAuth(authResponse) {
    const data = {
        token: authResponse.token,
        email: authResponse.email,
        fullName: authResponse.fullName,
        role: authResponse.role,
        expiresAt: Date.now() + authResponse.expiresInMs
    };
    localStorage.setItem(AUTH_KEY, JSON.stringify(data));
}

/**
 * Logout — clear token and redirect to login.
 */
function logout() {
    localStorage.removeItem(AUTH_KEY);
    window.location.href = '/login';
}

/**
 * Require auth — redirect to /login if not authenticated.
 * Call this at the top of protected pages.
 */
function requireAuth() {
    const auth = getAuth();
    if (!auth) {
        window.location.href = '/login';
        return null;
    }
    return auth;
}

/**
 * Require admin role — redirect to dashboard if not admin.
 */
function requireAdmin() {
    const auth = requireAuth();
    if (!auth) return null;
    if (auth.role !== 'ADMIN') {
        window.location.href = '/dashboard';
        return null;
    }
    return auth;
}

/**
 * Login form handler.
 */
async function handleLogin(event) {
    event.preventDefault();
    const form = event.target;
    const email = form.email.value.trim();
    const password = form.password.value;
    const errorBox = document.getElementById('login-error');
    const submitBtn = form.querySelector('button[type="submit"]');

    errorBox.style.display = 'none';
    submitBtn.disabled = true;
    submitBtn.innerHTML = '<span class="loader"></span> Logging in...';

    try {
        const response = await fetch('/api/auth/login', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ email, password })
        });

        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Login failed');
        }

        setAuth(data);
        window.location.href = '/dashboard';

    } catch (e) {
        errorBox.textContent = e.message;
        errorBox.style.display = 'block';
        submitBtn.disabled = false;
        submitBtn.innerHTML = 'Login';
    }
}

/**
 * Register form handler.
 */
async function handleRegister(event) {
    event.preventDefault();
    const form = event.target;
    const fullName = form.fullName.value.trim();
    const email = form.email.value.trim();
    const password = form.password.value;
    const errorBox = document.getElementById('register-error');
    const submitBtn = form.querySelector('button[type="submit"]');

    errorBox.style.display = 'none';
    submitBtn.disabled = true;
    submitBtn.innerHTML = '<span class="loader"></span> Creating account...';

    try {
        const response = await fetch('/api/auth/register', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ fullName, email, password })
        });

        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Registration failed');
        }

        setAuth(data);
        window.location.href = '/dashboard';

    } catch (e) {
        errorBox.textContent = e.message;
        errorBox.style.display = 'block';
        submitBtn.disabled = false;
        submitBtn.innerHTML = 'Create Account';
    }
}

/**
 * Switch between login/register tabs.
 */
function showTab(tabName) {
    document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
    document.querySelectorAll('.tab-content').forEach(c => c.style.display = 'none');
    document.querySelector(`[data-tab="${tabName}"]`).classList.add('active');
    document.getElementById(`${tabName}-content`).style.display = 'block';
}

/**
 * Render user info in navbar (used on protected pages).
 */
function renderUserInfo() {
    const auth = getAuth();
    if (!auth) return;

    const userInfo = document.getElementById('user-info');
    if (userInfo) {
        const roleBadge = auth.role === 'ADMIN'
            ? '<span class="role-badge admin">Admin</span>'
            : '<span class="role-badge">User</span>';

        userInfo.innerHTML = `
            ${roleBadge}
            <span>${auth.fullName}</span>
            <a href="#" onclick="logout(); return false;" style="color: var(--danger); font-weight: 500;">Logout</a>
        `;
    }

    // Show admin link if user is admin
    const adminLink = document.getElementById('admin-link');
    if (adminLink && auth.role === 'ADMIN') {
        adminLink.style.display = 'inline-block';
    }
}