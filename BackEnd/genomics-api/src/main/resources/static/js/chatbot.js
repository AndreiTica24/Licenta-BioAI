/* ============================================================================
   chatbot.js — Floating chat widget for genetics assistant
   ============================================================================ */

let chatHistory = [];

/**
 * Initialize chatbot widget — appends HTML to body and sets up listeners.
 */
function initChatbot() {
    // Only initialize if user is logged in
    if (!getAuth()) return;

    const widgetHtml = `
        <div id="chatbot-widget">
            <button id="chatbot-toggle" onclick="toggleChatbot()">
                <span id="chatbot-toggle-icon">💬</span>
            </button>
            <div id="chatbot-panel" style="display:none;">
                <div id="chatbot-header">
                    <div>
                        <strong>🧬 Genetics Assistant</strong>
                        <div style="font-size:11px; color:rgba(255,255,255,0.85);">
                            Ask about variants, genes, ClinVar...
                        </div>
                    </div>
                    <button onclick="toggleChatbot()" 
                            style="background:none; border:none; color:white; 
                                   font-size:20px; cursor:pointer; padding:0 4px;">×</button>
                </div>
                <div id="chatbot-messages">
                    <div class="chat-message bot">
                        Hi! I'm your genetics assistant. I can answer questions about:
                        <ul style="margin:8px 0 0 16px; padding:0; font-size:13px;">
                            <li>Variant classification (Pathogenic, VUS, Benign)</li>
                            <li>Genes and chromosomes</li>
                            <li>ClinVar, gnomAD, VEP</li>
                            <li>Bioinformatics terminology</li>
                        </ul>
                    </div>
                </div>
                <div id="chatbot-input-container">
                    <input type="text" id="chatbot-input" 
                           placeholder="Ask about a variant, gene, or term..." 
                           onkeydown="if(event.key==='Enter') sendChatMessage()">
                    <button onclick="sendChatMessage()" id="chatbot-send-btn">Send</button>
                </div>
            </div>
        </div>
    `;

    const container = document.createElement('div');
    container.innerHTML = widgetHtml;
    document.body.appendChild(container);
}

function toggleChatbot() {
    const panel = document.getElementById('chatbot-panel');
    const toggleIcon = document.getElementById('chatbot-toggle-icon');
    if (panel.style.display === 'none') {
        panel.style.display = 'flex';
        toggleIcon.textContent = '×';
        document.getElementById('chatbot-input').focus();
    } else {
        panel.style.display = 'none';
        toggleIcon.textContent = '💬';
    }
}

async function sendChatMessage() {
    const input = document.getElementById('chatbot-input');
    const sendBtn = document.getElementById('chatbot-send-btn');
    const question = input.value.trim();

    if (!question) return;

    appendChatMessage('user', question);
    input.value = '';
    input.disabled = true;
    sendBtn.disabled = true;

    // Show "typing" indicator
    const typingId = appendChatMessage('bot', '<span class="chat-typing">●●●</span>', true);

    try {
        const response = await apiPost('/api/chatbot/ask', { question });

        // Remove typing indicator
        const typingEl = document.getElementById(typingId);
        if (typingEl) typingEl.remove();

        // Format answer (basic markdown-like)
        const formattedAnswer = formatChatAnswer(response.answer);
        appendChatMessage('bot', formattedAnswer);

    } catch (e) {
        const typingEl = document.getElementById(typingId);
        if (typingEl) typingEl.remove();
        appendChatMessage('bot', '⚠️ Sorry, an error occurred: ' + e.message);
    }

    input.disabled = false;
    sendBtn.disabled = false;
    input.focus();
}

function appendChatMessage(role, content, isTyping = false) {
    const messagesDiv = document.getElementById('chatbot-messages');
    const messageDiv = document.createElement('div');
    const messageId = 'chat-msg-' + Date.now() + '-' + Math.random().toString(36).substring(7);
    messageDiv.id = messageId;
    messageDiv.className = 'chat-message ' + role;
    messageDiv.innerHTML = content;
    messagesDiv.appendChild(messageDiv);
    messagesDiv.scrollTop = messagesDiv.scrollHeight;
    return messageId;
}

function formatChatAnswer(text) {
    // Convert line breaks to <br>
    text = text.replace(/\n\n/g, '</p><p>').replace(/\n/g, '<br>');
    // Wrap in paragraph
    text = '<p>' + text + '</p>';
    // Bold **text**
    text = text.replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>');
    // Italic *text*
    text = text.replace(/(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)/g, '<em>$1</em>');
    return text;
}

// Auto-initialize when DOM is ready
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initChatbot);
} else {
    initChatbot();
}