/* fz docs chat assistant - vanilla JS, no framework. ES module.
 *
 * Modes:
 *   local  - WebLLM (https://esm.run/@mlc-ai/web-llm), runs in browser via WebGPU.
 *   groq   - BYOK call to api.groq.com (OpenAI-compatible /chat/completions).
 *   gemini - BYOK call to generativelanguage.googleapis.com.
 *
 * BYOK keys are stored in localStorage only and only sent to the chosen
 * provider's API endpoint - never to the docs origin.
 */

// Resolve asset paths relative to this script so the widget works whether
// the host page lives at /, /modules/, or any other subdirectory.
const SCRIPT_URL = new URL(import.meta.url);
const CHAT_ROOT = new URL("./", SCRIPT_URL);
const SYSTEM_PROMPT_URL = new URL("system-prompt.txt", CHAT_ROOT).href;

const LS_KEYS = {
  provider: "fz-chat-provider",
  groqKey: "fz-chat-groq-key",
  geminiKey: "fz-chat-gemini-key",
  openrouterKey: "fz-chat-openrouter-key",
  cerebrasKey: "fz-chat-cerebras-key",
  localModel: "fz-chat-local-model",
  consented: "fz-chat-local-consent",
};

const LOCAL_MODELS = {
  "Qwen2.5-0.5B-Instruct-q4f16_1-MLC": { label: "Qwen2.5 0.5B (~400 MB, fast)", size: "~400 MB" },
  "Phi-3.5-mini-instruct-q4f16_1-MLC": { label: "Phi-3.5 mini (~2.2 GB, better)", size: "~2.2 GB" },
};
const DEFAULT_LOCAL_MODEL = "Qwen2.5-0.5B-Instruct-q4f16_1-MLC";

const DISCLAIMER = "_AI-generated - verify against the API docs._";

let _systemPrompt = null;        // cached after first fetch
let _markedPromise = null;       // cached marked module load
let _purifyPromise = null;       // cached DOMPurify module load
let _webllmEngine = null;        // cached MLCEngine instance
let _webllmLoading = null;       // promise while engine loads

function ls(k, dflt = "") {
  try { return localStorage.getItem(k) ?? dflt; } catch (_) { return dflt; }
}
function lsSet(k, v) {
  try { localStorage.setItem(k, v); } catch (_) {}
}

async function loadSystemPrompt() {
  if (_systemPrompt) return _systemPrompt;
  const res = await fetch(SYSTEM_PROMPT_URL);
  if (!res.ok) throw new Error(`failed to fetch system-prompt.txt (${res.status})`);
  _systemPrompt = await res.text();
  return _systemPrompt;
}

async function loadMarked() {
  if (!_markedPromise) {
    _markedPromise = import("https://esm.run/marked").then(m => m.marked || m.default || m);
  }
  return _markedPromise;
}

async function loadPurify() {
  if (!_purifyPromise) {
    _purifyPromise = import("https://esm.run/dompurify").then(m => m.default || m);
  }
  return _purifyPromise;
}

// Render markdown safely. If DOMPurify failed to load, fall back to plain text
// (assigning textContent) so a bad provider response can never inject script tags.
function renderMarkdown(bubble, text, marked, purify) {
  if (marked && purify) {
    bubble.innerHTML = purify.sanitize(marked.parse(text || ""));
  } else if (marked) {
    // Defensive: refuse to render unsanitized HTML; fall back to plain text.
    bubble.textContent = text || "";
  } else {
    bubble.textContent = text || "";
  }
}

function el(tag, attrs = {}, ...children) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") node.className = v;
    else if (k === "html") node.innerHTML = v;
    else if (k.startsWith("on")) node.addEventListener(k.slice(2), v);
    else node.setAttribute(k, v);
  }
  for (const c of children) {
    if (c == null) continue;
    node.appendChild(typeof c === "string" ? document.createTextNode(c) : c);
  }
  return node;
}

// ----- Providers -----

async function callGroq(messages, key, onChunk, signal) {
  const res = await fetch("https://api.groq.com/openai/v1/chat/completions", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      "Authorization": `Bearer ${key}`,
    },
    body: JSON.stringify({
      model: "llama-3.1-8b-instant",
      messages,
      stream: false,
      temperature: 0.2,
      max_tokens: 700,
    }),
    signal,
  });
  if (!res.ok) {
    const t = await res.text();
    throw new Error(`Groq API error ${res.status}: ${t.slice(0, 240)}`);
  }
  const j = await res.json();
  const text = j.choices?.[0]?.message?.content || "";
  if (!text) throw new Error("Provider returned empty response");
  onChunk(text);
  return text;
}

async function callOpenRouter(messages, key, onChunk, signal) {
  // OpenRouter exposes many free-tier models (DeepSeek V3, Llama 3.3, etc.).
  // Uses OpenAI-compatible /chat/completions; model id includes ":free" suffix.
  const res = await fetch("https://openrouter.ai/api/v1/chat/completions", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      "Authorization": `Bearer ${key}`,
      "HTTP-Referer": location.origin,
      "X-Title": "fz docs assistant",
    },
    body: JSON.stringify({
      model: "deepseek/deepseek-chat-v3.1:free",
      messages,
      stream: false,
      temperature: 0.2,
      max_tokens: 700,
    }),
    signal,
  });
  if (!res.ok) {
    const t = await res.text();
    throw new Error(`OpenRouter API error ${res.status}: ${t.slice(0, 240)}`);
  }
  const j = await res.json();
  const text = j.choices?.[0]?.message?.content || "";
  if (!text) throw new Error("Provider returned empty response");
  onChunk(text);
  return text;
}

async function callCerebras(messages, key, onChunk, signal) {
  // Cerebras inference - extremely fast, generous free tier.
  // OpenAI-compatible /chat/completions.
  const res = await fetch("https://api.cerebras.ai/v1/chat/completions", {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      "Authorization": `Bearer ${key}`,
    },
    body: JSON.stringify({
      model: "llama3.1-8b",
      messages,
      stream: false,
      temperature: 0.2,
      max_tokens: 700,
    }),
    signal,
  });
  if (!res.ok) {
    const t = await res.text();
    throw new Error(`Cerebras API error ${res.status}: ${t.slice(0, 240)}`);
  }
  const j = await res.json();
  const text = j.choices?.[0]?.message?.content || "";
  if (!text) throw new Error("Provider returned empty response");
  onChunk(text);
  return text;
}

async function callGemini(messages, key, onChunk, signal) {
  // Convert OpenAI-style messages to Gemini contents.
  // System message becomes systemInstruction; the rest become contents.
  let systemInstruction = null;
  const contents = [];
  for (const m of messages) {
    if (m.role === "system") {
      systemInstruction = { parts: [{ text: m.content }] };
    } else {
      contents.push({
        role: m.role === "assistant" ? "model" : "user",
        parts: [{ text: m.content }],
      });
    }
  }
  // Pass key via x-goog-api-key header so it never appears in URL/Referer logs.
  const url = "https://generativelanguage.googleapis.com/v1beta/models/gemini-flash-latest:generateContent";
  const body = {
    contents,
    generationConfig: { temperature: 0.2, maxOutputTokens: 700 },
  };
  if (systemInstruction) body.systemInstruction = systemInstruction;
  const res = await fetch(url, {
    method: "POST",
    headers: {
      "Content-Type": "application/json",
      "x-goog-api-key": key,
    },
    body: JSON.stringify(body),
    signal,
  });
  if (!res.ok) {
    const t = await res.text();
    throw new Error(`Gemini API error ${res.status}: ${t.slice(0, 240)}`);
  }
  const j = await res.json();
  const text = j.candidates?.[0]?.content?.parts?.map(p => p.text).join("") || "";
  if (!text) throw new Error("Provider returned empty response");
  onChunk(text);
  return text;
}

async function loadWebLLM(modelId, statusFn) {
  if (_webllmEngine && _webllmEngine._fzModelId === modelId) return _webllmEngine;
  if (_webllmLoading) return _webllmLoading;
  _webllmLoading = (async () => {
    statusFn("Loading WebLLM runtime...");
    const webllm = await import("https://esm.run/@mlc-ai/web-llm");
    if (!navigator.gpu) {
      throw new Error("This browser has no WebGPU. Try Chrome or Edge, or switch to BYOK mode.");
    }
    // Bump the context window so the ~9k-token modules.json system prompt
    // fits. Default Qwen 2.5 / Phi 3.5 MLC bundles ship a conservative 4096
    // cap that overflows on first message. Also enable sliding window
    // attention so long chats can spill past the cap gracefully.
    const engine = await webllm.CreateMLCEngine(modelId, {
      initProgressCallback: (p) => {
        const pct = Math.round((p.progress || 0) * 100);
        statusFn(`Loading model: ${p.text || ""} ${pct}%`);
      },
    }, {
      context_window_size: 16384,
      sliding_window_size: -1,
      attention_sink_size: -1,
    });
    engine._fzModelId = modelId;
    _webllmEngine = engine;
    _webllmLoading = null;
    return engine;
  })();
  return _webllmLoading;
}

async function callLocal(messages, modelId, onChunk, statusFn, signal) {
  const engine = await loadWebLLM(modelId, statusFn);
  statusFn("Generating...");
  const completion = await engine.chat.completions.create({
    messages,
    stream: true,
    temperature: 0.2,
    max_tokens: 700,
  });
  let acc = "";
  for await (const chunk of completion) {
    if (signal?.aborted) {
      // WebLLM's iterator doesn't natively accept AbortSignal; bail out
      // cooperatively so the finally block in handleSubmit re-enables Send.
      throw new DOMException("Aborted", "AbortError");
    }
    const delta = chunk.choices?.[0]?.delta?.content || "";
    if (delta) {
      acc += delta;
      onChunk(acc);
    }
  }
  if (!acc) throw new Error("Provider returned empty response");
  return acc;
}

// ----- UI -----

function buildDialog() {
  // Header
  const settingsBtn = el("button", {
    class: "fz-chat-iconbtn",
    title: "Settings",
    "aria-label": "Settings",
    type: "button",
  }, "⚙");
  const closeBtn = el("button", {
    class: "fz-chat-iconbtn",
    title: "Close",
    "aria-label": "Close",
    type: "button",
  }, "✕");
  const header = el("div", { class: "fz-chat-header" },
    el("h2", {}, "Ask the docs"),
    settingsBtn,
    closeBtn,
  );

  // Settings panel (initially hidden). Local first because it's private and
  // free with no signup once WebGPU is available.
  const providerSel = el("select", { id: "fz-provider-select" },
    el("option", { value: "local" }, "Local (WebLLM, needs WebGPU)"),
    el("option", { value: "openrouter" }, "OpenRouter (free models)"),
    el("option", { value: "cerebras" }, "Cerebras (free, very fast)"),
    el("option", { value: "gemini" }, "Google Gemini (free tier)"),
    el("option", { value: "groq" }, "Groq (free)"),
  );
  const localModelSel = el("select", { id: "fz-local-model-select" });
  for (const [id, m] of Object.entries(LOCAL_MODELS)) {
    localModelSel.appendChild(el("option", { value: id }, m.label));
  }
  const groqKeyInput = el("input", {
    type: "password",
    id: "fz-groq-key",
    placeholder: "gsk_...",
    autocomplete: "off",
  });
  const openrouterKeyInput = el("input", {
    type: "password",
    id: "fz-openrouter-key",
    placeholder: "sk-or-...",
    autocomplete: "off",
  });
  const cerebrasKeyInput = el("input", {
    type: "password",
    id: "fz-cerebras-key",
    placeholder: "csk-...",
    autocomplete: "off",
  });
  const geminiKeyInput = el("input", {
    type: "password",
    id: "fz-gemini-key",
    placeholder: "AIza...",
    autocomplete: "off",
  });
  const localBlock = el("div", { class: "fz-provider-block fz-block-local" },
    el("label", {}, "Model"),
    localModelSel,
    el("p", { class: "fz-hint", html: 'First use downloads the model to your browser cache. Runs locally - no data leaves your device. Requires WebGPU (Chrome 113+ / Safari 17.4+ / Edge desktop). <a href="/webgpu-test.html" target="_blank" rel="noopener">Test your browser &rarr;</a>' }),
  );
  const groqBlock = el("div", { class: "fz-provider-block fz-block-groq" },
    el("label", {}, "Groq API key"),
    groqKeyInput,
    el("p", { class: "fz-hint", html: 'Stored in your browser only. <a href="https://console.groq.com/keys" target="_blank" rel="noopener">Get a Groq key &rarr;</a>'}),
  );
  const openrouterBlock = el("div", { class: "fz-provider-block fz-block-openrouter" },
    el("label", {}, "OpenRouter API key"),
    openrouterKeyInput,
    el("p", { class: "fz-hint", html: 'Stored in your browser only. Routes to free models like DeepSeek V3.1. <a href="https://openrouter.ai/keys" target="_blank" rel="noopener">Get an OpenRouter key &rarr;</a>'}),
  );
  const cerebrasBlock = el("div", { class: "fz-provider-block fz-block-cerebras" },
    el("label", {}, "Cerebras API key"),
    cerebrasKeyInput,
    el("p", { class: "fz-hint", html: 'Stored in your browser only. Fastest hosted inference (~2000 tok/s). <a href="https://cloud.cerebras.ai/" target="_blank" rel="noopener">Get a Cerebras key &rarr;</a>' }),
  );
  const geminiBlock = el("div", { class: "fz-provider-block fz-block-gemini" },
    el("label", {}, "Gemini API key"),
    geminiKeyInput,
    el("p", { class: "fz-hint", html: 'Stored in your browser only. <a href="https://aistudio.google.com/apikey" target="_blank" rel="noopener">Get a Gemini key &rarr;</a>' }),
  );
  const privacyNotice = el("p", { class: "fz-chat-notice fz-hidden" }, "");
  const settingsPanel = el("div", { class: "fz-chat-settings fz-hidden" },
    el("label", {}, "Provider"),
    providerSel,
    privacyNotice,
    groqBlock,
    openrouterBlock,
    cerebrasBlock,
    geminiBlock,
    localBlock,
  );

  // Consent / model-load button (shown for local mode before first load).
  const consentBlurb = el("div", { class: "fz-consent-blurb" },
    "First-time setup: clicking the button below downloads the selected model to your browser cache. Nothing leaves your device.");
  const loadBtn = el("button", { class: "fz-chat-loadbtn", type: "button" }, "Load model (~400 MB, runs locally)");
  const loadWrap = el("div", { class: "fz-load-wrap" }, consentBlurb, loadBtn);

  // Messages area
  const messages = el("div", { class: "fz-chat-messages", role: "log", "aria-live": "polite" });

  // Status bar
  const status = el("div", { class: "fz-status" }, "Ready.");

  // Input form
  const textarea = el("textarea", {
    placeholder: "Ask about fz...",
    rows: "1",
    "aria-label": "Your question",
  });
  const sendBtn = el("button", { type: "submit" }, "Send");
  const form = el("form", { class: "fz-chat-form" }, textarea, sendBtn);

  const root = el("div", { class: "fz-chat-root" },
    header, settingsPanel, loadWrap, messages, status, form);

  return {
    root,
    settingsBtn, closeBtn,
    settingsPanel,
    providerSel, localModelSel,
    groqKeyInput, openrouterKeyInput, cerebrasKeyInput, geminiKeyInput,
    localBlock, groqBlock, openrouterBlock, cerebrasBlock, geminiBlock,
    privacyNotice,
    consentBlurb, loadBtn, loadWrap,
    messages, status, form, textarea, sendBtn,
  };
}

function makeMessageNode(role, text, marked, purify) {
  const bubble = el("div", { class: "fz-bubble" });
  if (role === "assistant") {
    renderMarkdown(bubble, text, marked, purify);
  } else {
    bubble.textContent = text;
  }
  return el("div", { class: `fz-msg ${role}` },
    el("div", { class: "fz-role" }, role === "user" ? "You" : role === "assistant" ? "Assistant" : "System"),
    bubble,
  );
}

class ChatWidget {
  constructor() {
    const dlg = el("dialog", { id: "fz-chat-dialog" });
    const ui = buildDialog();
    dlg.appendChild(ui.root);
    document.body.appendChild(dlg);

    this.dlg = dlg;
    this.ui = ui;
    this.history = [];           // { role, content }[]
    this.marked = null;
    this.purify = null;
    this.systemPrompt = null;
    this.lastAssistantNode = null;
    this.activeAbort = null;     // AbortController for the in-flight request

    // Initial provider/model from localStorage. Default to Local (WebLLM) so
    // visitors with WebGPU get a free, private experience with no signup.
    // Users without WebGPU can switch to a hosted provider in Settings; the
    // load button also auto-redirects to OpenRouter on no-GPU.
    ui.providerSel.value = ls(LS_KEYS.provider, "local");
    ui.localModelSel.value = ls(LS_KEYS.localModel, DEFAULT_LOCAL_MODEL);
    ui.groqKeyInput.value = ls(LS_KEYS.groqKey, "");
    ui.openrouterKeyInput.value = ls(LS_KEYS.openrouterKey, "");
    ui.cerebrasKeyInput.value = ls(LS_KEYS.cerebrasKey, "");
    ui.geminiKeyInput.value = ls(LS_KEYS.geminiKey, "");

    this.refreshProviderBlocks();
    this.refreshLoadButton();
    this.refreshPrivacyNotice();

    // Wire events.
    ui.settingsBtn.addEventListener("click", () => {
      ui.settingsPanel.classList.toggle("fz-hidden");
    });
    ui.closeBtn.addEventListener("click", () => dlg.close());
    dlg.addEventListener("close", () => {
      // Cancel any in-flight provider call so the Send button isn't left
      // disabled forever if the user closes mid-stream.
      if (this.activeAbort) {
        try { this.activeAbort.abort(); } catch (_) {}
      }
    });

    ui.providerSel.addEventListener("change", () => {
      lsSet(LS_KEYS.provider, ui.providerSel.value);
      this.refreshProviderBlocks();
      this.refreshLoadButton();
      this.refreshPrivacyNotice();
    });
    ui.localModelSel.addEventListener("change", () => {
      lsSet(LS_KEYS.localModel, ui.localModelSel.value);
      this.refreshLoadButton();
    });
    ui.groqKeyInput.addEventListener("change", () => lsSet(LS_KEYS.groqKey, ui.groqKeyInput.value));
    ui.openrouterKeyInput.addEventListener("change", () => lsSet(LS_KEYS.openrouterKey, ui.openrouterKeyInput.value));
    ui.cerebrasKeyInput.addEventListener("change", () => lsSet(LS_KEYS.cerebrasKey, ui.cerebrasKeyInput.value));
    ui.geminiKeyInput.addEventListener("change", () => lsSet(LS_KEYS.geminiKey, ui.geminiKeyInput.value));

    ui.loadBtn.addEventListener("click", () => this.handleLoadClick());

    ui.form.addEventListener("submit", (ev) => {
      ev.preventDefault();
      this.handleSubmit();
    });
    ui.textarea.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter" && !ev.shiftKey) {
        ev.preventDefault();
        this.handleSubmit();
      }
    });
  }

  open() { if (!this.dlg.open) this.dlg.showModal(); }

  refreshProviderBlocks() {
    const p = this.ui.providerSel.value;
    this.ui.localBlock.style.display = p === "local" ? "" : "none";
    this.ui.groqBlock.style.display = p === "groq" ? "" : "none";
    this.ui.openrouterBlock.style.display = p === "openrouter" ? "" : "none";
    this.ui.cerebrasBlock.style.display = p === "cerebras" ? "" : "none";
    this.ui.geminiBlock.style.display = p === "gemini" ? "" : "none";
  }

  refreshLoadButton() {
    const p = this.ui.providerSel.value;
    if (p !== "local") {
      this.ui.loadWrap.style.display = "none";
      return;
    }
    const modelId = this.ui.localModelSel.value || DEFAULT_LOCAL_MODEL;
    const m = LOCAL_MODELS[modelId] || LOCAL_MODELS[DEFAULT_LOCAL_MODEL];
    if (_webllmEngine && _webllmEngine._fzModelId === modelId) {
      this.ui.loadWrap.style.display = "none";
      return;
    }
    this.ui.loadWrap.style.display = "";
    this.ui.loadBtn.textContent = `Load model (${m.size}, runs locally)`;
    this.ui.loadBtn.disabled = false;
  }

  refreshPrivacyNotice() {
    const p = this.ui.providerSel.value;
    const notice = this.ui.privacyNotice;
    const PROVIDER_NAMES = {
      groq: "Groq",
      openrouter: "OpenRouter",
      cerebras: "Cerebras",
      gemini: "Google Gemini",
    };
    if (PROVIDER_NAMES[p]) {
      notice.textContent = `Note: your question is sent to ${PROVIDER_NAMES[p]}. Your API key is stored in this browser only.`;
      notice.classList.remove("fz-hidden");
    } else {
      notice.classList.add("fz-hidden");
      notice.textContent = "";
    }
  }

  setStatus(text, isError = false) {
    this.ui.status.textContent = text;
    this.ui.status.classList.toggle("fz-error", !!isError);
  }

  appendMessage(role, text) {
    const node = makeMessageNode(role, text, this.marked, this.purify);
    this.ui.messages.appendChild(node);
    this.ui.messages.scrollTop = this.ui.messages.scrollHeight;
    return node;
  }

  updateAssistantStreaming(text) {
    if (!this.lastAssistantNode) {
      this.lastAssistantNode = this.appendMessage("assistant", text);
      return;
    }
    const bubble = this.lastAssistantNode.querySelector(".fz-bubble");
    renderMarkdown(bubble, text, this.marked, this.purify);
    this.ui.messages.scrollTop = this.ui.messages.scrollHeight;
  }

  async ensurePrereqs() {
    if (!this.systemPrompt) {
      this.setStatus("Loading docs catalog...");
      this.systemPrompt = await loadSystemPrompt();
    }
    if (!this.marked) {
      try { this.marked = await loadMarked(); } catch (_) { this.marked = null; }
    }
    if (!this.purify) {
      try { this.purify = await loadPurify(); } catch (_) { this.purify = null; }
    }
  }

  async handleLoadClick() {
    const modelId = this.ui.localModelSel.value || DEFAULT_LOCAL_MODEL;
    // Pre-flight WebGPU check so failures surface immediately and direct
    // the user to BYOK instead of failing deep in WebLLM's loader.
    if (!navigator.gpu) {
      this.setStatus("No WebGPU. Switching to OpenRouter. Run /webgpu-test.html for diagnostics.", true);
      this.ui.providerSel.value = "openrouter";
      this.ui.providerSel.dispatchEvent(new Event("change", { bubbles: true }));
      this.ui.settingsPanel.classList.remove("fz-hidden");
      return;
    }
    try {
      const adapter = await navigator.gpu.requestAdapter();
      if (!adapter) {
        this.setStatus("WebGPU API present but no adapter. See /webgpu-test.html for fixes.", true);
        this.ui.providerSel.value = "openrouter";
        this.ui.providerSel.dispatchEvent(new Event("change", { bubbles: true }));
        this.ui.settingsPanel.classList.remove("fz-hidden");
        return;
      }
    } catch (_) { /* fall through and let the loader surface the real error */ }

    this.ui.loadBtn.disabled = true;
    lsSet(LS_KEYS.consented, "1");
    try {
      await this.ensurePrereqs();
      await loadWebLLM(modelId, (s) => this.setStatus(s));
      this.setStatus("Model ready.");
      this.refreshLoadButton();
    } catch (err) {
      this.setStatus(`Load failed: ${err.message || err}`, true);
      this.ui.loadBtn.disabled = false;
    }
  }

  async handleSubmit() {
    const text = this.ui.textarea.value.trim();
    if (!text) return;
    const provider = this.ui.providerSel.value;

    // Validate provider readiness.
    const KEY_INPUTS = {
      groq:       this.ui.groqKeyInput,
      openrouter: this.ui.openrouterKeyInput,
      cerebras:   this.ui.cerebrasKeyInput,
      gemini:     this.ui.geminiKeyInput,
    };
    const KEY_NAMES = { groq: "Groq", openrouter: "OpenRouter", cerebras: "Cerebras", gemini: "Gemini" };
    if (KEY_INPUTS[provider] && !KEY_INPUTS[provider].value.trim()) {
      this.setStatus(`Add a ${KEY_NAMES[provider]} API key in Settings (gear icon).`, true);
      this.ui.settingsPanel.classList.remove("fz-hidden");
      return;
    }
    if (provider === "local") {
      const modelId = this.ui.localModelSel.value || DEFAULT_LOCAL_MODEL;
      if (!_webllmEngine || _webllmEngine._fzModelId !== modelId) {
        this.setStatus("Click 'Load model' first to download the local model.", true);
        return;
      }
    }

    this.ui.textarea.value = "";
    this.ui.sendBtn.disabled = true;
    this.appendMessage("user", text);
    this.history.push({ role: "user", content: text });
    this.lastAssistantNode = null;

    // Per-submit AbortController. The dialog 'close' listener aborts this so
    // the Send button is never left disabled if the user closes mid-stream.
    const abortCtrl = new AbortController();
    this.activeAbort = abortCtrl;

    try {
      await this.ensurePrereqs();
      const messages = [
        { role: "system", content: this.systemPrompt },
        ...this.history,
      ];
      this.setStatus("Thinking...");
      const onChunk = (txt) => this.updateAssistantStreaming(txt);
      let reply = "";
      if (provider === "local") {
        reply = await callLocal(messages, this.ui.localModelSel.value || DEFAULT_LOCAL_MODEL,
                                onChunk, (s) => this.setStatus(s), abortCtrl.signal);
      } else if (provider === "groq") {
        reply = await callGroq(messages, this.ui.groqKeyInput.value.trim(), onChunk, abortCtrl.signal);
      } else if (provider === "openrouter") {
        reply = await callOpenRouter(messages, this.ui.openrouterKeyInput.value.trim(), onChunk, abortCtrl.signal);
      } else if (provider === "cerebras") {
        reply = await callCerebras(messages, this.ui.cerebrasKeyInput.value.trim(), onChunk, abortCtrl.signal);
      } else if (provider === "gemini") {
        reply = await callGemini(messages, this.ui.geminiKeyInput.value.trim(), onChunk, abortCtrl.signal);
      } else {
        throw new Error(`unknown provider: ${provider}`);
      }
      // Append disclaimer if the model didn't echo it back already.
      if (reply && !/AI[- ]generated/i.test(reply)) {
        reply = reply + "\n\n" + DISCLAIMER;
        this.updateAssistantStreaming(reply);
      }
      this.history.push({ role: "assistant", content: reply });
      this.setStatus("Ready.");
    } catch (err) {
      const isAbort = err && (err.name === "AbortError" || /aborted/i.test(err.message || ""));
      if (isAbort) {
        // Silent on user-initiated abort; just clear the status.
        this.setStatus("Cancelled.");
      } else {
        const msg = err.message || String(err);
        this.appendMessage("system", `Error: ${msg}`);
        this.setStatus(`Error: ${msg}`, true);
      }
    } finally {
      // Always re-enable the Send button, no matter how we exited.
      this.ui.sendBtn.disabled = false;
      if (this.activeAbort === abortCtrl) this.activeAbort = null;
    }
  }
}

function mountFAB() {
  if (document.getElementById("fz-chat-fab")) return; // already mounted
  // Inject CSS.
  const cssHref = new URL("chat.css", CHAT_ROOT).href;
  if (!document.querySelector(`link[href="${cssHref}"]`)) {
    document.head.appendChild(el("link", { rel: "stylesheet", href: cssHref }));
  }
  let widget = null;
  const fab = el("button", {
    id: "fz-chat-fab",
    type: "button",
    "aria-label": "Open AI chat assistant",
    title: "Ask the docs (AI)",
  }, "✨ Ask AI");
  fab.addEventListener("click", () => {
    if (!widget) widget = new ChatWidget();
    widget.open();
  });
  document.body.appendChild(fab);
}

if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", mountFAB);
} else {
  mountFAB();
}

export { mountFAB };
