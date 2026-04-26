/* fz docs chat assistant - vanilla JS, no framework. ES module.
 *
 * Two paths, no API keys ever:
 *
 *   1. Local in-browser chat via WebLLM + WebGPU. Default Qwen 2.5 0.5B
 *      (~400 MB), upgradeable to Phi-3.5 mini. All inference on the
 *      user's machine; nothing leaves their browser.
 *
 *   2. "Open in your own AI account" — opens ChatGPT / Claude / Gemini /
 *      NotebookLM in a new tab, prefilled with a bootstrap prompt that
 *      points the assistant at our system-prompt.txt URL. The user is
 *      already logged in there, so no key is ever pasted into our chat.
 */

// Resolve asset paths relative to this script so the widget works whether
// the host page lives at /, /modules/, or any other subdirectory.
const SCRIPT_URL = new URL(import.meta.url);
const CHAT_ROOT = new URL("./", SCRIPT_URL);
const SYSTEM_PROMPT_URL = new URL("system-prompt.txt", CHAT_ROOT).href;

const LS_KEYS = {
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

// ----- Local provider (WebLLM) -----

async function loadWebLLM(modelId, statusFn) {
  if (_webllmEngine && _webllmEngine._fzModelId === modelId) return _webllmEngine;
  if (_webllmLoading) return _webllmLoading;
  _webllmLoading = (async () => {
    statusFn("Loading WebLLM runtime...");
    const webllm = await import("https://esm.run/@mlc-ai/web-llm");
    if (!navigator.gpu) {
      throw new Error("This browser has no WebGPU. Try Chrome or Edge desktop, or use 'Open in your own AI account' instead.");
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
  const closeBtn = el("button", {
    class: "fz-chat-iconbtn",
    title: "Close",
    "aria-label": "Close",
    type: "button",
  }, "✕");
  const header = el("div", { class: "fz-chat-header" },
    el("h2", {}, "Ask about FZ"),
    closeBtn,
  );

  // Settings panel (initially hidden). Two paths only:
  //   1. Local in-browser (WebLLM) — model selector below.
  //   2. Open in your own AI account — buttons below.
  // No API key inputs. The user is either using local or jumping out to a
  // service where they're already logged in.
  const localModelSel = el("select", { id: "fz-local-model-select" });
  for (const [id, m] of Object.entries(LOCAL_MODELS)) {
    localModelSel.appendChild(el("option", { value: id }, m.label));
  }
  const localBlock = el("div", { class: "fz-provider-block fz-block-local" },
    el("label", {}, "In-browser model (WebLLM)"),
    localModelSel,
    el("p", { class: "fz-hint", html: 'First use downloads the model to your browser cache. Runs locally - no data leaves your device. Requires WebGPU (Chrome 113+ / Safari 17.4+ / Edge desktop). <a href="/webgpu-test.html" target="_blank" rel="noopener">Test your browser &rarr;</a>' }),
  );

  // External AI links: open ChatGPT / Claude / Gemini in the user's logged-in
  // account, pre-filled with a bootstrap that asks the assistant to fetch our
  // system-prompt.txt for context. Clipboard fallback covers cases where the
  // LLM can't fetch URLs (or the docs are on localhost). Gemini routes to
  // Google AI Studio because gemini.google.com itself doesn't honour ?prompt=
  // natively (only via browser extension); aistudio.google.com does.
  const extChatGPT  = el("button", { class: "fz-ext-btn", type: "button", "data-target": "chatgpt" }, "ChatGPT");
  const extClaude   = el("button", { class: "fz-ext-btn", type: "button", "data-target": "claude"  }, "Claude");
  const extGemini   = el("button", { class: "fz-ext-btn", type: "button", "data-target": "gemini",
                                     title: "Opens Google AI Studio (gemini.google.com doesn't accept URL params)" }, "Gemini");
  const extCopyBtn  = el("button", { class: "fz-ext-copybtn", type: "button" }, "copy prompt to clipboard");
  const extStatus   = el("span", { class: "fz-ext-status" }, "");
  const externalBlock = el("div", { class: "fz-external-block fz-external-primary" },
    el("label", {}, "Option 1: Open in your AI account"),
    el("div", { class: "fz-external-row" }, extClaude, extChatGPT, extGemini),
    el("p", { class: "fz-hint" },
      "Opens a new tab in your logged-in account with a prompt that loads the FZ catalog as context. ",
      "Gemini routes to Google AI Studio. If the assistant can't fetch URLs (or you're on localhost), ",
      extCopyBtn, "."),
    el("div", { class: "fz-external-fallback" }, extStatus),
  );

  // The local-WebLLM section is Option 2 — collapsed by default
  // inside a <details> element. User clicks the disclosure to expand.
  const localDetails = el("details", { class: "fz-local-details" });
  const localSummary = el("summary", { class: "fz-local-summary" },
    "Option 2: Chat in this browser (local model, needs WebGPU)");
  localDetails.appendChild(localSummary);
  localDetails.appendChild(localBlock);

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
    header,
    externalBlock,         // primary path: open in user's own AI account
    localDetails,          // backup: WebLLM, collapsed by default
    loadWrap,              // load-model button (visible when local needs loading)
    messages,
    status,
    form,                  // textarea + Send (Send routes to local WebLLM)
  );

  return {
    root,
    closeBtn,
    localDetails,
    localModelSel, localBlock,
    extChatGPT, extClaude, extGemini, extCopyBtn, extStatus,
    consentBlurb, loadBtn, loadWrap,
    messages, status, form, textarea, sendBtn,
  };
}

// External-AI URL builders. Each returns the URL of the assistant's chat
// page pre-filled with a bootstrap prompt that points at our catalog.
//   - ChatGPT: ?q= (auto-submits as first message)
//   - Claude:  ?q= on /new (same)
//   - Gemini:  routes to Google AI Studio's ?prompt= on /prompts/new_chat
//              (gemini.google.com itself doesn't honour URL params natively).
function buildExternalUrl(target, promptUrl, question) {
  const intro =
    "I'm asking about FZ, a modular C++17 toolkit for error-bounded scientific data compression. " +
    `Please fetch this catalog and use it as authoritative context for every answer: ${promptUrl}` +
    " Only recommend ALGO_* names and module names that appear in the catalog; if the catalog doesn't cover something, say so.";
  const ask = question.trim() ? `\n\nMy question: ${question.trim()}` : "";
  const full = intro + ask;
  switch (target) {
    case "chatgpt":
      return "https://chatgpt.com/?q=" + encodeURIComponent(full);
    case "claude":
      return "https://claude.ai/new?q=" + encodeURIComponent(full);
    case "gemini":
      return "https://aistudio.google.com/prompts/new_chat?prompt=" + encodeURIComponent(full);
    default:
      return null;
  }
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

    ui.localModelSel.value = ls(LS_KEYS.localModel, DEFAULT_LOCAL_MODEL);
    this.refreshLoadButton();

    ui.closeBtn.addEventListener("click", () => dlg.close());
    dlg.addEventListener("close", () => {
      // Cancel any in-flight WebLLM call so the Send button isn't stuck
      // if the user closes mid-stream.
      if (this.activeAbort) {
        try { this.activeAbort.abort(); } catch (_) {}
      }
    });

    ui.localModelSel.addEventListener("change", () => {
      lsSet(LS_KEYS.localModel, ui.localModelSel.value);
      this.refreshLoadButton();
    });

    ui.loadBtn.addEventListener("click", () => this.handleLoadClick());

    // External-AI buttons: open the user's chosen account in a new tab.
    for (const btn of [ui.extChatGPT, ui.extClaude, ui.extGemini]) {
      btn.addEventListener("click", () => this.handleExternalClick(btn.dataset.target));
    }
    ui.extCopyBtn.addEventListener("click", () => this.handleCopyCatalog());

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

  refreshLoadButton() {
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

  // Open ChatGPT / Claude in the user's logged-in account, pre-filled with
  // a bootstrap that fetches our system prompt URL. SYSTEM_PROMPT_URL is
  // resolved at runtime from import.meta.url, so it auto-uses the deployed
  // host (localhost:9001 in dev, szcompressor.org/SZ3/... in prod).
  handleExternalClick(target) {
    const promptUrl = new URL(SYSTEM_PROMPT_URL).href;
    const question = (this.ui.textarea.value || "").trim();
    const url = buildExternalUrl(target, promptUrl, question);
    if (!url) return;
    window.open(url, "_blank", "noopener");
  }

  async handleCopyCatalog() {
    try {
      await this.ensurePrereqs();
      await navigator.clipboard.writeText(this.systemPrompt || "");
      const n = (this.systemPrompt || "").length;
      this.ui.extStatus.textContent = `Copied ${n.toLocaleString()} chars. Paste into your assistant chat.`;
    } catch (e) {
      this.ui.extStatus.textContent = "Copy failed: " + (e.message || e);
    }
    setTimeout(() => { this.ui.extStatus.textContent = ""; }, 6000);
  }

  async handleLoadClick() {
    const modelId = this.ui.localModelSel.value || DEFAULT_LOCAL_MODEL;
    // Pre-flight WebGPU check so failures surface immediately. The
    // external-AI buttons above the WebLLM section are the recommended
    // alternative when WebGPU isn't available.
    if (!navigator.gpu) {
      this.setStatus("This browser has no WebGPU. Use one of the 'Open in your AI account' buttons above instead, or run /webgpu-test.html for diagnostics.", true);
      return;
    }
    try {
      const adapter = await navigator.gpu.requestAdapter();
      if (!adapter) {
        this.setStatus("WebGPU API present but no compatible GPU adapter. Use 'Open in your AI account' above, or see /webgpu-test.html.", true);
        return;
      }
    } catch (_) { /* let WebLLM's loader surface the real error */ }

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

    // Send routes only to local WebLLM. External providers are reached via
    // the 'Open in your AI account' buttons above (they open a new tab).
    const modelId = this.ui.localModelSel.value || DEFAULT_LOCAL_MODEL;
    if (!_webllmEngine || _webllmEngine._fzModelId !== modelId) {
      this.setStatus("Click 'Load model' first (in the Local section above), or open the question in your AI account using the buttons above.", true);
      // Auto-expand the local <details> so the user sees the load button.
      if (this.ui.localDetails) this.ui.localDetails.open = true;
      return;
    }

    this.ui.textarea.value = "";
    this.ui.sendBtn.disabled = true;
    this.appendMessage("user", text);
    this.history.push({ role: "user", content: text });
    this.lastAssistantNode = null;

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
      let reply = await callLocal(messages, modelId,
                                  onChunk, (s) => this.setStatus(s), abortCtrl.signal);
      if (reply && !/AI[- ]generated/i.test(reply)) {
        reply = reply + "\n\n" + DISCLAIMER;
        this.updateAssistantStreaming(reply);
      }
      this.history.push({ role: "assistant", content: reply });
      this.setStatus("Ready.");
    } catch (err) {
      const isAbort = err && (err.name === "AbortError" || /aborted/i.test(err.message || ""));
      if (isAbort) {
        this.setStatus("Cancelled.");
      } else {
        const msg = err.message || String(err);
        this.appendMessage("system", `Error: ${msg}`);
        this.setStatus(`Error: ${msg}`, true);
      }
    } finally {
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
