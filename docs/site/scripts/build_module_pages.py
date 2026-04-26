#!/usr/bin/env python3
"""Generate one stub HTML page per fz module/algo entry for Pagefind to index.

Reads docs/modules.json + docs/module_page_template.html, writes one
<id>.html per entry plus an index.html into the output directory.

Fails loud if any referenced wiring_file or header path doesn't exist.

Use --validate-json to run an internal-consistency pass over modules.json
without writing any files (used in CI / pre-commit).
"""
import argparse
import html
import json
import sys
from pathlib import Path


# --- Sentinels for empty fields (fix for empty-array rendering) ---
EMPTY_LIST_HTML = "<li><em>(none listed)</em></li>"
EMPTY_KEYWORDS_HTML = '<span class="fz-meta"><em>(none)</em></span>'
EMPTY_REFS_HTML = '<li class="fz-meta"><em>References to be added during paper-mining pass.</em></li>'


def render_keywords(kws):
    if not kws:
        return EMPTY_KEYWORDS_HTML
    return " ".join(f"<span>{html.escape(k)}</span>" for k in kws)


def render_data_examples(items):
    if not items:
        return EMPTY_LIST_HTML
    return "".join(f"<li>{html.escape(x)}</li>" for x in items)


def render_pipeline(pipe):
    parts = []
    for stage in ("decomposition", "quantizer", "encoder", "lossless"):
        v = pipe.get(stage, "(none)")
        if v and v != "(none)":
            parts.append(f"{stage}: {html.escape(str(v))}")
        else:
            parts.append(f"{stage}: <em>(none)</em>")
    return "<br/>".join(parts)


def render_pipeline_notes(notes):
    """Optional prose annotation on the pipeline (P0 #2)."""
    if not notes:
        return ""
    return f'<p class="fz-meta"><em>{html.escape(notes)}</em></p>'


def render_references(refs):
    if not refs:
        return EMPTY_REFS_HTML
    out = []
    for r in refs:
        title = html.escape(r.get("title", "Untitled"))
        venue = html.escape(r.get("venue", ""))
        year = r.get("year", "")
        evidence = html.escape(r.get("evidence", "")) if r.get("evidence") else ""
        cite_bits = [title]
        if venue or year:
            cite_bits.append(f'<span class="year">{venue} {year}</span>'.strip())
        line = " &middot; ".join(cite_bits)
        if evidence:
            line += f'<span class="evidence">{evidence}</span>'
        out.append(f"<li>{line}</li>")
    return "".join(out)


def extract_doxygen_navbar(doxygen_index_html):
    """Pull out the Doxygen menu/topbar block from a generated HTML page.

    Returns a dict with three string members:
      head_extra: <script>/<link> tags to add inside <head> (jQuery, menu.js, ...)
      body_init:  inline initMenu(...) script + <div id="main-nav"> placeholder
                  to drop right after the title-area block.

    The path is one level up (catalog pages live under docs/html/modules/),
    so URLs that were relative to docs/html/ get prefixed with "../".
    """
    import re
    text = Path(doxygen_index_html).read_text()

    # Capture the <head> JS includes Doxygen relies on for the navbar/menu.
    # We grab from the first jquery script through dynsections (one block).
    head_match = re.search(
        r'(<script type="text/javascript" src="jquery\.js"></script>.*?'
        r'<script type="text/javascript" src="dynsections\.js"></script>)',
        text,
        re.DOTALL,
    )
    head_extra = head_match.group(1) if head_match else ""

    # Capture the <body> block that loads menudata.js + menu.js + calls
    # initMenu(...) and ends with the #main-nav placeholder.
    body_match = re.search(
        r'(<script type="text/javascript" src="menudata\.js"></script>.*?'
        r'<div id="main-nav"></div>)',
        text,
        re.DOTALL,
    )
    body_init = body_match.group(1) if body_match else '<div id="main-nav"></div>'

    # Rewrite the few asset paths so they resolve from docs/html/modules/
    # instead of docs/html/.  jquery.js, menudata.js, menu.js, dynsections.js
    # all live in the parent directory.
    for asset in ("jquery.js", "dynsections.js", "menudata.js", "menu.js"):
        head_extra = head_extra.replace(f'src="{asset}"', f'src="../{asset}"')
        body_init = body_init.replace(f'src="{asset}"', f'src="../{asset}"')

    # initMenu(relPath, ...) — its first arg is the prefix prepended to
    # every URL in menudata.js. From docs/html/modules/*.html the menu's
    # URLs (which are written relative to docs/html/) need a "../" prefix.
    # Match any combination of trailing args (search-engine on/off changes them).
    body_init = re.sub(
        r"initMenu\('',\s*",
        "initMenu('../',",
        body_init,
    )

    return {"head_extra": head_extra, "body_init": body_init}


def render_entry(entry, template, repo_root, navbar=None, fz_version=""):
    wiring_file = entry.get("wiring_file") or entry.get("header") or ""
    if wiring_file:
        full = repo_root / wiring_file
        if not full.exists():
            sys.exit(f"ERROR: wiring/header file does not exist for {entry['id']}: {wiring_file}")

    category = entry.get("category", "module")
    category_class = "algo" if category == "algorithm" else "module"
    title = entry.get("name") or entry["id"]

    # Wiring display: if no wiring_file, fall back to wiring_note prose.
    # The $wiring_block placeholder owns the whole <p> so we avoid
    # <code><em>...</em></code> when there's no real path.
    if wiring_file:
        wiring_block = f'<p class="fz-meta">Wired in <code>{html.escape(wiring_file)}</code></p>'
    elif entry.get("wiring_note"):
        wiring_block = f'<p class="fz-meta">{html.escape(entry["wiring_note"])}</p>'
    else:
        wiring_block = '<p class="fz-meta"><em>(no dedicated wiring file)</em></p>'

    # Pagefind ranking weights:
    #   - ALGO_* (algorithm) pages get h1=10.0 — these are the canonical
    #     user-facing entry points. They must outrank lower-level module
    #     pages on shared keywords ("interpolation", "svd", "zstd", etc.).
    #   - Niche/specialized algo variants (entry["niche_variant"] == true)
    #     drop to h1=6.0 so they don't outrank the canonical algo when the
    #     query phrase appears in both titles (e.g. "molecular dynamics"
    #     must land on ALGO_BIOMD, not ALGO_BIOMDXTC).
    #   - Module pages get h1=3.5 — they appear in results, but only win
    #     when the query is module-specific (e.g. "BitplaneEncoder").
    # Algorithms additionally get a body-weight bump on their leading
    # paragraph (canonical_intro below) so verbatim canonical search
    # terms are weighted strongly even outside the h1.
    if category == "algorithm":
        h1_weight = "6.0" if entry.get("niche_variant") else "10.0"
    else:
        h1_weight = "3.5"

    # Optional canonical-intro sentence: lead with verbatim search terms
    # so Pagefind's title+near-h1 proximity boost picks the *canonical*
    # algo over more specialized variants ("molecular dynamics" must
    # land on ALGO_BIOMD, not ALGO_BIOMDXTC; "zstd" must land on
    # ALGO_LOSSLESS, not Lossless_zstd; etc.). Body weight 3.0 means a
    # term in this paragraph counts ~3x compared to plain body text.
    canonical_terms = entry.get("canonical_terms") or []
    if canonical_terms and category == "algorithm":
        terms_str = ", ".join(html.escape(t) for t in canonical_terms)
        canonical_intro = (
            f'<p class="fz-canonical" data-pagefind-weight="3.0">'
            f'<strong>Canonical fz pipeline for:</strong> {terms_str}.</p>'
        )
    else:
        canonical_intro = ""

    subs = {
        "$title": html.escape(title),
        "$id": html.escape(entry["id"]),
        "$category": html.escape(category),
        "$category_class": category_class,
        "$best_for": html.escape(entry.get("best_for", "")),
        "$pipeline": render_pipeline(entry.get("pipeline", {})),
        "$pipeline_notes": render_pipeline_notes(entry.get("pipeline_notes", "")),
        "$wiring_block": wiring_block,
        "$wiring_file": html.escape(wiring_file),
        "$supported_dims": html.escape(", ".join(map(str, entry.get("supported_dims", [])))) or "(any)",
        "$supported_dtypes": html.escape(", ".join(entry.get("supported_dtypes", []))) or "(any)",
        "$error_modes": html.escape(", ".join(entry.get("error_modes", []))) or "(any)",
        "$data_examples": render_data_examples(entry.get("data_examples", [])),
        "$keywords": render_keywords(entry.get("keywords", [])),
        "$references": render_references(entry.get("references", [])),
        "$h1_weight": h1_weight,
        "$canonical_intro": canonical_intro,
        "$doxygen_navbar_head": (navbar or {}).get("head_extra", ""),
        "$doxygen_navbar_body": (navbar or {}).get("body_init", ""),
        "$projectnumber": html.escape(fz_version),
    }
    out = template
    for k in sorted(subs, key=len, reverse=True):
        out = out.replace(k, subs[k])
    return out


def render_index(algos, modules, navbar=None, fz_version=""):
    def section(title, entries):
        if not entries:
            return f"<h2>{html.escape(title)}</h2><p><em>(none)</em></p>"
        rows = []
        for e in entries:
            badge_class = "algo" if e.get("category") == "algorithm" else "module"
            best = html.escape(e.get("best_for", "")[:300])
            rows.append(
                f'<li><a href="{html.escape(e["id"])}.html"><strong>{html.escape(e.get("name", e["id"]))}</strong></a> '
                f'<span class="fz-badge {badge_class}">{html.escape(e.get("category", ""))}</span>'
                f'<br/><span class="fz-meta">{best}</span></li>'
            )
        return f'<h2>{html.escape(title)}</h2><ul class="fz-index">{"".join(rows)}</ul>'

    head_extra = (navbar or {}).get("head_extra", "")
    body_init = (navbar or {}).get("body_init", "")
    return f"""<!DOCTYPE html>
<html><head>
<meta charset="UTF-8"/>
<link rel="icon" type="image/png" sizes="32x32" href="../favicon.png"/>
<link rel="shortcut icon" href="../favicon.png"/>
<title>FZ Catalog: All algorithms and modules</title>
<link href="../tabs.css" rel="stylesheet" type="text/css"/>
<link href="../doxygen.css" rel="stylesheet" type="text/css"/>
<link href="../pagefind.css" rel="stylesheet" type="text/css"/>
<link href="../pagefind/pagefind-ui.css" rel="stylesheet"/>
<script src="../pagefind/pagefind-ui.js"></script>
{head_extra}
<script>
window.addEventListener('DOMContentLoaded', () => {{
    new PagefindUI({{
        element: "#pagefindsearch",
        showSubResults: true,
        translations: {{ placeholder: "" }}
    }});
    document.addEventListener('click', (e) => {{
        const search = document.getElementById('pagefindsearch');
        if (!search) return;
        const input = search.querySelector('.pagefind-ui__search-input');
        if (!input || search.contains(e.target)) return;
        if (input.value) {{
            input.value = '';
            input.dispatchEvent(new Event('input',  {{ bubbles: true }}));
            input.dispatchEvent(new Event('change', {{ bubbles: true }}));
        }}
    }});
}});
</script>
<!-- fz docs AI chat assistant: floating "Ask AI" button. -->
<script type="module" src="../assets/chat/chat.js"></script>
<style>
body {{ font-family: Roboto, sans-serif; }}
.fz-module-page {{ max-width: 980px; margin: 1.5em auto; padding: 0 1em; }}
.fz-badge {{ display: inline-block; background: #e8edf2; padding: 0.1em 0.5em; border-radius: 3px; font-size: 0.85em; }}
.fz-badge.algo {{ background: #d4edda; color: #155724; }}
.fz-badge.module {{ background: #fff3cd; color: #856404; }}
.fz-meta {{ color: #555; font-size: 0.92em; }}
.fz-index li {{ margin: 0.55em 0; list-style: none; padding-left: 0.5em; border-left: 3px solid #0366d6; }}
.fz-index {{ padding-left: 0; }}
</style>
</head>
<body>
<div id="top"><!-- doxygen requires this top div for the menubar -->
<div id="titlearea"><table cellspacing="0" cellpadding="0" style="width:100%"><tbody><tr style="height: 56px;">
<td id="projectalign" style="padding-left: 0.5em; width: 100%;">
  <div id="projectname">FZ Module Library&#160;<span id="projectnumber">{html.escape(fz_version)}</span></div>
  <div id="projectbrief">Modular toolkit for error-bounded scientific data compression</div>
</td>
<td class="search-bars-cell">
  <div class="search-bars">
    <div class="search-field search-field--fulltext">
      <span class="search-field__label">Search</span>
      <div id="pagefindsearch" class="search-field__mount"></div>
    </div>
  </div>
</td>
</tr></tbody></table></div>
{body_init}
</div><!-- /top -->
<div class="fz-module-page">
<p><a href="../index.html">&larr; Home</a></p>
<h1>Algorithm Catalog</h1>
<p>Browse the available compression algorithms and reusable pipeline modules in fz.
Use the search box above to find modules by data type, application domain, or technique.
For C++ symbol lookup (classes, functions, files) see the
<a href="../annotated.html">API reference</a>.</p>
{section("Algorithms", algos)}
{section("Pipeline modules", modules)}
</div>
</body></html>
"""


# ---------- internal-consistency validator (P1 #10) ----------

# Pipeline-stage values that should be ignored when checking module references.
# These are bracketed sentinels like "(internal)", "(any)", composite labels
# like "InterpolationDecomposition | BlockwiseDecomposition", or generic
# downstream markers.
_PIPELINE_NONREF_PREFIXES = ("(", "<")


def _pipeline_module_refs(pipe_dict):
    """Yield individual module IDs referenced by a pipeline entry.

    Splits on '|' and '+' so composite labels like 'A | B' or 'A + B' are
    decomposed into their components. Sentinels (parenthesized) are skipped.
    """
    for stage in ("decomposition", "quantizer", "encoder", "lossless"):
        v = pipe_dict.get(stage)
        if not v or not isinstance(v, str):
            continue
        v = v.strip()
        if not v or v.startswith(_PIPELINE_NONREF_PREFIXES):
            continue
        # Split composite "A | B" or "A + B" labels.
        for part in v.replace("+", "|").split("|"):
            tok = part.strip()
            # Drop trailing parenthetical detail like "X (internal)".
            paren = tok.find("(")
            if paren > 0:
                tok = tok[:paren].strip()
            if tok and not tok.startswith(_PIPELINE_NONREF_PREFIXES):
                yield tok


def validate(data):
    """Internal-consistency pass over modules.json.

    Returns (errors, warnings). Empty errors list means CI-pass.
    Warnings are advisory (e.g. module-shaped names not present in modules[]).
    """
    errs = []
    warns = []
    algos = data.get("algos", [])
    modules = data.get("modules", [])

    # Index by id; check for duplicates.
    by_id = {}
    for entry in algos + modules:
        eid = entry.get("id")
        if not eid:
            errs.append(f"entry missing 'id': {entry}")
            continue
        if eid in by_id:
            errs.append(f"duplicate id: {eid}")
        by_id[eid] = entry

    algo_ids = {a["id"] for a in algos if a.get("id")}
    module_ids = {m["id"] for m in modules if m.get("id")}

    # 1. used_by on modules must point at real algos. (HARD)
    for m in modules:
        for ref in m.get("used_by", []) or []:
            if ref not in algo_ids:
                errs.append(f"module {m.get('id')!r}: used_by references unknown algo {ref!r}")

    # 2. Pipeline stages on algos may reference module-shaped names that aren't
    #    indexed as standalone modules in modules[] (e.g. internal-only decomps
    #    like BlockwiseDecomposition). Emit as advisory warnings only -
    #    modules[] is curated to "distinctive" reusable modules, not every
    #    internal class. (SOFT)
    for a in algos:
        for ref in _pipeline_module_refs(a.get("pipeline", {}) or {}):
            if ref not in module_ids:
                suffixes = ("Decomposition", "Encoder", "Quantizer", "Predictor", "Compressor")
                if any(ref.endswith(s) for s in suffixes):
                    warns.append(
                        f"algo {a.get('id')!r}: pipeline references {ref!r} which is not in modules[]; "
                        f"if reusable, consider adding a module entry."
                    )

    # 3. Required fields per entry. (HARD)
    for entry in algos + modules:
        eid = entry.get("id", "<unknown>")
        for field in ("name", "category", "best_for"):
            if not entry.get(field):
                errs.append(f"entry {eid!r}: missing required field {field!r}")
        cat = entry.get("category")
        if cat not in {"algorithm", "decomposition", "quantizer", "encoder", "lossless", "module"}:
            errs.append(f"entry {eid!r}: unknown category {cat!r}")

    # 4. References should at least have title + venue + year if present. (HARD)
    for entry in algos + modules:
        for i, ref in enumerate(entry.get("references", []) or []):
            for f in ("title", "venue", "year"):
                if not ref.get(f):
                    errs.append(f"entry {entry.get('id')!r} reference[{i}]: missing {f!r}")

    return errs, warns


def main():
    p = argparse.ArgumentParser()
    p.add_argument("modules_json")
    p.add_argument("output_dir", nargs="?", default=None,
                   help="Output directory (omit when using --validate-json)")
    p.add_argument("--template", default=None,
                   help="Path to template HTML (default: docs/module_page_template.html)")
    p.add_argument("--repo-root", default=None,
                   help="Path to repo root (default: parent of modules_json's directory)")
    p.add_argument("--validate-json", action="store_true",
                   help="Run internal-consistency checks on modules.json and exit; do not write files.")
    args = p.parse_args()

    modules_json = Path(args.modules_json).resolve()
    # Default repo root: walk up from modules.json (docs/site/modules.json -> repo root).
    repo_root = Path(args.repo_root).resolve() if args.repo_root else modules_json.parent.parent.parent

    # Read FZ version from CMakeLists.txt so module pages render the same
    # "FZ Module Library 3.3.2" titlebar as Doxygen pages do.
    import os, re as _re
    fz_version = os.environ.get("FZ_VERSION", "")
    if not fz_version:
        try:
            cml = (repo_root / "CMakeLists.txt").read_text()
            m = _re.search(r"project\(SZ3\s+VERSION\s+([0-9.]+)\)", cml)
            if m: fz_version = m.group(1)
        except Exception:
            pass

    with open(modules_json) as f:
        data = json.load(f)

    if args.validate_json:
        errs, warns = validate(data)
        for w in warns:
            print(f"  warn: {w}", file=sys.stderr)
        if errs:
            print(f"VALIDATION FAILED ({len(errs)} errors):", file=sys.stderr)
            for e in errs:
                print(f"  - {e}", file=sys.stderr)
            sys.exit(2)
        print(f"OK: {len(data.get('algos', []))} algos + {len(data.get('modules', []))} modules; "
              f"no internal-consistency errors ({len(warns)} advisory warnings).")
        return

    if not args.output_dir:
        sys.exit("ERROR: output_dir is required unless --validate-json is set")
    output_dir = Path(args.output_dir).resolve()
    template_path = Path(args.template) if args.template else modules_json.parent / "module_page_template.html"

    if not template_path.exists():
        sys.exit(f"ERROR: template not found: {template_path}")

    template = template_path.read_text()
    output_dir.mkdir(parents=True, exist_ok=True)

    # Pull Doxygen's main-nav block out of the freshly generated index.html
    # so the catalog pages share the topbar (Main Page / Catalog / Classes /
    # Files / ...). If the Doxygen output isn't available yet, fall back to
    # an empty navbar — pages still render, just without the menu.
    navbar = None
    doxygen_index = output_dir.parent / "index.html"
    if doxygen_index.exists():
        try:
            navbar = extract_doxygen_navbar(doxygen_index)
        except Exception as ex:  # pragma: no cover
            print(f"warn: could not extract Doxygen navbar from {doxygen_index}: {ex}",
                  file=sys.stderr)
            navbar = None
    else:
        print(f"warn: Doxygen index not found at {doxygen_index}; "
              f"catalog pages will lack the topbar nav.", file=sys.stderr)

    algos = data.get("algos", [])
    modules = data.get("modules", [])

    for entry in algos + modules:
        rendered = render_entry(entry, template, repo_root, navbar=navbar, fz_version=fz_version)
        out_path = output_dir / f"{entry['id']}.html"
        out_path.write_text(rendered)

    index_path = output_dir / "index.html"
    index_path.write_text(render_index(algos, modules, navbar=navbar, fz_version=fz_version))

    total = len(algos) + len(modules)
    print(f"Wrote {total} module pages + index to {output_dir}")
    print(f"  algos: {len(algos)}, modules: {len(modules)}")


if __name__ == "__main__":
    main()
