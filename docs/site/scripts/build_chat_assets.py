#!/usr/bin/env python3
"""Build the chat assistant's system prompt.

Reads docs/modules.json and writes docs/html/assets/chat/system-prompt.txt
containing a compact JSON catalog wrapped in instructions.

The catalog is trimmed (heavyweight reference evidence dropped) to keep the
prompt small enough to fit comfortably in the 32K-context window of a
sub-1B local model. Roughly 12-18 KB after trimming.
"""
import argparse
import json
import sys
from pathlib import Path


SYSTEM_HEADER_TEMPLATE = """You are the fz documentation assistant. fz is an error-bounded scientific data compression toolkit (SZ3 / ZFP / SPERR / MGARD compositions). Answer ONLY using the JSON catalog below. Recommend specific ALGO_* values when the user describes their data. If the catalog does not cover the question, say so - do not invent modules, papers, or flags. Be concise (3-6 sentences). When you recommend an algorithm, name it (e.g. ALGO_BIOMD) and quote the relevant "best_for" snippet. End every response with: AI-generated - verify against the API docs.

Valid ALGO ids: {algo_allowlist}. Refuse to mention any other ALGO_*.

CATALOG:
"""


def trim_entry(e: dict) -> dict:
    """Return a trimmed copy of an algo/module entry suitable for the prompt.

    Drops verbose evidence prose from references (keeps title+venue+year),
    keeps the load-bearing best_for / pipeline / keywords / data_examples /
    error_modes / supported_dims / supported_dtypes fields.
    """
    out = {
        "id": e["id"],
        "name": e.get("name", e["id"]),
        "category": e.get("category", "module"),
        "best_for": e.get("best_for", ""),
    }
    pipe = e.get("pipeline")
    if pipe:
        out["pipeline"] = pipe
    if e.get("pipeline_notes"):
        out["pipeline_notes"] = e["pipeline_notes"]
    for f in ("supported_dims", "supported_dtypes", "error_modes",
              "keywords", "data_examples", "used_by", "canonical_terms"):
        v = e.get(f)
        if v:
            out[f] = v
    refs = e.get("references") or []
    if refs:
        out["references"] = [
            {k: r.get(k) for k in ("title", "venue", "year") if r.get(k)}
            for r in refs
        ]
    return out


def build_prompt(modules_json: Path) -> str:
    with open(modules_json) as f:
        data = json.load(f)
    trimmed = {
        "algos": [trim_entry(e) for e in data.get("algos", [])],
        "modules": [trim_entry(e) for e in data.get("modules", [])],
    }
    # Compact JSON: no whitespace, but keep keys sorted for diff-friendliness.
    body = json.dumps(trimmed, separators=(",", ":"), sort_keys=True)
    # Build the dynamic ALGO allowlist line so it stays in sync with modules.json.
    algo_ids = [
        a["id"] for a in data.get("algos", [])
        if isinstance(a.get("id"), str) and a["id"].startswith("ALGO_")
    ]
    allowlist = ", ".join(algo_ids) if algo_ids else "(none defined)"
    header = SYSTEM_HEADER_TEMPLATE.format(algo_allowlist=allowlist)
    return header + body + "\n"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("modules_json")
    p.add_argument("output_path",
                   help="Path to write system-prompt.txt "
                        "(typically docs/html/assets/chat/system-prompt.txt)")
    args = p.parse_args()

    src = Path(args.modules_json).resolve()
    dst = Path(args.output_path).resolve()
    if not src.exists():
        sys.exit(f"ERROR: modules.json not found: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    prompt = build_prompt(src)
    dst.write_text(prompt)
    kb = len(prompt.encode("utf-8")) / 1024.0
    print(f"Wrote {dst} ({kb:.1f} KB, {len(prompt)} chars)")


if __name__ == "__main__":
    main()
