#!/usr/bin/env python3
"""Check docs/site/modules.json against the tree.

The composition constraints in that file are the machine-readable form of facts the C++ tests
already enforce. This keeps the two from drifting: it verifies every entry points at a header that
exists, that the ids it cross-references resolve, and that the constraint vocabulary stays closed.
"""

import json
import os
import pathlib
import sys

REPO = pathlib.Path(__file__).resolve().parents[3]
MODULES = REPO / "docs" / "site" / "modules.json"

CATEGORIES = {"decomposition", "quantizer", "encoder", "lossless", "preprocessor", "predictor", "utility"}
COMPOSITION_KEYS = {
    "requires_state_num",
    "requires_construction_args",
    "requires_encoder",
    "pairs_with",
    "avoid_pairing",
    "out_range",
    "error_bound_regime",
    "expands_payload",
    "notes",
}
OUT_RANGE_VALUES = {"none (reports 0)", "signed (does not start at 0)"}


def main():
    doc = json.loads(MODULES.read_text())
    mods = doc["modules"]
    ids = {m["id"] for m in mods}
    errors = []

    if len(ids) != len(mods):
        errors.append("duplicate module ids")

    for m in mods:
        mid = m["id"]
        for field in ("id", "name", "category", "header", "best_for"):
            if not m.get(field):
                errors.append("%s: missing %s" % (mid, field))
        if m.get("category") not in CATEGORIES:
            errors.append("%s: unknown category %r" % (mid, m.get("category")))
        if not (REPO / m["header"]).is_file():
            errors.append("%s: header does not exist: %s" % (mid, m["header"]))

        comp = m.get("composition")
        if comp is None:
            continue
        unknown = set(comp) - COMPOSITION_KEYS
        if unknown:
            errors.append("%s: unknown composition keys %s" % (mid, sorted(unknown)))
        if "out_range" in comp and comp["out_range"] not in OUT_RANGE_VALUES:
            errors.append("%s: unknown out_range %r" % (mid, comp["out_range"]))
        for key in ("requires_encoder", "pairs_with"):
            for ref in comp.get(key, []):
                if ref not in ids:
                    errors.append("%s: %s references unknown module %r" % (mid, key, ref))
        for entry in comp.get("avoid_pairing", []):
            if entry.get("with") not in ids:
                errors.append("%s: avoid_pairing references unknown module %r" % (mid, entry.get("with")))
            if not entry.get("why"):
                errors.append("%s: avoid_pairing entry has no reason" % mid)

    for e in errors:
        print("modules.json: " + e)
    print("%d modules checked, %d with composition constraints, %d problems"
          % (len(mods), sum(1 for m in mods if "composition" in m), len(errors)))
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
