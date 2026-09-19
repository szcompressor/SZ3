#!/usr/bin/env python3
"""Every env key a workflow interpolates has to be declared in that same workflow.

GitHub interpolates an undeclared key to the empty string, so a step keeps running with an
argument missing. That is how splitting one workflow into three left `sz3 -3  -M  ` behind: the
compressor printed its usage text and exited non-zero, which reads as a compression failure.

    python3 tools/test/check_workflow_env.py [--workflow-dir .github/workflows]

Reads the workflow files as text rather than YAML so it needs nothing installed.
"""

import argparse
import pathlib
import re
import sys

# `${{ env.NAME }}`, with the whitespace GitHub allows on either side of the name.
REFERENCE = re.compile(r"\$\{\{\s*env\.([A-Za-z_][A-Za-z0-9_-]*)\s*\}\}")
# A key inside an `env:` mapping: `NAME: value`, where value may be empty or continue on later
# lines. An empty value is a finding too -- it reaches a program exactly as an undeclared key does.
ENV_KEY = re.compile(r"^(\s*)([A-Za-z_][A-Za-z0-9_-]*):(.*)$")
ENV_BLOCK = re.compile(r"^(\s*)env:\s*(#.*)?$")
# `echo "NAME=..." >> $GITHUB_ENV` and the PowerShell `"NAME=..." | Out-File ... GITHUB_ENV`, which
# declare a key for later steps without an env: block.
RUNTIME_KEY = re.compile(r"[\"']?([A-Za-z_][A-Za-z0-9_-]*)=")


def declared_keys(lines):
    """Names an env: mapping declares, and the subset of them whose value is empty."""
    names, empty = set(), set()
    block_indent = None
    key_indent = None
    for line in lines:
        stripped = line.rstrip("\n")
        if not stripped.strip() or stripped.lstrip().startswith("#"):
            continue
        indent = len(stripped) - len(stripped.lstrip())
        if block_indent is not None and indent <= block_indent:
            block_indent = key_indent = None
        if block_indent is not None:
            match = ENV_KEY.match(stripped)
            if match and (key_indent is None or len(match.group(1)) == key_indent):
                key_indent = len(match.group(1))
                names.add(match.group(2))
                value = match.group(3).strip()
                if value in ("", "''", '""'):
                    empty.add(match.group(2))
            continue
        block = ENV_BLOCK.match(stripped)
        if block:
            block_indent = len(block.group(1))
            key_indent = None
    return names, empty


def runtime_keys(text):
    """Names written to $GITHUB_ENV, which later steps may read."""
    names = set()
    for line in text.splitlines():
        if "GITHUB_ENV" not in line:
            continue
        names.update(match.group(1) for match in RUNTIME_KEY.finditer(line))
    return names - {"GITHUB_ENV"}


def check(path):
    text = path.read_text()
    declared, empty = declared_keys(text.splitlines())
    declared |= runtime_keys(text)
    problems = []
    for match in REFERENCE.finditer(text):
        name = match.group(1)
        line = text.count("\n", 0, match.start()) + 1
        if name not in declared:
            problems.append(f"{path}:{line}: env.{name} is used but never declared")
        elif name in empty:
            problems.append(f"{path}:{line}: env.{name} is declared with an empty value")
    return problems


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-dir", default=".github/workflows", type=pathlib.Path)
    args = parser.parse_args()

    paths = sorted(args.workflow_dir.glob("*.yml")) + sorted(args.workflow_dir.glob("*.yaml"))
    if not paths:
        print(f"no workflows under {args.workflow_dir}", file=sys.stderr)
        return 1

    problems = []
    for path in paths:
        problems += check(path)
    for problem in problems:
        print(problem, file=sys.stderr)
    print(f"checked {len(paths)} workflows, {len(problems)} undeclared or empty env references")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
