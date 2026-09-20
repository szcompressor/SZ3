#!/usr/bin/env python3
"""Every env key a workflow interpolates has to be declared in a scope that reaches it.

GitHub interpolates an undeclared key to the empty string, so a step keeps running with an
argument missing. That is how splitting one workflow into three left `sz3 -3  -M  ` behind: the
compressor printed its usage text and exited non-zero, which reads as a compression failure.

Scope counts for as much as spelling. A key declared in one job does not exist in another, a
step's own env: reaches that step alone, and a $GITHUB_ENV write reaches only the later steps of
the same job. Checking every reference against the union of every env: block in the file would
pass the one workflow this script exists to catch, so each scope is tracked separately.

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
MAPPING_KEY = re.compile(r"^(\s*)([A-Za-z_][A-Za-z0-9_-]*):")
STEPS_BLOCK = re.compile(r"^(\s*)steps:\s*(#.*)?$")
LIST_ITEM = re.compile(r"^(\s*)-(\s|$)")

EMPTY_VALUES = ("", "''", '""')


def line_scopes(lines):
    """For each line, the (job, step) it sits in -- the two scopes narrower than the workflow.

    `job` is the job id, or None for a line outside `jobs:`. `step` is the step's index within
    that job, or None for a line that belongs to the job itself (runs-on:, strategy:, ...).
    """
    scopes = []
    in_jobs = False
    job = job_indent = None
    step = step_indent = steps_indent = None
    for line in lines:
        stripped = line.rstrip("\n")
        body = stripped.strip()
        if body and not body.startswith("#"):
            indent = len(stripped) - len(stripped.lstrip())
            if indent == 0 and MAPPING_KEY.match(stripped):
                # A new top-level key closes every job and step under the previous one.
                in_jobs = body.split(":", 1)[0] == "jobs"
                job = job_indent = None
                step = step_indent = steps_indent = None
            elif in_jobs and job_indent is None:
                # The first key under `jobs:` fixes the indent the job ids sit at.
                match = MAPPING_KEY.match(stripped)
                if match:
                    job_indent = indent
                    job = match.group(2)
            elif in_jobs and indent == job_indent and MAPPING_KEY.match(stripped):
                job = MAPPING_KEY.match(stripped).group(2)
                step = step_indent = steps_indent = None
            elif in_jobs and job is not None:
                item = LIST_ITEM.match(stripped)
                if steps_indent is not None and indent <= steps_indent and not (
                    item and indent == step_indent
                ):
                    # A key at or above `steps:` is the job speaking again, not another step.
                    step = step_indent = steps_indent = None
                if steps_indent is None:
                    match = STEPS_BLOCK.match(stripped)
                    if match:
                        steps_indent = len(match.group(1))
                elif item and (step_indent is None or indent == step_indent):
                    step_indent = indent
                    step = 0 if step is None else step + 1
        scopes.append((job, step))
    return scopes


def declarations(lines):
    """Each key an env: mapping declares, as (index of the env: line, name, value)."""
    block_indent = key_indent = block_line = None
    for index, line in enumerate(lines):
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
                yield block_line, match.group(2), match.group(3).strip()
            continue
        block = ENV_BLOCK.match(stripped)
        if block:
            block_indent = len(block.group(1))
            key_indent = None
            block_line = index


def runtime_declarations(lines):
    """Names written to $GITHUB_ENV, as (line index, name). Only later steps of the job see them."""
    for index, line in enumerate(lines):
        if "GITHUB_ENV" not in line:
            continue
        for match in RUNTIME_KEY.finditer(line):
            if match.group(1) != "GITHUB_ENV":
                yield index, match.group(1)


def check(path):
    text = path.read_text()
    lines = text.splitlines()
    scopes = line_scopes(lines)

    workflow = {}
    jobs = {}  # job -> {name: value}
    steps = {}  # (job, step) -> {name: value}
    for index, name, value in declarations(lines):
        job, step = scopes[index]
        if job is None:
            target = workflow
        elif step is None:
            target = jobs.setdefault(job, {})
        else:
            target = steps.setdefault((job, step), {})
        target[name] = value

    # job -> name -> the earliest step that writes it; a write is visible from the next step on.
    runtime = {}
    for index, name in runtime_declarations(lines):
        job, step = scopes[index]
        if job is None or step is None:
            continue
        written = runtime.setdefault(job, {})
        written[name] = min(step, written.get(name, step))

    problems = []
    for match in REFERENCE.finditer(text):
        name = match.group(1)
        index = text.count("\n", 0, match.start())
        job, step = scopes[index]
        # GitHub's resolution order, innermost first. Anything outside this chain is a different
        # scope.
        chain = []
        if step is not None:
            chain.append(steps.get((job, step), {}))
            if runtime.get(job, {}).get(name, step) < step:
                chain.append({name: "<$GITHUB_ENV>"})
        if job is not None:
            chain.append(jobs.get(job, {}))
        chain.append(workflow)

        value = next((scope[name] for scope in chain if name in scope), None)
        if value is None:
            problems.append(f"{path}:{index + 1}: env.{name} is not declared in any scope that reaches it")
        elif value in EMPTY_VALUES:
            problems.append(f"{path}:{index + 1}: env.{name} is declared with an empty value")
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
