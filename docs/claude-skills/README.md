# fz skills for Claude Code

A small bundle of [Claude Code skills](https://docs.claude.com/en/docs/claude-code/skills) that teach Claude how to navigate, extend, and benchmark this repository.

## What's in here

| Skill | When Claude invokes it |
|---|---|
| `fz-overview` | "what algorithms does fz support?", "fz architecture", "what modules are in fz?" |
| `fz-compose-pipeline` | "compose a custom pipeline", "swap in a different encoder", "build a custom compressor from fz modules" |
| `fz-add-module` | "add a new encoder/quantizer/decomposition/lossless module" |
| `fz-add-algorithm` | "add a new ALGO_*", "register a new algorithm in the dispatcher" |
| `fz-bench-multibound` | "benchmark across multiple error bounds", "compare two algorithms head-to-head" |

## Install

Claude Code looks for skills in two locations:

- **`<project>/.claude/skills/`** — project-scoped, only loaded when Claude is run inside this repo
- **`~/.claude/skills/`** — user-scoped, loaded everywhere

Pick whichever scope fits. Both are gitignored by default; that's why this repo keeps the source-of-truth copies under `docs/claude-skills/` instead.

### Project-scoped (recommended for fz contributors)

From the repo root:

```bash
mkdir -p .claude/skills
for d in docs/claude-skills/*/; do ln -s "../../$d" ".claude/skills/$(basename "$d")"; done
```

(Or copy instead of symlink if you prefer to detach.)

### User-scoped (use these skills in any directory)

```bash
mkdir -p ~/.claude/skills
cp -r docs/claude-skills/* ~/.claude/skills/
```

## Verify

After installing, restart your Claude Code session (or run `/skills`) and confirm the five skills appear. Trigger one with a prompt like "what algorithms does fz support?" — Claude should cite `fz-overview` and answer from it.

## Editing

The canonical copies live under `docs/claude-skills/`. Edit there and re-symlink/copy. PRs against this directory are welcome.
