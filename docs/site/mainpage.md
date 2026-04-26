FZ Module Library {#mainpage}
=================

**FZ** is a modular toolkit for error-bounded scientific data compression.
The goal is to bring **all important lossy compression modules** into one
C++17 header-only library — including SZ3 (interpolation/Lorenzo), ZFP
(block transform), SPERR (wavelet/SPECK), MGARD (multigrid), and a growing
catalog of additional algorithms — so developers can compose new pipelines
from a shared pool of modules rather than picking one monolithic compressor.

@htmlonly
<div style="display: flex; gap: 18px; flex-wrap: wrap; margin: 1.5em 0 2em 0;">
  <a href="modules/index.html"
     style="flex: 1 1 320px; min-width: 280px; max-width: 460px;
            display: block; padding: 18px 22px;
            background: #f6f8fa; color: #283a5d; border: 1px solid #c4cfe5;
            border-radius: 8px; text-decoration: none;">
    <div style="font-size: 0.78em; letter-spacing: 0.08em; text-transform: uppercase; color: #6a7da3;">Browse</div>
    <div style="font-size: 1.4em; font-weight: 600; margin-top: 4px;">Algorithm Catalog &rarr;</div>
    <div style="margin-top: 8px; font-size: 0.92em;">
      Algorithms (SZ3, ZFP, SPERR, MGARD, MDZ, &hellip;) and reusable pipeline modules,
      indexed by data type, application domain, and technique.
    </div>
  </a>
  <a href="annotated.html"
     style="flex: 1 1 280px; min-width: 240px; max-width: 380px;
            display: block; padding: 18px 22px;
            background: #f6f8fa; color: #283a5d; border: 1px solid #c4cfe5;
            border-radius: 8px; text-decoration: none;">
    <div style="font-size: 0.78em; letter-spacing: 0.08em; text-transform: uppercase; color: #6a7da3;">Browse</div>
    <div style="font-size: 1.4em; font-weight: 600; margin-top: 4px;">API Reference &rarr;</div>
    <div style="margin-top: 8px; font-size: 0.92em;">
      Doxygen-generated class, namespace, and file reference for everything
      under <code>include/SZ3/</code>.
    </div>
  </a>
</div>
@endhtmlonly

Coding with AI agents
---------------------

If you're using Claude Code (or another agent-driven IDE assistant), point it at
[`CLAUDE.md`](https://github.com/szcompressor/SZ3/blob/fz/CLAUDE.md) at the
repo root for an LLM-tailored overview of the pipeline, the module catalog, and
the build. The repo also ships installable Claude skills under
[`docs/claude-skills/`](https://github.com/szcompressor/SZ3/tree/fz/docs/claude-skills)
covering the common workflows:

- `fz-overview` &mdash; pipeline + module map
- `fz-compose-pipeline` &mdash; wire existing modules into a new composition
- `fz-add-module` &mdash; add a brand-new Decomposition / Quantizer / Encoder / Lossless
- `fz-add-algorithm` &mdash; add a brand-new `ALGO_*` enum + dispatch + wiring
- `fz-bench-multibound` &mdash; multi-bound benchmark recipe

Symlink (or copy) any of those into `~/.claude/skills/` and your agent gets the
matching workflow on tap. Or just open the chat panel below and ask in natural
language.

Links
-----

[fzframework.org](https://fzframework.org/) &mdash; source repository, issue tracker,
and Python / HDF5 / Rust / Fortran bindings.
