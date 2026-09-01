#!/usr/bin/env python3
"""Check that every SZ3 header compiles on its own.

A header that only works because some other header happened to be included first is a
latent portability break: it compiles until a translation unit changes its include order,
or until a standard library stops providing something transitively.

Usage:
    python3 tools/test/check_headers.py [--build-dir DIR] [--jobs N] [--compiler CXX]

The build directory supplies the generated SZ3/version.hpp and the fetched Eigen; it
defaults to the first of build/, build-release/ that has one.
"""

import argparse
import concurrent.futures
import os
import pathlib
import subprocess
import sys
import tempfile

REPO = pathlib.Path(__file__).resolve().parents[2]
INCLUDE = REPO / "include"
SKIP_DIRS = {"thirdparty", "testing"}  # testing/ needs GoogleTest, which the library does not depend on


def find_include_dirs(build_dir):
    """Repo headers, the build tree's generated headers, and the third-party dirs it configured."""
    dirs = [INCLUDE]
    if build_dir:
        b = pathlib.Path(build_dir)
        if (b / "include").is_dir():
            dirs.append(b / "include")
        # The build tree records where CMake found Eigen and Zstd; reuse exactly those.
        cache = b / "CMakeCache.txt"
        if cache.is_file():
            for line in cache.read_text(errors="replace").splitlines():
                key, _, value = line.partition("=")
                if key.startswith("Eigen3_DIR"):
                    d = pathlib.Path(value.strip())
                    dirs += [d, d.parent / "eigen-src"]
                elif key.startswith(("ZSTD_CFLAGS", "ZSTD_INCLUDE_DIR")):
                    for token in value.split():
                        dirs.append(pathlib.Path(token[2:] if token.startswith("-I") else token))
        dirs += sorted(b.glob("_deps/*eigen*"))
    dirs += [pathlib.Path(p) for p in ("/usr/include/eigen3", "/usr/local/include/eigen3",
                                       "/opt/homebrew/include/eigen3", "/opt/homebrew/include",
                                       "/usr/local/include")]
    seen, out = set(), []
    for d in dirs:
        if d.is_dir() and d not in seen:
            seen.add(d)
            out.append(d)
    return out


def default_build_dir():
    for name in ("build", "build-release"):
        candidate = REPO / name
        if (candidate / "include" / "SZ3" / "version.hpp").is_file():
            return candidate
    return None


def headers():
    for path in sorted(INCLUDE.rglob("*.hpp")):
        rel = path.relative_to(INCLUDE)
        if SKIP_DIRS & set(rel.parts):
            continue
        yield rel


def check(rel, compiler, inc_flags, tmpdir):
    src = pathlib.Path(tmpdir) / (str(rel).replace(os.sep, "_") + ".cpp")
    src.write_text('#include "%s"\n' % rel.as_posix())
    proc = subprocess.run(
        [compiler, "-std=c++17", "-fsyntax-only"] + inc_flags + [str(src)],
        capture_output=True,
        text=True,
    )
    return rel, proc.returncode, proc.stderr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", default=None)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--compiler", default=os.environ.get("CXX", "c++"))
    args = ap.parse_args()

    build_dir = args.build_dir or default_build_dir()
    if build_dir is None:
        print("error: no build directory with a generated SZ3/version.hpp; configure one first")
        return 2
    inc_flags = ["-I%s" % d for d in find_include_dirs(build_dir)]

    all_headers = list(headers())
    failures = []
    with tempfile.TemporaryDirectory() as tmpdir:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(check, rel, args.compiler, inc_flags, tmpdir) for rel in all_headers]
            for future in concurrent.futures.as_completed(futures):
                rel, code, err = future.result()
                if code != 0:
                    failures.append((rel, err))

    for rel, err in sorted(failures):
        print("NOT SELF-CONTAINED: %s" % rel.as_posix())
        for line in err.splitlines():
            if "error:" in line:
                print("    " + line.strip())

    print("%d headers checked, %d not self-contained" % (len(all_headers), len(failures)))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
