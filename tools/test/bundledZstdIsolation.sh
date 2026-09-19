#!/bin/bash
# What a build with SZ3_USE_BUNDLED_ZSTD=ON is allowed to hand the rest of the machine: the
# vendored Zstd has to reach SZ3's consumers as a linkable target and as nothing else -- no header
# in the prefix, none reachable through the export, no symbol in a shared library.
#
#   tools/test/bundledZstdIsolation.sh <build-dir> <install-prefix>
#
# Assert which Zstd was used and what escaped, never just that the build exited 0. No set -e: every
# check below is explicit, and the total at the bottom catches a section that stopped early.
set -u
BUILD=$(cd "$1" && pwd)
PREFIX=$(cd "$2" && pwd)
LIBDIR=$PREFIX/lib
[ -d "$LIBDIR" ] || LIBDIR=$PREFIX/lib64
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
cd "$WORK" || exit 1

pass=0; fail=0; skip=0
ok()   { echo "PASS  $1"; pass=$((pass+1)); }
bad()  { echo "FAIL  $1"; shift; for l in "$@"; do echo "        $l"; done; fail=$((fail+1)); }
# Never counted as a pass: a suite that reports more checks than it ran is worse than no suite.
skipped() { echo "SKIP  $1"; skip=$((skip+1)); }

SZ3BIN=$BUILD/tools/sz3/sz3
[ -x "$SZ3BIN" ] || SZ3BIN=$BUILD/tools/sz3/sz3.exe
# Required, not skipped: a missing binary must never look like a passing check.
if [ ! -x "$SZ3BIN" ]; then
    echo "no sz3 binary under $BUILD/tools/sz3, so there is nothing to ask which Zstd it linked"
    exit 1
fi
# The one place a Mach-O build differs; everything else here reads files.
linked_libs() {
    if command -v ldd > /dev/null 2>&1; then ldd "$1" 2>&1
    elif command -v otool > /dev/null 2>&1; then otool -L "$1" 2>&1
    else echo ""; fi
}

echo "=== bundled Zstd: build $BUILD, prefix $PREFIX ==="

# ---------------------------------------------------------------- 1. it really is the bundled one
if [ -e "$LIBDIR/libsz3_zstd.a" ] || [ -e "$LIBDIR/sz3_zstd.lib" ]; then
    ok "bundled-zstd-archive-installed"
else
    bad "bundled-zstd-archive-installed" "no libsz3_zstd.a under $LIBDIR"
fi
linked_libs "$SZ3BIN" > linked.log
if grep -q 'libzstd\.so\|libzstd\.[0-9]*\.dylib\|libzstd\.dylib' linked.log; then
    bad "no-system-libzstd-linked" "asked for the bundled Zstd and linked the system one anyway" \
        "$(grep 'libzstd' linked.log)"
else
    ok "no-system-libzstd-linked"
fi

# ---------------------------------------------------------------- 2. what the export offers
TARGETS=$LIBDIR/cmake/SZ3/SZ3Targets.cmake
if [ -f "$TARGETS" ] && grep -q 'SZ3::zstd' "$TARGETS"; then
    ok "export-names-sz3-zstd"
else
    bad "export-names-sz3-zstd" "$TARGETS does not name SZ3::zstd, so a consumer has nothing to link"
fi
# An absolute path into this machine's build tree pins the export to the machine that produced it.
grep -h 'libsz3_zstd' "$LIBDIR"/cmake/SZ3/SZ3Targets*.cmake 2>/dev/null > vendored.log
if [ ! -s vendored.log ]; then
    # Not a pass: SZ3::zstd is an archive the export has to locate by path, so naming it nowhere
    # means a consumer resolves the name itself and this check had nothing to look at.
    bad "export-references-the-vendored-zstd-relocatably" \
        "no SZ3Targets*.cmake under $LIBDIR/cmake/SZ3 names libsz3_zstd at all"
elif grep -qv '_IMPORT_PREFIX' vendored.log; then
    bad "export-references-the-vendored-zstd-relocatably" \
        "referenced by an absolute path, so the export cannot move" "$(cat vendored.log)"
else
    ok "export-references-the-vendored-zstd-relocatably"
fi

# ---------------------------------------------------------------- 3. nothing leaked into the prefix
# These are the names the system Zstd's own package owns. Installing any of them makes SZ3
# uninstallable beside it, which is what Debian and Fedora reject a package for.
leaked=
for f in include/zstd.h include/zdict.h include/zbuff.h include/zstd_errors.h \
         lib/libzstd.a lib/libzstd.so include/SZ3/bundled_zstd/zstd.h; do
    [ -e "$PREFIX/$f" ] && leaked="$leaked $f"
done
if [ -n "$leaked" ]; then
    bad "no-zstd-file-the-system-package-owns" "installed:$leaked"
else
    ok "no-zstd-file-the-system-package-owns"
fi
find "$PREFIX" -name 'zstd.h' > anyheader.log 2>/dev/null
if [ -s anyheader.log ]; then
    bad "no-zstd-h-anywhere-in-the-prefix" "$(cat anyheader.log)"
else
    ok "no-zstd-h-anywhere-in-the-prefix"
fi

# A consumer compiles with every exported include directory on its command line, so a zstd.h
# sitting in one of them shadows the system Zstd's header for the whole translation unit.
grep -ho 'INTERFACE_INCLUDE_DIRECTORIES "[^"]*"' "$LIBDIR"/cmake/SZ3/*.cmake 2>/dev/null \
  | sed 's/.*"\(.*\)"/\1/' | tr ';' '\n' | sed "s|\${_IMPORT_PREFIX}|$PREFIX|" \
  | sed '/^$/d' | sort -u > incdirs.log
# Counted first: with no include directories to walk, the check below has nothing to look inside
# and would pass on an export that declares none.
if [ -s incdirs.log ]; then
    ok "export-declares-include-directories"
    sed 's/^/        /' incdirs.log
else
    bad "export-declares-include-directories" "no INTERFACE_INCLUDE_DIRECTORIES under $LIBDIR/cmake/SZ3"
fi
reachable=
while read -r d; do
    [ -e "$d/zstd.h" ] && reachable="$reachable $d"
done < incdirs.log
if [ -n "$reachable" ]; then
    bad "no-zstd-h-reachable-from-an-exported-include-dir" "a zstd.h in:$reachable"
else
    ok "no-zstd-h-reachable-from-an-exported-include-dir"
fi

# ---------------------------------------------------------------- 4. nothing leaked into the ABI
# A shared library that re-exports ZSTD_* interposes on every other Zstd in the process.
find "$LIBDIR" -maxdepth 1 \( -name '*.so' -o -name '*.so.*' -o -name '*.dylib' \) > shared.log 2>/dev/null
if ! command -v nm > /dev/null 2>&1; then
    skipped "no-shared-library-exports-zstd-symbols (no nm)"
elif [ ! -s shared.log ]; then
    # A named skip, not a pass: BUILD_SHARED_LIBS=OFF installs archives, and an archive carries the
    # symbols without exporting them. Counted, so the total below still adds up.
    skipped "no-shared-library-exports-zstd-symbols ($LIBDIR holds no shared library)"
else
    sed 's/^/        /' shared.log
    exporters=; readable=0
    while read -r lib; do
        nm -D --defined-only "$lib" > syms.log 2>/dev/null
        # Counted, because nm reporting nothing for every library reads exactly like a clean result.
        [ -s syms.log ] && readable=$((readable+1))
        grep -q ' T \(ZSTD_\|ZDICT_\|ZBUFF_\)' syms.log && exporters="$exporters $lib"
    done < shared.log
    if [ -n "$exporters" ]; then
        bad "no-shared-library-exports-zstd-symbols" "exports Zstd's symbols:$exporters"
    elif [ "$readable" = 0 ]; then
        bad "no-shared-library-exports-zstd-symbols" \
            "nm read a dynamic symbol table out of none of them, so nothing was inspected"
    else
        ok "no-shared-library-exports-zstd-symbols"
    fi
fi

echo
echo "  $pass passed, $fail failed, $skip skipped"

# Raise this with the check it comes with. A section that stops early otherwise shows only as a
# smaller number at the bottom that nobody compares.
EXPECTED=9
ran=$((pass + fail + skip))
if [ "$ran" -ne "$EXPECTED" ]; then
    echo "  the suite accounted for $ran checks, not $EXPECTED"
    exit 1
fi
exit $((fail > 0))
