#!/bin/bash
# The tokens of a link line that link an SZ3 library, one per line.
#
#   tools/test/sz3LinkTokens.sh <link.txt>
#   tools/test/sz3LinkTokens.sh --self-test
#
# Not "the tokens that contain sz3". On a GitHub runner the checkout is at
# /home/runner/work/SZ3/SZ3, because the repository is named SZ3, so the name of the project is
# in every absolute path on the line and -Wl,-rpath-link,<builddir>/lib matches it while linking
# nothing. A token links an SZ3 library if it is -l<name with sz3 in it>, or a path whose *file
# name* is an SZ3 library. -L, -I and -Wl,... name places to look, not libraries to link.
set -u

sz3_link_tokens() {
    tr '[:space:]' '\n' | awk '
        { t = $0; gsub(/^"+|"+$/, "", t) }
        t == ""   { next }
        t ~ /^-l/ { if (tolower(t) ~ /sz3/) print t; next }
        t ~ /^-/  { next }
        { n = tolower(t); sub(/.*\//, "", n)
          if (n ~ /sz3/ && n ~ /\.(so|a|dylib|dll|lib)(\.[0-9]+)*$/) print t }'
}

# Every check that reads this file is only as good as the line above, so prove it both ways.
# The tokens are real ones, off the link lines in the GROMACS workflow's own runs.
self_test() {
    bad=0
    while IFS='|' read -r tok want; do
        [ -n "$tok" ] || continue
        got=no
        [ -n "$(printf '%s\n' "$tok" | sz3_link_tokens)" ] && got=yes
        if [ "$got" != "$want" ]; then
            echo "sz3LinkTokens self-test: '$tok' -> $got, expected $want"
            bad=1
        fi
    done <<'CASES'
/home/runner/work/SZ3/SZ3/_prefix/lib/libhdf5sz3.so|yes
lib/libhdf5sz3.so|yes
lib/libhdf5sz3.so.3.3.3|yes
libSZ3c.a|yes
-lhdf5sz3|yes
-l:libhdf5sz3.so|yes
-Wl,-rpath-link,/home/runner/work/SZ3/SZ3/_m/found/b-OFF/lib|no
-Wl,-rpath,"$ORIGIN/../lib:/home/runner/work/SZ3/SZ3/_prefix/lib"|no
-L/home/runner/work/SZ3/SZ3/_prefix/lib|no
-I/home/runner/work/SZ3/SZ3/include|no
/home/runner/work/SZ3/SZ3/gbuild/lib/libgromacs.so.10|no
/home/runner/work/SZ3/SZ3/_m/found/b-OFF/lib|no
"CMakeFiles/h5md-test.dir/h5md.cpp.o"|no
/home/runner/micromamba/envs/h5/lib/libhdf5.so|no
CASES
    [ "$bad" = 0 ] || return 1
    echo "sz3LinkTokens self-test: ok"
}

case "${1:-}" in
--self-test)
    self_test
    ;;
"" | -h | --help)
    sed -n '2,12p' "$0"
    exit 1
    ;;
*)
    [ -f "$1" ] || { echo "no such link line: $1" >&2; exit 1; }
    sz3_link_tokens < "$1"
    ;;
esac
