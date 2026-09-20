#!/bin/bash
# The error sz3 -a reports, against the bound it was asked for.
#
#   tools/test/checkMaxError.sh <sz3 -a output> <tolerance>
#
# The parse is its own assertion: sed exits 0 when it matches nothing, so a renamed or missing line
# leaves an empty string that compares below the bound and reads as a compressor that stayed in it.
set -u
LOG=$1
TOL=$2

if [ ! -f "$LOG" ]; then
    echo "FAIL  max-absolute-error-was-reported"
    echo "        no such file: $LOG"
    exit 1
fi
if [ -z "$TOL" ]; then
    echo "FAIL  max-absolute-error-within-tolerance"
    echo "        no tolerance given, so there is nothing to compare against"
    exit 1
fi

ERR=$(grep -E 'Max absolute error =' "$LOG" \
      | sed -E 's/.*Max absolute error = *([0-9]+(\.[0-9]+)?([eE][+-]?[0-9]+)?).*/\1/')
if [ -z "$ERR" ]; then
    echo "FAIL  max-absolute-error-was-reported"
    echo "        no 'Max absolute error =' line in $LOG"
    echo "        last lines:"
    tail -10 "$LOG" | sed 's/^/          /'
    exit 1
fi
echo "PASS  max-absolute-error-was-reported"

if awk -v err="$ERR" -v tol="$TOL" 'BEGIN { exit (err > tol) }'; then
    echo "PASS  max-absolute-error-within-tolerance"
    echo "        $ERR <= $TOL"
else
    echo "FAIL  max-absolute-error-within-tolerance"
    echo "        $ERR exceeds the requested bound $TOL"
    exit 1
fi
