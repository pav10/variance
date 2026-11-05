#!/bin/bash
# varspec.sh D V0 OUT [NPROC]: all real-rooted monic integer polynomials of degree D whose
# roots have variance <= V0, with x -> +-x + k normalised (a1 in [-floor(D/2), 0]).
# Uses enum6.c in "all" mode (no filters).  Output lines: a1 a2 ... aD
D=$1; V0=$2; OUT=$3; NP=${4:-1}
BIN=${ENUMBIN:-/tmp/claude-0/bin}
SMAX=$(python3 -c "print($D*$V0 + 1e-9)")
python3 - $D $SMAX > $OUT.jobs <<'PY'
import sys
d=int(sys.argv[1]); smax=float(sys.argv[2])
for a1 in range(0, -(d//2)-1, -1):
    for a2 in range(-20*d, 10*d):
        S=(d-1)*a1*a1/d-2*a2
        if 0 < S <= smax: print(a1, a2)
PY
: > $OUT
cat $OUT.jobs | xargs -P $NP -n 2 sh -c "$BIN/enum_d$D \$0 \$1 $SMAX all 2>/dev/null" >> $OUT
rm -f $OUT.jobs
