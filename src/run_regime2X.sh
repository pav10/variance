#!/bin/bash
# regime 2, both rungs regular: degree-6 Case A enumeration with F1, F2 (downspread) but
# without the quadratic filter F3 (not valid for G when nu_2 < 2 nu).  Output: one poly per line.
NP=${1:-2}; OUT=results/regime2X; BIN=/tmp/claude-0/bin
python3 - > $OUT/jobs.txt <<'PY'
import math
d = 6
smax = d * ((d + 4) + math.sqrt((d + 4) ** 2 - 16)) / 4
jobs = []
for a1 in range(0, -(d // 2) - 1, -1):
    for a2 in range(-60, 30):
        S = (d - 1) * a1 * a1 / d - 2 * a2
        if 0 < S <= smax + 1e-9: jobs.append((S, a1, a2))
jobs.sort(reverse=True)
for S, a1, a2 in jobs: print(a1, a2)
PY
cat $OUT/jobs.txt | ENUM_NOF3=1 xargs -P $NP -n 2 sh -c "$BIN/enum6g_d6 \$0 \$1 > $OUT/out_\$0_\$1.txt 2> $OUT/log_\$0_\$1.txt"
cat $OUT/out_*.txt > $OUT/raw.txt
