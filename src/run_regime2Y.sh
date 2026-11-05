#!/bin/bash
# sextic Case A, a half-square rung: enum6g.c in mode GYL over all classes a1 in [-3,3]
# (a1 > 0 are the reflections x -> -x, covering the upper rung), nu <= 15.75 (S <= 94.5).
NP=${1:-3}; OUT=results/regime2Y; BIN=/tmp/claude-0/bin
python3 - > $OUT/jobs.txt <<'PY'
d = 6; smax = 94.5
jobs = []
for a1 in range(-3, 4):
    for a2 in range(-60, 30):
        S = (d - 1) * a1 * a1 / d - 2 * a2
        if 0 < S <= smax + 1e-9: jobs.append((S, a1, a2))
jobs.sort(reverse=True)
for S, a1, a2 in jobs: print(a1, a2)
PY
cat $OUT/jobs.txt | xargs -P $NP -n 2 sh -c "$BIN/enum6g_d6 \$0 \$1 94.500001 g GYL > $OUT/out_\$0_\$1.txt 2> $OUT/log_\$0_\$1.txt"
cat $OUT/out_*.txt > $OUT/raw.txt
