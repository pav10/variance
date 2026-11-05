#!/bin/bash
# validation in degree 4: GYL output vs brute force over all x_1 with S <= SM
D=4; SM=32.000001; OUT=/tmp/claude-0/gval
: > $OUT/gyl.txt; : > $OUT/all.txt
for a1 in -2 -1 0 1 2; do for a2 in $(seq -40 5); do
  S=$(python3 -c "print(($D-1)*$a1*$a1/$D - 2*$a2)")
  ok=$(python3 -c "print(1 if 0 < $S <= $SM else 0)")
  [ "$ok" = 1 ] || continue
  /tmp/claude-0/bin/enum6g_d4 $a1 $a2 $SM g GYL >> $OUT/gyl.txt 2>>$OUT/gyl.log
  /tmp/claude-0/bin/enum6g_d4 $a1 $a2 $SM all >> $OUT/all.txt 2>/dev/null
done; done
wc -l $OUT/gyl.txt $OUT/all.txt
