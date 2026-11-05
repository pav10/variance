#!/bin/bash
# usage: src/runcensus.sh DEG NMAX [NSH] [MAXLINES]
# runs the witness census on data/fields/degDEG.txt in NSH parallel gp processes
DEG=$1; NMAX=$2; NSH=${3:-4}; MAXL=${4:-0}
cd "$(dirname "$0")/.."
IN=data/fields/deg$DEG.txt
if [ "$MAXL" != "0" ]; then head -n $((MAXL+1)) $IN > /tmp/claude-0/in_$DEG.txt; IN=/tmp/claude-0/in_$DEG.txt; fi
mkdir -p results/census
for s in $(seq 0 $((NSH-1))); do
  OUT=results/census/deg${DEG}_N${NMAX}_part$s.txt; rm -f $OUT
  echo "TMAX=${TMAX:-0}; census(\"$IN\", \"$OUT\", $s, $NSH, $NMAX)" | gp -q -D parisizemax=4000000000 src/census.gp > results/census/deg${DEG}_part$s.log 2>&1 &
done
wait
cat results/census/deg${DEG}_N${NMAX}_part*.txt | sort -t'|' -k1,1n > results/census/deg${DEG}_N${NMAX}.txt
rm -f results/census/deg${DEG}_N${NMAX}_part*.txt
echo "done: $(wc -l < results/census/deg${DEG}_N${NMAX}.txt) fields"
cut -d'|' -f7 results/census/deg${DEG}_N${NMAX}.txt | sort | uniq -c
