#!/usr/bin/env bash
# Read-level Hi-C workflow on a subsampled .pairs file.
#
# The full library is ~10^9 read pairs; the point here is to show the pair
# classification and deduplication logic, not to rebuild a deep matrix.
set -euo pipefail
THREADS=${THREADS:-8}
PAIRS_IN=$1                       # e.g. data/hic/4DNFIYBWQG4A.pairs.gz
OUT=../results/hic; mkdir -p "$OUT"
CHROMSIZES=../refs/hg38.chrom.sizes
N=5000000

zcat "$PAIRS_IN" | head -n "$N" | bgzip -c > "$OUT/sub.pairs.gz"

# Sort by (chrom1, chrom2, pos1, pos2), then deduplicate. Optical/PCR duplicates
# in Hi-C are identified by the full 4-tuple of positions AND strands -- two
# distinct ligations can share one end, so deduplicating on a single end
# over-collapses real contacts.
pairtools sort --nproc "$THREADS" -o "$OUT/sorted.pairs.gz" "$OUT/sub.pairs.gz"
pairtools dedup --max-mismatch 1 --mark-dups \
  --output "$OUT/dedup.pairs.gz" \
  --output-stats "$OUT/dedup.stats.txt" \
  "$OUT/sorted.pairs.gz"

# Keep only unique-unique mappings on the main assembly.
pairtools select '(pair_type=="UU")' -o "$OUT/valid.pairs.gz" "$OUT/dedup.pairs.gz"

# cis fraction is the go/no-go metric for a Hi-C library: a low cis/trans ratio
# means ligations happened in dilute solution rather than in situ, i.e. the
# crosslinking or the in-nucleus step failed.
pairtools stats -o "$OUT/valid.stats.txt" "$OUT/valid.pairs.gz"
python3 - "$OUT/valid.stats.txt" <<'PY'
import sys
d = dict(l.split('\t')[:2] for l in open(sys.argv[1]) if '\t' in l)
cis, tot = float(d.get('cis', 0)), float(d.get('total_nodups', 1))
print(f"cis fraction: {cis/tot:.3f}   (< ~0.4 suggests random ligation)")
PY

cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
  "$CHROMSIZES":1000 "$OUT/valid.pairs.gz" "$OUT/sub.1000.cool"
cooler zoomify --balance -p "$THREADS" \
  -r 1000,2000,5000,10000,25000,50000,100000 "$OUT/sub.1000.cool"
echo "Wrote $OUT/sub.1000.mcool"
