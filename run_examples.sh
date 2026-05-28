#!/usr/bin/env bash
# Build the example analyses and open them (Mac/Linux). Needs Python 3.
cd "$(dirname "$0")" || exit 1
for c in pairwise_advanced nma_inconsistency dta cnma; do
  python3 run.py "configs/example_${c}.json" || { echo "Is Python 3 installed?"; exit 1; }
done
opener=""
command -v xdg-open >/dev/null 2>&1 && opener=xdg-open
command -v open >/dev/null 2>&1 && opener=open
if [ -n "$opener" ]; then
  for c in pairwise_advanced nma_inconsistency dta cnma; do "$opener" "output/example_${c}.html"; done
fi
echo "AACT example (needs a local AACT snapshot; export TRIALFORGE_AACT=/path first):"
echo "  python3 run.py configs/example_aact_finerenone.json"
echo "Make your own: copy a file in configs/, edit it, then: python3 run.py configs/your_file.json"
