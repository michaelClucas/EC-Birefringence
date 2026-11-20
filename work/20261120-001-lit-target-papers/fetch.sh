#!/bin/sh
# Re-fetch the sources read in this task into a scratch directory given as $1 (not committed).
# Usage: sh fetch.sh <outdir>   then run h2t.py / abs.py on the results via tools/pin.py.
set -e
out="$1"; mkdir -p "$out"
for id in 2509.13654 2510.07928 2601.13624 2512.21888 2202.13919 2011.11254 2205.13962; do
  curl -sL -A "research-client" -o "$out/$id.html" "https://arxiv.org/html/$id"
done
ids="2011.11254,2201.07682,2205.13962,2201.13347,2305.02268,2410.12089,2001.10465,2006.08061,1705.02523,2504.13154,2401.11973,2206.01635,2008.10334,2404.13771,2502.07654,2503.14452,2403.02096,2203.16567,2408.08063,1509.02461,1705.06387,1206.5546,2006.15982"
curl -sL -A "research-client" -o "$out/abs.xml" "https://export.arxiv.org/api/query?id_list=$ids&max_results=30"
