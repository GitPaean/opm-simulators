#!/usr/bin/env bash
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
flow=${1:-/home/kaib/OPM-devel/build/opm-simulators/bin/flow_comp}
output=${2:-/tmp/compositional-cases-0923-runs}

mkdir -p "$output"

for deck in "$here"/*/*.DATA; do
    name=$(basename "$deck" .DATA)
    mkdir -p "$output/$name"
    printf 'Running %s\n' "$name"
    "$flow" "$deck" "--output-dir=$output/$name" \
        > "$output/$name.log" 2>&1
done

python3 "$here/verify_rate_controls.py" "$output"
python3 "$here/verify_outputs.py" "$output" "${@:3}"

printf 'All 25 development cases passed the output and control checks. Logs and output: %s\n' "$output"
