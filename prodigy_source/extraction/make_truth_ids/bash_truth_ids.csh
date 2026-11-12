#!/bin/bash

# Generate compact NEXUS truth-ID catalogues for Step 1 and Step 2.
# Edit run_truth_ids.py to change resolution, snapshots, mass cut,
# selected steps, overwrite behaviour, or batch size.

script_dir="$(cd "$(dirname "$0")" && pwd)"
cd "$script_dir" || exit 1

python run_truth_ids.py
