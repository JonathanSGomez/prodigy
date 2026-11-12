# Coauthor review guide — PRODIGY (private staging)

This is the historical NEXUS source snapshot supplied for the PRODIGY manuscript. Its purpose is **scientific and software review**, not a public installation or production use.

## Suggested review order

1. **Scientific definitions:** compare `legacy_nexus/functions.py` and `legacy_nexus/metric.py` with the manuscript's definitions of nearest neighbours, Step 1/2 candidate membership, ground truth, meanPC, and the candidate-to-truth area criterion.
2. **Calibration modes:** inspect the flags and control flow in `legacy_nexus/main.py`, especially `best_combination`, `best_step`, `run_step1`, `run_step2`, `save_candidate_ids`, and snapshot/mass dependence.
3. **Reproducibility:** identify the minimal permitted input catalogue and provenance files needed to reproduce at least one $z$, mass cut, region and projection; confirm units and corresponding published metrics.
4. **Publication readiness:** decide which algorithms are intended as public functionality, which science-derived inputs can be redistributed, and who approves the license/copyright attribution.

## Specific scientific questions to resolve

- **Nearest-neighbour convention:** `cKDTree.query(points, k=nneighbors)` includes the source galaxy. Does the scientific definition consistently refer to the `n_n`-th *other* galaxy? Changing this convention would require rerunning/revalidating the calibration.
- **Correct fiducial per snapshot:** the historical code loads one `sn=sns[1]` combination in advance. Match the final Step 1 and Step 2 settings to redshift and stellar-mass cut independently.
- **Core membership:** verify whether Step 2 always returns a subset of the Step 1 detection, including the recursive-growth implementation.
- **Metric empty-set behaviour:** reconcile the sentinel `-1` path in `metric.py` with the manuscript's stated zero-purity/completeness convention.
- **Area constraint and units:** cross-check comoving and $h$ units and the exact ensemble averaging criterion.
- **Dispersions:** distinguish near-optimal configuration ranges ($\Delta\mathrm{meanPC}\le 0.02$) from statistical scatter across regions/projections. Verify the final values against the data used for the paper.

## Collaboration workflow

- Create a **GitHub Issue** for each distinct concern and quote the filename, function, and line range.
- Keep numerical-algorithm changes separate from documentation or portability changes.
- Prefer small pull requests with a regression comparison to the current paper outputs.
- Do **not** upload large simulation products, private local paths in new configuration files, SSH keys, or unpublished result archives.
- Confirm that every planned public distribution is approved by the relevant coauthors and data owners.

## What is *not* being claimed

This repository is not currently a pip-installable package; does not reproduce the final manuscript figures without additional private inputs; has not undergone independent scientific verification; and is not yet ready for Euclid/other survey catalogues.
