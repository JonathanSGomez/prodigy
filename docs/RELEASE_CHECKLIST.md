# Roadmap toward a public PRODIGY release

All items below are **pending** unless explicitly confirmed in a later commit.

## Scientific validation

- [ ] Freeze the calibration dataset and the exact published `z`/stellar-mass configurations.
- [ ] Resolve the nearest-neighbour self-inclusion convention without silently changing calibrated results.
- [ ] Reproduce one Step 1 and Step 2 result against the paper, including purity, completeness, meanPC, and areas.
- [ ] Verify the Step 2 subset/core interpretation and metric edge cases.
- [ ] Verify comoving units and the calibration aggregation/area constraint.

## Engineering

- [ ] Split pure detection from merger-tree ground truth and calibration, and introduce a documented CLI/API.
- [ ] Replace absolute workstation and cluster paths with configuration parameters.
- [ ] Replace destructive log cleanup in the historical runner.
- [ ] Provide redshift-/mass-dependent parameter tables with data provenance.
- [ ] Publish a legally redistributable small example with expected candidate IDs for both steps.
- [ ] Add deterministic tests, CI, version constraints, and cross-platform installation instructions.

## Governance and release

- [ ] Obtain collaborator approval and check third-party code/data redistribution rights.
- [ ] Confirm copyright attribution and BSD-3-Clause license choice with the full author team.
- [ ] Write `CITATION.cff` based on the accepted manuscript and authors' approved metadata.
- [ ] Create an initial release tag, archive it with Zenodo, and use the DOI in the A&A Code Availability statement.
- [ ] Change the GitHub repository visibility to public **only** after release approval.
