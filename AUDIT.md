# Audit of supplied NEXUS source (internal review)

Source: `nexus(2).zip` supplied by the authors. Inspected statically without running the simulation-dependent code. The files were parsed for Python syntax only.

## Mode controls in `prodigy_source/main.py`

- Lines 50–54: full search grid includes `ll1_values = range(100,8100,100)`, `nn_values=range(5,11)`, `alpha_values=np.arange(0.1,3.1,0.1)`, and `ll2p_values=range(1,101)`.
- Lines 56–71: `best_combination=True`, `best_step=1`, `run_step1=False`, `run_step2=True`. `best_combination=False` leaves full grid in place; `True` calls `get_best_combination()` in `functions.py`, which reads external best-configuration CSV files. With `best_step=1` and Step 2 active, the code explicitly expands `ll2p_values` to 1–100.
- CRITICAL: line 61 uses `sn=sns[1]` (the z~1 snapshot, 96) to select one global combination, while the processing loop subsequently iterates all available science snapshots. This does **not** automatically apply a distinct (nn, alpha, ll1, p2) per redshift as in Table 2 of the manuscript. Needs correction and regression tests.
- `save_candidate_ids=False` (line 84): by default only metric tables are written; scientific distribution should return candidate galaxy IDs/catalogues from both Step 1 and Step 2.

## Portability and data

- `data.py`, `main.py`, `run.py`, `functions.py` and truth-extraction scripts use private absolute paths `/data9/users/jsgomez/...` and `/home/jsgomez/...`.
- `data.py` loads `initial_tracks/*.txt` and AHF metadata; `metric.py` expects full progenitor nodes and other simulation truth inputs. These are unavailable to a generic observational user.
- The ZIP includes only an old short `nexus_best_meanPC_clean_summary.csv` covering Mstar1e9 and Mstar1e10, not the paper's full 12-row calibration table or Step 2 percentiles. Do not publish this as the definitive config.
- The filename and comments still use NEXUS; rename in a controlled refactor rather than blindly replacing all occurrences.

## Scientific validation checks before release

1. `functions.py`, `calcular_densidades_projection()`: `cKDTree.query(points,k=nneighbors)` includes the query point itself; `dists[:,nneighbors-1]` is consequently the (nn-1)-th other galaxy for distinct positions. Reconcile this implementation convention with the manuscript's description of the `n_n`-th nearest neighbour **before changing code**, since altering it might change all calibrations.
2. `metric.py`, `compute_match2()` returns -1 for purity and completeness when either candidate or truth set is empty; the manuscript states that for empty candidates P=C=0. Verify whether downstream aggregation filters sentinels before drawing conclusions.
3. Step 2 is grown on the full `x,y,dens,ID` arrays using an adaptive length estimated from Step 1 members. Check whether the resulting membership is always a true subset of Step 1 and whether the manuscript's wording should specify the exact relationship.
4. Confirm area units and conversion between ckpc/h, cMpc/h and h-free cMpc across algorithms and figures.
5. `run.py` uses `os.system('rm .../*')` to clear log directories. Replace with a safe portable launcher before broader use.
6. Verify that experimental plot modifications or manual calibration values do not appear in any advertised benchmark or automated example.

## Release architecture proposal

- `prodigy/`: reusable detector, no access to descendants or simulation-only truth.
- `calibration/`: metrics, merger-tree truth matching, grid search, and simulations-specific adapters.
- `config/`: 12 author-verified fiducial configurations (z=1,2,3,4; three mass cuts), plus metadata on units and calibration scope.
- `examples/`: small redistributable projected-galaxy catalogue, expected two-step candidate outputs and run instructions.
- `tests/`: exact numerical baselines for representative regions and projections.
- `docs/`: installation, parameter interpretation, input schema, output schema, limitations and citation.

Do not declare an executable public release until the code has passed regression tests against the manuscript results, and the authors have agreed on the license and redistribution permissions.
