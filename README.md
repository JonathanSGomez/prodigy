# PRODIGY

**PROtocluster Discovery and Identification algorithm calibrated on simulated GalaxY clusters**  
*Tracing galaxy protocluster regions with The Three Hundred Project*

> **Internal scientific review — source snapshot, not a public software release.**  
> This repository is currently private. The code below has not been validated as a standalone package or on observational catalogues. Do not cite this snapshot as a released version or redistribute its files without the authors' agreement.

PRODIGY is a two-stage framework for identifying galaxy protocluster candidates in projected galaxy distributions. **Step 1** selects the extended progenitor environment using nearest-neighbour density, an overdensity cut, and FoF-like connectivity. **Step 2** identifies a smaller, higher-confidence core through an adaptive linking length derived from galaxy separations inside the Step 1 selection. In the calibration, each step is evaluated against a different descendant-defined galaxy population from **The Three Hundred Project**.

The associated manuscript is *PRODIGY: Identifying the Progenitors of Present-day Galaxy Clusters* (Gómez, Ceverino, Hough & Yepes; in preparation). The manuscript considers $z=1,2,3,4$ and stellar-mass selections $M_\star>10^8,10^9,10^{10}\,M_\odot$.

## What is available here?

| Location | Purpose | Review status |
| --- | --- | --- |
| [`legacy_nexus/main.py`](legacy_nexus/main.py) | Historical Step 1/2 grid-search and metrics runner | Original source; needs mode separation |
| [`legacy_nexus/functions.py`](legacy_nexus/functions.py) | Projected density estimates, connectivity, parameter lookup | Original source; neighbour convention to validate |
| [`legacy_nexus/metric.py`](legacy_nexus/metric.py) | Simulation truth matching, purity/completeness, and areas | Original source; empty-set cases to validate |
| [`legacy_nexus/data.py`](legacy_nexus/data.py) | AHF catalogue and progenitor-data paths | Requires private paths and data |
| [`legacy_nexus/run.py`](legacy_nexus/run.py) | Historical batch launcher | Not portable; may delete logs |
| [`legacy_nexus/extraction/`](legacy_nexus/extraction/) | Supporting truth and truth-area extraction | Simulation-specific; review before reuse |
| [`legacy_nexus/initial_tracks/generate_centralhalosIDS.py`](legacy_nexus/initial_tracks/generate_centralhalosIDS.py) | Main-progenitor track preparation | Simulation-specific |
| [`AUDIT.md`](AUDIT.md) | Technical issues identified during static inspection | Prioritised review required |
| [`docs/REVIEW_GUIDE.md`](docs/REVIEW_GUIDE.md) | Coauthor review questions and process | Ready for comments |
| [`docs/RELEASE_CHECKLIST.md`](docs/RELEASE_CHECKLIST.md) | Requirements before public release | Not yet complete |
| [`docs/SOURCE_MANIFEST.md`](docs/SOURCE_MANIFEST.md) | Provenance and checksums of original scripts | For traceability |

### Deliberately **not** included

- Large `metric/results/` parameter-sweep tables and logs, including the omitted massive-tracer outputs.
- The Three Hundred simulation catalogues and merger-tree/halo products.
- Region truth tables, `initial_tracks/*.txt`, derived candidate IDs, and private local configurations.
- The manuscript's full, author-verified 12-row fiducial parameter file: **pending verification against the final calibrations**.

These are not necessary for browsing the source; redistributing simulation products requires separate authorisation.

## How `main.py` is currently configured

The historical runner uses `best_combination = False` to leave the full Step 1 and Step 2 parameter grids available, whereas `best_combination = True` reads preselected parameters from an **external** CSV via `get_best_combination()`.

**Important:** the supplied source uses `sn=sns[1]` to load that combination before looping over snapshots. It therefore does **not** yet automatically select the correct parameters for every redshift and stellar-mass threshold in the paper. The current defaults also enable Step 2 metrics but **do not** write candidate catalogues. Do not interpret this snapshot as a validated implementation of the final fiducial table.

## Setup and reproducibility status

The original scripts import Python, NumPy, Pandas, SciPy, Matplotlib and Shapely; see [`requirements.txt`](requirements.txt). Dependency versions, Python compatibility, scientific outputs, and execution on a clean computer remain **untested**. The scripts use internal absolute paths and simulation-derived ground truth. They **cannot currently be run from this repository alone**. Please review the source without launching `legacy_nexus/run.py` on a working machine; it contains log-directory cleanup code.

The intended public release will offer separately documented **detection**, **Step 1 calibration**, and **Step 2 calibration** modes, together with independently reproducible tests, a redistributable demonstration catalogue, and properly versioned configurations.

## For coauthors

Please begin with **[the coauthor review guide](docs/REVIEW_GUIDE.md)** and **[the scientific audit](AUDIT.md)**. GitHub Issues can be used for feedback with file and line references. The historical files are preserved unmodified in `legacy_nexus/`; changes to the calibrated algorithm should be discussed and regression-tested before being merged into a new package.

## Licensing and citation

A BSD 3-Clause [`LICENSE`](LICENSE) file was created when the private repository was initialised. **The team must confirm copyright attribution, coauthor approval and the redistribution status of included components before public release.** A paper DOI, software DOI/Zenodo archive, and `CITATION.cff` will be added only once the scientific release is approved and validated.
