import sys
sys.dont_write_bytecode = True

from data import *
from metric import *
import bz2
import time

# ============================================================
# NEXUS fast version v1.2
#
# Run as:
#     python main_v1.2.py 7k 100
#
# Main changes relative to main_v1.1.py:
#   - Candidate files are no longer written by default.
#   - Metrics are computed immediately after ids_step1/ids_step2.
#   - Truth catalogues are loaded once through metric.py.
#   - Output metric tables are written as:
#       /home/jsgomez/codes/nexus/metric/results/<resolution>/regionXXXX/
#           matches_snXXX_step1.dat
#           matches_snXXX_step2.dat
#
# This keeps the NEXUS candidate search logic almost unchanged, but avoids:
#   candidate IDs -> disk -> match_candidates.py -> metrics
# ============================================================

# -----------------------
# Region and snapshots
# -----------------------
region = int(sys.argv[2])
print(f"\n========== region {region:3d} ==========")

sns = [128, 96, 78, 65, 55]  # sn=128 is needed to define the reference density

#-------------------------------------------------------------------------------------------------
# Stellar-mass cut used for the truth catalogues, matching the metric/match_candidates.py default.
#-------------------------------------------------------------------------------------------------
z0_root_smass_cut   = 0.0
candidate_smass_cut = 1.0e8
truth_smass_cut     = candidate_smass_cut
host_index          = 0

print(f"candidate_smass_cut : {candidate_smass_cut:.3e}")
print(f"truth_smass_cut     : {truth_smass_cut:.3e}")

# -----------------------
# Parameter grid
# -----------------------
ll1_values   = range(100, 8100, 100)      # Mstar1e9, Mstar1e10, Mstar1e11: range(100, 8100, 100)  | Mstar1e8: range(100, 1100, 100)
nn_values    = range(5, 11)               # Mstar1e9, Mstar1e10, Mstar1e11: range(5, 11)           | Mstar1e8: range(9, 11)
alpha_values = np.arange(0.1, 3.1, 0.1)   # Mstar1e9, Mstar1e10, Mstar1e11: np.arange(0.1, 3, 0.1) | Mstar1e8: np.arange(0.1, 3.1, 0.1)

ll2p_values  = range(1, 101)   #[90]   #range(1, 101)

best_combination = True
best_step        = 1   # 1 o 2

if best_combination:
    nn_values, alpha_values, ll1_values, ll2p_values = get_best_combination( # si best_step = 1 te devuelve ll2p_values = [90].
        sn=sns[1],
        candidate_smass_cut=candidate_smass_cut,
        resolution=resolution,
        step=best_step,
    )

run_step1 = False
run_step2 = True

if (run_step2 == True) and (best_combination==True and best_step==1):
    ll2p_values  = range(1, 101)                                      #[90]   #range(1, 101)

# Projection convention:
#   1 -> XY plane, Z depth
#   2 -> XZ plane, Y depth
#   3 -> YZ plane, X depth
projection_values = [1, 2, 3]

# -----------------------
# New v1.2 controls
# -----------------------
# In v1.2 we do NOT write the candidate ID catalogues by default.
# Keep this switch for debugging or for future xy projection figures.
save_candidate_ids = False
save_metrics       = True  # True para corridas normales/calibración


# Metric output directory. This preserves the structure used by your plotting scripts.
def get_mass_label(mass_cut):
    exponent = int(np.log10(mass_cut))
    return f"Mstar1e{exponent}"

metric_out_base   = "/home/jsgomez/codes/nexus/metric/results"
mass_label        = get_mass_label(candidate_smass_cut)
metric_out_region = "%s/%s/%s/region%s" % (
    metric_out_base,
    resolution,
    mass_label,
    str(region).zfill(4)
)
os.makedirs(metric_out_region, exist_ok=True)

# -----------------------
# Density normalization per snapshot
# -----------------------
escala_log_por_snapshot = {
    "128": 11.0,
    "96" : 11.5,
    "78" : 12.0,
    "65" : 12.5,
    "55" : 13.0,
}


def load_clean_galaxies(region, sn, candidate_smass_cut):
    """
    Load AHF_halos once for a region/snapshot and select the galaxy sample used by NEXUS.
    The central progenitor is forced into the sample, as in the original code.
    """

    file_halos = "%s/region%s/%s.z%.3f.AHF_halos" % (
        dir_gizmo,
        region,
        str(sn).zfill(3),
        redshifts[sn]
    )

    if not os.path.exists(file_halos):
        print(f"     sn: {sn}  -> {file_halos} no encontrado")
        return None

    
    fmhires_cut = get_fmhires_cut(resolution)

    df = pd.read_csv(file_halos, sep="\t")

    idp = id_prog_centrals.loc[region-1, str(sn)]
    gal_p = df[df["#ID(1)"] == idp]

    if len(gal_p) == 0:
        print(f"     sn: {sn}  -> central progenitor {idp} no encontrado")
        return None

    gal_f = df[
               (df["M_star(65)"] > candidate_smass_cut) &
               (df["fMhires(38)"] > fmhires_cut)
              ]
    gal = pd.concat([gal_p, gal_f]).drop_duplicates().reset_index(drop=True)

    if resolution == "3k":
        Rvir = float(gal_p["Rvir(12)"].values[0])
    else:
        Rvir = float(gal_p["Rhalo(12)"].values[0])

    return {
        "gal"     : gal,
        "idp"     : int(idp),
        "ID"      : gal["#ID(1)"].values.astype(np.int64),
        "xc"      : float(gal_p["Xc(6)"].values[0]),
        "yc"      : float(gal_p["Yc(7)"].values[0]),
        "zc"      : float(gal_p["Zc(8)"].values[0]),
        "Rvir"    : Rvir,
    }


def compute_density_for_snapshot(gal, sn, nneighbors, projection):
    """
    Wrapper around the projected density estimator.

    projection = 1 -> XY plane, Z depth
    projection = 2 -> XZ plane, Y depth
    projection = 3 -> YZ plane, X depth
    """

    escala_log = escala_log_por_snapshot[str(sn)]
    x, y, dens = calcular_densidades_projection(gal, sn, escala_log, nneighbors, projection)
    dens = np.asarray(dens).reshape(-1)
    return x, y, dens


def get_projected_center(snapshot_data_sn, projection):
    """Return central projected coordinates for the chosen projection."""
    projection = int(projection)

    if projection == 1:
        return snapshot_data_sn["xc"], snapshot_data_sn["yc"]
    elif projection == 2:
        return snapshot_data_sn["xc"], snapshot_data_sn["zc"]
    elif projection == 3:
        return snapshot_data_sn["yc"], snapshot_data_sn["zc"]
    else:
        raise ValueError(f"projection no reconocida: {projection}. Usa 1=xy, 2=xz, 3=yz")


def compute_reference_overdensity(snapshot_data_128, density_data_128, projection):
    """
    Computes the z=0 reference overdensity around Rvir for one projection,
    exactly as in the original code but using the projected centre.
    """

    x, y, dens = density_data_128
    x_center, y_center = get_projected_center(snapshot_data_128, projection)
    Rvir = snapshot_data_128["Rvir"]

    r_centro = np.sqrt((x - x_center)**2 + (y - y_center)**2)
    p_shell = (r_centro > 0.9*Rvir) & (r_centro < 1.1*Rvir)

    if np.sum(p_shell) == 0:
        raise RuntimeError("No hay galaxias en el shell 0.9-1.1 Rvir para definir sobredensidad_ref (END)")

    return float(np.mean(dens[p_shell]))


def compute_candidate_area_from_indices(x, y, idx):
    """
    Candidate convex-hull area in (Mpc/h)^2.

    NEXUS x,y coordinates are in kpc/h, so we divide by 1e3 before computing area.
    """
    idx = np.asarray(idx, dtype=int)

    if len(idx) < 3:
        return 0.0

    return compute_convex_hull_area(x[idx] / 1.0e3, y[idx] / 1.0e3)


def grow_candidates(x, y, dens, ID, x_center, y_center, linking_length, dens_cut):
    """
    FoF-like growth around the central progenitor.

    Returns
    -------
    ids_selected : np.ndarray[int]
    idx_selected : np.ndarray[int]
    """

    n = len(ID)
    visited = np.zeros(n, dtype=bool)

    seed_x = np.array([x_center], dtype=float)
    seed_y = np.array([y_center], dtype=float)

    selected_idx = []

    while len(seed_x) > 0:
        new_idx = []

        for xs, ys in zip(seed_x, seed_y):
            available = np.where(~visited)[0]
            if len(available) == 0:
                continue

            dist = np.sqrt((x[available] - xs)**2 + (y[available] - ys)**2)
            ok = available[(dist < linking_length) & (dens[available] > dens_cut)]

            if len(ok) > 0:
                visited[ok] = True
                new_idx.extend(ok.tolist())

        if len(new_idx) == 0:
            break

        selected_idx.extend(new_idx)
        seed_x = x[new_idx]
        seed_y = y[new_idx]

    selected_idx = np.asarray(selected_idx, dtype=int)
    ids_selected = ID[selected_idx].astype(np.int64)

    return ids_selected, selected_idx

from scipy.spatial import cKDTree
import numpy as np

def grow_candidates_fast(x, y, dens, ID, x_center, y_center, linking_length, dens_cut):
    x = np.asarray(x)
    y = np.asarray(y)
    ID = np.asarray(ID, dtype=np.int64)
    dens = np.asarray(dens).reshape(-1)

    points = np.column_stack([x, y])
    tree = cKDTree(points)

    valid = dens > dens_cut

    seed_idx = tree.query_ball_point([x_center, y_center], r=linking_length)
    seed_idx = [i for i in seed_idx if valid[i]]

    visited = np.zeros(len(ID), dtype=bool)
    queue = list(seed_idx)

    for i in seed_idx:
        visited[i] = True

    head = 0
    while head < len(queue):
        i = queue[head]
        head += 1

        neigh = tree.query_ball_point(points[i], r=linking_length)

        for j in neigh:
            if valid[j] and not visited[j]:
                visited[j] = True
                queue.append(j)

    idx_candidates = np.where(visited)[0]
    ids_candidates = ID[idx_candidates]

    return ids_candidates, idx_candidates

def compute_ll2_from_step1(x, y, idx_step1, ll1, ll2_percentile):
    """
    Compute adaptive ll2 from the minimum distances among Step 1 candidates.
    If Step 1 has <= 1 candidate, return ll1 as a safe fallback.
    """

    if len(idx_step1) <= 1:
        return float(ll1)

    coords = np.column_stack([x[idx_step1], y[idx_step1]])

    dx = coords[:, 0][:, None] - coords[:, 0][None, :]
    dy = coords[:, 1][:, None] - coords[:, 1][None, :]
    dist = np.sqrt(dx*dx + dy*dy)
    np.fill_diagonal(dist, np.inf)

    dist_min = np.min(dist, axis=1)
    dist_min = dist_min[np.isfinite(dist_min)]

    if len(dist_min) == 0:
        return float(ll1)

    return float(np.percentile(dist_min, ll2_percentile))


def write_candidate_file(out_dir, sn, step, nn, alpha, ll1, ll2p, ids, dens_cut, ll2=None, projection=1):
    """
    Write candidate IDs directly into a .txt.bz2 file.

    In v1.2 this is not called unless save_candidate_ids=True.
    """

    os.makedirs(out_dir, exist_ok=True)

    fname = "%s/protocluster.sn%s.z%.3f.proj%d.nn%d.alpha%.2f.ll1%d.ll2p%d.txt.bz2" % (
        out_dir,
        str(sn).zfill(3),
        redshifts[sn],
        projection,
        nn,
        alpha,
        ll1,
        ll2p
    )

    with bz2.open(fname, "wt") as f:
        if step == 1:
            f.write(
                "# N galaxies = %d | ll1 = %.3f (kpc) | nn = %d | alpha = %.2f | dens_cut = %.6e\n" %
                (len(ids), ll1, nn, alpha, dens_cut)
            )
        else:
            f.write(
                "# N galaxies = %d | ll1 = %.3f (kpc) | ll2 = %.3f (kpc) | nn = %d | alpha = %.2f | dens_cut = %.6e\n" %
                (len(ids), ll1, ll2, nn, alpha, dens_cut)
            )

        for gid in ids:
            f.write(f"{int(gid)}\n")

    return fname


# ============================================================
# 1. Load all snapshots once for NEXUS candidate detection
# ============================================================
print("\nLoading snapshots once for NEXUS...")
snapshot_data = {}

for sn in sns:
    print(f"     loading sn={sn}")
    data_sn = load_clean_galaxies(region, sn, candidate_smass_cut)
    if data_sn is not None:
        snapshot_data[sn] = data_sn

if 128 not in snapshot_data:
    raise RuntimeError("sn=128 no esta disponible; no puedo calcular sobredensidad_ref")

valid_sns   = [sn for sn in sns if sn in snapshot_data]
science_sns = [sn for sn in valid_sns if sn != 128]

print(f"Snapshots validos: {valid_sns}")
print(f"Snapshots a procesar: {science_sns}")

# ============================================================
# 2. Load truth catalogues once for metrics
# ============================================================
ids_all_prog, ids_main_cluster, area_all_prog, area_main_cluster = get_progenitors(
    region=region,
    resolution=resolution,
    science_sns=science_sns,
    redshifts=redshifts,
    dir_gizmo=dir_gizmo,
    smass_cut=truth_smass_cut,
    z0_root_smass_cut=z0_root_smass_cut,
    host_index=host_index,
    projections=projection_values,
)

# Metric rows are accumulated by snapshot and step.
metric_rows = {}
for sn in science_sns:
    if run_step1:
        metric_rows[(sn, 1)] = []
    if run_step2:
        metric_rows[(sn, 2)] = []

# ============================================================
# 3. Main parameter sweep
# ============================================================
t0_all = time.time()

n_total_runs = 0

for projection in projection_values:

    print(f"\n================ projection = {projection} ({get_projection_label(projection)}) ================")

    for nn in nn_values:

        print(f"\n========== projection = {projection} | nn = {nn} ==========")

        # --------------------------------------------------------
        # Compute densities once per snapshot for this nn/projection
        # --------------------------------------------------------
        density_data = {}

        for sn in valid_sns:
            print(f"     computing densities: projection={projection}, sn={sn}, nn={nn}")
            gal = snapshot_data[sn]["gal"]
            density_data[sn] = compute_density_for_snapshot(gal, sn, nn, projection)

        # --------------------------------------------------------
        # Reference overdensity at z=0 for this nn/projection
        # --------------------------------------------------------
        sobredensidad_ref = compute_reference_overdensity(
            snapshot_data_128=snapshot_data[128],
            density_data_128=density_data[128],
            projection=projection
        )

        print(f"     sobredensidad_ref = {sobredensidad_ref:.6e}")

        # --------------------------------------------------------
        # Sweep alpha, snapshots, ll1, ll2p
        # --------------------------------------------------------
        for alpha in alpha_values:

            dens_cut = float(alpha * sobredensidad_ref)
            print(f"\n     alpha = {alpha:.2f} | dens_cut = {dens_cut:.6e}")

            for sn in science_sns:

                x, y, dens = density_data[sn]
                ID = snapshot_data[sn]["ID"]
                x_center, y_center = get_projected_center(snapshot_data[sn], projection)

                print(f"          projection={projection} sn={sn}")

                last_ids_step1 = None
                last_metric_step1 = None

                for ll1 in ll1_values:

                    # ----------------------
                    # Step 1
                    # ----------------------
                    ids_step1, idx_step1 = grow_candidates_fast(
                        x=x,
                        y=y,
                        dens=dens,
                        ID=ID,
                        x_center=x_center,
                        y_center=y_center,
                        linking_length=ll1,
                        dens_cut=dens_cut
                    )

                    area_candidates_step1 = compute_candidate_area_from_indices(x, y, idx_step1)

                    for ll2p in ll2p_values:

                        if run_step1:

                            if (
                                last_ids_step1 is not None and
                                len(ids_step1) == len(last_ids_step1) and
                                np.array_equal(ids_step1, last_ids_step1)
                            ):
                                row = last_metric_step1.copy()
                                row[4] = ll1
                                row[7] = ll2p

                            else:
                                row = make_metric_row(
                                    projection=projection,
                                    sn=sn,
                                    z=redshifts[sn],
                                    step=1,
                                    ll1=ll1,
                                    nn=nn,
                                    alpha=alpha,
                                    ll2p=ll2p,
                                    ids_candidates=ids_step1,
                                    ids_truth=ids_all_prog.get(sn, np.array([], dtype=np.int64)),
                                    area_candidates=area_candidates_step1,
                                    area_truth=area_all_prog[projection].get(sn, 0.0),
                                )

                                last_ids_step1 = ids_step1.copy()
                                last_metric_step1 = row.copy()

                            metric_rows[(sn, 1)].append(row)

                            if save_candidate_ids:
                                step1_dir = "%s/%s/%s/step1/region%s" % (
                                    dir_results,
                                    resolution,
                                    mass_label,
                                    str(region).zfill(4)
                                )
                                write_candidate_file(
                                    out_dir=step1_dir,
                                    sn=sn,
                                    step=1,
                                    nn=nn,
                                    alpha=alpha,
                                    ll1=ll1,
                                    ll2p=ll2p,
                                    ids=ids_step1,
                                    dens_cut=dens_cut,
                                    projection=projection
                                )

                        # ----------------------
                        # Step 2
                        # ----------------------
                        if run_step2:
                            ll2 = compute_ll2_from_step1(
                                x=x,
                                y=y,
                                idx_step1=idx_step1,
                                ll1=ll1,
                                ll2_percentile=ll2p
                            )

                            ids_step2, idx_step2 = grow_candidates_fast(
                                x=x,
                                y=y,
                                dens=dens,
                                ID=ID,
                                x_center=x_center,
                                y_center=y_center,
                                linking_length=ll2,
                                dens_cut=dens_cut
                            )

                            area_candidates_step2 = compute_candidate_area_from_indices(x, y, idx_step2)

                            metric_rows[(sn, 2)].append(
                                make_metric_row(
                                    projection=projection,
                                    sn=sn,
                                    z=redshifts[sn],
                                    step=2,
                                    ll1=ll1,
                                    nn=nn,
                                    alpha=alpha,
                                    ll2p=ll2p,
                                    ids_candidates=ids_step2,
                                    ids_truth=ids_main_cluster.get(sn, np.array([], dtype=np.int64)),
                                    area_candidates=area_candidates_step2,
                                    area_truth=area_main_cluster[projection].get(sn, 0.0),
                                )
                            )

                            if save_candidate_ids:
                                step2_dir = "%s/%s/%s/step2/region%s" % (
                                    dir_results,
                                    resolution,
                                    mass_label,
                                    str(region).zfill(4)
                                )
                                write_candidate_file(
                                    out_dir=step2_dir,
                                    sn=sn,
                                    step=2,
                                    nn=nn,
                                    alpha=alpha,
                                    ll1=ll1,
                                    ll2p=ll2p,
                                    ids=ids_step2,
                                    dens_cut=dens_cut,
                                    ll2=ll2,
                                    projection=projection
                                )

                        n_total_runs += 1

# ============================================================
# 4. Save metric tables
# ============================================================
for sn in science_sns:
    if save_metrics:
        if run_step1:
            save_metric_table(metric_rows[(sn, 1)], metric_out_region, sn, 1)
        if run_step2:
            save_metric_table(metric_rows[(sn, 2)], metric_out_region, sn, 2)
    else:
        print("Metric tables not saved because save_metrics = False")

print("\n========== SUMMARY ==========")
print(f"Region              : {region}")
print(f"Resolution          : {resolution}")
print(f"Valid snapshots     : {valid_sns}")
print(f"Science snapshots   : {science_sns}")
print(f"projection values   : {projection_values}")
print(f"ll1 values          : {len(list(ll1_values))}")
print(f"nn values           : {len(list(nn_values))}")
print(f"alpha values        : {len(alpha_values)}")
print(f"ll2p values         : {len(ll2p_values)}")
print(f"Combinations run    : {n_total_runs}")
print(f"Metric out dir      : {metric_out_region}")
print(f"Candidate IDs saved : {save_candidate_ids}")
print(f"Elapsed time [min]  : {(time.time()-t0_all)/60.0:.2f}")
print("END")
