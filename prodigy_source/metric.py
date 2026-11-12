#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

import bz2
import shutil


# ============================================================
# metric.py
#
# Helper functions to compute purity/completeness/F1 directly
# inside NEXUS, without first writing candidate IDs to disk.
#
# This module is intended to live next to main_v1.2.py and be
# imported with:
#
#     from metric import *
#
# ============================================================

def get_fmhires_cut(resolution):
    if resolution == "3k":
        return 0.996
    elif resolution == "7k":
        return 0.9995
    else:
        raise ValueError(f"resolution no reconocida: {resolution}")

def ensure_1d_int(arr):
    """Return arr as a 1D int64 array. Empty input -> empty int64 array."""
    arr = np.asarray(arr)
    if arr.size == 0:
        return np.array([], dtype=np.int64)
    return np.atleast_1d(arr).astype(np.int64)


# ============================================================
# Basic match statistics
# ============================================================
def compute_match(ids_candidates, ids_truth):
    """
    Compute N_candidates, N_truth, N_match, purity, completeness and F1.

    purity      = N_match / N_candidates
    completeness = N_match / N_truth
    F1          = 2 * purity * completeness / (purity + completeness)
    """

    ids_candidates = np.unique(ensure_1d_int(ids_candidates))
    ids_truth      = np.unique(ensure_1d_int(ids_truth))

    ids_match = np.intersect1d(ids_candidates, ids_truth, assume_unique=True)

    n_candidates = len(ids_candidates)
    n_truth      = len(ids_truth)
    n_match      = len(ids_match)

    purity = n_match / n_candidates if n_candidates > 0 else -1.0
    completeness = n_match / n_truth if n_truth > 0 else -1.0

    if purity + completeness > 0:
        f1 = 2.0 * purity * completeness / (purity + completeness)
    else:
        f1 = -1.0

    return n_candidates, n_truth, n_match, purity, completeness, f1

def compute_match2(ids_candidates, ids_truth):
    ids_candidates = ensure_1d_int(ids_candidates)
    ids_truth      = ensure_1d_int(ids_truth)

    n_candidates = len(ids_candidates)
    n_truth      = len(ids_truth)

    if n_candidates == 0 or n_truth == 0:
        return n_candidates, n_truth, 0, -1.0, -1.0, -1.0

    # ids_truth ya debería venir único desde get_all_progenitors_ids/get_main_cluster_ids.
    # ids_candidates también debería venir único desde grow_candidates_fast.
    p_match = np.isin(ids_candidates, ids_truth, assume_unique=True)
    n_match = int(np.sum(p_match))

    purity = n_match / n_candidates if n_candidates > 0 else -1.0
    completeness = n_match / n_truth if n_truth > 0 else -1.0

    if purity + completeness > 0:
        f1 = 2.0 * purity * completeness / (purity + completeness)
    else:
        f1 = -1.0

    return n_candidates, n_truth, n_match, purity, completeness, f1

def compute_convex_hull_area(x, y):
    """
    Convex-hull area in the units of x*y.

    If there are fewer than three unique projected points, the area is zero.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    p = np.isfinite(x) & np.isfinite(y)
    x = x[p]
    y = y[p]

    if len(x) < 3:
        return 0.0

    points = np.column_stack([x, y])
    points = np.unique(points, axis=0)

    if len(points) < 3:
        return 0.0

    try:
        hull = ConvexHull(points)
        return float(hull.volume)  # in 2D, volume is area
    except Exception:
        return 0.0


def get_projected_xy_from_props(props, projection):
    """
    Return projected x,y arrays from a properties dictionary.

    Properties are already in Mpc/h in metric.py.

    projection = 1 -> XY plane
    projection = 2 -> XZ plane
    projection = 3 -> YZ plane
    """
    projection = int(projection)

    if projection == 1:
        return props["x"], props["y"]
    elif projection == 2:
        return props["x"], props["z"]
    elif projection == 3:
        return props["y"], props["z"]
    else:
        raise ValueError(f"projection no reconocida: {projection}. Usa 1=xy, 2=xz, 3=yz")


def compute_truth_area(props_sn, ids_truth, projection):
    """Convex-hull area for a truth catalogue at one snapshot/projection."""
    ids_truth = ensure_1d_int(ids_truth)

    if len(ids_truth) == 0:
        return 0.0

    xproj, yproj = get_projected_xy_from_props(props_sn, projection)
    mask = np.isin(props_sn["id_halo"], ids_truth, assume_unique=False)

    return compute_convex_hull_area(xproj[mask], yproj[mask])


def make_metric_row(
    projection,
    sn,
    z,
    step,
    ll1,
    nn,
    alpha,
    ll2p,
    ids_candidates,
    ids_truth,
    area_candidates=-1.0,
    area_truth=-1.0,
):
    """Build one output row for matches_snXXX_stepY.dat.bz2."""

    n_cand, n_truth, n_match, purity, completeness, f1 = compute_match2(
        ids_candidates,
        ids_truth
    )

    area_candidates = float(area_candidates)
    area_truth      = float(area_truth)

    if area_truth > 0:
        area_ratio = area_candidates / area_truth
    else:
        area_ratio = -1.0

    return [
        int(projection), sn, z, step, ll1, nn, alpha, ll2p,
        n_cand, n_truth, n_match, purity, completeness, f1,
        area_candidates, area_truth, area_ratio
    ]


# ============================================================
# Output table writer
# ============================================================
def save_metric_table(rows, out_region, sn, step):
    """
    Save all metric rows for one snapshot and one step.

    Output columns:
        projection sn z step ll1 nn alpha ll2p N_candidates N_truth N_match
        purity completeness F1 area_candidates area_truth area_ratio
    """

    if len(rows) == 0:
        return None

    os.makedirs(out_region, exist_ok=True)

    df = pd.DataFrame(rows, columns=[
        "projection", "sn", "z", "step", "ll1", "nn", "alpha", "ll2p",
        "N_candidates", "N_truth", "N_match",
        "purity", "completeness", "F1",
        "area_candidates", "area_truth", "area_ratio"
    ])

    out_summary = f"{out_region}/matches_sn{str(sn).zfill(3)}_step{step}.dat"

    header = (
        f"{'projection':>10s} "
        f"{'sn':>4s} "
        f"{'z':>10s} "
        f"{'step':>5s} "
        f"{'ll1':>7s} "
        f"{'nn':>4s} "
        f"{'alpha':>8s} "
        f"{'ll2p':>6s} "
        f"{'N_candidates':>14s} "
        f"{'N_truth':>10s} "
        f"{'N_match':>10s} "
        f"{'purity':>12s} "
        f"{'completeness':>14s} "
        f"{'F1':>12s} "
        f"{'area_candidates':>18s} "
        f"{'area_truth':>14s} "
        f"{'area_ratio':>12s}"
    )

    np.savetxt(
        out_summary,
        df.values,
        fmt=[
            "%10d",     # projection
            "%4d",      # sn
            "%10.6f",   # z
            "%5d",      # step
            "%7d",      # ll1
            "%4d",      # nn
            "%8.2f",    # alpha
            "%6d",      # ll2p
            "%14d",     # N_candidates
            "%10d",     # N_truth
            "%10d",     # N_match
            "%12.8f",   # purity
            "%14.8f",   # completeness
            "%12.8f",   # F1
            "%18.8f",   # area_candidates
            "%14.8f",   # area_truth
            "%12.8f",   # area_ratio
        ],
        header=header,
        comments=""
    )

    with open(out_summary, "rb") as f_in:
        with bz2.open(out_summary + ".bz2", "wb", compresslevel=9) as f_out:
            shutil.copyfileobj(f_in, f_out)

    os.remove(out_summary)

    print(f"      saved metric table: {out_summary}.bz2")

    return out_summary + ".bz2"


# ============================================================
# GIZMO properties and truth catalogues
# ============================================================
def load_gizmo_properties(dir_gizmo, region, sn, redshift):
    """
    Load only the AHF_halos columns needed for truth catalogues.

    Parameters
    ----------
    dir_gizmo : str
        Resolution-specific GIZMO directory, e.g.
        /data9/.../music_planck/7k
    region : int
    sn : int
    redshift : float
    """

    props = {}

    file_halos = f"{dir_gizmo}/region{region}/{str(sn).zfill(3)}.z{redshift:.3f}.AHF_halos"

    (props["id_halo"], props["smass"], props["hmass"], props["id_hosthalo"],
     props["x"], props["y"], props["z"], props["r200"], props["fmhires"]) = np.loadtxt(
        file_halos,
        usecols=(0, 64, 3, 1, 5, 6, 7, 11, 37),
        unpack=True
    )

    props["id_halo"]     = props["id_halo"].astype(np.int64)
    props["id_hosthalo"] = props["id_hosthalo"].astype(np.int64)

    # kpc/h -> Mpc/h, consistent with the old metric code.
    for key in ["x", "y", "z", "r200"]:
        props[key] = props[key] / 1.0e3

    return props


def recenter_properties(props, x0, y0, z0):
    """Recenter positions relative to the z=0 main halo position."""
    props = props.copy()
    props["x"] = props["x"] - x0
    props["y"] = props["y"] - y0
    props["z"] = props["z"] - z0
    return props


def get_z0_cluster_members(props128, smass_cut_for_roots=0.0, host_index=0):
    """
    Return z=0 members inside R200 of the main halo.

    These IDs are used as root IDs for the full progenitor nodes file.
    The main halo itself is included.
    """

    host_id    = int(props128["id_halo"][host_index])
    xh, yh, zh = props128["x"][host_index], props128["y"][host_index], props128["z"][host_index]
    r200       = props128["r200"][host_index]

    p_sub = np.where(props128["id_hosthalo"] == host_id)[0]

    dist = np.sqrt(
        (props128["x"][p_sub] - xh)**2 +
        (props128["y"][p_sub] - yh)**2 +
        (props128["z"][p_sub] - zh)**2
    )
    p_sub = p_sub[dist < r200]

    p_members = np.unique(np.concatenate([[host_index], p_sub]))

    if smass_cut_for_roots > 0:
        p_members = p_members[props128["smass"][p_members] > smass_cut_for_roots]

    return props128["id_halo"][p_members].astype(np.int64), p_members


def load_full_nodes(nodes_file):
    """Load only the required columns from full progenitor nodes."""

    return np.loadtxt(
        nodes_file,
        usecols=(0, 1, 2, 4),
        dtype=[("root", "i8"), ("sn", "i8"), ("id", "i8"), ("smass", "f8")]
    )


def get_all_progenitors_ids(nodes, root_ids_z0, sn, smass_cut,
                            props_sn=None, fmhires_cut=None):
    """Load only the required columns from full progenitor nodes."""

    mask = (
        np.isin(nodes["root"], root_ids_z0) &
        (nodes["sn"] == sn) &
        (nodes["smass"] > smass_cut)
    )

    ids = np.unique(nodes["id"][mask].astype(np.int64))

    if props_sn is not None and fmhires_cut is not None:
        p_clean = props_sn["fmhires"] > fmhires_cut
        ids_clean = props_sn["id_halo"][p_clean].astype(np.int64)
        ids = np.intersect1d(ids, ids_clean)

    return ids


def load_old_main_progenitor_tracks(dir_gizmo, region):
    """
    Load old main-progenitor tracks.

    These are used only to define the Step 2 truth catalogue, i.e. the
    main-cluster/core around the main progenitor of the z=0 main halo.
    """

    base = f"{dir_gizmo}/region{region}"
    tracks = {}

    for prop in ["idhalo", "x", "y", "z", "r200", "smass"]:
        fname = f"{base}/128.z0.000.AHF_progenitors_{prop}"
        tracks[prop] = np.loadtxt(fname)
        if prop in ["x", "y", "z", "r200"]:
            tracks[prop] = tracks[prop] / 1.0e3

    tracks["idhalo"] = tracks["idhalo"].astype(np.int64)
    return tracks


def recenter_tracks(tracks, x0, y0, z0):
    """Recenter track positions relative to z=0 main halo position."""
    tracks = tracks.copy()
    tracks["x"] = tracks["x"] - x0
    tracks["y"] = tracks["y"] - y0
    tracks["z"] = tracks["z"] - z0
    return tracks


def get_main_cluster_ids(props_sn, tracks, host_index, sn, smass_cut,
                         fmhires_cut=None):
    """
    Return the core/main-cluster truth catalogue at snapshot sn.

    Definition: central main progenitor of the z=0 main halo + its subhalos
    inside its R200, applying smass_cut.
    """

    main_id = int(tracks["idhalo"][host_index, sn])

    if main_id <= 0:
        return np.array([], dtype=np.int64)

    xc = tracks["x"][host_index, sn]
    yc = tracks["y"][host_index, sn]
    zc = tracks["z"][host_index, sn]
    rc = tracks["r200"][host_index, sn]

    p_sub = np.where(props_sn["id_hosthalo"] == main_id)[0]
    p_cen = np.where(props_sn["id_halo"] == main_id)[0]
    p = np.unique(np.concatenate([p_cen, p_sub]))

    if len(p) == 0:
        return np.array([], dtype=np.int64)

    dist = np.sqrt(
        (props_sn["x"][p] - xc)**2 +
        (props_sn["y"][p] - yc)**2 +
        (props_sn["z"][p] - zc)**2
    )

    mask = (
        (dist < rc) &
        (props_sn["smass"][p] > smass_cut)
    )

    if fmhires_cut is not None:
        mask = mask & (props_sn["fmhires"][p] > fmhires_cut)

    p = p[mask]

    return np.unique(props_sn["id_halo"][p].astype(np.int64))

def get_progenitors(
    region,
    resolution,
    science_sns,
    redshifts,
    dir_gizmo,
    smass_cut=1.0e10,
    z0_root_smass_cut=0.0,
    host_index=0,
    projections=(1, 2, 3),
):
    """
    Load truth catalogues once and return dictionaries by snapshot.

    Returns
    -------
    ids_all_prog : dict[int, np.ndarray]
        ids_all_prog[sn] are all real progenitors of z=0 cluster members.

    ids_main_cluster : dict[int, np.ndarray]
        ids_main_cluster[sn] are the Step 2 core truth catalogues.
    """

    print("\nLoading truth catalogues for metrics...")

    # z=0 cluster members used as roots for all-progenitor truth.
    props128_raw = load_gizmo_properties(dir_gizmo, region, 128, redshifts[128])
    x0, y0, z0 = props128_raw["x"][host_index], props128_raw["y"][host_index], props128_raw["z"][host_index]
    props128 = recenter_properties(props128_raw, x0, y0, z0)

    root_ids_z0, _ = get_z0_cluster_members(
        props128,
        smass_cut_for_roots=z0_root_smass_cut,
        host_index=host_index
    )
    print(f"  z=0 cluster members used as roots: {len(root_ids_z0)}")

    # Full progenitor tree nodes.
    nodes_file = f"{dir_gizmo}/region{region}/128.AHF_full_progenitors_nodes.bz2"
    print(f"  loading full progenitors nodes: {nodes_file}")
    nodes = load_full_nodes(nodes_file)

    # Old tracks only for Step 2 core truth.
    print("  loading old main-progenitor tracks for Step 2 core")
    tracks = load_old_main_progenitor_tracks(dir_gizmo, region)
    tracks = recenter_tracks(tracks, x0, y0, z0)

    ids_all_prog = {}
    ids_main_cluster = {}
    area_all_prog = {int(proj): {} for proj in projections}
    area_main_cluster = {int(proj): {} for proj in projections}

    for sn in science_sns:
        redshift = redshifts[sn]
        fmhires_cut = get_fmhires_cut(resolution)

        props_sn_raw = load_gizmo_properties(dir_gizmo, region, sn, redshift)
        props_sn = recenter_properties(props_sn_raw, x0, y0, z0)

        ids_all_prog[sn] = get_all_progenitors_ids(
            nodes,
            root_ids_z0,
            sn,
            smass_cut,
            props_sn=props_sn,
            fmhires_cut=fmhires_cut
        )

        ids_main_cluster[sn] = get_main_cluster_ids(
            props_sn,
            tracks,
            host_index,
            sn,
            smass_cut,
            fmhires_cut=fmhires_cut
        )

        for proj in projections:
            proj = int(proj)
            area_all_prog[proj][sn] = compute_truth_area(
                props_sn,
                ids_all_prog[sn],
                projection=proj
            )
            area_main_cluster[proj][sn] = compute_truth_area(
                props_sn,
                ids_main_cluster[sn],
                projection=proj
            )

        print(
            f"  sn={sn:3d}: truth all={len(ids_all_prog[sn])}, "
            f"truth core={len(ids_main_cluster[sn])}"
        )

    return ids_all_prog, ids_main_cluster, area_all_prog, area_main_cluster
