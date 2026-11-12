import os
import re
import sys
import time
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from shapely.geometry import Polygon, Point
from scipy.spatial import cKDTree

def distancia(u, v, i, j):
    return np.hypot(u[j] - u[i], v[j] - v[i])

def g7(mdist, nneighbors):
    sept = []
    for fila in mdist:
        vec = np.sort(fila[fila > 0])
        #sept.append(vec[5] if len(vec) > 5 else np.nan)
        sept.append(vec[nneighbors-2] if len(vec) > nneighbors-2 else np.nan)
    return np.array(sept)

def calcular_densidades(gal, snap, escala_log, nneighbors):
    x  = gal['Xc(6)'].values
    y  = gal['Yc(7)'].values
    zc = gal['Zc(8)'].values
    D = (np.nanmax(zc)) - np.nanmin(zc)
    N = len(x)
    mdist = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i != j:
                mdist[i, j] = distancia(x, y, i, j)
    sept = g7(mdist, nneighbors)
    dens = np.zeros((N, 1))
    for i in range(N):
        if np.isfinite(sept[i]) and sept[i] > 0:
            #dens[i, 0] = 7 / (np.pi * sept[i]**2) * 10**(escala_log) / D
            dens[i, 0] = nneighbors / (np.pi * sept[i]**2) * 10**(escala_log) / D
    return x, y, dens

from scipy.spatial import cKDTree

def calcular_densidades2(gal, snap, escala_log, nneighbors):
    x  = gal['Xc(6)'].values
    y  = gal['Yc(7)'].values
    zc = gal['Zc(8)'].values

    D = np.nanmax(zc) - np.nanmin(zc)

    points = np.column_stack([x, y])
    tree = cKDTree(points)

    # Incluye la propia galaxia como distancia 0.
    # Por eso la distancia al vecino nneighbors corresponde a columna nneighbors-1.
    dists, _ = tree.query(points, k=nneighbors)

    sept = dists[:, nneighbors - 1]

    dens = np.zeros((len(x), 1))

    p = np.isfinite(sept) & (sept > 0) & np.isfinite(D) & (D > 0)

    dens[p, 0] = (
        nneighbors /
        (np.pi * sept[p]**2) *
        10**escala_log /
        D
    )

    return x, y, dens


# ============================================================
# Projection helpers for NEXUS
# ============================================================
def get_projection_columns(projection):
    """
    Return projected plane columns and line-of-sight column for AHF_halos pandas tables.

    projection = 1 -> XY plane, Z depth
    projection = 2 -> XZ plane, Y depth
    projection = 3 -> YZ plane, X depth
    """
    projection = int(projection)

    if projection == 1:
        return 'Xc(6)', 'Yc(7)', 'Zc(8)'
    elif projection == 2:
        return 'Xc(6)', 'Zc(8)', 'Yc(7)'
    elif projection == 3:
        return 'Yc(7)', 'Zc(8)', 'Xc(6)'
    else:
        raise ValueError(f"projection no reconocida: {projection}. Usa 1=xy, 2=xz, 3=yz")


def get_projection_label(projection):
    projection = int(projection)
    if projection == 1:
        return 'xy'
    elif projection == 2:
        return 'xz'
    elif projection == 3:
        return 'yz'
    else:
        raise ValueError(f"projection no reconocida: {projection}. Usa 1=xy, 2=xz, 3=yz")


def calcular_densidades_projection(gal, snap, escala_log, nneighbors, projection):
    """
    Fast projected density estimator for a chosen line of sight.

    projection = 1 -> neighbours in XY, depth Z
    projection = 2 -> neighbours in XZ, depth Y
    projection = 3 -> neighbours in YZ, depth X
    """
    xcol, ycol, zcol = get_projection_columns(projection)

    x  = gal[xcol].values
    y  = gal[ycol].values
    zc = gal[zcol].values

    D = np.nanmax(zc) - np.nanmin(zc)

    points = np.column_stack([x, y])
    tree = cKDTree(points)

    dists, _ = tree.query(points, k=nneighbors)
    sept = dists[:, nneighbors - 1]

    dens = np.zeros((len(x), 1))

    p = np.isfinite(sept) & (sept > 0) & np.isfinite(D) & (D > 0)

    dens[p, 0] = (
        nneighbors /
        (np.pi * sept[p]**2) *
        10**escala_log /
        D
    )

    return x, y, dens


def change_param(file, name, value, pattern, replacement):
    #print "setting the value of '"+name+"' to "+str(value)+" of the file:"+file
    with open(file, 'rt') as f: original_config = f.read()
    #print "configiguratio of "+file+":", original_config

    # Produce a copy from the original configuration
    # We put it in a dictionary we can modify it easily within set_var
    x = {'config': original_config + ''}

    # Replace variables in configuration so they contain specific values
    pattern     = pattern       % name                                                 #; print "pattern =", pattern
    replacement = replacement   % (name, str(value))
    x['config'] = re.sub(pattern, replacement, x['config'], flags=re.MULTILINE)        #; print x['config']

    # Write final configuration file with updated values
    with open(file, 'wt') as f: f.write(x['config'])

import glob
import time


def check_endfile(logfiles):
    """
    Cuenta cuántos logfiles contienen la palabra END.
    """
    n_ended = 0

    for logfile in logfiles:
        try:
            with open(logfile) as f:
                for line in f:
                    if "END" in line:
                        n_ended += 1
                        break
        except FileNotFoundError:
            pass

    return n_ended


def print_time(t):
    """
    Formatea tiempo en segundos a s, m s, o h m s.
    """
    t = int(t)

    if t < 60:
        return "time: %ds" % t

    elif t < 3600:
        minutes = int(t / 60)
        seconds = t - minutes * 60
        return "time: %dm %ds" % (minutes, seconds)

    else:
        hours   = int(t / 3600)
        minutes = int((t - hours * 3600) / 60)
        seconds = t - hours * 3600 - minutes * 60
        return "time: %dh %dm %ds" % (hours, minutes, seconds)


def progress_bar(dir_log, n_total, code_name, barsize=50, wait_time=1):
    """
    Espera hasta que n_total logfiles tengan END.
    """
    ti = time.time()

    logfiles = glob.glob("%s/*" % dir_log)
    n_ended = check_endfile(logfiles)

    while n_ended < n_total:

        time.sleep(wait_time)

        logfiles = glob.glob("%s/*" % dir_log)
        n_ended = check_endfile(logfiles)

        n_ended_normalized = int((n_ended / n_total) * barsize)
        inicial = "#" * n_ended_normalized
        final   = "-" * (barsize - n_ended_normalized)

        print(
            "[%s%s] | %3d de %3d | %s | codigo %s" %
            (inicial, final, n_ended, n_total, print_time(time.time() - ti), code_name),
            end="\r"
        )

    print("")

def get_best_combination(
    sn,
    candidate_smass_cut,
    resolution="7k",
    step=1,
    metric_best="meanPC",
    zone="allowed",
    projection="all",
    base_dir="/home/jsgomez/science/match_linkings",
):
    """
    Return best NEXUS parameters for Step 1 or Step 2.

    Step 1:
        reads best_step1_configs/best_step1_<res>_<Mstar>_snXXX_allproj.csv
        returns nn_values, alpha_values, ll1_values, ll2p_values

    Step 2:
        reads Step 1 best config plus
        best_step2_configs/best_step2_<res>_<Mstar>_snXXX_allproj.csv
        returns nn_values, alpha_values, ll1_values, ll2p_values

    Selection:
        metric_best == "meanPC"
        zone == "allowed"
        then choose the row with maximum meanPC.
    """

    import os
    import numpy as np
    import pandas as pd

    mass_label = "Mstar1e%d" % int(np.log10(candidate_smass_cut))
    sn_label = "sn%s" % str(int(sn)).zfill(3)

    def _read_best_file(step_number):
        subdir = "best_step%d_configs" % step_number
        fname = "best_step%d_%s_%s_%s_%sproj.csv" % (
            step_number,
            resolution,
            mass_label,
            sn_label,
            projection,
        )

        path = os.path.join(base_dir, subdir, fname)

        if not os.path.exists(path):
            raise FileNotFoundError("No existe best file:\n%s" % path)

        return pd.read_csv(path), path

    def _select_best_row(df, path):
        if "metric_best" not in df.columns:
            raise ValueError("No existe columna metric_best en:\n%s" % path)

        if "zone" not in df.columns:
            raise ValueError("No existe columna zone en:\n%s" % path)

        p = (
            (df["metric_best"].astype(str) == str(metric_best)) &
            (df["zone"].astype(str) == str(zone))
        )

        sub = df[p].copy()

        if len(sub) == 0:
            raise ValueError(
                "No encontre filas para metric_best=%s, zone=%s en:\n%s"
                % (metric_best, zone, path)
            )

        if metric_best not in sub.columns:
            raise ValueError(
                "La columna de metrica %s no existe en:\n%s"
                % (metric_best, path)
            )

        sub = sub[np.isfinite(sub[metric_best].astype(float))].copy()

        if len(sub) == 0:
            raise ValueError(
                "Todas las filas tienen %s no finito en:\n%s"
                % (metric_best, path)
            )

        idx = sub[metric_best].astype(float).idxmax()
        return sub.loc[idx]

    # -------------------------
    # Step 1 best parameters
    # -------------------------
    df1, file1 = _read_best_file(step_number=1)
    row1 = _select_best_row(df1, file1)

    nn = int(row1["nn"])
    alpha = float(row1["alpha"])

    if "ll1" in row1.index:
        ll1 = int(row1["ll1"])
    elif "ll1_mpc" in row1.index:
        ll1 = int(round(float(row1["ll1_mpc"]) * 1000.0))
    else:
        raise ValueError("No encontre ll1 ni ll1_mpc en:\n%s" % file1)

    ll2p = 90

    file2 = None
    row2 = None

    # -------------------------
    # Step 2 best parameter
    # -------------------------
    if int(step) == 2:
        df2, file2 = _read_best_file(step_number=2)
        row2 = _select_best_row(df2, file2)

        if "ll2p" not in row2.index:
            raise ValueError("No encontre columna ll2p en:\n%s" % file2)

        ll2p = int(row2["ll2p"])

    elif int(step) != 1:
        raise ValueError("step debe ser 1 o 2")

    print("")
    print("Using best NEXUS combination")
    print("  step       : %d" % int(step))
    print("  resolution : %s" % resolution)
    print("  mass_label : %s" % mass_label)
    print("  sn         : %d" % int(sn))
    print("  metric     : %s" % metric_best)
    print("  zone       : %s" % zone)
    print("  step1 file : %s" % file1)
    if file2 is not None:
        print("  step2 file : %s" % file2)
    print("  nn         : %d" % nn)
    print("  alpha      : %.3f" % alpha)
    print("  ll1        : %d kpc" % ll1)
    print("  ll2p       : %d" % ll2p)
    print("  Step1 %s   : %.5f" % (metric_best, float(row1[metric_best])))
    if row2 is not None:
        print("  Step2 %s   : %.5f" % (metric_best, float(row2[metric_best])))
    print("")

    return [nn], [alpha], [ll1], [ll2p]
