#!/usr/bin/env python3
"""
make_truth_areas.py

Compute projected convex-hull areas for previously saved NEXUS truth-ID
catalogues, without rebuilding the progenitor trees every time an area-vs-z
figure is made.

The calculation reproduces the truth-area definition used by NEXUS:

    1. Read the truth halo IDs at each available snapshot.
    2. Match those IDs to the corresponding AHF_halos catalogue.
    3. Convert AHF positions from ckpc/h to cMpc/h by dividing by 1000.
    4. Compute the 2D convex-hull area for the three projections:
           projection 1: XY
           projection 2: XZ
           projection 3: YZ
    5. Return area = 0 when fewer than three unique projected points exist.

A single compressed CSV is written per region and Step. It contains the
three projected areas plus their median and mean. Both native (cMpc/h)^2
and h-free cMpc^2 values are saved.

Example
-------
python make_truth_areas.py 7k \
    --regions 70 \
    --smass-cut 1.0e10 \
    --step 1 \
    --overwrite
"""

from __future__ import annotations

import argparse
import bz2
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, QhullError


TRUTH_FILE_RE = re.compile(
    r"^truth\.sn(?P<sn>\d+)\."
    r"z(?P<z>[-+0-9.eE]+)\."
    r"step(?P<step>\d+)\.txt\.bz2$"
)


def mass_label_from_cut(smass_cut: float) -> str:
    """Convert 1e10 -> Mstar1e10."""
    if smass_cut <= 0:
        raise ValueError("smass_cut debe ser positivo.")

    exponent = int(round(np.log10(float(smass_cut))))

    if not np.isclose(float(smass_cut), 10.0 ** exponent):
        raise ValueError(
            "Este flujo espera cortes de masa de la forma 10**N. "
            f"Valor recibido: {smass_cut:.6e}"
        )

    return f"Mstar1e{exponent}"


def default_regions(resolution: str) -> List[int]:
    """Default The Three Hundred region list."""
    resolution = str(resolution)

    if resolution == "7k":
        return [region for region in range(1, 325) if region != 107]

    if resolution == "3k":
        return list(range(1, 325))

    raise ValueError(
        f"resolution no reconocida: {resolution}. Usa 3k o 7k."
    )


def load_ids_bz2(path: Path) -> np.ndarray:
    """Read one integer halo ID per line from a bz2 file."""
    ids: List[int] = []

    with bz2.open(path, "rt") as stream:
        for line in stream:
            value = line.strip()

            if value:
                ids.append(int(value))

    if not ids:
        return np.array([], dtype=np.int64)

    return np.unique(np.asarray(ids, dtype=np.int64))


def discover_truth_catalogues(
    region_dir: Path,
    step: int,
) -> List[Dict[str, object]]:
    """
    Discover all truth-ID catalogues available in one region directory.
    """
    if not region_dir.is_dir():
        raise FileNotFoundError(region_dir)

    catalogues: List[Dict[str, object]] = []
    seen_snapshots: Dict[int, Path] = {}

    for path in sorted(region_dir.glob(f"truth.sn*.z*.step{step}.txt.bz2")):
        match = TRUTH_FILE_RE.match(path.name)

        if match is None:
            continue

        file_step = int(match.group("step"))

        if file_step != int(step):
            continue

        sn = int(match.group("sn"))
        redshift = float(match.group("z"))

        if sn in seen_snapshots:
            raise RuntimeError(
                f"Hay más de un fichero truth para sn={sn}:\n"
                f"  {seen_snapshots[sn]}\n"
                f"  {path}"
            )

        seen_snapshots[sn] = path

        catalogues.append(
            {
                "sn": sn,
                "z": redshift,
                "truth_file": path,
            }
        )

    # Lower snapshots correspond to earlier cosmic epochs. Sorting by sn
    # therefore naturally orders the catalogue from high z toward z=0.
    catalogues.sort(key=lambda item: int(item["sn"]))

    return catalogues


def resolve_halo_file(
    simulation_region_dir: Path,
    sn: int,
    redshift: float,
) -> Path:
    """
    Find the matching AHF_halos file.

    First use the exact z rounded to three decimals, matching the truth filename.
    If that exact path is unavailable, allow a unique snapshot-number fallback.
    """
    exact = (
        simulation_region_dir
        / f"{int(sn):03d}.z{float(redshift):.3f}.AHF_halos"
    )

    if exact.exists():
        return exact

    matches = sorted(
        simulation_region_dir.glob(
            f"{int(sn):03d}.z*.AHF_halos"
        )
    )

    if len(matches) == 1:
        return matches[0]

    if len(matches) == 0:
        raise FileNotFoundError(exact)

    raise RuntimeError(
        f"Se encontraron varios AHF_halos para sn={sn} en "
        f"{simulation_region_dir}:\n"
        + "\n".join(f"  {path}" for path in matches)
    )


def load_halo_positions(halo_file: Path) -> Dict[str, np.ndarray]:
    """
    Load halo IDs and coordinates from AHF_halos.

    NEXUS uses:
        ID -> column 0
        X  -> column 5
        Y  -> column 6
        Z  -> column 7

    Raw positions are divided by 1000, exactly as in NEXUS:
        ckpc/h -> cMpc/h
    """
    try:
        id_halo, x, y, z = np.loadtxt(
            halo_file,
            usecols=(0, 5, 6, 7),
            unpack=True,
        )
    except ValueError as error:
        raise ValueError(
            f"No se pudo leer el catálogo AHF_halos:\n{halo_file}"
        ) from error

    id_halo = np.atleast_1d(id_halo).astype(np.int64)

    return {
        "id_halo": id_halo,
        "x": np.atleast_1d(x).astype(float) / 1.0e3,
        "y": np.atleast_1d(y).astype(float) / 1.0e3,
        "z": np.atleast_1d(z).astype(float) / 1.0e3,
    }


def compute_convex_hull_area(
    x: Sequence[float],
    y: Sequence[float],
) -> float:
    """
    Reproduce NEXUS's projected convex-hull area definition.

    In two dimensions scipy.spatial.ConvexHull.volume is the enclosed area.
    """
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)

    finite = np.isfinite(x_arr) & np.isfinite(y_arr)
    points = np.column_stack(
        [x_arr[finite], y_arr[finite]]
    )

    if len(points) < 3:
        return 0.0

    points = np.unique(points, axis=0)

    if len(points) < 3:
        return 0.0

    try:
        hull = ConvexHull(points)
        return float(hull.volume)
    except (QhullError, ValueError):
        return 0.0


def projected_area_native(
    halo_data: Dict[str, np.ndarray],
    matched_mask: np.ndarray,
    projection: int,
) -> float:
    """
    Convex-hull area in native NEXUS units: (cMpc/h)^2.
    """
    projection = int(projection)

    if projection == 1:
        xproj = halo_data["x"]
        yproj = halo_data["y"]
    elif projection == 2:
        xproj = halo_data["x"]
        yproj = halo_data["z"]
    elif projection == 3:
        xproj = halo_data["y"]
        yproj = halo_data["z"]
    else:
        raise ValueError(
            f"projection no reconocida: {projection}. Usa 1, 2 o 3."
        )

    return compute_convex_hull_area(
        xproj[matched_mask],
        yproj[matched_mask],
    )


def compute_snapshot_row(
    region: int,
    resolution: str,
    mass_label: str,
    step: int,
    catalogue: Dict[str, object],
    simulation_region_dir: Path,
    little_h: float,
) -> Dict[str, object]:
    """
    Compute all three projected truth areas for one snapshot.
    """
    sn = int(catalogue["sn"])
    redshift = float(catalogue["z"])
    truth_file = Path(catalogue["truth_file"])

    ids_truth = load_ids_bz2(truth_file)

    halo_file = resolve_halo_file(
        simulation_region_dir=simulation_region_dir,
        sn=sn,
        redshift=redshift,
    )

    halo_data = load_halo_positions(halo_file)

    matched_mask = np.isin(
        halo_data["id_halo"],
        ids_truth,
        assume_unique=False,
    )

    matched_ids = np.unique(
        halo_data["id_halo"][matched_mask]
    )

    areas_native = {
        projection: projected_area_native(
            halo_data=halo_data,
            matched_mask=matched_mask,
            projection=projection,
        )
        for projection in (1, 2, 3)
    }

    native_values = np.asarray(
        [areas_native[1], areas_native[2], areas_native[3]],
        dtype=float,
    )

    # If coordinates are cMpc/h, an area in h-free cMpc^2 is obtained
    # by dividing by h^2.
    h2 = float(little_h) ** 2

    if h2 <= 0:
        raise ValueError("little_h debe ser positivo.")

    areas_cmpc2 = {
        projection: value / h2
        for projection, value in areas_native.items()
    }

    cmpc_values = np.asarray(
        [areas_cmpc2[1], areas_cmpc2[2], areas_cmpc2[3]],
        dtype=float,
    )

    return {
        "sn": sn,
        "z": redshift,
        "area_projection_median_cMpc2": float(np.median(cmpc_values)),
        "area_projection_mean_cMpc2": float(np.mean(cmpc_values)),
        "area_xy_cMpc2": float(areas_cmpc2[1]),
        "area_xz_cMpc2": float(areas_cmpc2[2]),
        "area_yz_cMpc2": float(areas_cmpc2[3]),
        "area_projection_median_cMpc_h2": float(np.median(native_values)),
        "area_projection_mean_cMpc_h2": float(np.mean(native_values)),
        "area_xy_cMpc_h2": float(areas_native[1]),
        "area_xz_cMpc_h2": float(areas_native[2]),
        "area_yz_cMpc_h2": float(areas_native[3]),
        "N_truth_ids": int(len(ids_truth)),
        "N_matched_halos": int(len(matched_ids)),
        "N_missing_ids": int(len(ids_truth) - len(matched_ids)),
        "region": int(region),
        "resolution": str(resolution),
        "mass_label": str(mass_label),
        "step": int(step),
        "truth_file": str(truth_file),
        "halo_file": str(halo_file),
    }


def make_truth_areas_for_region(
    region: int,
    resolution: str,
    smass_cut: float,
    step: int,
    truth_base: Path,
    simulation_base: Path,
    little_h: float,
    output_name: str | None,
    overwrite: bool,
) -> Path | None:
    """
    Compute one compressed truth-area table for one region.
    """
    mass_label = mass_label_from_cut(smass_cut)

    truth_region_dir = (
        truth_base
        / resolution
        / mass_label
        / f"step{int(step)}"
        / f"region{int(region):04d}"
    )

    simulation_region_dir = (
        simulation_base
        / resolution
        / f"region{int(region)}"
    )

    if output_name is None:
        output_name = f"truth_areas.step{int(step)}.csv.bz2"

    out_file = truth_region_dir / output_name

    if out_file.exists() and not overwrite:
        print(f"  existing output, skipped: {out_file}")
        return out_file

    catalogues = discover_truth_catalogues(
        region_dir=truth_region_dir,
        step=step,
    )

    if len(catalogues) == 0:
        print(
            f"  WARNING: no truth files found in {truth_region_dir}"
        )
        return None

    print("")
    print(f"========== region {int(region):4d} ==========")
    print(f"resolution       : {resolution}")
    print(f"mass_label       : {mass_label}")
    print(f"step             : {step}")
    print(f"N snapshots      : {len(catalogues)}")
    print(f"truth directory  : {truth_region_dir}")
    print(f"simulation dir   : {simulation_region_dir}")
    print(f"h                : {little_h:.6f}")

    rows: List[Dict[str, object]] = []

    for catalogue in catalogues:
        row = compute_snapshot_row(
            region=region,
            resolution=resolution,
            mass_label=mass_label,
            step=step,
            catalogue=catalogue,
            simulation_region_dir=simulation_region_dir,
            little_h=little_h,
        )

        rows.append(row)

        print(
            f"  sn={row['sn']:3d} z={row['z']:6.3f} "
            f"N={row['N_matched_halos']:5d} "
            f"Amed={row['area_projection_median_cMpc2']:12.5f} cMpc^2"
        )

    df = pd.DataFrame(rows)

    # Ensure a deterministic column order and high-z -> low-z evolution.
    columns = [
        "sn",
        "area_projection_median_cMpc2",
        "z",
        "area_projection_mean_cMpc2",
        "area_xy_cMpc2",
        "area_xz_cMpc2",
        "area_yz_cMpc2",
        "area_projection_median_cMpc_h2",
        "area_projection_mean_cMpc_h2",
        "area_xy_cMpc_h2",
        "area_xz_cMpc_h2",
        "area_yz_cMpc_h2",
        "N_truth_ids",
        "N_matched_halos",
        "N_missing_ids",
        "region",
        "resolution",
        "mass_label",
        "step",
        "truth_file",
        "halo_file",
    ]

    df = (
        df[columns]
        .sort_values(["sn"])
        .reset_index(drop=True)
    )

    truth_region_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    df.to_csv(
        out_file,
        index=False,
        compression="bz2",
    )

    print(f"saved: {out_file}")
    print("END")

    return out_file


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compute NEXUS truth convex-hull areas from compact truth-ID files."
        )
    )

    parser.add_argument(
        "resolution",
        choices=["3k", "7k"],
        help="Simulation resolution.",
    )

    parser.add_argument(
        "--regions",
        type=int,
        nargs="+",
        default=None,
        help="Regions to process. Default: all valid regions.",
    )

    parser.add_argument(
        "--smass-cut",
        type=float,
        default=1.0e10,
        help="Stellar-mass cut used by the truth catalogues.",
    )

    parser.add_argument(
        "--step",
        type=int,
        choices=[1, 2],
        default=1,
        help="Truth definition: 1=extended progenitors, 2=core.",
    )

    parser.add_argument(
        "--truth-base",
        type=Path,
        default=Path(
            "/data9/users/jsgomez/protocluster_truth"
        ),
        help="Base directory containing compact truth-ID catalogues.",
    )

    parser.add_argument(
        "--simulation-base",
        type=Path,
        default=Path(
            "/data9/users/jsgomez/simulations/"
            "TheThreeHundred/HYDROs/GIZMO/music_planck"
        ),
        help="Base directory containing resolution/region*/AHF_halos files.",
    )

    parser.add_argument(
        "--h",
        type=float,
        default=0.6777,
        help=(
            "Dimensionless Hubble parameter used to convert "
            "(cMpc/h)^2 to cMpc^2. Set --h 1 to keep native areas unchanged."
        ),
    )

    parser.add_argument(
        "--output-name",
        type=str,
        default=None,
        help=(
            "Output filename inside each region directory. "
            "Default: truth_areas.stepN.csv.bz2"
        ),
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite an existing output table.",
    )

    args = parser.parse_args()

    regions = (
        default_regions(args.resolution)
        if args.regions is None
        else args.regions
    )

    print("")
    print("========== make_truth_areas.py ==========")
    print(f"resolution      : {args.resolution}")
    print(f"N regions       : {len(regions)}")
    print(f"smass_cut       : {args.smass_cut:.3e}")
    print(f"step            : {args.step}")
    print(f"truth_base      : {args.truth_base}")
    print(f"simulation_base : {args.simulation_base}")
    print(f"h               : {args.h:.6f}")
    print("=========================================")

    n_saved = 0
    n_skipped = 0

    for region in regions:
        try:
            result = make_truth_areas_for_region(
                region=region,
                resolution=args.resolution,
                smass_cut=args.smass_cut,
                step=args.step,
                truth_base=args.truth_base,
                simulation_base=args.simulation_base,
                little_h=args.h,
                output_name=args.output_name,
                overwrite=args.overwrite,
            )

            if result is None:
                n_skipped += 1
            else:
                n_saved += 1

        except (FileNotFoundError, RuntimeError, ValueError) as error:
            n_skipped += 1
            print("")
            print(f"WARNING: region {region} skipped:")
            print(error)
            print("")

    print("")
    print("========== SUMMARY ==========")
    print(f"regions requested : {len(regions)}")
    print(f"tables available  : {n_saved}")
    print(f"regions skipped   : {n_skipped}")
    print("END")


if __name__ == "__main__":
    main()
