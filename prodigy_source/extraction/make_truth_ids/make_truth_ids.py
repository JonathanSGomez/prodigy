#!/usr/bin/env python3
"""
Generate compact NEXUS truth-ID catalogues for both physical definitions:

Step 1
------
All real progenitors of the z=0 cluster members.

Step 2
------
The core/main-cluster truth catalogue used by NEXUS:
the main progenitor of the z=0 central halo plus its subhalos inside
the main progenitor R200 at each snapshot.

The Step-1 and Step-2 implementations follow the definitions in
NEXUS metric.py.

Examples
--------
# Generate both Step 1 and Step 2 for region 70
python make_truth_ids.py 7k \
    --regions 70 \
    --snapshots 128 127 126 125 \
    --smass-cut 1.0e10 \
    --steps 1 2

# Generate only Step 2
python make_truth_ids.py 7k \
    --regions 70 \
    --snapshots 96 \
    --smass-cut 1.0e10 \
    --steps 2
"""

from __future__ import annotations

import argparse
import bz2
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, Sequence

import numpy as np


def read_redshifts(zpath: str | Path) -> np.ndarray:
    return np.absolute(
        np.loadtxt(
            Path(zpath) / "z_list.txt"
        )
    )


def get_fmhires_cut(resolution: str) -> float:
    if resolution == "3k":
        return 0.996

    if resolution == "7k":
        return 0.9995

    raise ValueError(
        f"resolution no reconocida: {resolution}"
    )


def mass_label_from_cut(smass_cut: float) -> str:
    if smass_cut <= 0:
        raise ValueError("smass_cut debe ser positivo.")

    exponent = int(round(np.log10(float(smass_cut))))

    if not np.isclose(
        float(smass_cut),
        10.0 ** exponent,
    ):
        raise ValueError(
            "Este flujo espera cortes de masa 10**N. "
            f"Valor recibido: {smass_cut:.6e}"
        )

    return f"Mstar1e{exponent}"


def default_regions(resolution: str) -> list[int]:
    if resolution == "7k":
        return [
            region
            for region in range(1, 325)
            if region != 107
        ]

    if resolution == "3k":
        return list(range(1, 325))

    raise ValueError(
        f"resolution no reconocida: {resolution}"
    )


def ensure_1d_int(values) -> np.ndarray:
    values = np.asarray(values)

    if values.size == 0:
        return np.array([], dtype=np.int64)

    return np.atleast_1d(values).astype(np.int64)


def load_gizmo_properties(
    dir_gizmo_resolution: str | Path,
    region: int,
    sn: int,
    redshift: float,
) -> Dict[str, np.ndarray]:
    """
    Load the AHF_halos columns used by the NEXUS truth definitions.
    """
    file_halos = (
        Path(dir_gizmo_resolution)
        / f"region{int(region)}"
        / f"{int(sn):03d}.z{float(redshift):.3f}.AHF_halos"
    )

    if not file_halos.exists():
        raise FileNotFoundError(file_halos)

    props: Dict[str, np.ndarray] = {}

    (
        props["id_halo"],
        props["smass"],
        props["hmass"],
        props["id_hosthalo"],
        props["x"],
        props["y"],
        props["z"],
        props["r200"],
        props["fmhires"],
    ) = np.loadtxt(
        file_halos,
        usecols=(0, 64, 3, 1, 5, 6, 7, 11, 37),
        unpack=True,
    )

    props["id_halo"] = (
        np.atleast_1d(props["id_halo"])
        .astype(np.int64)
    )

    props["id_hosthalo"] = (
        np.atleast_1d(props["id_hosthalo"])
        .astype(np.int64)
    )

    for key in [
        "smass",
        "hmass",
        "x",
        "y",
        "z",
        "r200",
        "fmhires",
    ]:
        props[key] = np.atleast_1d(
            props[key]
        ).astype(float)

    # AHF coordinates/radii: ckpc/h -> cMpc/h.
    for key in ["x", "y", "z", "r200"]:
        props[key] = props[key] / 1.0e3

    return props


def recenter_properties(
    props: Dict[str, np.ndarray],
    x0: float,
    y0: float,
    z0: float,
) -> Dict[str, np.ndarray]:
    """
    Recenter positions relative to the z=0 main halo, as in NEXUS.
    """
    out = dict(props)

    out["x"] = props["x"] - x0
    out["y"] = props["y"] - y0
    out["z"] = props["z"] - z0

    return out


# ============================================================
# Step 1: all real progenitors
# ============================================================
def get_z0_cluster_members(
    props128: Dict[str, np.ndarray],
    smass_cut_for_roots: float = 0.0,
    host_index: int = 0,
):
    """
    z=0 central halo plus its subhalos inside its R200.

    These IDs become the roots in the full-progenitor nodes file.
    """
    if host_index < 0 or host_index >= len(props128["id_halo"]):
        raise IndexError(
            f"host_index={host_index} fuera del catálogo z=0."
        )

    host_id = int(
        props128["id_halo"][host_index]
    )

    xh = props128["x"][host_index]
    yh = props128["y"][host_index]
    zh = props128["z"][host_index]
    r200 = props128["r200"][host_index]

    p_sub = np.where(
        props128["id_hosthalo"] == host_id
    )[0]

    dist = np.sqrt(
        (props128["x"][p_sub] - xh) ** 2
        + (props128["y"][p_sub] - yh) ** 2
        + (props128["z"][p_sub] - zh) ** 2
    )

    p_sub = p_sub[dist < r200]

    p_members = np.unique(
        np.concatenate(
            [
                np.array([host_index], dtype=int),
                p_sub,
            ]
        )
    )

    if smass_cut_for_roots > 0:
        p_members = p_members[
            props128["smass"][p_members]
            > smass_cut_for_roots
        ]

    return (
        props128["id_halo"][p_members]
        .astype(np.int64),
        p_members,
    )


def load_full_nodes(
    nodes_file: str | Path,
) -> np.ndarray:
    """
    Load the columns used by NEXUS from the full progenitor nodes file.
    """
    nodes_file = Path(nodes_file)

    if not nodes_file.exists():
        raise FileNotFoundError(nodes_file)

    nodes = np.loadtxt(
        nodes_file,
        usecols=(0, 1, 2, 4),
        dtype=[
            ("root", "i8"),
            ("sn", "i8"),
            ("id", "i8"),
            ("smass", "f8"),
        ],
    )

    return np.atleast_1d(nodes)


def get_all_progenitors_ids(
    nodes: np.ndarray,
    root_ids_z0: np.ndarray,
    sn: int,
    smass_cut: float,
    props_sn: Dict[str, np.ndarray] | None = None,
    fmhires_cut: float | None = None,
) -> np.ndarray:
    """
    Step-1 truth IDs, matching NEXUS get_all_progenitors_ids().
    """
    mask = (
        np.isin(
            nodes["root"],
            root_ids_z0,
        )
        & (nodes["sn"] == int(sn))
        & (nodes["smass"] > float(smass_cut))
    )

    ids = np.unique(
        nodes["id"][mask]
        .astype(np.int64)
    )

    if (
        props_sn is not None
        and fmhires_cut is not None
    ):
        clean = (
            props_sn["fmhires"]
            > float(fmhires_cut)
        )

        ids_clean = (
            props_sn["id_halo"][clean]
            .astype(np.int64)
        )

        ids = np.intersect1d(
            ids,
            ids_clean,
        )

    return ensure_1d_int(ids)


# ============================================================
# Step 2: main-progenitor core
# ============================================================
def _load_track_array(
    path: Path,
    integer: bool = False,
) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(path)

    values = np.loadtxt(path)

    # NEXUS indexes tracks as tracks[prop][host_index, sn].
    # Preserve a row dimension even if a file has a single track.
    values = np.asarray(values)

    if values.ndim == 1:
        values = values[np.newaxis, :]

    if integer:
        values = values.astype(np.int64)
    else:
        values = values.astype(float)

    return values


def load_old_main_progenitor_tracks(
    dir_gizmo_resolution: str | Path,
    region: int,
) -> Dict[str, np.ndarray]:
    """
    Load the old main-progenitor tracks used by NEXUS Step 2.

    Files:
        128.z0.000.AHF_progenitors_idhalo
        128.z0.000.AHF_progenitors_x
        128.z0.000.AHF_progenitors_y
        128.z0.000.AHF_progenitors_z
        128.z0.000.AHF_progenitors_r200
        128.z0.000.AHF_progenitors_smass
    """
    base = (
        Path(dir_gizmo_resolution)
        / f"region{int(region)}"
    )

    tracks: Dict[str, np.ndarray] = {}

    for prop in [
        "idhalo",
        "x",
        "y",
        "z",
        "r200",
        "smass",
    ]:
        path = (
            base
            / f"128.z0.000.AHF_progenitors_{prop}"
        )

        tracks[prop] = _load_track_array(
            path,
            integer=(prop == "idhalo"),
        )

        if prop in ["x", "y", "z", "r200"]:
            tracks[prop] = (
                tracks[prop] / 1.0e3
            )

    return tracks


def recenter_tracks(
    tracks: Dict[str, np.ndarray],
    x0: float,
    y0: float,
    z0: float,
) -> Dict[str, np.ndarray]:
    """
    Recenter Step-2 track coordinates exactly as in NEXUS.
    """
    out = dict(tracks)

    out["x"] = tracks["x"] - x0
    out["y"] = tracks["y"] - y0
    out["z"] = tracks["z"] - z0

    return out


def validate_track_index(
    tracks: Dict[str, np.ndarray],
    host_index: int,
    sn: int,
) -> None:
    shape = tracks["idhalo"].shape

    if host_index < 0 or host_index >= shape[0]:
        raise IndexError(
            f"host_index={host_index} fuera de tracks con shape={shape}."
        )

    if sn < 0 or sn >= shape[1]:
        raise IndexError(
            f"sn={sn} fuera de tracks con shape={shape}."
        )


def get_main_cluster_ids(
    props_sn: Dict[str, np.ndarray],
    tracks: Dict[str, np.ndarray],
    host_index: int,
    sn: int,
    smass_cut: float,
    fmhires_cut: float | None = None,
) -> np.ndarray:
    """
    Step-2 truth catalogue, matching NEXUS get_main_cluster_ids().

    Definition
    ----------
    Main progenitor of the z=0 main halo plus its subhalos inside
    the main progenitor R200, with the stellar-mass and high-resolution
    cuts applied.
    """
    validate_track_index(
        tracks=tracks,
        host_index=host_index,
        sn=sn,
    )

    main_id = int(
        tracks["idhalo"][host_index, sn]
    )

    if main_id <= 0:
        return np.array([], dtype=np.int64)

    xc = tracks["x"][host_index, sn]
    yc = tracks["y"][host_index, sn]
    zc = tracks["z"][host_index, sn]
    rc = tracks["r200"][host_index, sn]

    p_sub = np.where(
        props_sn["id_hosthalo"] == main_id
    )[0]

    p_cen = np.where(
        props_sn["id_halo"] == main_id
    )[0]

    p = np.unique(
        np.concatenate(
            [
                p_cen,
                p_sub,
            ]
        )
    )

    if len(p) == 0:
        return np.array([], dtype=np.int64)

    dist = np.sqrt(
        (props_sn["x"][p] - xc) ** 2
        + (props_sn["y"][p] - yc) ** 2
        + (props_sn["z"][p] - zc) ** 2
    )

    mask = (
        (dist < rc)
        & (
            props_sn["smass"][p]
            > float(smass_cut)
        )
    )

    if fmhires_cut is not None:
        mask &= (
            props_sn["fmhires"][p]
            > float(fmhires_cut)
        )

    p = p[mask]

    return np.unique(
        props_sn["id_halo"][p]
        .astype(np.int64)
    )


# ============================================================
# Output
# ============================================================
def truth_output_path(
    out_base: str | Path,
    resolution: str,
    mass_label: str,
    step: int,
    region: int,
    sn: int,
    redshift: float,
) -> Path:
    out_dir = (
        Path(out_base)
        / resolution
        / mass_label
        / f"step{int(step)}"
        / f"region{int(region):04d}"
    )

    return (
        out_dir
        / (
            f"truth.sn{int(sn):03d}."
            f"z{float(redshift):.3f}."
            f"step{int(step)}.txt.bz2"
        )
    )


def save_ids_bz2(
    ids,
    out_file: str | Path,
    overwrite: bool = False,
) -> bool:
    """
    Save IDs. Return True when written and False when skipped.
    """
    out_file = Path(out_file)

    if out_file.exists() and not overwrite:
        return False

    out_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    ids = ensure_1d_int(ids)

    with bz2.open(out_file, "wt") as stream:
        for halo_id in ids:
            stream.write(
                f"{int(halo_id)}\n"
            )

    return True


# ============================================================
# Region workflow
# ============================================================
def make_truth_ids_for_region(
    region: int,
    resolution: str,
    snapshots: Sequence[int],
    smass_cut: float,
    steps: Sequence[int],
    z0_root_smass_cut: float,
    host_index: int,
    dir_gizmo_base: str | Path,
    zpath: str | Path,
    out_base: str | Path,
    overwrite: bool = False,
) -> None:
    steps = sorted(
        set(int(step) for step in steps)
    )

    invalid_steps = [
        step
        for step in steps
        if step not in [1, 2]
    ]

    if invalid_steps:
        raise ValueError(
            f"steps no reconocidos: {invalid_steps}"
        )

    if len(steps) == 0:
        raise ValueError(
            "Debes solicitar Step 1, Step 2, o ambos."
        )

    redshifts = read_redshifts(zpath)

    dir_gizmo_resolution = (
        Path(dir_gizmo_base)
        / resolution
    )

    fmhires_cut = get_fmhires_cut(
        resolution
    )

    mass_label = mass_label_from_cut(
        smass_cut
    )

    print("")
    print(
        f"========== region {int(region):4d} =========="
    )
    print(f"resolution  : {resolution}")
    print(f"smass_cut   : {smass_cut:.3e}")
    print(f"steps       : {steps}")
    print(f"fmhires_cut : {fmhires_cut:.5f}")
    print(f"overwrite   : {overwrite}")

    # z=0 main halo position used by both NEXUS truth definitions.
    props128_raw = load_gizmo_properties(
        dir_gizmo_resolution=dir_gizmo_resolution,
        region=region,
        sn=128,
        redshift=redshifts[128],
    )

    if (
        host_index < 0
        or host_index >= len(
            props128_raw["id_halo"]
        )
    ):
        raise IndexError(
            f"host_index={host_index} fuera del catálogo z=0."
        )

    x0 = props128_raw["x"][host_index]
    y0 = props128_raw["y"][host_index]
    z0 = props128_raw["z"][host_index]

    props128 = recenter_properties(
        props128_raw,
        x0,
        y0,
        z0,
    )

    root_ids_z0 = None
    nodes = None
    tracks = None

    if 1 in steps:
        root_ids_z0, _ = get_z0_cluster_members(
            props128,
            smass_cut_for_roots=z0_root_smass_cut,
            host_index=host_index,
        )

        print(
            "z=0 cluster members used as Step-1 roots: "
            f"{len(root_ids_z0)}"
        )

        nodes_file = (
            dir_gizmo_resolution
            / f"region{int(region)}"
            / "128.AHF_full_progenitors_nodes.bz2"
        )

        print(
            f"loading Step-1 full progenitor nodes: "
            f"{nodes_file}"
        )

        nodes = load_full_nodes(
            nodes_file
        )

    if 2 in steps:
        print(
            "loading old main-progenitor tracks "
            "for Step-2 truth core"
        )

        tracks = load_old_main_progenitor_tracks(
            dir_gizmo_resolution=dir_gizmo_resolution,
            region=region,
        )

        tracks = recenter_tracks(
            tracks,
            x0,
            y0,
            z0,
        )

    for sn in snapshots:
        sn = int(sn)

        if sn < 0 or sn >= len(redshifts):
            raise IndexError(
                f"sn={sn} fuera de z_list con longitud {len(redshifts)}."
            )

        redshift = float(
            redshifts[sn]
        )

        props_sn_raw = load_gizmo_properties(
            dir_gizmo_resolution=dir_gizmo_resolution,
            region=region,
            sn=sn,
            redshift=redshift,
        )

        props_sn = recenter_properties(
            props_sn_raw,
            x0,
            y0,
            z0,
        )

        if 1 in steps:
            assert nodes is not None
            assert root_ids_z0 is not None

            ids_step1 = get_all_progenitors_ids(
                nodes=nodes,
                root_ids_z0=root_ids_z0,
                sn=sn,
                smass_cut=smass_cut,
                props_sn=props_sn,
                fmhires_cut=fmhires_cut,
            )

            out_step1 = truth_output_path(
                out_base=out_base,
                resolution=resolution,
                mass_label=mass_label,
                step=1,
                region=region,
                sn=sn,
                redshift=redshift,
            )

            written = save_ids_bz2(
                ids_step1,
                out_step1,
                overwrite=overwrite,
            )

            state = (
                "written"
                if written
                else "existing, skipped"
            )

            print(
                f"  sn={sn:3d} z={redshift:.3f} "
                f"Step 1 N_truth={len(ids_step1):6d} "
                f"[{state}] -> {out_step1}"
            )

        if 2 in steps:
            assert tracks is not None

            ids_step2 = get_main_cluster_ids(
                props_sn=props_sn,
                tracks=tracks,
                host_index=host_index,
                sn=sn,
                smass_cut=smass_cut,
                fmhires_cut=fmhires_cut,
            )

            out_step2 = truth_output_path(
                out_base=out_base,
                resolution=resolution,
                mass_label=mass_label,
                step=2,
                region=region,
                sn=sn,
                redshift=redshift,
            )

            written = save_ids_bz2(
                ids_step2,
                out_step2,
                overwrite=overwrite,
            )

            state = (
                "written"
                if written
                else "existing, skipped"
            )

            print(
                f"  sn={sn:3d} z={redshift:.3f} "
                f"Step 2 N_truth={len(ids_step2):6d} "
                f"[{state}] -> {out_step2}"
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Generate compact NEXUS truth IDs for "
            "Step 1, Step 2, or both."
        )
    )

    parser.add_argument(
        "resolution",
        choices=["3k", "7k"],
    )

    parser.add_argument(
        "--regions",
        type=int,
        nargs="+",
        default=None,
    )

    parser.add_argument(
        "--snapshots",
        type=int,
        nargs="+",
        default=[96],
    )

    parser.add_argument(
        "--smass-cut",
        type=float,
        default=1.0e10,
    )

    parser.add_argument(
        "--steps",
        type=int,
        nargs="+",
        choices=[1, 2],
        default=[1, 2],
        help="Truth definitions to write. Default: 1 2.",
    )

    parser.add_argument(
        "--z0-root-smass-cut",
        type=float,
        default=0.0,
        help=(
            "Optional z=0 root mass cut used only by Step 1."
        ),
    )

    parser.add_argument(
        "--host-index",
        type=int,
        default=0,
    )

    parser.add_argument(
        "--dir-gizmo-base",
        type=Path,
        default=Path(
            "/data9/users/jsgomez/simulations/"
            "TheThreeHundred/HYDROs/GIZMO/music_planck"
        ),
    )

    parser.add_argument(
        "--zpath",
        type=Path,
        default=Path(
            "/data9/users/jsgomez/simulations/TheThreeHundred"
        ),
    )

    parser.add_argument(
        "--out-base",
        type=Path,
        default=Path(
            "/data9/users/jsgomez/protocluster_truth"
        ),
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
    )

    args = parser.parse_args()

    regions = (
        default_regions(args.resolution)
        if args.regions is None
        else args.regions
    )

    print("")
    print("========== make_truth_ids.py ==========")
    print(f"resolution : {args.resolution}")
    print(f"N regions  : {len(regions)}")
    print(f"snapshots  : {args.snapshots}")
    print(f"smass_cut  : {args.smass_cut:.3e}")
    print(f"steps      : {args.steps}")
    print(f"out_base   : {args.out_base}")
    print(f"overwrite  : {args.overwrite}")
    print("=======================================")

    n_failed = 0

    for region in regions:
        try:
            make_truth_ids_for_region(
                region=region,
                resolution=args.resolution,
                snapshots=args.snapshots,
                smass_cut=args.smass_cut,
                steps=args.steps,
                z0_root_smass_cut=args.z0_root_smass_cut,
                host_index=args.host_index,
                dir_gizmo_base=args.dir_gizmo_base,
                zpath=args.zpath,
                out_base=args.out_base,
                overwrite=args.overwrite,
            )

        except Exception as error:
            n_failed += 1

            print("")
            print(
                f"ERROR: region {region} failed:"
            )
            print(error)
            print("")

    # Exact END line required by the progress-bar runner.
    print("END")

    return 1 if n_failed > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
