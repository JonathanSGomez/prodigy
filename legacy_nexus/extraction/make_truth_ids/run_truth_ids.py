#!/usr/bin/env python3
import sys
sys.dont_write_bytecode = True

import math
import subprocess
from pathlib import Path

from functions import progress_bar


# ============================================================
# Configuration
# ============================================================
resolution = "3k"

# Generate both definitions by default.
steps = [1, 2]

# Full range used previously.
snapshots = list( range(128, 20, -1) )

smass_cut = 1.0e9
z0_root_smass_cut = 0.0
host_index = 0

# Existing Step-1 files are skipped when False, while missing
# Step-2 files are generated.
overwrite = True

# One region = one independent Python process.
n_max_regions = 54 #108

script_dir = (
    Path(__file__)
    .resolve()
    .parent
)

code = (
    script_dir
    / "make_truth_ids.py"
)

mass_label = (
    "Mstar1e%d"
    % int(
        round(
            math.log10(
                smass_cut
            )
        )
    )
)

steps_label = "-".join(
    str(step)
    for step in steps
)

dir_log = (
    script_dir
    / "logfiles_truth"
    / resolution
    / mass_label
    / (
        "steps%s"
        % steps_label
    )
)

dir_log.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Region list
# ============================================================
if resolution == "7k":
    all_regions = [
        region
        for region in range(1, 325)
        if region != 107
    ]

elif resolution == "3k":
    all_regions = list(
        range(1, 325)
    )

else:
    raise ValueError(
        "resolution no reconocida: %s"
        % resolution
    )


subregions = [
    all_regions[
        i:i + n_max_regions
    ]
    for i in range(
        0,
        len(all_regions),
        n_max_regions,
    )
]


def clean_logs(directory: Path) -> None:
    for path in directory.glob(
        "*.log"
    ):
        path.unlink()


# ============================================================
# Run batches
# ============================================================
if not code.exists():
    raise FileNotFoundError(
        "No existe make_truth_ids.py:\n%s"
        % code
    )


for i_batch, regions in enumerate(
    subregions
):
    print("")
    print("====================================")
    print(
        "Batch %d / %d"
        % (
            i_batch + 1,
            len(subregions),
        )
    )
    print(
        "Regions:",
        regions[0],
        "to",
        regions[-1],
    )
    print(
        "N regions:",
        len(regions),
    )
    print(
        "resolution:",
        resolution,
    )
    print(
        "snapshots:",
        snapshots,
    )
    print(
        "smass_cut: %.3e"
        % smass_cut
    )
    print(
        "steps:",
        steps,
    )
    print(
        "overwrite:",
        overwrite,
    )
    print("====================================")
    print("")

    clean_logs(
        dir_log
    )

    processes = []
    log_handles = []

    for region in regions:
        logfile = (
            dir_log
            / (
                "truth.%s.%s.steps%s.region%s.log"
                % (
                    resolution,
                    mass_label,
                    steps_label,
                    str(region).zfill(4),
                )
            )
        )

        command = [
            sys.executable,
            str(code),
            resolution,
            "--regions",
            str(region),
            "--snapshots",
            *[
                str(sn)
                for sn in snapshots
            ],
            "--smass-cut",
            "%.8e" % smass_cut,
            "--steps",
            *[
                str(step)
                for step in steps
            ],
            "--z0-root-smass-cut",
            "%.8e"
            % z0_root_smass_cut,
            "--host-index",
            str(host_index),
        ]

        if overwrite:
            command.append(
                "--overwrite"
            )

        log_handle = open(
            logfile,
            "wt",
        )

        process = subprocess.Popen(
            command,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            cwd=str(script_dir),
        )

        processes.append(
            (
                region,
                process,
                logfile,
            )
        )

        log_handles.append(
            log_handle
        )

    progress_bar(
        dir_log=str(dir_log),
        n_total=len(regions),
        code_name="truth IDs Step 1+2",
        log_pattern=(
            "truth.%s.%s.steps%s.region*.log"
            % (
                resolution,
                mass_label,
                steps_label,
            )
        ),
    )

    failed_regions = []

    for (
        region,
        process,
        logfile,
    ) in processes:
        return_code = (
            process.wait()
        )

        if return_code != 0:
            failed_regions.append(
                region
            )

            print(
                "FAILED region %s. See %s"
                % (
                    str(region).zfill(4),
                    logfile,
                )
            )

    for log_handle in log_handles:
        log_handle.close()

    print("")
    print(
        "END batch %d"
        % (
            i_batch + 1
        )
    )

    if failed_regions:
        print(
            "Failed regions:",
            failed_regions,
        )
    else:
        print(
            "Failed regions: none"
        )

    print("")
