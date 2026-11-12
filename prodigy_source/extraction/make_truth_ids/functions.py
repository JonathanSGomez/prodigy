#!/usr/bin/env python3
import glob
import os
import time


def check_endfile(logfiles):
    """
    Count logs containing an exact END line.
    """
    n_ended = 0

    for logfile in logfiles:
        try:
            with open(logfile, "rt") as stream:
                for line in stream:
                    if line.strip() == "END":
                        n_ended += 1
                        break
        except FileNotFoundError:
            pass

    return n_ended


def print_time(seconds):
    seconds = int(seconds)

    if seconds < 60:
        return "time: %ds" % seconds

    if seconds < 3600:
        minutes = seconds // 60
        remaining = seconds - minutes * 60
        return "time: %dm %ds" % (
            minutes,
            remaining,
        )

    hours = seconds // 3600
    minutes = (
        seconds - hours * 3600
    ) // 60

    remaining = (
        seconds
        - hours * 3600
        - minutes * 60
    )

    return "time: %dh %dm %ds" % (
        hours,
        minutes,
        remaining,
    )


def progress_bar(
    dir_log,
    n_total,
    code_name,
    log_pattern="*.log",
    barsize=50,
    wait_time=1,
):
    """
    Wait until n_total logs contain an exact END line.
    """
    if n_total <= 0:
        return

    ti = time.time()

    while True:
        logfiles = glob.glob(
            os.path.join(
                dir_log,
                log_pattern,
            )
        )

        n_ended = check_endfile(
            logfiles
        )

        fraction = min(
            max(
                n_ended / float(n_total),
                0.0,
            ),
            1.0,
        )

        n_bar = int(
            fraction * barsize
        )

        initial = "#" * n_bar
        final = "-" * (
            barsize - n_bar
        )

        print(
            "[%s%s] | %3d de %3d | %s | codigo %s"
            % (
                initial,
                final,
                n_ended,
                n_total,
                print_time(
                    time.time() - ti
                ),
                code_name,
            ),
            end="\r",
            flush=True,
        )

        if n_ended >= n_total:
            break

        time.sleep(wait_time)

    print("")
