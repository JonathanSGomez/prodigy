import sys
sys.dont_write_bytecode = True

import os
from functions import progress_bar


resolution = "7k"

dir_log = "/home/jsgomez/codes/nexus/logfiles"
os.makedirs(dir_log, exist_ok=True)

n_max_regions = 108

# Para GIZMO 7k
to_remove = [107]
all_regions = [x for x in range(1, 325) if x not in to_remove]

subregions = [
    all_regions[i:i+n_max_regions]
    for i in range(0, len(all_regions), n_max_regions)
]

for i_batch, regions in enumerate(subregions):

    print("")
    print("====================================")
    print("Batch %d / %d" % (i_batch + 1, len(subregions)))
    print("Regions:", regions[0], "to", regions[-1])
    print("N regions:", len(regions))
    print("====================================")
    print("")

    # Limpiar logs anteriores
    if os.listdir(dir_log):
        os.system("rm %s/*" % dir_log)

    n_regions = 0

    for region in regions:

        logfile = "%s/nexus.%s.region%s.log" % (
            dir_log,
            resolution,
            str(region).zfill(4)
        )

        cmd = (
            "python main.py %s %d > %s 2>&1 &" %
            (resolution, region, logfile)
        )

        os.system(cmd)
        n_regions += 1

    progress_bar(dir_log, n_regions, "NEXUS")

    print("END batch %d" % (i_batch + 1))
    print("")
