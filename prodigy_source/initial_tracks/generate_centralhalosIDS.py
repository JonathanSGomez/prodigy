import os
import numpy as np

source = "gizmo"
resolution = "7k"

final_array = [] 
for region in range(1, 325):
    print("region: %s" % region)
    dir_gizmo  = "/data9/users/jsgomez/simulations/TheThreeHundred/HYDROs/GIZMO/music_planck/%s/region%s" % (resolution, region)

    if os.path.exists("%s/128.z0.000.AHF_progenitors_idhalo" % dir_gizmo):
        IDs = np.loadtxt("%s/128.z0.000.AHF_progenitors_idhalo" % dir_gizmo)
        IDs = np.array(IDs, dtype=int)

        centralhalosIDS = np.concatenate([ [region] , IDs[0, :][::-1] ]) #; print(len(centralhalosIDS))
    else:
        centralhalosIDS = np.array([-1] * 130)

    final_array.append(centralhalosIDS)

final_array = np.array(final_array, dtype=int)
#print(final_array)

if resolution == "3k": outfile = "centralhalosIDS_%s_%s_nonskipping.txt" % (source, resolution) # es mas completo que CRMratio2
if resolution == "7k": outfile = "centralhalosIDS_%s_%s_CRMratio2.txt"   % (source, resolution) # no existe el nonskipping
np.savetxt(outfile,
           final_array,
           fmt="%d", # %d17
           #header="regionID            128               127               126               125               124               123               122               121               120               119               118               117               116               115               114               113               112               111               110               109               108               107               106               105               104               103               102               101               100                99                98                97                96                95                94                93                92                91                90                89                88                87                86                85                84                83                82                81                80                79                78                77                76                75                74                73                72                71                70                69                68                67                66                65                64                63                62                61                60                59                58                57                56                55                54                53                52                51                50                49                48                47                46                45                44                43                42                41                40                39                38                37                36                35                34                33                32                31                30                29                28                27                26                25                24                23                22                21                20                19                18                17                16                15                14                13                12                11                10                9                 8                 7                 6                 5                 4                 3                 2                 1                 0",
           header = "regionID 128 127 126 125 124 123 122 121 120 119 118 117 116 115 114 113 112 111 110 109 108 107 106 105 104 103 102 101 100 99 98 97 96 95 94 93 92 91 90 89 88 87 86 85 84 83 82 81 80 79 78 77 76 75 74 73 72 71 70 69 68 67 66 65 64 63 62 61 60 59 58 57 56 55 54 53 52 51 50 49 48 47 46 45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0",
           comments=''
           )
