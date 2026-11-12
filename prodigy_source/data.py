from functions import *

resolution  = str(sys.argv[1])
dir_gizmo   = "/data9/users/jsgomez/simulations/TheThreeHundred/HYDROs/GIZMO/music_planck/%s" % resolution
dir_results = "/data9/users/jsgomez/protocluster_candidates"
zpath       = "/data9/users/jsgomez/simulations/TheThreeHundred"
redshifts   = np.absolute(np.loadtxt("%s/z_list.txt" % zpath))

#id_prog_centrals = pd.read_csv('initial_tracks/GIZMO_centralhalosIDS_noCRM.txt', sep='\t')

if resolution == "3k": id_prog_centrals = pd.read_csv('initial_tracks/centralhalosIDS_gizmo_3k_nonskipping.txt', sep=" ") # es mas completo que CRMratio2
if resolution == "7k": id_prog_centrals = pd.read_csv('initial_tracks/centralhalosIDS_gizmo_7k_CRMratio2.txt'  , sep=" ") # no existe el nonskipping

galaxies_dict    = {}
densities_dict   = {}
