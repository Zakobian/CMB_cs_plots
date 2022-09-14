from main_parser import parser
import plotter
import numpy as np
from BaseUniverse import universe

args = parser.parse_args()
universes=[universe(1),universe(0),universe(-1)]
cs2ms=[0.5]
cs2ps=[0.5]
print(cs2ps,cs2ms)
for experim,(cs2p,cs2m) in enumerate(zip(cs2ps,cs2ms)):
    plotter.plot_1(universes,np.array(args.nts),cs2m,cs2p,experim)
