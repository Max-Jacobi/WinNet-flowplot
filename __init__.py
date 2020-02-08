import numpy as np

Elnames = np.loadtxt("/home/mjacobi/Winnet/Winnet_scripts/my_scripts/data/elementlist", dtype=str)
St_N, St_Z = np.loadtxt('/home/mjacobi/Winnet/Winnet_scripts/my_scripts/data/stableiso.dat', unpack=True, usecols=(1, 2))
