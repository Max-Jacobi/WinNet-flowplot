import numpy as np
import os 

datapath = os.path.dirname(__file__) +'/data'
Elnames = np.loadtxt("{}/elementlist".format(datapath), dtype=str)
St_N, St_Z = np.loadtxt('{}/stableiso.dat'.format(datapath), unpack=True, usecols=(1, 2))
