import numpy as np
import os

datapath = os.path.dirname(__file__) + '/data'
Elnames = np.loadtxt("{}/elementlist".format(datapath), dtype=str)
St_N, St_Z = np.loadtxt('{}/stableiso.dat'.format(datapath), unpack=True, usecols=(1, 2))


def geName(N, Z):
    '''
    return name of isotpe from its neutron and proton numbers
    '''
    el = Elnames[Z]
    A = Z + N
    return el + str(A)


def getNZ(name):
    '''
    return neutron and proton number for the name of an isotpe
    '''
    name = name.lower()
    el = ''.join([let for let in name if let.isalpha()])
    A = int(''.join([let for let in name if let.isdigit()]))
    Z = int(np.argwhere(Elnames == el))
    N = A - Z
    return N, Z
