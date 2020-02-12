from .. import np
from ..isotope import Isotope
from .isotopecollection import IsotopeCollection


class Snapshot(IsotopeCollection):
    '''
    Get isotopes from snapshot and create a isotopecollection
    Input:
       path     - path to snapfile

    Atributes:
       path, num, time, temp, dens, isotopes
    '''

    def __init__(self, path):
        super(Snapshot, self).__init__()
        Ns, Zs, Ys = np.loadtxt(path, skiprows=3, usecols=(0, 1, 2), unpack=True)

        with open(path, 'r') as sf:
            sf.readline()
            header = sf.readline()
            self.time, self.temp, self.dens = np.array(header.split()).astype(float)

        self.isotopes = np.array([Isotope(Z=Z, N=N, Y=Y) for N, Z, Y in zip(Ns, Zs, Ys)])
        self.path = path
        self.num = int(path.split("_")[-1][:4])
