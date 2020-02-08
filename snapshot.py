from . import np
from .nucleus import nucleus


class snapshot(object):
    '''
    Input:
       path     - path to snapfile

    Atributes:
       path, num, time, temp, dens, nuclei
    '''

    def __init__(self, path):
        Ns, Zs, Ys = np.loadtxt(path, skiprows=3, usecols=(0, 1, 2), unpack=True)

        with open(path, 'r') as sf:
            sf.readline()
            header = sf.readline()
            self.time, self.temp, self.dens = np.array(header.split()).astype(float)

        self.nuclei = np.array([nucleus(Z=Z, N=N, Y=Y) for N, Z, Y in zip(Ns, Zs, Ys)])
        self.path = path
        self.num = int(path.split("_")[-1][:4])
