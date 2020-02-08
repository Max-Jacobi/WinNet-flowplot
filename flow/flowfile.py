import os
from .. import np
from ..nucleus import Nucleus
from .flux import Flow
from .flowcollection import FlowCollection


class FlowFile(FlowCollection):
    '''
    Get flows and nuclei from flowfile.

    Input:
       path        - path to FlowFile
       num         - number of flowfile
       ymin        - (optional) minimum abundance to be considered
    Attributes:
    - nuclei       - array of nucleus objects
    - flows        - array of absolute(!) flows (dY/dt)*dt
    Methods:
    - addNucleus
    - getNucleus
    - getMaxFlow
    - addition     - creates a new FlowFile instance
                     - flows are added together
                     - nucleus with higher abundance is kept
    '''

    def __init__(self, path, ymin=1e-10):
        super(FlowFile, self).__init__()

        self.path = path
        self.num = int(path.split("_")[-1][:-4])

        # nin, zin, yin, nout, zout, yout, fls = np.loadtxt(path, skiprows=3, unpack=True, )
        with open(path, 'r') as ff:
            header = ff.readline()
            if 'dt' in header:
                header = ff.readline()
                self.time, self.dt, self.temp, self.dens = np.array(header.split()).astype(float)
            else:
                header = ff.readline()
                self.time, self.temp, self.dens = np.array(header.split()).astype(float)
                self.dt = self._get_fake_dt()

        for ni, zi, yi, no, zo, yo, fl in np.loadtxt(path, skiprows=3):
            if (yi < ymin) or (fl <= 1e-99):
                continue
            if yo < ymin:
                yo = -np.inf
            self.addFlowFromZN(nin, zin, yin, nout, zout, yout, fl)

        self.sort()

        if len(self.flows) == 0:
            raise RuntimeError("No flows in {}".format(path))

    def _get_fake_dt(self):
        if self.num == 0:
            return 0
        base = '/'.join(path.split('/')[:-2])
        for file in os.listdir(base):
            if '.par' in file:
                parfile = '{}/{}'.format(base,file)
                with open(parfile) as pf:
                    for line in pf:
                        if 'snapshot_every' in line:
                            out_every = int(line.split('=')[-1].strip())
                            break
                break
            else:
                out_every = 1

        prev_path = '_'.join(self.path.split('_')[:-1]) + '_{:04d}.dat'.format(self.num-1)
        with open(prev_path, 'r') as ff:
            ff.readline()
            header = ff.readline()
            prev_time, _, _ = np.array(header.split()).astype(float)
        return (self.time - prev_time)/out_every

    def __repr__(self):
        return "FlowFile at {}: {} flows".format(self.path, len(self.flows))
