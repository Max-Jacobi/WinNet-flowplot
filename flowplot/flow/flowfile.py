import os
from .. import np
from .flowcollection import FlowCollection


class FlowFile(FlowCollection):
    '''
    Get flows and isotopes from flowfile and create a flowcollection

    Input:
       path        - path to FlowFile
       ymin        - (optional) minimum abundance to be considered
       flmin       - (optional) minimum flow to be considered (relative to bigest flow in file)
    Attributes:
    - isotopes     - array of isotope objects
    - flows        - array of absolute(!) flows (dY/dt)*dt
    - num          - number of flowfile
    '''

    def __init__(self, path, ymin=1e-20, flmin=1e-10):
        super(FlowFile, self).__init__()

        self.path = path
        self.num = int(path.split("_")[-1][:-4])

        with open(path, 'r') as ff:
            header = ff.readline()
            if 'dt' in header:
                header = ff.readline()
                self.time, self.dt, self.temp, self.dens = np.array(header.split()).astype(float)
            else:
                header = ff.readline()
                self.time, self.temp, self.dens = np.array(header.split()).astype(float)
                self.dt = self.get_fake_dt()

        nins, zins, yins, nouts, zouts, youts, fls = np.loadtxt(path, skiprows=3, unpack=True)
        mfl = max(fls) * flmin
        for nin, zin, yin, nout, zout, yout, fl in zip(nins, zins, yins, nouts, zouts, youts, fls):
            if fl < mfl:
                continue
            if yin < ymin:
                continue
            if yout < ymin:
                yout = -np.inf
            self.addFlowFromZN(nin, zin, yin, nout, zout, yout, fl)

        self.sort()

        if len(self.flows) == 0:
            raise RuntimeError("No flows in {}".format(path))

    def get_fake_dt(self):
        if self.num == 0:
            return 0
        base = '/'.join(self.path.split('/')[:-2])
        for file in os.listdir(base):
            if '.par' in file:
                parfile = '{}/{}'.format(base, file)
                with open(parfile) as pf:
                    for line in pf:
                        if 'snapshot_every' in line:
                            out_every = int(line.split('=')[-1].strip())
                            break
                break
            else:
                out_every = 1

        prev_path = '_'.join(self.path.split('_')[:-1]) + '_{:04d}.dat'.format(self.num-1)
        if not os.path.isfile(prev_path):
            return 1

        with open(prev_path, 'r') as ff:
            ff.readline()
            header = ff.readline()
            prev_time, _, _ = np.array(header.split()).astype(float)
        return (self.time - prev_time)/out_every

    def __repr__(self):
        return "FlowFile at {}: {} flows".format(self.path, len(self.flows))
