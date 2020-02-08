import os
from matplotlib.collections import PatchCollection
from matplotlib.patches import FancyArrow
from ... import np


def find_flows(path):
    '''
    Arguments:
    path    - path to flowdir
    returns (flowpaths, flowtimes, flowpaths)
    '''
    flowpaths = []
    flowtimes = []
    flowtemps = []
    for file in os.listdir(path):
        with open('{}/{}'.format(path, file), 'r') as ff:
            header = ff.readline()
            if 'dt' in header:
                header = ff.readline()
                time, dt, temp, dens = np.array(header.split()).astype(float)
            else:
                header = ff.readline()
                time, temp, dens = np.array(header.split()).astype(float)
        flowtimes.append(time)
        flowtemps.append(temp)
        flowpaths.append('{}/{}'.format(path, file)
        sort=np.argsort(flowtimes)[::-1]
        flowpaths=np.array(flowpaths)[sort]
        flowtimes=np.array(flowtimes)[sort]
        flowtemps=np.array(flowtemps)[sort]
    return flowpaths, flowtimes, flowpaths

def plot_flow(flowfile, **kwargs):



    def flow_to_arrow(flow):
        fl=flow.fl
        width=.2*self.FlowNorm(fl)

        ar=FancyArrow(flow.N0, fl ow.Z0, flow.dN, flow.dZ,
                      fc=self.FlowCmap(nor m(fl)),
                      ec='k',
                      width=width,
                      head_width=3*width,
                      head_length=2*width,
                      length_includes_head=True,
                     )
        return ar
