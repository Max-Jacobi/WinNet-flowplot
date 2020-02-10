import os
from matplotlib.collections import PatchCollection
from matplotlib.patches import FancyArrow
from .. import LogNorm, Normalize, cm
from ... import np


def plotFlowCollection(ax, flowcollection, **kwargs):
    '''
    kwargs:
    - flowrange
    - flownorm
    - flowcmap
    - ynorm
    - ycmap
    - ploty
    - plotflows
    '''

    MaxFlow = flowcollection.getMaxFlow()
    def fill_kwargs(kwargs):
        standart =  {'flowrange': 1,
                     'ynorm': LogNorm(1e-10, 1e-3, clip=True),
                     'flowcmap': cm.jet,
                     'ycmap': cm.jet,
                     'ploty': True,
                     'plotflows': True}

        for k,s in standart.items():
            if k not in kwargs.keys():
                kwargs[k] = s
                
        if 'flownorm' not in kwargs.keys():
            kwargs['flownorm'] = LogNorm(MaxFlow*10**(-kwargs['flowrange']), MaxFlow, clip=True)
            
        return kwargs

    kwargs = fill_kwargs(kwargs)
    FlowNorm = kwargs['flownorm']
    YNorm = kwargs['ynorm']
    FlowCmap = kwargs['flowcmap']
    YCmap = kwargs['ycmap']

    def FlowtoArrow(flow):
        width=.2*FlowNorm(flow.flow)
        ar = FancyArrow(flow.N0, flow.Z0, flow.dN, flow.dZ,
                        fc=FlowCmap(FlowNorm(flow.flow)),
                        ec='k',
                        width=width,
                        head_width=3*width,
                        head_length=2*width,
                        length_includes_head=True)
        return ar

    Xmin, Xmax, Ymin, Ymax = flowcollection.getBounds()

    if kwargs['ploty']:
        Yim = np.empty((Xmax+1, Ymax+1))
        Yim.fill(np.nan)
        for iso in flowcollection.isotopes:
            if iso.Y >= flowcollection.ymin:
                Yim[iso.N, iso.Z] = iso.Y
            ax.imshow(np.log10(Yim).T, origin='lower', cmap='jet', norm=Normalize(-10,-3))

    if kwargs['plotflows']:
        Arrows = []
        for fl in flowcollection.flows:
            if FlowNorm(fl.flow) <= 0:
                continue
            Arrows.append(FlowtoArrow(fl))
            Patches = PatchCollection(Arrows, match_original=True)
            ax.add_collection(Patches)
