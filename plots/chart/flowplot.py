from matplotlib.collections import PatchCollection
from matplotlib.patches import FancyArrow
from .. import LogNorm, cm


def plotFlowCollection(ax, flowcollection, **kwargs):
    '''
    kwargs:
    - flow_range:       2
    - norm:            LogNorm(MaxFlow*1e-flowrange, MaxFlow)
    - cmap:            'jet'
    - fixed_arrows:     False
    - kwargs for FancyArrow

    returns PatchCollection of arrows

    '''

    MaxFlow = flowcollection.getMaxFlow()
    frange = kwargs.pop('flow_range', 2)
    norm = kwargs.pop('norm', LogNorm(MaxFlow*10**(-frange), MaxFlow, clip=True))
    fixed = kwargs.pop('fixed_arrows', False)
    cmap = cm.get_cmap(kwargs.pop('cmap', 'jet'))

    ax.axis(flowcollection.getBounds())

    Arrows = []
    for fl in flowcollection.flows:
        if norm(fl.flow) <= 0:
            continue

        if fixed:
            width = .15
        else:
            width = .2*norm(fl.flow)
        st = {"fc": cmap(norm(fl.flow)),
              "ec": 'k',
              "width": width,
              "head_width": 3*width,
              "head_length": 2*width}
        nkwargs = kwargs.copy()
        for kv in st.items():
            nkwargs.setdefault(*kv)

        ar = FancyArrow(fl.N0, fl.Z0, fl.dN, fl.dZ,
                        length_includes_head=True,
                        **nkwargs)
        Arrows.append(ar)
    Patches = PatchCollection(Arrows, match_original=True)
    col = ax.add_collection(Patches, )
    return col
