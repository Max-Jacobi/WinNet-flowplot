from matplotlib.collections import PatchCollection
from matplotlib.patches import FancyArrow, Rectangle
from . import inset_axes
from .. import LogNorm, cm, plt
from ... import np


class FlowCollectionPlot(object):
    '''
    arguments:
    - ax             - matplotlib.Axes() object
    - flowcollection - FlowCollection object

    kwargs:
    - flow_range:      2
    - normalized:      True
    - norm:            LogNorm(1e-flow_range, 1) (*MaxFlow if not normalized)
    - cmap:            'jet'
    - scale_arrows:    True
    - kwargs for FancyArrow
    '''

    def __init__(self, ax, flowcollection, frange=2, normalized=True, scale_arrows=True, **kwargs):
        self.ax = ax
        self.flowcollection = flowcollection
        self.normalized = normalized
        self.frange = frange
        self.scale_arrows = scale_arrows
        self.MaxFlow = flowcollection.getMaxFlow()
        self.MinFlow = flowcollection.getMinFlow()

        if frange == 'all':
            if normalized:
                self.norm = kwargs.pop('norm', LogNorm(self.MinFlow/self.MaxFlow, 1, clip=True))
            else:
                self.norm = kwargs.pop('norm', LogNorm(self.MinFlow, self.MaxFlow, clip=True))
        else:
            if normalized:
                self.norm = kwargs.pop('norm', LogNorm(10**(-frange), 1, clip=True))
            else:
                self.norm = kwargs.pop('norm', LogNorm(self.MaxFlow*10**(-frange), self.MaxFlow, clip=True))
        self.cmap = cm.get_cmap(kwargs.pop('cmap', 'jet'))

        self.ax.axis(self.flowcollection.getBounds())
        kwargs.setdefault('lw', .3)

        Arrows = []
        Acol = []
        if normalized:
            flows = np.vectorize(lambda f: f.flow / self.MaxFlow)(self.flowcollection.flows)
        else:
            flows = np.vectorize(lambda f: f.flow)(self.flowcollection.flows)

        for fl, flow in zip(flows, self.flowcollection.flows):
            if frange != 'all' and self.norm(fl) <= 0:
                continue

            if scale_arrows:
                width = .2*self.norm(fl) + 0.01
            else:
                width = .15
            st = {"ec": 'k',
                  "width": width,
                  "head_width": 2.5*width,
                  "head_length": .5}
            nkwargs = kwargs.copy()
            for kv in st.items():
                nkwargs.setdefault(*kv)

            ar = FancyArrow(flow.N0, flow.Z0, flow.dN, flow.dZ,
                            length_includes_head=True,
                            **nkwargs)
            Arrows.append(ar)
            Acol.append(fl)
        Arrows = np.array(Arrows)
        Acol = np.array(Acol)
        self.Patches = PatchCollection(Arrows, norm=self.norm, cmap=self.cmap, match_original=True)
        self.Patches.set_array(Acol)
        self.ax.add_collection(self.Patches, )

    def addColorBar(self, xshift=-.06, yshift=.01, loc='lower right', **kwargs):
        '''
        adds a colorbar as inset at location loc
        xshift and yshift are bbox offset to anchor
        kwargs are passed to colorbar
        returns colorbar object
        '''
        if loc.split(' ')[1] == 'left':
            xshift += .035
        if loc.split(' ')[0] == 'upper':
            yshift -= .06
        self.cax = inset_axes(self.ax, width="3%", height="50%", loc=loc,
                              bbox_transform=self.ax.transAxes,
                              bbox_to_anchor=(xshift, yshift, 1, 1))
        self.cbar = plt.colorbar(self.Patches, cax=self.cax, **kwargs)
        if self.normalized:
            self.cax.set_title(r'$\frac{\Delta Y}{\rm{max}(\Delta Y)}$')
        else:
            self.cax.set_title(r'$\Delta Y$')

        return self.cbar
