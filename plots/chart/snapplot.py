from ... import np
from .. import plt, LogNorm
from . import inset_axes


class IsotopeCollectionPlot(object):
    '''
    kwargs: same as imshow
    - norm defaults to LogNorm(1e-10,10**(ceil(log10(MaxY))))
    - cmap defaults to 'jet'
    '''

    def __init__(self, ax, isotopecollection, grid=True, **kwargs):
        self.ax = ax
        self.isotopecollection = isotopecollection
        MaxY = isotopecollection.getMaxY()
        kwargs.setdefault('norm', LogNorm(1e-10, 10**(int((np.log10(MaxY)))), clip=True))
        kwargs.setdefault('cmap', 'jet')
        Xmin, Xmax, Ymin, Ymax = isotopecollection.getBounds()
        if grid:
            ax.vlines([x+.5 for x in range(Xmin-10, Xmax+10)], Ymin-10, Ymax+10, color='k', lw=.2)
            ax.hlines([y+.5 for y in range(Ymin-10, Ymax+10)], Xmin-10, Xmax+10, color='k', lw=.2)
        self.Yarray = np.empty((Xmax+1, Ymax+1))
        self.Yarray.fill(np.nan)
        for iso in isotopecollection.isotopes:
            if iso.Y >= isotopecollection.ymin:
                self.Yarray[iso.N, iso.Z] = iso.Y
        # self.Yarray = np.log10(self.Yarray)
        self.im = ax.imshow(self.Yarray.T, origin='lower', **kwargs)

    def addColorBar(self, xshift=-.06, yshift=0, loc='lower right', **kwargs):
        '''
        adds a colorbar as inset at location loc
        xshift and yshift are bbox to anchor offset
        kwargs are passed to colorbar
        returns colorbar object
        '''
        self.cax = inset_axes(self.ax, width="3%", height="50%", loc=loc,
                              bbox_transform=self.ax.transAxes,
                              bbox_to_anchor=(xshift, yshift, 1, 1))
        self.cbar = plt.colorbar(self.im, cax=self.cax, **kwargs)
        self.cax.set_title(r'$Y$')
        return self.cbar
