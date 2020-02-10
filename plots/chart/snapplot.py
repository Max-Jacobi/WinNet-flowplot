from ... import np


def plotIsotopeCollection(ax, isotopecollection, **kwargs):
    '''
    kwargs: same as imshow

    returns PatchCollection of arrows

    '''

    Xmin, Xmax, Ymin, Ymax = isotopecollection.getBounds()

    Yim = np.empty((Xmax+1, Ymax+1))
    Yim.fill(np.nan)
    for iso in isotopecollection.isotopes:
        if iso.Y >= isotopecollection.ymin:
            Yim[iso.N, iso.Z] = iso.Y
    im = ax.imshow(np.log10(Yim).T, origin='lower', **kwargs)
    return im
