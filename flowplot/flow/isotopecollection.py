from .. import np, getName
from ..isotope import Isotope


class IsotopeCollection(object):
    '''
    contains isotopes and their flows
    Attributes:
    - isotopes
    - flows
    Methods:
    - getIsotope
    - sort
    - getMaxY
    - getBounds
    Addition creates a new FlowCollection instance
       - flows are added together
       - isotopes with higher abundance are kept
    '''

    def __init__(self,  ymin=1e-10):
        self.ymin = ymin
        self.isotopes = np.array([])

    def addIsotope(self, **kwargs):
        '''
        Add a isotope object from either:
        - name
        - Z and N
        Optional:
        - Y
        '''
        iso = Isotope(**kwargs)
        self.isotopes = np.append(self.isotopes, iso)
        return iso

    def getIsotope(self, name=None, Z=None, N=None):
        '''
        Get a isotope object from either:
        - name
        - Z and N
        '''
        if name is None:
            if Z is None or N is None:
                raise ValueError("Give name, Z and N")
            name = getName(N, Z)
        ind = int(np.argwhere(self.isotopes.astype(str) == name.lower()))

        return self.isotopes[ind]

    def sort(self):
        '''
        sort isotopes by abundance
        '''
        abus = np.vectorize(lambda i: i.Y)(self.isotopes)
        self.isotopes = self.isotopes[np.argsort(abus)]

    def getMaxY(self):
        '''
        returns maximal abundance in isotopecollection
        '''
        return max(list(map(lambda n: n.Y, self.isotopes)))

    def getBounds(self, for_plot=False):
        '''
        returns (minN, minZ, maxN, maxZ)
        if for_plot is True add and subtract .5 so that ax.axis(IsotopeCollection().getBounds()) works
        '''
        maxZ = max(list(map(lambda n: n.Z, self.isotopes)))
        maxN = max(list(map(lambda n: n.N, self.isotopes)))
        minZ = min(list(map(lambda n: n.Z, self.isotopes)))
        minN = min(list(map(lambda n: n.N, self.isotopes)))
        if for_plot:
            return minN-.5, maxN+.5, minZ-.5, maxZ+.5
        else:
            return minN, maxN, minZ, maxZ

    def __add__(self, other):
        new = IsotopeCollection(ymin=min(self.ymin, other.ymin))

        isos1 = self.isotopes.astype(str)
        abs1 = np.vectorize(lambda n: n.Y)(self.isotopes)

        isos2 = other.isotopes.astype(str)
        abs2 = np.vectorize(lambda n: n.Y)(other.isotopes)

        for name, i1, i2 in zip(*np.intersect1d(isos1, isos2, assume_unique=True, return_indices=True)):
            new.addIsotope(name=name, Y=max(abs1[i1], abs2[i2]))
        for name in np.setxor1d(isos1, isos2, assume_unique=True):
            Y = np.concatenate((abs1[isos1 == name], abs2[isos2 == name]))[0]
            new.addIsotope(name=name, Y=Y)

        new.sort()
        return new

    def __repr__(self):
        return "IsotopeCollection: {} ".format(len(self.flows))
