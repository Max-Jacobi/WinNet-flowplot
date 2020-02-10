from .. import np
from ..isotope import Isotope


class IsotopeCollection(object):
    '''
    contains isotopes and their flows
    Attributes:
    - isotopes
    - flows
    Methods:
    - getMaxFlow
    - addFlowFromName
    - addFlowFromZN
    Addition creates a new FlowCollection instance
       - flows are added together
       - isotopes with higher abundance are kept
    '''

    def __init__(self,  ymin=1e-10):
        self.ymin = ymin
        self.isotopes = np.array([])
        self._known_isos = np.array([])


    def _addIsotope(self, **kwargs):
        '''
        Add a isotope object from either:
        - chk
        - name
        - Z and N
        Optional:
        - Y
        '''
        iso = Isotope(**kwargs)
        self._known_isos = np.append(self._known_isos, 1e3*iso.Z+iso.N)
        self.isotopes = np.append(self.isotopes, iso)
        return iso

    def _getIsotope(self, chk=None, name=None, Z=None, N=None):
        '''
        Get a isotope object from either:
        - chk
        - name
        - Z and N
        '''
        if name is not None:
            ind = int(np.argwhere(self.isotopes.astype(str) == name))
        elif chk is not None:
            ind = int(np.argwhere(self._known_isos == chk))
        elif Z is not none and N is not None:
            chk = 1e3*Z + N
            ind = int(np.argwhere(self._known_isos == chk))
        else:
            raise ValueError("Give name, Z and N, or checksum")
        return self.isotopes[ind]

    def sort(self):
        self.isotopes = self.isotopes[np.argsort(self.isotopes.astype(str))]

    def getMaxY(self):
        '''
        returns maximal flow inside flowfile
        '''
        return max(list(map(lambda n: n.Y, self.isotopes)))

    def getBounds(self):
        '''
        returns (minN, minZ, maxN, maxZ)
        '''
        maxZ = max(list(map(lambda n: n.Z, self.isotopes)))
        maxN = max(list(map(lambda n: n.N, self.isotopes)))
        minZ = min(list(map(lambda n: n.Z, self.isotopes)))
        minN = min(list(map(lambda n: n.N, self.isotopes)))
        return minN, maxN, minZ, maxZ

    def getIsotopesDict(self):
        keys = self.isotopes.astype(str)
        ab = map(lambda n: n.Y, self.isotopes)
        return dict(zip(keys, ab))

    def __add__(self, other):
        new = IsotopeCollection(ymin=min(self.ymin, other.ymin))

        isos1 = self.isotopes.astype(str)
        abs1 = np.array(list(map(lambda n: n.Y, self.isotopes)))

        isos2 = other.isotopes.astype(str)
        abs2 = np.array(list(map(lambda n: n.Y, other.isotopes)))

        for name, i1, i2 in zip(*np.intersect1d(isos1, isos2, assume_unique=True, return_indices=True)):
            new._addIsotope(name=name, Y=max(abs1[i1], abs2[i2]))
        for name in np.setxor1d(isos1, isos2, assume_unique=True):
            Y = np.concatenate((abs1[isos1 == name], abs2[isos2 == name]))[0]
            new._addIsotope(name=name, Y=Y)

        new.sort()
        return new

    def __repr__(self):
        return "IsotopeCollection: {} ".format(len(self.flows))
