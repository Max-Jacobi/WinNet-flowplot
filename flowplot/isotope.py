from . import Elnames, np, getNZ


class Isotope(object):

    """
    isotope: contains basic information of an isotope.
    For example the name, the name in the network, the amount of protons, neutrons and the mass number
    """

    def __init__(self, name=None, Z=None, N=None, Y=np.nan, chk=None):
        """
        Input either of:
          name       - name of the isotope
          Z          - proton number of isotope
          N          - neutron number of isotope
          chk        - Z*1e3+N
        Atributes:
          A, Z, N, Y,name, Name, el, El
        """

        if chk is not None:
            self.Z = int(chk//1e3)
            self.N = int(chk - self.Z*1e3)
            self.A = self.Z + self.N
            self._get_Name()
        elif name is not None:
            self.name = name.lower()
            self.Name = name[0].upper() + name[1:].lower()
            self.N, self.Z = getNZ(self.name)
        elif Z is not None and N is not None:
            self.Z = int(Z)
            self.N = int(N)
            self.A = int(Z + N)
            self._get_Name()
        else:
            raise(ValueError("Give either name or Z and N of isotope"))

        # Special cases
        if self.name == 'p':
            self.name = 'h1'
        if self.name == 'n' or self.name == 'neutrons':
            self.name = 'neutron'
        if self.name == 'd':
            self.name = 'h2'
        if self.name == 't':
            self.name = 'h3'

        self.Y = Y

    def _get_Name(self):
        self.el = Elnames[self.Z]
        self.El = self.el[0].upper() + self.el[1:]

        self.name = self.el + str(self.A)
        self.Name = self.El + str(self.A)

    def __repr__(self):
        repr = "Isotope: {}".format(self.Name)
        if not np.isnan(self.Y):
            repr += ": Abundance: {}".format(self.Y)
        return repr

    def __str__(self):
        return self.name
