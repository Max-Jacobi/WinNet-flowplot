from . import Elnames, np


class Isotope(object):

    """
    isotope: contains basic information of an isotope.
    For example the name, the name in the network, the amount of protons, neutrons and the mass number
    """

    def __init__(self, name=None, Z=None, N=None, Y=np.nan, chk=None):
        """
        Input either:
          name       - name of the isotope
        or:
          Z          - proton number of isotope
          N          - neutron number of isotope
        Atributes:
          A, Z, N, name, Name, el, El, ab
        """

        if chk is not None:
            self.Z = int(chk//1e3)
            self.N = int(chk - self.Z*1e3)
            self.A = self.Z + self.N
            self._get_Name()
        elif name is not None:
            self.name = name.lower()
            self.Name = name[0].upper() + name[1:].lower()
            self._get_ZN()
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

    def _get_ZN(self):
        A = ''.join(filter(lambda x: x.isdigit(), self.name))
        try:
            self.A = int(A)
        except ValueError:
            raise ValueError("Can't get A from name {}".format(self.name))

        self.el = ''.join(filter(lambda x: x.isalpha(), self.name))
        if len(self.el) not in [1, 2]:
            raise ValueError("Cant get element name from name {}".format(self.name))

        self.el = self.el.lower()
        self.El = self.el[0].upper() + self.el[1:]

        Z = np.argwhere(Elnames == self.el)
        try:
            self.Z = int(Z)
        except TypeError:
            print(self.el)
            print(np.argwhere(Elnames == self.el))
        except ValueError:
            raise ValueError("Can't get Z from element name {}".format(self.name))

        self.N = self.A - self.Z

    def __repr__(self):
        repr = "Isotope: {}".format(self.Name)
        if not np.isnan(self.Y):
            repr += ": Abundance: {}".format(self.Y)
        return repr

    def __str__(self):
        return self.name
