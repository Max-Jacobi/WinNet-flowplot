from .. import np
from .flux import Flow
from ..isotope import Isotope
from .isotopecollection import IsotopeCollection


class FlowCollection(IsotopeCollection):
    '''
    contains isotopes and their flows
    Attributes:
    - isotopes
    - flows
    Methods:
    - getMaxFlow
    Addition creates a new FlowCollection instance
       - flows are added together
       - isotopes with higher abundance are kept
    '''

    def __init__(self,  ymin=1e-10):
        super(FlowCollection, self).__init__()
        self.flows = np.array([])

    def _addFlowFromName(self, name, flow):
        '''
        for addition of flowcollection
        add Flow from name
        assumes Isotopes are in self.isotopes
        '''
        name_in, name_out = name.split('->')
        iso_in = self._getIsotope(name=name_in)
        iso_out = self._getIsotope(name=name_out)
        flow = Flow(iso_in, iso_out, flow)
        self.flows = np.append(self.flows, flow)
        iso_in.flow_out = np.append(iso_in.flow_out, flow)
        iso_out.flow_in = np.append(iso_out.flow_in, flow)
        return flow

    def _addFlowFromZN(self, Nin,  Zin, Yin, Nout, Zout, Yout, flow):
        '''
        for flowfile readin
        add Flow from Zin, Zout, Nin and Nout
        '''
        chk_in = 1e3*Zin + Nin
        chk_out = 1e3*Zout + Nout

        if chk_in in self._known_isos:
            iso_in = self._getIsotope(chk=chk_in)
        else:
            iso_in = self._addIsotope(chk=chk_in, Y=Yin)

        if chk_out in self._known_isos:
            iso_out = self._getIsotope(chk=chk_out)
        else:
            iso_out = self._addIsotope(chk=chk_out, Y=Yout)

        flow = Flow(iso_in, iso_out, flow)
        self.flows = np.append(self.flows, flow)
        iso_in.flow_out = np.append(iso_in.flow_out, flow)
        iso_out.flow_in = np.append(iso_out.flow_in, flow)
        return flow

    def _addIsotope(self, **kwargs):
        '''
        Add a isotope object from either:
        - chk
        - name
        - Z and N
        Optional:
        - Y
        '''
        iso = super(FlowCollection, self)._addIsotope(**kwargs)
        iso.flow_in = np.array([])
        iso.flow_out = np.array([])
        return iso

    def sort(self):
        super(FlowCollection, self).sort()
        self.flows = self.flows[np.argsort(self.flows.astype(str))]

    def getMaxFlow(self):
        '''
        returns maximal flow inside flowfile
        '''
        return max(list(map(lambda f: f.flow, self.flows)))

    def getFlowDict(self):
        keys = self.flows.astype(str)
        flows = map(lambda f: f.flow, self.flows)
        return dict(zip(keys, flows))

    def __add__(self, other):
        new = FlowCollection(ymin=min(self.ymin, other.ymin))

        isos1 = self.isotopes.astype(str)
        abs1 = np.array(list(map(lambda n: n.Y, self.isotopes)))
        flow_names1 = self.flows.astype(str)
        flows1 = np.array(list(map(lambda f: f.flow, self.flows)))

        isos2 = other.isotopes.astype(str)
        abs2 = np.array(list(map(lambda n: n.Y, other.isotopes)))
        flow_names2 = other.flows.astype(str)
        flows2 = np.array(list(map(lambda f: f.flow, other.flows)))

        for name, i1, i2 in zip(*np.intersect1d(isos1, isos2, assume_unique=True, return_indices=True)):
            new._addIsotope(name=name, Y=max(abs1[i1], abs2[i2]))
        for name in np.setxor1d(isos1, isos2, assume_unique=True):
            Y = np.concatenate((abs1[isos1 == name], abs2[isos2 == name]))[0]
            new._addIsotope(name=name, Y=Y)

        for name, i1, i2 in zip(*np.intersect1d(flow_names1, flow_names2, assume_unique=True, return_indices=True)):
            new._addFlowFromName(name, flows1[i1] + flows2[i2])
        for name in np.setxor1d(flow_names1, flow_names2, assume_unique=True):
            flow = np.concatenate((flows1[flow_names1 == name], flows2[flow_names2 == name]))[0]
            new._addFlowFromName(name, flow)

        new.sort()

        return new

    def __repr__(self):
        return "FlowCollection: {} flows".format(len(self.flows))
