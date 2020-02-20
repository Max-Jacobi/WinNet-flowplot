from .. import np
from .flux import Flow
from .isotopecollection import IsotopeCollection


class FlowCollection(IsotopeCollection):
    '''
    contains isotopes and their flows
    Attributes:
    - isotopes
    - flows
    Methods:
    - sort
    - getMaxFlow
    - getMinFlow
    - getFlowsTo
    Addition creates a new FlowCollection instance
       - flows are added together
       - keeps isotope with higher abundance
    '''

    def __init__(self,  ymin=1e-10):
        super(FlowCollection, self).__init__()
        self.flows = np.array([])

    def _addFlowFromName(self, name, flow):
        '''
        assumes Isotopes are in self.isotopes
        '''
        name_in, name_out = name.split('->')
        iso_in = self.getIsotope(name=name_in)
        iso_out = self.getIsotope(name=name_out)
        flow = Flow(iso_in, iso_out, flow)
        self.flows = np.append(self.flows, flow)
        iso_in.flow_out = np.append(iso_in.flow_out, flow)
        iso_out.flow_in = np.append(iso_out.flow_in, flow)
        return flow

    def _addFlowFromZN(self, Nin,  Zin, Yin, Nout, Zout, Yout, flow):
        '''
        add Flow from Zin, Zout, Nin and Nout
        checks for isotope in isotopes and adds it if it is not present
        '''
        if flow < 1e-99:
            return
        chk_in = 1e3*Zin + Nin
        chk_out = 1e3*Zout + Nout

        if chk_in in self._known_isos:
            iso_in = self.getIsotope(chk=chk_in)
        else:
            iso_in = self._addIsotope(chk=chk_in, Y=Yin)

        if chk_out in self._known_isos:
            iso_out = self.getIsotope(chk=chk_out)
        else:
            iso_out = self._addIsotope(chk=chk_out, Y=Yout)

        flow = Flow(iso_in, iso_out, flow)
        self.flows = np.append(self.flows, flow)
        iso_in.flow_out = np.append(iso_in.flow_out, flow)
        iso_out.flow_in = np.append(iso_out.flow_in, flow)
        return flow

    def _addIsotope(self, **kwargs):
        '''
        Add a isotope object from either of:
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
        '''
        sort flows by flow
        sort isotopes by abundance
        '''
        super(FlowCollection, self).sort()
        flows = np.vectorize(lambda f: f.flow)(self.flows)
        self.flows = self.flows[np.argsort(flows)]

    def getMaxFlow(self):
        '''
        returns maximal flow in flowfile
        '''
        return np.max(np.vectorize(lambda f: f.flow)(self.flows))

    def getMinFlow(self):
        '''
        returns minimal flow in flowfile
        '''
        return np.min(np.vectorize(lambda f: f.flow)(self.flows))

    def getFlowsTo(self, name, N):
        '''
        Returns a new flowcollection with all isotopes but only the flows
        that are connected to the isotope with name <name> by N steps.
        Flows are scaled to match their parent flow
        '''
        def prevFlows(flow):
            new_flows = []
            for fl in flow.iso_in.flow_in:
                if fl.flow < flow.flow*5e-2:
                    continue
                new_flows.append(Flow(fl.iso_in, fl.iso_out, fl.flow))
            sum_flows = sum(map(lambda f: f.flow, new_flows))
            for fl in new_flows:
                fl.flow *= flow.flow/sum_flows
            return np.array(new_flows)

        iso = self.getIsotope(name=name)
        subcol = FlowCollection()
        subcol.isotopes = self.isotopes
        subcol._known_isos = np.vectorize(lambda i: i.Z*1e3 + i.N)(subcol.isotopes)
        for fl in iso.flow_in:
            subcol._addFlowFromZN(Nin=fl.iso_in.N, Zin=fl.iso_in.Z, Yin=fl.iso_in.Y,
                                  Nout=fl.iso_out.N, Zout=fl.iso_out.Z, Yout=fl.iso_out.Y,
                                  flow=fl.flow)
        oldFlows = subcol.flows
        for nn in range(N):
            newFlows = []
            for fl in np.concatenate([prevFlows(fl) for fl in oldFlows]):
                if not np.any(str(fl) == subcol.flows.astype(str)):
                    newFlows.append(fl)
                    subcol._addFlowFromZN(Nin=fl.iso_in.N, Zin=fl.iso_in.Z, Yin=fl.iso_in.Y,
                                          Nout=fl.iso_out.N, Zout=fl.iso_out.Z, Yout=fl.iso_out.Y,
                                          flow=fl.flow)
            oldFlows = newFlows
        subcol.sort()
        return subcol

    def getFlowsFrom(self, name, N):
        '''
        Returns a new flowcollection with all isotopes but only the flows
        that are connected to the isotope with name <name> by N steps.
        Flows are scaled to match their parent flow
        '''
        def followingFlows(flow):
            new_flows = []
            for fl in flow.iso_out.flow_out:
                if fl.flow < flow.flow*5e-2:
                    continue
                new_flows.append(Flow(fl.iso_in, fl.iso_out, fl.flow))
            sum_flows = sum(map(lambda f: f.flow, new_flows))
            for fl in new_flows:
                fl.flow *= flow.flow/sum_flows
            return np.array(new_flows)

        iso = self.getIsotope(name=name)
        subcol = FlowCollection()
        subcol.isotopes = self.isotopes
        subcol._known_isos = np.vectorize(lambda i: i.Z*1e3 + i.N)(subcol.isotopes)
        for fl in iso.flow_out:
            subcol._addFlowFromZN(Nin=fl.iso_in.N, Zin=fl.iso_in.Z, Yin=fl.iso_in.Y,
                                  Nout=fl.iso_out.N, Zout=fl.iso_out.Z, Yout=fl.iso_out.Y,
                                  flow=fl.flow)
        oldFlows = subcol.flows
        for nn in range(N):
            newFlows = []
            for fl in np.concatenate([followingFlows(fl) for fl in oldFlows]):
                if not np.any(str(fl) == subcol.flows.astype(str)):
                    newFlows.append(fl)
                    subcol._addFlowFromZN(Nin=fl.iso_in.N, Zin=fl.iso_in.Z, Yin=fl.iso_in.Y,
                                          Nout=fl.iso_out.N, Zout=fl.iso_out.Z, Yout=fl.iso_out.Y,
                                          flow=fl.flow)
            oldFlows = newFlows
        subcol.sort()
        return subcol

    def __add__(self, other):
        new = FlowCollection(ymin=min(self.ymin, other.ymin))

        # get arrays of abundaces and flows as well as isotope and flow names
        isos1 = self.isotopes.astype(str)
        isos2 = other.isotopes.astype(str)
        abs1 = np.vectorize(lambda n: n.Y)(self.isotopes)
        abs2 = np.vectorize(lambda n: n.Y)(other.isotopes)
        flow_names1 = self.flows.astype(str)
        flow_names2 = other.flows.astype(str)
        flow_names1_reversed = np.vectorize(lambda f: '->'.join(f.split('->')[::-1]))(flow_names1)
        flow_names2_reversed = np.vectorize(lambda f: '->'.join(f.split('->')[::-1]))(flow_names2)
        flows1 = np.vectorize(lambda f: f.flow)(self.flows)
        flows2 = np.vectorize(lambda f: f.flow)(other.flows)

        # for isotopes in both sets add an isotope with max(Y1, Y2)
        for name, i1, i2 in zip(*np.intersect1d(isos1, isos2, assume_unique=True, return_indices=True)):
            new._addIsotope(name=name, Y=max(abs1[i1], abs2[i2]))
        # add the rest
        for name in np.setxor1d(isos1, isos2, assume_unique=True):
            Y = np.concatenate((abs1[isos1 == name], abs2[isos2 == name]))[0]
            new._addIsotope(name=name, Y=Y)

        # for flows in both sets add a flow with the sum of flows
        names, is1, is2 = np.intersect1d(flow_names1, flow_names2, assume_unique=True, return_indices=True)
        for name, flow in zip(names, flows1[is1]+flows2[is2]):
            new._addFlowFromName(name, flow)

        # if other contais the reverse of a flow in self but the flow in self is bigger
        # add the flow from self but substract the flow from other
        names, is1, is2 = np.intersect1d(flow_names1, flow_names2_reversed, assume_unique=True, return_indices=True)
        delflow = flows1[is1]-flows2[is2]
        mask = delflow > 0
        for name, flow in zip(names[mask], delflow[mask]):
            new._addFlowFromName(name, flow)

        # if other contais the reverse of a flow in self but the flow in other is bigger
        # add the flow from other but substract the flow from self
        mask = np.logical_not(mask)
        names = flow_names2[is2]
        for name, flow in zip(names[mask], -delflow[mask]):
            new._addFlowFromName(name, flow)

        # finaly for flows that are only self or other and not the reverse of anything just add the flow
        names = np.setxor1d(flow_names1, flow_names2, assume_unique=True)
        rev_names = np.concatenate([flow_names1_reversed, flow_names2_reversed])
        for name in names:
            if name not in rev_names:
                flow = np.concatenate((flows1[flow_names1 == name], flows2[flow_names2 == name]))[0]
                new._addFlowFromName(name, flow)

        new.sort()
        return new

    def __repr__(self):
        return "FlowCollection: {} flows".format(len(self.flows))
