from .. import np
from .flux import Flow


class FlowCollection(object):
    '''
    contains nuclei and their flows
    Attributes:
    - nuclei
    - flows
    Methods:
    - addNucleus
    - getNucleus
    - getMaxFlow
    - addition     - creates a new FlowCollection instance
                         - flows are added together
                         - keep nuclei with higher abundances
    '''

    def __init__(self,  ymin=1e-10):
        self.ymin = ymin
        self.flows = np.array([])
        self.nuclei = np.array([])
        self.known_nucs = np.array([])


    def addFlowFromName(self, name, flow):
        '''
        for addition of flowcollection
        add Flow from name 
        assumes Nuclei are in self.nuclei
        '''
        name_in, name_out = name.split('->')
        nuc_in = self.getNuccleus(name=name_in)
        nuc_out = self.getNuccleus(name=name_out)
        flow = Flow(nuc_in, nuc_out, flow)
        self.flows = np.append(self.flows, flow)
        nuc_in.flow_out = np.append(nuc_in.flow_out, flow)
        nuc_out.flow_in = np.append(nuc_out.flow_in, flow)
        return flow

    def addFlowFromZN(self, Nin,  Zin, Yin, Nout, Zout, Yout, flow):
        ''' for flowfile readin
        add Flow from Zin, Zout, Nin and Nout
        '''
        chk_in = 1e3*Zin + Nin
        chk_out = 1e3*Zout + Nout
        if chk_in in self.known_nucs:
            nuc_in = self.getNuccleus(chk=chk_in)
        else:
            nuc_in = self.addNucleus(chk=chk_in, Y=Y_in)
             
        if chk_out in self.known_nucs:
            nuc_out = self.getNuccleus(chk=chk_out)
        else:
            nuc_out = self.addNucleus(chk=chk_out, Y=Y_out)
        flow = Flow(nuc_in, nuc_out, flow)
        self.flows = np.append(self.flows, flow)
        nuc_in.flow_out = np.append(nuc_in.flow_out, flow)
        nuc_out.flow_in = np.append(nuc_out.flow_in, flow)
        return flow

    def __add__(self, other):
        new = FlowCollection(ymin=min(self.ymin, other.ymin))

        nucs1 = self.nuclei.astype(str)
        abs1 = map(lambda n: n.Y, self.nuclei)
        flow_names1 = self.flows.astype(str)
        flows1 = map(lambda f: f.flow, self.flows)

        nucs2 = other.nuclei.astype(str)
        abs2 = map(lambda n: n.Y, other.nuclei)
        flow_names2 = other.flows.astype(str)
        flows2 = map(lambda f: f.flow, other.flows)

        for name, i1, i2 in zip(*np.intersect1d(nucs1, nucs2, assume_unique=True, return_indices=True)):
            new.addNucleus(name=name, Y=max(abs1[i1], abs2[i2]))
        for name in np.setxor1d(nucs1, nucs2, assume_unique=True):
            Y = np.concatenate((abs1[nucs1 == name], abs2[nucs2 == name]))[0]        
            new.addNucleus(name=name, Y=Y)

        for name, i1, i2 in zip(*np.intersect1d(flow_names1, flow_names2, assume_unique=True, return_indices=True)):
            new.addFlowFromName(name, flows1[i1] + flows2[i2])
        for name in np.setxor1d(flow_names1, flow_names2, assume_unique=True):
            flow = np.concatenate((flows1[flow_names1 == name], flows2[flow_names2 == name]))[0]        
            new.addFlowFromName(name, flow)

        new.sort()
            
        return new

    def addNucleus(self, *kargs):
        '''
        Add a nucleus object from either:
        - chk 
        - name 
        - Z and N
        Optional:
        - Y
        '''
        nuc = Nucleus(**kwargs)
        self.known_nucs = np.append(self.known_nucs, chk)
        self.nuclei = np.append(self.nuclei, nuc)
        nuc.flow_in = np.array([])
        nuc.flow_out  = np.array([])
        return nuc

    def getNuccleus(self, chk=None, name=None, Z=None, N=None):
        '''
        Get a nucleus object from either:
        - chk 
        - name 
        - Z and N
        '''
        if chk is not None:
            ind = int(np.argwhere(self.known_nucs == chk))
        elif name is not None:
            ind = int(np.argwhere(self.nuclei.astype(int) == name))
        elif Z is not none and N is not None:
            chk = 1e3*Z + N
            ind = int(np.argwhere(self.known_nucs == chk))
        else:
            raise ValueError("Give name or checksum")
        return self.nuclei[ind]

    def sort(self, nuclei=True, flows=True):
        if nuclei:
            self.nuclei = self.nuclei[np.argsort(self.nuclei.astype(str))]
        if flows:
            self.flows = self.flows[np.argsort(self.flows.astype(str))]

    def getMaxFlow(self):
        return max(list(map(lambda f: f.flow, self.flows)))

    def getBounds(self):
        ''' returns minZ, minN, maxZ, maxN'''
        maxZ = max(list(map(lambda n: n.Z, self.nuclei)))
        maxN = max(list(map(lambda n: n.N, self.nuclei)))
        minZ = min(list(map(lambda n: n.Z, self.nuclei)))
        minN = min(list(map(lambda n: n.N, self.nuclei)))
        return minZ, minN, maxZ, maxN

    def getFlowDict(self):
        keys = self.flows.astype(str)
        flows = map(lambda f: f.flow, self.flows)
        return dict(zip(keys, flows))
        
    def getNucleiDict(self):
        keys = self.nuclei.astype(str)
        ab = map(lambda m: n.Y, self.nuclei)
        return dict(zip(keys, ab))

    def __repr__(self):
        return "flowcollection: {} flows".format(len(self.flows))

