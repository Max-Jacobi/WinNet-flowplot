class Flow(object):
    '''
    Input:
       nuc_in  - (nucleus object) target nucleus
       nuc_out - (nucleus object) product nucleus
    Attributes:
       nuc_in, nuc_out, dZ, dN, Z0, N0, flow
    '''

    def __init__(self, nuc_in, nuc_out, flow):
        self.Z0 = nuc_in.Z
        self.N0 = nuc_in.N
        self.dZ = nuc_out.Z - nuc_in.Z
        self.dN = nuc_out.N - nuc_in.N
        self.flow = flow
        self.nuc_in = nuc_in
        self.nuc_out = nuc_out

    def __repr__(self):
        return "flow: {}->{}".format(self.nuc_in.Name, self.nuc_out.Name)

    def __str__(self):
        return "{}->{}".format(self.nuc_in.Name, self.nuc_out.Name)
