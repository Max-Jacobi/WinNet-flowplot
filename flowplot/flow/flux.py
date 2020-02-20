class Flow(object):
    '''
    Input:
       nuc_in  - (nucleus object) target nucleus
       nuc_out - (nucleus object) product nucleus
       flow    - value of flow
    Attributes:
       nuc_in, nuc_out, dZ, dN, Z0, N0, flow
    '''

    def __init__(self, iso_in, iso_out, flow):
        self.Z0 = iso_in.Z
        self.N0 = iso_in.N
        self.dZ = iso_out.Z - iso_in.Z
        self.dN = iso_out.N - iso_in.N
        self.flow = flow
        self.iso_in = iso_in
        self.iso_out = iso_out

    def __repr__(self):
        return "flow: {}->{}".format(self.iso_in.Name, self.iso_out.Name)

    def __str__(self):
        return "{}->{}".format(self.iso_in.name, self.iso_out.name)
