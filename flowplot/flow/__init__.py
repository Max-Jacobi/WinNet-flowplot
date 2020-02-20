from .flowfile import FlowFile


def IntegrateFlows(paths, quite=False):
    '''
    Integrates all flowfiles in paths.
    Can take some time for large amounts of flows!
    Can take up a few 100mb of memory if alot of different flows are loaded!
    '''
    int_flows = FlowFile(paths[0])
    for ii, path in enumerate(paths):
        int_flows += FlowFile(path)
        if not quite:
            print("Progress: {:4.1f}%".format((ii+1)/len(paths)*100), end='\r')
    if not quite:
        print(F"Done: {len(int_flows.flows)} Flows on {len(int_flows.isotopes)} Nuclei")
    return int_flows
