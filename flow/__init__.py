from .flowfile import FlowFile
from .snapshot import Snapshot


def IntegrateFlows(paths, quite=False):
    '''
    Integrates all flowfiles with paths in paths.
    Can take some time for large amounts of flows!
    Can take up a few 100mb of memory if alot of different flows are loaded!
    '''
    int_flows = FlowFile(paths[0])
    for ii, path in enumerate(paths):
        int_flows += FlowFile(path)
        if not quite:
            print(f"Progress: {(ii+1)/len(paths)*100:5.1f}%", end='\r')
    if not quite:
        print(F"Done: {len(int_flows.flows)} Flows on {len(int_flows.isotopes)} Nuclei")
    return int_flows


def
