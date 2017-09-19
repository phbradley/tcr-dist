class DotDict(dict):
    def __getattr__(self, item):
        if item in self:
            return self[item]
        raise AttributeError

    def __setattr__(self, key, value):
        if key in self:
            self[key] = value
            return
        raise AttributeError

class TCRClone(DotDict):
    """Object that contains all info for a single TCR clone.
    As a DotDict attributes can be accessed using dot notation
    or standard dict key access."""

    alphaAA = ''
    betaAA = ''
    gammaAA = ''
    deltaAA = ''
    
    subjid = ''
    cloneid = ''
    epitope = ''

    def __init__(self, chain1, chain2, **kwargs):
        for k in chain1.keys():
            self[k] = chain1[k]
        for k in chain2.keys():
            self[k] = chain2[k]

        for k in kwargs.keys():
            self[k] = kwargs[k]

class TCRChain(DotDict):
    def __init__(self, **kwargs):
        for k in kwargs.keys():
            self[k] = kwargs[k]    