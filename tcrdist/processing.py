import pandas as pd
import numpy as np

__all__ = ['processNT',
           'readNT']

def processNT(chain, nuc, quals):
    """Process one nucleotide TCR sequence (any chain)

    Parameters
    ----------
    chain : str
        Values: alpha, beta, gamma, delta
    nuc : str
    quals : str
        Dot-separated list of integer scores

    Returns
    -------
    TCRChain object"""

    quals = np.array(quals.split('.'))




