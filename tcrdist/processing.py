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



def get_translation( seq, frame ):
    assert frame[0] in '+-'
    if frame[0] == '-': seq = logo_tools.reverse_complement( seq )
    offset = abs( int( frame ))-1
    assert offset in range(3)
    seq = seq[offset:].lower()
    naa = len(seq)/3
    protseq = ''
    codons = []
    for i in range(naa):
        codon = seq[3*i:3*i+3]
        codons.append( codon )
        if '#' in codon:
            protseq += '#'
        else:
            protseq += genetic_code.get( codon, 'X' )
    return protseq,codons
