import pandas as pd
import numpy as np
import itertools

from . import util
from .objects import DistanceParams
from . import tcr_distances

__all__ = ['computeVRegionDistances',
           'basicDistance']

def computeVRegionDistances(organism, params=None):
    if params is None:
        params = DistanceParams()

    d = tcr_distances.compute_all_v_region_distances(organism, params)
    return pd.DataFrame(d)

def basicDistance(chains, tcr1, tcr2, VRegionDists=None, params=None):
    if params is None:
        params = DistanceParams()
        
    assert tcr1.organism == tcr2.organism
    
    if VRegionDists is None:
        vrDists = tcr_distances.compute_all_v_region_distances(organism=tcr1.organism, params=params)
    
    dist=0.0
    if 'A' in chains:
        dist += min( (vrDists[x][y] for x in tcr1.va_reps for y in tcr2.va_reps))
        dist += tcr_distances.weighted_cdr3_distance(tcr1.cdr3a, tcr2.cdr3a, params)
    if 'B' in chains:
        dist += min((vrDists[x][y] for x in tcr1.vb_reps for y in tcr2.vb_reps))
        dist += tcr_distances.weighted_cdr3_distance(tcr1.cdr3b, tcr2.cdr3b, params)

    return params.scale_factor * dist

def computePairwiseDistances(chains, clonesDf, VRegionDists=None, params=None):
    if params is None:
        params = DistanceParams()
        
    if VRegionDists is None:
        VRegionDists = tcr_distances.compute_all_v_region_distances(organism=tcr1.organism, params=params)
    
    nTCRs = clonesDf.shape[0]
    pwDist = np.zeros((nTCRs, nTCRs))
    for i,j in itertools.product(range(nTCRs), range(nTCRs)):
        if i <= j:
            d = basicDistance(chains,
                              clonesDf.iloc[i],
                              clonesDf.iloc[j], 
                              VRegionDists=VRegionDists, params=params)
            pwDist[i,j] = d
            pwDist[j,i] = d
    pwDf = pd.DataFrame(pwDist, index=clonesDf.index, columns=clonesDf.index)
    return pwDf
    
    
def nearestNeighborDistance(chains, tcr1, referenceTCRs, pwDf, params=None):
    if params is None:
        params = DistanceParams()
    pass
    