import pandas as pd
import numpy as np
import itertools

from . import util
from .objects import DistanceParams
from . import tcr_distances
from .processing import getTCRID

__all__ = ['computeVRegionDistances',
           'basicDistance']

def computeVRegionDistances(organism, params=None):
    """Compute distances between all pairs of V-region genes for the
    specified organism. The resulting distance matrix (pd.DataFrame)
    can be (re)used for calculating the total distance between two clones.

    Params
    ------
    organism : str
        Values 'human' or 'mouse' are supported
    params : dict or td.objects.DistanceParams object
        Default parameters loaded by instantiating
        td.objects.DistanceParams or use None.
        Required keys are listed below.

    Returns
    -------
    pwDf : pd.DataFrame [n_regions x n_regions]
        DataFrame with identical index and columns, each containing
        labels of the V-gene regions

    
    Req. param keys
    ---------------
    weight_v_region
    gap_penalty_v_region
    distance_matrix"""

    if params is None:
        params = DistanceParams()

    """This returns a nested dict so index names and column names come from that"""
    d = tcr_distances.compute_all_v_region_distances(organism, params)
    pwDf = pd.DataFrame(d)
    return pwDf

def basicDistance(chains, tcr1, tcr2, VRegionDists=None, params=None):
    """Compute distance between two TCR clones. Distance can be based on one or
    both chains.

    Params
    ------
    chains : str
        E.g. use 'A' or 'B' or 'AB' for alpha-beta TCR analysis
    tcr1, tcr2 : DotDict or pd.Series
        Representation of TCR with dot notation enabled (req. set and get attribute)
        Can come from row of psDf or clonesDf or td.objects.TCRClone
        Requires the following attributes (some for each chain specified):
        organism, cdr3a, va_reps
    VRegionDists : pd.DataFrame
        Loaded via td.distances.computeVRegionDistances() or can be
        provided to avoid computation.
    params : dict or td.objects.DistanceParams object
        Default parameters loaded by instantiating td.objects.DistanceParams or use None.
        Required keys: distance_matrix, gap_penalty_cdr3_region, scale_factor

    Returns
    -------
    d : float
        Distance, possibly scaled by scale_factor"""

    if params is None:
        params = DistanceParams()
        
    assert tcr1.organism == tcr2.organism
    
    if VRegionDists is None:
        VRegionDists = tcr_distances.compute_all_v_region_distances(organism=tcr1.organism, params=params)
    
    dist=0.0
    if 'A' in chains:
        dist += min( (VRegionDists[x][y] for x in tcr1.va_reps.split(';') for y in tcr2.va_reps.split(';')))
        dist += tcr_distances.weighted_cdr3_distance(tcr1.cdr3a, tcr2.cdr3a, params)
    if 'B' in chains:
        dist += min((VRegionDists[x][y] for x in tcr1.vb_reps.split(';') for y in tcr2.vb_reps.split(';')))
        dist += tcr_distances.weighted_cdr3_distance(tcr1.cdr3b, tcr2.cdr3b, params)

    return params.scale_factor * dist

def computeBasicPWDistances(chains, clonesDf, VRegionDists=None, params=None):
    """Compute distance between all pairs of TCR clones in clonesDf.
    Uses td.distances.basicDistance() and organizes results in symmetric matrix.
    Can be called with clonesDf or psDf from td.processing.identifyClones() or
    td.processing.readPairedSequences()

    Reccomendation to use TCRID or CLONEID as index so that pwDf has identifiable
    columns and index.

    Example
    -------
    psDf = td.datasets.loadPSData('test_human_pairseqs')
    psDf = psDf.set_index('TCRID')
    pwDf = computeBasicPWDistances('AB', psDf):

    Params
    ------
    chains : str
        E.g. use 'A' or 'B' or 'AB' for alpha-beta TCR analysis
    clonesDf : pd.DataFrame
        Requires the following columns (some for each chain specified):
        organism, cdr3a, va_reps
    VRegionDists : pd.DataFrame
        Loaded via td.distances.computeVRegionDistances() or can be
        provided to avoid computation.
    params : dict or td.objects.DistanceParams object
        Default parameters loaded by instantiating td.objects.DistanceParams or use None.
        Required keys: distance_matrix, gap_penalty_cdr3_region, scale_factor

    Returns
    -------
    pwDf : pd.DataFrame
        Square, symmetric distance matrix with index of clonesDf as index and columns"""

    if params is None:
        params = DistanceParams()
        
    if VRegionDists is None:
        VRegionDists = tcr_distances.compute_all_v_region_distances(organism=clonesDf.organism.iloc[0],
                                                                    params=params)
    
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
    
    
def nearestNeighborDistance(chains, tcr, referenceTCRs, neighborPctile=10, neighborTopN=None, pwDf=None, weighted=False, VRegionDists=None, params=None):
    """Computes the nearest neighbor distance for a single TCR given
    a set of reference TCRs. The reference TCRs are all typically from a known group,
    typically an epitope-specific sort. The distance can be used for classification
    of new clones in the group.

    Must specify neighborPctile OR neighborTopN but not both.

    Initially looks for tcr and each reference TCR in pwDf using the TCR's identifier.
    Currently the identifier is the index from psDf or clonesDf that is used to compute pwDf.
    In the future this could be a hash that uniquely identifies a TCR clone.

    Params
    ------
    chains : str
        E.g. use 'A' or 'B' or 'AB' for alpha-beta TCR analysis
    tcr1 : DotDict or pd.Series
        Representation of TCR with dot notation enabled (req. set and get attribute)
        Can come from row of psDf or clonesDf or td.objects.TCRClone
        Requires the following attributes (some for each chain specified):
        organism, cdr3a, va_reps and an identifier which may be in pwDf
    referenceTCRs : pd.DataFrame
        Requires an index which matches pwDf and the following columns
        (some for each chain specified): organism, cdr3a, va_reps
    neighborPctile : float [0, 100]
        Percentile of nearest neighbors to include from referenceTCRs
    neighborTopN : int
        Number of nearest neighbors to include from referenceTCRs
    pwDf : pd.DataFrame
        Square, symmetric distance matrix which may contains the TCR or reference TCRs to
        save computation of distances.
    weighted : bool
        Compute a weighted nearest neighbor distance.
    VRegionDists : pd.DataFrame
        Loaded via td.distances.computeVRegionDistances() or can be
        provided to avoid computation.
    params : dict or td.objects.DistanceParams object
        Default parameters loaded by instantiating td.objects.DistanceParams or use None.
        Required keys: distance_matrix, gap_penalty_cdr3_region, scale_factor

    Returns
    -------
    d : float
        Mean of the distances to the nearest neighbors."""

    if params is None:
        params = DistanceParams()

    if VRegionDists is None:
        VRegionDists = tcr_distances.compute_all_v_region_distances(organism=tcr.organism,
                                                                    params=params)

    if (neighborPctile is None and neighborTopN is None) or (neighborPctile is not None and neighborTopN is not None):
        print('Must specify neighborPctile or neighborTopN, but not both.')

    if neighborTopN is None:
        neighborTopN = np.round(len(referenceTCRs) * (neighborPctile/100.))
    neighborTopN = int(neighborTopN)

    """Minimum is 1 neighbor and max is all neighbors"""
    neighborTopN = min(max(1, neighborTopN), len(neighborTopN))

    tcrID = getTCRID(tcr, chains)
   
    if pwDf is None or not tcrID in pwDf or not np.all([r in pwDf for r in referenceTCRs.index]):
        """Need to compute some or all distances first"""
        allD = np.zeros(len(referenceTCRs))
        for i,refID in enumerate(referenceTCRs.index):
            if pwDf is None or tcrID not in pwDf or refID not in pwDf:
                allD[i] = basicDistance(chains, tcr, referenceTCRs.loc[refID], VRegionDists=VRegionDists, params=params)
            else:
                allD[i] = pwDf.loc[tcrID, refID]
    else:
        allD = pwDf.loc[referenceTCRs.index, tcrID].values

    """Sort in place"""
    allD.sort()

    if not weighted:
        d = allD[:neighborTopN].mean()
    else:
        weights = 1. - np.arange(neighborTopN)/neighborTopN
        d = np.sum(allD[:neighborTopN] * weights) / weights.sum()

    return d

        
        




    
    