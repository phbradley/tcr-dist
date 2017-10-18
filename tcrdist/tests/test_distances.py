from __future__ import absolute_import, division, print_function
import os.path as op
import numpy as np
import pandas as pd
import itertools
import numpy.testing as npt
import pytest

import matplotlib
matplotlib.use('Agg')

import tcrdist as td

def exampleBetaChains():
    lines = """TRBV18*01 CASSPRHGISPLHF
    TRBV19*01,TRBV19*02 CASSPGGVTEAFF
    TRBV28*01 CASRGPGEPYEQYF
    TRBV4-1*01,TRBV4-1*02 CASSQAEGGLLSYEQYF
    TRBV6-6*01,TRBV6-6*02,TRBV6-6*03,TRBV6-6*04 CASMLGGPPQETQYF
    TRBV7-8*03 CASSHSWDSGTGELFF
    TRBV7-9*01,TRBV7-9*02,TRBV7-9*03 CASSSSGGAYNEQFF
    TRBV19*01,TRBV19*02 CASSIEGSGANVLTF
    TRBV7-8*01 CASSSSGGMNIQYF
    TRBV6-1*01 CASSEAPLLGGNEQYF
    TRBV29-1*01 CSVEARLF""".split('\n')

    tcrs = []
    for line in lines:
        l = line.split()
        tcr = td.objects.TCRChain()
        tcr.cdr3b = l[1]
        tcr.vb_genes = l[0].split(',')
        tcr.organism = 'human'
        tcr.vb_reps = ';'.join([td.util.get_rep(vb_gene, 'human') for vb_gene in tcr.vb_genes])
        tcrs.append(tcr)
    return tcrs
    
datasetsPath = op.join(td.__path__[0], 'datasets')

tempSkip = pytest.mark.skip(reason="Temporarily skipping for efficiency.")

#@tempSkip
def test_human_vregion():
    """Getting lots of NA distances..."""
    pwDf = td.distances.computeVRegionDistances('human')
    assert np.all(pwDf.index == pwDf.columns)

def test_mouse_vregion():
    pwDf = td.distances.computeVRegionDistances('mouse')
    assert np.all(pwDf.index == pwDf.columns)
    
def test_chain():
    tcrs = exampleBetaChains()
    pwmat = np.zeros((len(tcrs), len(tcrs)))
    for i,j in itertools.product(range(len(tcrs)), range(len(tcrs))):
        if j <= i:
            pwmat[i,j] = td.distances.basicDistance('B', tcrs[i], tcrs[j])
            pwmat[j,i] = pwmat[i,j]
    
def test_compute_paired_chain_distance():
    psDf = td.datasets.loadPSData('test_human_pairseqs')
    d = td.distances.basicDistance('AB', psDf.iloc[0], psDf.iloc[-1])

def test_compute_clones_distance():
    psDf = td.datasets.loadPSData('test_human_pairseqs')
    probDf = td.processing.computeProbs(psDf)
    psDf = psDf.join(probDf)
    clonesDf = td.processing.identifyClones(psDf)
    d = td.distances.basicDistance('AB', clonesDf.iloc[0], clonesDf.iloc[-1])
    
def test_compute_all_distances():
    psDf = td.datasets.loadPSData('test_human_pairseqs')
    probDf = td.processing.computeProbs(psDf)
    psDf = psDf.join(probDf)
    clonesDf = td.processing.identifyClones(psDf)
    clonesDf = clonesDf.set_index('CLONEID')
    pwDf = td.distances.computeBasicPWDistances('AB', clonesDf)
    assert pwDf.shape[0] == pwDf.shape[1]
    assert pwDf.shape[0] == clonesDf.shape[0]

def test_nbr_dist():
    psDf = td.datasets.loadPSData('test_human_pairseqs')
    psDf = psDf.set_index('TCRID')
    pwDf = td.distances.computeBasicPWDistances('AB', psDf)
    params = pd.objects.DistanceParams()
    VRegionDists = td.distance.computeVRegionDistances('human')
    
    d = nearestNeighborDistance(chains='AB',
                                tcr=psDf.iloc[0],
                                referenceTCRs=psDf,
                                neighborTopN=5,
                                pwDf=pwDf,
                                weighted=False,
                                VRegionDists=VRegionDists,
                                params=params)
    dw = nearestNeighborDistance(chains='AB',
                                tcr=psDf.iloc[0],
                                referenceTCRs=psDf,
                                neighborTopN=5,
                                pwDf=pwDf,
                                weighted=True,
                                VRegionDists=VRegionDists,
                                params=params)
    """Force recompute with pwDf = None"""
    d2 = nearestNeighborDistance(chains='AB',
                                tcr=psDf.iloc[0],
                                referenceTCRs=psDf,
                                neighborTopN=5,
                                pwDf=None,
                                weighted=False,
                                VRegionDists=VRegionDists,
                                params=params)
    assert d == d2

def test_embedding():
    psDf = td.datasets.loadPSData('test_human_pairseqs')
    psDf = psDf.set_index('TCRID')
    pwDf = td.distances.computeBasicPWDistances('AB', psDf)

    # xyDf = embedDistanceMatrix(pwDf, method='kpca', n_components=2)
    xyDf = td.embedding.plotEmbedding(pwDf, labels=psDf.epitope, method='kpca')