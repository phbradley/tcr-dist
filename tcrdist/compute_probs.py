import sys
import logging
import pandas as pd

from .util import get_top_genes
from . import tcr_sampler
from .translation import get_translation
from .amino_acids import amino_acids
from . import tcr_rearrangement

logger = logging.getLogger('compute_probs.py')

__all__ = ['filterOutRow',
           'samplerProb',
           'getMaskedSeqs',
           'rearrangementProb']
           
def filterOutRow(r, max_cdr3_length=30, allow_stop_codons=False, allow_X=False):
    """Assesses whether a single row from psDf should be filtered out"""
    chains = [s.split('_')[0] for s in r.index if s.find('status') > 0]
    
    fo = False
    for c in chains:
        fo = fo | r['%s_status' % c] != 'OK'
        fo = fo | 'UNK' in r['v%s_gene' % c]
        fo = fo | 'UNK' in r['j%s_gene' % c]
        fo = fo | 'TRa' in r['v%s_gene' % c]
        fo = fo | 'TRa' in r['j%s_gene' % c]
        fo = fo | len(r['cdr3%s_protseq' % c]) > max_cdr3_length
        
        for aa in r['cdr3%s_protseq' % c]:
            if not aa in amino_acids:
                assert aa in 'X*'
                if ( aa == '*' and not allow_stop_codons) or ( aa == 'X' and not allow_X ):
                    logger.debug('{} skipping: badseq: {} {}'.format(theid, cdr3a_protseq,cdr3b_protseq))
                    fo = True
                    break
    return fo

def samplerProb(r, chain):
    """Compute the probability for the CDR3 protein/nucleotide sequence.
    
    Parameter
    ---------
    r : pd.Series or TCRChain or TCRClone
        Row from the psDf or any other representation of a TCR chain
    chain : str
        Value 'a' or 'b' or 'A' or 'B'
    
    Returns
    -------
    prob_nucseq : float
        Probability based on the nucleotide sequence
    prob_protseq : float
        Probability based on the protein sequence"""
    
    c = chain.lower()
    
    if c == 'a':
        samplerFunc = tcr_sampler.alpha_cdr3_protseq_probability
    elif c == 'b':
        samplerFunc = tcr_sampler.beta_cdr3_protseq_probability
    
    prob_nucseq, new_cdr3_nucseq = samplerFunc(r.id,
                                               r.organism,
                                               r['v%s_gene'%c],
                                               r['j%s_gene'%c],
                                               cdr3_protseq='',
                                               cdr3_nucseq=r['cdr3%s_nucseq'%c],
                                               return_final_cdr3_nucseq=True)

    if new_cdr3_nucseq != r['cdr3%s_nucseq'%c]: ## note note note
        logger.warning('new_cdr3%s_nucseq: %s %s',c,len(new_cdr3_nucseq),new_cdr3_nucseq)
        logger.warning('old_cdr3%s_nucseq: %s %s',c,len(r['cdr3%s_nucseq'%c]), r['cdr3%s_nucseq'%c])
        new_cdr3_protseq = get_translation( new_cdr3_nucseq, '+1' )[0]
    else:
        new_cdr3_protseq = r['cdr3%s'%c]
        assert new_cdr3_protseq == get_translation( r['cdr3%s_nucseq'%c], '+1' )[0]
    
    prob_protseq = samplerFunc(r.id,
                               r.organism,
                               r['v%s_gene'%c],
                               r['j%s_gene'%c],
                               cdr3_protseq=new_cdr3_protseq)
    return prob_nucseq, prob_protseq

def getMaskedSeqs(r, chain):
    """Junction analysis
    
    Parameters
    ----------
    r : pd.Series or TCRChain
    chain : str
        Value 'a' or 'b' or 'A' or 'B'
    
    Returns
    -------
    cdr3_protseq_masked
    it
    cdr3_nuc_seq"""
    
    c = chain.lower()
    
    junction_results = tcr_sampler.analyze_junction( r.organism, r['v%s_gene'%c], r['j%s_gene'%c], r['cdr3%s'%c], r['cdr3%s_nucseq'%c] )

    cdr3_new_nucseq, cdr3_protseq_masked, cdr3_protseq_new_nucleotide_countstring, trims, inserts = junction_results

    """What are these checking for?"""
    #assert trims[1] + trims[2] + inserts[0] + inserts[1] + inserts[2] == 0
    #assert inserts[3] == len( cdr3_new_nucseq )

    it = '+%d-%d'%(sum(inserts[1:]),sum(trims))
    return cdr3_protseq_masked, it, cdr3_new_nucseq

def rearrangementProb(r, chain):
    """Compute rearrangement probability.
    
    Parameters
    ----------
    r : pd.Series or TCRChain
    chain : str
        Value 'a' or 'b' or 'A' or 'B'

    Returns
    -------
    v_rep_prob : float
    j_rep_prob : float"""
    
    c = chain.lower()
    
    v_countreps = r['v%s_countreps'%c].split(';')
    j_countreps = r['j%s_countreps'%c].split(';')
    v_rep_prob = max( [ tcr_rearrangement.all_countrep_pseudoprobs[r.organism][x] for x in v_countreps ] )
    j_rep_prob = max( [ tcr_rearrangement.all_countrep_pseudoprobs[r.organism][x] for x in j_countreps ] )

    return v_rep_prob, j_rep_prob  