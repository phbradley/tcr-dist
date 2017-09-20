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
           'rearrangementProb',
           'computeProbs']
           
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

def computeProbs(psDf, add_masked_seqs=True, filterOut=False, max_cdr3_length=30, allow_stop_codons=False, allow_X=False):
    """Compute probabilities for each TCR based on observed likelihood from an unpaired dataset
    and based on rearrangement likelihood
    
    Parameters
    ----------
    psDf : pd.DataFrame
        Paired-sequences dataset generated by td.processing.readPairedSequences()
        
    Returns
    -------
    psDf : pd.DataFrame
        Input dataset with new columns: a_protseq_prob, cdr3a_protseq_prob, va_rep_prob,
                                       ja_rep_prob, a_nucseq_prob, b_protseq_prob,
                                       cdr3b_protseq_prob, vb_rep_prob, jb_rep_prob, b_nucseq_prob

        optionally: cdr3a_protseq_masked, a_indels, cdr3a_new_nucseq, cdr3b_protseq_masked, b_indels, cdr3b_new_nucseq"""
    
    out = []
    for rowi, row in psDf.iterrows():
        """If iterrows is slow there are potentially ways to speed this up using psDf.apply()"""
        vals = {}
        vals['ind'] = rowi
        
        if filterOut:
            fo = filterOutRow(row,
                              max_cdr3_length=max_cdr3_length,
                              allow_stop_codons=allow_stop_codons,
                              allow_X=allow_X)
            if fo:
                """vals will be missing keys, which will be assigned Nan in outDf"""
                continue
            
        aprob_nucseq, aprob_protseq = samplerProb(row, 'a')
        va_rep_prob, ja_rep_prob = rearrangementProb(row, 'a')
        
        vals['a_protseq_prob'    ] = aprob_protseq * va_rep_prob * ja_rep_prob
        vals['cdr3a_protseq_prob'] = aprob_protseq
        vals['va_rep_prob'       ] = va_rep_prob
        vals['ja_rep_prob'       ] = ja_rep_prob
        vals['a_nucseq_prob'     ] = aprob_nucseq * va_rep_prob * ja_rep_prob
        
        bprob_nucseq, bprob_protseq = samplerProb(row, 'b')
        vb_rep_prob, jb_rep_prob = rearrangementProb(row, 'b')
        
        vals['b_protseq_prob'    ] = bprob_protseq * vb_rep_prob * jb_rep_prob
        vals['cdr3b_protseq_prob'] = bprob_protseq
        vals['vb_rep_prob'       ] = vb_rep_prob
        vals['jb_rep_prob'       ] = jb_rep_prob
        vals['b_nucseq_prob'     ] = bprob_nucseq * vb_rep_prob * jb_rep_prob
        
        if add_masked_seqs:
            cdr3a_protseq_masked, ita, cdr3a_new_nucseq = getMaskedSeqs(row, 'a')
            cdr3b_protseq_masked, itb, cdr3b_new_nucseq = getMaskedSeqs(row, 'b')

            vals[ 'cdr3a_protseq_masked'] = cdr3a_protseq_masked
            vals[ 'a_indels'] = ita
            vals[ 'cdr3a_new_nucseq' ] = cdr3a_new_nucseq
            vals[ 'cdr3b_protseq_masked'] = cdr3b_protseq_masked
            vals[ 'b_indels'] = itb
            vals[ 'cdr3b_new_nucseq' ] = cdr3b_new_nucseq
        out.append(vals)
    
    outDf = pd.DataFrame(out).set_index('ind')
    assert outDf.shape[0] == psDf.shape[0]
    return outDf
    