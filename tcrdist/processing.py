import pandas as pd
import numpy as np

from .blast import parse_unpaired_dna_sequence_blastn, get_qualstring
from .objects import TCRChain, TCRClone
from . import util

__all__ = ['processNT',
           'readPairedSequences']

def processNT(organism, chain, nuc, quals):
    """Process one nucleotide TCR sequence (any chain).

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

    ch = chain.lower()

    quals = np.array(quals.split('.')).astype(int)
    res = parse_unpaired_dna_sequence_blastn(organism, chain, nuc, info='',
                                            nocleanup=False, hide_nucseq=False,
                                            extended_cdr3=True,
                                            return_all_good_hits=True,
                                            max_bit_score_delta_for_good_hits=50)
    genes,evalues,status,all_good_hits_with_scores = res
    labels = ['v%s_gene','v%s_rep', 'v%s_mm', 'j%s_gene', 'j%s_rep', 'j%s_mm', 'cdr3%s_plus']
    tmp = {g:v for g,v in zip([lab % ch for lab in labels], genes)}
    tmp.update({'%s_evalue' % k.lower():evalues[k][0] for k in evalues.keys()})
    tmp.update({'%s_bitscore_gap' % k.lower():evalues[k][1] for k in evalues.keys()})

    tmp['%s_status' % ch] = 'OK' if not status else status
    tmp['%s_good_hits' % ch] = all_good_hits_with_scores

    tmp['cdr3%s' % ch],tmp['cdr3%s_nucseq' % ch] = tmp['cdr3%s_plus' % ch].split('-')
    tmp['cdr3%s_quals' % ch] = get_qualstring( tmp['cdr3%s_plus' % ch], nuc, quals )
    tmp['v%s_mismatches' % ch] = tmp['v%s_mm' % ch][0]
    tmp['j%s_mismatches' % ch] = tmp['j%s_mm' % ch][0]
    tmp['v%s_alignlen' % ch] = np.sum(tmp['v%s_mm' % ch])
    tmp['j%s_alignlen' % ch] = np.sum(tmp['j%s_mm' % ch])

    hits = tmp['%s_good_hits' % ch]
    if  hits and len(hits) == 2 and hits[0] and hits[1]:
        tmp['v%s_blast_hits' % ch] = ';'.join( '{}:{}'.format(x[0],x[1]) for x in hits[0] )
        tmp['j%s_blast_hits' % ch] = ';'.join( '{}:{}'.format(x[0],x[1]) for x in hits[1] )
        va_genes = util.get_top_genes( tmp['v%s_blast_hits' % ch] ) ## a set
        ja_genes = util.get_top_genes( tmp['j%s_blast_hits' % ch] )
        tmp['v%s_genes' % ch] = ';'.join( sorted( va_genes ) )
        tmp['j%s_genes' % ch] = ';'.join( sorted( ja_genes ) )
        tmp['v%s_reps' % ch]  = ';'.join( sorted( util.get_top_reps( tmp['v%s_blast_hits' % ch], organism ) ) )
        tmp['j%s_reps' % ch]  = ';'.join( sorted( util.get_top_reps( tmp['j%s_blast_hits' % ch], organism ) ) )
        tmp['v%s_countreps' % ch] = ';'.join( sorted( set( (util.get_mm1_rep_gene_for_counting(x,organism) for x in va_genes ))))
        tmp['j%s_countreps' % ch] = ';'.join( sorted( set( (util.get_mm1_rep_gene_for_counting(x,organism) for x in ja_genes ))))

    chain = TCRChain(**tmp)
    return chain

def readPairedSequences(organism, paired_seqs_file):
    """Read a TSV of paired-chain TCR data.
    Can also accomodate unpaired data or cells with more than two chains
    (e.g. gamma1 and gamma2)"""

    raw = pd.read_csv(paired_seqs_file, delimiter='\t')
    """E.g. ['a', 'b']"""
    chains = [s.split('_')[0] for s in raw.columns if s.find('nucseq') > 0]
    
    out = []
    for c in chains:
        out.append(raw.apply(lambda row: processNT(organism, c.upper(), row['%s_nucseq' % c], row['%s_quals' % c]).to_series(), axis=1))

    otherCols = [c for c in raw.columns if c.find('nucseq') == -1 and c.find('quals') == -1]
    out = [raw[otherCols]] + out
    psDf = pd.concat(out, axis=1)
    return psDf
    
def computeProbs(psDf, filter=False):
    pass




