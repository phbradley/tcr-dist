#import math
import sys
from basic import *
from util import get_top_genes
#import matplotlib
#if make_png: matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#import numpy as np
import tcr_sampler
from translation import get_translation
from amino_acids import amino_acids

new_probs = pipeline_params['new_probs']

if new_probs:
    print 'compute_probs: new_probs'
    import tcr_rearrangement_new as tcr_rearrangement ## all_rep_probs
else:
    print 'compute_probs: old_probs'
    import tcr_rearrangement



with Parser(locals()) as p:
    #p.str('args').unspecified_default().multiple().required()
    p.str('organism')
    p.str('infile').required()
    p.str('outfile').required()
    #p.int('min_mice').default(2)
    #p.float('float_arg')     # --float_arg 9.6
    #p.flag('plot')       # --flag_arg  (no argument passed)
    p.flag('verbose').shorthand('v')       # --flag_arg  (no argument passed)
    p.flag('allow_stop_codons')       # --flag_arg  (no argument passed)
    p.flag('allow_X')       # --flag_arg  (no argument passed)
    p.flag('clobber').shorthand('c')       # --flag_arg  (no argument passed)
    p.flag('add_masked_seqs')       # --flag_arg  (no argument passed)
    p.flag('filter')       # --flag_arg  (no argument passed)
    p.int('max_cdr3_length_for_filtering').default(30)       # --flag_arg  (no argument passed)
    p.flag('no_probabilities').described_as('Assign a probability of 1 to all TCRs.')
    #p.flag('find_exact_matches')       # --flag_arg  (no argument passed)
    #p.range('range_arg')     # --range_arg 1:2
    #p.multiword('multi_arg') # --multi_arg hello world
    #p.file('file_arg')       # --file_arg README.txt
    #p.directory('dir_arg')   # --dir_arg /tmp/
    #p.str('floatlist').cast(lambda x: [float(val) for val in x.split(',')])
    #p.multiword('intlist').cast(lambda x: [int(val) for val in x.split()])

assert add_masked_seqs

if exists(outfile):assert clobber

out = open(outfile,'w')

infields = []
outfields = []

for line in open(infile,'rU'):
    if not infields:
        if line[0] == '#':
            infields = line[1:-1].split('\t')
        else:
            infields = line[:-1].split('\t')

        outfields = infields[:]
        outfields.extend( ['a_protseq_prob','cdr3a_protseq_prob','va_rep_prob','ja_rep_prob','a_nucseq_prob',
                           'b_protseq_prob','cdr3b_protseq_prob','vb_rep_prob','jb_rep_prob','b_nucseq_prob' ] )

        if add_masked_seqs:
            outfields.extend( ['cdr3a_protseq_masked','a_indels','cdr3a_new_nucseq',
                               'cdr3b_protseq_masked','b_indels','cdr3b_new_nucseq' ] )

        out.write('\t'.join( outfields )+'\n' )
        continue
    assert infields

    l = parse_tsv_line( line[:-1], infields )

    if filter and 'status' in l and l['status'] != 'OK':continue

    theid = line.split("\t")[0]

    ## ALPHA
    va_gene = l['va_gene']
    ja_gene = l['ja_gene']
    vb_gene = l['vb_gene']
    jb_gene = l['jb_gene']
    cdr3a_protseq = l['cdr3a']
    cdr3a_nucseq  = l['cdr3a_nucseq']
    cdr3b_protseq = l['cdr3b']
    cdr3b_nucseq  = l['cdr3b_nucseq']

    if filter:
        if 'UNK' in va_gene+ja_gene or 'TRa' in va_gene+ja_gene: continue
        if 'UNK' in va_gene+ja_gene or 'TRa' in va_gene+ja_gene: continue
        if len(cdr3a_protseq)>max_cdr3_length_for_filtering:continue
        if len(cdr3b_protseq)>max_cdr3_length_for_filtering:continue

        ## check for stop codons
        skip_me = False
        for a in cdr3a_protseq + cdr3b_protseq:
            if a not in amino_acids:
                assert a in 'X*'
                if ( a == '*' and not allow_stop_codons) or ( a == 'X' and not allow_X ):
                    Log('{} skipping: badseq: {} {}'.format(theid, cdr3a_protseq,cdr3b_protseq))
                    skip_me = True
                    break

        if skip_me:
            continue

    ## probs are computed by reps
    va_reps = l['va_reps'].split(';')
    ja_reps = l['ja_reps'].split(';')
    va_countreps = l['va_countreps'].split(';')
    ja_countreps = l['ja_countreps'].split(';')
    va_cdr3_nucseq = tcr_sampler.get_v_cdr3_nucseq( organism, va_gene )
    ja_cdr3_nucseq = tcr_sampler.get_j_cdr3_nucseq( organism, ja_gene )
    va_cdr3_protseq,codons = get_translation( va_cdr3_nucseq, '+1' )
    ja_cdr3_protseq,codons = get_translation( ja_cdr3_nucseq, '+{}'.format(1+len(ja_cdr3_nucseq)%3))

    if no_probabilities or not tcr_rearrangement.probs_data_exist( organism,'A'):
        ##all probabilities will be set to 1 if this flag is set
        aprob_nucseq = 1
        aprob_protseq = 1
    else:
        aprob_nucseq,new_cdr3a_nucseq = tcr_sampler.alpha_cdr3_protseq_probability( theid, organism, va_gene, ja_gene,
                                                                                cdr3_protseq='',
                                                                                cdr3_nucseq=cdr3a_nucseq,  verbose=verbose,
                                                                                return_final_cdr3_nucseq=True )

        if new_cdr3a_nucseq != cdr3a_nucseq: ## note note note
            print 'new_cdr3a_nucseq:',len(new_cdr3a_nucseq),new_cdr3a_nucseq
            print 'old_cdr3a_nucseq:',len(cdr3a_nucseq),cdr3a_nucseq
            new_cdr3a_protseq = get_translation( new_cdr3a_nucseq, '+1' )[0]
        else:
            new_cdr3a_protseq = cdr3a_protseq[:]
            assert new_cdr3a_protseq == get_translation( cdr3a_nucseq, '+1' )[0]

        aprob_protseq = tcr_sampler.alpha_cdr3_protseq_probability( theid, organism, va_gene, ja_gene, new_cdr3a_protseq,
                                                                verbose=verbose )

    ## BETA
    vb_reps = l['vb_reps'].split(';')
    jb_reps = l['jb_reps'].split(';')
    vb_countreps = l['vb_countreps'].split(';')
    jb_countreps = l['jb_countreps'].split(';')
    vb_cdr3_nucseq = tcr_sampler.get_v_cdr3_nucseq( organism, vb_gene )
    jb_cdr3_nucseq = tcr_sampler.get_j_cdr3_nucseq( organism, jb_gene )
    vb_cdr3_protseq,codons = get_translation( vb_cdr3_nucseq, '+1' )
    jb_cdr3_protseq,codons = get_translation( jb_cdr3_nucseq, '+{}'.format(1+len(jb_cdr3_nucseq)%3))

    if no_probabilities or not tcr_rearrangement.probs_data_exist( organism,'B'):
        ##all probabilities will be set to 1 if this flag is set
        bprob_nucseq = 1
        bprob_protseq = 1
    else:
        bprob_nucseq, new_cdr3b_nucseq \
        = tcr_sampler.beta_cdr3_protseq_probability( theid, organism, vb_gene, jb_gene, cdr3_protseq='',
                                                     verbose=verbose, cdr3_nucseq=cdr3b_nucseq,
                                                     allow_early_nucseq_mismatches=True,
                                                     return_final_cdr3_nucseq=True )

        if new_cdr3b_nucseq != cdr3b_nucseq: ## note note note
            new_cdr3b_protseq = get_translation( new_cdr3b_nucseq, '+1' )[0]
        else:
            new_cdr3b_protseq = cdr3b_protseq[:]
            assert new_cdr3b_protseq == get_translation( cdr3b_nucseq, '+1' )[0]

        bprob_protseq = tcr_sampler.beta_cdr3_protseq_probability( theid, organism, vb_gene, jb_gene, new_cdr3b_protseq,
                                                               verbose=verbose )

    vals = dict(l)  #line.split('\t') + ['']*(len(outfields)-len(infields))

    if add_masked_seqs:
        ## junction analysis
        a_junction_results = tcr_sampler.analyze_junction( organism, va_gene, ja_gene, cdr3a_protseq, cdr3a_nucseq )
        b_junction_results = tcr_sampler.analyze_junction( organism, vb_gene, jb_gene, cdr3b_protseq, cdr3b_nucseq )

        cdr3a_new_nucseq, cdr3a_protseq_masked, cdr3a_protseq_new_nucleotide_countstring,a_trims,a_inserts \
            = a_junction_results
        cdr3b_new_nucseq, cdr3b_protseq_masked, cdr3b_protseq_new_nucleotide_countstring,b_trims,b_inserts \
            = b_junction_results

        # from tcr_sampler.py:
        # trims = ( v_trim, d0_trim, d1_trim, j_trim )
        # inserts = ( best_d_id, n_vd_insert, n_dj_insert, n_vj_insert )

        assert a_trims[1] + a_trims[2] + a_inserts[0] + a_inserts[1] + a_inserts[2] + b_inserts[3] == 0
        assert a_inserts[3] == len( cdr3a_new_nucseq )

        ita = '+%d-%d'%(sum(a_inserts[1:]),sum(a_trims))
        itb = '+%d-%d'%(sum(b_inserts[1:]),sum(b_trims))

        vals[ 'cdr3a_protseq_masked'] = cdr3a_protseq_masked
        vals[ 'a_indels'] = ita
        vals[ 'cdr3a_new_nucseq' ] = cdr3a_new_nucseq
        vals[ 'cdr3b_protseq_masked'] = cdr3b_protseq_masked
        vals[ 'b_indels'] = itb
        vals[ 'cdr3b_new_nucseq' ] = cdr3b_new_nucseq

    ## there's a little bit of a bias toward guys with more blast hits, ie shorted reads? since we take a max
    if no_probabilities or not ( tcr_rearrangement.probs_data_exist(organism,'A') and
                                 tcr_rearrangement.probs_data_exist(organism,'B') ):
        ##all probabilities will be set to 1 if this flag is set
        va_rep_prob = 1
        ja_rep_prob = 1
        vb_rep_prob = 1
        jb_rep_prob = 1
    elif new_probs:
        va_rep_prob = max( [ tcr_rearrangement.all_countrep_pseudoprobs[organism]['A']['V'][x] for x in va_countreps ] )
        ja_rep_prob = max( [ tcr_rearrangement.all_countrep_pseudoprobs[organism]['A']['J'][x] for x in ja_countreps ] )
        vb_rep_prob = max( [ tcr_rearrangement.all_countrep_pseudoprobs[organism]['B']['V'][x] for x in vb_countreps ] )
        jb_rep_prob = max( [ tcr_rearrangement.all_countrep_pseudoprobs[organism]['B']['J'][x] for x in jb_countreps ] )
    else:
        va_rep_prob = max( [ tcr_rearrangement.all_rep_probs[organism][x] for x in va_reps ] )
        ja_rep_prob = max( [ tcr_rearrangement.all_rep_probs[organism][x] for x in ja_reps ] )
        vb_rep_prob = max( [ tcr_rearrangement.all_rep_probs[organism][x] for x in vb_reps ] )
        jb_rep_prob = max( [ tcr_rearrangement.all_rep_probs[organism][x] for x in jb_reps ] )

    vals['a_protseq_prob'    ] = aprob_protseq * va_rep_prob * ja_rep_prob
    vals['cdr3a_protseq_prob'] = aprob_protseq
    vals['va_rep_prob'       ] = va_rep_prob
    vals['ja_rep_prob'       ] = ja_rep_prob
    vals['a_nucseq_prob'     ] = aprob_nucseq * va_rep_prob * ja_rep_prob

    vals['b_protseq_prob'    ] = bprob_protseq * vb_rep_prob * jb_rep_prob
    vals['cdr3b_protseq_prob'] = bprob_protseq
    vals['vb_rep_prob'       ] = vb_rep_prob
    vals['jb_rep_prob'       ] = jb_rep_prob
    vals['b_nucseq_prob'     ] = bprob_nucseq * vb_rep_prob * jb_rep_prob

    assert len(vals.keys()) == len(outfields)

    out.write( make_tsv_line( vals, outfields, '-' )+'\n' )
    out.flush()

out.close()
