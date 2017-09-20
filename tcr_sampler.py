from basic import *
import sys
import translation
from all_genes import all_genes, gap_character
from genetic_code import genetic_code, reverse_genetic_code
import logo_tools
import random

if pipeline_params['new_probs']:
    Log( 'tcr_sampler: new_probs' )
    import tcr_rearrangement_new as tcr_rearrangement
else:
    Log( 'tcr_sampler: old_probs' )
    import tcr_rearrangement


########################################################################################################################
default_mismatch_score_for_cdr3_nucseq_probabilities = -4 ## blast is -3
default_mismatch_score_for_junction_analysis = -4 ## blast is -3

def count_matches( a,b,mismatch_score=-3 ):
    assert a[0].lower() == a[0]
    #assert a[0].upper() == a[0]
    ## from the beginning
    match_score    = 1
    best_score=0
    score=0
    num_matches = 0
    for i in range(min(len(a),len(b))):
        if a[i] == b[i] or logo_tools.nuc_match_lower_case.get( (a[i],b[i]), False ):
            score += match_score
        else:
            score += mismatch_score
        if score >= best_score: ## gt OR EQUAL! take longer matched regions
            best_score = score
            num_matches = i+1
    return num_matches



def get_v_cdr3_nucseq( organism, v_gene, paranoid = False ):
    vg = all_genes[organism][v_gene]
    ab = vg.chain

    v_nucseq  = vg.nucseq
    v_nucseq_offset = vg.nucseq_offset
    v_nucseq = v_nucseq[ v_nucseq_offset: ]

    v_alseq = vg.alseq
    alseq_cpos = vg.cdr_columns[-1][0] -1
    #print organism, v_gene, old_v_alseq
    #print organism, v_gene, v_alseq
    numgaps = v_alseq[:alseq_cpos].count('.')
    v_cpos = alseq_cpos - numgaps
    v_nucseq = v_nucseq[3*v_cpos:] ## now v_nucseq starts with the 'C' codon


    # the following hacky adjustment is now incorporated in the dbfile
    # if organism == 'mouse':
    #     if v_gene == 'TRAV13D-1*01':
    #         #-----------------------------------
    #         #../../../tmp.blat:mismatch: V 6 imgt: a genome: t TRAV13D-1*01
    #         #tmp.whoah:whoah  6 act: t  98.7 exp: a   1.1 TRAV13D-1*01 TRAV13-1*01 620
    #         #tmp.whoah:whoah: expected: caaggtatcgtgt consensus: caaggtttcgtgt TRAV13D-1*01 620
    #         #tmp.3.whoah:whoah  6 act: t  97.4 exp: a   1.4 TRAV13D-1*01 TRAV13-1*01 642
    #         #tmp.3.whoah:whoah: expected: caaggtatcgtgt consensus: caaggtttcgtgt TRAV13D-1*01 642
    #         #tmp.la_mc.whoah:whoah  6 act: t  89.0 exp: a   7.0 TRAV13D-1*01 TRAV13-1*01 100
    #         #tmp.la_mc.whoah:whoah: expected: caaggtatcgtgt consensus: caaggtttcgtgt TRAV13D-1*01 100
    #         assert v_nucseq == 'tgtgctatggaac' ## CAM ## THIS WILL FAIL SINCE WE ADDED THIS TO THE DB...
    #         v_nucseq         = 'tgtgctttggaac' ## CAL


    return v_nucseq


def get_j_cdr3_nucseq( organism, j_gene, paranoid = False ):
    jg = all_genes[organism][j_gene]
    ab = jg.chain

    j_nucseq  = jg.nucseq
    j_nucseq_offset = jg.nucseq_offset

    ## goes up to (but not including) GXG
    num_genome_j_aas_in_loop = len( jg.cdrs[0].replace(gap_character,''))

    ## trim j_nucseq so that it extends up to the F/W position
    j_nucseq = j_nucseq[:3*num_genome_j_aas_in_loop + j_nucseq_offset]


    # the following hacky adjustments are now incorporated in the dbfile
    # if organism == 'mouse':
    #     if j_gene == 'TRAJ47*01':
    #         # -----------------------------------
    #         # ../../../tmp.blat:mismatch: J 2 imgt: c genome: g TRAJ47*01
    #         # ../../../tmp.blat:mismatch: J 24 imgt: g genome: t TRAJ47*01
    #         # tmp.whoah:whoah  2 act: g  81.9 exp: c   4.7 TRAJ47*01 TRAJ47*01 1412
    #         # tmp.whoah:whoah 24 act: t  82.7 exp: g  16.8 TRAJ47*01 TRAJ47*01 1412
    #         # tmp.whoah:whoah: expected: tgcactatgcaaacaagatgatctgt consensus: tggactatgcaaacaagatgatcttt TRAJ47*01 1412
    #         # tmp.3.whoah:whoah  2 act: g  81.6 exp: c   5.0 TRAJ47*01 TRAJ47*01 1362
    #         # tmp.3.whoah:whoah 24 act: t  82.7 exp: g  16.6 TRAJ47*01 TRAJ47*01 1362
    #         # tmp.3.whoah:whoah: expected: tgcactatgcaaacaagatgatctgt consensus: tggactatgcaaacaagatgatcttt TRAJ47*01 1362
    #         # tmp.la_mc.whoah:whoah  2 act: g  79.6 exp: c   5.3 TRAJ47*01 TRAJ47*01 113
    #         # tmp.la_mc.whoah:whoah 24 act: t  99.1 exp: g   0.9 TRAJ47*01 TRAJ47*01 113
    #         # tmp.la_mc.whoah:whoah: expected: tgcactatgcaaacaagatgatctgt consensus: tggactatgcaaacaagatgatcttt TRAJ47*01 113
    #         assert j_nucseq == 'tgcactatgcaaacaagatgatctgt' ## C at end
    #         j_nucseq         = 'tggactatgcaaacaagatgatcttt' ## F at end
    #     elif j_gene == 'TRAJ24*01':
    #         # -----------------------------------
    #         # ../../../tmp.blat:unaligned: J 0 TRAJ24*01
    #         # ../../../tmp.blat:unaligned: J 1 TRAJ24*01
    #         # ../../../tmp.blat:gapped: J 6 TRAJ24*01
    #         # tmp.whoah:whoah  2 act: c  60.3 exp: a  15.3 TRAJ24*01 TRAJ24*01 464
    #         # tmp.whoah:whoah  4 act: a  88.6 exp: c   2.8 TRAJ24*01 TRAJ24*01 464
    #         # tmp.whoah:whoah  5 act: c  93.3 exp: t   1.5 TRAJ24*01 TRAJ24*01 464
    #         # tmp.whoah:whoah  6 act: t  97.2 exp: g   1.1 TRAJ24*01 TRAJ24*01 464
    #         # tmp.whoah:whoah: expected: tgaactggccagtttggggaaactgcagttt consensus: gacaactgccagtttggggaaactgcagttt TRAJ24*01 464
    #         # tmp.3.whoah:whoah  2 act: c  60.8 exp: a  13.9 TRAJ24*01 TRAJ24*01 475
    #         # tmp.3.whoah:whoah  4 act: a  86.3 exp: c   4.2 TRAJ24*01 TRAJ24*01 475
    #         # tmp.3.whoah:whoah  5 act: c  94.5 exp: t   1.1 TRAJ24*01 TRAJ24*01 475
    #         # tmp.3.whoah:whoah  6 act: t  98.1 exp: g   1.1 TRAJ24*01 TRAJ24*01 475
    #         # tmp.3.whoah:whoah: expected: tgaactggccagtttggggaaactgcagttt consensus: gacaactgccagtttggggaaactgcagttt TRAJ24*01 475
    #         # tmp.la_mc.whoah:whoah  2 act: c  75.3 exp: a   4.3 TRAJ24*01 TRAJ24*01 93
    #         # tmp.la_mc.whoah:whoah  4 act: a  89.2 exp: c   2.2 TRAJ24*01 TRAJ24*01 93
    #         # tmp.la_mc.whoah:whoah  5 act: c  97.8 exp: t   1.1 TRAJ24*01 TRAJ24*01 93
    #         # tmp.la_mc.whoah:whoah  6 act: t  98.9 exp: g   0.0 TRAJ24*01 TRAJ24*01 93
    #         # tmp.la_mc.whoah:whoah: expected: tgaactggccagtttggggaaactgcagttt consensus: gacaactgccagtttggggaaactgcagttt TRAJ24*01 93
    #         assert j_nucseq == 'tgaactggccagtttggggaaactgcagttt'
    #         j_nucseq         = 'gacaactgccagtttggggaaactgcagttt'
    #         ## take the consensus
    #         ## given that there's an indel (and the alignment to the genome starts at j sequence position 3)
    #         ## it's hard to tell what to do at the beginning...



    return j_nucseq




def analyze_junction( organism, v_gene, j_gene, cdr3_protseq, cdr3_nucseq, force_d_id=0, return_cdr3_nucseq_src=False ):
    #print organism, v_gene, j_gene, cdr3_protseq, cdr3_nucseq
    assert v_gene.startswith('TR') #and v_gene[2] == j_gene[2]
    ab = all_genes[organism][v_gene].chain
    v_nucseq = get_v_cdr3_nucseq( organism, v_gene )
    j_nucseq = get_j_cdr3_nucseq( organism, j_gene )
    ## how far out do we match
    num_matched_v = count_matches( v_nucseq, cdr3_nucseq, default_mismatch_score_for_junction_analysis )

    num_matched_j = count_matches( ''.join( reversed( list( j_nucseq ) )),
                                   ''.join( reversed( list( cdr3_nucseq ))),
                                   default_mismatch_score_for_junction_analysis )


    if num_matched_v + num_matched_j > len(cdr3_nucseq):
        ## some overlap!
        extra = num_matched_v + num_matched_j - len(cdr3_nucseq )
        fake_v_trim = extra/2 ## now deterministic
        fake_j_trim = extra - fake_v_trim
        num_matched_v -= fake_v_trim
        num_matched_j -= fake_j_trim

    assert num_matched_v + num_matched_j <= len(cdr3_nucseq)

    if num_matched_v + num_matched_j == len(cdr3_nucseq):
        nseq = ''
    else:
        nseq = cdr3_nucseq[num_matched_v:len(cdr3_nucseq)-num_matched_j]

    ncount = [1] * len(cdr3_nucseq)
    cdr3_nucseq_src = ['N'] * len(cdr3_nucseq)
    for i in range(num_matched_v):
        ncount[i] = 0
        cdr3_nucseq_src[i] = 'V'
    for i in range(num_matched_j):
        ncount[-1-i] = 0
        cdr3_nucseq_src[-1-i] = 'J'

    assert num_matched_v + num_matched_j + len(nseq) == len(cdr3_nucseq)

    v_trim = len(v_nucseq)-num_matched_v
    j_trim = len(j_nucseq)-num_matched_j

    assert len(cdr3_nucseq) == len(v_nucseq) + len(nseq) + len(j_nucseq) - ( v_trim + j_trim )

    #d_info = ''
    n_vj_insert = 0
    n_vd_insert = 0
    n_dj_insert = 0
    d0_trim = 0
    d1_trim = 0
    best_d_id = 0

    if ab == 'A':
        n_vj_insert = len(nseq)

    elif ab == 'B':
        ## look for one of the d-gene segments
        max_overlap = 0
        for d_id, d_nucseq in tcr_rearrangement.all_trbd_nucseq[organism].iteritems():
            if force_d_id and d_id != force_d_id: continue
            for start in range(len(d_nucseq)):
                for stop in range(start,len(d_nucseq)):
                    overlap_seq = d_nucseq[start:stop+1]
                    if overlap_seq in nseq:
                        if len(overlap_seq)>max_overlap:
                            max_overlap = len(overlap_seq)
                            best_d_id = d_id
                            best_overlap_seq = overlap_seq
                            best_trim = (start,len(d_nucseq)-1-stop)

        if max_overlap: ## found a bit of d, although it might be bogus (eg 1 nt)
            pos = nseq.index( best_overlap_seq )
            for i in range(pos+num_matched_v,pos+num_matched_v+max_overlap):
                assert ncount[i] == 1
                ncount[i] = 0
                cdr3_nucseq_src[i] = 'D'
                nseq = nseq[:i-num_matched_v] + '+' + nseq[i+1-num_matched_v:]

            n_vd_insert = pos
            n_dj_insert = len(nseq) - pos - len(best_overlap_seq)
            d0_trim = best_trim[0]
            d1_trim = best_trim[1]

            expected_cdr3_nucseq_len = ( len(v_nucseq) + n_vd_insert +
                                         len(tcr_rearrangement.all_trbd_nucseq[organism][best_d_id]) + n_dj_insert +
                                         len(j_nucseq) -
                                         ( v_trim + d0_trim + d1_trim + j_trim ) )
            assert len(cdr3_nucseq) == expected_cdr3_nucseq_len


        else:
            best_d_id = 0
            n_vd_insert = 0
            n_dj_insert = 0
            d0_trim = 0
            d1_trim = 0




    if cdr3_protseq:
        assert 3*len(cdr3_protseq) == len(ncount)

        newseq = ''
        newseq_ncount = ''

        for i,a in enumerate(cdr3_protseq):
            nc = sum(ncount[3*i:3*i+3])
            newseq_ncount += `nc`
            if nc>1:
                newseq += a
            elif nc==1:
                newseq += a.lower()
            else:
                newseq += '-'

        ## output
        cdr3_protseq_masked = newseq[:]
        cdr3_protseq_new_nucleotide_countstring = newseq_ncount[:]
    else:## no protseq given (perhaps its an out of frame?)
        cdr3_protseq_masked = ''
        cdr3_protseq_new_nucleotide_countstring = ''

    new_nucseq = nseq[:]

    trims = ( v_trim, d0_trim, d1_trim, j_trim )
    inserts = ( best_d_id, n_vd_insert, n_dj_insert, n_vj_insert )

    ## new_nucseq spans the inserted nucleotide sequence and has '+' for D-nucleotides
    if return_cdr3_nucseq_src:
        return new_nucseq, cdr3_protseq_masked, cdr3_protseq_new_nucleotide_countstring, trims, inserts, cdr3_nucseq_src
    else:
        return new_nucseq, cdr3_protseq_masked, cdr3_protseq_new_nucleotide_countstring, trims, inserts


## what fraction of all the sequences that match nucseq at non-'n'-positions will code for protseq?
def get_coding_probability( nucseq, protseq ):
    assert len(nucseq) == 3*len(protseq)

    #probs = []
    total_prob = 1.0
    for i in range(len(protseq)):
        aa = protseq[i]
        ncodon = nucseq[3*i:3*i+3]
        ## what fraction of the possible codons here would actually code for the desired aa
        count=0
        for c in reverse_genetic_code[aa]:
            match = True
            for x,y in zip(c,ncodon):
                ## want to allow for the possibility of 'wskmyr' symbols in nucseq?
                if x!= y and y != 'n' and not logo_tools.nuc_match_lower_case.get( (x,y), False ):
                    match = False
                    break
            if match:
                count+=1
        if count==0:
            total_prob=0
            break
        num_n = ncodon.count('n')
        total_codons = 4**num_n
        prob = float(count)/total_codons
        #if num_n:
        #    print 'ncodon:',ncodon,aa,count,total_codons,prob
        #probs.append( prob )
        total_prob *= prob
    return total_prob



def alpha_cdr3_protseq_probability( theid, organism, v_gene, j_gene, cdr3_protseq,
                                    cdr3_nucseq = '', error_threshold = 0.05, verbose=False,
                                    allow_early_nucseq_mismatches = True,
                                    return_final_cdr3_nucseq = False ):
    nucleotide_match = ( cdr3_nucseq != '' )
    if nucleotide_match:
        assert not cdr3_protseq
        cdr3_protseq = translation.get_translation( cdr3_nucseq, '+1' )[0]
        assert len(cdr3_nucseq) == 3 * len(cdr3_protseq )

    ab = 'A'
    assert all_genes[organism][v_gene].chain == ab

    v_nucseq = get_v_cdr3_nucseq( organism, v_gene )
    j_nucseq = get_j_cdr3_nucseq( organism, j_gene )

    ## what is the largest amount of these nucseqs we could preserve and still get cdr3_protseq

    max_v_germline = 0
    len_v_nucseq = len(v_nucseq)
    max_j_germline = 0

    len_j_nucseq = len(j_nucseq)
    len_cdr3_protseq = len(cdr3_protseq)
    len_cdr3_nucseq = len(cdr3_nucseq)

    if nucleotide_match:
        if allow_early_nucseq_mismatches:
            mismatch_score = default_mismatch_score_for_cdr3_nucseq_probabilities
        else:
            mismatch_score = -100
        max_v_germline = count_matches( v_nucseq, cdr3_nucseq, mismatch_score )

        max_j_germline = count_matches( ''.join( reversed( list( j_nucseq ) )),
                                        ''.join( reversed( list( cdr3_nucseq ))),
                                        mismatch_score )

        if allow_early_nucseq_mismatches: ## obliterate the mismatches now
            max_v, max_j = max_v_germline, max_j_germline
            if max_v + max_j > len(cdr3_nucseq):
                ## some overlap!
                extra = max_v + max_j - len(cdr3_nucseq )
                #print 'TRIM extra',extra
                fake_v_trim = extra/2 ## now dterministic
                fake_j_trim = extra - fake_v_trim
                max_v -= fake_v_trim
                max_j -= fake_j_trim
            old_cdr3_nucseq = cdr3_nucseq[:]
            cdr3_nucseq = v_nucseq[:max_v] + \
                          cdr3_nucseq[ max_v : len_cdr3_nucseq-max_j ] + \
                          j_nucseq[len_j_nucseq-max_j:]
            if old_cdr3_nucseq != cdr3_nucseq:
                Log('{} early_cdr3a_nucseq_mismatch: {} {} before {} after {}'.format(theid, v_gene, j_gene,
                                                                                    old_cdr3_nucseq, cdr3_nucseq ) )
                assert len(cdr3_nucseq) == len(old_cdr3_nucseq)
    else:

        for i in range( len(v_nucseq)):
            i_aa = i/3 ## which aa do we code for?
            len_codon = (i%3) + 1
            if i_aa >= len(cdr3_protseq): break
            start = 3*i_aa
            codon = v_nucseq[ start:start+len_codon]
            target_aa = cdr3_protseq[ i_aa ]
            matched = False
            for c in reverse_genetic_code[target_aa]:
                if c.startswith(codon):
                    matched = True
            if verbose:
                print 'V',codon, target_aa, matched
            if matched:
                max_v_germline = i+1
            else:
                break

        ## how about J?
        for i in range( len_j_nucseq):
            i_aa = i/3 ## which aa do we code for?
            len_codon = (i%3) + 1
            if i_aa >= len(cdr3_protseq): break
            end   = len(j_nucseq)-3*i_aa
            codon = j_nucseq[max(0,end-len_codon):end]
            target_aa = cdr3_protseq[ len_cdr3_protseq-1-i_aa ]
            matched = False
            for c in reverse_genetic_code[target_aa]:
                if c.endswith(codon):
                    matched = True
            if verbose:
                print 'J',codon, target_aa, matched
            if matched:
                max_j_germline = i+1
            else:
                break


    min_insert = 3*len_cdr3_protseq - max_v_germline - max_j_germline
    if verbose:
        print 'max_v_germline:',max_v_germline, len_v_nucseq, v_nucseq, cdr3_nucseq

        print 'max_j_germline:',max_j_germline, len_j_nucseq, j_nucseq, cdr3_nucseq, \
            all_genes[organism][j_gene].protseq

        print 'min_insert:',min_insert,max_v_germline,max_j_germline

    total_prob = 0.0
    min_extra_trim = max(0,-1*min_insert)
    for extra_trim in range(min_extra_trim,100):
        old_total_prob = total_prob
        total_prob_this_trim = 0.0
        for extra_v_trim in range(0,extra_trim+1):
            extra_j_trim = extra_trim - extra_v_trim

            v_trim = len_v_nucseq - max_v_germline + extra_v_trim
            j_trim = len_j_nucseq - max_j_germline + extra_j_trim
            if v_trim > len_v_nucseq or j_trim > len_j_nucseq: continue

            n_insert = min_insert + extra_v_trim + extra_j_trim
            n_nucseq = v_nucseq[:len_v_nucseq-v_trim] + 'n'*n_insert + j_nucseq[j_trim:]

            assert len(n_nucseq) == 3*len_cdr3_protseq
            if nucleotide_match:
                coding_prob = 0.25 ** n_insert
            else:
                coding_prob = get_coding_probability( n_nucseq, cdr3_protseq )

            trim_prob = tcr_rearrangement.get_alpha_trim_probs( organism, v_trim, j_trim, n_insert )

            total_prob_this_trim += coding_prob * trim_prob
            total_prob += coding_prob * trim_prob

            if verbose:
                print 'coding_prob:',cdr3_protseq,v_trim,j_trim,n_insert,total_prob,coding_prob,trim_prob,n_nucseq

        if extra_trim>2 and total_prob_this_trim < error_threshold * old_total_prob:
            break
    if return_final_cdr3_nucseq:
        return total_prob, cdr3_nucseq
    else:
        return total_prob



def beta_cdr3_protseq_probability( theid, organism, v_gene, j_gene, cdr3_protseq,
                                   cdr3_nucseq = '', error_threshold = 0.05, verbose=False,
                                   allow_early_nucseq_mismatches = True,
                                   return_final_cdr3_nucseq = False ):
    nucleotide_match = ( cdr3_nucseq != '' )
    if nucleotide_match:
        assert not cdr3_protseq
        cdr3_protseq = translation.get_translation( cdr3_nucseq, '+1' )[0]
        assert len(cdr3_nucseq) == 3 * len(cdr3_protseq )

    ab = 'B'
    assert all_genes[organism][v_gene].chain == ab

    v_nucseq = get_v_cdr3_nucseq( organism, v_gene )
    j_nucseq = get_j_cdr3_nucseq( organism, j_gene )

    ## what is the largest amount of these nucseqs we could preserve and still get cdr3_protseq
    max_v_germline = 0
    max_j_germline = 0

    len_v_nucseq = len(v_nucseq)
    len_j_nucseq = len(j_nucseq)
    len_cdr3_nucseq = len(cdr3_nucseq)
    len_cdr3_protseq = len(cdr3_protseq)

    if nucleotide_match:
        if allow_early_nucseq_mismatches:
            mismatch_score = default_mismatch_score_for_cdr3_nucseq_probabilities
        else:
            mismatch_score = -100
        max_v_germline = count_matches( v_nucseq, cdr3_nucseq, mismatch_score )

        max_j_germline = count_matches( ''.join( reversed( list( j_nucseq ) )),
                                        ''.join( reversed( list( cdr3_nucseq ))),
                                        mismatch_score )

        if allow_early_nucseq_mismatches: ## obliterate the mismatches now
            max_v, max_j = max_v_germline, max_j_germline
            if max_v + max_j > len(cdr3_nucseq):
                ## some overlap!
                extra = max_v + max_j - len(cdr3_nucseq )
                #print 'TRIM extra',extra
                fake_v_trim = extra/2 ## now dterministic
                fake_j_trim = extra - fake_v_trim
                max_v -= fake_v_trim
                max_j -= fake_j_trim
            old_cdr3_nucseq = cdr3_nucseq[:]
            cdr3_nucseq = v_nucseq[:max_v] + \
                          cdr3_nucseq[ max_v : len_cdr3_nucseq-max_j ] + \
                          j_nucseq[len_j_nucseq-max_j:]
            if old_cdr3_nucseq != cdr3_nucseq:
                Log('{} early_cdr3a_nucseq_mismatch: before {} after {}'.format(theid, old_cdr3_nucseq, cdr3_nucseq ) )
                assert len(cdr3_nucseq) == len(old_cdr3_nucseq)

    else:
        ## V
        for i in range( len(v_nucseq)):
            i_aa = i/3 ## which aa do we code for?
            len_codon = (i%3) + 1
            if i_aa >= len(cdr3_protseq): break
            start = 3*i_aa
            codon = v_nucseq[ start:start+len_codon]
            target_aa = cdr3_protseq[ i_aa ]
            matched = False
            for c in reverse_genetic_code[target_aa]:
                if c.startswith(codon):
                    matched = True
            if verbose:
                print 'V',codon, target_aa, matched
            if matched:
                max_v_germline = i+1
            else:
                break

        ## J
        for i in range( len_j_nucseq):
            i_aa = i/3 ## which aa do we code for?
            len_codon = (i%3) + 1
            if i_aa >= len(cdr3_protseq): break
            end   = len(j_nucseq)-3*i_aa
            codon = j_nucseq[max(0,end-len_codon):end]
            target_aa = cdr3_protseq[ len_cdr3_protseq-1-i_aa ]
            matched = False
            for c in reverse_genetic_code[target_aa]:
                if c.endswith(codon):
                    matched = True
            if verbose:
                print 'J',codon, target_aa, matched
            if matched:
                max_j_germline = i+1
            else:
                break


    if verbose:
        print 'max_v_germline:',max_v_germline, len(v_nucseq)

    ## how about J?

    min_insert = 3*len_cdr3_protseq - max_v_germline - max_j_germline
    if verbose:
        print 'max_j_germline:',max_j_germline, len_j_nucseq,cdr3_protseq,\
            all_genes[organism][j_gene].protseq

        print 'min_insert:',min_insert,max_v_germline,max_j_germline

    if organism in ['human','mouse'] and j_gene[3] == 'B':
        trbj_index = int( j_gene[4] ) ## to decide which d genes to allow
        assert trbj_index in [1,2]
    else:
        ## no D/J compatibility check
        trbj_index=0

    total_prob = 0.0
    min_extra_trim = max(0,-1*min_insert)

    dids = tcr_rearrangement.all_trbd_nucseq[organism].keys()
    for extra_trim in range(min_extra_trim,100):
        old_total_prob = total_prob
        total_prob_this_trim = 0.0
        for extra_v_trim in range(0,extra_trim+1):
            extra_j_trim = extra_trim - extra_v_trim

            v_trim = len_v_nucseq - max_v_germline + extra_v_trim
            j_trim = len_j_nucseq - max_j_germline + extra_j_trim
            if v_trim > len_v_nucseq or j_trim > len_j_nucseq: continue

            n_insert = min_insert + extra_v_trim + extra_j_trim
            assert n_insert>=0 ## b/c of min_extra_trim
            total_prob_this_insert = 0.0

            ## now we are looking to fit part of the D gene into this middle region and still code for the right aas
            for did in dids:
                if trbj_index == 1:
                    if did == 1:
                        did_prob = 1.0
                    else:
                        continue
                else:
                    did_prob = 1.0/float(len(dids))
                d_nucseq = tcr_rearrangement.all_trbd_nucseq[organism][did]
                len_d_nucseq = len( d_nucseq )
                for d0_trim in range(len_d_nucseq+1):
                    for d1_trim in range(len_d_nucseq+1):
                        len_d_insert = len_d_nucseq - d0_trim - d1_trim
                        if len_d_insert < 0 or len_d_insert > n_insert: continue
                        #if len_d_insert == 0 and d1_trim: continue ## only hit this one once!
                        d_insert = d_nucseq[ d0_trim: len_d_nucseq-d1_trim]
                        num_n = n_insert - len_d_insert
                        for num_n_before_d in range(num_n+1):
                            num_n_after_d = num_n - num_n_before_d
                            assert num_n_after_d>=0

                            n_nucseq = ( v_nucseq[:len_v_nucseq-v_trim] +
                                         'n'*num_n_before_d + d_insert + 'n'*num_n_after_d +
                                         j_nucseq[j_trim:] )

                            assert len(n_nucseq) == 3*len_cdr3_protseq

                            trim_prob = tcr_rearrangement.get_beta_trim_probs( organism, did,
                                                                               v_trim, d0_trim, d1_trim, j_trim,
                                                                               num_n_before_d, num_n_after_d )
                            if not trim_prob: continue

                            if nucleotide_match:
                                assert len(n_nucseq) == len_cdr3_nucseq
                                matched = True
                                #print n_nucseq, cdr3_nucseq
                                for a,b in zip( n_nucseq, cdr3_nucseq ):
                                    if a!=b and a!= 'n':
                                        matched=False
                                if matched:
                                    coding_prob = 0.25 ** num_n
                                else:
                                    coding_prob = 0.0
                            else:
                                coding_prob = get_coding_probability( n_nucseq, cdr3_protseq )
                            prob = did_prob * coding_prob * trim_prob

                            total_prob_this_insert += prob ## just for status output
                            total_prob_this_trim += prob
                            total_prob += prob

                            if verbose and coding_prob:
                                print 'coding_prob:',cdr3_protseq,"trims:",v_trim,d0_trim,d1_trim,j_trim,\
                                    "inserts:",num_n_before_d,num_n_after_d,\
                                    "d_insert:",d_insert,\
                                    "total_prob:",total_prob,"prob:",prob,"coding_prob:",coding_prob,\
                                    "trim_prob:",trim_prob,n_nucseq


            if verbose:
                print 'n_insert:',n_insert,extra_v_trim,extra_j_trim,'total_prob:',total_prob,\
                    'total_prob_this_insert:',total_prob_this_insert


        if extra_trim>2 and total_prob_this_trim < error_threshold * old_total_prob:
            break

    if return_final_cdr3_nucseq:
        return total_prob, cdr3_nucseq
    else:
        return total_prob

def setup_random_sampling_list( probs ):
    #print 'setup_random_sampling_list:',probs
    total_prob = 0.0
    l = []
    for k,p in probs.iteritems():
        total_prob += p
        l.append( ( total_prob,k ) )
    assert abs( 1-total_prob )<1e-3 ## since we normalized already...
    return l

def sample_from_random_sampling_list( l ):
    #print 'sample_from_random_sampling_list:',l
    f = random.random()
    for (prob,k) in l:
        if f<=prob:
            return k
    return l[-1][1]

def sample_alpha_sequences( organism, nsamples, v_gene, j_gene, force_aa_length = 0,
                            in_frame_only = True, no_stop_codons = True,
                            max_tries = 100000000,
                            include_annotation= False ):
    ab = 'A'
    #organism = 'mouse'
    bases = 'acgt'

    trim_probs = tcr_rearrangement.all_trim_probs[ organism ]

    v_nucseq = get_v_cdr3_nucseq( organism, v_gene )
    j_nucseq = get_j_cdr3_nucseq( organism, j_gene )

    v_nucseq_len = len( v_nucseq )
    j_nucseq_len = len( j_nucseq )

    max_v_trim = min( 15, v_nucseq_len -3 )
    max_j_trim = min( 15, j_nucseq_len -3 )

    seqs = []

    nsampled = 0

    v_trim_probsl = setup_random_sampling_list( trim_probs[ 'A_v_trim' ] )
    j_trim_probsl = setup_random_sampling_list( trim_probs[ 'A_j_trim' ] )
    vj_insert_probsl = setup_random_sampling_list( trim_probs[ 'A_vj_insert' ] )

    ntries=0
    while nsampled < nsamples:
        ntries += 1
        if ntries%100000==0:
            Log('sample_alpha_sequences: ntries: {} {} {} {} force_aa_length {}'\
                .format(ntries,nsampled,v_gene,j_gene,force_aa_length))
        if ntries>max_tries:
            break

        vtrim = sample_from_random_sampling_list( v_trim_probsl )
        jtrim = sample_from_random_sampling_list( j_trim_probsl )
        n_insert = sample_from_random_sampling_list( vj_insert_probsl )

        if vtrim > max_v_trim or jtrim > max_j_trim: continue

        ## check if in frame?
        cdr3_len = v_nucseq_len + j_nucseq_len + n_insert - ( vtrim + jtrim )

        if in_frame_only and (cdr3_len%3 != 0): continue
        if force_aa_length and force_aa_length != cdr3_len/3: continue

        vj_insert = ''
        for i in range(n_insert): vj_insert += random.choice( bases )

        cdr3_seq = v_nucseq[:v_nucseq_len-vtrim] + vj_insert + j_nucseq[jtrim:]

        assert len(cdr3_seq) == cdr3_len

        ## check for stop codons?
        protseq = ''
        for i in range(cdr3_len/3):
            protseq += genetic_code[ cdr3_seq[ 3*i : 3*i+3 ] ]
        if '*' in protseq and no_stop_codons:
            continue

        if include_annotation:
            cdr3_annotation = ( 'V'*(v_nucseq_len-vtrim) +
                                'N'*n_insert +
                                'J'*(j_nucseq_len-jtrim) )
            assert len(cdr3_seq) == len(cdr3_annotation)
            seqs.append( ( cdr3_seq, protseq, cdr3_annotation ) )
        else:
            seqs.append( ( cdr3_seq, protseq ) )
        nsampled += 1
    return seqs


def sample_beta_sequences( organism, nsamples, v_gene, j_gene, force_aa_length = 0,
                           in_frame_only = True,
                           no_stop_codons = True,
                           max_dj_insert = 10,
                           max_tries = 100000000,
                           include_annotation = False ):
    trim_probs = tcr_rearrangement.all_trim_probs[ organism ]

    #organism = 'mouse'
    bases = 'acgt'

    v_nucseq = get_v_cdr3_nucseq( organism, v_gene )
    j_nucseq = get_j_cdr3_nucseq( organism, j_gene )
    v_nucseq_len = len( v_nucseq )
    j_nucseq_len = len( j_nucseq )

    max_v_trim = min( 15, v_nucseq_len -3 )
    max_j_trim = min( 15, j_nucseq_len -3 )

    v_trim_probsl = setup_random_sampling_list( trim_probs[ 'B_v_trim' ] )
    j_trim_probsl = setup_random_sampling_list( trim_probs[ 'B_j_trim' ] )
    vd_insert_probsl = setup_random_sampling_list( trim_probs[ 'B_vd_insert' ] )
    dj_insert_probsl = setup_random_sampling_list( trim_probs[ 'B_dj_insert' ] )

    dids = tcr_rearrangement.all_trbd_nucseq[organism].keys()

    d_trim_probsl = dict( zip( dids, [ setup_random_sampling_list( trim_probs['B_D{}_d01_trim'.format(x)] )
                                       for x in dids] ) )
    # d_trim_probsl = { 1: setup_random_sampling_list( trim_probs['B_D1_d01_trim'] ),
    #                   2: setup_random_sampling_list( trim_probs['B_D2_d01_trim'] ) }

    jno = 0 # no D filtering
    if organism in ['human','mouse'] and j_gene[2] == 'B':
        jno = int( j_gene[4] )
        assert jno in [1,2]
    if jno == 1:
        allowed_dgenes = [1]
    else:
        allowed_dgenes = dids[:]

    seqs = []
    nsampled = 0
    ntries=0

    while nsampled < nsamples:
        ntries += 1
        if ntries%100000==0:
            Log('sample_beta_sequences: ntries: {} {} {} {} force_aa_length {}'\
                .format(ntries,nsampled,v_gene,j_gene,force_aa_length))

        if ntries>max_tries:
            break
        ## pick d segment
        dgene = random.choice( allowed_dgenes )
        d_nucseq = tcr_rearrangement.all_trbd_nucseq[organism][dgene]
        d_nucseq_len = len(d_nucseq)

        vtrim = sample_from_random_sampling_list( v_trim_probsl )
        jtrim = sample_from_random_sampling_list( j_trim_probsl )
        n_vd_insert = sample_from_random_sampling_list( vd_insert_probsl )
        n_dj_insert = sample_from_random_sampling_list( dj_insert_probsl )

        n_d0_trim,n_d1_trim = sample_from_random_sampling_list( d_trim_probsl[dgene] )

        if n_dj_insert > max_dj_insert: continue ## some of these are probably bogus (sequencing errors)
        if vtrim > max_v_trim or jtrim > max_j_trim or (n_d0_trim + n_d1_trim) > d_nucseq_len: continue

        ## check if in frame?
        cdr3_len = v_nucseq_len + j_nucseq_len + d_nucseq_len + n_vd_insert + n_dj_insert - \
                   ( vtrim+jtrim+n_d0_trim+n_d1_trim )

        if in_frame_only and (cdr3_len%3 != 0): continue
        if force_aa_length and force_aa_length != cdr3_len/3: continue

        vd_insert = ''
        for i in range(n_vd_insert): vd_insert += random.choice( bases )
        dj_insert = ''
        for i in range(n_dj_insert): dj_insert += random.choice( bases )

        ## choose inserted bases
        d_seq = d_nucseq[ n_d0_trim : d_nucseq_len-n_d1_trim]

        cdr3_seq = v_nucseq[:v_nucseq_len-vtrim] + vd_insert + d_seq + dj_insert + j_nucseq[jtrim:]


        assert len(cdr3_seq) == cdr3_len

        ## check for stop codons?
        protseq = ''
        for i in range(cdr3_len/3):
            protseq += genetic_code[ cdr3_seq[ 3*i : 3*i+3 ] ]
        if '*' in protseq and no_stop_codons:
            continue

        if include_annotation:
            cdr3_annotation = ( 'V'*(v_nucseq_len-vtrim) +
                                'N'*n_vd_insert +
                                'D'*len(d_seq)+
                                'N'*n_dj_insert+
                                'J'*(j_nucseq_len-jtrim) )
            assert len(cdr3_seq) == len(cdr3_annotation)
            seqs.append( ( cdr3_seq, protseq, cdr3_annotation ) )
        else:
            seqs.append( ( cdr3_seq, protseq ) )
        nsampled += 1
    return seqs

def sample_tcr_sequences( organism, nsamples, v_gene, j_gene,
                          force_aa_length = 0,
                          in_frame_only = True,
                          no_stop_codons = True,
                          max_tries = 100000000,
                          include_annotation = False ):
    ab = all_genes[organism][v_gene].chain
    assert ab in 'AB'
    if ab == 'A':
        return sample_alpha_sequences( organism, nsamples, v_gene, j_gene, force_aa_length = force_aa_length,
                                       in_frame_only = in_frame_only, no_stop_codons = no_stop_codons,
                                       max_tries = max_tries, include_annotation = include_annotation )
    else:
        return sample_beta_sequences( organism, nsamples, v_gene, j_gene, force_aa_length = force_aa_length,
                                      in_frame_only = in_frame_only, no_stop_codons = no_stop_codons,
                                      max_tries = max_tries, include_annotation = include_annotation )

def add_masked_CDR3_sequences_to_tcr_dict( organism, vals ):
    ## this code is mostly taken from compute_probs.py
    va_gene = vals['va_gene']
    ja_gene = vals['ja_gene']
    vb_gene = vals['vb_gene']
    jb_gene = vals['jb_gene']
    cdr3a_protseq = vals['cdr3a']
    cdr3a_nucseq  = vals['cdr3a_nucseq']
    cdr3b_protseq = vals['cdr3b']
    cdr3b_nucseq  = vals['cdr3b_nucseq']

    a_junction_results = analyze_junction( organism, va_gene, ja_gene, cdr3a_protseq, cdr3a_nucseq )
    b_junction_results = analyze_junction( organism, vb_gene, jb_gene, cdr3b_protseq, cdr3b_nucseq )

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


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

