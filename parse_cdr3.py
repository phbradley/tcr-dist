from basic import *
from all_genes import all_genes, gap_character


def get_cdr3_and_j_match_counts( organism, ab, qseq, j_gene, min_min_j_matchlen = 3,
                                 extended_cdr3 = False,
                                 max_missing_aas_at_cdr3_cterm = 2 ): # new (was 0)
    #fasta = all_fasta[organism]
    jg = all_genes[organism][j_gene]

    errors = []

    ## qseq starts at CA...
    assert qseq[0] == 'C'

    num_genome_j_aas_in_loop = len(jg.cdrs[0].replace(gap_character,''))-2
    if extended_cdr3: num_genome_j_aas_in_loop += 2 ## up to but not including GXG

    ## history: was only for alpha
    aseq = qseq[:] ## starts at the C position

    ja_gene = j_gene
    #assert ja_gene in fasta
    ja_seq = jg.protseq #fasta[ ja_gene ]

    min_j_matchlen = min_min_j_matchlen+3

    while min_j_matchlen >= min_min_j_matchlen:
        ntrim =0
        while ntrim+min_j_matchlen<len(ja_seq) and ja_seq[ntrim:ntrim+min_j_matchlen] not in aseq:
            ntrim += 1

        jatag = ja_seq[ntrim:ntrim+min_j_matchlen]
        if jatag in aseq:
            break
        else:
            min_j_matchlen -= 1

    #print 'min_j_matchlen:',min_j_matchlen,'jatag:',jatag,'ntrim:',ntrim,'ja_seq:',ja_seq,'qseq',qseq

    if jatag not in aseq:
        Log(`( 'whoah',ab,aseq,ja_seq )`)
        errors.append( 'j{}tag_not_in_aseq'.format(ab) )
        return '-',[100,0],errors
    elif ja_seq.count( jatag ) != 1:
        Log(`( 'whoah2',ab,aseq,ja_seq )`)
        errors.append( 'multiple_j{}tag_in_jseq'.format(ab) )
        return '-',[100,0],errors
    else:
        pos = aseq.find( jatag )
        looplen = pos - ntrim + num_genome_j_aas_in_loop
        if not extended_cdr3:
            aseq = aseq[3:]
            looplen -= 3 ## dont count CAX
        if len(aseq)<looplen:
            num_missing = looplen - len(aseq )
            if num_missing > max_missing_aas_at_cdr3_cterm:
                Log(`( 'short',ab,aseq,ja_seq )`)
                errors.append( ab+'seq_too_short' )
                return '-',[100,0],errors ## early return
            suffix = ja_seq[num_genome_j_aas_in_loop-num_missing:num_genome_j_aas_in_loop]
            ## NOTE -- changin qseq, aseq
            print 'max_missing_aas_at_cdr3_cterm:',max_missing_aas_at_cdr3_cterm,'num_missing:',num_missing
            aseq += suffix
            qseq += suffix
        cdrseq = aseq[:looplen ]

    ## now count mismatches in the J gene, beyond the cdrseq
    j_seq = jg.protseq #fasta[ j_gene ] ## not sure why we do this again (old legacy code)
    if qseq.count( cdrseq ) > 1:
        Log('multiple cdrseq occurrences %s %s'%(qseq,cdrseq))
        errors.append('multiple_cdrseq_occ')
        return '-',[100,0],errors
    assert qseq.count(cdrseq) == 1
    start_counting_qseq = qseq.find(cdrseq)+len(cdrseq)
    start_counting_jseq = num_genome_j_aas_in_loop
    j_match_counts = [0,0]
    #assert extended_cdr3 ## otherwise I think this count is not right?
    #print 'here',start_counting_qseq,start_counting_jseq,len(qseq)
    for qpos in range( start_counting_qseq, len(qseq)):
        jpos = start_counting_jseq + (qpos-start_counting_qseq)
        #print 'here',qpos,jpos
        if jpos>= len(j_seq): break
        if qseq[qpos] == j_seq[jpos]:
            j_match_counts[1] += 1
        else:
            j_match_counts[0] += 1

    return cdrseq, j_match_counts,errors




def parse_cdr3( organism, ab, qseq, v_gene, j_gene, q2v_align,
                extended_cdr3 = False, max_missing_aas_at_cdr3_cterm = 2 ):
    ## v_align is a mapping from 0-indexed qseq positions to 0-indexed v_gene protseq positions
    #fasta = all_fasta[ organism ]
    #align_fasta = all_align_fasta[ organism ]
    vg = all_genes[organism][v_gene]

    errors = []

    ## what is the C position in this v gene?
    v_seq = vg.protseq #fasta[ v_gene ]
    v_alseq = vg.alseq #align_fasta[ v_gene ]
    assert v_seq == v_alseq.replace(gap_character,'')

    alseq_cpos = vg.cdr_columns[-1][0] - 1 ## now 0-indexed
    #alseq_cpos = alseq_C_pos[organism][ab] - 1 ## now 0-indexed
    numgaps = v_alseq[:alseq_cpos].count(gap_character)

    cpos = alseq_cpos - numgaps ## 0-indexed
    cpos_match = -1

    v_match_counts = [0,0]

    qseq_len = len(qseq)
    for (qpos,vpos) in sorted( q2v_align.iteritems() ):
        #print 'q2v-align:',qpos, vpos, cpos
        if qpos == len(qseq):
            continue ## from a partial codon at the end
        if vpos == cpos:
            cpos_match = qpos
        elif vpos <= cpos:
            ## only count v mismatches here
            if qseq[qpos] == v_seq[vpos]:
                v_match_counts[1] += 1
            else:
                v_match_counts[0] += 1

    if cpos_match<0 or qseq[ cpos_match ] != 'C':
        ## problemo
        Log('failed to find blast match to C position')
        errors.append('no_V{}_Cpos_blastmatch'.format(ab))
        return '-',[100,0],[100,0],errors

    cdrseq, j_match_counts, other_errors \
        = get_cdr3_and_j_match_counts( organism, ab, qseq[ cpos_match: ], j_gene,
                                       extended_cdr3 = extended_cdr3,
                                       max_missing_aas_at_cdr3_cterm = max_missing_aas_at_cdr3_cterm )

    return cdrseq, v_match_counts, j_match_counts, errors+other_errors

