from basic import *
from amino_acids import amino_acids
from tcr_distances_blosum import blosum
from paths import path_to_db
import translation

cdrs_sep = ';'
gap_character = '.'

all_genes = {}

class TCR_Gene:
    def __init__( self, l ):
        self.id = l['id']
        self.organism = l['organism']
        self.chain = l['chain']
        self.region = l['region']
        self.nucseq = l['nucseq']
        self.alseq = l['aligned_protseq']
        self.cdrs = l['cdrs'].split(cdrs_sep) if l['cdrs'] else []
        ## these are still 1-indexed !!!!!!!!!!!!!!
        self.cdr_columns = [ map( int,x.split('-')) for x in l['cdr_columns'].split(cdrs_sep) ] if self.cdrs else []
        frame = l['frame']
        assert frame in ['+1','+2','+3','1','2','3']
        self.nucseq_offset = int( frame[-1] )-1 ## 0, 1 or 2 (0-indexed for python)
        self.protseq = translation.get_translation( self.nucseq, frame )[0]
        assert self.protseq == self.alseq.replace(gap_character,'')
        # sanity check
        if self.cdrs:
            #print self.cdrs
            #print [ self.alseq[ x[0]-1 : x[1] ] for x in self.cdr_columns ]
            assert self.cdrs == [ self.alseq[ x[0]-1 : x[1] ] for x in self.cdr_columns ]

## need to make this a little more configurable (cmdline??)

db_file = path_to_db+'/'+pipeline_params['db_file']
assert exists(db_file)

lines = parse_tsv_file( db_file )

for l in lines:
    g = TCR_Gene( l )
    if g.organism not in all_genes:
        all_genes[g.organism] = {} # map from id to TCR_Gene objects
    all_genes[g.organism][g.id] = g


verbose = ( __name__ == '__main__' )

for organism,genes in all_genes.iteritems():

    for ab in 'AB':
        org_merged_loopseqs = {}
        for id,g in genes.iteritems():
            if g.chain == ab and g.region == 'V':
                loopseqs = g.cdrs[:-1] ## exclude CDR3 Nterm
                org_merged_loopseqs[id] = ' '.join( loopseqs )

        all_loopseq_nbrs = {}
        all_loopseq_nbrs_mm1 = {}
        for id1,seq1 in org_merged_loopseqs.iteritems():
            g1 = genes[id1]
            cpos = g1.cdr_columns[-1][0] - 1 #0-indexed
            alseq1 = g1.alseq
            minlen = cpos+1
            assert len(alseq1) >= minlen
            if alseq1[cpos] != 'C':
                print 'funny cpos:',id1,alseq1,g1.cdrs[-1]

            all_loopseq_nbrs[id1] = []
            all_loopseq_nbrs_mm1[id1] = []
            for id2,seq2 in org_merged_loopseqs.iteritems():
                g2 = genes[id2]
                alseq2 = g2.alseq
                assert len(alseq2) >= minlen
                assert len(seq1) == len(seq2)
                if seq1 == seq2:
                    all_loopseq_nbrs[id1].append( id2 )
                    all_loopseq_nbrs_mm1[id1].append( id2 )
                    continue

                ## count mismatches between these two, maybe count as an "_mm1" nbr
                loop_mismatches = 0
                loop_mismatches_cdrx = 0
                loop_mismatch_seqs =[]
                spaces=0
                for a,b in zip( seq1,seq2):
                    if a==' ':
                        spaces+=1
                        continue
                    if a!= b:
                        if a in '*.' or b in '*.':
                            loop_mismatches += 10
                            break
                        else:
                            assert a in amino_acids and b in amino_acids
                            if spaces<=1:
                                loop_mismatches += 1
                                loop_mismatch_seqs.append( ( a,b ) )
                            else:
                                assert spaces==2
                                loop_mismatches_cdrx += 1
                            if loop_mismatches>1:
                                break
                if loop_mismatches <=1:
                    all_mismatches = 0
                    for a,b in zip( alseq1[:cpos+2],alseq2[:cpos+2]):
                        if a!= b:
                            if a in '*.' or b in '*.':
                                all_mismatches += 10
                            else:
                                assert a in amino_acids and b in amino_acids
                                all_mismatches += 1
                    #dist = tcr_distances.blosum_sequence_distance( seq1, seq2, gap_penalty=10 )
                    if loop_mismatches<=1 and loop_mismatches + loop_mismatches_cdrx <= 2 and all_mismatches<=10:
                        if loop_mismatches == 1:
                            blscore= blosum[(loop_mismatch_seqs[0][0],loop_mismatch_seqs[0][1])]
                        else:
                            blscore = 100
                        if blscore>=1:
                            all_loopseq_nbrs_mm1[id1].append( id2 )
                            if loop_mismatches>0 and verbose:
                                mmstring = ','.join(['%s/%s'%(x[0],x[1]) for x in loop_mismatch_seqs])
                                gene1 = id1[:id1.index('*')]
                                gene2 = id2[:id2.index('*')]
                                if gene1 != gene2:
                                    print 'v_mismatches:',organism,mmstring,blscore,id1,id2,\
                                        loop_mismatches,loop_mismatches_cdrx,all_mismatches,seq1
                                    print 'v_mismatches:',organism,mmstring,blscore,id1,id2,\
                                        loop_mismatches,loop_mismatches_cdrx,all_mismatches,seq2


        for id in all_loopseq_nbrs:
            rep = min( all_loopseq_nbrs[id] )
            assert org_merged_loopseqs[id] == org_merged_loopseqs[ rep ]
            genes[id].rep = rep
            if verbose:
                print 'vrep %s %15s %15s %s'%(organism, id, rep, org_merged_loopseqs[id])


        ## merge mm1 nbrs to guarantee transitivity
        while True:
            new_nbrs = False
            for id1 in all_loopseq_nbrs_mm1:
                new_id1_nbrs = False
                for id2 in all_loopseq_nbrs_mm1[id1]:
                    for id3 in all_loopseq_nbrs_mm1[id2]:
                        if id3 not in all_loopseq_nbrs_mm1[id1]:
                            all_loopseq_nbrs_mm1[id1].append( id3 )
                            if verbose:
                                print 'new_nbr:',id1,'<--->',id2,'<--->',id3
                            new_id1_nbrs = True
                            break
                    if new_id1_nbrs:
                        break
                if new_id1_nbrs:
                    new_nbrs = True
            if verbose:
                print 'new_nbrs:',ab,organism,new_nbrs
            if not new_nbrs:
                break

        for id in all_loopseq_nbrs_mm1:
            rep = min( all_loopseq_nbrs_mm1[id] )
            genes[id].mm1_rep = rep
            if verbose:
                print 'mm1vrep %s %15s %15s %s'%(organism, id, rep,org_merged_loopseqs[id])


    ## setup Jseq reps
    for ab in 'AB':
        jloopseqs = {}
        for id,g in genes.iteritems():
            if g.chain == ab and g.region == 'J':
                num = len( g.cdrs[0].replace( gap_character, '' ) )
                jloopseq = g.protseq[:num+3] ## go all the way up to and including the GXG
                jloopseqs[id] = jloopseq
        all_jloopseq_nbrs = {}
        for id1,seq1 in jloopseqs.iteritems():
            all_jloopseq_nbrs[id1] = []
            for id2,seq2 in jloopseqs.iteritems():
                if seq1 == seq2:
                    all_jloopseq_nbrs[id1].append( id2 )
        for id in all_jloopseq_nbrs:
            rep = min( all_jloopseq_nbrs[id] )
            genes[id].rep = rep
            genes[id].mm1_rep = rep # just so we have an mm1_rep field defined...
            assert jloopseqs[id] == jloopseqs[ rep ]
            if verbose:
                print 'jrep %s %15s %15s %15s'%(organism, id, rep, jloopseqs[id])



def get_cdr3_and_j_match_counts( organism, ab, qseq, j_gene, min_min_j_matchlen = 3,
                                 extended_cdr3 = False ):
    #fasta = all_fasta[organism]
    jg = all_genes[organism][j_gene]

    errors = []

    ## qseq starts at CA...
    assert qseq[0] == 'C'

    num_genome_j_positions_in_loop = len(jg.cdrs[0].replace(gap_character,''))-2
    #num_genome_j_positions_in_loop = all_num_genome_j_positions_in_loop[organism][ab][j_gene]

    if extended_cdr3: num_genome_j_positions_in_loop += 2 ## up to but not including GXG

    ## history: was only for alpha
    aseq = qseq[:] ## starts at the C position

    ja_gene = j_gene
    assert ja_gene in fasta
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
        looplen = pos - ntrim + num_genome_j_positions_in_loop
        if not extended_cdr3:
            aseq = aseq[3:]
            looplen -= 3 ## dont count CAX
        if len(aseq)<looplen:
            Log(`( 'short',ab,aseq,ja_seq )`)
            errors.append( ab+'seq_too_short' )
            return '-',[100,0],errors

        cdrseq = aseq[:looplen ]

    ## now count mismatches in the J gene, beyond the cdrseq
    j_seq = jg.protseq #fasta[ j_gene ] ## not sure why we do this again (old legacy code)
    if qseq.count( cdrseq ) > 1:
        Log('multiple cdrseq occurrences %s %s'%(qseq,cdrseq))
        errors.append('multiple_cdrseq_occ')
        return '-',[100,0],errors
    assert qseq.count(cdrseq) == 1
    start_counting_qseq = qseq.find(cdrseq)+len(cdrseq)
    start_counting_jseq = num_genome_j_positions_in_loop
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




def parse_cdr3( organism, ab, qseq, v_gene, j_gene, q2v_align, extended_cdr3 = False ):
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

    cdrseq, j_match_counts, other_errors = get_cdr3_and_j_match_counts( organism, ab, qseq[ cpos_match: ], j_gene,
                                                                        extended_cdr3 = extended_cdr3 )

    return cdrseq, v_match_counts, j_match_counts, errors+other_errors


if __name__ == '__main__':
    ## show J alignments
    pass
