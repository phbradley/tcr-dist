from basic import *
from amino_acids import amino_acids
from tcr_distances_blosum import blosum
from paths import path_to_db

## these are indexed by organism
## these are potentially used outside this python file #########################
pb_cdrs = {}
all_align_fasta = {}
all_fasta = {}
all_num_genome_j_positions_in_loop = {}
all_loopseq_representative = {}
all_loopseq_representative_mm1 = {}
all_jseq_representative = {}
all_core_positions = {}
all_merged_loopseqs = {}

gap_character = '.'

#################################################################################


fasta_dir = path_to_db+'/fasta/'

verbose = ( __name__ == '__main__' )

## these go out of date when we remove gaps from the alignment, below
##
## look at imgt cdrs
## http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
## imgt_cdr_positions = [ ( 27, 38 ), ( 56, 65 ) ] ## 1-indexed
## these are 1-indexed !!!
##
## note that the TRAV mouse alignment seems to be shifted by 1 relative to IMGT for FR1 and FR2 and by 2 for C104 (->106)
## looks like the insertion happens around TRAV-alpha alignment position 86
## the other three agree at anchor positions C23, W41, C104
##
pb_cdr_positions   = { 'mouse': {'A': [ ( 28, 39 ), ( 57, 66 ), (82, 88) ],
                                 'B': [ ( 27, 38 ), ( 56, 65 ), (81, 86) ] },
                       'human': {'A': [ ( 27, 38 ), ( 56, 65 ), (81, 86) ],
                                 'B': [ ( 27, 38 ), ( 56, 65 ), (81, 86) ] },
                       }


alseq_C_pos = { 'mouse':{'A':106,'B':104},
                'human':{'A':104,'B':104} } ## 1-indexed

## 1-indexed:
extra_alignment_columns = { 'mouse':{'A':[9,86],'B':[] }, ## 1-indexed
                            'human':{'A':[],'B':[] } }


core_positions_generic_1indexed = [
    21, 23, 25,   ##  23 is C
    39, 41,       ##  41 is W
    53, 54, 55,
    78,           ##            maybe also 80?
    89,           ##  89 is L
    102, 103, 104 ## 104 is C
]


default_num_positions_after_GXG = { 'mouse': {'A':7, 'B':6 },
                                    'human': {'A':7, 'B':6 }}


## I did this by hand for the initial db version. Need to come up with a better way in the future...
##
all_funny_jseq = """
mouse TRAJ19*01               IYRGFHKFSSGIESKHNVSP KFSSGIE
mouse TRAJ20*01                 SGNYKLGVESVTMMSVRA KLGVESV
mouse TRAJ25*01                 RTKVSSVFGTWRRLLVKP SSVFGTW
mouse TRAJ29*01                NSGSRELVLGREARLSMIE ELVLGRE
mouse TRAJ3*01               EFSYSSKLIFGAETKLRNPPY KLIFGAE
mouse TRAJ35*01              QTGFASALTFGSGTKVIPCLP ALTFGSG
mouse TRAJ41*01                 VSNTSSMLAEAPHYWSHP SSMLAEA
mouse TRAJ45*02                 NTGGADRLTFGKGTQLII RLTFGKG
mouse TRAJ49*01               NTGYQNFYFGKGTSLTVIPS NFYFGKG
mouse TRAJ59*01               LLKREDKATFATGGYEAEED KATFATG
mouse TRBJ1-6*01                 SYNSPLYFAAGTRLTVT PLYFAAG
mouse TRBJ2-6*01                   ALALTDWQPIEQPMR ALTDWQP
human TRAJ16*01                FSDGQKLLFARGTMLKVDL KLLFARG
human TRAJ59*01                  KEGNRKFTFGMGTQVRV KFTFGMG
human TRAJ61*01                YRVNRKLTFGANTRGIMKL KLTFGAN
human TRBJ2-2P*01                  LRGAAGRLGGGLLVL GAAGRLG
"""
funny_jseq = {}
for line in all_funny_jseq.split('\n'):
    l = line.split()
    if len(l) == 4:
        organism = l[0]
        ab = l[1][2]
        if organism not in funny_jseq: funny_jseq[organism] = {}
        if ab not in funny_jseq[organism]:funny_jseq[organism][ab] = []
        funny_jseq[organism][ab].append( l[3] )



for organism in [ 'mouse','human' ]:

    ## read the TR-V alignments
    align_file = '{}imgt_{}_TR_protein_sequences_with_gaps.fasta'.format( fasta_dir, organism )
    assert exists(align_file)
    align_fasta = {}
    tr_prefixes = ['TRBV','TRAV']
    for line in open( align_file,'r'):
        if line[0] == '>':
            id = line.split('|')[1]
            if id[:4] in tr_prefixes: assert id not in align_fasta
            align_fasta[id] = ''
        else:
            align_fasta[id] += line.split()[0]

    for id in align_fasta:
        if id[3] == 'V' and id[:4] in tr_prefixes:
            cpos = alseq_C_pos[organism][id[2]]
            if len(align_fasta[id])<cpos:
                if verbose:
                    print 'short alseq:',organism,id,cpos,len(align_fasta[id])
                align_fasta[id] += gap_character*(cpos - len(align_fasta[id]))

            if align_fasta[id][cpos-1] != 'C' and verbose:
                print 'bad cpos',id, align_fasta[id][cpos-1]

    fastafile = '{}imgt_{}_TR_protein_sequences.fasta'.format( fasta_dir, organism )
    assert exists(fastafile)

    ## read the fasta file
    fasta = {}
    tr_prefixes = ['TRBV','TRBJ','TRAV','TRAJ']
    for line in open( fastafile,'r'):
        if line[0] == '>':
            id = line.split('|')[1]
            if id[:4] in tr_prefixes: assert id not in fasta
            fasta[id] = ''
        else:
            fasta[id] += line.split()[0]

    ##

    if True: ## setup num_genome_j_positions_in_loop
        all_num_genome_j_positions_in_loop[organism] = {}
        for ab in 'AB':
            all_num_genome_j_positions_in_loop[organism][ab] = {}
            ids = [ x for x in align_fasta if x[2] == ab and x[3] == 'J' ] ## TRxV ids
            ids.sort()
            L = max( len(align_fasta[x]) for x in ids )

            default_suffixlen = default_num_positions_after_GXG[ organism ][ab]
            for id in ids:
                suffixlen = default_suffixlen
                alseq = align_fasta[id]
                if '*' in alseq:
                    all_num_genome_j_positions_in_loop[organism][ab][id] = len(alseq) - suffixlen - 5
                    continue
                if not ( alseq[-suffixlen-1] == 'G' and alseq[-suffixlen-3] == 'G'):
                    suffixlen = 0
                    for word in funny_jseq[organism][ab]:
                        if word in alseq:
                            assert not suffixlen
                            suffixlen = len(alseq) - alseq.find(word) - len(word)
                assert suffixlen
                starred = False
                num_spaces = L - len(alseq)
                if suffixlen != default_suffixlen or alseq[-suffixlen-1] != 'G' or alseq[-suffixlen-3] != 'G':
                    starred = True
                    num_spaces += ( suffixlen - default_suffixlen )
                ## there are 2 residues between the loop and the GXG
                all_num_genome_j_positions_in_loop[organism][ab][id] = len(alseq) - suffixlen - 5
                if verbose:
                    print '%s %-20s %s%s   %s'%(organism,id,' '*num_spaces,alseq,'**********'*starred)

    if True: ## show alignment columns in plain text ################################### and debugging
        col_len = 15000000
        for ab in 'AB':
            cpos = alseq_C_pos[organism][ab] - 1 ## 0-indexed
            L = cpos+2
            ids = [ x for x in align_fasta if x[2] == ab and x[3] == 'V' and len(align_fasta[x])>=L ] ## TRxV ids
            ids.sort()

            for pos in range(L):
                cdrtag = '|'
                for start,stop in pb_cdr_positions[organism][ab]:
                    if pos+1>=start and pos+1<=stop:
                        cdrtag = '+'
                col = ''.join( [align_fasta[x][pos] for x in ids[:col_len] ] )
                if pos == cpos:
                    assert col.count('C') > (95*len(col))/100
                if verbose:
                    print '%s TR%sV %s %4d %s'%(organism,ab,cdrtag,pos+1,col)




    pb_cdrs[organism] = {}
    for id,alseq in align_fasta.iteritems():
        if id.startswith('TRBV') or id.startswith('TRAV'):
            pb_cdrs[organism][id] = []
            # cdr1 = align_fasta[id][pb_cdr_positions[0][0]-1 : pb_cdr_positions[0][1] ]
            # cdr2 = align_fasta[id][pb_cdr_positions[1][0]-1 : pb_cdr_positions[1][1] ]
            # cdrX = align_fasta[id][pb_cdr_positions[2][0]-1 : pb_cdr_positions[2][1] ]
            #print '%-15s %s %s %s'%(id,cdr1,cdr2,cdrX)
            for start,stop in pb_cdr_positions[organism][id[2]]:
                assert stop<=len(alseq)
                pos1 = start-1
                pos2 = stop-1
                alseq_loop = alseq[pos1:pos2+1]
                numgaps = alseq[:pos1].count(gap_character)
                loop_start = pos1-numgaps
                loop_stop = loop_start + len(alseq_loop) - 1 - alseq_loop.count(gap_character)
                a = ''.join( alseq_loop.split(gap_character) )
                b = fasta[id][loop_start:loop_stop+1]
                #print id, alseq_loop, b, loop_start, loop_stop#, cdr1, cdr2, cdrX
                pb_cdrs[organism][id].append( ( b, loop_start, loop_stop ) ) ## 0-indexed start and stop, inclusive
                assert a == b


    ## setup reps
    all_loopseq_representative[organism] = {}
    all_loopseq_representative_mm1[organism] = {}
    all_merged_loopseqs[ organism ] = {}

    for ab in 'AB':
        org_merged_loopseqs = {}
        for id,alseq in align_fasta.iteritems():
            if id[2] == ab and id[3] == 'V':
                loopseqs = []
                for start,stop in pb_cdr_positions[organism][ab]: ## start,stop are 1-indexed
                    pos1 = start-1
                    pos2 = stop-1
                    alseq_loop = alseq[pos1:pos2+1]
                    loopseqs.append( alseq_loop )
                org_merged_loopseqs[id] = ' '.join( loopseqs )
        all_loopseq_nbrs = {}
        all_loopseq_nbrs_mm1 = {}
        for id1,seq1 in org_merged_loopseqs.iteritems():
            cpos = alseq_C_pos[organism][ab] - 1 ## 0-indexed
            alseq1 = align_fasta[id1]
            minlen = cpos+1
            if len(alseq1)<minlen:
                alseq1 = align_fasta[id1] + 'X'*( minlen-len(alseq1))
                if verbose:
                    print 'short_align:',id1,len(alseq1),minlen,alseq1

            all_loopseq_nbrs[id1] = []
            all_loopseq_nbrs_mm1[id1] = []
            for id2,seq2 in org_merged_loopseqs.iteritems():
                alseq2 = align_fasta[id2]
                if len(alseq2)<minlen:
                    alseq2 = align_fasta[id2] + 'X'*( minlen-len(alseq2))
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
            all_loopseq_representative[organism][id] = min( all_loopseq_nbrs[id] )
            assert org_merged_loopseqs[id] == org_merged_loopseqs[ all_loopseq_representative[organism][id] ]
            if verbose:
                print 'vrep %s %15s %15s %s'%(organism, id, all_loopseq_representative[organism][id],
                                              org_merged_loopseqs[id])
            all_merged_loopseqs[ organism ][ id ] = org_merged_loopseqs[id][:]


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
                                print 'new_nbr:',id1,id2,id3
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
            all_loopseq_representative_mm1[organism][id] = min( all_loopseq_nbrs_mm1[id] )
            if verbose:
                print 'mm1vrep %s %15s %15s %s'%(organism, id, all_loopseq_representative_mm1[organism][id],
                                                 org_merged_loopseqs[id])


    ## setup Jseq reps
    all_jseq_representative[organism] = {}

    for ab in 'AB':
        jloopseqs = {}
        for id,jseq in fasta.iteritems():
            if id[2] == ab and id[3] == 'J':
                num = all_num_genome_j_positions_in_loop[organism][ab][id]
                jloopseq = jseq[:num+5] ## go all the way up to and including the GXG
                jloopseqs[id] = jloopseq
        all_jloopseq_nbrs = {}
        for id1,seq1 in jloopseqs.iteritems():
            all_jloopseq_nbrs[id1] = []
            for id2,seq2 in jloopseqs.iteritems():
                #assert len(seq1) == len(seq2)
                if seq1 == seq2:
                    all_jloopseq_nbrs[id1].append( id2 )
        for id in all_jloopseq_nbrs:
            all_jseq_representative[organism][id] = min( all_jloopseq_nbrs[id] )
            assert jloopseqs[id] == jloopseqs[ all_jseq_representative[organism][id] ]
            if verbose:
                print 'jrep %s %15s %15s %15s'%(organism, id, all_jseq_representative[organism][id],
                                                jloopseqs[id])

    ## setup core positions
    all_core_positions[ organism ] = {}
    for id,alseq in align_fasta.iteritems():
        if id[2] in 'AB' and id[3] == 'V':
            ab = id[2]
            poslist = []
            poslist_seq = ''
            align_pos_mapping = {}
            for i in range(1,len(alseq)+1):
                if i in extra_alignment_columns[organism][ab]:continue
                last_mapped_pos = 0
                if align_pos_mapping: last_mapped_pos = max( align_pos_mapping.keys() )
                align_pos_mapping[ last_mapped_pos + 1 ] = i
                if False and verbose:
                    print 'align_pos_mapping:',organism,ab,last_mapped_pos+1,i

            for generic_align_pos1 in core_positions_generic_1indexed:
                align_pos = align_pos_mapping[ generic_align_pos1 ] - 1 ## now 0-indexed, shifted for extra columns
                numgaps = alseq[:align_pos].count(gap_character)
                if alseq[align_pos] in amino_acids:
                    poslist.append( align_pos - numgaps )
                else:
                    poslist.append( -1 )
                poslist_seq += alseq[ align_pos ]

            if verbose:
                print 'core_poslist_seq:',organism,poslist_seq,id

            all_core_positions[organism][id] = poslist

    all_align_fasta[ organism ] = align_fasta
    all_fasta[ organism ] = fasta

    #exit()


## return a list of tuples, each tuple is start and stop positions (0-indexed) of the loop
##
def parse_core_positions( organism, ab, qseq, v_hit ):
    fasta = all_fasta[organism]
    align_fasta = all_align_fasta[organism]

    v_gene = v_hit.hit_id
    v_seq = fasta[ v_gene ]
    v_alseq = align_fasta[ v_gene ]
    assert v_seq == ''.join( v_alseq.split(gap_character))


    core_positions = all_core_positions[ organism ][ v_gene ]

    qcore_positions = [-1]*len(core_positions)
    mismatches = 0

    for qpos, ( hpos, haa ) in v_hit.q2hmap.iteritems():
        if haa in amino_acids:
            assert haa == v_seq[hpos]
            if hpos in core_positions:
                qcore_positions[ core_positions.index( hpos ) ] = qpos
                if qseq[ qpos ] != haa:
                    mismatches += 1

    return qcore_positions, mismatches


## return a list of tuples, each tuple is start and stop positions (0-indexed) of the loop
##
def parse_other_cdrs( organism, ab, qseq, v_hit ):
    fasta = all_fasta[organism]
    align_fasta = all_align_fasta[organism]

    v_gene = v_hit.hit_id
    v_seq = fasta[ v_gene ]
    v_alseq = align_fasta[ v_gene ]
    assert v_seq == ''.join( v_alseq.split(gap_character))

    cdrs = pb_cdrs[organism][ v_gene ]

    assert len(cdrs) == 3 ## CDR1, CDR2, and CDRX

    qseq_cdrs = []
    mismatches = 0
    for loopseq,start,stop in cdrs:
        assert loopseq == v_seq[start:stop+1]
        ## what aligns to this region in v_hit
        q_loop_positions = []
        for qpos, ( hpos, haa ) in v_hit.q2hmap.iteritems():
            if haa in amino_acids:
                assert haa == v_seq[hpos]
                if hpos >= start and hpos <= stop:
                    q_loop_positions.append( qpos )
        if q_loop_positions:
            qstart = min( q_loop_positions )
            qstop  = max( q_loop_positions )
            qseq_cdrs.append( ( qseq[qstart:qstop+1], qstart, qstop ) )
            for i in range(qstart,qstop+1):
                if i not in v_hit.q2hmap or v_hit.q2hmap[i][1] != qseq[i]: mismatches += 1

        else:
            qseq_cdrs.append( ('-',-1,-1))
    return qseq_cdrs, mismatches




def get_cdr3_and_j_match_counts( organism, ab, qseq, j_gene, min_min_j_matchlen = 3,
                                 extended_cdr3 = False ):
    fasta = all_fasta[organism]

    errors = []

    ## qseq starts at CA...
    assert qseq[0] == 'C'

    num_genome_j_positions_in_loop = all_num_genome_j_positions_in_loop[organism][ab][j_gene]

    if extended_cdr3: num_genome_j_positions_in_loop += 2 ## up to but not including GXG

    ## history: was only for alpha
    aseq = qseq[:] ## starts at the C position

    ja_gene = j_gene
    assert ja_gene in fasta
    ja_seq = fasta[ ja_gene ]

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
    j_seq = fasta[ j_gene ]
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
    fasta = all_fasta[ organism ]
    align_fasta = all_align_fasta[ organism ]

    errors = []

    ## what is the C position in this v gene?
    v_seq = fasta[ v_gene ]
    v_alseq = align_fasta[ v_gene ]
    assert v_seq == ''.join( v_alseq.split(gap_character))

    alseq_cpos = alseq_C_pos[organism][ab] - 1 ## now 0-indexed
    numgaps = v_alseq[:alseq_cpos].count(gap_character)

    cpos = alseq_C_pos[organism][ab] - numgaps - 1 ## 0-indexed
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
    for org in all_num_genome_j_positions_in_loop:
        for ab in all_num_genome_j_positions_in_loop[org]:
            for id in all_num_genome_j_positions_in_loop[org][ab]:
                num = all_num_genome_j_positions_in_loop[org][ab][id]
                jseq = all_fasta[org][id]
                print 'jseq %10s %15s %10s %s'%(org,id,jseq[:num],jseq[num:])

