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


fasta_dir = path_to_db+'/old/fasta/'

verbose = False
#verbose = ( __name__ == '__main__' )

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


## these aren't used right now but could be added later
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



if __name__ == '__main__':

    ## silly but the nucseqs for D regions do not seem to be loaded by the other code in this file!
    dlines = """>K02545|TRBD1*01|Homo sapiens|F|D-REGION|82..93|12 nt|1| | | | |12+0=12| | |
gggacagggggc
>X02987|TRBD2*01|Homo sapiens|F|D-REGION|140..155|16 nt|1| | | | |16+0=16| | |
gggactagcggggggg
>M14159|TRBD2*02|Homo sapiens|F|D-REGION|569..584|16 nt|1| | | | |16+0=16| | |
gggactagcgggaggg
>X00933|TRBD1*01|Mus musculus|F|D-REGION|156..167|12 nt|1| | | | |12+0=12| | |
gggacagggggc
>X00934|TRBD2*01|Mus musculus|F|D-REGION|170..183|14 nt|1| | | | |14+0=14| | |
gggactgggggggc""".split('\n')
    d_fasta = {}
    for line in dlines:
        if line[0] == '>':
            l = line[1:].split('|')
            id = l[1]
            organism = {'Homo sapiens':'human', 'Mus musculus':'mouse'}[l[2]]
            if organism not in d_fasta:
                d_fasta[organism] = {}
        else:
            d_fasta[organism][id] = line



    ## this next part is taken from the old version of read_sanger_data
    #from genetic_code import genetic_code
    #import logo_tools
    from translation import get_translation

    all_offsets = {}
    all_fasta_sd = {}

    for organism in ['mouse','human']:
        all_offsets[organism] = {}
        all_fasta_sd[organism] = {}
        prot = 'protein'
        nuc = 'nucleotide'
        for ab in 'AB':
            all_offsets[organism][ab] ={}
            all_fasta_sd[organism][ab] ={}
            for vj in 'VJ':
                all_offsets[organism][ab][vj] ={}
                all_fasta_sd[organism][ab][vj] ={}
                myfasta = {}
                all_fastafile = {}
                for np in [prot,nuc]:
                    myfasta[np] ={}
                    fastafile = path_to_db+'/old/fasta/imgt_%s_TR_%s_sequences.fasta.TR%s%s.fasta'\
                                %( organism, np, ab, vj )
                    all_fastafile[np] = fastafile
                    assert exists( fastafile )

                    id = ''
                    for line in open( fastafile,'r'):
                        if line[0] == '>':
                            id = line[1:-1]
                            myfasta[np][id] = ''
                        else:
                            assert id
                            myfasta[np][id] += line[:-1]
                all_fasta_sd[organism][ab][vj] = myfasta
                for id in myfasta[prot]:
                    assert id in myfasta[nuc]
                    pseq = myfasta[prot][id]
                    nseq = myfasta[nuc][id]
                    myframe = -1
                    for i in range(3):
                        tseq = get_translation( nseq, '+{}'.format(i+1) )[0]
                        if pseq in tseq:
                            myframe = i + 3*tseq.index(pseq)
                    assert myframe >= 0
                    num_after = len(nseq) - 3*len(pseq) - myframe
                    all_offsets[organism][ab][vj][id] = ( myframe, num_after )




    ## make a single tsv file with the following fields
    ##
    ## id organism region
    ## chain is A or B -- where A means alpha-like (VJ recombining) and B means beta-like (VDJ recombining)
    ## region is V D J
    ##
    ## cdrs: comma-separated list of protein sequences for cdr regions
    ##
    outfields = "id organism chain region nucseq frame aligned_protseq cdr_columns cdrs".split()

    cdrs_sep = ';'

    outfile = 'db/alphabeta_db.tsv'
    out = open(outfile,'w')
    out.write('\t'.join( outfields )+'\n' )

    for organism in all_fasta: ## protein sequences
        for chain in 'AB':
            ## what ids are relevant here?
            v_ids = sorted( [x for x in all_merged_loopseqs[organism] if x[2] == chain ] )
            j_ids = sorted( [x for x in all_jseq_representative[organism] if x[2] == chain ] )

            # some sanity checks
            for id in v_ids + j_ids:
                assert id in all_fasta[organism]

            d_ids = []
            for id in all_fasta[organism]:
                if id[2] != chain: continue
                region = id[3]
                assert region in 'VDJC'
                if region in 'VJ':
                    assert id in v_ids or id in j_ids
                elif region == 'D':
                    d_ids.append( id )
            d_ids.sort()

            max_alseq_length = max( ( len(all_align_fasta[organism][x]) for x in v_ids ) )

            for id in v_ids:
                region = 'V'
                protseq = all_fasta[organism][id]
                alseq = all_align_fasta[organism][id]
                alseq += gap_character*( max_alseq_length - len(alseq))
                nucseq = all_fasta_sd[organism][chain][region][nuc][id]
                protseq2 = all_fasta_sd[organism][chain][region][prot][id]
                assert protseq == protseq2
                assert protseq == ''.join( alseq.split(gap_character))


                offset,num_after = all_offsets[organism][chain][region][id]
                frame = offset+1 ## 1-indexed frame
                protseq3 = get_translation( nucseq, '+{}'.format(frame))[0]
                #print organism,id,protseq,protseq3,nucseq
                assert protseq == protseq3

                ## get the cdrs
                loopseq = all_merged_loopseqs[ organism ][ id ].split()
                assert len(loopseq) == 3

                alseq_cpos = alseq_C_pos[organism][chain]-1 ## 0-indexed
                numgaps = alseq[:alseq_cpos].count(gap_character)
                cpos = alseq_cpos - numgaps
                cdr3_Nterm_protseq = protseq[cpos:] if cpos < len(protseq) else ''

                if not cdr3_Nterm_protseq or cdr3_Nterm_protseq[0] != 'C':
                    print 'funny CDR3:',organism, id, alseq


                Ntrunc = ( loopseq[0][0] == gap_character )
                if Ntrunc: ## N term truncation?
                    print 'Ntrunc:',organism,id,alseq,loopseq

                Ctrunc = not cdr3_Nterm_protseq

                if organism == 'mouse' and id == 'TRAV13D-1*01':
                    #-----------------------------------
                    #../../../tmp.blat:mismatch: V 6 imgt: a genome: t TRAV13D-1*01
                    #tmp.whoah:whoah  6 act: t  98.7 exp: a   1.1 TRAV13D-1*01 TRAV13-1*01 620
                    #tmp.whoah:whoah: expected: caaggtatcgtgt consensus: caaggtttcgtgt TRAV13D-1*01 620
                    #tmp.3.whoah:whoah  6 act: t  97.4 exp: a   1.4 TRAV13D-1*01 TRAV13-1*01 642
                    #tmp.3.whoah:whoah: expected: caaggtatcgtgt consensus: caaggtttcgtgt TRAV13D-1*01 642
                    #tmp.la_mc.whoah:whoah  6 act: t  89.0 exp: a   7.0 TRAV13D-1*01 TRAV13-1*01 100
                    #tmp.la_mc.whoah:whoah: expected: caaggtatcgtgt consensus: caaggtttcgtgt TRAV13D-1*01 100
                    old_cdr3_nucseq, new_cdr3_nucseq = 'tgtgctatggaac', 'tgtgctttggaac' ### CAM.. vs CAL..
                    assert nucseq.endswith( old_cdr3_nucseq )
                    nucseq = nucseq[: -1*len(old_cdr3_nucseq)] + new_cdr3_nucseq
                    cdr3_Nterm_protseq = get_translation( new_cdr3_nucseq, '+1' )[0]
                    old_protseq = get_translation( old_cdr3_nucseq, '+1' )[0]
                    print 'update V region:',id,new_cdr3_nucseq,cdr3_Nterm_protseq
                    assert old_protseq in alseq
                    alseq = alseq.replace( old_protseq, cdr3_Nterm_protseq )

                cdr_columns = pb_cdr_positions[organism][chain] + [[alseq_C_pos[organism][chain], max_alseq_length ]]

                cdrs_new = [ alseq[x[0]-1:x[1]] for x in cdr_columns ]

                assert loopseq == cdrs_new[:3]
                if cdr3_Nterm_protseq:
                    assert cdr3_Nterm_protseq == cdrs_new[3][:len(cdr3_Nterm_protseq)]

                outl = { 'id': id,
                         'organism': organism,
                         'functional': 0 if Ntrunc or Ctrunc else 1,
                         'chain': chain,
                         'region': region,
                         'nucseq': nucseq,
                         'aligned_protseq': alseq,
                         'frame': '+{}'.format( offset+1 ), ## convention for frame is 1-indexed
                         'cdr_columns':cdrs_sep.join( '{}-{}'.format(x[0],x[1]) for x in cdr_columns ),
                         'cdrs': cdrs_sep.join( cdrs_new ),
                         }

                out.write( make_tsv_line( outl, outfields )+'\n' )


            splits = []

            for id in j_ids:
                protseq = all_fasta[organism][id]
                # add 2 since we go all the way up to (but not including) the GXG
                num_genome_positions_in_CDR3 = all_num_genome_j_positions_in_loop[organism][chain][id] + 2
                assert num_genome_positions_in_CDR3 < len(protseq)
                splits.append( [ num_genome_positions_in_CDR3, len(protseq) - num_genome_positions_in_CDR3 ] )

            max_in  = max( ( x[0] for x in splits) )
            max_out = max( ( x[1] for x in splits) )

            for id in j_ids:
                region = 'J'
                protseq = all_fasta[organism][id]
                nucseq = all_fasta_sd[organism][chain][region][nuc][id]
                assert protseq == all_fasta_sd[organism][chain][region][prot][id]

                offset,num_after = all_offsets[organism][chain][region][id]

                num_genome_positions_in_CDR3 = all_num_genome_j_positions_in_loop[organism][chain][id] + 2


                if organism == 'mouse' and id == 'TRAJ47*01':
                    # -----------------------------------
                    # ../../../tmp.blat:mismatch: J 2 imgt: c genome: g TRAJ47*01
                    # ../../../tmp.blat:mismatch: J 24 imgt: g genome: t TRAJ47*01
                    # tmp.whoah:whoah  2 act: g  81.9 exp: c   4.7 TRAJ47*01 TRAJ47*01 1412
                    # tmp.whoah:whoah 24 act: t  82.7 exp: g  16.8 TRAJ47*01 TRAJ47*01 1412
                    # tmp.whoah:whoah: expected: tgcactatgcaaacaagatgatctgt consensus: tggactatgcaaacaagatgatcttt TRAJ47*01 1412
                    # tmp.3.whoah:whoah  2 act: g  81.6 exp: c   5.0 TRAJ47*01 TRAJ47*01 1362
                    # tmp.3.whoah:whoah 24 act: t  82.7 exp: g  16.6 TRAJ47*01 TRAJ47*01 1362
                    # tmp.3.whoah:whoah: expected: tgcactatgcaaacaagatgatctgt consensus: tggactatgcaaacaagatgatcttt TRAJ47*01 1362
                    # tmp.la_mc.whoah:whoah  2 act: g  79.6 exp: c   5.3 TRAJ47*01 TRAJ47*01 113
                    # tmp.la_mc.whoah:whoah 24 act: t  99.1 exp: g   0.9 TRAJ47*01 TRAJ47*01 113
                    # tmp.la_mc.whoah:whoah: expected: tgcactatgcaaacaagatgatctgt consensus: tggactatgcaaacaagatgatcttt TRAJ47*01 113
                    old_cdr3_nucseq = 'tgcactatgcaaacaagatgatctgt' ## C at end
                    new_cdr3_nucseq = 'tggactatgcaaacaagatgatcttt' ## F at end, SAME LENGTH!
                    num_cdr3_nucs = offset + 3*num_genome_positions_in_CDR3
                    assert nucseq[:num_cdr3_nucs] == old_cdr3_nucseq
                    old_protseq = get_translation( nucseq, '+{}'.format( offset+1 ) )[0]
                    assert old_protseq == protseq
                    nucseq = new_cdr3_nucseq + nucseq[ num_cdr3_nucs: ]
                    protseq = get_translation( nucseq, '+{}'.format( offset+1 ) )[0]

                elif organism == 'mouse' and id == 'TRAJ24*01':
                    # -----------------------------------
                    # ../../../tmp.blat:unaligned: J 0 TRAJ24*01
                    # ../../../tmp.blat:unaligned: J 1 TRAJ24*01
                    # ../../../tmp.blat:gapped: J 6 TRAJ24*01
                    # tmp.whoah:whoah  2 act: c  60.3 exp: a  15.3 TRAJ24*01 TRAJ24*01 464
                    # tmp.whoah:whoah  4 act: a  88.6 exp: c   2.8 TRAJ24*01 TRAJ24*01 464
                    # tmp.whoah:whoah  5 act: c  93.3 exp: t   1.5 TRAJ24*01 TRAJ24*01 464
                    # tmp.whoah:whoah  6 act: t  97.2 exp: g   1.1 TRAJ24*01 TRAJ24*01 464
                    # tmp.whoah:whoah: expected: tgaactggccagtttggggaaactgcagttt consensus: gacaactgccagtttggggaaactgcagttt TRAJ24*01 464
                    # tmp.3.whoah:whoah  2 act: c  60.8 exp: a  13.9 TRAJ24*01 TRAJ24*01 475
                    # tmp.3.whoah:whoah  4 act: a  86.3 exp: c   4.2 TRAJ24*01 TRAJ24*01 475
                    # tmp.3.whoah:whoah  5 act: c  94.5 exp: t   1.1 TRAJ24*01 TRAJ24*01 475
                    # tmp.3.whoah:whoah  6 act: t  98.1 exp: g   1.1 TRAJ24*01 TRAJ24*01 475
                    # tmp.3.whoah:whoah: expected: tgaactggccagtttggggaaactgcagttt consensus: gacaactgccagtttggggaaactgcagttt TRAJ24*01 475
                    # tmp.la_mc.whoah:whoah  2 act: c  75.3 exp: a   4.3 TRAJ24*01 TRAJ24*01 93
                    # tmp.la_mc.whoah:whoah  4 act: a  89.2 exp: c   2.2 TRAJ24*01 TRAJ24*01 93
                    # tmp.la_mc.whoah:whoah  5 act: c  97.8 exp: t   1.1 TRAJ24*01 TRAJ24*01 93
                    # tmp.la_mc.whoah:whoah  6 act: t  98.9 exp: g   0.0 TRAJ24*01 TRAJ24*01 93
                    # tmp.la_mc.whoah:whoah: expected: tgaactggccagtttggggaaactgcagttt consensus: gacaactgccagtttggggaaactgcagttt TRAJ24*01 93
                    old_cdr3_nucseq = 'tgaactggccagtttggggaaactgcagttt'
                    new_cdr3_nucseq = 'gacaactgccagtttggggaaactgcagttt' ## SAME LENGTH!
                    ## take the consensus
                    ## given that there's an indel (and the alignment to the genome starts at j sequence position 3)
                    ## it's hard to tell what to do at the beginning...
                    num_cdr3_nucs = offset + 3*num_genome_positions_in_CDR3
                    assert nucseq[:num_cdr3_nucs] == old_cdr3_nucseq
                    old_protseq = get_translation( nucseq, '+{}'.format( offset+1 ) )[0]
                    assert old_protseq == protseq
                    nucseq = new_cdr3_nucseq + nucseq[ num_cdr3_nucs: ]
                    protseq = get_translation( nucseq, '+{}'.format( offset+1 ) )[0]

                prefix = gap_character * ( max_in - num_genome_positions_in_CDR3 )
                suffix = gap_character * ( max_out - ( len(protseq) - num_genome_positions_in_CDR3 ) )
                alseq = prefix + protseq + suffix
                cdr_columns = [ [1,max_in] ]
                cdrs = [ alseq[:max_in] ]

                outl = { 'id': id,
                         'organism': organism,
                         'functional': 0 if '*' in protseq else 1,
                         'chain': chain,
                         'region': region,
                         'nucseq': nucseq,
                         'frame': '+{}'.format( offset+1 ), ## convention for frame is 1-indexed
                         'aligned_protseq':alseq,
                         'cdr_columns':'{}-{}'.format( cdr_columns[0][0], cdr_columns[0][1] ),
                         'cdrs': cdrs[0]
                         }

                out.write( make_tsv_line( outl, outfields )+'\n' )

            d_protseqs = [ get_translation( x, '+1' )[0] for x in d_fasta[organism].values() ]
            max_d_protseq_len = max( ( len(x) for x in d_protseqs ) )

            # D genes
            for id in d_ids:
                if id[2] != chain: continue
                region = 'D'
                frame = '+1'
                nucseq = d_fasta[organism][id]
                protseq = get_translation(nucseq,frame)[0]
                alseq = protseq+gap_character*( max_d_protseq_len-len(protseq))

                outl = { 'id': id,
                         'organism': organism,
                         'functional': 1,
                         'chain': chain,
                         'region': region,
                         'nucseq': nucseq,
                         'frame': frame,
                         'aligned_protseq': alseq,
                         'cdr_columns':'',
                         'cdrs': ''
                         }

                out.write( make_tsv_line( outl, outfields )+'\n' )

    out.close()
