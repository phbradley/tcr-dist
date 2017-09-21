import re
from basic import *
from amino_acids import amino_acids
from tcr_distances_blosum import blosum
from paths import path_to_db
from translation import get_translation

gap_character = '.'

verbose = False
#verbose = ( __name__ == '__main__' )

## look at imgt cdrs
## http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html
## imgt_cdr_positions = [ ( 27, 38 ), ( 56, 65 ) ] ## 1-indexed
## these are 1-indexed !!!
##
## note that the TRAV mouse alignment seems to be shifted by 1 relative to IMGT for FR1 and FR2 and by 2 for C104 (->106)
## looks like the insertion happens around TRAV-alpha alignment position 86
## the other three agree at anchor positions C23, W41, C104
##
## indexed by id[2]
pb_cdr_positions   = { 'mouse': {'A': [ ( 28, 39 ), ( 57, 66 ), (82, 88) ],
                                 'D': [ ( 28, 39 ), ( 57, 66 ), (82, 88) ], #since we shift
                                 'B': [ ( 27, 38 ), ( 56, 65 ), (81, 86) ],
                                 'G': [ ( 27, 38 ), ( 56, 65 ), (81, 86) ] },
                       'human': {'A': [ ( 27, 38 ), ( 56, 65 ), (81, 86) ],
                                 'B': [ ( 27, 38 ), ( 56, 65 ), (81, 86) ],
                                 'G': [ ( 27, 38 ), ( 56, 65 ), (81, 86) ],
                                 'D': [ ( 27, 38 ), ( 56, 65 ), (81, 86) ] },
}


CWC_positions = { 'mouse': { 'A': [24,42,106], 'D': [24,42,106], 'B': [23,41,104], 'G': [23,41,104] }, # use id[2]
                  'human': { 'A': [23,41,104], 'D': [23,41,104], 'B': [23,41,104], 'G': [23,41,104] } } #D-shift


## 1-indexed:
extra_alignment_columns = { 'mouse':{'A':[9,86],'B':[],'G':[],'D':[9,86] }, ## 1-indexed
                            'human':{'A':[],'B':[],'G':[],'D':[] } }

# core_positions_generic_1indexed = [
#     21, 23, 25,   ##  23 is C
#     39, 41,       ##  41 is W
#     53, 54, 55,
#     78,           ##            maybe also 80?
#     89,           ##  89 is L
#     102, 103, 104 ## 104 is C
# ]


outfields = "id organism chain region nucseq frame aligned_protseq cdr_columns cdrs".split()

cdrs_sep = ';'

outfile = 'db/gammadelta_db.tsv'
out = open(outfile,'w')
out.write('\t'.join( outfields )+'\n' )


for organism in [ 'mouse','human' ]:


    ## the alignments:
    big_fasta = {}
    fastadir = '/home/pbradley/tcr_scripts/db/genedb_090617/'
    all_functional = {}
    for fasta_tag, big_fasta_file in \
        [ ['align',fastadir+'IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP'],
          ['prot',fastadir+'IMGTGENEDB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP'],
          ['nuc',fastadir+'IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP']]:
        assert exists( big_fasta_file )

        desired_species = {'human':'Homo sapiens','mouse':'Mus musculus'}[ organism ]
        desired_regions = ['V-REGION','J-REGION','D-REGION']
        desired_prefix = 'TR'
        fasta = {}
        for line in open( big_fasta_file,'rU'):
            if line[0] == '>':
                id = ''
                l = line[1:-1].split('|')
                species = l[2]
                functional = l[3]
                region = l[4]
                prefix = l[1][:2]
                if desired_species in species and region in desired_regions and prefix == desired_prefix:
                    id = l[1]
                    assert region[0] == id[3]
                    fasta[id] = ''
                    if id in all_functional:
                        assert functional == all_functional[id]
                    else:
                        all_functional[id] = functional
                else:
                    id=''
                    if False and desired_species in species and region not in desired_regions:
                        print 'not region:',region,fasta_tag
                    if False and desired_species in species and prefix != 'TR':
                        print 'not prefix:',prefix,fasta_tag
            else:
                if id:
                    fasta[id] += line.split()[0]

        big_fasta[fasta_tag] = fasta
        print 'num_ids:',fasta_tag,len(fasta.keys())

    align_fasta = big_fasta['align']
    prot_fasta = big_fasta['prot']
    nuc_fasta = big_fasta['nuc']

    for chain in 'AB':
        ## get relevant V regions
        v_ids = []
        for id in align_fasta:
            if id[3] == 'V':
                if chain == 'A' and id[2] == 'G' or chain == 'B' and id[2] in 'AD':
                    v_ids.append( id )
        maxlen = max( ( len(align_fasta[x]) for x in v_ids ) )
        for id in sorted(v_ids):
            abdg = id[2]
            alseq = align_fasta[id]
            alseq += gap_character*( maxlen - len(alseq))

            ## for mouse delta, adjust alignment to match mouse alpha
            if organism=='mouse' and abdg=='D':
                pos1=extra_alignment_columns[organism]['A'][0] - 1
                pos2=extra_alignment_columns[organism]['A'][1] - 1
                alseq = alseq[:pos1]+gap_character+alseq[pos1:]
                alseq = alseq[:pos2]+gap_character+alseq[pos2:]
                assert alseq.endswith(gap_character+gap_character)
                alseq = alseq[:-2]

            assert len(alseq) == maxlen

            cwc = ''.join( ( alseq[x-1] for x in CWC_positions[organism][abdg] ) )
            if cwc != 'CWC' and len(cwc.replace(gap_character,''))==3:
                print 'funny CWC:',cwc,alseq, organism, id, all_functional[id]

            extraseq = ''.join( ( alseq[x-1] for x in extra_alignment_columns[organism][abdg] ) )
            if extraseq and extraseq.replace(gap_character,'') :
                print 'extra:',extraseq, organism, id, all_functional[id]


            protseq = prot_fasta[id]
            nucseq = nuc_fasta[id]

            assert protseq == alseq.replace(gap_character,'')

            #print ' ',protseq
            myframe = -1
            for i in range(3):
                tseq = get_translation( nucseq, '+{}'.format(i+1) )[0]
                #print i,tseq
                if protseq in tseq:
                    myframe = i + 3*tseq.index(protseq)
            #print id, myframe
            if myframe==-1:
                print 'bad frame:',id, myframe,protseq ### NOTE SKIPPING THIS ONE WITH A BAD PROTEIN SEQUENCE
                continue
            assert myframe >= 0 and myframe<3

            cpos = CWC_positions[organism][abdg][-1] # 1-indexed
            cdr_columns = pb_cdr_positions[organism][abdg] + [[cpos,maxlen]] ## all 1-indexed
            cdrs = [ alseq[x[0]-1:x[1]] for x in cdr_columns ]

            region = 'V'

            outl = { 'id': id,
                     'organism': organism,
                     #'functional': 0 if Ntrunc or Ctrunc else 1,
                     'chain': chain,
                     'region': region,
                     'nucseq': nucseq,
                     'aligned_protseq': alseq,
                     'frame': '+{}'.format( myframe+1 ), ## convention for frame is 1-indexed
                     'cdr_columns':cdrs_sep.join( '{}-{}'.format(x[0],x[1]) for x in cdr_columns ),
                     'cdrs': cdrs_sep.join( cdrs ),
                     }

            out.write( make_tsv_line( outl, outfields )+'\n' )



        ## now the J regions
        j_ids = []
        for id in align_fasta:
            if id[3] == 'J':
                if chain == 'A' and id[2] == 'G' or chain == 'B' and id[2] == 'D':
                    j_ids.append( id )


        bounds = {}
        for id in j_ids:
            jseq = prot_fasta[id]
            #print 'jseq:',organism, chain, jseq, id
            m = re.search( 'F[AG].G', jseq )
            assert m
            txt = m.group(0)
            assert jseq.count(txt)==1
            num_in = jseq.index(txt)+1 # in the CDR3
            num_out = len(jseq) - num_in
            bounds[id] = [ num_in, num_out ]

        maxin = max( ( x[0] for x in bounds.values()))
        maxout = max( ( x[1] for x in bounds.values()))


        #maxlen = max( ( len(align_fasta[x]) for x in v_ids ) )
        for id in j_ids:
            jseq = prot_fasta[id]
            num_in,num_out = bounds[id]
            alseq = gap_character*(maxin -num_in) + jseq + gap_character*(maxout-num_out)
            print 'jseq:',organism, chain, alseq, id

            protseq = prot_fasta[id]
            nucseq = nuc_fasta[id]

            assert protseq == alseq.replace(gap_character,'')

            #print ' ',protseq
            myframe = -1
            for i in range(3):
                tseq = get_translation( nucseq, '+{}'.format(i+1) )[0]
                #print i,tseq
                if protseq in tseq:
                    myframe = i + 3*tseq.index(protseq)
            #print id, myframe
            if myframe==-1:
                print 'bad frame:',id, myframe,protseq ### NOTE SKIPPING THIS ONE WITH A BAD PROTEIN SEQUENCE
                continue
            assert myframe >= 0 and myframe<3

            cdr_columns = [[1,maxin]]
            cdrs = [ alseq[:maxin]]

            region = 'J'

            outl = { 'id': id,
                     'organism': organism,
                     #'functional': 0 if Ntrunc or Ctrunc else 1,
                     'chain': chain,
                     'region': region,
                     'nucseq': nucseq,
                     'aligned_protseq': alseq,
                     'frame': '+{}'.format( myframe+1 ), ## convention for frame is 1-indexed
                     'cdr_columns':cdrs_sep.join( '{}-{}'.format(x[0],x[1]) for x in cdr_columns ),
                     'cdrs': cdrs_sep.join( cdrs ),
                     }

            out.write( make_tsv_line( outl, outfields )+'\n' )



        if chain == 'B':
            ## now the D regions
            d_ids = []
            for id in align_fasta:
                if id[2:4] == 'DD':
                    d_ids.append( id )


            maxlen = max( ( len(prot_fasta[x]) for x in d_ids ) )


            #maxlen = max( ( len(align_fasta[x]) for x in v_ids ) )
            for id in d_ids:
                protseq = prot_fasta[id]
                alseq = protseq + gap_character*(maxlen-len(protseq))
                nucseq = nuc_fasta[id]

                assert protseq == alseq.replace(gap_character,'')

                #print ' ',protseq
                myframe = -1
                for i in range(3):
                    tseq = get_translation( nucseq, '+{}'.format(i+1) )[0]
                    #print i,tseq
                    if protseq in tseq:
                        myframe = i + 3*tseq.index(protseq)
                #print id, myframe
                if myframe==-1:
                    print 'bad frame:',id, myframe,protseq ### NOTE SKIPPING THIS ONE WITH A BAD PROTEIN SEQUENCE
                    continue
                assert myframe >= 0 and myframe<3

                region = 'D'

                outl = { 'id': id,
                         'organism': organism,
                         #'functional': 0 if Ntrunc or Ctrunc else 1,
                         'chain': chain,
                         'region': region,
                         'nucseq': nucseq,
                         'aligned_protseq': alseq,
                         'frame': '+{}'.format( myframe+1 ), ## convention for frame is 1-indexed
                         'cdr_columns':'',
                         'cdrs': ''
                         }

                out.write( make_tsv_line( outl, outfields )+'\n' )
out.close()
