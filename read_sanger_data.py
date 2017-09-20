from basic import *
import logo_tools
import sys
import random
import blast
from genetic_code import genetic_code
from all_genes import all_genes, gap_character
import parse_cdr3
from paths import path_to_blast_executables
from translation import get_translation

AB = 'AB'

def get_blast_nucseq_database( organism, chain, region ):
    """The nucleotide sequence database"""
    assert chain in 'AB' and region in 'VJ'
    db_dir = path_to_current_db_files()
    assert exists(db_dir)
    blast_dir = db_dir+'/blast_dbs'
    if not exists(blast_dir): mkdir( blast_dir )
    fastafile = blast_dir+'/nucseq_{}_{}_{}.fasta'.format( organism, chain, region )
    if not exists( fastafile ):
        ids = sorted( ( id for (id,g) in all_genes[organism].iteritems() if g.chain== chain and g.region==region ) )
        out = open(fastafile,'w')
        for id in ids:
            out.write('>{}\n{}\n'.format( id, all_genes[organism][id].nucseq ) )
        out.close()
    ## now make the blast db files
    dbfiles = glob(fastafile+'.n*')
    if not dbfiles:
        cmd = '{}/formatdb -p F -i {}'.format( path_to_blast_executables, fastafile )
        Log(cmd)
        system(cmd)
        dbfiles = glob(fastafile+'.n*')
        if not dbfiles:
            Log('blast db creation failed!')
            exit()
    return fastafile

all_offsets = {}
all_fasta = {}

prot = 'protein'
nuc = 'nucleotide'


for organism in all_genes:
    all_offsets[organism] = {}
    all_fasta[organism] = {}
    for ab in AB:
        all_offsets[organism][ab] ={}
        all_fasta[organism][ab] ={}
        for vj in 'VJ':
            ids = sorted( ( id for id,g in all_genes[organism].iteritems() if g.chain == ab and g.region == vj ) )

            all_offsets[organism][ab][vj] ={}
            all_fasta[organism][ab][vj] ={prot:{},nuc:{}}
            for id in ids:
                g = all_genes[organism][id]
                all_fasta[organism][ab][vj][nuc][id] = g.nucseq
                all_fasta[organism][ab][vj][prot][id] = g.protseq
                all_offsets[organism][ab][vj][id] = g.nucseq_offset

def reverse_q2hmap( qseq, hseq, hit ):
    assert hit.h_strand == -1
    q2hmap = {}
    for qpos,(hpos,na) in hit.q2hmap.iteritems():
        if hpos>=0:
            new_qpos = len(qseq)-qpos-1
            new_na = logo_tools.base_partner[ na ]
            assert new_na == hseq[ hpos ]
            q2hmap[new_qpos] = (hpos,new_na)
    return q2hmap

## should be just a single query sequence
def get_all_hits_with_evalues_and_scores( blastfile ):
    query = ''
    in_hits = False
    hits = []
    for line in open( blastfile,'r'):
        if line.startswith('Query='):
            assert not query
            query= line.split()[1]
        elif line.startswith('Sequences producing'):
            in_hits = True
        elif line.startswith('>'):
            assert in_hits
            break
        elif in_hits:
            l = line.split()
            if len(l) >= 3:
                evalue = l[-1]
                if evalue.startswith('e'): evalue = '1'+evalue
                bitscore = int( l[-2] )
                hitid = l[0]
                if '|' in hitid: hitid = hitid.split('|')[1]
                hits.append( ( hitid, bitscore, float(evalue) ) )
                #print 'new_hit:',hits[-1]
    return hits


##
##    returns genes,evalues,status ## status is a list, may be empty
##
def parse_unpaired_dna_sequence_blastn( organism, ab, blast_seq, info,
                                        verbose, nocleanup, hide_nucseq,
                                        extended_cdr3,
                                        return_all_good_hits = False,
                                        max_bit_score_delta_for_good_hits = 50,
                                        max_missing_aas_at_cdr3_cterm = 2 ):

    ## make this a little more unique
    blast_tmpfile = 'tmp%d%s%s%f%s.fa'%(len(blast_seq),organism,ab,random.random(),blast_seq[:3])
    #print 'blast_tmpfile:',blast_tmpfile
    #assert not exists(blast_tmpfile)

    genes =  ( 'UNK', 'UNK', [100,0], 'UNK', 'UNK', [100,0], '-' )

    status = []

    evalues = {'V'+ab:(1,0),'J'+ab:(1,0)}

    all_good_hits_with_scores = [[],[]]

    if verbose:
        print 'blast_seq:',info,ab,blast_seq

    if len(blast_seq) <= 20:
        status.append('short_{}_blast_seq_{}'.format(ab,len(blast_seq)))
    else:

        out = open(blast_tmpfile,'w')
        out.write('>tmp\n%s\n'%blast_seq)
        out.close()

        ## now blast against V and J
        top_hits = []

        for ivj,vj in enumerate('VJ'):
            dbfile = get_blast_nucseq_database( organism, ab, vj ) # also ensures that it exists
            assert exists(dbfile)
            blastall_exe = path_to_blast_executables+'/blastall'
            assert exists( blastall_exe )
            cmd = '%s -F F -p blastn -i %s -d %s -v 100 -b 1 -o %s.blast'\
                  %( blastall_exe, blast_tmpfile, dbfile, blast_tmpfile )
            #print cmd
            system(cmd)

            if verbose:
                print 'blast:',info,ab,vj, '='*50
                print ''.join( open(blast_tmpfile+'.blast','r').readlines())
                print '='*80

            ## try parsing the results
            evalue_threshold = 1e-1
            identity_threshold = 20
            hits = blast.parse_blast_alignments( blast_tmpfile+'.blast', evalue_threshold, identity_threshold )
            hits_scores = get_all_hits_with_evalues_and_scores( blast_tmpfile+'.blast' ) ## id,bitscore,evalue
            if hits and hits[ hits.keys()[0]]:

                top_hit = hits[ hits.keys()[0] ][0]
                top_id, top_bit_score, top_evalue = hits_scores[0]

                all_good_hits_with_scores[ivj] \
                    = [ x for x in hits_scores if top_bit_score-x[1] <= max_bit_score_delta_for_good_hits ]

                assert top_hit.hit_id == top_id
                ## figure out the score gap to the next non-equivalen
                bit_score_gap = top_bit_score
                top_rep = all_genes[ organism ][top_id].rep
                for ( id, bit_score, evalue ) in hits_scores[1:]:
                    if all_genes[organism][id].rep != top_rep:
                        bit_score_gap = top_bit_score - bit_score
                        break
                evalues[vj+ab] = ( top_hit.evalue, bit_score_gap )
                top_hits.append( top_hit )
            else:
                status.append('no_{}{}_blast_hits'.format(vj,ab))

        if len(top_hits) == 2: ## hits in both v and j

            v_hit = top_hits[0]
            j_hit = top_hits[1]

            v_gene = v_hit.hit_id
            j_gene = j_hit.hit_id
            v_rep = all_genes[organism][v_gene].rep
            j_rep = all_genes[organism][j_gene].rep

            v_nucseq = all_fasta[organism][ab]['V'][nuc][v_hit.hit_id]
            j_nucseq = all_fasta[organism][ab]['J'][nuc][j_hit.hit_id]
            v_protseq = all_fasta[organism][ab]['V'][prot][v_hit.hit_id]

            ## this might fail if these guys are pseudo-genes...
            ## so filter out the non-aa-matching j genes...
            ##
            v_hitseq_frame = all_offsets[organism][ab]['V'][ v_hit.hit_id ]
            j_hitseq_frame = all_offsets[organism][ab]['J'][ j_hit.hit_id ]

            ## tricky if the hits are on different strands!
            ##
            assert v_hit.q_strand == 1 ## I think this is the blastn convention...
            assert j_hit.q_strand == 1


            if v_hit.h_strand != j_hit.h_strand:
                Log(`( 'ERR V/J strand mismatch:',v_hit.h_strand, v_hit.evalue, j_hit.h_strand, j_hit.evalue)`)
                genes = ( v_gene.replace('TRAV','TRaV' ).replace('TRBV','TRbV'),
                          v_rep .replace('TRAV','TRaV' ).replace('TRBV','TRbV'), [100,0],
                          j_gene.replace('TRAJ','TRaJ' ).replace('TRBJ','TRbJ'),
                          j_rep .replace('TRAJ','TRaJ' ).replace('TRBJ','TRbJ'), [100,0], '-' )
                status.append('vj_{}_strand_mismatch'.format(ab))
            else:

                v_q2hmap = v_hit.q2hmap
                j_q2hmap = j_hit.q2hmap

                if v_hit.h_strand == -1:
                    ## switch stuff around...
                    ## have to mess with the alignment too
                    v_q2hmap = reverse_q2hmap( blast_seq, v_nucseq, v_hit )
                    j_q2hmap = reverse_q2hmap( blast_seq, j_nucseq, j_hit )

                    blast_seq = logo_tools.reverse_complement( blast_seq )
                    if verbose:
                        print 'reverse-comp blast_seq:',ab


                q_vframes = {}
                for qpos,(vpos,vna) in v_q2hmap.iteritems():
                    if vpos>=0:
                        f = ( qpos - vpos + v_hitseq_frame )%3
                        q_vframes[f] = q_vframes.get(f,0)+1
                q_vframe = max( [(count,x) for x,count in q_vframes.iteritems() ] )[1]

                q_jframes = {}
                for qpos,(jpos,jna) in j_q2hmap.iteritems():
                    if jpos>=0:
                        f = ( qpos - jpos + j_hitseq_frame )%3
                        q_jframes[f] = q_jframes.get(f,0)+1

                q_jframe = max( [(count,x) for x,count in q_jframes.iteritems() ] )[1]



                #q_frame_vstart = ( v_hitseq_frame + v_hit.q_start - v_hit.h_start )%3
                #q_frame_jstart = ( j_hitseq_frame + j_hit.q_start - j_hit.h_start )%3

                ## construct a protein sequence alignment between translation of blast_seq and
                q2v_align = {}
                for qpos,(vpos,vna) in sorted( v_q2hmap.iteritems() ):
                    if vpos>=0:
                        f = ( qpos - vpos + v_hitseq_frame )%3
                        if f != q_vframe: continue
                        v_protpos = ( vpos - v_hitseq_frame )/3
                        q_protpos = ( qpos - q_vframe )/3
                        if q_protpos in q2v_align:
                            if q2v_align[ q_protpos ] != v_protpos:
                                Log('indel?? {} {} {}'.format(organism,ab,info))

                        q2v_align[ q_protpos ] = v_protpos
                        ## this could be aligning a position that's not actually in the translated protein
                        ## sequence if there are 1 or 2 nucleotides at the end...



                if q_vframe != q_jframe: ## out of frame
                    Log(`( 'ERR frame mismatch:',q_vframe, v_hit.evalue, q_jframe, j_hit.evalue)`)
                    if verbose:
                        print 'frame mismatch',q_vframe,q_jframe
                    # genes = ( v_gene.replace('TRAV','TRaV' ).replace('TRBV','TRbV'),
                    #           v_rep .replace('TRAV','TRaV' ).replace('TRBV','TRbV'), [100,0],
                    #           j_gene.replace('TRAJ','TRaJ' ).replace('TRBJ','TRbJ'),
                    #           j_rep .replace('TRAJ','TRaJ' ).replace('TRBJ','TRbJ'), [100,0], '-' )
                    status.append('vj_{}_frame_mismatch'.format(ab))

                    ## fiddle with blast_seq
                    ## for each 'extra' nucleotide inserted between v and j, add two '#' characters after the nucleotide
                    last_v_align_pos = max( v_q2hmap.keys() )
                    first_j_align_pos = min( j_q2hmap.keys() )

                    ## add some '#' characters to blast_seq to get V and J back into frame
                    num_to_insert = (q_vframe - q_jframe)%3
                    insertpos = max(last_v_align_pos+1, (last_v_align_pos+first_j_align_pos)/2)

                    blast_seq = blast_seq[:insertpos] + '#'*num_to_insert + blast_seq[insertpos:]

                    # num_inserted_nucleotides = (q_jframe - q_vframe)%3
                    # new_blast_seq = blast_seq[:last_q_align_pos+1]
                    # extra_seq = blast_seq[last_q_align_pos+1:]
                    # for i in range(num_inserted_nucleotides):
                    #     new_blast_seq += extra_seq[0] + '##'
                    #     extra_seq = extra_seq[1:]
                    # new_blast_seq += extra_seq
                    # blast_seq = new_blast_seq[:]


                qseq, codons = get_translation( blast_seq, '+%d'%(q_vframe+1) )
                cdr3,v_mm,j_mm,errors = parse_cdr3.parse_cdr3( organism, ab, qseq, v_hit.hit_id, j_hit.hit_id,
                                                               q2v_align, extended_cdr3 = extended_cdr3,
                                                               max_missing_aas_at_cdr3_cterm=max_missing_aas_at_cdr3_cterm )
                if verbose:
                    print 'cdr3:',ab,cdr3,cdr3 in qseq,'q_vframe:',q_vframe

                status.extend(errors)

                if cdr3 != '-':
                    ## the cdr3 sequence should be contained in qseq, unless qseq was missing 1-2 rsds at cterm
                    if not hide_nucseq:
                        if cdr3 in qseq: ## the old way, without any missing C-term rsds of CDR3
                            offset = qseq.find(cdr3)
                            cdr3_codons = codons[ offset:offset+len(cdr3) ]
                            cdr3 += '-'+''.join(cdr3_codons)
                        else:
                            num_missing_cterm_aas = 1
                            while num_missing_cterm_aas < max_missing_aas_at_cdr3_cterm and \
                                  cdr3[:-1*num_missing_cterm_aas] not in qseq:
                                num_missing_cterm_aas += 1
                            assert cdr3[:-1*num_missing_cterm_aas] in qseq
                            ## this is a nuisance...
                            assert extended_cdr3 # it's the new default anyhow
                            jg = all_genes[organism][j_hit.hit_id]
                            j_nucseq = jg.nucseq
                            j_cdr3len = len( jg.cdrs[0].replace(gap_character,''))
                            j_cdr3_nucseq = jg.nucseq[: jg.nucseq_offset + 3*j_cdr3len ]
                            missing_nucseq = j_cdr3_nucseq[-3*num_missing_cterm_aas:]
                            offset = qseq.find(cdr3[:-1*num_missing_cterm_aas])
                            cdr3_codons = codons[ offset:offset+len(cdr3)-num_missing_cterm_aas ]
                            cdr3_nucseq = ''.join(cdr3_codons)+missing_nucseq
                            assert len(cdr3_nucseq) == 3 * len(cdr3)
                            assert get_translation( cdr3_nucseq,'+1')[0] == cdr3
                            cdr3 += '-'+cdr3_nucseq

                        # if verbose:
                        #     cdr3_nucseq = ''.join( cdr3_codons ).upper()
                        #     nucseq_startpos = 3*offset + q_vframe
                        #     alt_nucseq = blast_seq[ nucseq_startpos:nucseq_startpos+len(cdr3_nucseq) ]
                        #     rc1 = logo_tools.reverse_complement(cdr3_nucseq)
                        #     rc2 = logo_tools.reverse_complement(blast_seq)
                        #     print 'cdr3_nucseq',ab,offset,cdr3_nucseq,cdr3_nucseq in blast_seq,\
                        #         blast_seq.index(cdr3_nucseq),alt_nucseq,rc1 in rc2
                if '#' in blast_seq: ## sign of out-of-frame
                    v_gene = v_gene.replace('TRAV','TRaV' ).replace('TRBV','TRbV')
                    v_rep  = v_rep .replace('TRAV','TRaV' ).replace('TRBV','TRbV')
                    j_gene = j_gene.replace('TRAJ','TRaJ' ).replace('TRBJ','TRbJ')
                    j_rep  = j_rep .replace('TRAJ','TRaJ' ).replace('TRBJ','TRbJ')
                    protseq,nucseq = cdr3.split('-')
                    if protseq and nucseq:
                        if protseq.count('#') != nucseq.count('#'):
                            assert nucseq.count('#') == 2
                            assert protseq.count('#') == 1
                            protseq = protseq.replace('#','##')
                            cdr3 = '{}-{}'.format(protseq,nucseq)

                genes = ( v_gene, v_rep, v_mm, j_gene, j_rep, j_mm, cdr3 )

                if cdr3 != "-":
                    cdr3aa = cdr3.split("-")[0]
                    if len(cdr3aa) < 5:
                        status.append('cdr3{}_len_too_short'.format(ab))

    if not nocleanup:
        files = glob(blast_tmpfile+'*')
        for file in files:
            remove( file )

    assert len(genes) == 7

    if return_all_good_hits:
        return genes,evalues,status,all_good_hits_with_scores ## status is a list, maybe be empty
    else:
        return genes,evalues,status ## status is a list, maybe be empty




def parse_paired_dna_sequences( program, organism, aseq, bseq, info = '',
                                verbose = False, nocleanup = False, hide_nucseq = False,
                                extended_cdr3 = False, return_all_good_hits = False ):
    if program == 'blastn':
        if return_all_good_hits:
            a_genes, a_evalues, a_status, a_hits \
                = parse_unpaired_dna_sequence_blastn( organism, 'A', aseq, info,
                                                      verbose, nocleanup, hide_nucseq,
                                                      extended_cdr3, return_all_good_hits=True )

            b_genes, b_evalues, b_status, b_hits \
                = parse_unpaired_dna_sequence_blastn( organism, 'B', bseq, info,
                                                      verbose, nocleanup, hide_nucseq,
                                                      extended_cdr3, return_all_good_hits=True )
        else:
            a_genes, a_evalues, a_status = parse_unpaired_dna_sequence_blastn( organism, 'A', aseq, info,
                                                                               verbose, nocleanup, hide_nucseq,
                                                                               extended_cdr3, return_all_good_hits )

            b_genes, b_evalues, b_status = parse_unpaired_dna_sequence_blastn( organism, 'B', bseq, info,
                                                                               verbose, nocleanup, hide_nucseq,
                                                                               extended_cdr3, return_all_good_hits )
    elif program == 'blastx':
        ## can recover this code from early versions of the repository if need be...
        print 'blastx parsing not supported at the moment'
        Log( 'blastx parsing not supported at the moment' )
        exit()
        # assert not return_all_good_hits
        # a_genes, a_evalues, a_status = parse_unpaired_dna_sequence_blastx( organism, 'A', aseq, info,
        #                                                                    verbose, nocleanup, hide_nucseq,
        #                                                                    extended_cdr3 )

        # b_genes, b_evalues, b_status = parse_unpaired_dna_sequence_blastx( organism, 'B', bseq, info,
        #                                                                    verbose, nocleanup, hide_nucseq,
        #                                                                    extended_cdr3 )
    else:
        Log('bad program: '+program)
        assert False

    assert len(a_genes) == 7 ## v_gene, v_rep, v_mm, j_gene, j_rep, j_mm, cdr3
    assert len(b_genes) == 7 ## ditto

    ab_genes = {'A': a_genes, 'B': b_genes }
    evalues = {}
    for tag,evalue in a_evalues.iteritems():
        evalues[tag] = evalue
    for tag,evalue in b_evalues.iteritems():
        evalues[tag] = evalue

    return ab_genes, evalues, a_status+b_status, [a_hits,b_hits]



if __name__ == '__main__':
    ## make all the dbs
    for organism in all_genes:
        for ab in 'AB':
            for vj in 'VJ':
                dbfile = get_blast_nucseq_database( organism, ab, vj ) # also ensures that it exists
                print 'dbfile:',dbfile
