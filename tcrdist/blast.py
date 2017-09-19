from glob import glob
import tempfile
import logging
import os
import os.path as op

from .paths import path_to_db, path_to_blast_executables
from .all_genes import all_genes
from . import cdr3s_human
from . import logo_tools
from .translation import get_translation

"""Both these imports are for calling two different versions of parse_cdr3??"""
from . import cdr3s_human
from . import parse_cdr3

logger = logging.getLogger('blast.py')

def getFASTAOffsets(all_genes):
    all_offsets = {}
    all_fasta = {}

    for organism in all_genes:
        all_offsets[organism] = {}
        all_fasta[organism] = {}
        for ab in 'AB':
            all_offsets[organism][ab] ={}
            all_fasta[organism][ab] ={}
            for vj in 'VJ':
                ids = sorted( ( id for id,g in all_genes[organism].iteritems() if g.chain == ab and g.region == vj ) )

                all_offsets[organism][ab][vj] ={}
                all_fasta[organism][ab][vj] ={'protein':{},'nucleotide':{}}
                for id in ids:
                    g = all_genes[organism][id]
                    all_fasta[organism][ab][vj]['nucleotide'][id] = g.nucseq
                    all_fasta[organism][ab][vj]['protein'][id] = g.protseq
                    all_offsets[organism][ab][vj][id] = g.nucseq_offset
    return all_offsets, all_fasta

all_offsets, all_fasta = getFASTAOffsets(all_genes)

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

def get_blast_nucseq_database( organism, chain, region ):
    """The nucleotide sequence database"""
    assert chain in 'AB' and region in 'VJ'
    db_dir = path_to_db
    assert op.exists(db_dir)
    blast_dir = db_dir+'/blast_dbs'
    if not op.exists(blast_dir): os.mkdir( blast_dir )
    fastafile = blast_dir+'/nucseq_{}_{}_{}.fasta'.format( organism, chain, region )
    if not op.exists( fastafile ):
        ids = sorted( ( id for (id,g) in all_genes[organism].iteritems() if g.chain== chain and g.region==region ) )
        with open(fastafile,'w') as out:
            for id in ids:
                out.write('>{}\n{}\n'.format( id, all_genes[organism][id].nucseq ) )

    ## now make the blast db files
    dbfiles = glob(fastafile+'.n*')
    if not dbfiles:
        cmd = '{}/formatdb -p F -i {}'.format( path_to_blast_executables, fastafile )
        logger.debug(cmd)
        os.system(cmd)
        dbfiles = glob(fastafile+'.n*')
        if not dbfiles:
            logger.error('blast db creation failed!')
            exit()
    return fastafile

def parse_unpaired_dna_sequence_blastn( organism, ab, blast_seq, info,
                                        nocleanup, hide_nucseq,
                                        extended_cdr3,
                                        return_all_good_hits = False,
                                        max_bit_score_delta_for_good_hits = 50 ):

    if not op.exists('./tmp'):
        os.mkdir('./tmp')
    handle, blast_tmpfile = tempfile.mkstemp(suffix='.fa', prefix='blasttmp', dir='./tmp')
    os.close(handle)

    genes =  ( 'UNK', 'UNK', [100,0], 'UNK', 'UNK', [100,0], '-' )

    status = []

    evalues = {'V'+ab:(1,0),'J'+ab:(1,0)}

    all_good_hits_with_scores = [[],[]]

    logger.debug('blast_seq: %s, %s, %s',info,ab,blast_seq)

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
            assert op.exists(dbfile)
            blastall_exe = path_to_blast_executables+'/blastall'
            assert op.exists( blastall_exe )
            cmd = '%s -F F -p blastn -i %s -d %s -v 100 -b 1 -o %s.blast'\
                  %( blastall_exe, blast_tmpfile, dbfile, blast_tmpfile )
            logger.debug('blast cmd: %s', cmd)
            os.system(cmd)

            logger.debug('blast: %s, %s, %s',info,ab,vj)
            logger.debug(''.join( open(blast_tmpfile+'.blast','r').readlines()))

            ## try parsing the results
            evalue_threshold = 1e-1
            identity_threshold = 20
            hits = parse_blast_alignments( blast_tmpfile+'.blast', evalue_threshold, identity_threshold )
            hits_scores = get_all_hits_with_evalues_and_scores( blast_tmpfile+'.blast' ) ## id,bitscore,evalue
            if hits and hits[ hits.keys()[0]]:

                top_hit = hits[ hits.keys()[0] ][0]
                top_id, top_bit_score, top_evalue = hits_scores[0]

                all_good_hits_with_scores[ivj] \
                    = [ x for x in hits_scores if top_bit_score-x[1] <= max_bit_score_delta_for_good_hits ]

                assert top_hit.hit_id == top_id
                ## figure out the score gap to the next non-equivalen
                bit_score_gap = top_bit_score
                if vj == 'V':
                    old_rep_map = cdr3s_human.all_loopseq_representative[organism]
                else:
                    old_rep_map = cdr3s_human.all_jseq_representative[organism]
                old_top_rep = old_rep_map[top_id]
                top_rep = all_genes[ organism ][top_id].rep
                assert old_top_rep == top_rep
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
            old_v_rep = cdr3s_human.all_loopseq_representative[ organism ][ v_gene ]
            old_j_rep = cdr3s_human.all_jseq_representative[ organism ][ j_gene ]
            v_rep = all_genes[organism][v_gene].rep
            j_rep = all_genes[organism][j_gene].rep
            assert old_v_rep == v_rep and old_j_rep == j_rep

            v_nucseq = all_fasta[organism][ab]['V']['nucleotide'][v_hit.hit_id]
            j_nucseq = all_fasta[organism][ab]['J']['nucleotide'][j_hit.hit_id]
            v_protseq = all_fasta[organism][ab]['V']['protein'][v_hit.hit_id]

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
                logger.error('ERR V/J strand mismatch: %s, %s, %s, %s',v_hit.h_strand, v_hit.evalue, j_hit.h_strand, j_hit.evalue)
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
                    logger.debug('reverse-comp blast_seq: %s', ab)


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
                                logger.error('indel?? {} {} {}'.format(organism,ab,info))

                        q2v_align[ q_protpos ] = v_protpos
                        ## this could be aligning a position that's not actually in the translated protein
                        ## sequence if there are 1 or 2 nucleotides at the end...



                if q_vframe != q_jframe: ## out of frame
                    logger.error('ERR frame mismatch: %s, %s, %s, %s',q_vframe, v_hit.evalue, q_jframe, j_hit.evalue)
                    logger.debug('frame mismatch: %s, %s',q_vframe,q_jframe)
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
                old_cdr3,v_mm,j_mm,errors = cdr3s_human.parse_cdr3( organism, ab, qseq, v_hit.hit_id, j_hit.hit_id,
                                                                    q2v_align, extended_cdr3 = extended_cdr3 )
                cdr3,v_mm,j_mm,errors = parse_cdr3.parse_cdr3( organism, ab, qseq, v_hit.hit_id, j_hit.hit_id,
                                                               q2v_align, extended_cdr3 = extended_cdr3 )
                assert cdr3 == old_cdr3

                logger.debug('cdr3: %s, %s, %s, %s, %s',ab,cdr3,cdr3 in qseq,'q_vframe:',q_vframe)

                status.extend(errors)

                if cdr3 in qseq:
                    if not hide_nucseq:
                        offset = qseq.find(cdr3)
                        cdr3_codons = codons[ offset:offset+len(cdr3) ]
                        cdr3 += '-'+''.join(cdr3_codons)
                        '''    
                        cdr3_nucseq = ''.join( cdr3_codons ).upper()
                        nucseq_startpos = 3*offset + q_vframe
                        alt_nucseq = blast_seq[ nucseq_startpos:nucseq_startpos+len(cdr3_nucseq) ]
                        rc1 = logo_tools.reverse_complement(cdr3_nucseq)
                        rc2 = logo_tools.reverse_complement(blast_seq)
                        logger.debug('cdr3_nucseq: %s, %s, %s, %s, %s, %s, %s',ab,offset,cdr3_nucseq,cdr3_nucseq in blast_seq,\
                            blast_seq.index(cdr3_nucseq),alt_nucseq,rc1 in rc2)
                        '''
                else:
                    assert cdr3 == '-'
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
            os.remove( file )

    assert len(genes) == 7

    if return_all_good_hits:
        return genes,evalues,status,all_good_hits_with_scores ## status is a list, may be empty
    else:
        return genes,evalues,status ## status is a list, may be empty

### most of this will only work for Blast  --  not psi-blast
def is_new_query_id_line( line ):
    return line and line.split() and line.split()[0] == 'Query='

def is_new_hit_id_line( line ):
    return line and line[0] == '>'

def is_new_match_line( line ):
    return line and line[:3] == ' Sc'

def is_query_alignment_line( line ):
    return line and line.split() and line.split()[0] == 'Query:'

class BlastMatch:
    """ Reads a single blast match from a set of input lines"""

    def get_line( s ):
        if not s.lines: return ''
        line = s.lines[0]
        del s.lines[0]
        return line


    def __init__( s, lines, query_id, hit_id ):
        s.valid = False
        s.query_id = query_id
        s.hit_id = hit_id
        s.lines = lines[:]

        line = s.get_line()
        assert is_new_match_line( line )

        assert line.split()[5][:6] == 'Expect' and line.split()[6] == '='
        evalue = line.split()[7]
        if evalue[-1] == ',': evalue = evalue[:-1]
        #assert evalue[-1] == ','
        if evalue[0] == 'e': evalue = '1'+evalue
        s.evalue = float( evalue )

        line = s.get_line()
        ident = line.split()[3]
        assert ident[0] == '('
        if ident[-3:] == '%),':
            ident = int( ident[1:-3] )
        else:
            assert ident[-2:] == '%)'
            ident = int( ident[1:-2] ) ## pb added this line 10/15/15 (!)
        s.identities = ident

        is_blastx_output = False
        is_blastn_output = False
        s.frame = 'NA'
        q_strand = 1
        h_strand = 1
        while line and not is_query_alignment_line( line ):
            line = s.get_line()
            if line.startswith(' Frame'):
                #print line[:-1]
                s.frame = line.split()[2]
                is_blastx_output = True
            if line.startswith(' Strand'):
                is_blastn_output = True
                #print line[:-1]
                l = line.split()
                assert l[3] == '/'
                assert len(l) == 5
                if l[2] == 'Plus':
                    q_strand = 1
                else:
                    assert l[2] == 'Minus'
                    q_strand = -1
                if l[4] == 'Plus':
                    h_strand = 1
                else:
                    assert l[4] == 'Minus'
                    h_strand = -1

        ## loop through possible several lines of alignment
        qstarts = []
        qstops = []
        qseq = ''
        hstarts = []
        hstops = []
        hseq = ''
        middleseq= ''
        while line:
            assert is_query_alignment_line( line )

            qstarts.append( int(line.split()[1])-1) ## 0-indexed
            qstops.append( int(line.split()[3])-1 )
            new_qseq = line.split()[2]
            qseq = qseq + new_qseq
            middleseq_offset = line.index( new_qseq )
            line = s.get_line()
            middleseq += line[ middleseq_offset: middleseq_offset+len(new_qseq)]
            line = s.get_line()
            hstarts.append( int(line.split()[1])-1)
            hstops.append( int(line.split()[3])-1)
            hseq = hseq + line.split()[2]

            ## advance to next query alignment line
            while line and not is_query_alignment_line( line ): line = s.get_line()

        qstart = q_strand * min( [q_strand*x for x in qstarts ] )
        qstop  = q_strand * max( [q_strand*x for x in qstops ] )
        hstart = h_strand * min( [h_strand*x for x in hstarts ] )
        hstop  = h_strand * max( [h_strand*x for x in hstops ] )

        if not is_blastx_output:
            assert q_strand*(qstop-qstart) == len( ''.join( qseq.split('-') ) ) - 1
            assert h_strand*(hstop-hstart) == len( ''.join( hseq.split('-') ) ) - 1

        assert len(hseq) == len(qseq)

        #print hit_id
        #print qseq
        #print hseq

        ## success
        s.h_start = hstart ## 0-indexed
        s.h_stop = hstop
        s.h_strand = h_strand
        s.h_align = hseq
        s.q_start = qstart
        s.q_stop = qstop
        s.q_strand = q_strand
        s.q_align = qseq
        s.middleseq = middleseq

        q2hmap = {}
        ia = 0
        ib = 0
        gaps = '.-'
        for i,a in enumerate( qseq ):
            b = hseq[i]
            if a not in gaps:
                if b in gaps: q2hmap[ qstart+ia ] = (-1,'-')
                else:         q2hmap[ qstart+ia ] = (hstart+ib,b)
            if a not in gaps: ia += q_strand
            if b not in gaps: ib += h_strand
        s.q2hmap = q2hmap ## 0-indexed numbering wrt to fullseq
        s.valid = True

def parse_blast_alignments( blastfile, evalue_threshold, identity_threshold ):
    ## THIS WILL NOT WORK FOR PSIBLAST

    data = open( blastfile, 'r' )

    line = data.readline()
    while line and not is_new_query_id_line( line ): line = data.readline()

    hits = {}
    while line:
        assert is_new_query_id_line( line )

        query_id = line.split()[1]
        hits[ query_id ] = []

        ## read to the start of the alignments
        line = data.readline()
        while line and not is_new_hit_id_line( line ) and not is_new_query_id_line( line ): line = data.readline()

        if not is_new_hit_id_line( line ): continue ## no hits

        while line:
            assert is_new_hit_id_line( line )

            hit_id = line.split()[0][1:]

            ## read to the first match for this hit
            while line and not is_new_match_line( line ): line = data.readline()

            while line:
                assert is_new_match_line( line )

                hitlines = [ line ]

                line = data.readline()
                while line and \
                        not is_new_match_line( line ) and \
                        not is_new_hit_id_line( line ) and \
                        not is_new_query_id_line( line ):
                    hitlines.append( line )
                    line = data.readline()

                ## now in a new match
                m = BlastMatch( hitlines, query_id, hit_id )

                if m.evalue <= evalue_threshold and m.identities >= identity_threshold:
                    hits[ query_id ].append( m )

                if not is_new_match_line( line ): break
            if not is_new_hit_id_line( line ): break
        if not is_new_query_id_line( line ):
            assert not line
            break
    data.close()

    return hits

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

def get_qualstring( cdr3seqtag, nucseq_in, quals_in ):
    nucseq = nucseq_in[:].upper()
    quals = quals_in[:]
    if '-' not in cdr3seqtag: return ''
    if len( cdr3seqtag.split('-') )!= 2: return ''
    cdr3nucseq = cdr3seqtag.split('-')[1].upper()
    if not cdr3nucseq: return ''
    if cdr3nucseq not in nucseq:
        nucseq = logo_tools.reverse_complement(nucseq)
        quals[::-1]
    if nucseq.count( cdr3nucseq ) == 0:
        # print cdr3nucseq
        # print nucseq
        # print nucseq_in
        return ''
    elif nucseq.count( cdr3nucseq ) == 1:
        pos = nucseq.index( cdr3nucseq)
    else:
        ## multiple occurrences of cdr3nucseq in nucseq, choose the one with the best quality
        max_min_qual = 0
        bestpos = 0
        offset = 0
        while nucseq.count( cdr3nucseq, offset ):
            pos = nucseq.find( cdr3nucseq, offset )
            min_qual = min( quals[pos:pos+len(cdr3nucseq)] )
            logger.debug('multiple matches: %s %s %s',pos,min_qual,max_min_qual)
            if min_qual>max_min_qual:
                max_min_qual = min_qual
                bestpos = pos
            offset = pos+1
        pos = bestpos
    return '.'.join( [ '%d'%x for x in quals[pos:pos+len(cdr3nucseq)] ] )