

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

