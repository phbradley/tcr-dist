from basic import *
from all_genes import all_genes, gap_character
from amino_acids import amino_acids
from tcr_distances_blosum import blosum, bsd4

class DistanceParams:
    def __init__(self, config_string=None ):
        self.gap_penalty_v_region = 4
        self.gap_penalty_cdr3_region = 12 # same as gap_penalty_v_region=4 since weight_cdr3_region=3 is not applied
        self.weight_v_region = 1
        self.weight_cdr3_region = 3
        self.distance_matrix = bsd4
        self.align_cdr3s = False
        self.trim_cdr3s = True
        self.scale_factor = 1.0
        if config_string:
            l = config_string.split(',')
            for tag,val in [x.split(':') for x in l ]:
                if tag == 'gap_penalty_cdr3_region':
                    self.gap_penalty_cdr3_region = float(val)
                elif tag == 'gap_penalty_v_region':
                    self.gap_penalty_v_region = float(val)
                elif tag == 'weight_cdr3_region':
                    self.weight_cdr3_region = float(val)
                elif tag == 'weight_v_region':
                    self.weight_v_region = float(val)
                elif tag == 'scale_factor':
                    self.scale_factor = float(val)
                elif tag == 'align_cdr3s':
                    assert val in ['True','False']
                    self.align_cdr3s = ( val == 'True' )
                elif tag == 'trim_cdr3s':
                    assert val in ['True','False']
                    self.trim_cdr3s = ( val == 'True' )
                else:
                    print 'unrecognized tag:',tag
                    assert False
            print 'config_string: {} self: {}'.format( config_string, self )

    def __str__(self):
        return 'DistanceParams: gap_penalty_v_region= {} gap_penalty_cdr3_region= {} weight_v_region= {} weight_cdr3_region= {} align_cdr3s= {} trim_cdr3s= {}'\
            .format( self.gap_penalty_v_region, self.gap_penalty_cdr3_region,
                     self.weight_v_region, self.weight_cdr3_region,
                     self.align_cdr3s, self.trim_cdr3s )

default_distance_params = DistanceParams()

def blosum_character_distance( a, b, gap_penalty, params ):
    if a== gap_character and b == gap_character:
        return 0
    elif a == '*' and b == '*':
        return 0
    elif a == gap_character or b == gap_character or a=='*' or b=='*':
        return gap_penalty
    else:
        # assert a in amino_acids
        # assert b in amino_acids
        # maxval = min( blosum[(a,a)], blosum[(b,b)] )
        # return maxval - blosum[(a,b)]
        return params.distance_matrix[ (a,b) ]

def blosum_sequence_distance( aseq, bseq, gap_penalty, params ):
    assert len(aseq) == len(bseq)
    dist = 0.0
    for a,b in zip(aseq,bseq):
        if a == ' ':
            assert b== ' '
        else:
            dist += blosum_character_distance( a, b, gap_penalty, params )
    return dist

def align_cdr3s( a, b, gap_character ):
    if len(a) == len(b):
        return (a[:],b[:])

    if len(a)<len(b): ## s0 is the shorter sequence
        s0,s1 = a,b
    else:
        s0,s1 = b,a

    lendiff = len(s1)-len(s0)

    best_score=-1000
    best_gappos=0 # in case len(s0) == 1

    # the gap comes after s0[gappos]

    for gappos in range(len(s0)-1):
        score=0
        for i in range(gappos+1):
            score += blosum[ (s0[i],s1[i]) ]
        for i in range(gappos+1,len(s0)):
            score += blosum[ (s0[i],s1[i+lendiff]) ]
        if score>best_score:
            best_score = score
            best_gappos = gappos
    ## insert the gap
    s0 = s0[:best_gappos+1] + gap_character*lendiff + s0[best_gappos+1:]

    assert len(s0) == len(s1)

    if len(a)<len(b): ## s0 is the shorter sequence
        return ( s0, s1 )
    else:
        return ( s1, s0 )

## align
##
##   shortseq[        ntrim: gappos   ] with longseq[ ntrim: gappos ] and
##   shortseq[ -1*remainder: -1*ctrim ] with longseq[ -1*remainder: -1*ctrim ]
##
## but be careful about negative indexing if ctrim is 0
##
## the gap comes after position (gappos-1) ie there are gappos amino acids before the gap
##
##
## DOES NOT INCLUDE THE GAP PENALTY
##
def sequence_distance_with_gappos( shortseq, longseq, gappos, params ):
    ntrim = 3 if params.trim_cdr3s else 0
    ctrim = 2 if params.trim_cdr3s else 0
    remainder = len(shortseq)-gappos
    dist = 0.0
    count =0
    if ntrim < gappos:
        for i in range(ntrim,gappos):
            #print i,shortseq[i],longseq[i],params.distance_matrix[(shortseq[i],longseq[i])]
            dist += params.distance_matrix[ (shortseq[i], longseq[i] ) ]
            count += 1
    #print 'sequence_distance_with_gappos1:',gappos,ntrim,ctrim,remainder,dist
    if ctrim < remainder:
        for i in range(ctrim, remainder):
            #print -1-i,shortseq[-1-i],longseq[-1-i],params.distance_matrix[(shortseq[-1-i],longseq[-1-i])]
            dist += params.distance_matrix[ (shortseq[-1-i], longseq[-1-i] ) ]
            count += 1
    #print 'sequence_distance_with_gappos2:',gappos,ntrim,ctrim,remainder,dist
    return dist,count


def weighted_cdr3_distance( seq1, seq2, params ):
    shortseq,longseq = (seq1,seq2) if len(seq1)<=len(seq2) else (seq2,seq1)

    ## try different positions of the gap
    lenshort = len(shortseq)
    lenlong = len(longseq)
    lendiff = lenlong - lenshort

#    assert lenshort>3 ##JCC testing
    assert lenshort > 1##JCC testing
    assert lendiff>=0
    if params.trim_cdr3s:
        assert lenshort > 3+2 ## something to align... NOTE: Minimum length of cdr3 protein carried into clones file is currently set in the read_sanger_data.py script!

    if not params.align_cdr3s:
        ## if we are not aligning, use a fixed gap position relative to the start of the CDR3
        ## that reflects the typically longer and more variable-length contributions to
        ## the CDR3 from the J than from the V. For a normal-length
        ## CDR3 this would be after the Cys+5 position (ie, gappos = 6; align 6 rsds on N-terminal side of CDR3).
        ## Use an earlier gappos if lenshort is less than 11.
        ##
        gappos = min( 6, 3 + (lenshort-5)/2 )
        best_dist,count = sequence_distance_with_gappos( shortseq, longseq, gappos, params )

    else:
        ## the CYS and the first G of the GXG are 'aligned' in the beta sheet
        ## the alignment seems to continue through roughly CYS+4
        ## ie it's hard to see how we could have an 'insertion' within that region
        ## gappos=1 would be a insertion after CYS
        ## gappos=5 would be a insertion after CYS+4 (5 rsds before the gap)
        ## the full cdr3 ends at the position before the first G
        ## so gappos of len(shortseq)-1 would be gap right before the 'G'
        ## shifting this back by 4 would be analogous to what we do on the other strand, ie len(shortseq)-1-4
        min_gappos = 5
        max_gappos = len(shortseq)-1-4
        while min_gappos>max_gappos:
            min_gappos -= 1
            max_gappos += 1
        for gappos in range( min_gappos, max_gappos+1 ):
            dist, count = sequence_distance_with_gappos( shortseq, longseq, gappos, params )
            if gappos>min_gappos:
                assert count==best_count
            if gappos == min_gappos or dist < best_dist:
                best_dist = dist
                best_gappos = gappos
                best_count = count
        #print 'align1:',shortseq[:best_gappos] + '-'*lendiff + shortseq[best_gappos:], best_gappos, best_dist
        #print 'align2:',longseq, best_gappos, best_dist


    ## Note that weight_cdr3_region is not applied to the gap penalty
    ##
    return  params.weight_cdr3_region * best_dist + lendiff * params.gap_penalty_cdr3_region

def compute_all_v_region_distances( organism, params ):
    rep_dists = {}
    for chain in 'AB': # don't compute inter-chain distances
        repseqs = []
        for id,g in all_genes[organism].iteritems():
            if g.chain == chain and g.region == 'V':
                merged_loopseq = ' '.join( g.cdrs[:-1])
                repseqs.append( ( id, merged_loopseq  ) )
                rep_dists[ id ] = {}

        for r1,s1 in repseqs:
            for r2,s2 in repseqs:
                #if r1[2] != r2[2]: continue
                rep_dists[r1][r2] = params.weight_v_region * \
                                    blosum_sequence_distance( s1, s2, params.gap_penalty_v_region, params )

    return rep_dists

def compute_distance(t1,t2,chains,rep_dists,distance_params): #    t1/2 = [ va_reps, vb_reps, l['cdr3a'], l['cdr3b'] ]
    dist=0.0
    if 'A' in chains:
        dist += min( ( rep_dists[x][y] for x in t1[0] for y in t2[0] ) ) +\
                weighted_cdr3_distance( t1[2], t2[2], distance_params )
    if 'B' in chains:
        dist += min( ( rep_dists[x][y] for x in t1[1] for y in t2[1] ) ) +\
                weighted_cdr3_distance( t1[3], t2[3], distance_params )
    return distance_params.scale_factor * dist


def compute_auc( l0, l1, sign_factor=1 ):
    ## l0 are the true positives, l1 are the false positives
    ## if sign_factor==1 then lower scores are better, otherwise it's the opposite
    ##
    if not l0:
        return 0.0, [0,1], [0,0]
    elif not l1:
        return 1.0, [0,0,1], [0,1,1]

    l = [ (sign_factor*x,0) for x in l0 ] + [ (sign_factor*x,-1) for x in l1 ] ## in ties, take the false positive first
    l.sort()

    xvals = []
    yvals = []

    counts = [0,0]
    totals = [len(l0),len(l1)]

    area=0.0
    width = 1.0/totals[1]
    for ( score, neg_tcr_class ) in l:
        tcr_class = -1*neg_tcr_class
        counts[ tcr_class ] += 1
        xval = float( counts[1] ) / totals[1]
        yval = float( counts[0] ) / totals[0]
        xvals.append( xval )
        yvals.append( yval )
        if tcr_class==1: area += yval * width

    return area,xvals,yvals

def get_rank( val, l ): ## does not require that the list l is sorted
    num_lower = 0
    num_upper = 0

    epsilon = 1e-6

    lower_neighbor = val-10000
    upper_neighbor = val+10000

    for x in l:
        if x<val-epsilon:
            num_lower += 1
            lower_neighbor = max( lower_neighbor, x )
        elif x>val+epsilon:
            num_upper += 1
            upper_neighbor = min( upper_neighbor, x )

    total = len(l)
    num_equal = total - num_lower - num_upper
    assert num_equal >=0

    if num_upper == 0:
        return 100.0

    elif num_lower == 0:
        return 0.0

    else:
        assert upper_neighbor>lower_neighbor
        interp = (val-lower_neighbor)/(upper_neighbor-lower_neighbor)

        #if num_equal>0:print 'num_equal:',num_equal

        interp_num_lower = num_lower + interp * ( 1 + num_equal )

        return (100.0*interp_num_lower)/total



## negative nbrdist_percentile means take exactly -nbrdist_percentile topn
def sort_and_compute_nbrdist_from_distances( l, nbrdist_percentile, dont_sort=False ):
    if not dont_sort: l.sort()
    assert l[0]<=l[-1]
    if nbrdist_percentile<0:
        n = max( 1, min(len(l), -1*nbrdist_percentile ) )
    else:
        n = max(1, ( nbrdist_percentile * len(l) )/100 )
    return sum( l[:n])/float(n)

## negative nbrdist_percentile means take exactly -nbrdist_percentile topn
def sort_and_compute_weighted_nbrdist_from_distances( l, nbrdist_percentile, dont_sort=False ):
    if not dont_sort: l.sort()
    assert l[0]<=l[-1]
    if nbrdist_percentile<0:
        n = max( 1, min(len(l), -1*nbrdist_percentile ) )
    else:
        n = max(1, ( nbrdist_percentile * len(l) )/100 )

    total_wt = 0.0
    nbrdist=0.0
    for i,val in enumerate( l[:n] ):
        wt = 1.0 - float(i)/n
        total_wt += wt
        nbrdist += wt * val

    return nbrdist / total_wt



# if __name__ == '__main__': ## hacking
#     a,b = "DVGYKL DPAGNTGKL".split()
#     #a,b = "GEGSNNRI GYNTNTGKL".split()
#     #a,b = "GDRYAQGL GDVDYAQGL".split()
#     print align_cdr3s( a,b,'.')
#     exit()


if __name__ == '__main__':
    # generate an input file for tcr-dist calculation in C++
    #
    params = DistanceParams()

    for aa in amino_acids:
        print 'AAdist',aa,
        for bb in amino_acids:
            print '{:.3f}'.format( bsd4[(aa,bb)] ),
        print

    rep_dists = compute_all_v_region_distances( 'human', params )

    ## this part only works with the classic db (getting chain from id[2] is bad for gammadelta)
    vb_genes = [ x for x in rep_dists.keys() if x[2] == 'B' ]
    vb_genes.sort()

    print 'num_v_genes',len(vb_genes)

    for v1 in vb_genes:
        print 'Vdist',v1,
        for v2 in vb_genes:
            print '{:.3f}'.format(rep_dists[v1][v2]),
        print
