import util
import tcr_distances
import copy

class DistanceTCR:
    ## cdr3s are expected to run from conserved Cys to F/Y/W position right before the 'GXG' motif (inclusive of C+F)
    def __init__( self, organism, va_genes, vb_genes, cdr3a, cdr3b ):
        self.va_reps = frozenset( ( util.get_rep( x, organism ) for x in va_genes ) )
        self.vb_reps = frozenset( ( util.get_rep( x, organism ) for x in vb_genes ) )
        self.cdr3a = cdr3a[:]
        self.cdr3b = cdr3b[:]

class DistanceCalculator:
    def __init__( self, organism, params = None ):
        self.organism = organism[:]
        if params is None:
            self.params = tcr_distances.DistanceParams()
        else:
            self.params = copy.deepcopy( params )
        self.rep_dists = tcr_distances.compute_all_v_region_distances( organism, self.params )

    ## set chains == 'A' (or 'B') for single-chain alpha (or beta) distances
    ## tcr1 and tcr2 are DistanceTCR objects
    def distance( self, tcr1, tcr2, chains='AB' ):
        dist=0.0
        if 'A' in chains:
            dist += min( ( self.rep_dists[x][y] for x in tcr1.va_reps for y in tcr2.va_reps ) ) +\
                    tcr_distances.weighted_cdr3_distance( tcr1.cdr3a, tcr2.cdr3a, self.params )
        if 'B' in chains:
            dist += min( ( self.rep_dists[x][y] for x in tcr1.vb_reps for y in tcr2.vb_reps ) ) +\
                    tcr_distances.weighted_cdr3_distance( tcr1.cdr3b, tcr2.cdr3b, self.params )

        return self.params.scale_factor * dist

if __name__ == '__main__':
    lines = """TRBV18*01 CASSPRHGISPLHF
TRBV19*01,TRBV19*02 CASSPGGVTEAFF
TRBV28*01 CASRGPGEPYEQYF
TRBV4-1*01,TRBV4-1*02 CASSQAEGGLLSYEQYF
TRBV6-6*01,TRBV6-6*02,TRBV6-6*03,TRBV6-6*04 CASMLGGPPQETQYF
TRBV7-8*03 CASSHSWDSGTGELFF
TRBV7-9*01,TRBV7-9*02,TRBV7-9*03 CASSSSGGAYNEQFF
TRBV19*01,TRBV19*02 CASSIEGSGANVLTF
TRBV7-8*01 CASSSSGGMNIQYF""".split('\n')


    lines = """TRBV6-1*01 CASSEAPLLGGNEQYF
TRBV29-1*01 CSVEARLF""".split('\n')

    organism = 'human'

    tcrs = []
    for line in lines:
        l = line.split()
        vb_genes = l[0].split(',')
        cdr3b = l[1]
        ## bogus place-holders
        va_genes = []
        cdr3a = ''
        tcrs.append( DistanceTCR( organism, va_genes, vb_genes, cdr3a, cdr3b ) )

    calculator = DistanceCalculator( organism )

    for i,t1 in enumerate(tcrs):
        for j,t2 in enumerate(tcrs):
            if j<=i:continue
            dist = calculator.distance( t1, t2, chains = 'B' )
            print t1.vb_reps, t1.cdr3b, t2.vb_reps, t2.cdr3b, dist
            #if dist<100:

