from basic import *
import util
from scipy import stats
from tcr_distances import get_rank
from all_genes import all_genes

import numpy as np

distance_threshold_default = pipeline_params['distance_threshold_25']


with Parser(locals()) as p:
    p.str('clones_file').required()
    p.multiword('epitopes').cast(lambda x:x.split())
    p.str('organism').required()
    p.str('outfile_prefix')
    p.int('nbrdist_percentile').default(10)
    p.flag('verbose')
    p.flag('show')
    p.float('distance_threshold').default( distance_threshold_default )
    p.flag('unweighted_nbrdist')

wtd_nbrdist = not unweighted_nbrdist #silly

if outfile_prefix is None:
    outfile_prefix = clones_file[:-4]

outlogfile = '{}_sharing.log'.format( clones_file[:-4] )
print 'making:',outlogfile
outlog =open( outlogfile,'w')

fake_chains = util.detect_fake_chains( clones_file )

import matplotlib
if not show: matplotlib.use('Agg')
import matplotlib.pyplot as plt

def confidence_interval( k, N, interval_fraction ):

    def get_prob( logp,k=k,N=N ):
        return stats.binom.pmf(k,N,math.exp(logp) )

    p0 = float(k)/N
    logp0 = math.log(p0)

    stepsize = 0.01

    #prob0 = stats.binom.pmf( k, N, p0 )
    #fwdcum = prob0
    #revcum = prob0

    upper = logp0
    lower = logp0

    upperprob = get_prob(upper)
    lowerprob = get_prob(lower)

    cumulative = 0.0

    nextlowerprob = get_prob( lower-stepsize )
    nextupperprob = get_prob( upper+stepsize )

    bounds = []
    while True:
        last_cumulative = cumulative
        if nextlowerprob > nextupperprob:
            cumulative += 0.5 * ( lowerprob + nextlowerprob )
            lower -= stepsize
            lowerprob = nextlowerprob
            nextlowerprob = get_prob( lower-stepsize )
        else:
            cumulative += 0.5 * ( upperprob + nextupperprob )
            upper += stepsize
            upperprob = nextupperprob
            nextupperprob = get_prob( upper+stepsize )
        #print last_cumulative, cumulative, lower,upper
        bounds.append( ( cumulative, lower, upper ) )
        if cumulative-last_cumulative < 0.001 * last_cumulative:
            break
    total = cumulative

    for (cum,lower,upper) in bounds:
        if cum/total > interval_fraction:
            break

    interval = ( math.exp(lower), math.exp(upper) )
    #print 'frac',cum/total, len(bounds), p0, interval

    return p0, interval, cum/total


#confidence_interval( 20, 10000, 0.95 )



#exit()




all_tcrs = {}
all_info = []

infields = []

all_protprobs = {}
all_nucprobs = {}
for ab in ['A','B','AB']:
    all_protprobs[ab] = []
    all_nucprobs[ab] = []

clones_file_with_nbrdists = '{}_nbrdists.tsv'.format( clones_file[:-4] )
assert exists( clones_file_with_nbrdists )

for line in open( clones_file_with_nbrdists,'r'):
    if not infields:
        if line[0] == '#':
            infields = line[1:-1].split('\t')
        else:
            infields = line[:-1].split('\t')
        continue
    assert infields

    l = parse_tsv_line( line[:-1], infields )

    mouse = l['subject']
    epitope = l['epitope']

    if epitopes and epitope not in epitopes: continue

    ## we probably should be using the 'va_genes' and 'ja_genes' info in the tsv line
    ## which comes from the clone finding process... oh, well. stick with this for now.
    ##
    va_genes = set( l['va_genes'].split(';') )
    ja_genes = set( l['ja_genes'].split(';') )
    vb_genes = set( l['vb_genes'].split(';') )
    jb_genes = set( l['jb_genes'].split(';') )

    va_reps = set(( all_genes[organism][x].rep for x in va_genes ))
    ja_reps = set(( all_genes[organism][x].rep for x in ja_genes ))
    vb_reps = set(( all_genes[organism][x].rep for x in vb_genes ))
    jb_reps = set(( all_genes[organism][x].rep for x in jb_genes ))

    protprob = { 'A': float(l['a_protseq_prob']),
                 'B': float(l['b_protseq_prob']),
                 'AB': float(l['a_protseq_prob']) * float( l['b_protseq_prob'] ) }

    nucprob  = { 'A': float(l[ 'a_nucseq_prob']),
                 'B': float(l[ 'b_nucseq_prob']),
                 'AB': float(l[ 'a_nucseq_prob']) * float( l[ 'b_nucseq_prob'] ) }

    for ab in protprob:
        all_protprobs[ab].append( protprob[ab] )
        all_nucprobs [ab].append(  nucprob[ab] )

    #clone_id = l['clone_id']
    line_index = len(all_info)
    all_info.append( dict(l))

    tcr = [ va_reps, ja_reps, vb_reps, jb_reps, l['cdr3a'], l['cdr3b'], line_index, protprob ]

    nuctcr = [ va_genes, ja_genes, vb_genes, jb_genes, l['cdr3a_nucseq'], l['cdr3b_nucseq'], line_index, nucprob ]


    if epitope not in all_tcrs:
        all_tcrs[epitope] = {}
    if mouse not in all_tcrs[epitope]:
        all_tcrs[epitope][mouse] = []
    all_tcrs[epitope][mouse].append( [tcr,nuctcr] )


def same_tcr( a,b,chains,comparison_mode):
    if comparison_mode == 2:## distance-based
        global all_chain_dists
        dist = all_chain_dists[ chains ][ a[6] ][ b[6] ]
        return ( dist <= distance_threshold * len(chains) )

    else:
        ## comparison_mode==0: full comparison
        ## comparison_mode==1: ignore CDR3s, ie only look at gene segments
        if 'A' in chains and ( (comparison_mode==0 and a[4] != b[4]) or \
                               a[0].isdisjoint(b[0]) or a[1].isdisjoint(b[1]) ): return False
        if 'B' in chains and ( (comparison_mode==0 and a[5] != b[5]) or \
                               a[2].isdisjoint(b[2]) or a[3].isdisjoint(b[3]) ): return False
    return True


smallest_nonzero_nucprob = {}
smallest_nonzero_protprob = {}

for ab in all_nucprobs:
    smallest_nonzero_nucprob[ab] = min( ( x for x in all_nucprobs[ab] if x>0 ) )
    smallest_nonzero_protprob[ab] = min( ( x for x in all_protprobs[ab] if x>0 ) )


## update the protprobs
for epitope in all_tcrs:
    for mouse,tcrs in all_tcrs[epitope].iteritems():
        for ptcr,ntcr in tcrs:
            for chains in ['A','B','AB']:
                if ptcr[-1][chains]==0: ptcr[-1][chains] = smallest_nonzero_protprob[chains]
                if ntcr[-1][chains]==0: ntcr[-1][chains] = smallest_nonzero_nucprob [chains]


## look at clonality vs probs or rank score, and sharing vs probs/rank/clonality
## so need to cluster tcrs to identify shared ones

#for skip_epitope in ['NONE'] + all_tcrs.keys():
if True:
    skip_epitope =''

    ## now adding a third list for _c arrays: is_clonal=2 if clone_size==1 and clone_size==max_clone_size
    ##
    nucprobs_c    = [ [], [], [] ] ## indexed by is-clonal
    protprobs_c   = [ [], [], [] ] ## indexed by is-clonal
    protprobs_s   = [ [], [] ] ## indexed by is-shared
    protprob_ranks_c   = [ [], [], [] ] ## indexed by is-clonal
    protprob_ranks_s   = [ [], [] ] ## indexed by is-shared
    nbrdist_rank_scores_s = [ [], [] ]
    nbrdist_rank_scores_c = [ [], [], [] ]
    nbrdist_scores_s = [ [], [] ]
    nbrdist_scores_c = [ [], [], [] ]

    rank_suffix = '_nbrdist{}rank'.format(nbrdist_percentile)
    nbrdist_suffix = '_nbrdist'+str(nbrdist_percentile)

    if wtd_nbrdist:
        rank_suffix = '_wtd'+rank_suffix
        nbrdist_suffix = '_wtd'+nbrdist_suffix


    table_c = [ [0,0], [0,0] ]
    table_s = [ [0,0], [0,0] ]

    all_nbrdist_rank_scores = []
    all_nbrdist_scores = []
    all_protprobs = []
    all_protprob_ranks = []
    all_nucprobs = []

    for epitope in all_tcrs:
        if epitope==skip_epitope: continue
        #
        epitope_protprobs = []
        for mouse in all_tcrs[epitope]:
            epitope_protprobs.extend( [ math.log10( ptcr[-1]['AB'] ) for ptcr,ntcr in all_tcrs[epitope][mouse] ] )

        #print epitope, type(epitope_protprobs[0]), epitope_protprobs[:3]

        for mouse,tcrs in all_tcrs[epitope].iteritems():
            max_clone_size = max( ( int( all_info[x[0][6]]['clone_size'] ) for x in tcrs ) )
            for ptcr,ntcr in tcrs:
                ## is this clone shared across other mice?
                info = all_info[ ptcr[6] ]
                clone_size = int( info['clone_size'] )

                ## define is_clonal and is_shared
                #is_clonal = 2 if ( clone_size==1 and clone_size==max_clone_size) else 1 if (clone_size>1 ) else 0
                is_clonal = 1 if (clone_size>1 ) else 0
                is_shared = 0
                for other_mouse, other_tcrs in all_tcrs[epitope].iteritems():
                    if other_mouse==mouse: continue
                    for ptcr2,ntcr2 in other_tcrs:
                        if same_tcr( ptcr, ptcr2, 'AB', comparison_mode=0 ):
                            is_shared = 1
                            break
                    if is_shared: break

                ## now the scores:
                nucprob  = math.log10( ntcr[-1]['AB'] )
                protprob = math.log10( ptcr[-1]['AB'] )
                pp_rank = get_rank( protprob, epitope_protprobs )

                nbrdist_rank_score = float( info['{}_AB{}'.format(epitope,rank_suffix)] )
                nbrdist_score = float( info['{}_AB{}'.format(epitope,nbrdist_suffix)] )

                nucprobs_c      [ is_clonal ].append( nucprob )
                protprobs_c     [ is_clonal ].append( protprob )
                nbrdist_rank_scores_c   [ is_clonal ].append( nbrdist_rank_score )
                nbrdist_scores_c[ is_clonal ].append( nbrdist_score )
                protprob_ranks_c[ is_clonal ].append( pp_rank )

                protprobs_s     [ is_shared ].append( protprob )
                nbrdist_rank_scores_s   [ is_shared ].append( nbrdist_rank_score )
                nbrdist_scores_s[ is_shared ].append( nbrdist_score )
                protprob_ranks_s[ is_shared ].append( pp_rank )

                if is_clonal<2:
                    table_c[ is_clonal ][ is_shared ] += 1
                    table_s[ is_shared ][ is_clonal ] += 1

                all_protprobs.append( protprob )
                all_nucprobs.append( nucprob )
                all_nbrdist_rank_scores.append( nbrdist_rank_score )
                all_nbrdist_scores.append( nbrdist_score )
                all_protprob_ranks.append( pp_rank )



    assert len(protprob_ranks_c[0]) == len(protprobs_c[0])
    assert len(protprob_ranks_c[1]) == len(protprobs_c[1])
    assert len(protprob_ranks_s[0]) == len(protprobs_s[0])
    assert len(protprob_ranks_s[1]) == len(protprobs_s[1])

    ## look for correlations between real-valued guys
    tagl = [ ( 'protprob', all_protprobs ),
             ( 'nucprob', all_nucprobs ),
             ( 'protprob_rank', all_protprob_ranks ),
             ( 'nbrdist_score', all_nbrdist_scores ),
             ( 'nbrdist_rank_score', all_nbrdist_rank_scores ) ]
    for xtag,xvals in tagl:
        for ytag,yvals in tagl:
            if ytag <= xtag: continue
            slope, intercept, r_value, p_value, std_err = stats.linregress(xvals,yvals)
            rho, p_value2 = stats.spearmanr( xvals,yvals)
            tau, p_value3 = stats.kendalltau( xvals,yvals)
            print 'overall_correlations: R: {:.4f} {:3g} rho: {:.4f} {:.3g} tau: {:.4f} {:.3g} {} {}'\
                .format( r_value, p_value,
                         rho, p_value2,
                         tau, p_value3,
                         xtag, ytag )



    #print table
    #print len(nucprobs[0]), len(nucprobs[1]), nucprobs[0][0]
    #print len(nbrdist_rank_scores[0]), len(nbrdist_rank_scores[1]), nbrdist_rank_scores[0][0]


    plt.figure(1,figsize=(14,6))

    nrows = 1
    ncols = 7
    plotno=0

    all_plots = [ [ 'is shared', 'protprob',protprobs_s],
                  [ 'is clonal', 'protprob',protprobs_c],
                  [ 'is clonal',  'nucprob', nucprobs_c],
                  [ 'is shared', 'nbrdist_rank_score',nbrdist_rank_scores_s],
                  [ 'is clonal', 'nbrdist_rank_score',nbrdist_rank_scores_c] ]


    for xlabel,ylabel,dats in all_plots:
        plotno += 1
        plt.subplot(nrows,ncols,plotno)

        plt.boxplot( dats )
        if len(dats[0])==0 or len(dats[1])==0:
            p=1
        else:
            t, p = stats.ttest_ind( dats[0], dats[1] )
        plt.title('P-value\n{:.3g}'.format( p ) )
        mn0 = sum(dats[0])/len(dats[0]) if dats[0] else 0.0
        mn1 = sum(dats[1])/len(dats[1]) if dats[1] else 0.0
        print 'overall_P {:9.3g} mn0 {:7.3f} mn1 {:7.3f} {} vs {}'\
            .format( p, mn0, mn1, '_'.join( xlabel.split()), '_'.join( ylabel.split()) )

        plt.xticks( [1,2], ['no\n({})'.format(len(dats[0])),'yes\n({})'.format(len(dats[1]))], fontsize=8 )
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)


    ## now show sharing, clonal frequencies
    for xlabel,ylabel,table in [ [ 'is shared','P( is clonal )',table_s ],
                                 [ 'is clonal','P( is shared )',table_c ] ]:

        oddsratio, p_table = stats.fisher_exact( table )
        plotno += 1
        plt.subplot(nrows,ncols,plotno)

        plist =[]
        totals = []
        for ii in range(2):
            total = table[ii][0] + table[ii][1]
            plist.append( 0 if table[ii][1] == 0 else float( table[ ii ][1] )/ total )
            totals.append( total )

        plt.bar( [0,1], [plist[0],plist[1]] )
        plt.xticks( [0.4,1.4], ['no\n({})'.format(totals[0]),'yes\n({})'.format(totals[1])], fontsize=8 )
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title('P-value\n{:.3g}'.format(p_table))
        print 'overall_P {:.3g} {} vs {}'.format( p_table, '_'.join( xlabel.split()), '_'.join( ylabel.split()) )


    plt.subplots_adjust(bottom = 0.1, top=0.9, right=0.98, left=0.07, wspace=0.65 )

    filetag = '{}_nbrdist{}'.format( '_wtd' if wtd_nbrdist else '', str(nbrdist_percentile) )
    pngfile = '{}_sharing_and_clonality{}.png'.format(outfile_prefix,filetag)
    print 'making',pngfile
    plt.savefig(pngfile)
    util.readme(pngfile,"""These plots explore the relationship between clonality and sharing of TCRs across mice for the same epitope. For the purpose of
    this analysis a TCR is "clonal" if it has a clone_size of at least 2 and is "shared" if it is seen in more than one mouse (ie subject). protprob and nucprob are the
    amino acid and nucleotide (respectively) generation probabilities under a very simple model of the rearrangement process. nbrdist_rank_score is a measure of
    repertoire sampling density nearby a given TCR: we compute an average distance to a TCRs nearest neighbors ("{}") and then percentile this over the repertoire to
    get a normalized nbr-distance measure that goes from 0 (many nearby TCRs in the repertoire) to 100 (very few).
    """.format(rank_suffix))

    if epitopes is None:
        epitopes = all_tcrs.keys()[:]
        epitopes.sort()

    plt.figure(2,figsize=(14,4*len(epitopes)))


    nrows = len(epitopes)
    ncols = 7

    plotno=0

    for desired_epitope in epitopes:

        ##
        nucprobs_c    = [ [], [] ] ## indexed by is-clonal
        protprobs_c   = [ [], [] ] ## indexed by is-clonal
        protprobs_s   = [ [], [] ] ## indexed by is-shared
        nbrdist_rank_scores_s = [ [], [] ]
        nbrdist_rank_scores_c = [ [], [] ]
        nbrdist_scores_s = [ [], [] ]
        nbrdist_scores_c = [ [], [] ]

        rank_suffix = '_nbrdist{}rank'.format(nbrdist_percentile)
        nbrdist_suffix = '_nbrdist'+str(nbrdist_percentile)

        if wtd_nbrdist:
            rank_suffix = '_wtd'+rank_suffix
            nbrdist_suffix = '_wtd'+nbrdist_suffix

        table_c = [ [0,0], [0,0] ]
        table_s = [ [0,0], [0,0] ]

        all_nbrdist_rank_scores = []
        all_protprobs = []
        all_nucprobs = []

        for epitope in all_tcrs:
            if epitope!=desired_epitope: continue
            #
            for mouse,tcrs in all_tcrs[epitope].iteritems():
                for ptcr,ntcr in tcrs:
                    ## is this clone shared across other mice?
                    info = all_info[ ptcr[6] ]
                    clone_size = int( info['clone_size'] )

                    ## define is_clonal and is_shared
                    is_clonal = 1 if ( clone_size>1 ) else 0
                    is_shared = 0
                    for other_mouse, other_tcrs in all_tcrs[epitope].iteritems():
                        if other_mouse==mouse: continue
                        for ptcr2,ntcr2 in other_tcrs:
                            if same_tcr( ptcr, ptcr2, 'AB', comparison_mode=0 ):
                                is_shared = 1
                                break
                        if is_shared: break

                    ## now the scores:
                    #nucprob = ptcr[-1]['AB']
                    nucprob = ntcr[-1]['AB']
                    protprob = ptcr[-1]['AB']

                    nbrdist_rank_score = float( info['{}_AB{}'.format(epitope,rank_suffix)] )
                    nbrdist_score = float( info['{}_AB{}'.format(epitope,nbrdist_suffix)] )

                    nucprobs_c   [ is_clonal ].append( math.log10(nucprob) )
                    protprobs_c  [ is_clonal ].append( math.log10(protprob) )
                    nbrdist_rank_scores_c[ is_clonal ].append( nbrdist_rank_score )
                    nbrdist_scores_c[ is_clonal ].append( nbrdist_score )

                    protprobs_s  [ is_shared ].append( math.log10(protprob) )
                    nbrdist_rank_scores_s[ is_shared ].append( nbrdist_rank_score )
                    nbrdist_scores_s[ is_shared ].append( nbrdist_score )

                    table_c[ is_clonal ][ is_shared ] += 1
                    table_s[ is_shared ][ is_clonal ] += 1

                    all_protprobs.append( math.log10(protprob))
                    all_nucprobs.append( math.log10(nucprob))
                    all_nbrdist_rank_scores.append( nbrdist_rank_score )

        ## look for correlations between real-valued guys
        tagl = [ ( 'protprob', all_protprobs ),
                 ( 'nucprob', all_nucprobs ),
                 ('nbrdist_rank_score', all_nbrdist_rank_scores ) ]
        for xtag,xvals in tagl:
            for ytag,yvals in tagl:
                if ytag <= xtag: continue
                slope, intercept, r_value, p_value, std_err = stats.linregress(xvals,yvals)
                rho, p_value2 = stats.spearmanr( xvals,yvals)
                tau, p_value3 = stats.kendalltau( xvals,yvals)
                print 'epitope_correlations: R: {:.4f} {:3g} rho: {:.4f} {:.3g} tau: {:.4f} {:.3g} {} {} {}'\
                    .format( r_value, p_value,
                             rho, p_value2,
                             tau, p_value3,
                             xtag, ytag, desired_epitope )



        for xlabel,ylabel,dats in [ [ 'is shared', 'protprob',protprobs_s],
                                    [ 'is clonal', 'protprob',protprobs_c],
                                    [ 'is clonal',  'nucprob', nucprobs_c],
                                    [ 'is shared', 'nbrdist_rank_score',nbrdist_rank_scores_s],
                                    [ 'is clonal', 'nbrdist_rank_score',nbrdist_rank_scores_c] ]:
                                    # [ 'is shared', 'nbrdist_score',nbrdist_scores_s],
                                    # [ 'is clonal', 'nbrdist_score',nbrdist_scores_c] ]:
            plotno += 1
            plt.subplot(nrows,ncols,plotno)

            plt.boxplot( dats )
            if len(dats[0])==0 or len(dats[1])==0:
                p=1
            else:
                t, p = stats.ttest_ind( dats[0], dats[1] )
            plt.xticks( [1,2], ['no\n({})'.format(len(dats[0])),'yes\n({})'.format(len(dats[1]))], fontsize=8 )
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title('{} P-value\n{:.3g}'.format( desired_epitope, p ) )
            print 'epitope_P {:.3g} {} vs {} {}'.format( p, '_'.join( xlabel.split()), '_'.join( ylabel.split()),
                                                         desired_epitope )

        ## now show sharing, clonal frequencies
        for xlabel,ylabel,table in [ [ 'is shared','P( is clonal )',table_s ],
                                     [ 'is clonal','P( is shared )',table_c ] ]:

            oddsratio, p_table = stats.fisher_exact( table )
            plotno += 1
            plt.subplot(nrows,ncols,plotno)

            plist =[]
            totals = []
            for ii in range(2):
                total = table[ii][0] + table[ii][1]
                plist.append( 0 if table[ii][1] == 0 else float( table[ ii ][1] )/ total )
                totals.append( total )

            plt.bar( [0,1], [plist[0],plist[1]] )
            plt.xticks( [0.4,1.4], ['no\n({})'.format(totals[0]),'yes\n({})'.format(totals[1])], fontsize=8 )
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title('{} P-value\n{:.3g}'.format(desired_epitope,p_table))
            print 'epitope_P {:.3g} {} vs {} {}'.format( p_table, '_'.join( xlabel.split()), '_'.join( ylabel.split()),
                                                         desired_epitope )


    plt.subplots_adjust(bottom = 0.1, top=0.9, right=0.98, left=0.07, wspace=0.65, hspace=0.4)

    filetag = '{}_nbrdist{}'.format( '_wtd' if wtd_nbrdist else '', str(nbrdist_percentile) )
    pngfile = '{}_sharing_and_clonality_by_epitope{}.png'.format(outfile_prefix,filetag)
    print 'making',pngfile
    plt.savefig(pngfile)
    util.readme(pngfile,"""Same plots as above, but now broken down by epitope.
    """)




    if show:
        plt.show()
#exit()




## analyze clonality
print "Mouse\tEpitope\tTotal_tcrs_this_mouse\tSame_pairs_this_mouse\tTotal_pairs_this_mouse\tp_this_mouse\tClone_sizes"
print "total_tcrs_this_mouse = sum(clone_sizes) and clone_sizes = [ int(all_info[x[0][6]] ['clone_size']) for x in tcrs ]"
print "same_pairs_this_mouse = sum( ( x*(x-1) for x in clone_sizes ) )"
print "total_pairs_this_mouse = total_tcrs_this_mouse * (total_tcrs_this_mouse-1)"
print "p_this_mouse = float( same_pairs_this_mouse )/total_pairs_this_mouse"

for epitope in all_tcrs:
    total_pairs = 0
    same_pairs = 0
    mice_infostrings = []
    for mouse,tcrs in all_tcrs[epitope].iteritems():
        clone_sizes = [ int(all_info[x[0][6]] ['clone_size']) for x in tcrs ]
        total_tcrs_this_mouse = sum(clone_sizes)
        same_pairs_this_mouse = sum( ( x*(x-1) for x in clone_sizes ) )
        total_pairs_this_mouse = total_tcrs_this_mouse * (total_tcrs_this_mouse-1)
        same_pairs  += same_pairs_this_mouse
        total_pairs += total_pairs_this_mouse
        if total_tcrs_this_mouse>1:
            p_this_mouse = float( same_pairs_this_mouse )/total_pairs_this_mouse
            mice_infostrings.append( '{:.6f},{}'.format(p_this_mouse,total_tcrs_this_mouse))
        else:
            p_this_mouse = 0
        print mouse, epitope, total_tcrs_this_mouse, same_pairs_this_mouse, total_pairs_this_mouse, p_this_mouse, clone_sizes
    p = float( same_pairs )/total_pairs

    #inv_p = 1.0/p if p!=0 else 1000.0
    if same_pairs:
        p0,interval,fraction = confidence_interval( same_pairs, total_pairs, 0.95 )
    else:
        interval = [p,p]
    def safe_inverse(p):
        return 1.0/p if p else 0.0
    outlog.write( 'clone_diversity: {:5s} {:9.3f} {:9.3f} {:9.3f} {:9.3f} {:9.3f} {}\n'\
                  .format( epitope,
                           safe_inverse(p), safe_inverse(interval[1]), safe_inverse(interval[0]),
                           p, 1.0-p,
                           ';'.join(mice_infostrings)))

##load distance matrix
all_chain_dists = {}

total_lines = len(all_info)

for chains in ['A','B','AB']:
    distfile = '{}_{}.dist'.format(clones_file[:-4],chains)
    Log('reading '+distfile)
    assert exists(distfile)
    all_dists = []
    for line in open( distfile,'r'):
        l = line.split()
        clone_id = l[0]
        index = len(all_dists)
        assert all_info[ index ]['clone_id'] == clone_id
        dists = [ float(x) for x in l[1:] ]
        assert total_lines == len(dists)
        all_dists.append( dists )
    all_chain_dists[chains] = all_dists


## look at epitope diversity using a gaussian-weighted "overlap" measure
for epitope1 in all_tcrs:
    for epitope2 in all_tcrs:
        if epitope2<epitope1: continue
        for same_mouse in [True,False]:
            ## if same_mouse==True and epitope1==epitope2 it means that we include same-mouse distances
            ## if same_mouse==True and epitope1!=epitope2 it means that we only use same-mouse distances
            for chains in ['A','B','AB']:
                dist_sdev = distance_threshold * len(chains)

                overlap_sum = 0.0
                total_sum = 0.0

                for m1,tcrs1 in all_tcrs[epitope1].iteritems():
                    for m2,tcrs2 in all_tcrs[epitope2].iteritems():
                        if epitope1==epitope2:
                            if m2<m1:continue # dont count pairs twice
                            if (not same_mouse) and m1==m2: continue
                        else:
                            if same_mouse != (m1==m2): continue

                        same_mouse_same_epitope = ( (m1==m2) and ( epitope1==epitope2 ) )

                        for it1,t1 in enumerate(tcrs1):
                            for it2,t2 in enumerate(tcrs2):
                                if same_mouse_same_epitope and it2<=it1: continue
                                dist = all_chain_dists[ chains ][ t1[0][6] ][ t2[0][6] ]
                                overlap_sum += math.exp( -1.0 * (dist/dist_sdev)**2 )
                                total_sum += 1.0
                if total_sum:
                    p0 = overlap_sum / total_sum
                    diversity = 1.0/p0
                else:
                    p0 = 0.0
                    diversity = 0.0

                N1 = sum( ( len(x) for x in all_tcrs[epitope1].values() ) )
                N2 = sum( ( len(x) for x in all_tcrs[epitope2].values() ) )

                outlog.write( 'GAUSSDIV SM{:d} SE{:d} {:{}s} {:{}s} {:2s} div: {:9.1f} overlap: {:.3f} total: {:.1f} N1: {} N2: {}\n'\
                              .format( same_mouse, epitope1 == epitope2,
                                       epitope1, max((len(x) for x in all_tcrs)),
                                       epitope2, max((len(x) for x in all_tcrs)),
                                       chains, diversity, overlap_sum, total_sum,
                                       N1, N2 ) )
#exit()

## look at epitope diversity using a gaussian-weighted shannon's entropy measure
for epitope in all_tcrs:
    for chains in ['A','B','AB']:
        dist_sdev = distance_threshold * len(chains)

        entropy = 0.0
        count = 0
        for m1,tcrs1 in all_tcrs[epitope].iteritems():
            for it1,t1 in enumerate(tcrs1):
                my_pval_sum = 0.0
                my_pval_norm = 0.0

                for m2,tcrs2 in all_tcrs[epitope].iteritems():
                    for it2,t2 in enumerate(tcrs2):
                        if m1==m2 and it1==it2: continue ## no self-distance
                        dist = all_chain_dists[ chains ][ t1[0][6] ][ t2[0][6] ]
                        my_pval_sum += math.exp( -1.0 * (dist/dist_sdev)**2 )
                        my_pval_norm += 1

                entropy -= math.log( my_pval_sum/my_pval_norm )
                count += 1
        entropy /= count
        diversity = 2**entropy

        outlog.write( 'GAUSSDIVSHANNON {:{}s} {:2s} div: {:12.3f} entropy: {:.3f} N: {}\n'\
                      .format( epitope, max((len(x) for x in all_tcrs)),
                               chains, diversity, entropy, count ))

## look at epitope diversity using wtdnbrdist scores
nbrdist_suffix = '_wtd_nbrdist'+str(nbrdist_percentile)

for epitope in all_tcrs:
    for ab in ['A','B','AB']:
        avg_nbrdist = 0.0
        total = 0
        for mouse,tcrs in all_tcrs[epitope].iteritems():
            for ptcr,ntcr in tcrs:
                ## is this clone shared across other mice?
                info = all_info[ ptcr[6] ]
                nbrdist_score = float( info['{}_{}{}'.format(epitope,ab,nbrdist_suffix)] )
                avg_nbrdist += nbrdist_score
                total += 1
        avg_nbrdist /= total
        outlog.write( 'avg_nbrdist: {} {} {:.3f} {} {}\n'.format( epitope, ab, avg_nbrdist, total, nbrdist_suffix ))




for comparison_mode in range(3): ## normal, genes only, distance-based

    for same_mouse in [True,False]:
        for ii_nuc in [0,1]:
            Log('comparison_mode: {} same_mouse: {} ii_nuc: {}'.format(comparison_mode,same_mouse,ii_nuc))
            for epitope1 in all_tcrs:
                for epitope2 in all_tcrs:
                    if epitope2<epitope1: continue
                    if same_mouse and epitope1==epitope2: continue

                    all_div = {}

                    for chains in ['A','B','AB']:
                        if chains in fake_chains: continue

                        overlaps = []

                        for m1,tcrs1 in all_tcrs[epitope1].iteritems():
                            for m2,tcrs2 in all_tcrs[epitope2].iteritems():
                                if epitope1==epitope2 and m2<=m1:continue
                                if same_mouse != ( m1==m2 ): continue

                                overlap=0
                                for t1 in tcrs1:
                                    for t2 in tcrs2:
                                        if same_tcr(t1[ii_nuc],t2[ii_nuc],chains,comparison_mode):
                                            overlap += 1

                                overlaps.append( ( overlap, len(tcrs1)*len(tcrs2) ) )

                        overlaps.sort()
                        overlaps.reverse()

                        total_pairs = len(overlaps)
                        overlapping_pairs = len( [x for x in overlaps if x[0]>0 ] )

                        if overlapping_pairs:
                            overlap_count = sum( ( x[0] for x in overlaps ) )
                            total_count = sum( ( x[1] for x in overlaps ) )

                            p0, interval, fraction = confidence_interval( overlap_count, total_count, 0.95 )
                            diversity = 1.0/p0
                            diversity_upper = 1.0/interval[0]
                            diversity_lower = 1.0/interval[1]
                        else: ## no overlap at all
                            diversity = 0.
                            diversity_upper = 0
                            diversity_lower = 0

                        all_div[chains] = diversity

                        if chains == 'AB':
                            est = all_div.get('A',1.0)*all_div.get('B',1.0)
                            estimate = ' est_unpaired: {:9.1f} ratio: {:9.3f} '.format( est, 0 if diversity==0 else est/diversity)
                        else:
                            estimate = ' '

                        top3_string = ' '.join( ['{}'.format(x) for x,y in overlaps[:3]] )
                        #top3_string = ' '.join( ['{}/{}'.format(x,y) for x,y in overlaps[:3]] )
                        outlog.write('{} CM{:d} SM{:d} SE{:d} {:{}s} {:{}s} {:2s} div: {:9.1f} {:9.1f} {:9.1f} {}pairs: {:4d} o_pairs: {:4d} top3: {}\n'\
                                     .format( 'NA' if ii_nuc else 'AA',
                                              comparison_mode,
                                              same_mouse, epitope1 == epitope2,
                                              epitope1, max((len(x) for x in all_tcrs)),
                                              epitope2, max((len(x) for x in all_tcrs)),
                                              chains, diversity, diversity_lower, diversity_upper,
                                              estimate,
                                              total_pairs, overlapping_pairs, top3_string ) )
outlog.close()
