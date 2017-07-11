from basic import *
import html_colors
import tcr_distances
import util

with Parser(locals()) as p:
    #p.str('args').unspecified_default().multiple().required()
    p.str('clones_file').required() ##
    p.int('nbrdist_percentile').default(10)
    #p.float('float_arg')     # --float_arg 9.6
    p.flag('show')       # --flag_arg  (no argument passed)
    p.flag('skip_controls')       # --flag_arg  (no argument passed)
    p.flag('simple_correlations')       # --flag_arg  (no argument passed)
    #p.range('range_arg')     # --range_arg 1:2
    #p.multiword('multi_arg') # --multi_arg hello world
    #p.file('file_arg')       # --file_arg README.txt
    #p.directory('dir_arg')   # --dir_arg /tmp/
    #p.str('floatlist').cast(lambda x: [float(val) for val in x.split(',')])
    p.multiword('nbrdist_percentiles').cast(lambda x: [int(val) for val in x.split()]).default("10 25")

import matplotlib
if not show: matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

clones_file_with_nbrdists = '{}_nbrdists.tsv'.format(clones_file[:-4])
assert exists( clones_file_with_nbrdists )

fake_chains = util.detect_fake_chains( clones_file )

##
all_tcrs = {}
all_info = []
infields = []
Log('reading '+clones_file_with_nbrdists)
for line in open( clones_file_with_nbrdists,'r'):
    if not infields:
        if line[0] == '#':
            infields = line[1:-1].split('\t')
        else:
            infields = line[:-1].split('\t')
        continue
    assert infields

    l = parse_tsv_line( line[:-1], infields )

    epitope = l['epitope']

    if epitope not in all_tcrs: all_tcrs[epitope] = []
    all_tcrs[epitope].append( l )
    all_info.append( l )


## let's look at covariation between the different nbrdist metrics
nrows = 2
ncols = 3

if skip_controls:
    epitopes = [ x for x in all_tcrs if 'NEG_CNTRL' not in x and '_DP' not in x and 'EMPTY' not in x and 'NA_' not in x]
else:
    epitopes = all_tcrs.keys()[:]

epitopes.sort()
max_epitope_len = max( ( len(x) for x in epitopes ) )
suffixes = [ '_wtd_nbrdist'+str(nbrdist_percentile) ]
#suffixes = [ '_nbrdist'+str(nbrdist_percentile), '_wtd_nbrdist'+str(nbrdist_percentile) ]

if simple_correlations:
    suffixes_corr = ['_wtd_nbrdist'+str(nbrdist_percentile)]
    all_chains_corr = ['AB']
    figsize = (12,12)
else:
    suffixes_corr = suffixes
    all_chains_corr = ['A','B','AB']
    figsize = (12,8)

nrows = len(suffixes_corr)
ncols = len(all_chains_corr)

from scipy.stats import pearsonr
plt.figure(4,figsize=figsize)
plotno=0


Log('computing nbrdist correlations')
for suffix in suffixes_corr:
    for chains in all_chains_corr:
        plotno += 1
        plt.subplot(nrows,ncols,plotno)

        all_vals = []
        for ep in epitopes:
            tag = '{}_{}{}'.format( ep, chains, suffix )
            all_vals.append( [ float( x[tag] ) for x in all_info ] )

        if chains in fake_chains:
            continue

        A = np.zeros( ( len(epitopes),len(epitopes) ) )
        D = np.zeros( ( len(epitopes),len(epitopes) ) )
        for i_ep1,ep1 in enumerate(epitopes):
            A[i_ep1][i_ep1] = 1.0
            vals1 = all_vals[i_ep1]
            for i_ep2,ep2 in enumerate(epitopes):
                if i_ep2 <= i_ep1: continue
                vals2 = all_vals[i_ep2]

                r,p = pearsonr( vals1, vals2 )
                A[i_ep1][i_ep2] = r
                A[i_ep2][i_ep1] = r
                D[i_ep1][i_ep2] = 1.0-r
                D[i_ep2][i_ep1] = 1.0-r
                if simple_correlations:
                    print 'RANK_CORR {:7.3f} {} {} {}'.format( r, chains, ep1, ep2 )

        epitopes, leaves = util.tree_sort( epitopes, D )

        A2 = np.zeros( ( len(epitopes),len(epitopes) ) )
        for i in range(len(epitopes)):
            for j in range(len(epitopes)):
                A2[i,j] = A[ leaves[i], leaves[j] ]


        plt.imshow( A2, interpolation='nearest', #aspect = aspect,
                    vmin=0.4, vmax=1.0 )

        if chains == 'AB':
            A_AB = A2

        if suffix==suffixes_corr[-1]:
            plt.xticks( range(len(epitopes)), epitopes, rotation='vertical' )
        else:
            plt.xticks( [], [] )

        if False and chains==all_chains_corr[0]: ## no yticks since the order is changing
            plt.yticks( range(len(epitopes)), epitopes )
        else:
            plt.yticks( [],[] )

        if suffix==suffixes_corr[0]: plt.title(chains)
        if chains==all_chains_corr[0]:
            assert suffix == '_wtd_nbrdist'+str(nbrdist_percentile) ## hack
            ylabel = 'wtd_nbrdist{}'.format( nbrdist_percentile )
            plt.ylabel(ylabel)


plt.suptitle('epitope-epitope nbrdist-score correlations (colorscale=0.4 to 1.0)')
plt.subplots_adjust( left=0.15, bottom=0.15, top=0.925, right = 0.925, hspace=0.02, wspace=0.02 )
pngfile = '{}_epitope_correlations_{}.png'.format(clones_file[:-4],nbrdist_percentile )
print 'making:',pngfile
plt.savefig(pngfile)
util.readme(pngfile,"""These heat maps show correlations between nbrdist scores for different repertoire datasets.
The idea is that, given one reference epitope-specific repertoire, we can assign a number to any TCR in any of the
different epitope-specific repertoires which is the average distance
to the nearest TCRs in the reference repertoire. Two different repertoires give us two different sets of numbers that we
can compare
over the entire merged set of TCRs and ask how well they correlate. Two very similar repertoires should give similar
nbrdist scores for most TCRs in the big dataset, whereas very different repertoires will have uncorrelated nbrdist
scores. What's plotted is just the linear Pearson correlation coefficient between the nbrdist scores assigned by
each pair of reference repertoires, colored from blue (values less than 0.4) to dark red (1.0 -- perfect correlation).
""")



A = np.zeros( ( len(epitopes),len(epitopes) ) )
D = np.zeros( ( len(epitopes),len(epitopes) ) )
Log('analyzing epitope-epitope mutual nbrdist scores')
for i_ep1,ep1 in enumerate(epitopes):
    for i_ep2,ep2 in enumerate(epitopes):
        if i_ep2<i_ep1: continue
        for chains in ['A','B','AB']:
            dists = []
            ns = []
            for r in range(2):
                e1 = ep1 if r==0 else ep2
                e2 = ep2 if r==0 else ep1

                ## use rank score for this
                tag1 = '{}_{}_wtd_nbrdist{}rank'.format( e1, chains, nbrdist_percentile )
                tag2 = '{}_{}_wtd_nbrdist{}rank'.format( e2, chains, nbrdist_percentile )

                ## rank score wrt the other epitope (e2) minus rank score for self epitope (e1)
                vals = [ float( x[tag2] ) for x in all_tcrs[e1] ]
                #vals = [ float( x[tag2] ) - float( x[tag1] ) for x in all_tcrs[e1] ]

                #vals.sort()
                #n = len(vals)/4
                #sumval = sum( vals[ :n ] )/float(n)
                ## straight-up average
                sumval = sum( vals )/float(len(vals))
                dists.append( sumval )
                ns.append( len(vals) )
                print 'ORDERED {:7.1f} {:3d} {:{}s} {:{}s} {} {}'\
                    .format( dists[-1], ns[-1], e1, max_epitope_len, e2, max_epitope_len,
                             chains, suffix )
            avgval = sum(dists)/2
            print 'AVG {:7.1f} {:3d} {:{}s} {:{}s} {} {}'\
                .format( avgval, min(ns), ep1, max_epitope_len, ep2, max_epitope_len,
                         chains, suffix )
            A[i_ep1][i_ep2] = 100 - avgval ## so red is small dists
            A[i_ep2][i_ep1] = 100 - avgval ## avgval will be mostly in range of 50 (same rep) to 100 (very different)

            dist = max( 0, avgval-50.0 ) if i_ep1 != i_ep2 else 0.0
            D[i_ep1][i_ep2] = dist
            D[i_ep2][i_ep1] = dist

## reorder the epitopes
epitopes,leaves = util.tree_sort( epitopes, D )

A2 = np.zeros( ( len(epitopes),len(epitopes) ) )
for i in range(len(epitopes)):
    for j in range(len(epitopes)):
        A2[i,j] = A[ leaves[i], leaves[j] ]

## let's imshow the epitope distances
plt.figure(5,figsize=(10,10))

plt.imshow( A2, interpolation='nearest', vmin=0, vmax=50 ) ## since we took 100-val

plt.xticks( range(len(epitopes)), epitopes, rotation='vertical' )
plt.yticks( range(len(epitopes)), epitopes )

plt.title('epitope-epitope average mutual NNdist rank scores (50=red, 100=blue)')
plt.subplots_adjust( left=0.2, bottom=0.2, top=0.95, right = 0.95 )

pngfile = '{}_epitope_epitope_avg_nbrdist_rank_scores.png'.format(clones_file[:-4] )
print 'making:',pngfile
plt.savefig(pngfile)
util.readme(pngfile,"""This heat map shows a distance measure between repertoires based on the NNdistance rank score: for any pair of repertoires, we can compute
the NNdistance rank scores of the TCRs in the first repertoire with respect to the second as a reference, and also the NNdistance rank scores of the TCRs in the second repertoire
now using the first as the reference. We can average these numbers to get an average 'mutual' NNdistance rank score for a pair of repertoires. Lower means closer, and
is colored red here.
""")
#exit()


# plt.figure(7,figsize=(10,10))
# plt.imshow( A_AB, interpolation='nearest', vmin=0.0, vmax=1.0 )
# plt.xticks( range(len(epitopes)), epitopes, rotation='vertical' )
# plt.yticks( range(len(epitopes)), epitopes )
# plt.subplots_adjust( left=0.2, bottom=0.2, top=0.95, right = 0.95 )
# plt.suptitle('epitope-epitope raw-rank-score correlations (colorscale=0.0 to 1.0)')
# pngfile = '{}_epitope_correlations_AB_{}.png'.format(clones_file[:-4],nbrdist_percentile )
# print 'making:',pngfile
# plt.savefig(pngfile)
# util.readme(pngfile,"""Epitope correlations, AB chains, colored 0.0 to 1.0
# """)


## read distances
all_chain_dists = {}
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
        all_dists.append( dists )
    all_chain_dists[chains] = all_dists




## lets first look at epitope similarities

A = np.zeros( ( len(epitopes),len(epitopes) ) )
B = np.zeros( ( len(epitopes),len(epitopes) ) )
Log('analyzing epitope-epitope tcr distances')
for i_ep1,ep1 in enumerate(epitopes):
    for i_ep2,ep2 in enumerate(epitopes):
        if i_ep2<i_ep1: continue
        for chains in ['AB']:
            ## first compute based on all the distances
            ee_dists = []
            dist_threshold = 35*len(chains) ## not scaled down by 0.01
            small_dist_count = 0
            for i1,l1 in enumerate(all_info):
                if l1['epitope'] != ep1: continue
                dists = all_chain_dists[chains][i1]
                for i2,l2 in enumerate(all_info):
                    if l2['epitope'] != ep2: continue
                    if i1 != i2:
                        ee_dists.append(dists[i2])
                        if dists[i2] <= dist_threshold:
                            small_dist_count += ( dist_threshold - dists[i2] ) / dist_threshold
            ee_dists.sort()


            wtd_nbrdist10 = tcr_distances.sort_and_compute_weighted_nbrdist_from_distances( ee_dists, 10, dont_sort=True )

            if chains == 'AB': ## for imshow plot of distances
                A[i_ep1][i_ep2] = -1 * wtd_nbrdist10
                A[i_ep2][i_ep1] = -1 * wtd_nbrdist10
                B[i_ep1][i_ep2] = ( float(small_dist_count)/len(ee_dists) )**(1.0/4)
                B[i_ep2][i_ep1] = ( float(small_dist_count)/len(ee_dists) )**(1.0/4)

            # if False:
            #     nbrdists = []
            #     wtd_nbrdists = []
            #     for percentile in [5,10,25]:
            #         nbrdists.append( tcr_distances.sort_and_compute_nbrdist_from_distances( ee_dists, percentile,
            #                                                                         dont_sort=True ) )
            #         wtd_nbrdists.append( tcr_distances.sort_and_compute_weighted_nbrdist_from_distances( ee_dists, percentile,
            #                                                                                      dont_sort=True ) )
            #     if ep1<=ep2:
            #         print 'dists {:2s} {:{}s} {:{}s} {} {}'\
            #             .format( chains, ep1, max_epitope_len, ep2, max_epitope_len,
            #                      ' '.join( ( '%7.2f'%x for x in nbrdists ) ),
            #                      ' '.join( ( '%7.2f'%x for x in wtd_nbrdists ) ) )


            ## try computing avg nbrdist for ep1 vs ep2 , vice versa, and subtract from this self-avg-nbrdists...

            #continue
## let's imshow the epitope distances
plt.figure(2,figsize=(10,10))

#aspect = 'auto'

plt.imshow( A, interpolation='nearest', #aspect = aspect,
            vmin=-250, vmax=0 )

plt.xticks( range(len(epitopes)), epitopes, rotation='vertical' )
plt.yticks( range(len(epitopes)), epitopes )
#plt.yticks( range(len(entropy_keys_single)), entropy_keys_single )

#plt.title('gene entropies (colorscale: {:.1f}-{:.1f})'\
#          .format(min_entropy_for_colorscale,max_entropy_for_colorscale))

plt.title('epitope-epitope repertoire distance distribution summary statistic (red is smaller distances)')
plt.subplots_adjust( left=0.2, bottom=0.2, top=0.95, right = 0.95 )

pngfile = '{}_epitope_distances.png'.format(clones_file[:-4] )
print 'making:',pngfile
plt.savefig(pngfile)
util.readme(pngfile,"""Another way of comparing two repertoires/datasets is just to look at all the distances between a TCR in
one repertoire and a TCR in another, If they are similar repertoires, these numbers will be smaller than if the
repertoires are very different. What's plotted here is the average of the top 10 percent of inter-repertoire TCR distances,
with a color scale running from red (low distances, similar repertoires) to blue (high distances, divergent repertoires).
""")



## let's imshow the epitope distances
plt.figure(3,figsize=(10,10))

#aspect = 'auto'

plt.imshow( B, interpolation='nearest' ), #aspect = aspect,
#vmin=0, vmax=0 )

plt.xticks( range(len(epitopes)), epitopes, rotation='vertical' )
plt.yticks( range(len(epitopes)), epitopes )
#plt.yticks( range(len(entropy_keys_single)), entropy_keys_single )

#plt.title('gene entropies (colorscale: {:.1f}-{:.1f})'\
#          .format(min_entropy_for_colorscale,max_entropy_for_colorscale))

plt.title('epitope-epitope repertoire distance fraction below threshold')
plt.subplots_adjust( left=0.2, bottom=0.2, top=0.95, right = 0.95 )

pngfile = '{}_epitope_small_distances_fraction.png'.format(clones_file[:-4] )
print 'making:',pngfile
plt.savefig(pngfile)
# not being used right now
#util.readme(pngfile,"""
#""")


plt.figure(1,figsize=(12,8))

ncols = 3
nrows = len(suffixes)

epitope_colors = dict( zip( epitopes, html_colors.get_rank_colors_no_lights( len(epitopes)) ) )

from scipy.stats import gaussian_kde

all_min_max = {'A':[-50,250],'B':[-50,250],'AB':[-25,500]}

plotno = 0
for ii_suffix,suffix in enumerate(suffixes):
    for ii_chains, chains in enumerate( ['A','B','AB'] ):
        plotno += 1
        plt.subplot( nrows, ncols, plotno )

        if chains in fake_chains: continue

        mn,mx = all_min_max[chains ]

        sortl = []
        for epitope in epitopes:
            tcrs = all_tcrs[epitope]
            tag = '{}_{}{}'.format( epitope, chains, suffix )
            vals = [ float( x[tag] ) for x in tcrs ]
            avgval = sum(vals)/len(vals)
            sortl.append( ( avgval, epitope ) )
        sortl.sort()

        for (avgval,epitope) in sortl:
            tcrs = all_tcrs[epitope]
            tag = '{}_{}{}'.format( epitope, chains, suffix )
            vals = [ float( x[tag] ) for x in tcrs ]
            avgval = sum(vals)/len(vals)

            #mn = min(vals)
            #mx = max(vals)

            density = gaussian_kde( vals )
            xs = np.linspace( mn, mx, 100 )
            ys = density(xs)

            plt.plot( xs, ys, c=epitope_colors[epitope],label='{:5.1f} {}'.format(avgval,epitope) )
            #ymax = 1.1 * max(ys0)
            #plt.title('{} {} AUC= {}'.format(epitope, scoretag_for_output, int(100*area)))
            #plt.ylim( (0,ymax))

        plt.xlim( (mn,mx))
        plt.subplot(nrows,ncols,plotno)
        plt.legend(fontsize=8,frameon=False)
        plt.xlabel( '{}_{}'.format( chains, suffixes[ (plotno-1)/ncols ] ) )

plt.subplots_adjust(left=0.06,right=0.97,bottom=0.08,top=0.95)
plt.suptitle('nbr-distance distributions, lower is more clustered; legend shows average values')

pngfile = '{}_nbrdist_distributions.png'.format( clones_file[:-4] )
print 'making:',pngfile
plt.savefig(pngfile)
# not being used right now
#util.readme(pngfile,"""
#""")

###

for nbrdist_perc in nbrdist_percentiles:

    plt.clf()

    ncols = 3
    nrows = 2

    suffix = '_wtd_nbrdist'+str(nbrdist_perc)


    epitope_colors = dict( zip( epitopes, html_colors.get_rank_colors_no_lights( len(epitopes)) ) )

    from scipy.stats import gaussian_kde

    all_min_max = {'A':[-50,250],'B':[-50,250],'AB':[-25,500]}

    plotno = 0
    for ii_cdf in range(2):
        for ii_chains, chains in enumerate( ['A','B','AB'] ):
            plotno += 1
            plt.subplot( nrows, ncols, plotno )

            if chains in fake_chains: continue

            mn,mx = all_min_max[chains ]

            sortl = []
            realmn,realmx=500,0
            for epitope in epitopes:
                tcrs = all_tcrs[epitope]
                tag = '{}_{}{}'.format( epitope, chains, suffix )
                vals = [ float( x[tag] ) for x in tcrs ]
                avgval = sum(vals)/len(vals)
                sortl.append( ( avgval, epitope ) )
                realmn = min( realmn, min(vals))
                realmx = max( realmx, max(vals))
            sortl.sort()

            if ii_cdf==0:
                for (avgval,epitope) in sortl:
                    tcrs = all_tcrs[epitope]
                    tag = '{}_{}{}'.format( epitope, chains, suffix )
                    vals = [ float( x[tag] ) for x in tcrs ]
                    avgval = sum(vals)/len(vals)

                    #mn = min(vals)
                    #mx = max(vals)

                    density = gaussian_kde( vals )
                    xs = np.linspace( mn, mx, 100 )
                    ys = density(xs)

                    plt.plot( xs, ys, c=epitope_colors[epitope],label='{:5.1f} {}'.format(avgval,epitope) )
                    #ymax = 1.1 * max(ys0)
                    #plt.title('{} {} AUC= {}'.format(epitope, scoretag_for_output, int(100*area)))
                    #plt.ylim( (0,ymax))

                plt.xlim( (mn,mx))
                if ii_chains==0: plt.ylabel( 'PDF' )
                plt.legend(fontsize=8,frameon=False)
            else:
                assert ii_cdf==1
                for (avgval,epitope) in sortl:
                    tcrs = all_tcrs[epitope]
                    tag = '{}_{}{}'.format( epitope, chains, suffix )
                    xvals = [ float( x[tag] ) for x in tcrs ]
                    xvals.sort()

                    N = len(xvals)
                    yvals = [ float(x+1)/N for x in range(N) ]

                    yvals = [0.0]+yvals
                    xvals = [0.0]+xvals

                    avgval = sum(xvals)/N

                    plt.plot( xvals, yvals, c=epitope_colors[epitope],label='{:5.1f} {}'.format(avgval,epitope) )
                    #ymax = 1.1 * max(ys0)
                    #plt.title('{} {} AUC= {}'.format(epitope, scoretag_for_output, int(100*area)))
                    #plt.ylim( (0,ymax))

                plt.xlim( (-5,realmx+5))
                if ii_chains==0: plt.ylabel( 'CDF' )
                plt.legend(fontsize=8,frameon=False,loc='best')
            plt.xlabel( '{}_{}'.format( chains, suffix ) )

    plt.subplots_adjust(left=0.06,right=0.97,bottom=0.08,top=0.95)
    plt.suptitle('nbr-distance distributions, lower is more clustered; legend shows average values')

    pngfile = '{}_nbrdist{}_distributions_w_cdf.png'.format( clones_file[:-4], nbrdist_perc )
    print 'making:',pngfile
    plt.savefig(pngfile)
    util.readme(pngfile,"""These plots show distributions of the nbr-distance (aka NNdistance) score for the different
epitopes in the dataset. The score assigned to a given TCR is just the weighted average of the distance
to the nearest {}th percentile
of the other TCRs in the repertoire. Lower values indicate more clustering nearby the TCR in question. The three
columns are based on (left to right) alpha-chain-only, beta-chain-only, and paired-chains distances.
The top row shows the smoothed frequency distributions of nbr-distance score, and the bottom row
shows the cumulative distributions: the height of the curve above a given nbrdistance value is the fraction
of the repertoire with scores less than or equal to that value (so a more clustered repertoire will have a curve that rises
higher earlier than a less clustered repertoire). The legends give the average value of
nbr-distance over the entire repertoire. Lower values correspond to more clustered / less diverse repertoires. These
averages are also the areas above the CDF curves, and are the diversity metrics plotted in the middle panels of the diversity
summary plots.
    """.format(nbrdist_perc))

if show:
    plt.show()




