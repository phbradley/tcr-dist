from basic import *
import html_colors
import util

with Parser(locals()) as p:
    p.str('clones_file').required()
    p.str('organism').required()
    p.str('outfile_prefix')
    p.flag('show')       # --flag_arg  (no argument passed)
    p.flag('paper_figs')       # --flag_arg  (no argument passed)
    p.flag('include_gene_matching')       # --flag_arg  (no argument passed)
    p.multiword('skip_epitopes').cast(lambda x:x.split())
    p.multiword('epitopes').cast(lambda x:x.split())


logfile = clones_file[:-4]+'_sharing.log'
auc_logfile= clones_file[:-4]+'_random_aucs.log'

fake_chains = util.detect_fake_chains( clones_file )

if not exists( logfile ):
    print 'Sorry, you need to run analyze_overlap_compute_simpsons.py before this script'
    print 'If successful, it will generate',logfile
    exit()

if not exists( auc_logfile ):
    print 'Sorry, you need to run random_tcr_distances.py and read_random_tcr_distances.py before this script'
    print 'If successful, they will generate',auc_logfile
    exit()

cmdline_epitopes = epitopes

if not outfile_prefix:
    outfile_prefix = logfile[:-4]

import matplotlib
if not show:
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
#import numpy as np

assert exists(logfile)

epitope_diversity = {} ## store overall diversity measures

## order for the plots?
##
## first diversity: alpha, beta, alpha-beta
##
## then cross-reactivity:
##
##


plt.figure(1,figsize=(12,30))

nrows = 5
ncols = 1



if include_gene_matching:
    colorkey = "(red=full-chain matching, green=gene-segment matching, blue=distance-threshold matching)"
else:
    colorkey = "(red=standard, blue=distance-matching, green=TCRdiv)"

if True:
    ## let's plot simpsons diversity

    all_dats = {}
    all_mice_dats = {}
    for line in open(logfile,'r'):
        l = line.split()
        if l[0] != 'clone_diversity:':continue
        epitope = l[1]

        if skip_epitopes and epitope in skip_epitopes: continue
        if cmdline_epitopes and epitope not in cmdline_epitopes: continue

        div,divlo,divhi,p = [float(x) for x in l[2:6] ]
        if divlo == 0 and divhi==0:
            assert p==0

        all_dats[epitope] = (div,divlo,divhi)

        info_strings = [x.split(',') for x in l[-1].split(';')]
        all_mice_dats[epitope] = [ (float(x[0]), int(x[1])) for x in info_strings ]


    for r in range(1,3):


        ax = plt.subplot(nrows,ncols,r) ## amino acid diversity
        if r== 2:
            ax.set_yscale( "log" )

        ax.yaxis.grid(True,which='major')

        counter=0

        xvals = []
        yvals = []
        locs_labels = []

        epitopes = all_dats.keys()

        l = [ ( all_dats[x][0],x ) for x in epitopes ]
        l.sort()

        epitopes = [x[1] for x in l ] ## now sorted

        #allyvals = []
        for epitope in epitopes:
            counter += 1


            if r==1:
                if all_dats[epitope][1] == 0:
                    div,divlo,divhi = 1.,1.,1.
                else:
                    div,divlo,divhi = [(1.0-1.0/x) for x in all_dats[epitope] ]
            else:
                if all_dats[epitope][1]==0: continue
                div,divlo,divhi = all_dats[epitope]

            #allyvals.extend( all_dats[epitope] )

            line_style = '-'
            point_style = 'o'
            color = 'r'

            ## a line from divlo to divhi
            plt.plot( [counter,counter], [divlo,divhi], line_style, color = color )

            plt.scatter( [counter], [div], marker=point_style, color = color )
            locs_labels.append( ( counter, epitope ) )

            if r==1:
                if epitope in all_mice_dats:
                    yvals = []
                    xvals = []
                    for (p,total_tcrs_this_mouse) in all_mice_dats[epitope]:
                        yvals.append(1.0-p)
                        xvals.append(counter)
                        plt.text( counter, 1.0-p, '  {}'.format(total_tcrs_this_mouse), fontsize='xx-small',va='center' )
                    plt.scatter( xvals, yvals, marker='_', c='r' )

                plt.ylim((-0.05,1.05))

        plt.xticks( [x[0] for x in locs_labels], [x[1] for x in locs_labels], rotation='vertical', fontsize=10 )
        plt.xlim((0,counter+1))
        if r==1:
            plt.title("1-simpson's diversity of clone-size distributions")
        else:
            plt.title("inverse simpson's diversity of clone-size distributions")

if True:
    ## let's plot simpsons diversity

    all_dats = {}
    for line in open(logfile,'r'):
        l = line.split()
        if l[0] != 'AA':continue
        assert l[7] == 'div:'
        aana = l[0]
        comparison_mode = int( l[1][2] )
        same_mouse = int( l[2][2] )
        same_epitope = int( l[3][2] )
        chains = l[6]
        #ichains = all_chains.index(chains )
        epitope,ep2 = l[4:6]

        if skip_epitopes and epitope in skip_epitopes: continue
        if cmdline_epitopes and epitope not in cmdline_epitopes: continue
        if comparison_mode==1 and not include_gene_matching: continue

        if same_mouse or not same_epitope: continue

        div,divlo,divhi = [float(x) for x in l[8:11] ]

        if div==0:continue

        if chains not in all_dats: all_dats[chains] = {}

        if epitope not in all_dats[chains]: all_dats[chains][epitope] = {}

        top3 = ','.join( line[ line.index('top3:')+5:-1].split())
        if not top3: top3='-'

        all_dats[chains][epitope][comparison_mode] = (div,divlo,divhi,top3)

    ## read the "green" data-- gaussian simpsons
    for line in open(logfile,'r'):
        l = line.split()
        if l[0] != 'GAUSSDIV': continue

        comparison_mode = 1 ## hack-- replace gene segment matching
        same_mouse = int( l[1][2] )
        same_epitope = int( l[2][2] )
        chains = l[5]
        #ichains = all_chains.index(chains )
        epitope,ep2 = l[3:5]
        div = float(l[7])

        if skip_epitopes and epitope in skip_epitopes: continue
        if cmdline_epitopes and epitope not in cmdline_epitopes: continue
        #if comparison_mode==1 and not include_gene_matching: continue

        if not same_mouse or not same_epitope: continue

        if chains not in all_dats: all_dats[chains] = {}
        if epitope not in all_dats[chains]: all_dats[chains][epitope] = {}
        all_dats[chains][epitope][comparison_mode] = (div,div,div,'')

        ## use this as a summary statistic for diversity
        if chains not in epitope_diversity: epitope_diversity[chains] = {}
        if epitope not in epitope_diversity[chains]: epitope_diversity[chains][epitope] = [0,0,0,0]
        epitope_diversity[chains][epitope][0] = div


    ax = plt.subplot(nrows,ncols,3) ## amino acid diversity
    ax.set_yscale( "log" )
    ax.yaxis.grid(True,which='major')


    ## what order to plot? increasing div at comparison_mode0 from left to right


    counter=0

    xvals = []
    yvals = []
    locs_labels = []

    for chains in ['A','B','AB']:
        counter += 1
        if chains not in all_dats: continue
        epitopes = all_dats[chains].keys()

        comparison_mode_for_sorting = 1 ## identity based on gene segments
        comparison_mode_for_sorting = 2 ## based on small distances


        l = []
        for epitope in epitopes:
            if comparison_mode_for_sorting in all_dats[chains][epitope]:
                sortval = all_dats[chains][epitope][comparison_mode_for_sorting][0]
            else:
                sortval = min( ( dats[0] for dats in all_dats[chains][epitope].values() ) )
            l.append( ( sortval, epitope ) )

        #l = [ ( all_dats[chains][x][comparison_mode_for_sorting][0],x ) for x in epitopes ]
        l.sort()

        epitopes = [x[1] for x in l ] ## now sorted


        for epitope in epitopes:
            counter += 1

            for comparison_mode in [0,2,1]:
                if comparison_mode not in all_dats[chains][epitope].keys(): continue

                div,divlo,divhi,top3 = all_dats[chains][epitope][comparison_mode]

                line_style = ['-','--',':'][ comparison_mode ]
                point_style = ['o','s','D'][ comparison_mode ]
                color = 'rgb'[ comparison_mode ]

                ## a line from divlo to divhi
                if divlo<divhi: plt.plot( [counter,counter], [divlo,divhi], line_style, color = color )

                plt.scatter( [counter], [div], marker=point_style, color = color )
            locs_labels.append( ( counter, '{} ({})'.format(epitope,chains) ) )


    plt.xticks( [x[0] for x in locs_labels], [x[1] for x in locs_labels], rotation='vertical', fontsize=10 )
    plt.xlim((1,counter+1))
    plt.title("inv-simpson's sharing-based repertoire diversity "+colorkey)

if True: ####################################################################### now plot sharing

    for desired_same_mouse in [0,1]:

        all_dats = {}
        for line in open(logfile,'r'):
            l = line.split()
            if l[0] != 'AA':continue
            assert l[7] == 'div:'
            aana = l[0]
            comparison_mode = int( l[1][2] )
            same_mouse = int( l[2][2] )
            same_epitope = int( l[3][2] )
            chains = l[6]
            #ichains = all_chains.index(chains )
            ep1,ep2 = l[4:6]

            if skip_epitopes and (ep1 in skip_epitopes or ep2 in skip_epitopes): continue
            if cmdline_epitopes and (ep1 not in cmdline_epitopes or ep2 not in cmdline_epitopes): continue
            if comparison_mode==1 and not include_gene_matching: continue

            if same_epitope or ( same_mouse != desired_same_mouse ): continue

            if chains not in all_dats: all_dats[chains] = {}

            epitope = (ep1,ep2)

            if float(l[8]) == 0:continue

            div,divlo,divhi = [1.0/float(x) for x in l[8:11] ]

            top3 = ','.join( line[ line.index('top3:')+5:-1].split())
            if not top3: top3='-'

            if epitope not in all_dats[chains]: all_dats[chains][epitope] = {}
            all_dats[chains][epitope][comparison_mode] = (div,divlo,divhi,top3)



        ## read the "green" data-- gaussian simpsons
        for line in open(logfile,'r'):
            l = line.split()
            if l[0] != 'GAUSSDIV': continue

            comparison_mode = 1 ## hack-- replace gene segment matching
            same_mouse = int( l[1][2] )
            same_epitope = int( l[2][2] )
            chains = l[5]
            #ichains = all_chains.index(chains )
            ep1,ep2 = l[3:5]
            div = float(l[7])
            if div==0:continue
            div = 1.0/div

            epitope=(ep1,ep2)

            if skip_epitopes and epitope in skip_epitopes: continue
            if cmdline_epitopes and epitope not in cmdline_epitopes: continue
            #if comparison_mode==1 and not include_gene_matching: continue

            if same_epitope or same_mouse!=desired_same_mouse: continue

            if chains not in all_dats: all_dats[chains] = {}
            if epitope not in all_dats[chains]: all_dats[chains][epitope] = {}
            all_dats[chains][epitope][comparison_mode] = (div,div,div,'')







        ax = plt.subplot(nrows,ncols,4+desired_same_mouse) ## amino acid diversity
        ax.set_yscale( "log" )
        ax.yaxis.grid(True,which='major')


        ## what order to plot? increasing div at comparison_mode0 from left to right


        counter=0

        xvals = []
        yvals = []
        locs_labels = []

        for chains in ['A','B','AB']:
            counter += 1
            if chains not in all_dats: continue
            epitopes = [ x for x,y in all_dats[chains].iteritems() if y.keys() != [1] ]

            #comparison_mode_for_sorting = 1 ## identity based on gene segments
            comparison_mode_for_sorting = 2 ## based on small distances

            l = []
            for epitope in epitopes:
                if comparison_mode_for_sorting in all_dats[chains][epitope]:
                    sortval = all_dats[chains][epitope][comparison_mode_for_sorting][0]
                else:
                    sortval = min( ( dats[0] for dats in all_dats[chains][epitope].values() ) )
                l.append( ( sortval, epitope ) )
            l.sort()
            l.reverse()

            epitopes = [x[1] for x in l ] ## now sorted


            for epitope in epitopes:
                counter += 1

                for comparison_mode in [0,2,1]:
                    if comparison_mode not in all_dats[chains][epitope].keys(): continue

                    div,divlo,divhi,top3 = all_dats[chains][epitope][comparison_mode]

                    line_style = ['-','--',':'][ comparison_mode ]
                    point_style = ['o','s','D'][ comparison_mode ]
                    color = 'rgb'[ comparison_mode ]

                    ## a line from divlo to divhi
                    plt.plot( [counter,counter], [divlo,divhi], line_style, color = color )

                    plt.scatter( [counter], [div], marker=point_style, color = color )
                locs_labels.append( ( counter, '{}-{} ({})'.format(epitope[0],epitope[1],chains) ) )


        plt.xticks( [x[0] for x in locs_labels], [x[1] for x in locs_labels], rotation='vertical', fontsize=10 )
        plt.xlim((1,counter+1))
        if desired_same_mouse == 0:
            plt.title("cross-reactive probability, different mice (same colors as above)")
        else:
            plt.title("cross-reactive probability, same mice (same colors as above)")


plt.subplots_adjust(bottom=0.1,top=0.97,left=0.05,right=0.97,hspace=0.3)

pngfile = '{}.png'.format(outfile_prefix)
print 'making:',pngfile
plt.savefig(pngfile)
util.readme(pngfile,"""These next five plots reflect different notions of sharing or repetition in the repertoire (ie, seeing "the same" TCR more than once).
The top two plots corresponds to repetition within a single mouse (ie, clonality), and give two representations of Simpson's measure: 1-Simpson's in the top
plot and 1/Simpson's in the bottom. The top plot shows 1-Simpson's for each of the individual mice (labeled by #reads in the mouse) together with an averaged value over
the entire repertoire (red disk, mice are weighted based on number of reads).
<br>
The middle plot uses the Simpson's diversity framework to analyze sharing across mice for the same epitope. The three colors correspond to three notions of
sharing: red points-- seeing the exact same TCR (alpha or beta or both chains); blue points-- seeing two TCRs within a (small) distance of one another; green points-- a
Gaussian-smoothed version of the blue points ("TCRdiv", see explanatory text for the very first figure on this page).
<br>
The bottom two plots look at cross-reactivity: sharing of TCRs between epitopes, either between different mice (plot#4) or within the same mouse (plot#5). The colors
are the same as in the middle plot.
""")

# util.make_readme( pngfile, """

# These four plots

# """ )

## read the other diversity summary statistic
for line in open(logfile,'r'):
    if line.startswith('GAUSSDIVSHANNON'):
        l = line.split()
        epitope,chains = l[1:3]
        shannon_diversity = float( l[4] )
        if skip_epitopes and epitope in skip_epitopes: continue
        if cmdline_epitopes and epitope not in cmdline_epitopes: continue
        epitope_diversity[chains][epitope][1] = shannon_diversity
    elif line.startswith('avg_nbrdist:'):
        l = line.split()
        epitope,chains = l[1:3]
        avg_nbrdist = float( l[3] )
        if skip_epitopes and epitope in skip_epitopes: continue
        if cmdline_epitopes and epitope not in cmdline_epitopes: continue
        epitope_diversity[chains][epitope][2] = avg_nbrdist

## load the random values for avg_nbrdist
rand_divs = {}
for chains in ['A','B','AB']: rand_divs[chains] = [[],[],[],[]]
randfile = '{}/{}_rand_divs_new.txt'.format(path_to_current_db_files(),organism)
if not exists(randfile):
    Log('WARNING:: plot_sharing.py: missing random divergences file: {}'.format(randfile))
else:
    for line in open(randfile,'r'):
        l = line.split()
        if line.startswith("GAUSSDIV SM1 SE1"):
            chains = l[5]
            div = float(l[7])
            rand_divs[chains][0].append( div )
        elif line.startswith("GAUSSDIVSHANNON"):
            chains = l[2]
            div = float( l[4] )
            rand_divs[chains][1].append( div )
        elif line.startswith('avg_nbrdist:'):
            fake_epitope,chains = l[1:3]
            avg_nbrdist = float( l[3] )
            rand_divs[chains][2].append( avg_nbrdist )

## let's try to get epitope diversity information from the aucs for discrimination from random tcrs
desired_nbrdist_tag_suffix = 'wtd_nbrdist10'
desired_nbrdist_label = 'nbrdist10p'
for line in open( auc_logfile,'r'):
    if line.startswith('auc_random '):
        l = line.split()
        auc = float( l[1] )
        chains = l[4]
        epitope = l[5]
        nbrdist_tag_suffix = l[6]
        if skip_epitopes and epitope in skip_epitopes: continue
        if cmdline_epitopes and epitope not in cmdline_epitopes: continue
        if nbrdist_tag_suffix == desired_nbrdist_tag_suffix:
            epitope_diversity[chains][epitope][3] = -1*auc ## sort from least to most diverse



plt.figure(2,figsize=(12,12))

## now make a diversity summary bar-plot figure
nrows = 3 ## A,B,AB
#ncols = 3
ncols = 3

epitopes = epitope_diversity['A'].keys()
#epitopes.sort()

plotno=0

for chains in ['A','B','AB']:
    if chains in fake_chains: continue
    #for ( divindex, divtype ) in [ ( 0, 'inv_simpson_gaussdist' ), (1, 'shannon_gaussdist'), (2, 'avg_nbrdist') ]:
    for ( divindex, divtype ) in [ ( 0, 'TCRdiv' ),
                                   ( 2, 'avg_nbrdist10p aka NN-distance'),
                                   ( 3, 'AUROC_random {}'.format(desired_nbrdist_label)) ]:

        l = [ ( epitope_diversity[chains][x][divindex],x) for x in epitopes ]
        l.sort()

        if paper_figs:
            for div,epitope in l:
                print 'epitope_diversity:',chains,epitope,'_'.join(divtype.split()),div

        plotno+=1
        ax = plt.subplot(nrows,ncols,plotno )

        if divindex<=1: ax.set_yscale( "log" )
        ax.yaxis.grid(True,which='major')

        if divindex==3: # special case
            aucs = [-100*x[0] for x in l] ## now go from 0 to 100
            heights = [x-100 for x in aucs]
            bottoms = [100]*len(l)
            lefts = range(len(l))
            plt.bar( lefts, heights, bottom=bottoms )
        else:
            plt.bar( range(len(l)), [x[0] for x in l ] )
        plt.xticks( [x+0.4 for x in range(len(l))], [x[1] for x in l ], rotation='vertical', fontsize=10 )
        plt.title('{} {}'.format(chains, divtype))

        rdivs = rand_divs[chains][divindex]
        if rdivs:
            rdivs.sort()
            for val in rdivs:
                plt.plot( [0,len(l)], [val,val], ':r' )

        locs,labels = plt.yticks()
        minval = min( ( x[0] for x in l ) )
        maxval = max( ( x[0] for x in l ) )
        ymn,ymx = plt.ylim()
        if divindex==2: ## only makes sense for linear plot
            locsep = locs[1]-locs[0]
            max_lowerloc = max( ( x for x in locs if x < minval-0.5*locsep ) )
            ymn = max( max_lowerloc, ymn )
        if divindex<=1: ## trim off some randoms?
            if max(locs) > maxval * 100:
                min_upperloc = min( ( x for x in locs if x > maxval*10 ) )
                ymx = min(ymx,min_upperloc)
            ##
        if divindex==3:
            ymn,ymx = ymx,ymn
        plt.ylim( ymn,ymx)

plt.suptitle( 'repertoire diversity summary' )

plt.subplots_adjust(bottom=0.1,top=0.93,left=0.05,right=0.97,hspace=0.35)

pngfile = '{}_diversity.png'.format(outfile_prefix)
print 'making:',pngfile
plt.savefig(pngfile)
util.readme(pngfile,"""These bar plots represent three ways of looking at the total diversity of an epitope-specific
repertoire, measured over single-chains (top two rows) and full receptors (bottom row). The left column shows the TCRdiv
diversity measure (a variant of the (inverse of) Simpson's Diversity Index that accounts for receptor similarity
as well as identity). The middle column shows the average NNdistance aka nbrdistance over all the TCR clones in the
repertoire and captures the average distance from TCR clones to their nearest neighbors in the repertoire.
The right column is based on comparison between the epitope-specific TCRs and random TCRs, specifically
looking at the "AUROC" measure. AUROC measures the area under the ROC curve, and captures the trade-off between true and
false positives as the value of a quality or sorting metric is varied from best to worst. An AUROC of 100 (perfect
separation) means that all the true positives occur before any of the false positives. An AUROC of 50 (random)
means that the positives and negatives are equally distributed with respect to the sorting measure. The sort used here
is the "nbrdist10P" which measures the weighted average distance to the closest tenth of the reference repertoire.
The AUROCs
are plotted with a flipped y-axis for easy comparison to the other columns,
the idea being that a more diverse repertoire will be harder to reliably separate
from the random background and hence will have smaller AUROC values (ie, lower AUROC == higher diversity).

""")



if show:
    plt.show()
