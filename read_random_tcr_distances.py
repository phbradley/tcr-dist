from basic import *
import numpy as np
import parse_tsv
from scipy.stats import gaussian_kde
import util
#from tcr_distances_blosum
import tcr_distances
import html_colors
from all_genes import all_genes

with Parser(locals()) as p:
    #p.str('organism').required()
    p.str('clones_file').required()
    p.str('file_tag').default('')
    p.str('outfile_prefix')
    p.str('organism').required()
    p.float('nbrdist_rescale').default(1.0)     # --float_arg 9.6
    p.str('single_nbrdist_tag_for_plotting').default('wtd_nbrdist10')
    p.flag('include_non_wtd')       # --flag_arg  (no argument passed)
    p.flag('show')       # --flag_arg  (no argument passed)
    p.flag('paper_supp')       # --flag_arg  (no argument passed)
    p.flag('exclude_xr')       # --flag_arg  (no argument passed)

import matplotlib
matplotlib.rcParams['mathtext.default'] = 'regular'
if not show: matplotlib.use('Agg')
import matplotlib.pyplot as plt

def greek_seg( seg ):
    return seg.replace( 'A', r'$\alpha$' ).replace( 'B', r'$\beta$' )

fake_chains = util.detect_fake_chains( clones_file )

## don't exclude xr guys just for discrimination vs random
##


ABs = ['A','B','AB']

## read random distances
random_nbrdists_file = '{}_random_nbrdists{}.tsv'.format( clones_file[:-4],file_tag )
assert exists(random_nbrdists_file)

outlogfile = '{}_random_aucs.log'.format( clones_file[:-4] )
print 'making:',outlogfile
outlog =open( outlogfile,'w')



if outfile_prefix is None:
    outfile_prefix = random_nbrdists_file[:-4]

infields = []
rand_nbrdists = []
for line in open( random_nbrdists_file,'r'):
    if not infields:
        infields = line[:-1].split('\t')
        rand_nbrdist_tags = [ x for x in infields if ( 'rank' not in x and 'nbrdist' in x ) ]
    else:
        l = parse_tsv_line( line[:-1], infields )
        nbrdists = dict( ( ( x, l[x]) for x in rand_nbrdist_tags ) )
        rand_nbrdists.append( nbrdists )

clones_file_with_nbrdists = '{}_nbrdists.tsv'.format( clones_file[:-4] )
header = open( clones_file_with_nbrdists,'r').readline()[:-1].split('\t')

nbrdist_tags = [ x for x in header if x in rand_nbrdist_tags and ( 'wtd' in x or include_non_wtd ) ]
nbrdist_tags.sort()
num_nbrdist_tags = len(nbrdist_tags)

Log('parsing {} for {} nbrdist_tags'.format( clones_file_with_nbrdists, num_nbrdist_tags ) )

tcr_fields = nbrdist_tags + [ 'va_genes','vb_genes','cdr3a','cdr3b']
all_tcrs = parse_tsv.parse_tsv_file( clones_file_with_nbrdists, ['epitope'], tcr_fields )

## look for cross-reactive tcrs
for e in all_tcrs:
    new_tcrs = []
    for l in all_tcrs[e]:
        va_reps = frozenset( ( all_genes[organism][x].rep for x in l[-4].split(';') ) )
        vb_reps = frozenset( ( all_genes[organism][x].rep for x in l[-3].split(';') ) )
        new_tcrs.append( [ (nbrdist_rescale * float(x)) for x in  l[:num_nbrdist_tags] ] +
                         [ va_reps, vb_reps, l[-2], l[-1], False ] ) ## add X-reactive flag
    all_tcrs[e] = new_tcrs

va_genes_index = num_nbrdist_tags
vb_genes_index = num_nbrdist_tags+1
cdr3a_index = num_nbrdist_tags+2
cdr3b_index = num_nbrdist_tags+3

def same_tcr( a, b ):
    return ( a[cdr3a_index] == b[cdr3a_index] and
             a[cdr3b_index] == b[cdr3b_index] and
             (not a[va_genes_index].isdisjoint( b[va_genes_index])) and
             (not a[vb_genes_index].isdisjoint( b[vb_genes_index]) ) )

if exclude_xr:
    print 'finding X-reactive'

    for e1 in all_tcrs:
        for e2 in all_tcrs:
            if e1==e2: continue
            for a in all_tcrs[e1]:
                for b in all_tcrs[e2]:
                    if same_tcr(a,b):
                        print 'X-react:',e1,a[cdr3a_index],a[cdr3b_index],e2,b[cdr3a_index],b[cdr3b_index]
                        a[-1] = True
                        b[-1] = True



epitopes = all_tcrs.keys()
epitopes.sort()

all_aucs_random = {}
all_aucs_others = {}

nbrdist_tags_for_plotting = []
for epitope in epitopes:
    for chains in ABs:
        prefix = '{}_{}_'.format(epitope,chains)
        for ii_nbrdist_tag, nbrdist_tag in enumerate( nbrdist_tags ):
            if not nbrdist_tag.startswith(prefix): continue
            nbrdist_tag_suffix = nbrdist_tag[len(prefix):]
            nbrdist_tag_suffix_w_chains = chains+'_'+nbrdist_tag[len(prefix):]
            if nbrdist_tag_suffix not in all_aucs_others:
                all_aucs_others[nbrdist_tag_suffix] = []
                all_aucs_random[nbrdist_tag_suffix] = []
                nbrdist_tags_for_plotting.append( nbrdist_tag_suffix )
            if nbrdist_tag_suffix_w_chains not in all_aucs_others:
                all_aucs_others[nbrdist_tag_suffix_w_chains] = []
                all_aucs_random[nbrdist_tag_suffix_w_chains] = []


## sort the nbrdist_tags_for_plotting
l = []
for tag in nbrdist_tags_for_plotting:
    if 'nbrdist' in tag:
        num = int( tag[ tag.index('nbrdist')+7:] )
        if num<0:
            l.append( ( -1*num, tag, 'nbrdist{}'.format(-1*num)) )
        else:
            l.append( ( 100*num,tag, 'nbrdist{}P'.format(num)))
    else:
        assert 'dens' in tag
        sdev = int( tag[ tag.index('dens')+4:] )
        l.append( ( 10000*sdev, tag, 'nbrdens{}'.format(sdev)))

l.sort()

nbrdist_tags_for_plotting = [x[1] for x in l] ## now sorted from most focused to most averaged
nbrdist_labels_for_plotting = [x[2] for x in l] ## now sorted from most focused to most averaged
#print 'nbrdist_tags_for_plotting:',nbrdist_tags_for_plotting
#print 'nbrdist_labels_for_plotting:',nbrdist_tags_for_plotting


nrows = len(epitopes)
ncols = len(nbrdist_tags_for_plotting)

nrows_single = len(epitopes)
ncols_single = 6

top_margin_inches = 1.25
bottom_margin_inches = 0.5
left_margin_inches = 1.0
right_margin_inches = 1.0
plot_height_inches = 2.0 * nrows
plot_width_inches  = 3.0 * ncols

fig_height = top_margin_inches + plot_height_inches + bottom_margin_inches
fig_width = left_margin_inches + plot_width_inches + right_margin_inches
fig_width_summary = 12.0

top_margin = float( plot_height_inches + bottom_margin_inches ) / fig_height
bottom_margin = float( bottom_margin_inches ) / fig_height
left_margin = float( left_margin_inches ) / fig_width
right_margin = float( left_margin_inches + plot_width_inches ) / fig_width

assert single_nbrdist_tag_for_plotting in nbrdist_tags_for_plotting

## these hold dats for the single_nbrdist_tag_for_plotting nbrdist_tag
save_epitope_nbrdists = {}
save_epitope_rocs = {}

for ii_epitope, epitope in enumerate(epitopes):
    other_epitopes = [ x for x in epitopes if x != epitope ]
    for ii_chains,chains in enumerate( ABs ):
        prefix = '{}_{}_'.format(epitope,chains)
        for ii_nbrdist_tag, nbrdist_tag in enumerate( nbrdist_tags ):
            if not nbrdist_tag.startswith(prefix): continue
            nbrdist_tag_suffix = nbrdist_tag[len(prefix):]
            nbrdist_tag_suffix_w_chains = chains+'_'+nbrdist_tag[len(prefix):]
            if nbrdist_tag_suffix not in all_aucs_others:
                all_aucs_others[nbrdist_tag_suffix] = []
                all_aucs_random[nbrdist_tag_suffix] = []
            if nbrdist_tag_suffix_w_chains not in all_aucs_others:
                all_aucs_others[nbrdist_tag_suffix_w_chains] = []
                all_aucs_random[nbrdist_tag_suffix_w_chains] = []

            positive_nbrdists = [ x[ ii_nbrdist_tag ] for x in all_tcrs[ epitope ] if not x[-1] ]
            negative_nbrdists_random = [ nbrdist_rescale*float( x[ nbrdist_tag ] ) for x in rand_nbrdists ]

            sign_factor = -1.0 if 'dens' in nbrdist_tag else 1.0

            auc_random, xvals_random, yvals_random = tcr_distances.compute_auc( positive_nbrdists, negative_nbrdists_random,
                                                                                sign_factor = sign_factor )

            outlog.write( 'auc_random {:7.3f} {:6d} {:6d} {:2s} {} {}\n'\
                              .format( auc_random, len(positive_nbrdists), len(negative_nbrdists_random),
                                       chains, epitope, nbrdist_tag_suffix ) )
            all_aucs_random[nbrdist_tag_suffix].append( auc_random )
            all_aucs_random[nbrdist_tag_suffix_w_chains].append( auc_random )



            ## now versus TCRs from other epitopes
            negative_nbrdists_others = [ x[ ii_nbrdist_tag ] for e2 in other_epitopes for x in all_tcrs[e2] if not x[-1] ]

            auc_others, xvals_others, yvals_others = tcr_distances.compute_auc( positive_nbrdists, negative_nbrdists_others,
                                                                                sign_factor = sign_factor )

            outlog.write( 'auc_others {:7.3f} {:6d} {:6d} {:2s} {} {}\n'\
                          .format( auc_others, len(positive_nbrdists), len(negative_nbrdists_others),
                                   chains, epitope, nbrdist_tag_suffix ) )
            all_aucs_others[nbrdist_tag_suffix].append( auc_others )
            all_aucs_others[nbrdist_tag_suffix_w_chains].append( auc_others )

            ## plotting ############################################################################################
            # which plotno?
            assert nbrdist_tag_suffix in nbrdist_tags_for_plotting
            ii_nbrdist_tag_suffix = nbrdist_tags_for_plotting.index( nbrdist_tag_suffix )
            plotno = ii_epitope * ncols + ii_nbrdist_tag_suffix + 1

            ## the kde-smoothed distributions
            plt.figure(2*ii_chains+1,figsize=(fig_width,fig_height))
            plt.subplot( nrows,ncols,plotno )

            mn = min( positive_nbrdists + negative_nbrdists_random + negative_nbrdists_others )
            mx = max( positive_nbrdists + negative_nbrdists_random + negative_nbrdists_others )

            if chains in fake_chains:
                continue ## the next line will fail if all positive_nbrdists are 0
            positive_density = gaussian_kde( positive_nbrdists )
            negative_density_random = gaussian_kde( negative_nbrdists_random )
            if other_epitopes: negative_density_others = gaussian_kde( negative_nbrdists_others )

            #line_style = ['--',':','-'][ ii_chains ]
            line_style = '-'

            xs = np.linspace( mn, mx, 100 )
            ys0 = positive_density(xs)
            if other_epitopes: ys1 = negative_density_others(xs)
            ys2 = negative_density_random(xs)
            plt.plot( xs, ys0, line_style, c='r' )
            if other_epitopes: plt.plot( xs, ys1, line_style, c='g' )
            plt.plot( xs, ys2, line_style, c='b' )

            plt.title('{} {}'.format(epitope, nbrdist_labels_for_plotting[ ii_nbrdist_tag_suffix ] ) )
            plt.suptitle('Smoothed {} {} distributions\nRed= epitope-specific TCRs\ngreen= TCRs from other epitopes\nblue= random TCRs'.format(chains, 'nbrdens' if 'dens' in nbrdist_tag else 'nbrdist' ),y=1.0)
            plt.subplots_adjust( hspace = 0.3, bottom = bottom_margin, right = right_margin,
                                 top = top_margin, left = left_margin )

            if 'dens' in nbrdist_tag:
                ymn,ymx = plt.ylim()
                ymx = min( ymx, 4*max(ys0) )
                plt.ylim( (ymn,ymx) )

            plt.figure(2*ii_chains+2,figsize=(fig_width,fig_height))
            plt.subplot( nrows,ncols,plotno )

            #line_style = ['--',':','-'][ ii_chains ]
            line_style = '-'

            if other_epitopes:
                plt.plot( xvals_others, yvals_others, line_style, c='g',label='others {}'.format(int(100*auc_others) ))
            plt.plot( xvals_random, yvals_random, line_style, c='b',label='random {}'.format(int(100*auc_random) ))
            plt.title('{} {}'.format(epitope, nbrdist_labels_for_plotting[ ii_nbrdist_tag_suffix ] ) )

            plt.suptitle('{} ROC curves (AUROCs in legend)\nRed= epitope-specific TCRs\ngreen= TCRs from other epitopes\nblue= random TCRs'.format(chains),y=1.0)
            plt.subplots_adjust( hspace = 0.3, bottom = bottom_margin, right = right_margin, top = top_margin,
                                 left = left_margin )
            plt.legend(loc='lower right',fontsize='small')


            if nbrdist_tag_suffix == single_nbrdist_tag_for_plotting:
                save_epitope_rocs[ epitope ] = ( xvals_random, yvals_random, auc_random )
                save_epitope_nbrdists[ epitope ] = positive_nbrdists

                ## now put everything into a single plot
                plt.figure(7,figsize=(fig_width_summary,fig_height))
                plotno = ii_epitope * ncols_single + ii_chains + 1
                plt.subplot( nrows_single, ncols_single, plotno )

                ## density:
                plt.plot( xs, ys0, line_style, c='r' )
                if other_epitopes: plt.plot( xs, ys1, line_style, c='g' )
                plt.plot( xs, ys2, line_style, c='b' )
                plt.yticks([],[])
                if paper_supp:
                    locs,labels = plt.xticks()
                    print 'locs,labels:',locs,labels
                    newlocs,newlabels = [],[]
                    for loc in locs:
                        if abs(loc/100.0 - int(loc/100.0))<1e-3:
                            newlocs.append(loc)
                            newlabels.append(str(int(loc+1e-3)))
                    plt.xticks(newlocs,newlabels)
                plt.title('{} {}'.format(epitope,greek_seg(chains)))

                if 'dens' in nbrdist_tag:
                    ymn,ymx = plt.ylim()
                    ymx = min( ymx, 4*max(ys0) )
                    plt.ylim( (ymn,ymx) )
                plotno = ii_epitope * ncols_single + ii_chains + 4
                plt.subplot( nrows_single, ncols_single, plotno )

                ## ROC:
                if other_epitopes:
                    plt.plot( xvals_others, yvals_others, line_style, c='g',label='others {}'.format(int(100*auc_others) ))
                plt.plot( xvals_random, yvals_random, line_style, c='b',label='random {}'.format(int(100*auc_random) ))
                plt.title('{} {}'.format(epitope,greek_seg(chains)))
                plt.legend(loc='lower right',fontsize='small')

                if paper_supp: ## special ticks
                    plt.xticks( [0,.2,.4,.6,.8,1.0], ['0','.2','.4','.6','.8','1'] )
                    plt.yticks( [0,.2,.4,.6,.8,1.0], ['0','.2','.4','.6','.8','1'] )

                ## suptitle
                if paper_supp:
                    title = 'Smoothed NN-distance histograms (left columns) and ROC curves (right columns, AUROCs given in legend)\nRed= epitope-specific TCRs, green= TCRs from other epitopes, blue= random background TCRs'
                    plt.suptitle( title, y=( top_margin + 4.0 )/5.0 )
                else:
                    plt.suptitle('Smoothed {} histograms (left cols) and ROC curves (right cols, AUROCs in legend)\nRed= epitope-specific TCRs, green= TCRs from other epitopes, blue= random TCRs'\
                                 .format( nbrdist_labels_for_plotting[ ii_nbrdist_tag_suffix ] ), y=1.0)

                plt.subplots_adjust( hspace = 0.3, bottom = bottom_margin, right = 0.98, top = top_margin, left = 0.05 )


for nbrdist_tag_suffix in sorted( all_aucs_others.keys()):
    l_random = all_aucs_random[nbrdist_tag_suffix]
    l_others = all_aucs_others[nbrdist_tag_suffix]
    outlog.write( 'avg_auc_random {:7.3f} {}\n'.format( sum( l_random ) / len( l_random ), nbrdist_tag_suffix ) )
    outlog.write( 'avg_auc_others {:7.3f} {}\n'.format( sum( l_others ) / len( l_others ), nbrdist_tag_suffix ) )

outlog.close()


for ii_roc in range(2):
    for ii_chains,chains in enumerate(ABs):
        figno = 2*ii_chains + ii_roc + 1
        figtype = 'roc' if ii_roc else 'nbrdists'
        pngfile = '{}_{}_{}.png'.format(outfile_prefix, figtype, chains )
        print 'making:',pngfile
        plt.figure(figno)
        plt.savefig(pngfile)
        if ii_roc==0:
            util.readme(pngfile, """KDE-smoothed histograms of different {} nbrdist measures, comparing each epitope-specific repertoire (red)
            to TCRs from the other repertoires (green) and to random TCRs (blue).<br><br>
            """.format(chains))
        else:
            util.readme( pngfile, """ROC curves showing true- (y-axis) versus false-positives (x-axis) as the sorting metric increases. Legends give the area under the curve (AUROC) values, ranging from 100 (perfect discrimination, curve goes straight
            up then straight across), to 50 (random), to 0 (total failure, all negatives come before all positives, curve goes straight across then straight up).<br><br>
            """.format(chains))


figno = 7
pngfile = '{}_summary.png'.format(outfile_prefix)
plt.figure(figno)
print 'making:',pngfile
plt.savefig(pngfile)
util.readme( pngfile, """KDE-smoothed nbrdist10P histograms for alpha (col 1), beta (col 2), and alpha-beta (col 3). ROC curves for
the same in columns 4-6.<br><br>""")


## make a figure showing all the nbrdist distributions (left col) and all the ROC curves
figno = 10 ## hope this is larger than the others
plt.figure(figno,figsize=(12,8))

plt.subplot(121) ## the nbrdist distributions
mn = 0
mx = max( ( max(x) for x in save_epitope_nbrdists.values() ) )


for epitope,color in zip( epitopes, html_colors.get_rank_colors_no_lights(len(epitopes))):
    nbrdists = save_epitope_nbrdists[epitope]
    positive_density = gaussian_kde( nbrdists )

    xs = np.linspace( mn, mx, 100 )
    ys0 = positive_density(xs)

    plt.plot( xs,ys0,c=color,label='{} ({:.1f})'.format(epitope,sum(nbrdists)/len(nbrdists)) )

plt.legend(fontsize=9,frameon=False,loc='best')

plt.subplot(122) ## now the ROC curves


for epitope,color in zip( epitopes, html_colors.get_rank_colors_no_lights(len(epitopes))):
    xvals,yvals,auc = save_epitope_rocs[epitope]
    plt.plot( xvals,yvals,c=color,label='{} ({:.3f})'.format(epitope,auc))

plt.xlim((0,1.0))
plt.ylim((0,1.0))
plt.legend(fontsize=9,frameon=False,loc='best')

pngfile = '{}_nbrdist_roc_superpositions.png'.format(outfile_prefix)
print 'making:',pngfile
plt.savefig(pngfile)
util.readme( pngfile, """
Superimposed KDE-smoothed NNdistance ({}) distributions (left) and ROC curves (right) for all epitopes (paired-chain analyses).<br><br>
""".format(single_nbrdist_tag_for_plotting))


if show:
    plt.show()

