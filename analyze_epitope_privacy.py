import sys
from basic import *
import tcr_distances
import parse_tsv
import numpy as np
from scipy.cluster import hierarchy
from scipy.spatial import distance
import util
import html_colors
from all_genes import all_genes


with Parser(locals()) as p:
    #p.str('args').unspecified_default().multiple().required()
    p.str('organism').required()
    p.str('clones_file').required()
    p.int('nrepeat').default(100)
    p.int('num_random_trees').default(3)
    p.int('nbrdist_percentile').default(25)
    p.float('tree_height_inches').default(3.0)
    p.float('max_leaf_font_size').default(10.0)
    p.float('min_leaf_font_size').default(6.0)
    p.multiword('epitopes').cast(lambda x:x.split())
    p.multiword('all_chains').cast(lambda x:x.split()).default("A B AB")
    p.flag('scale_bars_in_key')       # --flag_arg  (no argument passed)
    p.flag('no_mouse_labels')       # --flag_arg  (no argument passed)
    p.flag('paper_figs')       # --flag_arg  (no argument passed)
    p.flag('constant_seed')       # --flag_arg  (no argument passed)
    p.str('distance_params')
    p.str('outfile_prefix')

if constant_seed: random.seed(1)

distance_params = tcr_distances.DistanceParams( config_string = distance_params )

if outfile_prefix is None:
    outfile_prefix = clones_file[:-4]

import matplotlib
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#recompute_nbrdists = True

clones_file_with_nbrdists = '{}_nbrdists.tsv'.format(clones_file[:-4])
assert exists( clones_file_with_nbrdists )

## read the epitope-specific TCRs
all_tcrs = parse_tsv.parse_tsv_file( clones_file, ['epitope','subject'], ['va_genes','vb_genes','cdr3a','cdr3b'],
                                     save_l = True ) ## last element will be the full parse_tsv_line dictionary

if not epitopes:
    epitopes = all_tcrs.keys()[:]
    epitopes.sort()


Log('reading {}'.format(clones_file_with_nbrdists))

nbrdist_tag_suffix = '_wtd_nbrdist{}'.format(nbrdist_percentile)
nbrdist_tags = [ x+'_'+y+nbrdist_tag_suffix for x in epitopes for y in all_chains ]

all_nbrdists = parse_tsv.parse_tsv_file( clones_file_with_nbrdists, ['epitope','subject'], nbrdist_tags )

## BEGIN stolen from compute_distances.py
Log( 'precomputing v-region distances')
rep_dists = tcr_distances.compute_all_v_region_distances( organism, distance_params )
Log( 'done precomputing v-region distances' )




def compute_mouse_distances_fast( reorder, mice, mouse_indices, D ):
    inter_mouse_distances = []
    intra_mouse_distances = []
    inter_mouse_distance_averages = []
    intra_mouse_distance_averages = []

    for m1 in mice:
        for m2 in mice:
            if m1>m2: continue
            avg_dist = 0.0
            count=0
            for i1 in mouse_indices[m1]:
                for i2 in mouse_indices[m2]:
                    dist = D[ reorder[i1], reorder[i2] ]
                    if m1==m2:
                        if i1>=i2: continue
                        intra_mouse_distances.append( dist )
                    else:
                        inter_mouse_distances.append( dist )
                    avg_dist += dist
                    count += 1
            if not count: count=1
            if m1==m2:
                intra_mouse_distance_averages.append( avg_dist/count )
            else:
                inter_mouse_distance_averages.append( avg_dist/count )
    return intra_mouse_distances, inter_mouse_distances, intra_mouse_distance_averages, inter_mouse_distance_averages


def get_mouse_average_nbrdists( reorder, mice, mouse_indices, nbrdists ): ## returns a list ordered same as mice
    avgs = []
    for m in mice:
        avgs.append( sum( ( nbrdists[ reorder[ x ] ] for x in mouse_indices[m] ) ) / len( mouse_indices[m] ) )
    return avgs

def get_Z( val, mean, sdev ):
    return (val-mean)/sdev if sdev else 0.0

nrows = len(epitopes)
ncols = 1+num_random_trees

fig_height = nrows * tree_height_inches + 1.5
fig_width = 12

bottom_margin = 0.5 / fig_height
top_margin = ( fig_height-1.0 ) / fig_height

#print 'fig_height:',fig_height

plt.figure(1,figsize=(fig_width,fig_height))
plt.subplots_adjust(left=0.4,top=top_margin, bottom=bottom_margin,hspace=0.3)

plotno=0
figcounter = 3 ## reserve figures 1 and 2 for the old plots

epitope_zps = {}
epitope_zp2s = {}

for epitope in epitopes:
    tcrs = []
    tcr_infos = []
    mouse_indices = {}
    num_mouse_tcrs = {}

    mice = all_tcrs[epitope].keys()[:]
    mice.sort()

    for mouse in mice:
        mouse_indices[ mouse ] = []
        num_mouse_tcrs[ mouse ] = len(all_tcrs[epitope][mouse] )
        assert num_mouse_tcrs[ mouse ] == len( all_nbrdists[epitope][mouse] )
        for l in all_tcrs[epitope][mouse]:
            ## ( va_reps, vb_reps, cdr3a, cdr3b )
            index = len( tcrs )
            mouse_indices[mouse].append( index )
            va_reps = frozenset( ( all_genes[organism][x].rep for x in l[0].split(';') ) )
            vb_reps = frozenset( ( all_genes[organism][x].rep for x in l[1].split(';') ) )
            tcrs.append( ( va_reps, vb_reps, l[2], l[3] ) )
            assert len(l) == 5
            dict_from_parse_tsv_line = l[4]
            tcr_infos.append( dict_from_parse_tsv_line )

    util.assign_label_reps_and_colors_based_on_most_common_genes_in_repertoire( tcr_infos, organism )
    ## now stored as va_label_rep, jb_label_rep


    num_tcrs = len( tcrs )
    num_mice = len( mice )

    if num_mice<=1:
        Log( `( 'only one mouse?',epitope,mouse_indices.keys() )` )
        continue

    for chains in all_chains:

        my_nbrdist_index = nbrdist_tags.index( epitope+'_'+chains+nbrdist_tag_suffix )
        nbrdists = []
        for mouse in mice:
            for nbrdistl in all_nbrdists[epitope][mouse]:
                nbrdists.append( float( nbrdistl[ my_nbrdist_index ] ) )

        Log('START filling distance matrix {} {}'.format(epitope,chains))

        D = np.zeros( ( num_tcrs, num_tcrs ) )
        for i1,t1 in enumerate( tcrs ):
            for i2,t2 in enumerate( tcrs ):
                if i2<=i1: continue
                dist = tcr_distances.compute_distance( tcrs[ i1 ], tcrs[ i2 ], chains, rep_dists, distance_params )
                D[i1,i2]= dist
                D[i2,i1]= dist
        Log('DONE filling distance matrix {} {}'.format(epitope,chains))

        reorder = range(len(tcrs))

        intras,inters,real_intra_avs,real_inter_avs = compute_mouse_distances_fast( reorder, mice, mouse_indices, D )

        save_inters = [ real_inter_avs ]


        assert len(intras) + len(inters) == ( num_tcrs * (num_tcrs-1) ) /2

        avg_inter = sum(inters)/len(inters)
        avg_intra = sum(intras)/len(intras)

        rand_inters = []
        rand_intras = []
        rand_inter_avs = []
        rand_intra_avs = []

        save_rand_inters = []

        real_avg_nbrdists = get_mouse_average_nbrdists( reorder, mice, mouse_indices, nbrdists )

        all_rand_avg_nbrdists = []

        for r in range(nrepeat):
            random.shuffle( reorder )
            intras,inters,intra_avs,inter_avs = compute_mouse_distances_fast( reorder, mice, mouse_indices, D )

            if r<num_random_trees:
                save_inters.append( inter_avs )

            rand_intras.append( sum(intras)/len(intras))
            rand_inters.append( sum(inters)/len(inters))
            rand_intra_avs.append( intra_avs )
            rand_inter_avs.append( inter_avs )
            all_rand_avg_nbrdists.append( get_mouse_average_nbrdists( reorder, mice, mouse_indices, nbrdists ) )

        mean_intra,sdev_intra = get_mean_and_sdev( rand_intras )
        mean_inter,sdev_inter = get_mean_and_sdev( rand_inters )

        nlower1 = len( [ x for x in rand_intras if x < avg_intra ] )
        nlower2 = len( [ x for x in rand_inters if x < avg_inter ] )

        z_intra = ( avg_intra-mean_intra)/sdev_intra
        z_inter = ( avg_inter-mean_inter)/sdev_inter
        assert abs( z_intra + z_inter )<1e-3

        p_intra = (1.0*nlower1)/nrepeat
        print 'rep   Z: {:9.3f} P: {:9.6f} intra {:9.3f} inter {:9.3f} ntcrs: {:3d} nmice: {:3d} {} {}'\
            .format( z_intra, p_intra, avg_intra, avg_inter,
                     num_tcrs, num_mice,
                     chains, epitope )

        if chains == 'AB':
            epitope_zps[epitope] = ( z_intra, p_intra )

        ## look for mice with surprisingly low intras
        mouse_zp = {}
        for ii,real_av in enumerate( real_intra_avs ):
            rand_avs = [ x[ii] for x in rand_intra_avs ]
            mean,sdev = get_mean_and_sdev( rand_avs )
            nlower = len( [ x for x in rand_avs if x < real_av ] )
            Z = get_Z( real_av, mean, sdev )
            p = ( 1.0 * nlower ) / nrepeat if sdev else .5
            mouse_zp[ mice[ii] ] = (Z,p)
            print 'mouse Z: {:9.3f} P: {:9.6f} {:2d} {} {} {}'\
                .format( Z, p,
                         len(mouse_indices[mice[ii]]), mice[ii] , chains, epitope )

        ## look at mouse nbrdists
        rand_Zsums = [0.0]*nrepeat
        Zsum = 0.0
        mouse_zp2 = {}
        for ii, mouse in enumerate(mice):
            real_av = real_avg_nbrdists[ii]
            rand_avs = [x[ii] for x in all_rand_avg_nbrdists ]
            mean,sdev = get_mean_and_sdev( rand_avs )
            nlower = len( [ x for x in rand_avs if x < real_av ] )
            nhigher = len( [ x for x in rand_avs if x > real_av ] )
            Z = get_Z( real_av, mean, sdev )
            plower = ( 1.0 * nlower ) / nrepeat if sdev else .5
            phigher = ( 1.0 * nhigher ) / nrepeat if sdev else .5
            mouse_zp2[ mouse ] = (Z,min(plower,phigher))
            print 'mouse-nbrdist Z: {:9.3f} P: {:9.6f} {:2d} {} {} {}'\
                .format( Z, min( plower, phigher ),
                         len(mouse_indices[mice[ii]]), mice[ii] , chains, epitope )
            Zsum += abs(Z)

            for jj,rand_av in enumerate( rand_avs ):
                randZ = get_Z( rand_av, mean, sdev )
                rand_Zsums[jj] += abs(randZ)

        ## compare Zsum to rand_Zsums
        mean,sdev = get_mean_and_sdev( rand_Zsums )
        ZZ = get_Z( Zsum, mean, sdev )
        nhigher = len( [ x for x in rand_Zsums if x > Zsum ] )
        Zp = ( 1.0 * nhigher ) / nrepeat if sdev else .5
        if chains == 'AB':
            epitope_zp2s[epitope] = ( ZZ, Zp )


        print 'rep-nbrdist   Z: {:9.3f} P: {:9.6f} ntcrs: {:3d} nmice: {:3d} {} {}'\
            .format( ZZ, Zp, num_tcrs, num_mice, chains, epitope )


        mouse_pairs = []
        for m1 in mice:
            for m2 in mice:
                if m1>=m2: continue
                mouse_pairs.append( ( m1,m2) )
        assert len(mouse_pairs) == len( real_inter_avs )

        for ii,real_av in enumerate( real_inter_avs ):
            rand_avs = [ x[ii] for x in rand_inter_avs ]
            mean,sdev = get_mean_and_sdev( rand_avs )
            nlower = len( [ x for x in rand_avs if x < real_av ] )
            mouse1,mouse2 = mouse_pairs[ii]
            print 'mice  Z: {:9.3f} P: {:9.6f} {:2d} {:2d} {} {} {} {}'\
                .format( get_Z( real_av, mean, sdev ), ( 1.0 * nlower ) / nrepeat if sdev else .5,
                         len(mouse_indices[mouse1]), len(mouse_indices[mouse2]),
                         mouse1, mouse2, chains, epitope )


        if chains == 'AB': ## make some plots

            label_height_inches = float( tree_height_inches ) / num_mice
            leaf_font_size = max( min_leaf_font_size, min( max_leaf_font_size,
                                                           int( floor( 0.5+label_height_inches * 72.0 * 0.85 ) ) ) )


            for ii in range(1+num_random_trees):
                inters = save_inters[ii]
                mouse_D = np.zeros( ( num_mice,num_mice ) )

                #print len(inters), len(mouse_pairs)
                assert len(inters) == len(mouse_pairs)

                for i1,m1 in enumerate(mice):
                    for i2,m2 in enumerate(mice):
                        if m1>=m2: continue
                        mouse_D[i1,i2] = inters[0]
                        mouse_D[i2,i1] = inters[0]
                        del inters[0]
                assert not inters

                y = distance.squareform( mouse_D, checks=True )
                assert len(y) == ( num_mice*(num_mice-1) )/2

                Z = hierarchy.average( y )
                #Z = hierarchy.average( mouse_D )
                c,coph_dists = hierarchy.cophenet(Z,y)

                # leaves = hierarchy.leaves_list( Z )
                # tree_ordered_mice = [ mice[x] for x in leaves ]
                # #print 'leaves:',leaves
                print 'coph:',epitope,ii,c
                #print Z[:3]

                plotno += 1
                ax = plt.subplot( nrows, ncols, plotno )
                def get_stars_from_pvalue( p ):
                    if p<=0.001:
                        return '***'
                    elif p<=0.01:
                        return '**'
                    elif p<=0.05:
                        return '*'
                    else:
                        return ''
                if ii==0:
                    save_Z = Z
                    labels = [ '{} ({})  {:3s} Z: {:4.1f} {:4.1f} P: {:5.3f} {:5.3f}'\
                               .format( x, num_mouse_tcrs[x],
                                        get_stars_from_pvalue( min([mouse_zp[x][1], 1.0-mouse_zp[x][1], mouse_zp2[x][1] ] ) ),
                                        mouse_zp[x][0], mouse_zp2[x][0],
                                        min(mouse_zp[x][1],1.0-mouse_zp[x][1]), mouse_zp2[x][1] )
                               for x in mice ]
                    save_labels = labels
                else:
                    labels = None
                hierarchy.dendrogram( Z, ax=ax, orientation='right',labels = labels, leaf_font_size=leaf_font_size )

                #if ii==0:
                #    plt.ylabel(epitope)

                ## change to -1*z_intra so higher values reflect greater heterogeneity
                ## TODO: check sign of Z scores for individual subjects...
                info = '{}  Z: {:4.1f}  P: {:5.3f}'.format( epitope, -1*z_intra, min(p_intra,1.-p_intra) )
                if plotno <= ncols and ii>0:
                    plt.title('rand{}'.format(ii))
                elif ii==0:
                    plt.title(info)


            ############################################################################################################
            ## now we are going to secretly make a single figure for this epitope
            figcounter += 1

            fig_width = 2. if paper_figs else 12.
            fig_height = 6. if paper_figs else 12.

            fig = plt.figure(figcounter,figsize=(fig_width,fig_height))

            if paper_figs:
                tree_yfraction = 0.8
                tree_xfraction = 0.5
                ## [left,bottom,width,height]
                ax = fig.add_axes( [ 0.0, 0.0, 1.0 - tree_xfraction, tree_yfraction ] )
            else:
                ax = plt.subplot( 131 ) ## left plot, will have all the mouse-bars

            Z = save_Z ## recover the tree for the observed distances
            leaves = hierarchy.leaves_list( Z )
            tree_ordered_mice = [ mice[x] for x in leaves ]

            ## we need to figure out how to color all the different genes
            rep_colors = {}
            reptypes = []
            ordered_reps_list = {}
            for ab in 'ab':
                for vj in 'vj':
                    reptype = vj+ab
                    reptypes.append( reptype )
                    counts = {}
                    for info in tcr_infos:
                        rep = info['{}_label_rep'.format(reptype)]
                        counts[rep] = counts.get(rep,0)+1
                    l = [ ( counts[x],x) for x in counts ]
                    l.sort()
                    l.reverse()
                    colors = html_colors.get_rank_colors_no_lights( len(l) )
                    ordered_reps_list[reptype] = [ (x[0],x[1],y) for x,y in zip(l,colors)] #(count,rep,color)
                    # for (count,rep),color in zip( l, colors ):
                    #     rep_colors[rep] = color

            ## each mouse is going to get four bars: va ja vb jb
            ## each bar for mouse i is going to run in y from (i+pad) to (i+1-pad)
            pad = 0.1
            mouse_bars_height = 1.0-2*pad
            bar_width = 0.9 if paper_figs else 0.8
            bars = []
            ylocs = []
            ylabels = []
            for ii, mouse in enumerate( tree_ordered_mice ):
                ylocs.append( ii+0.5 )
                ylabels.append( '' if no_mouse_labels else save_labels[ mice.index(mouse) ] )
                infos = [ tcr_infos[x] for x in mouse_indices[mouse] ]
                total = float( len(infos))
                for jj,reptype in enumerate( reptypes ):
                    mouse_reps = [ x[reptype+'_label_rep'] for x in infos ]

                    bar_bottom = ii+1-pad ## the TOP for the top bar, we shift it down before we append
                    bar_left = jj-bar_width/2.0

                    for (repcount,rep,color) in ordered_reps_list[ reptype ]:
                        frac = mouse_reps.count(rep) / total
                        if frac:
                            bar_height = frac * mouse_bars_height
                            bar_bottom -= bar_height
                            bars.append( ( bar_left, bar_height, bar_bottom, color ) )
                    assert abs( bar_bottom - (ii+pad))<1e-3

            ## add some labels
            if paper_figs:
                for jj,reptype in enumerate( reptypes ):
                    top_of_top_bar = len(tree_ordered_mice)-pad
                    bar_left = jj-bar_width/2.0

                    label = { 'va': r'V$\alpha$',
                              'ja': r'J$\alpha$',
                              'vb': r'V$\beta$',
                              'jb': r'J$\beta$' }[ reptype ]

                    font_family="Droid Sans Mono"
                    fontdict = {'family': font_family, 'size': 12 }
                    plt.text( bar_left+bar_width/2., top_of_top_bar + pad, label,
                              va='bottom', ha='center', fontdict=fontdict )


            plt.bar( width  = bar_width,
                     left   = [x[0] for x in bars ],
                     height = [x[1] for x in bars ],
                     bottom = [x[2] for x in bars ],
                     color  = [x[3] for x in bars ],
                     edgecolor = 'none' )

            plt.ylim( (0.0,len(mice) ) )

            if paper_figs:
                plt.axis('off')
            else:
                plt.yticks( ylocs, ylabels )
                plt.xticks( range(4), reptypes )


            print 'bars xlim:',plt.xlim()

            ## now add the full-repertoire bars at the top if paper_figs
            if paper_figs:
                ax = fig.add_axes( [ 0.0, tree_yfraction, 1.0 - tree_xfraction, 1.0 - tree_yfraction ] )
                bars = []
                for jj,reptype in enumerate( reptypes ):

                    bar_bottom = 1-pad ## the TOP for the top bar, we shift it down before we append
                    bar_left = jj-bar_width/2.0

                    total = len(tcr_infos)
                    for (repcount,rep,color) in ordered_reps_list[ reptype ]:
                        frac = float(repcount) / total
                        assert frac
                        bar_height = frac * mouse_bars_height
                        bar_bottom -= bar_height
                        bars.append( ( bar_left, bar_height, bar_bottom, color ) )


                    if reptype == 'jb': ## add a label on the right
                        font_family="Droid Sans Mono"
                        fontdict = {'family': font_family, 'size': 12 }
                        plt.text( bar_left+1.3*bar_width, bar_bottom + 0.5 * mouse_bars_height,
                                  'All\nsubjects', va='center', ha='left', fontdict=fontdict )


                plt.bar( width  = bar_width,
                         left   = [x[0] for x in bars ],
                         height = [x[1] for x in bars ],
                         bottom = [x[2] for x in bars ],
                         color  = [x[3] for x in bars ],
                         edgecolor = 'none' )
                plt.ylim((-0.2,1))
                plt.axis('off')


            ## now the dendrogram
            if paper_figs:
                ax = fig.add_axes( [ 1.0 -tree_xfraction, 0.0, tree_xfraction, tree_yfraction ] )
            else:
                ax = plt.subplot( 132 ) ## middle plot, will have the dendrogram

            hierarchy.dendrogram( Z, ax=ax, orientation='right',labels = ['']*len(mice), link_color_func = None )
            plt.axis('off')

            if not paper_figs:
                ## now we're going to add a third column which gives the names of the genes
                plt.subplot(133)
                plot_height = fig_height*0.85 ## not used if paper_figs

                max_reps_to_show = 30
                bottoms = []
                heights = []
                colors = []
                left = -1*bar_width/2.
                for ii,reptype in enumerate(reptypes):
                    numreps = min( max_reps_to_show, len( ordered_reps_list[ reptype ] ) )
                    inches_per_rep = plot_height / ( 4.0 * numreps )
                    rep_font_size = min( 10.0, max( 6.0, int( floor( 0.5+ inches_per_rep * 72.0 * 0.85 ) ) ) )
                    print 'rep_font_size:',reptype,rep_font_size
                    topcount = ordered_reps_list[ reptype ][0][0]
                    pad = 0.05
                    reptype_bottom = 3-ii + pad
                    reptype_top = reptype_bottom +1. - 2*pad
                    rep_step = ( reptype_top - reptype_bottom ) / float(numreps)
                    ## in this loop we will at least draw the labels, even if we don't do the bars...
                    for jj,(count,rep,color) in enumerate( ordered_reps_list[ reptype ][:numreps] ):
                        middle = reptype_top - jj*rep_step - rep_step/2.
                        height = rep_step * float(count)/topcount
                        height = 0.9*rep_step
                        if not scale_bars_in_key:
                            bottoms.append( middle - height/2.)
                            heights.append( height )
                            colors.append( color )
                        num_tcrs_this_epitope = len(tcr_infos )
                        label = '{:4.1f}% {}'.format( ( 100.0*count)/num_tcrs_this_epitope, rep )
                        plt.text( left+1.1 * bar_width, middle, label, ha='left', va='center', color=color,
                                  fontsize = rep_font_size )
                    if numreps<len(ordered_reps_list[reptype] ):
                        num_missed = len(ordered_reps_list[reptype]) - numreps
                        plt.text( left+1.1 * bar_width, reptype_bottom - rep_step/2.,
                                  'and {} more'.format(num_missed),
                                  ha='left', va='center', color=color,
                                  fontsize = rep_font_size )

                    if scale_bars_in_key:
                        ## now make the bars themselves to scale their heights by rep counts
                        bar_top = reptype_top # this will change
                        total = len(tcr_infos)
                        for (count,rep,color) in ordered_reps_list[ reptype ]:
                            height = ( reptype_top - reptype_bottom ) * float(count) / total
                            bar_top -= height
                            bottoms.append( bar_top ) ## we've already shifted it down
                            heights.append( height )
                            colors.append( color )


                plt.bar( width  = [bar_width]*len(heights),
                         left = [left]*len(heights), height = heights, bottom = bottoms, color = colors,
                         edgecolor = 'none' )


                plt.axis('off')

            if not paper_figs:
                plt.subplots_adjust( left=0.1 if no_mouse_labels else 0.4,top=0.92, bottom=0.08, right = 0.9, wspace = 0.03 )
                plt.suptitle( 'mouse tree for {} showing gene composition of mouse repertoires'.format(epitope))

            ## now save the figure
            pngfile = '{}_{}_subject_tree.png'.format( outfile_prefix, epitope )
            print 'making:',pngfile
            plt.savefig( pngfile )
            util.readme( pngfile, """
This is a hierarchical clustering dendrogram of the subjects for the {} repertoire. The distance between a pair of subjects is defined to be the average distance between receptors in one subject and receptors in the other. The gene frequency composition of each subject is shown in the four stacks of colored bars in the middle (to the left of the tree, which is drawn with thin blue lines). To the right of the tree is a key for the gene segment coloring schemes which also shows the frequencies of each gene in the combined repertoire.
""".format(epitope))

            plt.close()


if plotno==0:
    Log('no multi-subject epitopes??')
    exit()

if paper_figs: ## just the individual tree images
    exit()

plt.suptitle('mouse trees based on avg TCR-TCR distance between mice',x=0.5,y=1.0)
pngfile = '{}_subject_trees.png'.format( outfile_prefix )
print 'making:',pngfile
plt.savefig( pngfile )
util.readme( pngfile, """This figure shows hierachical clustering (average linkage) trees of the different
mice (subjects), where the distance between two mice is defined as the average distance between TCRs in one
mouse and TCRs in the other. This could help us to interpret any high repertoire-heterogeneity Z-scores in the
"_subject_heterogeneity.png" bar plot earlier on, for example if it turns out that there are distinct sub-clusters of mice
with similar repertoires where the sub-clusters are quite different from one another.
See the earlier plot's readme for the explanation of how we assign a Z-score
and P-value to each repertoire reflecting the degree of repertoire heterogeneity across different mice. These
are the values shown after the epitope name above each row of trees. For this plot, more positive Z-scores mean
the intra-mouse distances are smaller than the inter-mouse distances (ie, more heterogeneity),
while negative Z-scores mean the opposite.
The first column of trees are the observed ones, while the next three are based on datasets in which the TCRs are
randomly reshuffled between mice (preserving each mouse's TCR count).
<br>
To provide a little more info about each mouse's individual repertoire, I looked at the intra-mouse distances for
TCRs from that mouse versus distances between randomly sampled sets of the same number of TCRs and assigned a Z-score
and P-value to each mouse (shown after the mouse names labeling each leaf in the trees). A negative Z-score means that
the mouse's TCRs are less diverse than average, while a positive Z-score means that the mouse's TCRs are more diverse
(more distant from one another). The P-value is a rough measure of how surprising these deviations from randomly
chosen TCRs actually are. The number of TCRs for each mouse is shown in parentheses after the mouse name. And *'s are
used to flag significant P-values (0.05, 0.01, 0.001).
""")


## now make a summary plot showing the epitope Z-scores and p-values
plt.figure(2,figsize=(14,8))
l = [ ( -1*epitope_zps[x][0], x ) for x in sorted( epitope_zps.keys()) ]
l.sort()

lefts = range( len(l) )
heights = [x[0] for x in l ]

plt.bar( lefts, heights, width=0.8 )
plt.subplots_adjust(bottom=0.25)
plt.xticks( [x+0.4 for x in lefts], ['{} P: {:.3f}'.format( x[1], epitope_zps[x[1]][1] ) for x in l], rotation='vertical', fontsize=9 )
plt.ylabel('Z-score\n(abs. values greater than 2 or 3 start to look significant')
plt.suptitle('Z-scores for intra- versus inter-mouse distances\nLarger Z means intra-mouse distances are smaller than inter-mouse distances\ni.e. greater heterogeneity across mice')
pngfile = '{}_subject_heterogeneity.png'.format( outfile_prefix )
print 'making:',pngfile
plt.savefig( pngfile )
util.readme( pngfile, """This summary plot is aimed at answering the question: are some mice (or humans) sampling from a different repertoire than others?
There is a companion plot later on (_subject_trees.png) which adds a little more detail. The idea behind the analysis is to compute
distances between all TCRs in the dataset, and then split those distances into two sets: intra-mouse distances (distances between
TCRs from the same mouse), and inter-mouse distances (distances between TCRs from different mice). If mice are sampling from different
repertoires then the intra-mouse distances might be expected to be smaller than the inter-mouse distances. To assess significance of an
observed difference, I use randomization: the TCRs are randomly re-assigned to the mice 1000 times (preserving the size of each mouse's
individual repertoire), and each time the intra- and inter-mouse
distances are calculated and the mean values of each are stored. To assess the significance of the actually observed mean intra- and inter-
mouse distances, I compute the mean and standard deviation of the 1000 shuffled means and use those to assign a Z-score to the observed
values. Larger absolute values of the Z-score suggest greater significance. To get an approximate P-value I keep track of how often the
random mean intra-mouse distance is less than the observed mean intra-mouse distance, and divide that number by 1000 (the number of random
reshufflings). This bar plot shows the Z-scores; the P-values are listed below the epitopes at the bottom.
""")



## now make a summary plot showing the epitope Z-scores and p-values
plt.figure(3,figsize=(14,8))
l = [ ( epitope_zp2s[x][0], x ) for x in sorted( epitope_zp2s.keys()) ]
l.sort()

lefts = range( len(l) )
heights = [x[0] for x in l ]

plt.bar( lefts, heights, width=0.8 )
plt.subplots_adjust(bottom=0.25)
plt.xticks( [x+0.4 for x in lefts], ['{} P: {:.3f}'.format( x[1], epitope_zp2s[x[1]][1] ) for x in l], rotation='vertical', fontsize=9 )
plt.suptitle("Z-scores for mouse NNdistance rank score distributions\nHigher Z means more tendency for 'head-like' mice and 'tail-like' mice")
pngfile = '{}_mouse_nbrdist_rank_score_heterogeneity.png'.format( outfile_prefix )
print 'making:',pngfile
plt.savefig( pngfile )
util.readme( pngfile, """Phil write something here""")






