from basic import *
import numpy as np
from sklearn.decomposition import KernelPCA
import html_colors
import parse_tsv
import util
from operator import add

with Parser(locals()) as p:
    p.str('organism').required()
    p.str('clones_file').required()
    p.str('pngfile_prefix')
    p.int('max_labels').default(5)
    p.float('distance_scale_factor').default(0.01)
    p.float('paper_figs_dpi').default(100.)
    p.float('Dmax')
    p.float('minval')
    p.float('maxval')
    p.flag('show')
    p.flag('paper_figs')
    p.flag('vertical') ## each epitope is a vertical stack of scatter plots
    p.flag('showmotifs')
    p.multiword('epitopes').cast(lambda x: x.split())

if pngfile_prefix is None:
    pngfile_prefix = clones_file[:-4]

import matplotlib
matplotlib.rcParams['mathtext.default'] = 'regular'
if not show: matplotlib.use('Agg')
import matplotlib.pyplot as plt

greek_ab = {'a':r'$\alpha$', 'b':r'$\beta$'}

all_tcrs = parse_tsv.parse_tsv_file( clones_file, ['epitope'], ['clone_id'], True )

if showmotifs:
    motifs_file = clones_file[:-4]+'_motifs.tsv'
    if not exists( motifs_file ):
        showmotifs = False
    else:
        assert exists(motifs_file)

        all_motifs = parse_tsv.parse_tsv_file( motifs_file, ['epitope','chain'],
                                               ['id','showmotif', 'chi_squared',
                                                'matches_with_nbrs', 'matches_with_nbrs_consensus', 'expected_fraction',
                                                'cluster_number', 'is_cluster_center','cluster_consensus' ] )

if epitopes is None:
    epitopes = all_tcrs.keys()[:]
    epitopes.sort()

## first load all the distance matrices so we can use a consistent Dmax scaling (make landscapes comparable)
all_Ds = []
for epitope in epitopes:

    tcrs = all_tcrs[epitope]
    tcr_infos = [x[1] for x in all_tcrs[epitope]]
    num_tcrs = len(tcr_infos)

    ## read distances
    distfile = '{}_AB_{}.dist'.format(clones_file[:-4],epitope)

    all_dists = []
    for line in open( distfile,'r'):
        l = line.split()
        clone_id = l[0]
        index = len(all_dists)
        assert tcr_infos[ index ]['clone_id'] == clone_id
        dists = [ distance_scale_factor*float(x) for x in l[1:] ]
        assert len(dists) == num_tcrs
        all_dists.append( dists )

    D = np.matrix(all_dists)
    print epitope, D.shape, 'D.max()=',D.max()
    all_Ds.append( D )

max_Dmax = max( ( D.max() for D in all_Ds ) )
all_Dmax = [ D.max() for D in all_Ds ]
med_Dmax = get_median( all_Dmax )

print 'max_Dmax:',max_Dmax,'med_Dmax:',med_Dmax,'cmdline_Dmax:',Dmax

if Dmax is None:
    Dmax = med_Dmax ## Note: using median value of D.max()

## now do kpca
all_xys = []
jc_all_xys = {}
for epitope,D in zip(epitopes,all_Ds):
    old_Dmax = D.max()
    D = np.minimum( D, np.full( D.shape, Dmax ) )
    print 'true_Dmax:',old_Dmax,'using_Dmax:',Dmax,epitope

    pca = KernelPCA(kernel='precomputed')
    gram = 1 - ( D / Dmax )
    xy = pca.fit_transform(gram)

    xs = xy[:,0]
    ys = xy[:,1]

    all_xys.append( xy )
    jc_all_xys[epitope] = xy
all_vals = reduce( add, [list(xy[:,0]) for xy in all_xys ] ) + \
           reduce( add, [list(xy[:,1]) for xy in all_xys ] )


if minval is None:
    minval, maxval = min(all_vals), max(all_vals)
else:
    print 'minval {} maxval {}'.format(min(all_vals),max(all_vals))
    assert minval <= min(all_vals)
    assert maxval >= max(all_vals)

print 'minval {} maxval {}'.format(minval,maxval)


span = maxval-minval
minval -= 0.03*span
maxval += 0.03*span
span = maxval-minval

ticks = [ minval + x*span/10 for x in range(1,10) ]
#yticks = [ minval + x*span/10 for x in range(1,10) ]

if vertical:
    assert not showmotifs # still have to mod plotno calc below
    ncols = len(epitopes)
    nrows = 4 + 2*showmotifs
else:
    nrows = len(epitopes)
    ncols = 4 + 2*showmotifs

#plotno = 0


plot_width_inches = 2.0
left_margin_inches = 0.5
right_margin_inches = 0.5
top_margin_inches = 0.5
bottom_margin_inches = 0.5
if vertical:
    horizontal_spacer = 0.5
    vertical_spacer = 0.04
else:
    if paper_figs:
        vertical_spacer = 0.04
        horizontal_spacer = 0.04
        left_margin_inches = 0.05
        right_margin_inches = 0.05
        top_margin_inches = 0.05
        bottom_margin_inches = 0.05
    else:
        vertical_spacer = 0.5
        horizontal_spacer = 0.04

fig_height = top_margin_inches + bottom_margin_inches + nrows * plot_width_inches + (nrows-1)*vertical_spacer
fig_width = left_margin_inches + right_margin_inches + ncols * plot_width_inches + (ncols-1)*horizontal_spacer

bottom_margin = bottom_margin_inches / fig_height
top_margin = (fig_height-top_margin_inches)/fig_height

left_margin = left_margin_inches / fig_width
right_margin = ( fig_width-right_margin_inches)/fig_width

hspace = vertical_spacer*(nrows-1) / fig_height
wspace = horizontal_spacer*(ncols-1) / fig_width

plt.figure(1,figsize=(fig_width,fig_height))




def setup_gridl( xy, minval, maxval, step, nx, ny, max_labels ):

    grid_counts = {}
    for i in range(num_tcrs):
        x,y = xy[i][0], xy[i][1]
        xg = int( floor( (x-minval)/step ) )
        yg = int( floor( (y-minval)/step ) )
        g= (xg,yg)
        grid_counts[g] = grid_counts.get(g,0)+1

    gridl = []
    max_count=0

    while len(gridl)<max_labels:
        gridl = []
        for xg in range(ngrid-nx+1):
            for yg in range(ngrid-ny+1):
                count=0
                for xxg in range(xg,xg+nx):
                    for yyg in range(yg,yg+ny):
                        count += grid_counts.get( (xxg,yyg),0 )
                if count<=max_count:
                    #print 'safe:',count,epitope,xg,yg
                    gridl.append( ( count, (xg,yg) ) )

        max_count += 1

    return gridl


def add_new_label_and_update_gridl( xy, minval, maxval, step, nx, ny, label_text, label_color, label_indices, gridl,
                                    fontfamily = None ):
    if not gridl:
        return ## can't do anything

    ## find the best grid for this guy
    min_avgdist = 1e6
    for ii_gridl, (count,(xg,yg)) in enumerate(gridl):
        x = minval + xg*step + nx*step*0.5
        y = minval + yg*step + ny*step*0.5
        avgdist = 0.
        for i in label_indices:
            dist = sqrt( ( xy[i][0]-x)**2 + (xy[i][1]-y)**2 )
            avgdist += dist
        avgdist /= len(label_indices)
        if avgdist<min_avgdist:
            min_avgdist = avgdist
            best_ii_gridl = ii_gridl

    ## now draw a label at this point
    count,(xg,yg) = gridl[ best_ii_gridl ]

    ## remove gridl entries that would overlap with this new label
    for ii in range(len(gridl)-1,-1,-1):
        ## should we delete this guy?
        count2,(xg2,yg2) = gridl[ii]
        if abs(xg2-xg)<nx and abs(yg2-yg)<ny:
            del gridl[ii]

    x0 = minval + xg*step
    y0 = minval + yg*step
    if xg==0: x0 += step*0.3 ## scoot over a little bit
    if fontfamily:
        fontdict = {'family': fontfamily, 'size': 6 }
        plt.text( x0, y0, label_text, color=label_color, va = 'bottom', ha='left', fontdict=fontdict )
    else:
        plt.text( x0, y0, label_text, color=label_color, va = 'bottom', ha='left', fontsize=6 )

jcmaxlen = 0

kPCAset = set()
for ii_epitope, epitope in enumerate( epitopes ):
    plt.figure(1,figsize=(fig_width,fig_height))
    xy = all_xys[ ii_epitope ]

    tcrs = all_tcrs[epitope]
    tcr_infos = [x[1] for x in all_tcrs[epitope]]
    num_tcrs = len(tcr_infos)

    util.assign_label_reps_and_colors_based_on_most_common_genes_in_repertoire( tcr_infos, organism )
    ## now stored as va_label_rep, jb_label_rep

    for jj_reptype, reptype in enumerate( ['va','ja','vb','jb'] ):
        repcounts = {}
        tcr_reps = []
        tcr_colors = []
        rep_colors = {}
        for l in tcr_infos:
            rep = l[reptype+'_label_rep']
            color = l[reptype+'_label_rep_color']
            repcounts[rep] = repcounts.get(rep,0)+1
            tcr_reps.append( rep )
            tcr_colors.append( color )
            rep_colors[rep] = color
        repl = [ (y,x) for x,y in repcounts.iteritems()]
        repl.sort()
        repl.reverse()

        #plotno += 1
        if vertical:
            plotno = jj_reptype*ncols + ii_epitope+1
        else:
            plotno = ii_epitope*ncols + jj_reptype+1
        axes = plt.subplot(nrows,ncols,plotno)
        #axes.set_aspect('equal')
        plt.scatter( xy[:,0], xy[:,1], s=10, c=tcr_colors, edgecolors='none' )
        plt.xticks(ticks,[])
        plt.yticks(ticks,[])


        ## can we figure out where to put labels?
        ## look for rows of grid points that are empty
        ##
        ngrid = 20
        step = (maxval-minval)/ngrid
        nx = 4
        ny = 1

        gridl = setup_gridl( xy, minval, maxval, step, nx, ny, max_labels )

        ## now add labels for the top few reps
        mincount_for_labels = float(num_tcrs)/25 ## 4%
        for count,rep in repl[:max_labels]:
            if count<mincount_for_labels:break
            if not gridl:break

            rep_indices = [ x for x in range(num_tcrs) if tcr_reps[x] == rep ]

            add_new_label_and_update_gridl( xy, minval, maxval, step, nx, ny, rep, rep_colors[rep], rep_indices, gridl )

        plt.xlim((minval,maxval))
        plt.ylim((minval,maxval))

        if paper_figs:
            textpad = 0.05* (maxval-minval)
            plt.text( minval +textpad, maxval-textpad, epitope, ha='left',va='top', fontsize=14 )
            plt.text( maxval -textpad, maxval-textpad, '{}{} colors'.format( reptype[0].upper(), greek_ab[ reptype[1]]),
                      ha = 'right', va='top', fontsize=11 )
        else:
            plt.title('{}  {}{}'.format(epitope, reptype[0].upper(), greek_ab[ reptype[1] ] ) )


    if showmotifs:
        if epitope not in all_motifs: continue
        for jj_motifchain, motifchain in enumerate( 'AB' ):
            if motifchain not in all_motifs[epitope]: continue

            motifids = [] ## in decreasing order of chi-squared
            motifid2info = {}
            motif_cluster_centers = {} ## map from cluster_number to motifid
            motif_cluster_members = {} ## map from cluster_number to motifids

            for ( motifid, showmotif, chi_squared, matches_with_nbrs, matches_with_nbrs_consensus, expected_fraction,
                  cluster_number, is_cluster_center,cluster_consensus ) in all_motifs[epitope][motifchain]:
                #if float(expected_fraction)>max_expected_fraction: continue

                motifids.append( motifid )

                cnum = int( cluster_number )
                info = { 'cluster_number':cnum,
                         'matches_with_nbrs':matches_with_nbrs,
                         'showmotif':showmotif,
                         'consensus':matches_with_nbrs_consensus,
                         'cluster_consensus':cluster_consensus,
                         'chi_squared':float(chi_squared) }


                motifid2info[ motifid ] = info
                if cnum not in motif_cluster_members:
                    motif_cluster_members[cnum] = []
                motif_cluster_members[cnum].append( motifid )
                if int(is_cluster_center):
                    motif_cluster_centers[ cnum ] = motifid


            #print motif_cluster_centers
            assert max( motif_cluster_centers.keys() ) == len(motif_cluster_centers.keys())-1 ## 0-indexed cluster #s

            l_clusters = [ ( max( ( motifid2info[m]['chi_squared'] for m in members ) ), cnum )
                  for cnum,members in motif_cluster_members.iteritems() ]
            l_clusters.sort()
            l_clusters.reverse() ## decreasing order of top chi_squared for each cluster
            #print l_clusters

            motif_cluster_colors = dict( zip( [ x[1] for x in l_clusters ],
                                              html_colors.get_rank_colors_no_lights( len(motif_cluster_centers) ) ) )

            plotno = ii_epitope*ncols + 4 + jj_motifchain + 1
            assert not vertical
            axes = plt.subplot(nrows,ncols,plotno)
            #axes.set_aspect('equal')
            ## first show them all light grey
            xs = xy[:,0]
            ys = xy[:,1]

            light_gray = '#D3D3D3'
            tcr_colors = [light_gray]*num_tcrs


            ngrid = 20
            step = (maxval-minval)/ngrid
            nx = 6 ## little bit bigger
            ny = 1

            gridl = setup_gridl( xy, minval, maxval, step, nx, ny, max_labels )

            motif_cluster_indices = {}
            for cnum in motif_cluster_colors:
                motif_cluster_indices[cnum] = set()


            all_matched_ids = set()
            for motifid in reversed( motifids ):
                info = motifid2info[ motifid ]
                cluster_number = info['cluster_number'] #0-indexed
                color = motif_cluster_colors[ cluster_number ]
                #showmotif = motifl[0]
                matched_ids = frozenset( info['matches_with_nbrs'].split(',') ) ## including nbrs
                all_matched_ids |= matched_ids
                indices = frozenset( [ x for x in range(num_tcrs) if tcr_infos[x]['clone_id'] in matched_ids ] )
                motif_cluster_indices[cluster_number] |= indices
                for ind in indices:
                    tcr_colors[ind] = color

            for (chi_squared,cnum) in l_clusters[:max_labels]:
                members = motif_cluster_members[ cnum ]
                indices = motif_cluster_indices[ cnum ]
                color = motif_cluster_colors[cnum]
                if len(indices)<mincount_for_labels: continue
                num_visible = tcr_colors.count(color)
                if num_visible<mincount_for_labels: continue
                top_chi_squared, top_motifid = max( ( (motifid2info[m]['chi_squared'],m) for m in members ) )
                assert abs(top_chi_squared-chi_squared)<1e-3

                #label = motifid2info[top_motifid]['consensus']
                #label = motifid2info[ motif_cluster_centers[ cnum ] ]['consensus']
                label = motifid2info[ motif_cluster_centers[ cnum ] ]['cluster_consensus']

                if gridl:
                    add_new_label_and_update_gridl( xy, minval, maxval, step, nx, ny, label, color,
                                                    indices, gridl, fontfamily='monospace' )

            assert len(tcr_colors) == len(xs)

            ## now do all the plotting
            plt.scatter( [ xs[i] for i in range(num_tcrs) if tcr_colors[i] == light_gray ],
                         [ ys[i] for i in range(num_tcrs) if tcr_colors[i] == light_gray ],
                         s=10, c=light_gray, edgecolors='none' )
            for motifid in reversed( motifids ):
                info = motifid2info[ motifid ]
                cluster_number = info['cluster_number'] #0-indexed
                color = motif_cluster_colors[ cluster_number ]
                plt.scatter( [ xs[i] for i in range(num_tcrs) if tcr_colors[i] == color ],
                             [ ys[i] for i in range(num_tcrs) if tcr_colors[i] == color ],
                             s=10, c=color, edgecolors='none' )

            for ji in range(num_tcrs): #JCC--adding ability to output kPC info
                templi = []
                for jx in tcr_infos[ji]:
                    templi.append(tcr_infos[ji][jx].strip())
                tempstr = tcrs[ji][0] + "\t" +  tcr_infos[ji]["epitope"] + "\t" +  str(xs[ji]) + "\t" +  str(ys[ji]) + "\t" +  tcr_colors[ji] + "\t" + "\t".join(str(x) for x in jc_all_xys[epitope][ji])
                jcmaxlen = max(jcmaxlen, len(tempstr.split("\t")))
                kPCAset.add(tempstr) #--JCC
            print '{} {} matched: {} {:.6f}'.format(epitope,motifchain,len(all_matched_ids),
                                                    float(len(all_matched_ids))/num_tcrs)
            plt.xticks(ticks,[])
            plt.yticks(ticks,[])

            plt.xlim((minval,maxval))
            plt.ylim((minval,maxval))

            plt.title('{} {}-motifs'.format(epitope,motifchain.upper()))


    if paper_figs:
        ## make a special secret figure just for this epitope
        margin=0.1
        sep=0.07
        pwidth=2.
        if vertical:
            height=2*margin+3*sep+4*pwidth
            width=2*margin+pwidth
        else:
            width=2*margin+3*sep+4*pwidth
            height=2*margin+pwidth
        plt.figure(2,figsize=(width,height))
        plt.clf()
        for jj_reptype, reptype in enumerate( ['va','ja','vb','jb'] ):
            repcounts = {}
            tcr_reps = []
            tcr_colors = []
            rep_colors = {}
            for l in tcr_infos:
                rep = l[reptype+'_label_rep']
                color = l[reptype+'_label_rep_color']
                repcounts[rep] = repcounts.get(rep,0)+1
                tcr_reps.append( rep )
                tcr_colors.append( color )
                rep_colors[rep] = color
            repl = [ (y,x) for x,y in repcounts.iteritems()]
            repl.sort()
            repl.reverse()

            nr,nc = (4,1) if vertical else (1,4)
            plt.subplot(nr,nc,1+jj_reptype)
            #axes.set_aspect('equal')
            textpad = 0.03*(maxval-minval)
            plt.scatter( xy[:,0], xy[:,1], s=10, c=tcr_colors, edgecolors='none' )
            if reptype=='va': ## epitope label
                plt.text( maxval-1.5*textpad, maxval-1.5*textpad, epitope,
                          va='top', ha='right', fontsize=14 )
            plt.text( maxval-textpad, minval+textpad,
                      '{}{} colors'.format( reptype[0].upper(), greek_ab[ reptype[1] ] ),
                      va='bottom', ha='right', fontsize=8 )
            # plt.text( (maxval+minval)/2., minval+textpad,
            #           '{}{} colors'.format( reptype[0].upper(), greek_ab[ reptype[1] ] ),
            #           va='bottom', ha='center', fontsize=8 )
            plt.xticks(ticks,[])
            plt.yticks(ticks,[])


            ## can we figure out where to put labels?
            ## look for rows of grid points that are empty
            ##
            ngrid = 20
            step = (maxval-minval)/ngrid
            nx = 4
            ny = 1

            gridl = setup_gridl( xy, minval, maxval, step, nx, ny, max_labels )

            ## now add labels for the top few reps
            mincount_for_labels = float(num_tcrs)/25 ## 4%
            for count,rep in repl[:max_labels]:
                if count<mincount_for_labels:break
                if not gridl:break

                rep_indices = [ x for x in range(num_tcrs) if tcr_reps[x] == rep ]

                add_new_label_and_update_gridl( xy, minval, maxval, step, nx, ny, rep, rep_colors[rep], rep_indices, gridl )

            plt.xlim((minval,maxval))
            plt.ylim((minval,maxval))

            #plt.title('{}  {}'.format(epitope,reptype.upper()))
        plt.subplots_adjust( hspace=(3*sep)/height, wspace=sep/width, left=margin/width, right=(width-margin)/width,
                             bottom=margin/height, top=(height-margin)/height )
        pngfile = '{}_{}_kpca.png'.format(pngfile_prefix,epitope)
        print 'making:',pngfile
        plt.savefig(pngfile,dpi=paper_figs_dpi)
        # ## make a special secret figure just for this epitope
        # margin=0.1
        # sep=0.1
        # pwidth=2.
        # width=2*margin+sep+2*pwidth
        # plt.figure(2,figsize=(width,width))
        # plt.clf()
        # for jj_reptype, reptype in enumerate( ['va','ja','vb','jb'] ):
        #     repcounts = {}
        #     tcr_reps = []
        #     tcr_colors = []
        #     rep_colors = {}
        #     for l in tcr_infos:
        #         rep = l[reptype+'_label_rep']
        #         color = l[reptype+'_label_rep_color']
        #         repcounts[rep] = repcounts.get(rep,0)+1
        #         tcr_reps.append( rep )
        #         tcr_colors.append( color )
        #         rep_colors[rep] = color
        #     repl = [ (y,x) for x,y in repcounts.iteritems()]
        #     repl.sort()
        #     repl.reverse()

        #     plt.subplot(2,2,1+jj_reptype)
        #     #axes.set_aspect('equal')
        #     plt.scatter( xy[:,0], xy[:,1], s=10, c=tcr_colors, edgecolors='none' )
        #     plt.xticks(ticks,[])
        #     plt.yticks(ticks,[])


        #     ## can we figure out where to put labels?
        #     ## look for rows of grid points that are empty
        #     ##
        #     ngrid = 20
        #     step = (maxval-minval)/ngrid
        #     nx = 4
        #     ny = 1

        #     gridl = setup_gridl( xy, minval, maxval, step, nx, ny, max_labels )

        #     ## now add labels for the top few reps
        #     mincount_for_labels = float(num_tcrs)/25 ## 4%
        #     for count,rep in repl[:max_labels]:
        #         if count<mincount_for_labels:break
        #         if not gridl:break

        #         rep_indices = [ x for x in range(num_tcrs) if tcr_reps[x] == rep ]

        #         add_new_label_and_update_gridl( xy, minval, maxval, step, nx, ny, rep, rep_colors[rep], rep_indices, gridl )

        #     plt.xlim((minval,maxval))
        #     plt.ylim((minval,maxval))

        #     #plt.title('{}  {}'.format(epitope,reptype.upper()))
        # space, lower, upper = sep/width, margin/width, (width-margin)/width
        # plt.subplots_adjust( hspace=space, wspace=space, left=lower, right=upper, bottom=lower, top=upper )
        # pngfile = '{}_{}_kpca.png'.format(pngfile_prefix,epitope)
        # print 'making:',pngfile
        # plt.savefig(pngfile,dpi=paper_figs_dpi)

plt.figure(1,figsize=(fig_width,fig_height))

plt.subplots_adjust( hspace=hspace, wspace = wspace, left=left_margin, right = right_margin,
                     bottom= bottom_margin, top=top_margin )
    #plt.suptitle('epitope={}   2D kernal-PCA projection'.format(epitope),size='large')

pngfile = pngfile_prefix + '_kpca.png'
print 'making',pngfile
plt.savefig(pngfile,dpi=paper_figs_dpi)

util.readme( pngfile, """
Kernel Principle Components Analysis (kPCA) 2D projection plots for the repertoires. Each row is a repertoire,
and each point in the plots corresponds to a single TCR clone, with the points arranged so as to keep nearby
TCRs (as measured by TCRdist) nearby in 2 dimensions. The four different panels are the same 2D projection
colored by gene usage for the four different segments (left to right: Va,Ja,Vb,Jb)""")

if show:
    plt.show()

print "kPCA results:"
print "clone.id" + "\t" + "epitope" + "\t" + "XS" + "\t" + "YS" + "\t" + "Color" + "\t" + "\t".join([("kPC" + str(i)) for i in range(jcmaxlen-5)])
for l in kPCAset:
    print l
