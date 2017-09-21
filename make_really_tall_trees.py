from basic import *
import score_trees_devel
import svg_basic
import numpy as np
import util
import html_colors
import scipy.stats
import copy
import random
from operator import add
import tcr_sampler
#from mannwhitneyu import mannwhitneyu as mannwhitneyu_exact #too slow


with Parser(locals()) as p:
    p.str('clones_file').required()
    p.str('organism').required()
    p.float('distance_scale_factor').default(0.01)
    p.float('extra_color_schemes_none_score').shorthand('none_score')
    p.flag('dont_trim_labels')
    p.flag('constant_seed')
    p.multiword('only_epitopes').cast(lambda x:x.split())
    p.multiword('extra_color_schemes').shorthand('colors').cast(lambda x:x.split())

if constant_seed: random.seed(1)
fake_chains = util.detect_fake_chains( clones_file )

probs_cs            = 'probs'
sharing_cs          = 'sharing'
cross_reactivity_cs = 'cross_reactivity'
clonality_cs        = 'clonality'
min_other_nbrdist_cs   = 'min_other_nbrdist'

color_scheme_explanations = {}


color_schemes = [ probs_cs, sharing_cs, cross_reactivity_cs, clonality_cs, min_other_nbrdist_cs ]

extra_color_scheme_prefix = 'FC_'

if extra_color_schemes:
    for tsvtag in extra_color_schemes:
        color_schemes.append( extra_color_scheme_prefix + tsvtag )

gap_character = '-' ## different from some other places
min_cluster_size = 1

#cluster_radius = {'AB':1.0, 'A':0.5, 'B':0.5}
distance_threshold_25_scaled = distance_scale_factor * pipeline_params[ 'distance_threshold_25' ]

cluster_radius = {
    'A' : distance_threshold_25_scaled*2, ## was 0.5
    'B' : distance_threshold_25_scaled*2, ## was 0.5
    'AB': distance_threshold_25_scaled*4  ## was 1.0
}


tree_width = 750
ymargin = 30 ## right now we dont want top text to get cut off
xmargin = 10

text_column_separation = 4
labels_tree_separation = 10

#labels_width = 700##approx
#tree_x0 = labels_width + xmargin + xpad

#total_svg_width = tree_x0 + tree_width + xmargin

branch_width_fraction = 0.1

log10_of_zero = -100
def get_safe_log10(f ):
    if f==0: return log10_of_zero
    else: return math.log10(f)



def pad_to_middle( s, num ): ## with spaces
    if len(s) >= num:return s
    extra = num-len(s)
    before = extra/2
    after = extra-before
    return ' '*before + s + ' '*after

def get_primary_number( gene_name ):
    tmp = gene_name[:]
    while tmp and not tmp[0].isdigit():
        tmp = tmp[1:]
    if not tmp:
        ## for example, if gene_name=='TRGJP'
        return 0
    assert tmp[0].isdigit()
    if tmp.isdigit():
        return int(tmp)
    else:
        tmp2 = ''
        while tmp[0].isdigit():
            tmp2 += tmp[0]
            tmp = tmp[1:]
        return int(tmp2)



## based these next two functions from ../distances.new.py
# def tree_leaves( tree ):
#     if tree[0] == tree[1]:
#         return [ tree[0] ]
#     else:
#         return tree_leaves( tree[0] ) + tree_leaves( tree[1] )

def tree_splits_ttest( edge_pvals, tree, other_leaves, leaf_scores, leaf_names, info='', pvalue_threshold = 1e-2 ):

    if tree[0] == tree[1]:
        return
    for i in range(2):
        a_leaves = score_trees_devel.Node_members( tree[i] )
        b_leaves = score_trees_devel.Node_members( tree[(i+1)%2] ) + other_leaves
        a_scores = reduce( add, [ leaf_scores[x] for x in a_leaves ] )
        b_scores = reduce( add, [ leaf_scores[x] for x in b_leaves ] )
        if a_scores and b_scores:
            amean = sum(a_scores)/len(a_scores)
            bmean = sum(b_scores)/len(b_scores)
            if amean != bmean:
                t, t_pvalue1 = scipy.stats.ttest_ind( a_scores, b_scores, equal_var = True  )
                t, t_pvalue2 = scipy.stats.ttest_ind( a_scores, b_scores, equal_var = False )
                u, u_pvalue = scipy.stats.mannwhitneyu( a_scores, b_scores )

                maxp = max( [t_pvalue1, t_pvalue2, u_pvalue] )
                if maxp < pvalue_threshold:

                    # if len(a_scores)<5 or len(b_scores)<5: ## try the exact calculation
                    #     print 'calc exact:',len(a_scores),len(b_scores)
                    #     u2, u_pvalue_exact = mannwhitneyu_exact( a_scores, b_scores )
                    #     print 'exact MWU:',u,u2,u_pvalue, u_pvalue_exact
                    print 'pvalue: {:.3e} {:.3e} {:.3e} {:.3e} {} a: {} {} {:.2f} b: {} {} {:.2f}'\
                        .format( maxp, t_pvalue1, t_pvalue2, u_pvalue, info,
                                 len(a_leaves), leaf_names[a_leaves[0]], amean,
                                 len(b_leaves), leaf_names[b_leaves[0]], bmean )
                    k = tuple( sorted( a_leaves ) )
                    assert k not in edge_pvals
                    symbol = '-' if amean < bmean else '+'
                    edge_pvals[ k ] = [ maxp, t_pvalue1, t_pvalue2, u_pvalue, symbol ]

        ## recurse
        tree_splits_ttest( edge_pvals, tree[i], b_leaves, leaf_scores, leaf_names, info )


def label_pval_edges( cmds, edge_pvals, subtree, plotting_info ):
    sizes, node_position, Transform, canvas_tree_w_factor, canvas_tree_min_rmsd = plotting_info

    if score_trees_devel.IsALeaf( subtree ):
        return
    else:
        big_rmsd = subtree[2]
        for ii in range(2):
            iitree = subtree[ii]
            little_rmsd = iitree[2]
            #assert little_rmsd<=big_rmsd
            assert little_rmsd <= big_rmsd+1e-3
            if little_rmsd > big_rmsd:
                print 'WHOAH:',little_rmsd,big_rmsd
            leaves = tuple( sorted( score_trees_devel.Node_members(iitree) ) )
            if leaves in edge_pvals:
                pvals = edge_pvals[leaves]
                symbol = pvals[-1]

                center= score_trees_devel.Center( iitree, node_position, sizes, use_sizes_as_weights=True )
                size = score_trees_devel.Size( iitree, sizes )

                line_width = max(1,int(floor(0.5+ size*canvas_tree_w_factor )))
                box_x0 = Transform(max(canvas_tree_min_rmsd,little_rmsd)) ; box_x1 = Transform(big_rmsd)
                sep = 3

                cmds.append( svg_basic.make_text( '{:.0E} {}'.format( pvals[0], symbol ),
                                                  ( box_x0+sep, center-line_width/2-sep),
                                                  10, font_family="Droid Sans Mono" ) )

            label_pval_edges( cmds, edge_pvals, iitree, plotting_info )


## little class to store some info
class TCR:
    def __init__( self, l ):
        self.subject = l['subject']
        self.epitope = l['epitope']
        self.cdr3a = l['cdr3a']
        self.cdr3b = l['cdr3b']
        if 'cdr3a_protseq_masked' not in l or 'cdr3b_protseq_masked' not in l:
            tcr_sampler.add_masked_CDR3_sequences_to_tcr_dict( organism, l )
        self.cdr3a_masked = l['cdr3a_protseq_masked']
        self.cdr3b_masked = l['cdr3b_protseq_masked']
        self.a_indels = l['a_indels']
        self.b_indels = l['b_indels']
        self.clone_id = l['clone_id']
        self.clone_size = int( l['clone_size'] )
        self.info = copy.deepcopy(l)

        #genes= []
        self.reps_for_sharing = []
        self.reps_for_counting = []
        for ab in 'ab':
            for vj in 'vj':
                #hits = l['{}{}_blast_hits'.format(vj,ab)]
                genes = set( l['{}{}_genes'.format(vj,ab)].split(';') )
                self.reps_for_sharing.append( util.reps_from_genes( genes, organism=organism, mm1=False ) )
                #self.reps_for_counting.append( util.countreps_from_genes( genes, organism=organism ) )
        return

## parse the clones_file ############################################################################################3

all_tcrs = {}
#all_color_scores = {}

infields = []
clones_file_with_nbrdists = '{}_nbrdists.tsv'.format(clones_file[:-4])
assert exists( clones_file_with_nbrdists )

Log('parsing {}'.format(clones_file_with_nbrdists))
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
    tcr = TCR( l )


    ## figure out the color scores, one for each chains
    ## color by probs
    tcr.color_scores = {}
    tcr.color_scores[probs_cs] = { 'A': get_safe_log10( float(l['a_protseq_prob']) ),
                                   'B': get_safe_log10( float(l['b_protseq_prob']) ),
                                   'AB': get_safe_log10( float(l['a_protseq_prob']) * float( l['b_protseq_prob'] ) ) }

    #clone_size = float(l['clone_size'])

    tcr.color_scores[clonality_cs] = { 'A': tcr.clone_size, 'B': tcr.clone_size, 'AB': tcr.clone_size }

    ## look at rank scores for other epitopes
    #suffix = '_rank25'
    suffix = '_wtd_nbrdist10rank'
    tcr.color_scores[ min_other_nbrdist_cs ] = {}
    for ab in ['A','B','AB']:
        other_ranks = []
        for tag,val in l.iteritems():
            suf = '_{}{}'.format(ab,suffix)
            if tag.endswith(suf):
                ep = tag[:-1*len(suf)]
                if ep != epitope:
                    other_ranks.append( ( float(val), ep ) )
        if other_ranks:
            tcr.color_scores[ min_other_nbrdist_cs ][ ab ] = min( other_ranks ) ## (val,other_ep)
        else:
            tcr.color_scores[ min_other_nbrdist_cs ][ ab ] = (0.0,'NA')

    ## extra color scores
    if extra_color_schemes:
        for tsvtag in extra_color_schemes:
            scheme = extra_color_scheme_prefix + tsvtag
            tcr.color_scores[scheme] = {}
            for ab in ['A','B','AB']:
                tcr.color_scores[scheme][ab] = float( l[tsvtag] )

    if epitope not in all_tcrs:
        all_tcrs[epitope] = []

    all_tcrs[ epitope ].append( tcr )



epitopes = sorted( all_tcrs.keys()[:] )


# color_score_range = {'A':None, 'B':None, 'AB':None }
color_score_range_probs_cs = {}
color_score_range_probs_cs[ 'A'] = (  -9.0,  -5.0 )
color_score_range_probs_cs[ 'B'] = (  -9.5,  -5.5 )
color_score_range_probs_cs['AB'] = ( -16.0, -12.0 )

def same_tcr( t1, t2, chains ): ## t1 = (subject,genes,reps,cdr3a,cdr3b,...)
    #t1_reps = t1[2] ## [va,ja,vb,jb]
    #t2_reps = t2[2]
    if 'A' in chains:
        if ( t1.reps_for_sharing[0].isdisjoint( t2.reps_for_sharing[0] ) or
             t1.reps_for_sharing[1].isdisjoint( t2.reps_for_sharing[1] ) or
             t1.cdr3a != t2.cdr3a ):
            return False
    if 'B' in chains:
        if ( t1.reps_for_sharing[2].isdisjoint( t2.reps_for_sharing[2] ) or
             t1.reps_for_sharing[3].isdisjoint( t2.reps_for_sharing[3] ) or
             t1.cdr3b != t2.cdr3b ):
            return False
    return True

for epitope in epitopes:
    if only_epitopes and epitope not in only_epitopes: continue

    tcrs = all_tcrs[epitope]
    infos = [x.info for x in tcrs]
    ## this fills xx_label_rep and xx_label_rep_color in each dict in the infos list
    util.assign_label_reps_and_colors_based_on_most_common_genes_in_repertoire( infos, organism )
    rep_colors = {}
    for tcr, info in zip(tcrs,infos):
        tcr.single_reps = []
        for ab in 'ab':
            for vj in 'vj':
                rep = info[ '{}{}_label_rep'.format(vj,ab) ]
                color = info[ '{}{}_label_rep_color'.format(vj,ab) ]
                tcr.single_reps.append( rep )
                rep_colors[ rep ] = color

    epitope_mice = list( set( [ x.subject for x in tcrs ] ) )
    epitope_mice.sort() ## will use these to label the tree


    ## let's compute sharing
    Log('compute sharing '+epitope)
    for tcr in tcrs:
        tcr.color_scores[ cross_reactivity_cs ] = {}
        tcr.color_scores[ sharing_cs ] = {}

        for ab in ['A','B','AB']:

            ## look for sharing with other epitopes or other mice (same epitope)
            counts_xr = {}
            mice = [[], [] ]
            for ep2 in epitopes:
                for tcr2 in all_tcrs[ ep2 ]:
                    if same_tcr( tcr, tcr2, ab ):
                        if epitope != ep2:
                            counts_xr[ ep2 ] = counts_xr.get( ep2,0)+1
                        mice[ epitope == ep2 ].append( tcr2.subject )

            if counts_xr:
                other_epitopes = counts_xr.keys()[:]
                other_epitopes.sort()
            else:
                other_epitopes = []
            tcr.color_scores[ cross_reactivity_cs ][ab] = ( 1+len(other_epitopes), other_epitopes )#(Nepitopes,other-eps)
            tcr.color_scores[ sharing_cs ][ab] = [ len(set(mice[1])), len(set(mice[0]+mice[1])) ]


    for ab in ['A','B','AB']:
        if ab in fake_chains: continue

        radius = cluster_radius[ab]
        distfile = '{}_{}_{}.dist'.format( clones_file[:-4], ab, epitope )
        assert exists(distfile)
        Log('reading '+distfile)

        N=0
        all_nbrs = []
        all_dists = []
        for line in open( distfile,'r'):
            l = line.split()
            clone_id = l[0]
            assert tcrs[ len(all_nbrs) ].clone_id == clone_id
            dists = [ distance_scale_factor*float(x) for x in l[1:] ]
            if not N:
                N = len(dists)
            else:
                assert N == len(dists)

            nbrs = []
            for ii,d in enumerate(dists):
                if d <= radius:
                    nbrs.append( ii )
            all_dists.append( dists )
            all_nbrs.append( nbrs )

        deleted = [False]*N

        centers = []
        all_members = []

        Log('clustering {} tcrs'.format(len(tcrs)))
        while True:
            clusterno = len(centers)

            best_nbr_count =0
            for i in range(N):
                if deleted[i]: continue
                nbr_count = 0
                for nbr in all_nbrs[i]:
                    if not deleted[nbr]:
                        nbr_count+=1
                if nbr_count > best_nbr_count:
                    best_nbr_count = nbr_count
                    center = i

            if best_nbr_count < min_cluster_size:
                break

            centers.append( center )
            members = [center]

            deleted[center] = True
            for nbr in all_nbrs[center]:
                if not deleted[nbr]:
                    deleted[nbr] = True
                    members.append( nbr )

            assert len(members) == best_nbr_count
            all_members.append( frozenset(members) )

        num_clusters = len(centers)
        num_tcrs = len(tcrs)

        ## I think this will give a better ordering of the TCRs along the tree
        ## order from 1....N by going through the members of the clusters, largest to smallest
        old2new_index = {}
        new2old_index = {}
        last_index=-1
        for members in all_members:
            for member in members:
                last_index += 1
                old2new_index[ member ] = last_index
                new2old_index[ last_index ] = member
        assert len(old2new_index) == num_tcrs

        ## how much vertical space will we need?
        ##
        label_fontsize = 10
        tree_height = label_fontsize * len(tcrs)
        total_svg_height = tree_height + 2*ymargin


        max_clone_count=max( ( x.clone_size for x in tcrs ) )

        def clonality_fraction( clone_size ):
            global max_clone_count
            if max_clone_count==1: return 0.0
            exponent = 1.0/3
            mx = max_clone_count**exponent
            cs = clone_size**exponent
            return ( cs-1.0)/(mx-1.0)

        def clonality_color( clone_size ):
            return svg_basic.rgb_from_fraction( clonality_fraction( clone_size ) )

        #max_num_mice_same_epitope = max( ( x.color_scores[ sharing_cs ][ab][0] for x in tcrs ) )
        max_num_mice_all_epitopes = max( ( x.color_scores[ sharing_cs ][ab][1] for x in tcrs ) )
        max_num_epitopes = max( ( x.color_scores[ cross_reactivity_cs ][ab][0] for x in tcrs ) )

        def num_mice_all_epitopes_color( num_mice_all_epitopes ):
            if max_num_mice_all_epitopes==1:
                return svg_basic.rgb_from_fraction( 0.0 )
            else:
                return svg_basic.rgb_from_fraction( float(num_mice_all_epitopes-1)/(max_num_mice_all_epitopes-1))

        def num_epitopes_color( num_epitopes ):
            if max_num_epitopes == 1:
                return svg_basic.rgb_from_fraction( 0.0 )
            else:
                return svg_basic.rgb_from_fraction( float(num_epitopes-1)/(max_num_epitopes-1) )

        tree = None

        for color_scheme in color_schemes:

            ## let's fill out an array of color scores
            my_color_scores = [ x.color_scores[ color_scheme ][ab] for x in tcrs ]

            my_color_scores_labels = ['']*num_tcrs ## will go into the label text, at the end

            color_score_range = None

            if color_scheme == probs_cs:
                ## adjust scores of low-prob guys
                min_good_score = min( [x for x in my_color_scores if x!=log10_of_zero] )
                my_color_scores_floats = [ max(min_good_score,x) for x in my_color_scores ]

                color_score_range = color_score_range_probs_cs[ab]
            elif color_scheme == clonality_cs:
                my_color_scores_floats = [ clonality_fraction(float(x)) for x in my_color_scores ]
                my_color_scores_labels = [ '{:2d}'.format(x) for x in my_color_scores ]

            elif color_scheme == cross_reactivity_cs:
                my_color_scores_floats = [float(x[0]) for x in my_color_scores ]
                my_color_scores_labels = [ ' '.join(x[1]) for x in my_color_scores ]

            elif color_scheme == sharing_cs:
                my_color_scores_floats = [ float(x[0]) for x in my_color_scores ] # N-mice this epitope
                my_color_scores_labels = [ '{}'.format(x[0]) for x in my_color_scores ]

            elif color_scheme == min_other_nbrdist_cs:
                my_color_scores_floats = [-1*float(x[0]) for x in my_color_scores ]
                my_color_scores_labels = [ '{:3d} {}'.format( int(floor(0.5+x[0])), x[1] ) for x in my_color_scores ]

                color_score_range = (-100,0)

            elif color_scheme.startswith(extra_color_scheme_prefix):
                my_color_scores_floats = [ None if x == extra_color_schemes_none_score else x for x in my_color_scores ]
                my_color_scores_labels = ['NA' if x == extra_color_schemes_none_score else '{:.1f}'.format(x)
                                          for x in my_color_scores ]

            if my_color_scores_floats.count(None) == len(my_color_scores_floats):
                print 'skipping empty color scheme:',color_scheme,epitope
                continue ## no scores for this guy


            mn_score_color = min( [ x for x in my_color_scores_floats if x!=None ] )
            mx_score_color = max( [ x for x in my_color_scores_floats if x!=None ] )

            def get_tcr_score_color( index ):
                score = my_color_scores_floats[index]
                if score==None:return 'black'
                if color_score_range:
                    mn,mx = color_score_range
                else:
                    mn,mx = mn_score_color, mx_score_color
                if mx==mn: mx=mn+1
                return svg_basic.rgb_from_fraction( max(0.0,min(1.0, ( score-mn )/(mx-mn ) ) ) )

            ## now get some info together for plotting
            all_center_dists = {}
            all_scores = []
            sizes = []
            names = []
            infos = [] ## for the color score correlations

            for new_index in range(num_tcrs):
                old_index = new2old_index[ new_index ]
                names.append( '' )
                sizes.append( 1 )
                infos.append( '{} {}'.format( tcrs[old_index].cdr3a[3:-2], tcrs[old_index].cdr3b[3:-2] ))
                color_score = my_color_scores_floats[ old_index ]
                if color_score==None:
                    all_scores.append( [] )
                else:
                    all_scores.append( [color_score] )
                for other_new_index in range(num_tcrs):
                    other_old_index = new2old_index[ other_new_index ]
                    dist = all_dists[ old_index ][ other_old_index ]
                    all_center_dists[ (new_index,other_new_index) ] = dist
                    all_center_dists[ (other_new_index,new_index) ] = dist

            percentile = -1

            Log('Make_tree')
            if not tree:
                tree = score_trees_devel.Make_tree_new( all_center_dists, len(names),
                                                        score_trees_devel.Update_distance_matrix_AL,
                                                        all_scores, score_trees_devel.CallAverageScore(percentile) )
            else:
                tree = score_trees_devel.Copy_tree_update_scores( tree, all_scores,
                                                                  score_trees_devel.CallAverageScore(percentile))

            ## look for branches with high/low scores
            edge_pvals = {}
            tree_splits_ttest( edge_pvals, tree, [], all_scores, infos, epitope+'_'+ab+"_"+color_scheme )


            ## the x-values dont matter here
            ## but the y-values do
            tree_p0 = [10, ymargin ]
            tree_p1 = [1000, tree_height+ymargin ]


            ## node_position tells us where the different clusters are located, vertically
            ##
            Log('Canvas_tree 1st time')
            tmp_plotter = svg_basic.SVG_tree_plotter()
            node_position,Transform,canvas_tree_min_rmsd, canvas_tree_w_factor \
                = score_trees_devel.Canvas_tree( tree, names, sizes, tree_p0, tree_p1, branch_width_fraction,
                                                 tmp_plotter, label_internal_nodes = False,
                                                 score_range_for_coloring = color_score_range )

            cmds = []
            cmds.append( svg_basic.make_text( '{} {} #tcrs: {}  {} colors'.format(epitope,ab,len(tcrs),color_scheme),
                                              ( 10, ymargin ),
                                              30, font_family="Droid Sans Mono" ) ) ## label the epitope

            ## now let's add some text for each tcr
            num_columns = 6 + 2*len(ab)
            text_columns = []
            for i in range( num_columns ):
                text_columns.append( [] )
            header = ['']*num_columns

            for old_index, tcr in enumerate( tcrs ):
                ## 1. have some text on the left with the cdr3s, masked cdr3s and indles
                ## 2. clonality text
                ## 3. sharing text -- yeah
                ## 4. gene segments text, actually each gene segment separate,maybe since they are different colors
                ## 5. color label
                ##
                ## which mouse/subject
                icol = 0
                text_columns[ icol ].append( ( '{:2d}{:2s}'.format(epitope_mice.index(tcr.subject)+1,tcr.subject[:2]),
                                               'black' ) )
                header[icol] = 'M#'
                icol += 1


                text0 = ''
                if 'A' in ab:
                    text0 += '{} {} {}'\
                        .format( pad_to_middle( tcr.cdr3a if dont_trim_labels else tcr.cdr3a[3:-2], 15 ),
                                 pad_to_middle( tcr.cdr3a_masked if dont_trim_labels else tcr.cdr3a_masked[3:-2], 15 ),
                                 pad_to_middle( tcr.a_indels, 6 ) )
                if 'B' in ab:
                    if text0: text0 += ' '
                    text0 += '{} {} {}'\
                        .format( pad_to_middle( tcr.cdr3b if dont_trim_labels else tcr.cdr3b[3:-2], 15 ),
                                 pad_to_middle( tcr.cdr3b_masked if dont_trim_labels else tcr.cdr3b_masked[3:-2], 15 ),
                                 pad_to_middle( tcr.b_indels, 6 ) )

                text_columns[ icol ].append( ( text0, 'black' ) )
                icol += 1

                ## clonality
                header[icol] = ' C'
                text_columns[ icol ].append( ( '{:2d}'.format(tcr.clone_size), clonality_color( tcr.clone_size ) ) )
                icol += 1

                ## sharing
                header[icol] = ' M '
                num_mice_this_epitope, num_mice_all_epitopes = tcr.color_scores[ sharing_cs ][ab]
                text_columns[ icol ].append( ( '{},{}'.format( num_mice_this_epitope, num_mice_all_epitopes ),
                                               num_mice_all_epitopes_color( num_mice_all_epitopes ) ) )
                icol += 1

                ## cross-reactivity:
                header[icol] = 'E'
                num_epitopes = tcr.color_scores[ cross_reactivity_cs ][ab][0] ## includes this one
                text_columns[ icol ].append( ( str(num_epitopes), num_epitopes_color( num_epitopes ) ) )
                icol += 1

                ## gene segments
                reps = tcr.single_reps
                if 'A' in ab:
                    assert reps[0].startswith('TR') and reps[1].startswith('TR')
                    text_columns[ icol ].append( ( '{}{:02d}'.format( reps[0][2:4], get_primary_number( reps[0] ) ),
                                                   rep_colors[reps[0]] ) )
                    icol += 1
                    text_columns[ icol ].append( ( '{}{:02d}'.format( reps[1][2:4], get_primary_number( reps[1] ) ),
                                                   rep_colors[reps[1]] ) )
                    icol += 1
                    # text_columns[ icol ].append( ( 'AV{:02d}'.format( get_primary_number( reps[0] ) ), rep_colors[reps[0]] ) ) ; icol += 1
                    # text_columns[ icol ].append( ( 'AJ{:02d}'.format( get_primary_number( reps[1] ) ), rep_colors[reps[1]] ) ) ; icol += 1
                if 'B' in ab:
                    text_columns[ icol ].append( ( '{}{:02d}'.format( reps[2][2:4], get_primary_number( reps[2] ) ),
                                                   rep_colors[reps[2]] ) )
                    icol += 1
                    text_columns[ icol ].append( ( '{}{:02d}'.format( reps[3][2:4], get_primary_number( reps[3] ) ),
                                                   rep_colors[reps[3]] ) )
                    icol += 1
                    # text_columns[ icol ].append( ( 'BV{:02d}'.format( get_primary_number( reps[2] ) ), rep_colors[reps[2]] ) ) ; icol += 1
                    # text_columns[ icol ].append( ( 'BJ{}'.format( reps[3][4:])                      , rep_colors[reps[3]] ) ) ; icol += 1

                text_columns[ icol ].append( ( my_color_scores_labels[old_index], get_tcr_score_color( old_index ) ) ) ; icol += 1
                assert icol == num_columns


            ## now go through and figure out how wide each of the text columns is
            x_offset = xmargin
            for col,header_tag in zip( text_columns, header ):
                assert len(col) == num_tcrs
                maxlen = max((len(x[0]) for x in col ))
                if not maxlen: continue
                for old_index, ( text,color ) in enumerate(col):
                    new_index = old2new_index[ old_index ]
                    ypos = node_position[ new_index ]
                    lower_left = [ x_offset, ypos+0.5*label_fontsize*0.75 ]
                    cmds.append( svg_basic.make_text( text, lower_left, label_fontsize, color=color ) )
                if header_tag:
                    max_ypos = max( node_position.values() )
                    lower_left = [ x_offset, max_ypos+2.0*label_fontsize*0.75 ]
                    cmds.append( svg_basic.make_text( header_tag, lower_left, label_fontsize, color='black' ) )


                x_offset += text_column_separation + label_fontsize * 0.6 * maxlen


            ## how wide should the tree be?

            tree_p0 = [x_offset + labels_tree_separation, ymargin ]
            tree_p1 = [tree_p0[0] + tree_width, tree_height+ymargin ]

            total_svg_width = tree_p1[0] + xmargin

            ## node_position tells us where the different clusters are located, vertically
            ##

            plotter = svg_basic.SVG_tree_plotter()
            Log('Canvas_tree 2nd time')
            node_position,Transform,canvas_tree_min_rmsd, canvas_tree_w_factor = \
                score_trees_devel.Canvas_tree( tree, names, sizes, tree_p0, tree_p1, branch_width_fraction,
                                               plotter, label_internal_nodes = False,
                                               score_range_for_coloring = color_score_range )

            cmds.extend( plotter.cmds )


            ## label pvals
            if edge_pvals:
                plotting_info=(sizes, node_position, Transform, canvas_tree_w_factor, canvas_tree_min_rmsd)
                label_pval_edges( cmds, edge_pvals, tree, plotting_info )


            ## now we make an svg file
            prefix = '{}_tree_{}_{}_{}'.format(clones_file[:-4],ab,epitope,color_scheme)
            print 'create: {}.png'.format(prefix)

            svg_basic.create_file( cmds, total_svg_width, total_svg_height, prefix+'.svg', create_png=True)

