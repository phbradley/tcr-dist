from basic import *
import score_trees_devel
import svg_basic
from amino_acids import amino_acids
from tcr_distances import align_cdr3s
#from tcr_distances_blosum import blosum
import numpy as np
import util
import tcr_sampler ## for analyze_junction



with Parser(locals()) as p:
    #p.str('args').unspecified_default().multiple().required()
    p.str('clones_file').required()
    p.str('organism').required()
    p.str('outfile_prefix')
    p.float('distance_scale_factor').default(0.01) ## for sensible scale bars
    p.float('xmargin').default(10)
    p.float('ymargin').default(0)
    p.float('title_shift').default(40) ## only used if paper_figs is True
    p.str('color_scheme').default('probs')
    #p.int('int_arg').shorthand('i')
    #p.float('float_arg')     # --float_arg 9.6
    p.flag('verbose')       # --flag_arg  (no argument passed)
    p.flag('constant_seed')       # --flag_arg  (no argument passed)
    p.int('random_seed')       # --flag_arg  (no argument passed)
    p.flag('hacking')       # --flag_arg  (no argument passed)
    p.flag('paper_figs')       # --flag_arg  (no argument passed)
    p.flag('junction_bars')       # --flag_arg  (no argument passed)
    #p.range('range_arg')     # --range_arg 1:2
    #p.multiword('multi_arg') # --multi_arg hello world
    #p.file('file_arg')       # --file_arg README.txt
    #p.directory('dir_arg')   # --dir_arg /tmp/
    #p.str('floatlist').cast(lambda x: [float(val) for val in x.split(',')])
    #p.multiword('intlist').cast(lambda x: [int(val) for val in x.split()])
    p.multiword('ABs').cast(lambda x: x.split())
    p.multiword('epitopes').cast(lambda x: x.split())

if outfile_prefix is None:
    outfile_prefix = clones_file[:-4]

if constant_seed:
    random.seed(1)

if random_seed != None:
    print 'random_seed:',random_seed
    random.seed(random_seed)

if paper_figs and ABs == None:
    ABs = ['AB']

if ABs == None:
    ABs = ['A','B','AB']

fake_chains = util.detect_fake_chains( clones_file )
for ch in fake_chains:
    if ch in ABs:
        del ABs[ ABs.index(ch)]

gene_logo_name_trim = 2 if 'gammadelta' in pipeline_params['db_file'] else 4


font_family = "Droid Sans Mono"

greek_alpha = '&#x3b1;'
greek_beta  = '&#x3b2;'


junction_bars_color = { 'V':  'black',
                        'N1': 'red',
                        'N':  'red',
                        'N2': 'red',
                        'D':  'blue',
                        'J':  'gray' }

junction_bars_color = { 'V':  'black',
                        'N1': 'red',
                        'N':  'red',
                        'N2': 'red',
                        'D':  'blue',
                        'J':  'gray' }

junction_bars_color = { 'V':  'silver',
                        'N1': 'red',
                        'N':  'red',
                        'N2': 'red',
                        'D':  'black',
                        'J':  'dimgray' }


gap_character = '-' ## different from some other places
min_cluster_size = 1
min_cluster_size_for_glyphs = 5
min_cluster_fraction_for_glyphs = 0.03

max_covered_fraction_for_glyphs = 0.75

max_glyphs = 10

max_tcrs_for_trees = 300 ## if more, sub-sample
#max_tcrs_for_trees = 300 if paper_figs else 200


## this is a little silly, but it allows us to keep roughly the same clustering threshold when we update the distance measure
##
distance_threshold_25_scaled = distance_scale_factor * pipeline_params[ 'distance_threshold_25' ]

cluster_radius = {
    'A' : distance_threshold_25_scaled*2, ## was 0.5
    'B' : distance_threshold_25_scaled*2, ## was 0.5
    'AB': distance_threshold_25_scaled*4  ## was 1.0
}

vj_logo_width = 200
ypad = 300
xpad = 30
#xmargin = 10 #cmdline option now
xmargin_right = 10
glyphs_to_tree_spacer = 10
glyph_size_text_width = 40
pwmplusgaps_width = 770
tree_width = 1000
tree_height = 2000
pwm_height = 100
junction_bars_height = 35.0 if junction_bars else 0.0
ab_glyphs_spacer = 15


if True: #hacking
    fac = 0.5 if paper_figs else 0.333
    tree_width *= fac
    tree_height *= fac
    pwm_height *= fac
    junction_bars_height *= fac
    xpad *= fac
    ypad *= fac

    pwmplusgaps_width *= fac
    vj_logo_width *= fac

junction_bars_ypad = 3

glyph_height = pwm_height + junction_bars_height
boxpad = 2

if paper_figs: ypad = 30

## xmargin is the true margin and also the space between glyphs and tree

glyph_region_width = 2 * vj_logo_width + pwmplusgaps_width + 2 * xpad

tree_x0 = xmargin + glyph_size_text_width + glyph_region_width + glyphs_to_tree_spacer

total_svg_width = tree_x0 + tree_width + xmargin_right

branch_width_fraction = 0.3

log10_of_zero = -100



## parse the clones_file ############################################################################################3

all_tcr_infos = parse_tsv_file( clones_file, ['epitope'], [], True )

if not epitopes:
    epitopes = all_tcr_infos.keys()
    epitopes.sort()

all_tcrs = {}
all_color_scores = {}

for epitope in epitopes:
    infos = all_tcr_infos[epitope]
    util.assign_label_reps_and_colors_based_on_most_common_genes_in_repertoire( infos, organism )

    all_tcrs[epitope] = []
    all_color_scores[epitope] = []

    for l in infos:

        mouse = l['subject']
        epitope = l['epitope']
        cdr3a = l['cdr3a']
        cdr3b = l['cdr3b']

        ## note that we are using mm1 reps here that also dont have allele info
        va_rep = l['va_label_rep']
        ja_rep = l['ja_label_rep']
        vb_rep = l['vb_label_rep']
        jb_rep = l['jb_label_rep']

        if junction_bars:
            a_junction_results = tcr_sampler.analyze_junction( organism, l['va_gene'], l['ja_gene'],
                                                               cdr3a, l['cdr3a_nucseq'], return_cdr3_nucseq_src=True )
            b_junction_results = tcr_sampler.analyze_junction( organism, l['vb_gene'], l['jb_gene'],
                                                               cdr3b, l['cdr3b_nucseq'], return_cdr3_nucseq_src=True )

            cdr3a_new_nucseq, cdr3a_protseq_masked, cdr3a_protseq_new_nucleotide_countstring,\
                a_trims, a_inserts, cdr3a_nucseq_src = a_junction_results
            cdr3b_new_nucseq, cdr3b_protseq_masked, cdr3b_protseq_new_nucleotide_countstring,\
                b_trims, b_inserts, cdr3b_nucseq_src = b_junction_results
            ## try to distinguish between N before D and N after D
            for i in range(len(cdr3b_nucseq_src)):
                if cdr3b_nucseq_src[i] == 'N':
                    if cdr3b_nucseq_src[:i].count('D')==0:
                        cdr3b_nucseq_src[i] = 'N1'
                    else:
                        cdr3b_nucseq_src[i] = 'N2'
        else:
            cdr3a_nucseq_src = ['V']*(3*len(cdr3a)) ## hack, unused
            cdr3b_nucseq_src = ['V']*(3*len(cdr3b))

        assert len(cdr3a_nucseq_src) == 3*len(cdr3a)
        assert len(cdr3b_nucseq_src) == 3*len(cdr3b)

        all_tcrs[ epitope ].append( ( mouse, va_rep, ja_rep, vb_rep, jb_rep, cdr3a, cdr3b,
                                      cdr3a_nucseq_src, cdr3b_nucseq_src, l['clone_id'] ) )

        ## figure out the color scores, one for each chains
        if color_scheme == 'probs':
            ## color by probs
            def get_safe_log10(f ):
                if f==0: return log10_of_zero
                else: return math.log10(f)

            color_scores = { 'A': get_safe_log10( float(l['a_protseq_prob']) ),
                             'B': get_safe_log10( float(l['b_protseq_prob']) ),
                             'AB': get_safe_log10( float(l['a_protseq_prob']) * float( l['b_protseq_prob'] ) ) }

        all_color_scores[epitope].append( color_scores )





total_y_offset = ymargin
all_cmds = { 'A':[], 'B':[], 'AB':[] }

color_score_range = {'A':None, 'B':None, 'AB':None }

if color_scheme == 'probs':
    if organism == 'mouse':
        ## from ../distances.new.py
        color_score_range[ 'A'] = (  -9.0,  -5.0 )
        color_score_range[ 'B'] = (  -9.5,  -5.5 )
        color_score_range['AB'] = ( -16.0, -12.0 )
    else:
        assert organism == 'human'
        ## from ../distances.new.py
        ## try shifting down by 1 for a/b and 2 for ab
        color_score_range[ 'A'] = (  -10.0,  -6.0 )
        color_score_range[ 'B'] = (  -10.5,  -6.5 )
        color_score_range['AB'] = ( -20.0, -14.0 )



for epitope in epitopes:
    if hacking and epitope != 'M158': continue
    tcrs = all_tcrs[epitope]
    assert len(tcrs) == len(all_color_scores[epitope] )

    ## setup a coloring scheme for reps
    rep_colors = {}
    for info in all_tcr_infos[epitope]:
        for vj in 'vj':
            for ab in 'ab':
                rep   = info[vj+ab+'_label_rep']
                color =  info[vj+ab+'_label_rep_color']
                rep_colors[rep] = color


    total_tcrs = len(tcrs)

    ab_centers = {}
    ab_all_members = {}
    #if epitope != 'NP':continue
    #for ab in ['AB']:
    for ab in ABs:
        if hacking and ab != 'A': continue


        if color_scheme == 'probs':
            ## adjust scores of low-prob guys
            min_good_score = 0
            for color_scores in all_color_scores[epitope]:
                if color_scores[ab] != log10_of_zero:
                    min_good_score = min( min_good_score, color_scores[ab] )
            print 'min_good_score:',epitope,ab,min_good_score
            for color_scores in all_color_scores[epitope]:
                if color_scores[ab] == log10_of_zero:
                    color_scores[ab] = min_good_score

        radius = cluster_radius[ab]
        distfile = '{}_{}_{}.dist'.format( clones_file[:-4], ab, epitope )
        assert exists(distfile)


        N=0
        all_nbrs = []
        all_dists = []
        for line in open( distfile,'r'):
            l = line.split()
            clone_id = l[0]
            assert tcrs[ len(all_nbrs) ][-1] == clone_id
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

            #print epitope, radius, len(centers), best_nbr_count, uniqdlines[ center ][2:8]

            deleted[center] = True
            for nbr in all_nbrs[center]:
                if not deleted[nbr]:
                    deleted[nbr] = True
                    members.append( nbr )

            assert len(members) == best_nbr_count
            all_members.append( frozenset(members) )

        ab_centers[ab] = centers
        ab_all_members[ab] = all_members

        ## possibly subsample
        tree_indices = range(len(tcrs))
        if len(tree_indices) > max_tcrs_for_trees:
            tree_indices = random.sample( tree_indices, max_tcrs_for_trees )

        ## now get some info together for plotting
        all_center_dists = {}
        all_scores = []
        sizes = []
        names = []

        real_cluster_number2fake_cluster_number = {}
        fake_cluster_number2real_cluster_number = {}
        fake_ic=-1
        real_sizes = []
        for ic,center in enumerate(centers):
            ok_members = []
            for m in all_members[ic]:
                if m in tree_indices:
                    ok_members.append(m)
            if not ok_members:
                real_cluster_number2fake_cluster_number[ic] = -1
                continue
            fake_ic += 1
            real_cluster_number2fake_cluster_number[ic] = fake_ic
            fake_cluster_number2real_cluster_number[fake_ic] = ic
            assert fake_ic == len(sizes)
            size = len(ok_members)
            real_size = len(all_members[ic])
            all_scores.append( [ all_color_scores[epitope][x][ab] for x in ok_members ] )
            sizes.append( size )
            real_sizes.append( real_size )
            # if False and real_size>=my_min_cluster_size_for_glyphs:
            #     names.append( '{}'.format(real_size) )
            # else:
            names.append('')

        ## now do cluster center distances now that we know which ones are in the plot
        for ic,center in enumerate(centers):
            fake_ic = real_cluster_number2fake_cluster_number[ic]
            if fake_ic<0: continue
            for jc,other_center in enumerate(centers):
                fake_jc = real_cluster_number2fake_cluster_number[jc]
                if fake_jc<0: continue
                all_center_dists[ (fake_ic,fake_jc) ] = all_dists[center][other_center]



        print 'num_tcrs:',len(tcrs),'num_clusters:',len(centers),'fake_num_tcrs',sum(sizes),\
            'fake_num_clusters:',len(sizes)

        percentile = -1

        tree = score_trees_devel.Make_tree( all_center_dists, len(names),
                                            score_trees_devel.Update_distance_matrix_AL,
                                            all_scores, percentile )


        plotter = svg_basic.SVG_tree_plotter()

        tree_p0 = [tree_x0, total_y_offset + ypad]
        tree_p1 = [tree_x0 + tree_width, total_y_offset + ypad + tree_height ]



        ## node_position tells us where the different clusters are located, vertically
        ##
        node_position,Transform,canvas_tree_min_rmsd, canvas_tree_w_factor = \
            score_trees_devel.Canvas_tree( tree, names, sizes, tree_p0, tree_p1, branch_width_fraction,
                                           plotter, label_internal_nodes = False,
                                           score_range_for_coloring = color_score_range[ab] )

        max_rmsd_for_glyphs = 3.0*radius

        my_min_cluster_size_for_glyphs = max( min_cluster_size_for_glyphs,
                                              int(floor(0.5 + min_cluster_fraction_for_glyphs*total_tcrs)))

        while True:


            ## get all internal horizontal edges (ie subtrees) that are merged at an rmsd below a threshold
            ## and have at least my_min_cluster_size_for_glyphs
            def get_good_edges( subtree, node_position, sizes, real_sizes ):
                if score_trees_devel.IsALeaf( subtree ):
                    return []
                else:
                    big_rmsd = subtree[2]
                    edges = []
                    for ii in range(2):
                        iitree = subtree[ii]
                        little_rmsd = iitree[2]
                        center= score_trees_devel.Center( iitree, node_position, sizes, use_sizes_as_weights=True )
                        real_size = score_trees_devel.Size( iitree, real_sizes )
                        if little_rmsd <= max_rmsd_for_glyphs and real_size >= my_min_cluster_size_for_glyphs:
                            ## this tree is OK
                            edges.append( tuple( ( little_rmsd, big_rmsd, center, real_size,
                                                   tuple( sorted( score_trees_devel.Node_members(iitree) ) ) ) ) )
                        edges.extend( get_good_edges( iitree, node_position, sizes, real_sizes ) )
                return edges


            good_edges = get_good_edges( tree, node_position, sizes, real_sizes )

            # if True:
            #     for edge in good_edges:
            #         print 'good:',epitope, ab, edge
                #exit()

            ## figure out which edges we should draw glyphs for
            ## take all the true clusters
            ## try to maximize coverage?
            ##
            glyph_cmds = []
            glyph_edges = []
            while len(glyph_edges) < max_glyphs and len(glyph_edges)<len(good_edges):

                ## figure out who is covered
                covered = set()
                for e in glyph_edges:
                    for fake_ic in e[4]:
                        covered.add( fake_ic )

                new_edge = ()
                best_score = -1e6
                for e in good_edges:
                    if e in glyph_edges: continue
                    (little_rmsd,big_rmsd,center,real_size,clusters) = e
                    new_size = real_size
                    for fake_ic in clusters:
                        if fake_ic in covered:
                            new_size -= real_sizes[fake_ic]
                    if new_size<my_min_cluster_size_for_glyphs: continue
                    if (real_size-new_size) > max_covered_fraction_for_glyphs * real_size: continue

                    if little_rmsd<1e-3: ## true cluster
                        sortscore = 1000 + new_size
                    else:
                        sortscore = new_size - ( little_rmsd / max_rmsd_for_glyphs ) * total_tcrs
                    if sortscore > best_score:
                        new_edge = e[:]
                        best_score = sortscore
                if not new_edge: break
                #print 'new_edge:',new_edge
                glyph_edges.append( new_edge )

                ## let's make a box around this edge in the tree
                (little_rmsd,big_rmsd,center,real_size,clusters) = new_edge
                fake_size = sum( [sizes[x] for x in clusters] )
                line_width = max(1,int(floor(0.5+ fake_size*canvas_tree_w_factor )))
                box_x0 = Transform(max(canvas_tree_min_rmsd,little_rmsd)) ; box_x1 = Transform(big_rmsd)
                assert little_rmsd<big_rmsd
                glyph_cmds.append( svg_basic.rectangle( ( box_x0, center - line_width/2.0 ),
                                                        ( box_x1, center + line_width/2.0 ),
                                                        'none', 'black', stroke_width=3, dashed=True ) )
                glyph_cmds.append( svg_basic.make_text( '%d'%real_size, ( box_x1+4, center + 0.75*15.0/2) , 15, font_family=font_family))

            if len( glyph_edges ) < max_glyphs and my_min_cluster_size_for_glyphs>1:
                my_min_cluster_size_for_glyphs -= 1
            else:
                break

        cmds = glyph_cmds[:]

        glyph_location = dict( [ (e,e[2]) for e in glyph_edges ] )

        min_glyph_loc = tree_p0[1] + glyph_height/2
        max_glyph_loc = tree_p1[1] - glyph_height/2

        if False and paper_figs:
            ## scrunch everything down by a factor of 0.75
            scale = 0.68
            def new_loc( loc ):
                return max_glyph_loc - scale * ( max_glyph_loc - loc )
            min_glyph_loc = new_loc( min_glyph_loc )
            for e in glyph_location:
                glyph_location[e] = new_loc( glyph_location[e] )


        while True: ## keep looping until we are bump-free
            l = [(y,x) for x,y in glyph_location.iteritems() ]

            l.sort()
            bump = False
            stepsize = 5 ## pixels
            for (loc1,edge1),(loc2,edge2) in zip( l[:-1], l[1:] ):
                sep = loc2-loc1
                if sep < 1.25*glyph_height:
                    bump = True
                    if glyph_location[edge1]-stepsize > min_glyph_loc:
                        glyph_location[edge1] -= stepsize
                    if glyph_location[edge2]+stepsize < max_glyph_loc:
                        glyph_location[edge2] += stepsize
            if not bump:
                break

        #exit()

        ## now let's draw some cluster summary glyphs next to the cluster tree

        ## tree starts at total_y_offset+ypad
        if paper_figs:
            cmds.append( svg_basic.make_text( '{} ({} TCRs)'.format(epitope,total_tcrs),
                                              ( xmargin + glyph_size_text_width, total_y_offset+ypad+title_shift ),
                                              50, font_family=font_family ) ) ## label the epitope
            # cmds.append( svg_basic.make_text( '{} N={}'.format(epitope,total_tcrs),
            #                                   ( 0.5* ( tree_p0[0]+tree_p1[0] ) - 50, tree_p0[1] - 25 ),
            #                                   30, font_family=font_family ) ) ## label the epitope
        else:

            cmds.append( svg_basic.make_text( '{} {} #tcrs: {}'.format(epitope,ab,total_tcrs),
                                              ( 10, total_y_offset+ypad ) , 30,
                                              font_family=font_family ) ) ## label the epitope


        ## write out the glyph-cluster sizes
        for glyph_edge in glyph_edges:
            real_size = glyph_edge[3]
            loc = glyph_location[glyph_edge]
            ## silly xloc, was hard-coded to 5
            cmds.append( svg_basic.make_text( '%3d'%real_size, ( xmargin-5, loc+pwm_height/4.) , 20, font_family=font_family))



        for ii_ab2, ab2 in enumerate(ab):
            if len(ab) == 1:
                ab_fraction = 1.0
            else:
                ab_fraction = float(glyph_region_width-ab_glyphs_spacer) / ( 2*glyph_region_width )

            my_vj_logo_width     = ab_fraction * vj_logo_width
            my_pwmplusgaps_width = ab_fraction * pwmplusgaps_width
            my_xpad              = ab_fraction * xpad

            if ab2== 'A':
                v_index, j_index, cdr3_index, nucseq_src_index = 1, 2, 5, 7
                junction_bars_order = ['V','N','J']
                greek_letter = greek_alpha
            else:
                assert ab2=='B'
                v_index, j_index, cdr3_index, nucseq_src_index = 3, 4, 6, 8
                junction_bars_order = ['V','N1','D','N2','J']
                greek_letter = greek_beta

            ## the old way
            # for ic, (center,members) in enumerate( zip( centers, all_members ) ):
            #     fake_ic = real_cluster_number2fake_cluster_number[ic]
            #     if fake_ic<0: continue

            #     size = len(members)
            #     if size < my_min_cluster_size_for_glyphs: continue

            for glyph_edge in glyph_edges:
                (little_rmsd, big_rmsd, y_center, real_size, clusters ) = glyph_edge
                if hacking and real_size != 73:
                    continue
                ## which cluster center is the best choice
                min_rmsd = 1000
                best_fake_ic = -1
                members = []
                #clusters.sort() ## prefer center of larger clusters; clusters already sorted now
                if len(clusters)>1:assert clusters[0]<clusters[-1] and clusters[0]<clusters[1]
                for fake_ic in clusters:
                    ic = fake_cluster_number2real_cluster_number[ fake_ic ]
                    assert real_sizes[fake_ic] == len(all_members[ ic ])
                    members.extend( all_members[ ic ] )
                    rmsd=0.0
                    for fake_jc in clusters:
                        jc = fake_cluster_number2real_cluster_number[ fake_jc ]
                        rmsd += all_center_dists[(fake_ic,fake_jc)]
                        assert  all_center_dists[(fake_ic,fake_jc)] == all_dists[ centers[ic] ][ centers[jc] ]#sanity
                    if rmsd<min_rmsd-1e-3: ## take first in case of ties
                        min_rmsd = rmsd
                        best_fake_ic = fake_ic

                assert len(members) == real_size

                distl = []
                for m1 in members:
                    avgdis = 0.
                    maxdis = 0.
                    for m2 in members:
                        avgdis += all_dists[m1][m2]
                        maxdis = max(maxdis, all_dists[m1][m2] )
                    avgdis/=len(members)
                    distl.append( ( avgdis, m1 ) )
                    #distl.append( ( maxdis, m1 ) )
                distl.sort()
                distl_dict = dict( ( (y,x) for x,y in distl ) ) ## for debugging

                best_ic = fake_cluster_number2real_cluster_number[ best_fake_ic ]
                center = centers[ best_ic ] ## not used any more
                center = distl[0][1]
                #print 'best center:',center,best_ic,best_fake_ic,min_rmsd,clusters


                ## count v,j gene reps
                v_count = {}
                j_count = {}
                for rep in [tcrs[x][v_index] for x in members ]: v_count[rep] = v_count.get(rep,0)+1
                for rep in [tcrs[x][j_index] for x in members ]: j_count[rep] = j_count.get(rep,0)+1

                center_cdr3 = tcrs[center][ cdr3_index ][3:-2]
                if verbose:
                    print 'center_cdr3: {} {} {} {:15s} {:15s} {:9.3f} {:2d} {}'\
                        .format( epitope,ab,ab2,
                                 tcrs[center][v_index],tcrs[center][j_index],
                                 distl_dict[center],
                                 len(center_cdr3), center_cdr3 )
                L = len(center_cdr3)

                pwm = {}
                junction_pwm = {}
                gap_count = {}
                for i in range(L):
                    pwm[i] = dict(zip(amino_acids+[gap_character],[0]*21))
                    gap_count[i]=0
                for i in range(3*L):
                    junction_pwm[i] = dict( zip( junction_bars_order+[gap_character],
                                                 [0.]*(1+len(junction_bars_order))))

                for member in members:
                    member_cdr3 = tcrs[member][cdr3_index][3:-2]
                    member_junction = tcrs[member][nucseq_src_index][9:-6] ## a list
                    assert len(member_junction) == 3*len(member_cdr3)
                    a,b = align_cdr3s( center_cdr3, member_cdr3, gap_character )
                    if verbose:
                        print 'member_cdr3: {} {} {} {:15s} {:15s} {:9.3f} {:2d} {} {} {}'\
                            .format( epitope,ab,ab2,
                                     tcrs[member][v_index],tcrs[member][j_index],
                                     distl_dict[member],
                                     len( member_cdr3 ), member_cdr3, a, b )


                    for i in range(len(a)):
                        if a[i] == gap_character:
                            if i and a[i-1]==gap_character:continue
                            gap_count[i-1] += 1

                        else: ## b[i] could be gap_character or an amino_acid
                            pwmpos = i - a[:i].count(gap_character)
                            pwm[pwmpos][ b[i] ] += 1
                            for j in range(3*pwmpos,3*pwmpos+3):
                                if b[i] == gap_character:
                                    junction_pwm[j]['-'] += 1
                                else:
                                    bpos = i-b[:i].count(gap_character)
                                    assert b[i] == member_cdr3[bpos] ## sanity
                                    bsrc = member_junction[3*bpos + j-3*pwmpos ]
                                    junction_pwm[j][bsrc] += 1

                ## normalize the pwms
                for i in range(L):
                    tot = float(sum(pwm[i].values()))
                    for aa in pwm[i]:
                        pwm[i][aa] /= tot
                for i in range(3*L):
                    tot = float(sum(junction_pwm[i].values()))
                    for aa in junction_pwm[i]:
                        junction_pwm[i][aa] /= tot

                ## now we want to make a pwm
                total_gaps = sum(gap_count.values())


                column_width = my_pwmplusgaps_width / ( L + float(total_gaps)/len(members) )

                y0 = glyph_location[glyph_edge] - glyph_height/2.
                y1 = y0 + pwm_height
                # y0 = node_position[fake_ic] - pwm_height/2.
                # y1 = node_position[fake_ic] + pwm_height/2.

                ## make a v-gene logo
                vl = [(y,x[gene_logo_name_trim:],rep_colors[x]) for x,y in v_count.iteritems()]
                jl = [(y,x[gene_logo_name_trim:],rep_colors[x]) for x,y in j_count.iteritems()]
                #vl = [(y,x[4:x.index('*')]) for x,y in v_count.iteritems()]
                #jl = [(y,x[4:x.index('*')]) for x,y in j_count.iteritems()]

                single_glyph_width = 2* my_vj_logo_width + my_pwmplusgaps_width + 2*my_xpad
                x0 = xmargin + glyph_size_text_width + ii_ab2 * ( single_glyph_width + ab_glyphs_spacer )

                cmds.append( svg_basic.make_stack( (x0,y0), (x0 + my_vj_logo_width,y1), vl ) )

                ## a box around the V-logo
                cmds.append( svg_basic.rectangle( ( x0-boxpad,y0-boxpad), (x0 + my_vj_logo_width+boxpad,y1+boxpad ),
                                                  'none', 'black', stroke_width=1 ) )

                if junction_bars: ## label V-logo down below
                    text = 'V'+greek_letter
                    fontsize = junction_bars_height*0.9
                    p0 = [ x0 + 0.5* my_vj_logo_width - 0.6*fontsize, y1+junction_bars_height ]
                    cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )



                x0 += my_vj_logo_width + my_xpad

                cmds.append( svg_basic.rectangle( ( x0-boxpad,y0-boxpad), (x0 + my_pwmplusgaps_width+boxpad,y1+boxpad ),
                                                  'none', 'black', stroke_width=1 ) )



                #print ic,L,size,y0,y1

                ## now show each column
                prev_gap_column_width = 0.0
                jb_rights = [] #debugging
                for pos in range(L):
                    ## first the column of aas
                    colpwm={}
                    colpwm[0] = pwm[pos]
                    if verbose:
                        print 'colpwm:',pos,pwm[pos]
                    cmds.append( svg_basic.protein_logo( (x0,y0), (x0+column_width,y1), colpwm ) )

                    save_x0 = x0 ## for junction_bars

                    x0 += column_width

                    ## any gaps?
                    if gap_count[pos]:
                        gap_column_width = float( column_width * gap_count[pos] ) / len(members)
                        cmds.append( svg_basic.text_in_box( (x0,y0), (x0+gap_column_width,y1), gap_character, 'black' ) )

                        x0+= gap_column_width
                    else:
                        gap_column_width = 0.0


                    if junction_bars:
                        junction_bar_width = ( column_width + gap_column_width/2. + prev_gap_column_width/2. )/3.
                        junction_bar_x0 = save_x0 - prev_gap_column_width/2.
                        print 'left:',junction_bar_x0,'right:',junction_bar_x0 + 3.*junction_bar_width,\
                            'prev_gap_column_width:',prev_gap_column_width,'gap_column_width:',gap_column_width,\
                            'save_x0:',save_x0,'column_width:',column_width,'junction_bar_width:',junction_bar_width
                        if jb_rights:
                            assert abs( junction_bar_x0 - jb_rights[-1] )<1e-3
                        jb_rights.append( junction_bar_x0 + 3.*junction_bar_width )

                        y2 = y1+junction_bars_height
                        for j in range(3):
                            col = junction_pwm[3*pos+j]
                            lcol = [ ( col[x],x) for x in junction_bars_order ]
                            # lcol = [ (y,x) for x,y in col.iteritems()]
                            # lcol.sort()
                            # lcol.reverse()
                            y1shift = y1+ junction_bars_ypad
                            ## largest at the top
                            for frac,a in lcol:
                                if a==gap_character: continue
                                y1shift_next = y1shift + frac * junction_bars_height
                                color = junction_bars_color[ a ]
                                p0 = [ junction_bar_x0+ j   *junction_bar_width, y1shift]
                                p1 = [ junction_bar_x0+(j+1)*junction_bar_width, y1shift_next ]
                                cmds.append( svg_basic.rectangle( p0, p1, fill=color, stroke=color ) )
                                y1shift = y1shift_next


                    prev_gap_column_width = gap_column_width


                x0 += my_xpad

                ## now the J-logo
                cmds.append( svg_basic.make_stack( (x0,y0), (x0+my_vj_logo_width,y1), jl ) )

                cmds.append( svg_basic.rectangle( ( x0-boxpad,y0-boxpad), (x0 + my_vj_logo_width+boxpad,y1+boxpad ),
                                                  'none', 'black', stroke_width=1 ) )

                if junction_bars: ## label V-logo down below
                    text = 'J'+greek_letter
                    fontsize = junction_bars_height * 0.9
                    p0 = [ x0 + 0.5* my_vj_logo_width - 0.6*fontsize, y1+junction_bars_height ]
                    cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )




        # epitope_tree_cmds.extend(cmds)
        # epitope_tree_width = max( epitope_tree_width, total_svg_width )
        # epitope_tree_height += 2200

        all_cmds[ ab ].extend( plotter.cmds )
        all_cmds[ ab ].extend( cmds )

    total_y_offset += ypad + tree_height ## increment once per epitope


for ab in all_cmds:
    ## now make the svg file
    prefix = '{}_tall_tree_{}'.format(outfile_prefix,ab)

    svg_height = total_y_offset+ypad+ymargin
    svg_basic.create_file( all_cmds[ab], total_svg_width, svg_height, prefix+'.svg', create_png=True )

    if ab == 'AB':
        util.readme( prefix+'.png',"""
        These are TCRdist clustering trees for the different repertoires, with distances calculated over both chains. Down below this in the .html output are the clustering trees for distances calculated over each chain individually. To make these trees, the repertoire is first clustered with a fixed distance threshold (a TCRdist of {:.2f} for single chain distances and {:.2f} for paired alpha+beta chain distances) using a simple greedy approach that iteratively finds the TCR with the greatest number of neighbors within that distance, adds that cluster center and its neighbors to the list of clusters and deletes them from the repertoire, repeating until all the receptors have been clustered. These clusters are the leaves of the tree, and they are joined together for visualization purposes by an average-linkage hierarchical clustering approach that uses the matrix of distances between the cluster centers. The vertical thickness of the leaves and branches is proportional to the number of TCR clones represented by those branches. Repertoires with more than 300 receptors are subsampled to 300 after clustering so that the number of leaves in the tree doesn't become too large. Trees with all the TCR clones represented can be found by following the links toward the top of the .html output.<br><br>

        TCR logos are shown to the left of the tree for a representative subset of the branches (enclosed in dashed boxes labeled with their size; logos and boxed branches come in the same vertical order). Each logo panel shows the V- (left) and J- (right) gene frequencies in 'logo' format (height scaled by frequency, most frequent gene at the top; the IMGT gene names are trimmed to remove the leading TRAV/TRAJ/TRBV/TRBJ). In the middle is a CDR3 amino acid sequence logo. The colored bars below the CDR3 logo summarize the inferred rearrangement history of the grouped receptors by showing the nucleotide source, colored as follows: V region, light gray; J region, dark gray; D region, black; N insertions, red. THe number of TCRs contributing to the logo is shown to the left and should match the number next to the corresponding boxed branch of the tree.<br><br>

        This analysis is conducted at the level of clonotypes -- each expanded clone is condensed to a single receptor sequence for the purpose of clustering.

        """.format( cluster_radius['A'] / distance_scale_factor,
                    cluster_radius[ab] / distance_scale_factor ))
