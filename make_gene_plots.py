from basic import *
import html_colors
import svg_basic
import util
import numpy as np



with Parser(locals()) as p:
    p.str('clones_file').required()
    p.multiword('epitopes').cast(lambda x:x.split())
    p.multiword('force_pairing_order').cast(lambda x:x.split()).described_as('Left-right order for segments in gene usage diagrams')
    p.str('organism').required()
    p.str('outfile_prefix')
    p.str('countrep_enrichments_file')
    p.str('defs_id_prefix').default('')
    p.float('min_ami_for_colorscale').default(0.114) # from shuffling experiments
    p.float('max_ami_for_colorscale').default(0.5)
    p.float('min_entropy_for_colorscale').default(0.0)
    p.float('max_entropy_for_colorscale').default(5.0)
    p.float('min_jsd_for_colorscale').default(0.02259) ## from background vs background comparisons
    p.float('max_jsd_for_colorscale').default(0.0)
    p.float('min_gene_frequency_for_labels').default(0.05)
    p.float('vj_pairings_left_margin').default(50)
    p.float('vj_pairings_top_margin').default(50)
    p.flag('use_color_gradients')
    p.flag('reverse_gradients')
    p.flag('no_pairing_text')
    p.flag('paper_figs')
    p.flag('paper_supp')
    p.flag('consistentfigcolors')
    p.set_help_prefix("""
    This script makes a set of plots that illustrate gene segment usage in the dataset. The default prefix for plot output is the name of the clones file with the .tsv trimmed off.

    Plot descriptions:

    <prefix>_cdr3lens.png: Visualizes the length distribution of the CDR3 segments, colored by gene usage

    <prefix>_gene_segment_pies.png: Shows the gene segment usage for each repertoire as pie plots.

    <prefix>_gene_entropies_and_mi.png: Heat maps of gene usage distribution entropies, differences from background, and mutual information (reflecting covariation between gene usage in different segments)

    <prefix>_vj_pairings.png: 'Chord' style diagrams showing gene usage and covariation in graphical format.

    """)



#if paper_supp:
#    paper_figs = True ## NOTE

if not countrep_enrichments_file and exists( clones_file[:-4]+'_gene_probs.tsv' ):
    countrep_enrichments_file = clones_file[:-4]+'_gene_probs.tsv'
    print 'countrep_enrichments_file:',countrep_enrichments_file


if not outfile_prefix:
    outfile_prefix = clones_file[:-4]

#import numpy as np

segtypes = segtypes_uppercase[:] ## local use


pval_threshold_for_plotting_gene_correlations = 1e-2
pval_threshold_for_svg_correlations = 1e-6

#num_tcrs_to_choose_randomly = 100
#num_random_repeats = 2 ## I guess this wasn't useful (see ../cluster_dists_clusters.py)

greek_alpha = '&#x3b1;'
greek_beta  = '&#x3b2;'

segtype2greek_label = { 'VA':'V'+greek_alpha, 'JA':'J'+greek_alpha,
                        'VB':'V'+greek_beta , 'JB':'J'+greek_beta }



## load epitope jsd values
epitope_jsds = {}
jsd_tsvfile = clones_file[:-4] + '_JS_divergence.tsv'
if not exists( jsd_tsvfile ):
    print 'Sorry, you need to run analyze_gene_frequencies.py before running make_gene_plots.py'
    exit()

lines = parse_tsv_file( jsd_tsvfile, [], ['epitope'] + [x+'_jsd_normed' for x in segtypes_lowercase] )
for line in lines:
    epitope = line[0]
    vals = map(float,line[1:])
    epitope_jsds[epitope] = {}
    assert len(vals)== len(segtypes)
    for segtype,val in zip( segtypes, vals ):
        epitope_jsds[ epitope ][ segtype ] = val

epitope_entropies = {}
epitope_mis = {}
epitope_correlations = {}
epitope_correlations_svg = {}
epitope_repcounts = {}
epitope_repcounts_by_len = {}

min_cdr3len = 100
max_cdr3len = 0

all_tcrs = parse_tsv_file( clones_file, ['epitope'], [], True )

gradient_id_counter = 0
## returns id, cmd
def linear_gradient_cmd( x1, y1, x2, y2, offsets, colors, spread_method="pad" ):

    global gradient_id_counter
    global defs_id_prefix
    gradient_id_counter += 1
    id = '{}lingrad{:d}'.format(defs_id_prefix,gradient_id_counter)

    stoplines = ''
    assert len(offsets) == len(colors)

    for offset,color in zip( offsets,colors):
        stoplines += """
      <stop offset="{:.1f}%"   stop-color="{}" stop-opacity="1"/>
""".format( offset, color )

    cmd = """
  <defs>
    <linearGradient id="{}"
                    x1="{:.1f}%" y1="{:.1f}%"
                    x2="{:.1f}%" y2="{:.1f}%"
                    spreadMethod="{}">
{}
    </linearGradient>
  </defs>
""".format( id, x1, y1, x2, y2, spread_method, stoplines )

    return id, cmd


def roundlo(x): return int(floor(x))
def roundhi(x): return int(floor(1.0+x-1e-6))


## pixels ####################################################
left_margin = vj_pairings_left_margin #default is 50
right_margin = 50
top_margin = vj_pairings_top_margin #default is 50
bottom_margin = 50
yspacer = 50

flat_band = 50
final_flat_band = flat_band if use_color_gradients else 2.5*flat_band
middle_band = 400
slope_weight = 100
pairing_svg_y_offset = top_margin
pairing_svg_cmds = []
path_def_counter = 0
make_enrichment_glyphs = ( countrep_enrichments_file != None )
if make_enrichment_glyphs:
    all_countrep_enrichment = parse_tsv_file( countrep_enrichments_file, [ 'epitope','gene' ], ['jsd_prob_enrich'] )


if not epitopes:
    epitopes = all_tcrs.keys()[:]
    epitopes.sort()

for epitope in epitopes:

    ## this fills in *_label_rep fields in all_tcrs dictionary
    util.assign_label_reps_and_colors_based_on_most_common_genes_in_repertoire( all_tcrs[epitope], organism )

    epitope_entropies[epitope] = {}
    epitope_mis[epitope] = {}
    epitope_correlations[epitope] = []
    epitope_repcounts[epitope] = {}
    epitope_correlations_svg[epitope] = {}

    tcrs = []
    for fulltcr in all_tcrs[epitope]:
        tcrs.append( ( fulltcr['va_label_rep'], fulltcr['ja_label_rep'],
                       fulltcr['vb_label_rep'], fulltcr['jb_label_rep'],
                       len(fulltcr['cdr3a']), len(fulltcr['cdr3b'] ) ) ) # not subtracting 5 any more


    repcounts = {}
    repcounts2 = {}

    repcounts_by_len = {}

    for i,r in enumerate(segtypes):
        repcounts[r] = {}
        repcounts_by_len[r] = {}
        for s in segtypes[i+1:]:
            repcounts2[(r,s)] = {}

    rep_index = dict(zip(segtypes,range(len(segtypes))))

    for tcr in tcrs:
        assert len(tcr) == 6
        for r in segtypes:
            rep = tcr[ rep_index[r] ]
            repcounts[r][rep] = repcounts[r].get(rep,0)+1

            assert r[1] in 'AB'
            cdr3len = tcr[4] if r[1]=='A' else tcr[5]

            min_cdr3len = min(min_cdr3len,cdr3len)
            max_cdr3len = max(max_cdr3len,cdr3len)

            if cdr3len not in repcounts_by_len[r]:
                repcounts_by_len[r][cdr3len] = {}
            repcounts_by_len[r][cdr3len][rep] = repcounts_by_len[r][cdr3len].get(rep,0)+1

        for rs in repcounts2:
            rep = (tcr[ rep_index[rs[0]]], tcr[ rep_index[rs[1]]] )
            repcounts2[rs][rep] = repcounts2[rs].get(rep,0)+1

    for r in segtypes:
        for s in segtypes:
            rs=(r,s)
            if rs in repcounts2:
                for rep1 in repcounts[r]:
                    for rep2 in repcounts[s]:
                        rep=(rep1,rep2)
                        if rep not in repcounts2[rs]:
                            repcounts2[rs][rep] = 0



    epitope_repcounts[epitope] = dict( repcounts )
    epitope_repcounts_by_len[epitope] = dict( repcounts_by_len )

    N = len(tcrs)

    ## compute entropies, mutual informations
    for r in segtypes:
        entropy=0
        for rep,count in repcounts[r].iteritems():
            prob=float(count)/N
            entropy -= prob * math.log(prob,2)
        print 'ENT {:4s} {} entropy: {:7.3f} entropy_pow2: {:7.3f} N: {:6d}'.format(epitope,r,entropy,2**entropy,N)
        epitope_entropies[epitope][r] = entropy


    from sklearn.metrics import adjusted_mutual_info_score
    from scipy.stats import hypergeom

    all_ab_amis = []
    all_amis = {}

    for rs in repcounts2:
        ab_pairing = ( rs[0][1] != rs[1][1] )
        cluster_pairing = ( rs[0][0] == 'C' or rs[1][0] == 'C' )

        mi=0.0
        entropy=0
        for (rep1,rep2),count in repcounts2[rs].iteritems():
            pxy = float(count)/N
            if pxy>0: entropy -= pxy*math.log(pxy,2)
            count1 = repcounts[rs[0]][rep1]
            count2 = repcounts[rs[1]][rep2]
            px = float(count1)/N
            py = float(count2)/N
            if pxy>0: mi += pxy * math.log( (pxy/ (px*py)), 2 )

            ## lets look at the significance of this overlap
            expected = px * py * N
            pval = 1

            if count > expected:
                ## compute hypergeometric distn prob
                max_possible_overlap = min(count1,count2)
                x = np.arange(0,max_possible_overlap+1)
                cdf = hypergeom.cdf( x, N, count1, count2 ) ## cdf is accumulated prob <= val
                sf  = hypergeom.sf( x, N, count1, count2 )
                pval = sf[count-1] ## now greater than or equal to count
                if pval<1e-3:
                    print 'PVAL: {:4s} {:12.3e} {}-{} {:15s} {:15s} overlap: {:4d} expect: {:7.1f} count1: {:4d} count2: {:4d} '\
                        .format(epitope,pval,rs[0],rs[1],str(rep1),str(rep2),count,expected,count1,count2)
                #exit()

                if pval<pval_threshold_for_svg_correlations:
                    #print 'svg pval!',rep1,rep2,pval
                    epitope_correlations_svg[epitope][(rep1,rep2)] = ( pval, count/expected )
                    epitope_correlations_svg[epitope][(rep2,rep1)] = ( pval, count/expected )

            if count < expected:
                ## compute hypergeometric distn prob
                max_possible_overlap = min(count1,count2)
                x = np.arange(0,max_possible_overlap+1)
                cdf = hypergeom.cdf( x, N, count1, count2 ) ## cdf is accumulated prob <= val
                sf  = hypergeom.sf( x, N, count1, count2 )
                pval = cdf[count] ## less than or equal to count
                if pval<1e-3:
                    print 'PVAL: {:4s} {:12.3e} {}-{} {:15s} {:15s} overlap: {:4d} expect: {:7.1f} count1: {:4d} count2: {:4d} '\
                        .format(epitope,pval,rs[0],rs[1],str(rep1),str(rep2),count,expected,count1,count2)
                #exit()
                if pval<pval_threshold_for_svg_correlations:
                    #print 'svg pval!',rep1,rep2,pval
                    epitope_correlations_svg[epitope][(rep1,rep2)] = ( pval, count/expected )
                    epitope_correlations_svg[epitope][(rep2,rep1)] = ( pval, count/expected )


            if ab_pairing and (not cluster_pairing) and pval<pval_threshold_for_plotting_gene_correlations:
                if count==0:
                    logenrich = math.log(  0.25 / expected, 2 )
                else:
                    logenrich = math.log( count / expected, 2 )

                epitope_correlations[epitope].append ( ( logenrich, -1*math.log( pval,10 ), rs, rep1, rep2 ) )


        ## compute an adjusted mutual information score
        labels0 = []
        labels1 = []

        tcr_labels0 = []
        tcr_labels1 = []

        for tcr in tcrs:
            l0 = tcr[ rep_index[ rs[0] ] ]
            l1 = tcr[ rep_index[ rs[1] ] ]
            if l0 not in labels0: labels0.append( l0 )
            if l1 not in labels1: labels1.append( l1 )
            tcr_labels0.append( labels0.index(l0) )
            tcr_labels1.append( labels1.index(l1) )

        ami = adjusted_mutual_info_score( tcr_labels0, tcr_labels1 )


        if ab_pairing:
            all_ab_amis.append( ( ami, rs ) )

        all_amis[ (rs[0],rs[1]) ] = ami
        all_amis[ (rs[1],rs[0]) ] = ami

        print 'MI {:4s} {}-{} MI: {:7.3f} AMI: {:7.3f} MI_pow2 {:7.3f} entropy: {:7.3f} entropy_pow2: {:7.3f}'\
            .format(epitope,rs[0],rs[1],mi,ami,2**mi,entropy,2**entropy)

        epitope_entropies[epitope][rs] = entropy
        epitope_mis[epitope][rs] = (mi,ami)



    all_ab_amis.sort()
    all_ab_amis.reverse()

    top_pairing = all_ab_amis[0]
    print 'top ab pairing:',top_pairing

    middle_alpha = top_pairing[1][0]
    middle_beta  = top_pairing[1][1]
    assert middle_alpha in ['VA','JA']
    assert middle_beta in ['VB','JB']
    other_alpha = 'JA' if middle_alpha=='VA' else 'VA'
    other_beta  = 'JB' if middle_beta =='VB' else 'VB'


    ypixel_scale = max(1,int( 0.5 + 600.0/len(tcrs) ) )

    if paper_figs:
        ypixel_scale = 600.0/len(tcrs)

    elif paper_supp:
        ypixel_scale = 900.0/len(tcrs)

    #hacking
    #slope_weight = 1
    #ypixel_scale = 1


    pairing_svg_width = left_margin + right_margin + 3*(flat_band+middle_band) + final_flat_band

    if force_pairing_order:
        assert len(force_pairing_order) == 4
        reps = force_pairing_order[:]
    else:
        reps = [ other_alpha, middle_alpha, middle_beta, other_beta ]

    ff='Droid Sans Mono'

    if paper_figs:
        epitope_fontsize = 60
        midpoint = left_margin + 2*flat_band + 1.5*middle_band
        pairing_svg_cmds.append( svg_basic.make_text( '{}'.format(epitope,len(tcrs) ),
                                                      [midpoint-0.5*0.6*epitope_fontsize*len(epitope),
                                                       pairing_svg_y_offset+epitope_fontsize-20], epitope_fontsize,
                                                      font_family=ff ) )
    elif paper_supp:
        epitope_fontsize = 30
        midpoint = left_margin + 2*flat_band + 1.5*middle_band
        pairing_svg_cmds.append( svg_basic.make_text( '{}'.format(epitope,len(tcrs) ),
                                                      [midpoint-0.5*0.6*epitope_fontsize*len(epitope),
                                                       pairing_svg_y_offset+epitope_fontsize-20], epitope_fontsize,
                                                      font_family=ff ) )
    else:
        pairing_svg_cmds.append( svg_basic.make_text( '{} num_clones= {} ({}x y-pixel scale)'.format(epitope,len(tcrs),ypixel_scale),
                                                      [left_margin, pairing_svg_y_offset+20], 20, font_family=ff ) )



    correlation_fontsize = 16. if paper_supp else 14.
    correlation_fontheight = correlation_fontsize*0.75

    for ii in range(3):
        correlation_paths = []
        r0 = reps[ii]
        r1 = reps[ii+1]
        ami = all_amis[ (r0,r1)]

        x0 = left_margin + ii*( flat_band + middle_band )

        if paper_figs or paper_supp:
            text = segtype2greek_label[ r0 ]
            fontsize = 40. if paper_figs else 20.
            xtext = x0+0.5*flat_band-0.5*0.6*fontsize*2
            ytext = pairing_svg_y_offset+yspacer-6
            ## hacking
            ytext -= 6
            if ii==0:
                xtext += 8
            pairing_svg_cmds.append( svg_basic.make_text( text, [ xtext, ytext ], fontsize, font_family=ff ) )
            if ii==2: ## add the final column label
                text = segtype2greek_label[ r1 ]
                xtext = x0+1.5*flat_band-0.5*0.6*fontsize*2+middle_band
                xtext -= 8
                pairing_svg_cmds.append( svg_basic.make_text( text, [ xtext, ytext ], fontsize, font_family=ff ) )
        else:
            pairing_svg_cmds.append( svg_basic.make_text( r0, [x0+5, pairing_svg_y_offset+yspacer-3],
                                                          20, font_family=ff ) )
            if ii==2:
                pairing_svg_cmds.append( svg_basic.make_text( r1, [x0+flat_band+middle_band,
                                                                   pairing_svg_y_offset+yspacer-3], 20, font_family=ff))
        if not paper_figs:
            pairing_svg_cmds.append( svg_basic.make_text( '(AMI: {:.2f})'.format(ami),
                                                          [x0+flat_band+middle_band/2.5, pairing_svg_y_offset+yspacer-5],
                                                          12, font_family=ff  ))


        vl = [ (y,x) for x,y in repcounts[r0].iteritems() ]
        jl = [ (y,x) for x,y in repcounts[r1].iteritems() ]

        vl.sort() ; vl.reverse()
        jl.sort() ; jl.reverse()

        vcolors = dict(zip( [x[1] for x in vl], html_colors.get_rank_colors_no_lights( len(vl) ) ) )
        jcolors = dict(zip( [x[1] for x in jl], html_colors.get_rank_colors_no_lights( len(jl) ) ) )

        reps2tcrs = {}

        for t in tcrs:
            vj = ( t[ rep_index[r0] ], t[ rep_index[r1] ] )
            if vj not in reps2tcrs:reps2tcrs[vj] = []
            reps2tcrs[vj].append( t )

        ## on the left, the V-segments, ordered by counts
        ## on the right, J-segments, ordered by counts

        ## need to assign a vertical range to each v/j segment
        ## start with, one pixel per tcr
        ##


        jcounts = {}

        yleft=yspacer+pairing_svg_y_offset
        for vcount,v in vl:
            y0_right = yspacer+pairing_svg_y_offset
            vcolor = vcolors[v]
            for jcount,j in jl:
                vj=(v,j)
                jcolor = jcolors[j]

                num_tcrs = len(reps2tcrs.get(vj,[]))
                num_tcrs_scaled = num_tcrs * ypixel_scale

                if True:
                    stroke_width = roundhi(num_tcrs_scaled)

                    ## ok make a spline
                    yright = y0_right + jcounts.get(j,0)*ypixel_scale

                    #line/spline points
                    j_flat_band = flat_band if ii<2 else final_flat_band
                    points = [ (roundlo(x0), yleft + 0.5*num_tcrs_scaled ),
                               (x0+flat_band, yleft+0.5*num_tcrs_scaled ),
                               (roundhi(x0+flat_band+middle_band), yright+0.5*num_tcrs_scaled ),
                               (roundhi(x0+flat_band+middle_band+j_flat_band), yright+0.5*num_tcrs_scaled ) ]


                    path1_cmds = 'M {} {} L {} {} M {} {} C {} {}, {} {}, {} {}'\
                        .format( points[0][0], points[0][1], ## start of v-line
                                 points[1][0], points[1][1], ## end point of v-line
                                 points[1][0], points[1][1],
                                 points[1][0] +slope_weight, points[1][1], ## control for spline start
                                 points[2][0] -slope_weight, points[2][1], ## control for spline end
                                 points[2][0], points[2][1] )

                    if num_tcrs:
                        if use_color_gradients:
                            path1a_cmds = 'M {} {} L {} {}'\
                                .format( points[0][0], points[0][1],  ## start of v-line
                                         points[1][0], points[1][1] ) ## end point of v-line
                            pairing_svg_cmds.append( '<path d="{}" stroke="{}" stroke-width="{}" fill="none"/>'\
                                                     .format(path1a_cmds,vcolor, stroke_width ) )
                            ## define the gradient
                            path1b_cmds = 'M {} {} C {} {}, {} {}, {} {}'\
                                .format( points[1][0], points[1][1],
                                         points[1][0] +slope_weight, points[1][1], ## control for spline start
                                         points[2][0] -slope_weight, points[2][1], ## control for spline end
                                         points[2][0], points[2][1] )
                            #v_line_rhs_fraction = float(flat_band) / (flat_band + middle_band )
                            offsets = [0, 25.0, 75.0, 100]
                            #offsets = [0, 45.0, 55.0, 100]
                            #offsets = [0, 90.0, 99.0, 100]
                            if reverse_gradients:
                                colors = [jcolor, jcolor, vcolor, vcolor]
                            else:
                                colors = [vcolor, vcolor, jcolor, jcolor]
                            gradient_id, gradient_cmd = linear_gradient_cmd( 0, 0, 100, 0, offsets, colors )
                            pairing_svg_cmds.append( gradient_cmd )
                            pairing_svg_cmds.append( '<path d="{}" stroke="url(#{})" stroke-width="{}" fill="none"/>'\
                                                     .format(path1b_cmds, gradient_id, stroke_width ) )
                        else:
                            pairing_svg_cmds.append( '<path d="{}" stroke="{}" stroke-width="{}" fill="none"/>'\
                                                     .format(path1_cmds,vcolor, stroke_width ) )

                        if ii==2: ## add the right-most flat band
                            path2_cmds = 'M {} {} L {} {}'\
                                .format( points[2][0], points[2][1], ## start of j-line
                                         points[3][0], points[3][1] ) ## end of j-line

                            pairing_svg_cmds.append( '<path d="{}" stroke="{}" stroke-width="{}" fill="none"/>'\
                                                     .format(path2_cmds,jcolor, stroke_width) )


                    if vj in epitope_correlations_svg[epitope] and not paper_figs:
                        #print 'vj has correlations:',vj,epitope_correlations_svg[epitope][vj]
                        if not num_tcrs:
                            #print 'make dotted line!',vj
                            if not ( no_pairing_text or paper_figs ):
                                if paper_supp:
                                    assert use_color_gradients
                                    ## define the gradient
                                    path1b_cmds = 'M {} {} C {} {}, {} {}, {} {}'\
                                        .format( points[1][0], points[1][1],
                                                 points[1][0] +slope_weight, points[1][1], ## control for spline start
                                                 points[2][0] -slope_weight, points[2][1], ## control for spline end
                                                 points[2][0], points[2][1] )
                                    #v_line_rhs_fraction = float(flat_band) / (flat_band + middle_band )
                                    offsets = [0, 25.0, 75.0, 100]
                                    #offsets = [0, 45.0, 55.0, 100]
                                    #offsets = [0, 90.0, 99.0, 100]
                                    colors = [vcolor, vcolor, jcolor, jcolor]
                                    gradient_id, gradient_cmd = linear_gradient_cmd( 0, 0, 100, 0, offsets, colors )
                                    pairing_svg_cmds.append( gradient_cmd )
                                    pairing_svg_cmds.append( '<path d="{}" stroke="url(#{})" stroke-width="2" stroke-dasharray="5,5" fill="none"/>'\
                                                             .format(path1b_cmds, gradient_id, stroke_width ) )

                                else:
                                    dotted_path_cmds = 'M {} {} L {} {} M {} {} C {} {}, {} {}, {} {}'\
                                        .format( points[0][0], points[0][1], ## start of v-line
                                                 points[1][0], points[1][1], ## end point of v-line
                                                 points[1][0], points[1][1],
                                                 points[1][0] +slope_weight, points[1][1], ## control for spline start
                                                 points[2][0] -slope_weight, points[2][1], ## control for spline end
                                                 points[2][0], points[2][1] )
                                    pairing_svg_cmds.append( '<path d="{}" stroke="{}" stroke-width="2" stroke-dasharray="5,5" fill="none"/>'\
                                                             .format(path1_cmds,vcolor ) )

                        ## new way, just use regular text elements
                        ## pretend that the spline is actually a straight line between these points
                        swf=0.4
                        yshift = correlation_fontheight*0.5
                        p0 = ( points[1][0]+slope_weight*swf, points[1][1]+yshift )
                        p1 = ( points[2][0]-slope_weight*swf, points[2][1]+yshift )

                        dx = p1[0]-p0[0]
                        dy = p1[1]-p0[1]

                        ## so, what is the rotation we'll need?
                        rotangle = math.atan2(dy,dx) * ( 180.0 / math.pi )

                        step = 0.05
                        lower_left  = [ p0[0] + step*dx, p0[1] + step*dy ]

                        step = 0.95
                        lower_right = [ p0[0] + step*dx, p0[1] + step*dy ]

                        pval,enrich = epitope_correlations_svg[epitope][vj]

                        ## write some curved text
                        if enrich==0:
                            msg = '0x ({:.0E})'.format(pval)
                        elif enrich<0.1:
                            msg = '{:.2f}x ({:.0E})'.format(enrich,pval)
                        else:
                            msg = '{:.1f}x ({:.0E})'.format(enrich,pval)

                        fill1,fill2 = 'black','black'
                        if vcolor=='black':
                            fill1 = 'gold'
                        if jcolor=='black' and use_color_gradients:
                            fill2 = 'gold'


                        cmd1 = '<text x="{:.3f}" y="{:.3f}" font-size="{}" font-family="{}" fill="{}" transform="rotate({:.3f},{:.3f},{:.3f})" >{}</text>\n'\
                            .format( lower_left[0], lower_left[1], correlation_fontsize, ff, fill1,
                                     rotangle, lower_left[0], lower_left[1], msg )

                        cmd2 = '<text text-anchor="end" x="{:.3f}" y="{:.3f}" font-size="{}" font-family="{}" fill="{}" transform="rotate({:.3f},{:.3f},{:.3f})" >{}</text>\n'\
                            .format( lower_right[0], lower_right[1], correlation_fontsize, ff, fill2,
                                     rotangle, lower_right[0], lower_right[1], msg )

                        correlation_paths.append( ( pval, (cmd1, cmd2) ) )
                        #print 'corr cmd1:',vj,cmd1
                        #print 'corr cmd2:',vj,cmd2


                    yleft += num_tcrs_scaled
                    jcounts[j] = jcounts.get(j,0)+num_tcrs
                y0_right += jcount * ypixel_scale


        ## try doing the p-val paths
        correlation_paths.sort()
        correlation_paths.reverse() ## go in decreasing order of p-val so the most significant are on top

        ## now write the text
        for (pval, cmds ) in correlation_paths:
            pairing_svg_cmds.extend( cmds )

        ## let's label the alleles in the left stack (and right stack if ii==2)
        fontsize = 40 if paper_figs else 20.0 if paper_supp else 20
        fontheight = 0.75*fontsize
        fontwidth =  0.6 *fontsize

        min_height_for_labels = fontheight+1

        for jj,(r,ll,repcolors) in enumerate( [ (r0,vl,vcolors),(r1,jl,jcolors)]  ):
            if ii<2 and jj>0:continue


            ## label in white?
            x = x0 + jj*(flat_band+middle_band)
            ystart = yspacer+pairing_svg_y_offset

            for ( count,rep) in ll:
                if count*ypixel_scale < min_height_for_labels: break
                #ystop    = ystart + count*ypixel_scale
                midpoint =  ystart + count*ypixel_scale*0.5
                text = rep[2:]

                lower_left = [ x+2, midpoint+fontheight/2.0 ]

                my_flat_band = final_flat_band if ii==2 and jj==1 else flat_band
                bgcolor = repcolors[rep]
                textcolor = 'black' if ((paper_figs or paper_supp) and bgcolor!= 'black') else 'white'
                textcolor = 'black' if bgcolor!= 'black' else 'white'

                if True or paper_figs or paper_supp: ## center the text, unless on either side...
                    text_width = fontwidth*len(text)
                    lower_left_ha = {'left'  : lower_left,
                                     'right' : [ x+my_flat_band-text_width, midpoint+fontheight/2.0 ],
                                     'center': [ x+0.5*my_flat_band-0.5*text_width, midpoint+fontheight/2.0 ]}
                    if jj==0 and ii==0: ## left-most guy
                        ha = 'left'
                    elif jj==1 and ii==2: ## right-most guy
                        ha = 'right'
                    else:
                        ha = 'center'
                    pairing_svg_cmds.append( svg_basic.make_text( text, lower_left_ha[ha], fontsize, color=textcolor,
                                                                  font_family=ff))


                elif (True or jj==1) and fontwidth*len(text)>my_flat_band: # right-most set, dont want to over-run
                    myfontsize=int(0.5+(my_flat_band-4)/(len(text)*0.6))
                    pairing_svg_cmds.append( svg_basic.make_text( text, lower_left, myfontsize, color=textcolor,
                                                                  font_family=ff))
                else:
                    pairing_svg_cmds.append( svg_basic.make_text( text, lower_left, fontsize, color=textcolor,
                                                                  font_family=ff))



                ## add an enrichment glyph?
                if make_enrichment_glyphs:
                    enrich = float( all_countrep_enrichment[ epitope ][ rep ][0][0] )
                    if enrich>=2. or enrich<=0.5:
                        ## add a glyph
                        if paper_supp or paper_figs:
                            arrow_length = 1.35 * min_height_for_labels
                            arrow_width = 3.5
                        else:
                            arrow_length = 1.35 * min_height_for_labels
                            arrow_width = 1.5
                        #arrow_length = min_height_for_labels
                        #arrow_width = 2.5
                        eg_sep = 14.0
                        if 'A' in r:
                            center = [ lower_left_ha[ha][0] + text_width + eg_sep, midpoint ]
                        else:
                            #print rep
                            assert 'B' in r
                            center = [ lower_left_ha[ha][0] - eg_sep, midpoint ]

                        pairing_svg_cmds += svg_basic.enrichment_glyph_cmds( center, arrow_length, arrow_width,
                                                                             enrich )


                ystart += count*ypixel_scale



    pairing_svg_y_offset += 2*yspacer + len(tcrs)*ypixel_scale


if no_pairing_text:
    tmpcmds = pairing_svg_cmds[:]
    pairing_svg_cmds = []
    for cmd in tmpcmds:
        if '<text' in cmd:
            print 'skip:',cmd
        else:
            pairing_svg_cmds.append( cmd )


## make svg file
svgfile = '{}_vj_pairings.svg'.format( outfile_prefix)
print 'making',svgfile
bg_color = None if paper_figs else 'white'
svg_basic.create_file( pairing_svg_cmds, pairing_svg_width, pairing_svg_y_offset+bottom_margin,
                       svgfile, create_png = True, background_color = bg_color )

#exit() #hacking

if paper_supp:
    exit()


util.readme(svgfile[:-3]+'png',"""These diagrams depict the gene-segment pairing structure of the datasets. The four
genes are arrayed left to right with the alphas on the left and the betas on the right. Below each gene-type label (eg "VA")
is a color-stack showing all the TCR clones and how they break down into the different genes for that gene-type. Each clone
is devoted a constant vertical height in pixels indicated in the text at the top (N pixels in "Nx y-pixel scale"). The curved
segments joining neighboring gene-stacks show how the two gene distributions pair up, with the thickness of the segments
corresponding to the number of clones having those two segments (scaled by the indicated y-pixel scale). Significant gene-gene
pairings (positive or negative correlations with a P-value less than 1e-6) are labeled at the beginning and ending of the
corresponding segments. Gene-gene pairings which are not observed and for which this under-representation is significant
are indicated by dashed segments with P-value labels. Enrichments (depletions) of gene segments relative to
background are shown for all labeled genes by up (down) arrows where the number of arrowheads reflects the base-2
logarithm of the fold change, rounded down (one arrowhead means 2 <= fold change < 4,
two arrowheads means 4 <= fold change < 8, and so on).
<br>
<br>
The left-right ordering of the segment types is chosen so that VA and JA are on the left, VB and JB are on the right,
and the alpha-beta pairing with the largest adjusted mutual information is in the middle.
""")

## make some imshow plots

make_png = True
import matplotlib
if make_png: matplotlib.use('Agg')

import matplotlib.pyplot as plt


#epitopes = epitope_entropies.keys()
#epitopes.sort()

num_epitopes =len(epitopes)


######################################################################################
## from here below we are using a different epitopes order....
##
##
## first let's get a sensible epitope order: compute kl divergence between gene segment frequency distributions

epitope_divergences = np.zeros( (len(epitopes),len(epitopes)))

for segtype in segtypes:

    for i,ep1 in enumerate( epitopes ):
        icounts = epitope_repcounts[ep1][segtype]
        itot = sum( icounts.values() )
        for j,ep2 in enumerate( epitopes ):
            if j<=i: continue
            jcounts = epitope_repcounts[ep2][segtype]
            jtot = sum( jcounts.values() )
            js_div = 0.0
            for k in set( icounts.keys() + jcounts.keys() ):
                p = float( icounts.get(k,0) ) / itot
                q = float( jcounts.get(k,0) ) / jtot
                m = 0.5 * ( p + q )
                if p: js_div += 0.5 * p * math.log(p/m,2)
                if q: js_div += 0.5 * q * math.log(q/m,2)
            epitope_divergences[i,j] += js_div
            epitope_divergences[j,i] += js_div


for i,ep1 in enumerate( epitopes ):
    for j,ep2 in enumerate( epitopes ):
        if j<=i:continue
        print 'epitope_divergences: {:9.3f} {} {}'.format( epitope_divergences[i,j], ep1, ep2 )


if len(epitopes)>1:
    ## let's use scipy/matplotlib
    from scipy.cluster import hierarchy
    from scipy.spatial import distance

    y = distance.squareform( epitope_divergences, checks=True )

    assert len(y) == ( len(epitopes)*(len(epitopes)-1) )/2
    Z = hierarchy.average( y )
    c,coph_dists = hierarchy.cophenet(Z,y)

    leaves = hierarchy.leaves_list( Z )
    print 'old epitopes:',epitopes
    print 'leaves:',leaves
    epitopes = [ epitopes[x] for x in leaves ]
    print 'new epitopes:',epitopes
    print 'coph:',c


######################################################################################
######################################################################################

## make bar charts of cdr3 length split by rep

ncols = len(segtypes)

nrows = num_epitopes

top_margin_inches = 0.5
bottom_margin_inches = 0.25
plot_height_inches = 2.0 * nrows

fig_height = top_margin_inches + plot_height_inches + bottom_margin_inches
fig_width = 2.0 * ncols + 0.75

top_margin = float( plot_height_inches + bottom_margin_inches ) / fig_height
bottom_margin = float( bottom_margin_inches ) / fig_height


plt.figure(0,figsize=(fig_width,fig_height))

plotno=0
for epitope in epitopes:
    for segtype in segtypes:
        all_counts = epitope_repcounts[epitope][segtype]
        len_counts = epitope_repcounts_by_len[epitope][segtype]

        plotno += 1
        plt.subplot(nrows,ncols,plotno)

        all_l = [ (y,x) for x,y in all_counts.iteritems() ]
        all_l.sort()
        all_l.reverse()

        all_total = sum((x[0] for x in all_l))

        reps_sorted = [x[1] for x in all_l]

        rep_colors = dict( zip( reps_sorted, html_colors.get_rank_colors_no_lights(len(reps_sorted))) )

        min_count_for_label = 0.05 * all_total

        rep_labels = dict( [ (x[1], x[1][2:] if x[0] >= min_count_for_label else '' ) for x in all_l ] )

        total = sum(all_counts.values())

        lefts=[]
        heights=[]
        bottoms=[]
        colors = []
        #labels = []
        bar_reps = []
        #widths = []

        for le in range(min_cdr3len,max_cdr3len+1):
            ## total bar height
            counts = len_counts.get(le,{})
            total_this_len = sum( counts.values() )
            frac = float( total_this_len )/total

            l = [ (y,x) for x,y in counts.iteritems() ]
            l.sort() # smallest to largest
            #l.reverse()

            height_total=0
            for count,rep in l:
                height = float(count)/total
                lefts.append( le-0.4 )
                bottoms.append( height_total )
                colors.append( rep_colors[rep] )
                # if rep_labels[rep] not in labels:
                #     labels.append( rep_labels[rep])
                # else:
                #     labels.append( '' )
                heights.append( height )
                bar_reps.append( rep )
                height_total += height
        bars = plt.bar( lefts, heights, width=0.8, bottom=bottoms, color=colors,edgecolor='none' )
        assert len(bars) == len(bar_reps)
        legend_bars = [ bars[ bar_reps.index(x) ] for x in reps_sorted if rep_labels[x] ]
        legend_labels = [           rep_labels[x] for x in reps_sorted if rep_labels[x] ]

        #plt.legend()
        plt.legend( legend_bars, legend_labels, fontsize=7, markerscale=0.5, loc='upper left', handlelength=1,
                    frameon=False )

        locs,labels = plt.yticks()
        newlocs = []
        newlabels = []
        for loc in locs:
            if int(100*loc)%10==0:
                newlocs.append( loc )
                num = int(100*loc)/10
                newlabels.append( '.{}'.format(num))

        plt.yticks(newlocs,newlabels)


        if epitope==epitopes[0]:
            plt.title('{}'.format(segtype))#,segtype))


        if segtype==segtypes[-1]:
            x=0.15+plt.xlim()[1]
            y=sum(plt.ylim())/2.
            plt.text(x,y,epitope,fontdict={'fontsize':14})#horizontalalignment='right',verticalalignment='center')
            #plt.text(x,y,segtype,fontdict={'fontsize':14},rotation=270)#horizontalalignment='right',verticalalignment='center')


        # if epitope==epitopes[0]:
        #     angle = math.pi
        #     x=1.5*radius*math.cos( angle )
        #     y=1.5*radius*math.sin( angle )
        #     plt.text(x,y,segtype,fontdict={'fontsize':8},horizontalalignment='right',verticalalignment='center')

        # if segtype == segtypes[-1]:
        #     angle = 3*math.pi/2
        #     x=1.5*radius*math.cos( angle )
        #     y=1.5*radius*math.sin( angle )
        #     plt.text(x,y,epitope,fontdict={'fontsize':8},horizontalalignment='center',verticalalignment='top')


plt.subplots_adjust(left=0.05,right=0.9,bottom=bottom_margin,top=top_margin )

pngfile = '{}_cdr3lens.png'.format(outfile_prefix)
print 'making',pngfile
plt.savefig(pngfile)
util.readme(pngfile,"""These bar plots show the cdr3-length distributions for each epitope, colored by
gene segment. Each epitope is a single row. The left two columns show CDR3-alpha length distributions,
colored by V-segment in the first column and J-segment in the second column. The right two columns depict the
CDR3-beta length distributions in the same manner. Segments comprising at least 5 percent of the epitope's
dataset are labeled.
""")


#plt.show()
##


######################################################################################
######################################################################################

## dynamic figure sizing!
nrows = len(segtypes)
ncols = len(epitopes)

preferred_plot_width = 12.0
preferred_plot_height = 12.0

preferred_cell_size = max( 2.0, min( preferred_plot_height/nrows, preferred_plot_width/ncols ) )

plot_width = ncols * preferred_cell_size
plot_height = nrows * preferred_cell_size

fontsize_labels = 8.
fontsize_names = 12.

for repeat in range(3):
    if plot_width <= 1.2 * preferred_plot_width and plot_height <= 1.2 * preferred_plot_height: break

    if plot_width / preferred_plot_width > plot_height / preferred_plot_height: ## too wide
        plot_width *= 0.75
        plot_height *= 0.9
        fontsize_labels *= 0.9

    else: ## too tall
        plot_height *= 0.75
        plot_width *= 0.9
        fontsize_labels *= 0.9

fontsize_labels = max(5,int(floor(0.5+fontsize_labels)))

epitope_labels = dict( zip( epitopes, ( '{} ({})'.format(x,len(all_tcrs[x])) for x in epitopes ) ) )


fudge = 1.2
bottom_spacer = 0.3 # inches
left_margin_inches   = 0.25
bottom_tree_height_inches = 2.5
bottom_margin_inches = fudge * max( ( len(epitope_labels[x]) for x in epitopes ) ) * \
                       0.75 * fontsize_names / 72.0 + bottom_spacer + bottom_tree_height_inches

## allow for labels hanging off
top_margin_inches = 0.5
right_margin_inches = 0.5

fig_width  =   left_margin_inches + plot_width  + right_margin_inches
fig_height = bottom_margin_inches + plot_height +   top_margin_inches

top_margin    = float( bottom_margin_inches + plot_height ) / fig_height
bottom_margin = float( bottom_margin_inches ) / fig_height
left_margin   = float( left_margin_inches ) / fig_width
right_margin  = float( left_margin_inches + plot_width ) / fig_width

bottom_tree_height_fraction = float( bottom_tree_height_inches ) / fig_height


fig = plt.figure(3,figsize=(fig_width,fig_height))

## make pie charts for all the epitopes:  _gene_segment_pies.png


plotno=0
repslist = []
for segtype in segtypes:
    for epitope in epitopes:
        counts = epitope_repcounts[epitope][segtype]
        plotno += 1
        plt.subplot(nrows,ncols,plotno)

        l = [ (y,x) for x,y in counts.iteritems() ]
        l.sort()
        l.reverse()

        reps_sorted = [x[1] for x in l]
        for tempr in reps_sorted:
            if tempr not in repslist:
                repslist.append(tempr)
        if consistentfigcolors:
            rep_colors = dict( zip( repslist, html_colors.get_rank_colors_no_lights(len(repslist)))) ##This will keep colors the same across pies
        else:
            rep_colors = dict( zip( reps_sorted, html_colors.get_rank_colors_no_lights(len(reps_sorted))) )

        total = sum(counts.values())

        wedges, texts = plt.pie( [x[0] for x in l] )
        for (w,(count,rep)) in zip(wedges,l):
            w.set_facecolor(rep_colors[rep])
            w.set_edgecolor('none')
            #assert abs(w.r-1)<1e-3
            ## maybe add a label
            frac = float(count)/total
            if frac > min_gene_frequency_for_labels:##
                thresh = 0.3*w.r
                angle = math.pi * 0.5 * ( w.theta1 + w.theta2 ) / 180.0
                x=1.05*w.r*math.cos( angle )
                y=1.05*w.r*math.sin( angle )
                ha = 'left'   if x>thresh else ( 'center' if x>-thresh else 'right' )
                va = 'bottom' if y>thresh else ( 'center' if y>-thresh else 'top' )
                plt.text(x,y,'{}'.format(rep[2:]),color=rep_colors[rep],
                         fontdict={'fontsize':fontsize_labels},horizontalalignment=ha,verticalalignment=va)

        radius = wedges[0].r
        # if False and epitope==epitopes[0]:
        #     angle = math.pi
        #     x=1.5*radius*math.cos( angle )
        #     y=1.5*radius*math.sin( angle )
        #     plt.text(x,y,segtype,fontdict={'fontsize':8},horizontalalignment='right',verticalalignment='center')

        if segtype == segtypes[-1]:
            angle = 3*math.pi/2
            x=1.5*radius*math.cos( angle )
            y=1.5*radius*math.sin( angle )
            plt.text(x,y,epitope_labels[epitope],fontdict={'fontsize':fontsize_names},rotation='vertical',
                     horizontalalignment='center',verticalalignment='top')


plt.subplots_adjust(left=left_margin,right=right_margin,bottom=bottom_margin,top=top_margin )

if len(epitopes)>1:
    ## now let's add in the tree figure if we can
    ax = fig.add_axes( [ left_margin, 0.0, right_margin - left_margin, bottom_tree_height_fraction] )
    hierarchy.dendrogram( Z, ax=ax, orientation='bottom' )
    ax.axis('off')

pngfile = '{}_gene_segment_pies.png'.format(outfile_prefix)
print 'making',pngfile
plt.savefig(pngfile)
util.readme(pngfile,"""These pie charts depict the gene segment composition of the epitope-specific datasets,
with each column corresponding to an epitope (labeled at the bottom, with total number of clones given in
parentheses). Genes making up at least 5 percent of the population are labeled. The rows correspond to
V-alpha, J-alpha, V-beta, and J-beta distributions. At the bottom, the epitopes are clustered based on these gene-frequencies
using the Jensen-Shannon divergence measure.
""")

#exit()
##

#print 'HACK!','{}_gene_segment_pies.png'.format(outfile_prefix)
#exit()

######################################################################################

plt.figure(1,figsize=(8,12))

entropy_keys_single = [ 'VA','JA','VB','JB' ]
#entropy_keys_single = [ 'VA','JA','CA','VB','JB','CB' ]

A = np.zeros( ( len(epitopes), len(entropy_keys_single) ) )

for i,e in enumerate(epitopes):
    for j,k in enumerate(entropy_keys_single):
        A[i,j] = epitope_entropies[e][k]

A = A.transpose()

plt.subplot(311)

aspect = 'auto'

plt.imshow( A, aspect = aspect, interpolation='nearest',
            vmin=min_entropy_for_colorscale, vmax=max_entropy_for_colorscale )

plt.xticks( range(len(epitopes)), epitopes, rotation='vertical' )
plt.yticks( range(len(entropy_keys_single)), entropy_keys_single )

plt.title('gene entropies (colorscale: {:.2f}-{:.2f})'\
          .format(min_entropy_for_colorscale,max_entropy_for_colorscale))





###################################################################################################
A = np.zeros( ( len(epitopes), len(entropy_keys_single) ) )

for i,e in enumerate(epitopes):
    for j,k in enumerate(entropy_keys_single):
        A[i,j] = epitope_jsds[e][k]

A = A.transpose()

plt.subplot(312)

aspect = 'auto'

if max_jsd_for_colorscale==0.0: #not set on cmdline
    min_jsd_for_colorscale = 0.0
    max_jsd_for_colorscale = np.amax(A)

plt.imshow( A, aspect = aspect, interpolation='nearest',
            vmin=min_jsd_for_colorscale, vmax=max_jsd_for_colorscale )

plt.xticks( range(len(epitopes)), epitopes, rotation='vertical' )
plt.yticks( range(len(entropy_keys_single)), entropy_keys_single )

plt.title('gene J-S divergence to background (colorscale: {:.3f}-{:.3f})'\
          .format(min_jsd_for_colorscale,max_jsd_for_colorscale))




###################################################################################################

entropy_keys_double = [] ## this has to agree with the setup of repcounts2 up above (and I think it does)
for i,a in enumerate(entropy_keys_single):
    for j,b in enumerate(entropy_keys_single):
        if j<=i:continue
        if a[1] != b[1]: continue
        entropy_keys_double.append( ( a,b) )
for i,a in enumerate(entropy_keys_single):
    for j,b in enumerate(entropy_keys_single):
        if j<=i:continue
        if a[1] == b[1]: continue
        entropy_keys_double.append( ( a,b) )

A = np.zeros( ( len(epitopes), len(entropy_keys_double) ) )

for i,e in enumerate(epitopes):
    for j,k in enumerate(entropy_keys_double):
        A[i,j] = epitope_mis[e][k][1] ## adjusted_mi

A = A.transpose()

plt.subplot(313)

plt.imshow( A, aspect = aspect, interpolation='nearest',
            vmin = min_ami_for_colorscale, vmax=max_ami_for_colorscale )

plt.xticks( range(len(epitopes)), epitopes, rotation='vertical' )
plt.yticks( range(len(entropy_keys_double)), ['{}-{}'.format(x[0],x[1]) for x in entropy_keys_double ] )

plt.title('gene-gene adjusted mutual information (colorscale: {:.3f}-{:.3f})'\
          .format(min_ami_for_colorscale,max_ami_for_colorscale))
plt.subplots_adjust( hspace=0.35, top=0.95 )
pngfile = '{}_gene_entropies_and_mi.png'.format(outfile_prefix)
print 'making',pngfile
plt.savefig(pngfile)
if paper_figs: plt.savefig(pngfile+'.svg')
util.readme(pngfile,"""These three plots contain information on the gene-segment distributions and
correlations. <br><br>
The top plot (gene entropies) shows the total Shannon entropy of the gene segment
probability distributions for each epitope (columns) and the four gene types (rows: VA,JA,VB,JB). Red (higher
entropy) means more disorder (total frequency divided more evenly among larger numbers of genes) while
blue (less entropy) reflects stronger preference for one or a few genes. Segments like V-alpha with a lot of
possible genes will tend to have higher entropies than segments like J-beta with only a few.
<br><br>
The middle plot (gene relative entropies) shows the difference in gene frequency distributions between the observed, epitope-specific
repertoires and control or background distributions derived from large, non-epitope-specific datasets. The number being plotted
is the Jensen-Shannon divergence between the background and observed distributions, divided by the mean Shannon entropy
of the background and observed distributions.
<br><br>
The bottom plot (gene-gene adjusted mutual information) shows the mutual information between the different
gene segment distributions, with red indicating that knowledge of the identity of one gene (for example V-alpha) gives
more information about the identity of another (e.g., V-beta), and blue suggesting that the corresponding gene-pair
for that row are approximately independent for the epitope corresponding to that column. The actual number ranges
from 0 (completely independent) to 1 (completely dependent) and is an adjusted form of the mutual information that
takes into account the numbers of counts (computed using the scikit-learn routine sklearn.metrics.adjusted_mutual_info_score)
""")


#############################################################################################3

plt.figure(2,figsize=(15,12))

## scatter plot of all

colors = 'rgbcmyk'

for e,color in zip(epitopes,colors[:len(epitopes)] ):
    xvals = [x[0] for x in epitope_correlations[e]]
    yvals = [x[1] for x in epitope_correlations[e]]

    sizes = [ max(1,min(40,x[1])) for x in epitope_correlations[e]]

    plt.scatter(xvals, yvals, s=sizes,c=color,edgecolors='none')

    plotted = []

    min_text_distance = 1.0

    sigl = [ ( x[1]/5 + abs(x[0] ), x) for x in epitope_correlations[e] ]
    sigl.sort()
    sigl.reverse()

    for sig, (log2enrich,neglog10prob,rs,rep1,rep2) in sigl:
        x=log2enrich
        y=neglog10prob/4
        if neglog10prob>8 and sig > 4:
            tooclose = False
            for u,v in plotted:
                if sqrt((x-u)**2+(y-v)**2) < min_text_distance:
                    tooclose=True
                    break
            if not tooclose:
                plotted.append( ( x,y) )
                plt.text( log2enrich, neglog10prob, '{} {} {}'\
                          .format( e,
                                   rep1[2:rep1.index('*')] if '*' in rep1 else rep1[2:],
                                   rep2[2:rep2.index('*')] if '*' in rep2 else rep2[2:] ),
                          fontsize=8)

xmn,xmx = plt.xlim()
ymn,ymx = plt.ylim()
plt.plot( [0,0], [ymn,ymx],c='gray' )
plt.plot( [xmn,xmx], [0,0],c='gray' )
plt.xlim( (xmn,xmx) )
plt.ylim( (ymn,ymx) )
plt.xlabel('log2enrich (0-counts set to 0.25)')
plt.ylabel('-1*log10(P-value)')

plt.title('alpha-beta gene correlations')
plt.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95 )

pngfile = '{}_gene_gene_correlations.png'.format(outfile_prefix)
print 'making',pngfile
plt.savefig(pngfile)
## this is not currently being used
# util.readme(pngfile,"""
# """)
