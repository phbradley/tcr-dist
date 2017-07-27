from basic import *
import parse_tsv
import scipy
from amino_acids import HP, GES, KD, aa_charge, amino_acids
from operator import add
import html_colors
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import util

with Parser(locals()) as p:
    p.str('clones_file').required()
    p.str('organism')

fake_chains = util.detect_fake_chains( clones_file )

table1_file_prefix = clones_file[:-4]+'_summary_table'
table2_file_prefix = clones_file[:-4]+'_CDR3_table'

## what are the different columns

header1 = [ ( 'epitope','s' ),
            ( 'num_individuals', 'd' ),
            ( 'num_parsed_reads', 'd' ),
            ( 'num_clones', 'd' ),
            ( 'clonality', '.3f' ),
            ( 'TCRdiv-a', '.1f' ),
            ( 'TCRdiv-b', '.1f' ),
            ( 'TCRdiv-ab', '.1f' ),
            ( 'Pshare-a', '.3E' ),
            ( 'Pshare-b', '.3E' ),
            ( 'Pshare-ab', '.3E' ),
            ( 'discrimination (AUROC)', '.3f' ),
            ( 'heterogeneity (Z)', '.1f' )
]

header2 = [ ( 'epitope','s' ),
            ( 'N', 'd' ), ## terse output for the CDR3 table
            ( 'a_len', '.1f' ),
            ( 'a_charge', '.1f' ),
            ( 'a_hydro1', '.1f' ),
            ( 'a_hydro2', '.2f' ),
            ( 'b_len', '.1f' ),
            ( 'b_charge', '.1f' ),
            ( 'b_hydro1', '.1f' ),
            ( 'b_hydro2', '.2f' ),
            ( 'ab_len', '.1f' ),
            ( 'ab_charge', '.1f' ),
            ( 'ab_hydro1', '.1f' ),
            ( 'ab_hydro2', '.2f' ),
]

header_tags = [x[0] for x in header1 + header2]

footnotes = { 'clonality': "Clonality measured using 1.0-Simpson's Diversity Index, weighted-average across subjects",
              'TCRdiv-a': "Inverse-Simpson's style measure of repertoire diversity using Gaussian-smoothed distance matching (alpha)",
              'TCRdiv-b': "Inverse-Simpson's style measure of repertoire diversity using Gaussian-smoothed distance matching (beta)",
              'TCRdiv-ab': "Inverse-Simpson's style measure of repertoire diversity using Gaussian-smoothed distance matching (both chains)",
              'Pshare-a': "Probability that a clone from one subject has the same alpha-chain AAs as a clone from another",
              'Pshare-b': "Probability that a clone from one subject has the same beta-chain AAs as a clone from another",
              'Pshare-ab': "Probability that a clone from one subject has the same full-chain AAs as a clone from another",
              'discrimination (AUROC)': "AUROC for discriminating epitope-specific TCRs from background set using NN-distance score",
              'heterogeneity (Z)': "Z-score of difference between inter- and intra-subject distances (large and positive suggests subjects may have 'different' repertoires)",
              'N': "num_clones",
}


for ab in ['a','b','ab']:
    footnotes[ab+'_len'] = "Mean CDR3-{} length".format(ab)
    footnotes[ab+'_charge'] = "Mean CDR3-{} total charge".format(ab)
    footnotes[ab+'_hydro1'] = "Mean CDR3-{} total hydrophobicity (GES scale)".format(ab)
    #footnotes[ab+'_hydro2'] = "Mean total hydrophobicity (KD scale)"
    footnotes[ab+'_hydro2'] = "Mean CDR3-{} total hydrophobicity (HP scale)".format(ab)

for tag in footnotes: ## no typos
    assert tag in header_tags

all_dats = {}
def add_dat( epitope, tag, val ):
    global all_dats
    assert tag in header_tags
    all_dats[ epitope ][ tag ] = val


## parse the clones file
all_tcrs = parse_tsv.parse_tsv_file( clones_file, ['epitope','subject'], ['cdr3a','cdr3b','clone_size'] )
epitopes = all_tcrs.keys()[:]
epitopes.sort()

def get_charge( cdr3 ):
    return sum( ( aa_charge.get(x,0.0) for x in cdr3 ) )

def get_hp1( cdr3 ):
    return sum( ( -1*GES.get(x,0.0) for x in cdr3 ) )

def get_hp2( cdr3 ):
    return sum( ( HP.get(x,0.0) for x in cdr3 ) )
    #return sum( ( KD.get(x,0.0) for x in cdr3 ) )


all_scores = {}
for epitope in epitopes:
    all_dats[epitope] = {}
    all_scores[epitope] = {}
    add_dat(epitope, 'epitope', epitope )
    mice = all_tcrs[epitope].keys()
    add_dat(epitope, 'num_individuals', len(mice) )
    tcrs = reduce( add, all_tcrs[epitope].values() )
    tcrs = [ [x[0], x[1], x[0]+x[1], int(x[2]) ] for x in tcrs ]
    add_dat(epitope, 'num_clones', len(tcrs) )
    add_dat(epitope, 'N', len(tcrs) )
    add_dat(epitope, 'num_parsed_reads', sum( ( x[3] for x in tcrs ) ) )

    for ii,ab in enumerate(['a','b','ab']):
        cdrs = [x[ii] for x in tcrs]
        lens = [len(x) for x in cdrs ]
        add_dat(epitope, '{}_len'.format(ab), get_mean_and_sdev( lens )[0] )
        add_dat(epitope, '{}_charge'.format(ab), get_mean_and_sdev( [ get_charge(x) for x in cdrs ] )[0] )
        add_dat(epitope, '{}_hydro1'.format(ab), get_mean_and_sdev( [    get_hp1(x) for x in cdrs ] )[0] )
        add_dat(epitope, '{}_hydro2'.format(ab), get_mean_and_sdev( [    get_hp2(x) for x in cdrs ] )[0] )

        all_scores[epitope]['{}_len'.format(ab)] = lens
        all_scores[epitope]['{}_charge'.format(ab)] = [ get_charge(x) for x in cdrs ]
        all_scores[epitope]['{}_hydro1'.format(ab)] = [ get_hp1(x) for x in cdrs ]
        all_scores[epitope]['{}_hydro2'.format(ab)] = [ get_hp2(x) for x in cdrs ]

## start getting these things

## read all the heterogeneity Z scores

logfile = '{}_aep.log'.format(clones_file)
assert exists( logfile )
for line in open( logfile,'r'):
    l = line.split()
    if not l:
        continue
    if l[0] == 'rep':
        epitope = l[-1]
        chains = l[-2]
        if chains == 'AB':
            add_dat( epitope ,  'heterogeneity (Z)', -1 * float( l[2] ) ) ## now higher means more heterogeneity (duh)

## read all sharing-type scores

logfile = '{}_sharing.log'.format(clones_file[:-4])
assert exists( logfile )

for line in open( logfile,'r'):
    l = line.split()
    if not l:
        continue
    elif l[0] == 'clone_diversity:':
        epitope = l[1]
        p = float( l[5] )
        add_dat( epitope ,  'clonality', 1.0 - p )
    elif line.startswith("GAUSSDIV SM1 SE1"):
        epitope = l[3]
        chains = l[5]
        div = float(l[7])
        add_dat(epitope, 'TCRdiv-'+chains.lower(), div )
        # if chains == 'AB':
        #     add_dat(epitope, 'diversity', div )
    elif line.startswith('avg_nbrdist:'):
        l = line.split()
        epitope,chains = l[1:3]
        if chains == 'AB':
            #add_dat(epitope, 'avg_nbrdist', float( l[3] ) )
            pass
    elif line.startswith('AA CM0 SM0 SE1'):
        l = line.split()
        assert l[7] == 'div:'
        assert l[4] == l[5]
        epitope = l[4]
        chains = l[6]
        div = float( l[9] )
        p_sharing = 1.0/div if div else 0.0
        tag = 'Pshare-'+chains.lower()
        add_dat(epitope, tag, p_sharing )

## fill in 0 sharing -- lines may not be getting written out
for epitope in all_dats:
    for chains in ['a','b','ab']:
        tag= 'Pshare-'+chains
        if tag not in all_dats[epitope]:
            print 'missing:',epitope,tag
            add_dat( epitope, tag, 0.0 )


## load auc random
desired_nbrdist_tag_suffix = 'wtd_nbrdist10'
#desired_nbrdist_label = 'nbrdist10p'
logfile= clones_file[:-4]+'_random_aucs.log'
assert exists( logfile )

for line in open( logfile,'r'):
    if line.startswith('auc_random '):
        l = line.split()
        auc = float( l[1] )
        chains = l[4]
        epitope = l[5]
        nbrdist_tag_suffix = l[6]
        if chains == 'AB' and nbrdist_tag_suffix == desired_nbrdist_tag_suffix:
            add_dat( epitope, 'discrimination (AUROC)', auc )


for header, table_file_prefix, table_id, caption_text in [ ( header1, table1_file_prefix, 'summaryTable', 'Summary information on the dataset' ),
                                                    ( header2, table2_file_prefix, 'cdr3Table', 'Summary information on the CDR3 loops' ) ]:


    ## now write out the table
    table_file = table_file_prefix+'.html'
    out = open( table_file, 'w' )
    out.write("""<table id="{}" style="width:100%" class="tablesorter">
    <caption>{}</caption>
    <thead>
    <tr>
    """.format( table_id, caption_text) )

    footnote_tags = []

    for tag,fmt in header:
        if tag in footnotes:
            footnote_tags.append( tag )
            out.write('<th>{}<sup>{}</sup></th>\n'.format(tag,len(footnote_tags)))
        else:
            out.write('<th>{}</th>\n'.format(tag))
    out.write("""
    </tr>
    </thead>
    """)
    if footnote_tags:
        out.write('<tfoot>\n')

        for ii,tag in enumerate(footnote_tags):
            out.write('<tr><th colspan="{}" style="text-align:left"><sup>{}</sup>{}</th></tr>\n'\
                      .format( len(header), ii+1, footnotes[tag] ) )

        out.write("</tfoot>\n")

    out.write("<tbody>\n")

    for epitope in epitopes:
        out.write('<tr>\n')
        for tag,fmt in header:
            if tag not in all_dats[epitope]:
                out.write('<td>N/A</td>\n')
            else:
                out.write('<td>{:{}}</td>\n'.format( all_dats[epitope][tag], fmt ) )
        out.write('</tr>\n')
    out.write("""
    </tbody>
    </table>
    """)

    out.close()

    ## now make a tsv version
    table_file = table_file_prefix+'.tsv'
    out = open( table_file, 'w' )

    ## header
    out.write('\t'.join( [ x[0] for x in header ] )+'\n')

    for epitope in epitopes:
        vals = []
        for tag,fmt in header:
            if tag not in all_dats[epitope]:
                val = 'N/A'
            else:
                val = '{:{}}'.format( all_dats[epitope][tag], fmt )
            vals.append( val )
        out.write('\t'.join( vals )+'\n' )
    out.close()


## make distribution plots
from scipy.stats import gaussian_kde
pngfile = clones_file[:-4]+'_cdr3_distributions.png'
util.readme( pngfile, """
Distributions of CDR3 properties for the different epitopes""" )

colors = html_colors.get_rank_colors_no_lights(len(epitopes))

scoretag_suffixes = ['_len','_charge','_hydro1','_hydro2']

nrows = 3
ncols = len(scoretag_suffixes)

plt.figure(1,figsize=(12,12))
plotno=0
for ab in ['a','b','ab']:
    if ab.upper() in fake_chains: continue
    for suf in scoretag_suffixes:
        plotno += 1
        plt.subplot(nrows,ncols,plotno)
        scoretag = ab+suf
        allvals = reduce( add, [ all_scores[x][scoretag] for x in epitopes ] )
        mn,mx = min(allvals),max(allvals)

        for epitope,color in zip( epitopes, colors ):
            vals = all_scores[epitope][scoretag]
            if 'hydro' in suf:
                density = gaussian_kde( vals )
                xs = np.linspace( mn, mx, 100 )
                ys = density(xs)

            else:
                count = {}
                for val in vals:
                    count[val] = count.get( val,0)+1
                imn = int(floor(mn+0.5))
                imx = int(floor(mx+0.5))
                xs = range(imn,imx+1)
                ys = [ float(count.get(x,0))/len(vals) for x in xs]

            plt.plot( xs, ys, c=color,label=epitope)
        rn = mx-mn
        plt.xlim( (mn-rn/10,mx+rn/10))
        plt.legend(fontsize=6,frameon=False,loc='best')
        plt.title(scoretag)

print 'making:',pngfile
plt.savefig(pngfile)



