from basic import *

## what gets checked to decide if we re-run certain steps?
##
## _*.dist
## _random_nbrdists.tsv
## _cdr3_motifs_{epitope}.log
## _tree_AB_*png
##


MAX_MOTIFS_TIME_IN_SECONDS = 24 * 60 * 60 ## a day

with Parser(locals()) as p:
    p.str('pair_seqs_file').described_as('Name of the pair_seqs file (input option #1)')
    p.str('parsed_seqs_file').described_as('Name of the parsed_seqs file (input option #2)')
    p.str('clones_file').described_as('Name of the clones file (input option #3)')
    p.flag('only_parsed_seqs')
    p.flag('find_cdr3_motifs_in_parallel').described_as('Will spawn multiple simultaneous CDR3 motif finding jobs')
    p.flag('only_clones')
    p.flag('constant_seed')
    p.str('extra_make_really_tall_trees_args').shorthand('mrtt_args')
    p.str('extra_find_clones_args').shorthand('fc_args')
    p.str('organism').required()
    p.str('webdir').described_as('Location where the index.html summary output file will be generated. Default is <clones_file>_web/')
    p.str('distance_params')
    p.int('min_quality_for_singletons').default(20).described_as('Minimum CDR3 region quality score for singleton clones')
    p.float('seed_threshold_for_motifs')
    p.flag('make_fake_quals').described_as("(for --pair_seqs_file I/O, arg is passed to read_pair_seqs.py) Create a fictitious quality string for each nucleotide sequence")
    p.flag('make_fake_ids').described_as("(for --pair_seqs_file I/O, arg passed to read_pair_seqs.py) Create an id for each line based on index in the pair_seqs file")
    p.flag('make_fake_alpha').described_as("Create a fictitious alpha chain sequence")
    p.flag('make_fake_beta').described_as("Create a fictitious beta chain sequence")
    p.flag('force')
    p.flag('borderline_motifs').described_as("Option to relax the CDR3 motif finding significance thresholds; may be helpful for small datasets")
    p.flag('webstatus')
    p.flag('dry_run')
    p.flag('intrasubject_nbrdists').described_as('Include TCRs from the same subject when computing the nbrdist (aka NNdistance) score')
    p.flag('consistentfigcolors')
    p.flag('no_probabilities').described_as('Assign a probability of 1 to all TCRs.')
    p.set_help_prefix("""

################
INPUT
################

This script will run a pipeline of analysis tools starting from three possible input filetypes:

    pair_seqs_file: A .tsv (tab-separated values) file with info on alpha and beta chain
       sequences and quality scores, epitope, subject, and id.
       Required fields: id epitope subject a_nucseq b_nucseq a_quals b_quals
       For further details, run "python read_pair_seqs.py -h"
       In particular,
          * if you don't have quality score info you can add --make_fake_quals and fictitious
            quality scores will be created.
          * if you don't have ids, you can add --make_fake_ids and default ids will be created based
            on position in the file.


    parsed_seqs_file: A processed sequence file with V and J genes assigned, CDR3s parsed, etc. Will be produced
       as an intermediate when you start from a pair_seqs_file and could be then reused to re-run downstream
       analyses. Could also be generated from the output of other tools (eg MIXCR, conversion scripts yet to come).

    clones_file: A processed version of the parsed_seqs_file in which per-TCR probabilities have been assigned and
       clones identified. One line per clone. Most of the scripts in the pipeline take a clones_file as input.

    Use the corresponding command line option to point to the type of file you have.


#################
OUTPUT
#################

    Running the pipeline will generate a slew of different files. The best place to start will be the html output
    which can be found in the file

    <clones_file>_web/index.html

    where <clones_file> is the name of the clones_file produced by or input to the pipeline. To change the directory
    where index.html will live, use the --webdir command line option.


#################
NOTES
#################

    The pipeline has been tested on datasets as large as 5000 sequences, but it does start to get slow.
    Eventually we would like to code some of the time-intensive steps in C/C++/other. If this would help
    you, let me know!

    CDR3 motif finding in particular is slow for larger datasets. If you have multiple cores available,
    consider using the option:

       --find_cdr3_motifs_in_parallel

    Also, be aware of the --min_quality_for_singletons flag which can cause some (singleton) TCR clones
    to be filtered out if they have bad sequence read quality scores. The default is 20.


#################
SUPPORT/FEEDBACK
#################

    This analysis pipeline is a work in progress. Please direct questions/suggestions to

    pbradley@fredhutch.org

    Thank you!


""")

## imports are slow if we only wanted the --help output, so do these now
try:
    import numpy
except:
    print '[ERROR] Failed to import the python module scipy-- is it installed? I really need it.'
    exit()

try:
    import scipy
except:
    print '[ERROR] Failed to import the python module scipy-- is it installed? I really need it.'
    exit()

try:
    import matplotlib
except:
    print '[ERROR] Failed to import the python module matplotlib-- is it installed? I really need it.'
    exit()

try:
    import sklearn
except:
    print """
=============================================================================
[ERROR] failed to import the python module sklearn (scikit-learn)
[ERROR] Some analyses (kernelPCA plots, adjusted_mutual_information) will fail
[ERROR] Take a look at http://scikit-learn.org/stable/install.html
=============================================================================
"""

import time
import os
import subprocess
import sys
from paths import path_to_scripts, path_to_tablesorter_files


if pair_seqs_file:
    assert not parsed_seqs_file
    if not os.path.isfile(pair_seqs_file):
        print "Error: file " + pair_seqs_file + " does not exist."
        sys.exit()
    else:
        #checkinput(pair_seqs_file)
        #pair_seqs_file = pair_seqs_file[:-4]+ "_cleanedinput.tsv"
        parsed_seqs_file = pair_seqs_file[:-4]+'_parsed_seqs.tsv'

if parsed_seqs_file:
    assert not clones_file
    if not pair_seqs_file:
        if not os.path.isfile(parsed_seqs_file):
            print "Error: file " + parsed_seqs_file + " does not exist."
            sys.exit()
        #else:
            #checkinput(parsed_seqs_file)
            #parsed_seqs_file = parsed_seqs_file[:-4]+ "_cleanedinput.tsv"
    probs_file = parsed_seqs_file[:-4]+'_probs.tsv'
    clones_file = '{}_mq{}_clones.tsv'.format( probs_file[:-4], min_quality_for_singletons )

if distance_params:
    distance_params_args = ' --distance_params {} '.format( distance_params )
else:
    distance_params_args = ' '

if constant_seed:
    constant_seed_args = ' --constant_seed '
else:
    constant_seed_args = ' '

if not webdir:
    webdir = '{}_web'.format(clones_file[:-4])

if webdir.endswith('/'):webdir = webdir[:-1]

if not exists(webdir) and not only_parsed_seqs and not only_clones:
    os.mkdir(webdir)


if not only_parsed_seqs and not only_clones:
    files = glob(path_to_tablesorter_files+'/*')
    for file in files:
        #print 'copying tablesorter file:',file
        system('cp {} {}'.format( file, webdir ) )

webfile = '{}/index.html'.format(webdir)

if not ( only_clones or only_parsed_seqs ):
    print '\nWill generate summary output file: {}\n'.format(webfile)

webdir_contains_input_files = ( os.path.dirname(os.path.normpath(os.path.realpath( webfile ))) ==
                                os.path.dirname(os.path.normpath(os.path.realpath( clones_file ))) )

if webstatus:
    out = open(webfile,'w')
    out.write("""<!doctype html>
<title>Running</title>
<h1>Analysis is in progress, sorry for the delay, try reloading from time to time</h1>
""")
    out.close()

all_logfiles = []
all_errfiles = []

def run(cmd):
    if not dry_run:
        print cmd
        if webstatus: ## we want a continuously updating index.html file
            outwebstatus = open(webfile,'a')
            outwebstatus.write('<h2>Running:</h2>\n{}<br><br>\n'.format(cmd))
            outwebstatus.close()
        system(cmd)
        cmdl = cmd.split()
        if len(cmdl)>=4 and cmdl[-4] == '>' and cmdl[-2] == '2>':
            logfile = cmdl[-3]
            errfile = cmdl[-1]
            all_logfiles.append( logfile )
            all_errfiles.append( errfile )
            if webstatus:
                errlines = '<br>'.join( popen('tail '+errfile).readlines() )
                outwebstatus = open(webfile,'a')
                outwebstatus.write('<i>Last few stderr lines from run:</i><br><br>{}<br><br>\n'\
                                   .format(errlines))
                outwebstatus.close()



if pair_seqs_file and ( force or not exists( parsed_seqs_file ) ):
    cmd = 'python {}/read_pair_seqs.py {} {} {} {} --organism {} --infile {} --outfile {} -c > {}.log 2> {}.err'\
          .format( path_to_scripts,
                   ' --make_fake_ids ' if make_fake_ids else '',
                   ' --make_fake_quals ' if make_fake_quals else '',
                   ' --make_fake_alpha ' if make_fake_alpha else '',
                   ' --make_fake_beta ' if make_fake_beta else '',
                   organism, pair_seqs_file, parsed_seqs_file, parsed_seqs_file, parsed_seqs_file )
    run( cmd )

if only_parsed_seqs:
    exit()

if parsed_seqs_file and ( force or not exists( clones_file ) ):
    if no_probabilities:
        noprobsarg = "--no_probabilities"
    else:
        noprobsarg = " "
    ## compute probs
    if force or not exists( probs_file ):
        cmd = 'python {}/compute_probs.py --organism {}  --infile {} --outfile {} {}  -c --filter --add_masked_seqs > {}.log 2> {}.err'\
              .format( path_to_scripts, organism, parsed_seqs_file, probs_file, noprobsarg, probs_file, probs_file )
        run(cmd)

    ## find the clones
    if force or not exists( clones_file ) or extra_find_clones_args:
        cmd = 'python {}/find_clones.py {} --organism {}  --infile {} --outfile {}  -c --min_quality_for_singletons {} > {}.log 2> {}.err'\
            .format( path_to_scripts, extra_find_clones_args if extra_find_clones_args else ' ',
                     organism, probs_file, clones_file, min_quality_for_singletons, clones_file, clones_file )
        run(cmd)

assert exists( clones_file )

if only_clones:
    exit()

all_clones = parse_tsv_file( clones_file, ['epitope','subject'], ['cdr3a'], False )

epitopes = all_clones.keys()[:]
epitopes.sort()


## make a mouse table

cmd = 'python {}/make_mouse_table.py --clones_file {} > {}_mmt.log 2> {}_mmt.err'\
      .format( path_to_scripts, clones_file, clones_file, clones_file )
run(cmd)


## precompute some info on gene frequencies
cmd = 'python {}/analyze_gene_frequencies.py --organism {}  --clones_file {} > {}_agf.log 2> {}_agf.err'\
      .format( path_to_scripts, organism, clones_file, clones_file, clones_file )
run(cmd)


## make gene plots (entropy, relentropy, ami, covariation, pie charts of gene usage) and VJ pairings
cmd = 'python {}/make_gene_plots.py {} --organism {}  --clones_file {} --use_color_gradients > {}_mgp.log 2> {}_mgp.err'\
    .format( path_to_scripts, ' --consistentfigcolors '*consistentfigcolors,
             organism, clones_file, clones_file, clones_file )
run(cmd)


## compute distances
distfiles = glob('{}_*.dist'.format(clones_file[:-4]))

cmd = 'python {}/compute_distances.py {} {} --organism {} --clones_file {} > {}_cd.log 2> {}_cd.err'\
    .format( path_to_scripts, distance_params_args, ' --intrasubject_nbrdists '*intrasubject_nbrdists,
             organism, clones_file, clones_file, clones_file )
if force or not distfiles: run(cmd)


## plot nbr-distance histograms
cmd = 'python {}/plot_nbrdist_distributions.py --clones_file {} --nbrdist_percentiles 5 10 25 > {}_pnd.log 2> {}_pnd.err'\
    .format( path_to_scripts, clones_file, clones_file, clones_file )
run(cmd)

## compare with random tcrs
random_nbrdists_file = '{}_random_nbrdists.tsv'.format(clones_file[:-4] )
if not exists( random_nbrdists_file ):
    cmd = 'python {}/random_tcr_distances.py {} {} --organism {} --clones_file {} > {}_rtd.log 2> {}_rtd.err'\
          .format( path_to_scripts, constant_seed_args, distance_params_args, organism,
                   clones_file, clones_file, clones_file )
    run(cmd)

## now read the output of the random nbrdists
#assert exists( random_nbrdists_file ) ## tmp hacking

cmd = 'python {}/read_random_tcr_distances.py --organism {} --clones_file {} > {}_rrtd.log 2> {}_rrtd.err'\
      .format( path_to_scripts, organism, clones_file, clones_file, clones_file )
run(cmd)


## analyze overlap
#logfile = '{}_sharing.log'.format(clones_file[:-4])
if force or True: #not exists( logfile ):
    cmd = 'python {}/analyze_overlap_compute_simpsons.py --organism {} --clones_file {} > {}_aocs.log 2> {}_aocs.err'\
          .format( path_to_scripts, organism, clones_file, clones_file, clones_file )
    run(cmd)

## make overlap plot
cmd = 'python {}/plot_sharing.py --organism {} --clones_file {} > {}_ps.log 2> {}_ps.err'\
    .format( path_to_scripts, organism, clones_file, clones_file, clones_file )
run(cmd)


## make tall trees
cmd = 'python {}/make_tall_trees.py {} --organism {} --clones_file {} --junction_bars > {}_mtt.log 2> {}_mtt.err'\
      .format( path_to_scripts, constant_seed_args, organism, clones_file, clones_file, clones_file )
run(cmd)

## analyze intra-subject privacy
cmd = 'python {}/analyze_epitope_privacy.py {} {} --organism {} --clones_file {} --all_chains AB --nrepeat 1000 --tree_height_inches 5.0 --nbrdist_percentile 10 > {}_aep.log 2> {}_aep.err'\
      .format( path_to_scripts, constant_seed_args, distance_params_args, organism,
               clones_file, clones_file, clones_file )
run(cmd)

## find motifs #################################################
#motifs_files = glob('{}_cdr3_motifs_*log'.format(clones_file[:-4] ) )

max_ng_lines = 5000000 ## mouse beta has ~2 million; both human have way more

if 1: #force or not motifs_files:
    nsamples = 25
    max_motif_len = 100
    #min_count=10
    fixlen = False
    nofilter = False


    all_procs = {} ## only used if we are finding motifs in parallel
    for ep in epitopes:
        outfile = '{}_cdr3_motifs_{}.log'.format( clones_file[:-4], ep )
        errfile = '{}_cdr3_motifs_{}.err'.format( clones_file[:-4], ep )

        if exists( outfile ) and not force: continue ############### NOTE NOTE NOTE

        num_clones = sum( ( len(all_clones[ep][x]) for x in all_clones[ep].keys() ) )
        if seed_threshold_for_motifs:
            my_seed_threshold_for_motifs = seed_threshold_for_motifs
        elif num_clones<200:
            my_seed_threshold_for_motifs = 10
        else:
            my_seed_threshold_for_motifs = None ## use the default

        min_min_count = 5  ## otherwise we will have trouble making the chi-squared threshold of 75.0 in read_motifs for the top motif
        max_min_count = 10 ## even for larger sets... is this OK?

        min_count = min( max_min_count, max( min_min_count, num_clones/10 ) )

        print 'num_clones:',ep,num_clones,'my_seed_threshold_for_motifs:',my_seed_threshold_for_motifs,'min_count:',min_count

        if borderline_motifs:
            min_count = 3
            my_seed_threshold_for_motifs = 5.

        extra_args = ' --chi_squared_threshold_for_seeds {} '.format( my_seed_threshold_for_motifs ) \
                     if my_seed_threshold_for_motifs else ''


        cmd = 'python {}/find_cdr3_motifs.py {} --organism {} {} --clones_file {} --min_count {} --epitopes {} {} {} --verbose --big --nsamples {} --max_motif_len {} --max_ng_lines {} > {} 2> {}'\
            .format( path_to_scripts, constant_seed_args, organism, extra_args, clones_file, min_count, ep,
                     ' --nofilter '*nofilter,
                     ' --force_random_len '*fixlen,
                     nsamples, max_motif_len, max_ng_lines,
                     outfile, errfile )

        if find_cdr3_motifs_in_parallel:
            print 'SPAWN: cmd'
            proc = subprocess.Popen( cmd, shell=True )
            all_procs[ep] = proc
        else:
            run(cmd)

    if find_cdr3_motifs_in_parallel:
        total_time = 0 ; sleepseconds = 10
        while total_time < MAX_MOTIFS_TIME_IN_SECONDS:
            time.sleep(sleepseconds)
            total_time += sleepseconds
            all_done = True
            for ep,proc in all_procs.iteritems():
                retval = proc.poll()
                print total_time, ep, retval
                if retval==None: all_done = False

            if all_done: break

        for ep,proc in all_procs.iteritems():
            retval = proc.poll()
            if retval==None:
                ## didn't finish!!!
                print 'ACK motif finding didnt finish!!!',ep
                proc.kill()


## make motifs summary

extra_args = ' --min_chi_squared 30 --min_top_chi_squared 30 ' if borderline_motifs else ''

cmd = 'python {}/read_motifs.py {} {} --junction_bars --max_ng_lines {} --organism {} --clones_file {} > {}_rm.log 2> {}_rm.err'\
      .format( path_to_scripts, extra_args, constant_seed_args, max_ng_lines, organism,
               clones_file, clones_file, clones_file )
run(cmd)

## make kpca landscape plots
cmd = 'python {}/make_kpca_plots.py --organism {} --clones_file {} --showmotifs > {}_kpca.log 2> {}_kpca.err'\
    .format( path_to_scripts, organism, clones_file, clones_file, clones_file )
run(cmd)


## make a summary table
cmd = 'python {}/make_summary_table.py --clones_file {} > {}_mst.log 2> {}_mst.err'\
      .format( path_to_scripts, clones_file, clones_file, clones_file )
run(cmd)

summary_table_file = clones_file[:-4]+'_summary_table.html'
cdr3_table_file = clones_file[:-4]+'_CDR3_table.html'
assert exists( summary_table_file ) and exists( cdr3_table_file )


## make a bunch of trees
tree_files = glob('{}_tree_AB_*png'.format(clones_file[:-4]))
if force or not tree_files or extra_make_really_tall_trees_args:
    cmd = 'python {}/make_really_tall_trees.py {} {} --organism {} --clones_file {} > {}_mrtt.log 2> {}_mrtt.err'\
        .format( path_to_scripts, constant_seed_args,
                 extra_make_really_tall_trees_args if extra_make_really_tall_trees_args else ' ',
                 organism, clones_file, clones_file, clones_file )
    run(cmd)


## now make some more pages, one for each color scheme a/b/ab combination, and maybe also one for each epitope
all_tree_files = {}
color_schemes = []

for epitope in epitopes:
    tree_files = glob('{}_tree_[AB]_{}_*png'.format(clones_file[:-4],epitope)) +  glob('{}_tree_AB_{}_*png'.format(clones_file[:-4],epitope))

    ## actually lets glob by epitope... easier for epitopes w underscores...
    for file in tree_files:
        suffix = file[ len(clones_file)-4:-4] ## also trim .png
        print 'tree_file:',file
        print 'suffix:',suffix
        assert suffix.startswith('_tree_')
        suffix = suffix[6:]
        ab = suffix.split('_')[0]
        suffix = suffix[len(ab)+1:]
        assert suffix.startswith(epitope)
        myfile = True
        for other_ep in epitopes:
            if other_ep != epitope and other_ep.startswith(epitope) and suffix.startswith(other_ep):
                myfile = False
        if not myfile:
            continue
        suffix = suffix[len(epitope)+1:]
        color_scheme = suffix
        if color_scheme not in color_schemes: color_schemes.append( color_scheme )
        all_tree_files[ (ab,epitope,color_scheme)] = file

color_schemes.sort()

## make some tree webpages



tree_pages = []
subject_tree_suffixes = []
for epitope in epitopes:
    subject_tree_suffixes.append( '_{}_subject_tree.png'.format( epitope ) )
    for ab in ['A','B','AB']:
        pagename = '{} trees {}'.format( epitope, ab )
        pages = {}
        for frm in ['png','svg']:
            treefile = '{}/trees_{}_{}_{}.html'.format(webdir,epitope,ab,frm)
            pages[frm] = treefile

            out = open(treefile,'w')

            out.write("""<!doctype html>
            <title>{} {}</title>
            <div class="page">
            """.format(clones_file.split('/')[-1][:6], pagename ) )

            for cs in color_schemes:
                file = all_tree_files.get( (ab,epitope,cs), None )
                if file:
                    if frm == 'png':
                        out.write('{}\n<br>\n'.format(file))
                        if not webdir_contains_input_files: run('cp {} {}'.format(file,webdir))
                        out.write('<img src="{}" />\n'.format( file.split('/')[-1] ) )
                        out.write('<br>\n')
                    else:
                        assert frm == 'svg'
                        svgfile = file[:-3]+'svg'
                        out.write('{}\n<br>\n'.format(svgfile))
                        out.writelines(open(svgfile,'r').readlines())
                        out.write('<br>\n')
            out.write('</div>\n')
            out.close()
        tree_pages.append( [ pagename, pages ] )

for cs in color_schemes:
    for ab in ['A','B','AB']:
        pagename = '{} trees {}'.format( cs,ab)
        pages = {}
        for frm in ['png','svg']:
            treefile = '{}/trees_{}_{}_{}.html'.format(webdir,cs,ab,frm)
            pages[frm] = treefile

            out = open(treefile,'w')

            out.write("""<!doctype html>
            <title>{} {}</title>
            <div class="page">
            """.format(clones_file.split('/')[-1][:6], pagename ))

            for epitope in epitopes:
                file = all_tree_files.get( (ab,epitope,cs), None )
                if file:
                    if frm == 'png':
                        out.write('{}\n<br>\n'.format(file))
                        if not webdir_contains_input_files: run('cp {} {}'.format(file,webdir))
                        out.write('<img src="{}" />\n'.format( file.split('/')[-1] ) )
                    else:
                        assert frm == 'svg'
                        svgfile = file[:-3]+'svg'
                        out.write('{}\n<br>\n'.format(svgfile))
                        out.writelines(open(svgfile,'r').readlines())
                        out.write('<br>\n')
                out.write('<br>\n')
            out.write('</div>\n')
            out.close()
        tree_pages.append( [ pagename, pages ] )


# _random_nbrdists_nbrdists_AB.dens.png
# _random_nbrdists_roc_AB.dens.png
# _random_nbrdists_summary.dens.png
# _mouse_nbrdist_rank_score_heterogeneity.png


pngfile_suffixes = """
_cdr3_distributions.png
_kpca.png
_sharing_diversity.png
_subject_heterogeneity.png
_subject_table.png
_gene_segment_pies.png
_cdr3lens.png
_gene_entropies_and_mi.png
_nbrdist10_distributions_w_cdf.png
_nbrdist25_distributions_w_cdf.png
_random_nbrdists_nbrdists_AB.png
_random_nbrdists_roc_AB.png
_random_nbrdists_summary.png
_random_nbrdists_nbrdist_roc_superpositions.png
_sharing_and_clonality_wtd_nbrdist10.png
_sharing_and_clonality_by_epitope_wtd_nbrdist10.png
_epitope_correlations_10.png
_epitope_epitope_avg_nbrdist_rank_scores.png
_epitope_distances.png
_sharing.png
_subject_trees.png
_vj_pairings.png
_motif_summary_ngbig.png
_tall_tree_AB.png
_tall_tree_A.png,_tall_tree_B.png
""".split('\n') + subject_tree_suffixes

# _random_nbrdists_nbrdists_A.png
# _random_nbrdists_nbrdists_B.png
# _random_nbrdists_roc_A.png
# _random_nbrdists_roc_B.png
#_gene_gene_correlations.png

if webstatus:
    system('cp {} {}.old.html'.format(webfile,webfile))


out = open(webfile,'w')

table_style = """
<style>
table {
    font-family: arial, sans-serif;
    border-collapse: collapse;
    width: 100%;
}

td, th {
    border: 1px solid #dddddd;
    text-align: center;
    padding: 8px;
}

tr:nth-child(even) {
    background-color: #dddddd;
}


table.tablesorter {
        font-family:arial;
        background-color: #CDCDCD;
        margin:10px 0pt 15px;
        font-size: 10pt;
        width: 100%;
        text-align: left;
}

table.tablesorter thead tr th, table.tablesorter tfoot tr th {
        background-color: #e6EEEE;
        border: 1px solid #FFF;
        font-size: 10pt;
        padding: 4px;
}

table.tablesorter thead tr .header {
        background-image: url(bg.gif);
        background-repeat: no-repeat;
        background-position: center right;
        cursor: pointer;
}
table.tablesorter tbody td {
        color: #3D3D3D;
        padding: 4px;
        background-color: #FFF;
        vertical-align: top;
}
table.tablesorter tbody tr.odd td {
        background-color:#F0F0F6;
}
table.tablesorter thead tr .headerSortUp {
        background-image: url(asc.gif);
}
table.tablesorter thead tr .headerSortDown {
        background-image: url(desc.gif);
}
table.tablesorter thead tr .headerSortDown, table.tablesorter thead tr .headerSortUp {
background-color: #8dbdd8;
}


</style>
"""

tablesorter_js_includes = """
<script type="text/javascript" src="jquery-latest.js"></script>
<script type="text/javascript" src="jquery.tablesorter.js"></script>
"""

## the parser for scientific notation below came from an answer to the question:
## http://stackoverflow.com/questions/4126206/javascript-parsefloat-1-23e-7-gives-1-23e-7-when-need-0-000000123
##

tablesorter_js_code = """
<script>
	$(document).ready(function()
  {
  $("#summaryTable").tablesorter();
  $("#cdr3Table").tablesorter();
  }
	);

	// add parser through the tablesorter addParser method
	$.tablesorter.addParser({
	// set a unique id
	id: 'scinot',
	is: function(s) {
  return /[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?/.test(s);
	},
	format: function(s) {
  return $.tablesorter.formatFloat(s);
	},
	type: 'numeric'
	});

</script>

"""

out.write("""<!doctype html>
<html>
<head>
<title>{}</title>
{}
{}
</head>
<body>

<h1>Overview:</h1>

This page contains the results of TCR repertoire analyses conducted with the tcr-dist pipeline.
For a general overview of the
analysis methods, please consult the reference(s) mentioned in the "CITING" section of the README.txt file in the
tcr-dist repository. Individual plots should be preceded by some explanatory text below.

<br>
<br>

<h1>Command:</h1>
{}
<br>

<h1>Summary table: (click on headers to sort)</h1>
{}
<br>

<h1>CDR3 table: (click on headers to sort)</h1>
{}
<br>

{}

""".format( clones_file.split('/')[-1][:-4], table_style, tablesorter_js_includes,
            ' '.join(sys.argv),
            ''.join( open( summary_table_file,'r').readlines()),
            ''.join( open( cdr3_table_file,'r').readlines()),
            tablesorter_js_code ) )

if webstatus:
    out.write('<h1>Files:</h1>\n')
    if pair_seqs_file:
        ## maybe copy
        if not webdir_contains_input_files:
            system('cp {} {}'.format(pair_seqs_file, webdir ) )
        out.write('<a href="{}">pair_seqs file</a><br>\n'\
                  .format( os.path.basename(pair_seqs_file) ))

    if parsed_seqs_file:
        ## maybe copy
        if not webdir_contains_input_files:
            system('cp {} {}'.format(parsed_seqs_file, webdir ) )
        out.write('<a href="{}">parsed_seqs file</a><br>\n'\
                  .format( os.path.basename(parsed_seqs_file) ))

    ## definitely want clones file
    if not webdir_contains_input_files:
        system('cp {} {}'.format(clones_file, webdir ) )
    out.write('<a href="{}">clones file</a><br>\n'\
              .format( os.path.basename(clones_file) ))

out.write("""
<h1>TCR Clustering trees</h1>
The following links point to pages containing detailed TCRdist hierarchical clustering trees,
grouped by epitope (the first set) or coloring scheme (the second set).
<br>
<br>
"clonality" trees are colored by the size of each TCR clonotype (number of members).
<br>
<br>
"cross_reactivity" trees are colored by occurrence of TCRs in other epitope datasets.
<br>
<br>
"min_other_nbrdist" trees are colored by the lowest (percentiled) nbrdist score for an epitope
other than the source epitope for each receptor, to get an idea of which other repertoires a
given receptor is similar to.
<br>
<br>
"probs" trees are colored by TCR generation probability according to a very simple model of the
rearrangement process.
<br>
<br>
"sharing" trees are colored by the number of subjects a given TCR occurs in.
<br>
<br>
The trees are annotated to the left with the subject (truncated but prefixed with a unique integer), CDR3 amino acid
sequence, non-germline CDR3 amino acids, number of insertions and deletions, clone size, number of subjects
in which the TCR occurs (same repertoire, all repertoires), number of repertoires in which the TCR occurs,
gene families, and the numerical value of the coloring parameter (and any associated string-tag where relevant) if
not already present in the preceding info.
<br>
<br>
<br>
""")

for pagename, pages in tree_pages:
    out.write('{} <a href="{}">(png)</a> <a href="{}">(svg)</a><br>\n'\
              .format( pagename, pages['png'].split('/')[-1], pages['svg'].split('/')[-1] ))

for sufs in pngfile_suffixes:
    if not sufs:continue
    files= [ clones_file[:-4] + x for x in sufs.split(',') ]
    out.write('<br><br><h1>{}</h1>\n'.format(sufs))
    for file in files:
        readme_file = file+'.readme'
        if exists( readme_file):
            readme_text = ''.join( open(readme_file,'r').readlines())
            out.write('{}<br>\n'.format(readme_text))
        if not exists(file):
            svgfile = file[:-4]+'.svg'
            if exists( svgfile ): ## use svg version instead
                out.write("<br><br>CONVERSION TO .png FAILED, USING .svg VERSION!!!<br>\n")
                out.writelines( open( svgfile,'r').readlines())
            else:
                print 'missing:',file
                out.write("<br><br><br>The image file {} is missing. If there is just a single subject, then some of the subject_tree and subject_heterogeneity analyses don't pertain. Or if this is the motifs summary it may be that there were no motifs found. Try grepping for 'Error' in the files: <clones_file>*.err\n".format(file))
        else:
            ## can't use 'run' command since that messes with the webfile
            if not webdir_contains_input_files: system('cp {} {}'.format(file,webdir))
            out.write('<img src="{}" />\n'.format( file.split('/')[-1] ) )
    out.write('<br>\n')

out.write("""
</body>
</html>
""")

out.close()

