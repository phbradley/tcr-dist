from basic import *
from amino_acids import amino_acids
import re
import logo_tools
import svg_basic
import tcr_sampler
from all_genes import all_genes



if 'alphabeta' in pipeline_params['db_file']: ### alphabeta tcrs
    all_ng_logfiles = {'mouse':{}, 'human':{}}

    ## tmp1 is from SRP004475 runs SRR085977 and SRR673913 (see ~/csjobs/tcrseq/input/job1.list)
    ## SRP abstract: In order to see the extent of sharing of TCR-beta repertoire between different individual,
    ##   three samples of CD4 CD8 thymocyte from three B6 mice are sequenced
    ## SRR085977 and SRR673913, both from SRX031459
    ##
    ## tmp3 is from SRP010815 runs SRR407172 and SRR407173
    ## SRP abstract: By sequencing tens of millions of TCRa chain transcripts from naive mouse CD8+ T cells, we observed a hugely diverse repertoire, comprising nearly all possible TRAV-TRAJ combinations. Our findings are not compatible with sequential coordinate gene recombination, but rather with a model in which contraction and DNA looping in the TCRad locus provide equal access to TRAV and TRAJ gene segments, similar to that demonstrated for IgH gene recombination Overall design: High-throughput sequencing of entire TCRa repertoire from C57Bl/6 mice
    ## citation: Genolet R, Stevenson BJ, Farinelli L, Osteras M et al. Highly diverse TCRa chain repertoire of pre-immune CD8+ T cells reveals new insights in gene recombination. EMBO J 2012 Apr 4;31(7):1666-78. PMID: 22373576
    ##

    all_ng_logfiles['mouse']['A'] = '/home/pbradley/csdat/tcr/paul_thomas/scjobs_tcrseq_tmp3.read_sra_matches.log'
    all_ng_logfiles['mouse']['B'] = '/home/pbradley/csjobs/tcrseq/tmp1.read_sra_matches.log'

    ## unpaired sequence data from adaptive, howie study
    all_ng_logfiles['human']['A'] = '/home/pbradley/csdat/tcr/paul_thomas/raw_data/howie/tmp8.read_sra_matches.log.subjectZ_gdna_data_TCRA'
    all_ng_logfiles['human']['B'] = '/home/pbradley/csdat/tcr/paul_thomas/raw_data/howie/tmp7.read_sra_matches.log.subjectZ_gdna_data_TCRB'


elif 'gammadelta' in pipeline_params['db_file']: # gamma-deltas
    ## these are sketchy since gammadelta gene usage is so tissue specific!
    all_ng_logfiles = { 'human':{}}
    all_ng_logfiles = {
        'human': {
            # from blood and bone marrow
            'A': '/home/pbradley/csdat/tcr/paul_thomas/raw_data/gammadelta/tmp59.read_sra_matches.log',
            # from blood
            'B': '/home/pbradley/csdat/tcr/paul_thomas/raw_data/gammadelta/tmp57.read_sra_matches.log'
            }
    }

# for org in all_ng_logfiles:
#     for ab in 'AB':
#         filename = all_ng_logfiles[org][ab]
#         print filename
#         assert exists(filename)


for organism in all_ng_logfiles:

    ng_logfiles = all_ng_logfiles[organism]

    for ab,ng_logfile in ng_logfiles.iteritems():
        min_v_score = 10
        min_j_score = 8
        counter=0
        seen = {}
        num_chains=0

        ## NOTE
        outfile = '{}/redo_new_nextgen_chains_{}_{}.tsv'.format( path_to_current_db_files(), organism, ab )
        print 'making:',outfile

        out = open(outfile,'w')
        out.write('\t'.join( ['v_reps','j_reps','cdr3','cdr3_nucseq'] )+'\n')

        for line in open(ng_logfile,'r'):
            counter+=1
            if not counter%1000000:Log(`counter`+' '+`num_chains`+' '+ng_logfile)

            l = line.split()
            if l[0] != 'GENES': continue
            v_gene, v_rep, num_v_reps, j_gene, j_rep, num_j_reps, cdr3, v_score, j_score = l[1:10]

            v_genes = l[-4].split(',')
            j_genes = l[-3].split(',')

            v_reps = frozenset( all_genes[organism][x].rep for x in v_genes )
            j_reps = frozenset( all_genes[organism][x].rep for x in j_genes )

            cdr3_nucseq = l[-2]   ## different line length for alpha vs beta

            if all_genes[organism][v_gene].chain != ab:
                print 'whoah funny line:',line[:-1]
                continue

            if int(v_score) <min_v_score or int(j_score) < min_j_score: continue

            if cdr3[0] != 'C' or len(cdr3)>25:
                continue

            skip_me = False
            for v_gene in v_genes:
                for j_gene in j_genes:
                    if v_gene in seen:
                        if j_gene in seen[v_gene]:
                            if cdr3_nucseq in seen[v_gene][j_gene]:
                                skip_me = True
                                break
            if skip_me:
                continue

            for v_gene in v_genes:
                for j_gene in j_genes:
                    if v_gene in seen:
                        if j_gene not in seen[v_gene]:
                            seen[v_gene][j_gene] = set()
                    else:
                        seen[v_gene] = {j_gene:set()}

                    seen[v_gene][j_gene].add( cdr3_nucseq)

            ## now add to the different places
            out.write('\t'.join( [','.join(v_reps),','.join(j_reps),cdr3,cdr3_nucseq] )+'\n')

            num_chains += 1

        Log('read {} {}-chains from {}'.format(num_chains,ab,ng_logfile))

        out.close()
