from basic import *
import time

with Parser(locals()) as p:
    p.flag('from_pair_seqs').described_as('Restart from the nucleotide sequences (default is to start from parsed clones file)')
    p.flag('multicore').described_as('Allow the pipeline to start multiple processes for motif finding; also run mouse and human analyses simultaneously')
    p.set_help_prefix( "\nThis script will re-run the analysis from the primary citation; it will take a few hours. Use --multicore to speed things up if you have multiple cores to burn\n" )



## These are the distance parameters used in the original publication in Nature.
##
## Since then, we've seen evidence that a gap penalty in the CDR3 region of 12 performs slightly better.
## So that is the new default; here the scale_factor term accounts for the difference in TCRdist values
## (lower overall w/ classic gap penalty of 8)
##
distance_params='gap_penalty_v_region:4,gap_penalty_cdr3_region:8,weight_cdr3_region:3,align_cdr3s:False,trim_cdr3s:True,scale_factor:1.0439137134052388'


for organism in ['mouse','human']:

    if from_pair_seqs:

        pair_seqs_file = 'datasets/{}_pairseqs_v1.tsv'.format(organism)
        assert exists( pair_seqs_file )

        extra_args = ' --find_cdr3_motifs_in_parallel ' if multicore else ''
        cmd_suffix = ' &' if multicore else ''

        cmd = 'nice python run_basic_analysis.py {} --distance_params {} --organism {} --pair_seqs_file {} > {}.log 2> {}.err {}'\
            .format( extra_args, distance_params, organism, pair_seqs_file, pair_seqs_file, pair_seqs_file, cmd_suffix )
        print cmd
        system(cmd)
        time.sleep(1.0) ## short pause
    else:

        clones_file = 'datasets/{}_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv'.format(organism)
        assert exists( clones_file )

        extra_args = ' --find_cdr3_motifs_in_parallel ' if multicore else ''
        cmd_suffix = ' &' if multicore else ''

        cmd = 'nice python run_basic_analysis.py {} --distance_params {} --organism {} --clones_file {} > {}.log 2> {}.err {}'\
            .format( extra_args, distance_params, organism, clones_file, clones_file, clones_file, cmd_suffix )
        print cmd
        system(cmd)
        time.sleep(1.0) ## short pause

