from basic import *

with Parser(locals()) as p:
    p.flag('from_pair_seqs').described_as('Restart from the nucleotide sequences (default is to start from parsed clones file)')


for organism in ['mouse','human']:

    if from_pairseqs:

        pair_seqs_file = 'datasets/{}_pairseqs_v1.tsv'.format(organism)
        assert exists( pair_seqs_file )

        cmd = 'nice python run_basic_analysis.py --organism {} --pair_seqs_file {} --find_cdr3_motifs_in_parallel > {}.log 2> {}.err &'\
            .format( organism, pair_seqs_file, pair_seqs_file, pair_seqs_file )
        print cmd
        system(cmd)
    else:

        clones_file = 'datasets/{}_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv'.format(organism)
        assert exists( clones_file )

        cmd = 'nice python run_basic_analysis.py --organism {} --clones_file {} --find_cdr3_motifs_in_parallel > {}.log 2> {}.err &'\
            .format( organism, clones_file, clones_file, clones_file )
        print cmd
        system(cmd)

