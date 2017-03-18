from basic import *

for organism in ['mouse','human']:

    pair_seqs_file = 'datasets/{}_pairseqs_v1.tsv'.format(organism)
    assert exists( pair_seqs_file )

    cmd = 'nice python run_basic_analysis.py --organism {} --pair_seqs_file {} --find_cdr3_motifs_in_parallel > {}.log 2> {}.err &'\
        .format( organism, pair_seqs_file, pair_seqs_file, pair_seqs_file )
    print cmd
    system(cmd)

