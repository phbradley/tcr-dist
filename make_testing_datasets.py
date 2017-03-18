from basic import *

for filetag, max_epitopes, max_subjects, max_tcrs_per_subject in \
    [ [ 'test_tiny', 1, 1, 15 ],
      [ 'test_small', 3, 3, 15 ] ]:

    fields = 'id	epitope	subject	a_nucseq	b_nucseq	a_quals	b_quals'.split()

    for organism in ['mouse','human']:
        oldfile = './datasets/{}_pairseqs_v1.tsv'.format(organism)
        assert exists( oldfile )

        all_tcrs = parse_tsv_file( oldfile, ['epitope','subject'], [], True )

        epitopes = sorted( all_tcrs.keys() )[:max_epitopes]
        newfile = './testing/{}_{}_pairseqs_v1.tsv'.format(filetag, organism)

        print 'making',newfile

        out = open(newfile,'w')
        out.write( '\t'.join( fields )+'\n' )

        for epitope in epitopes:
            epitope_tcrs = all_tcrs[epitope]
            subjects = sorted( epitope_tcrs.keys() )[:max_subjects]
            for subject in subjects:
                tcrs = epitope_tcrs[subject]
                for outl in tcrs[:max_tcrs_per_subject]:
                    out.write(make_tsv_line( outl, fields )+'\n' )
        out.close()

        cmd = 'python run_basic_analysis.py --organism {} --pair_seqs_file {} > {}.log 2> {}.err &'\
              .format( organism, newfile, newfile, newfile )
        print cmd
        system(cmd)


