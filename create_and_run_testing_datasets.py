from basic import *
import time

if not isdir('testing/'):
    mkdir('testing/')

cmd = 'rm testing/* testing/*web/*'
print cmd
system(cmd)



#examples = [ [ 'test_tiny', 1, 1, 15 ],
#             [ 'test_small', 3, 3, 15 ] ]

examples = [ [ 'test_small', 3, 3, 15 ] ]

#examples = [ [ 'test_tiny', 1, 1, 15 ] ]

for filetag, max_epitopes, max_subjects, max_tcrs_per_subject in examples:

    fields = 'id	epitope	subject	a_nucseq	b_nucseq	a_quals	b_quals'.split()

    for organism in ['mouse','human']:
        ## this is the dataset from the paper
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


        ## use intrasubject_nbrdists here because these mini-repertoires may contain only a single subject
        ## if there's only one subject, then we can't compute a nbrdist if we exclude intra-subject distances
        cmd = 'python run_basic_analysis.py --constant_seed --intrasubject_nbrdists --organism {} --pair_seqs_file {} > {}.log 2> {}.err &'\
              .format( organism, newfile, newfile, newfile )
        print cmd
        system(cmd)
        time.sleep(1.0) ## short pause


