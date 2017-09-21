from basic import *
from all_genes import all_genes


with Parser(locals()) as p:
    p.int('maxlines').default(10000000)
    p.flag('show_vj_counts')

dbfile = pipeline_params['db_file']

if 'alphabeta' in dbfile:

    all_logfiles = {
        'mouse':
        { 'A':
          [ '/home/pbradley/csdat/tcr/paul_thomas/raw_data/sra/SRX1074287_SRX1074288_TCR_LA_MC_PCR/tmp.read_sra_matches.alpha.log',
            '/home/pbradley/csdat/tcr/paul_thomas/scjobs_tcrseq_tmp3.read_sra_matches.log' ],
          'B':
          [ '/home/pbradley/csjobs/tcrseq/tmp1.read_sra_matches.log',
            '/home/pbradley/csdat/tcr/paul_thomas/raw_data/sra/SRX1074287_SRX1074288_TCR_LA_MC_PCR/tmp.read_sra_matches.beta.log',
            '/home/pbradley/csdat/tcr/paul_thomas/raw_data/sra/friedman_SRP015131_TCR_beta_repertoire/tmp.read_sra_matches.log' ],
      },
        'human':
        { 'B':
          [ '/home/pbradley/scjobs/tcrseq/tmp7.read_sra_matches.log.subjectX_gdna_data_TCRB',
            '/home/pbradley/scjobs/tcrseq/tmp7.read_sra_matches.log.subjectY_gdna_data_TCRB',
            '/home/pbradley/scjobs/tcrseq/tmp7.read_sra_matches.log.experiment5_gdna_data_BreastPBMC2_TCRB',
            '/home/pbradley/scjobs/tcrseq/tmp7.read_sra_matches.log.experiment5_gdna_data_KidneyPBMC3_TCRB' ],
          'A':
          [ '/home/pbradley/scjobs/tcrseq/tmp8.read_sra_matches.log.subjectX_gdna_data_TCRA',
            '/home/pbradley/scjobs/tcrseq/tmp8.read_sra_matches.log.subjectY_gdna_data_TCRA',
            '/home/pbradley/scjobs/tcrseq/tmp8.read_sra_matches.log.experiment5_gdna_data_BreastPBMC2_TCRA',
            '/home/pbradley/scjobs/tcrseq/tmp8.read_sra_matches.log.experiment5_gdna_data_KidneyPBMC3_TCRA', ],
      }
    }

elif 'gammadelta' in dbfile:
    all_logfiles = {
        'human':
        { 'B':
          [ '/home/pbradley/scjobs/tcrseq/tmp57.read_sra_matches.log' ],
          'A':
          [ '/home/pbradley/scjobs/tcrseq/tmp59.read_sra_matches.log' ]
      }
    }

else: # nothing else yet
    exit(1)


for organism in all_logfiles:
    outfile = '{}/redo_nextgen_tuple_counts_v2_{}_max10M.log'.format( path_to_current_db_files(), organism )
    out = open(outfile,'w')

    for chain in all_logfiles[organism]:
        print 'making', outfile
        for logfile in all_logfiles[organism][chain]:
            assert exists(logfile)

            counter=0
            v_counts = {}
            j_counts = {}
            vj_counts = {}
            srcfile=''
            seen = set()
            for line in open(logfile,'r'):
                l = line.split()
                nucseq = l[-2]
                if nucseq in seen: continue
                seen.add( nucseq )
                srcfile = l[-1]
                counter+=1
                if counter>maxlines:break
                if counter%100000==0:
                    Log(`(srcfile,counter,organism,chain,len(v_counts.keys()),len(j_counts.keys()),logfile)`)
                try:
                    # v_genes = [ x for x in l[-4].split(',')
                    #             if all_genes[organism][x].chain == chain and all_genes[organism][x].region == 'V' ]
                    # j_genes = [ x for x in l[-3].split(',')
                    #             if all_genes[organism][x].chain == chain and all_genes[organism][x].region == 'J' ]
                    #j_genes = l[-3].split(',')
                    v_genes = tuple( sorted( set( all_genes[organism][x].count_rep for x in l[-4].split(',') ) ) )
                    j_genes = tuple( sorted( set( all_genes[organism][x].count_rep for x in l[-3].split(',') ) ) )
                except:
                    print 'bad line:',l[-4:],line
                    continue

                v_counts[v_genes] = v_counts.get(v_genes,0)+1
                j_counts[j_genes] = j_counts.get(j_genes,0)+1
                vj_genes = ( v_genes, j_genes )
                vj_counts[ vj_genes ] = vj_counts.get( vj_genes,0 )+1


            l = [ (y,x) for x,y in v_counts.iteritems() ]
            l.sort()
            l.reverse()
            for count,genes in l:
                out.write( 'TUPLE_COUNT {} {} {} {}\n'.format( count, ','.join( genes ), organism, logfile ) )

            l = [ (y,x) for x,y in j_counts.iteritems() ]
            l.sort()
            l.reverse()
            for count,genes in l:
                out.write( 'TUPLE_COUNT {} {} {} {}\n'.format( count, ','.join( genes ), organism, logfile ) )

            if show_vj_counts:
                l = [ (y,x) for x,y in vj_counts.iteritems() ]
                l.sort()
                l.reverse()
                for count,(v_genes,j_genes) in l:
                    out.write( 'TUPLE_COUNT {} {} {} {} {}\n'\
                               .format( count, ','.join( v_genes ), ','.join( j_genes ), organism,
                                        logfile ) )

    out.close()

        #exit()
