from basic import *
#import read_sanger_data
import tcr_sampler
import tcr_rearrangement_new
import logo_tools
from all_genes import all_genes, gap_character
from translation import get_translation

#import matplotlib
#if make_png: matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#import numpy as np

with Parser(locals()) as p:
    #p.str('args').unspecified_default().multiple().required()
    p.str('logfile')
    p.str('listfile')
    p.str('chain').required() #either A or B
    p.str('organism').required()
    p.int('min_cdr3_len').default(3+1+2) ## 'CAX' + 1 + 'XF'
    p.int('min_v_score').default(10)
    p.int('min_j_score').default(8)
    p.int('mismatch_score_for_correcting_cdr3_seqs').default(-6)
    p.flag('dump_out_of_frame')       #
    p.flag('woof')                    # show OOF info for out-of-frames
    p.flag('correct_cdr3_seqs')       #
    p.multiword('logfiles').cast(lambda x:x.split())

if not logfiles:
    logfiles = []
if logfile:
    logfiles.append( logfile )
if listfile:
    logfiles.extend( [x[:-1] for x in open(listfile,'r') ] )

assert chain in 'AB'

bases = 'acgt'

allowed_v_genes = frozenset( ( id for id,g in all_genes[organism].iteritems() if g.chain == chain and g.region=='V' ))
allowed_j_genes = frozenset( ( id for id,g in all_genes[organism].iteritems() if g.chain == chain and g.region=='J' ))


for logfile in logfiles:
    for line in open( logfile,'r'):
        if 'v_offset:' not in line: continue
        if 'UNK' in line: continue
        #print line
        l = line.split()
        assert l[-2] == 'filename:' # new format
        assert l[-4] == 'nucleic:'
        fastq_file = l[-1]
        nucseq = l[-3]
        assert l[-10] == 'tsvfile_line_number:'
        seqid = 'L'+l[-9]
        vpos = l.index('vseg:')
        jpos = l.index('jseg:')
        v_gene = l[vpos+1]
        j_gene = l[jpos+1]
        v_score,v_offset,v_0start,v_0stop,v_mismatches,v_reversed = map(int,[l[vpos+3+2*x] for x in range(6) ] )
        j_score,j_offset,j_0start,j_0stop,j_mismatches,j_reversed = map(int,[l[jpos+3+2*x] for x in range(6) ] )
        vpos2 = l.index('all_v_hits:')
        jpos2 = l.index('all_j_hits:')
        all_v_genes = l[vpos2+1].split(',')
        all_j_genes = l[jpos2+1].split(',')
        if v_score < min_v_score or j_score < min_j_score: continue
        if v_reversed != j_reversed: continue

        if v_0stop >= j_0start: continue

        if v_gene not in allowed_v_genes or j_gene not in allowed_j_genes:
            continue

        ## if we parsed an adaptive TCRAD file we may have alpha j-regions, for example
        all_v_genes = [ x for x in all_v_genes if x in allowed_v_genes ]
        all_j_genes = [ x for x in all_j_genes if x in allowed_j_genes ]

        if v_reversed:
            nucseq = logo_tools.reverse_complement(nucseq)

        # get the cdr3 sequence
        vg = all_genes[organism][v_gene]
        jg = all_genes[organism][j_gene]
        v_protseq = vg.protseq
        j_protseq = jg.protseq
        v_nucseq  = vg.nucseq
        j_nucseq  = jg.nucseq

        v_nucseq_offset = vg.nucseq_offset
        j_nucseq_offset = jg.nucseq_offset

        assert v_nucseq_offset in range(3)
        assert j_nucseq_offset in range(3)

        ## get the full matched sequence
        ## extend or trim on either end so that it goes from 'C' to 'F' before GXG
        ## are the frames of V match and J match in agreement?
        v_alseq = vg.alseq
        alseq_cpos = vg.cdr_columns[-1][0] - 1 ## 0-indexed
        numgaps = v_alseq[:alseq_cpos].count('.')
        v_cpos_protseq = alseq_cpos - numgaps
        v_cpos_nucseq = 3*v_cpos_protseq + v_nucseq_offset ## this is where the nucleotide cdr3 should start

        num_genome_j_aas_in_loop = len(jg.cdrs[0].replace(gap_character,'')) ## to GXG
        j_gpos_nucseq = 3*num_genome_j_aas_in_loop + j_nucseq_offset ## start of the 'G' codon

        #v_match_begin = v_0start + v_offset
        v_match_end = v_0stop + v_offset

        #j_match_begin = v_0start + j_offset
        #j_match_end = v_0stop + j_offset

        if v_match_end < v_cpos_nucseq:
            ## want at least a little overlap with the V cdr3 stretch
            #Log('early V match')
            continue

        ## define the cdr3 nucseq
        q_cpos_nucseq = v_cpos_nucseq - v_offset
        q_gpos_nucseq = j_gpos_nucseq - j_offset

        if q_cpos_nucseq >= q_gpos_nucseq:
            #Log('wrong order')
            continue

        prefix = ''
        suffix = ''
        if q_cpos_nucseq <0:
            ## need to add some junk from the V nucseq
            num = -1*q_cpos_nucseq
            prefix = v_nucseq[ v_cpos_protseq : v_cpos_protseq + num ]
        if q_gpos_nucseq >len(nucseq):
            num = q_gpos_nucseq - len(nucseq)
            suffix = j_nucseq[ j_gpos_nucseq-num:j_gpos_nucseq]

        cdr3_nucseq = prefix + nucseq[ max(0,q_cpos_nucseq) : min(len(nucseq),q_gpos_nucseq) ] + suffix

        v_rep = vg.rep
        j_rep = jg.rep

        v_reps = set( [ all_genes[organism][x].rep for x in all_v_genes ] )
        j_reps = set( [ all_genes[organism][x].rep for x in all_j_genes ] )

        if dump_out_of_frame:
            if len(cdr3_nucseq)%3:
                print nucseq
            continue

        if len(cdr3_nucseq)%3:
            if woof:
                print 'OOF {} {} {} {} {:d} {:d} {} {} {} {}:{}:{}'\
                    .format( v_gene, v_rep, j_gene, j_rep,
                             v_score, j_score,
                             ','.join(all_v_genes),
                             ','.join(all_j_genes),
                             cdr3_nucseq,
                             logfile, fastq_file, seqid
                    )
            continue
        if len(cdr3_nucseq)/3 < min_cdr3_len:
            continue

        ## in frame
        cdr3_protseq, codons = get_translation( cdr3_nucseq, '+1' )

        if '*' in cdr3_protseq or 'X' in cdr3_protseq:
            continue

        original_cdr3_nucseq = cdr3_nucseq[:]
        original_cdr3_protseq = cdr3_protseq[:]

        if chain == 'A':
            if correct_cdr3_seqs:
                tmp_results = tcr_sampler.analyze_junction\
                              ( organism, v_gene, j_gene, cdr3_protseq, cdr3_nucseq,
                                return_corrected_cdr3_seqs = True,
                                mismatch_score = mismatch_score_for_correcting_cdr3_seqs )
                corrected_cdr3_nucseq, corrected_cdr3_protseq = list(tmp_results)[-2:]

                if corrected_cdr3_nucseq != cdr3_nucseq:
                    print 'fixing early sequence error',cdr3_nucseq,'==>',corrected_cdr3_nucseq
                    cdr3_nucseq = corrected_cdr3_nucseq[:]
                    cdr3_protseq = corrected_cdr3_protseq[:]

            junction_results = tcr_sampler.analyze_junction( organism, v_gene, j_gene, cdr3_protseq, cdr3_nucseq )
            new_nucseq, cdr3_protseq_masked, cdr3_protseq_new_nucleotide_countstring, trims, inserts \
                = junction_results

            ( v_trim, d0_trim, d1_trim, j_trim ) = trims
            ( d_id, n_vd_insert, n_dj_insert, n_vj_insert ) = inserts
            if not new_nucseq: new_nucseq = '-'

            assert d0_trim + d1_trim + n_vd_insert + n_dj_insert == 0
            junction_info = '%s %s -%d -%d +%d'%( new_nucseq, cdr3_protseq_masked, v_trim, j_trim, n_vj_insert )
        else:
            junction_infos = []

            possible_d_ids = tcr_rearrangement_new.all_trbd_nucseq[organism].keys()[:] ; possible_d_ids.sort()
            for force_d_id in possible_d_ids:
                ## fix sequences if necessary ########################
                if correct_cdr3_seqs and force_d_id == possible_d_ids[0]: ## just the first time
                    tmp_results = tcr_sampler.analyze_junction\
                                  ( organism, v_gene, j_gene, cdr3_protseq, cdr3_nucseq,
                                    force_d_id = force_d_id,
                                    return_corrected_cdr3_seqs = True,
                                    mismatch_score = mismatch_score_for_correcting_cdr3_seqs )
                    corrected_cdr3_nucseq, corrected_cdr3_protseq = list(tmp_results)[-2:]

                    if corrected_cdr3_nucseq != cdr3_nucseq:
                        print 'fixing early sequence error',cdr3_nucseq,'==>',corrected_cdr3_nucseq
                        cdr3_nucseq = corrected_cdr3_nucseq[:]
                        cdr3_protseq = corrected_cdr3_protseq[:]

                ## now analyze
                junction_results = tcr_sampler.analyze_junction( organism, v_gene, j_gene, cdr3_protseq, cdr3_nucseq,
                                                                 force_d_id = force_d_id )
                new_nucseq, cdr3_protseq_masked, cdr3_protseq_new_nucleotide_countstring, trims, inserts\
                    = junction_results

                ( v_trim, d0_trim, d1_trim, j_trim ) = trims
                ( d_id, n_vd_insert, n_dj_insert, n_vj_insert ) = inserts
                if not new_nucseq: new_nucseq = '-'

                assert n_vj_insert == 0
                if d_id:
                    assert d_id == force_d_id
                    num_d_nucleotides = len( tcr_rearrangement_new.all_trbd_nucseq[organism][d_id] ) - d0_trim - d1_trim
                else:
                    num_d_nucleotides = 0
                junction_info = '%s %s D%d %d -%d -%d -%d -%d +%d +%d'\
                                %( new_nucseq, cdr3_protseq_masked, d_id, num_d_nucleotides, v_trim,
                                   d0_trim, d1_trim, j_trim,
                                   n_vd_insert, n_dj_insert )
                junction_infos.append( ( num_d_nucleotides, junction_info ) )
            junction_infos.sort()
            junction_infos.reverse() ## put the one with the largest num_d_nucleotides in the front
            junction_info = ' '.join([x[1] for x in junction_infos])


        original_info = '{},{},{}'.format('yesmut' if original_cdr3_nucseq != cdr3_nucseq else 'nomut',
                                          original_cdr3_nucseq,
                                          original_cdr3_protseq) if correct_cdr3_seqs else ''

        print 'GENES %s %s %d %s %s %d %s %d %d %d %d %s %s %s %s %s %s:%s:%s'\
            %( v_gene, v_rep, len(v_reps),
               j_gene, j_rep, len(j_reps),
               cdr3_protseq, v_score, j_score, v_0start, j_0stop, junction_info,
               original_info,
               ','.join(all_v_genes),
               ','.join(all_j_genes),
               cdr3_nucseq,
               logfile, fastq_file, seqid
            )
        #except:
        #    print 'ERR bad line:',line





