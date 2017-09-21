#from phil import *
from basic import *
import matplotlib
import numpy as np
import tcr_sampler
import tcr_rearrangement_new
#import cdr3_properties
from scipy.stats import poisson
import sys
import util
from amino_acids import amino_acids
from all_genes import all_genes


with Parser(locals()) as p:
    #p.str('args').unspecified_default().multiple().required()
    p.str('logfile').required()
    p.str('organism').required()
    p.str('chain').required().described_as("Either 'A' or 'B'")
    p.str('outfile')
    #p.int('min_v_score').default(10)
    #p.int('min_j_score').default(8)
    p.int('min_count').default(100)
    p.int('max_lines')
    p.int('subsample_lines')
    p.int('default_xmax').default(15)
    #p.int('int_arg').shorthand('i')
    #p.float('float_arg')     # --float_arg 9.6
    p.int('max_v_trim').default(20)
    p.int('max_j_trim').default(20)
    p.int('max_vj_insert').default(35)
    p.int('min_d_nucs_for_imotifs').default(5)
    p.int('imotif_len').default(6)
    p.flag('check_nucseqs')       # --flag_arg  (no argument passed)
    p.flag('skip_raw_data')       # --flag_arg  (no argument passed)
    p.flag('gene_frequencies')       # --flag_arg  (no argument passed)
    p.flag('hacking')       # --flag_arg  (no argument passed)
    #p.flag('force_good_did')       # --flag_arg  (no argument passed)
    p.flag('allow_repeats')       # --flag_arg  (no argument passed)
    p.flag('uniq_by_nucseq')       # --flag_arg  (no argument passed)
    p.flag('uniq_by_nucseq_and_tsv_file')       # --flag_arg  (no argument passed)
    p.flag('show_poisson')       # --flag_arg  (no argument passed)
    p.flag('dump_probs')       # --flag_arg  (no argument passed)
    p.flag('make_png')       # --flag_arg  (no argument passed)
    p.flag('imotifs')       # --flag_arg  (no argument passed)

force_good_did = True

if dump_probs:
    gene_frequencies= True

if make_png: matplotlib.use('Agg')
import matplotlib.pyplot as plt

if hacking:
    for mu in [0.1,0.5,1,2,4]:
        rv = poisson(mu)
        mean, var, skew, kurt = rv.stats(moments='mvsk')
        print mu,'mean, var, skew, kurt:',mean, var, skew, kurt
        xvals = np.arange(25)
        #xvals = np.arange(rv.ppf(0.01), rv.ppf(0.99 ))
        yvals = rv.pmf(xvals)
        plt.plot( xvals, yvals, label=`mu` )
    plt.legend()
    plt.show()
    exit()


possible_d_ids = tcr_rearrangement_new.all_trbd_nucseq[organism].keys()
num_d_ids = len(possible_d_ids)

bases = 'acgt'

all_new_nucseqs = {}
if imotifs:
    if chain == 'A':
        fake_did = 1
        all_new_nucseqs = {fake_did:[]}
    else:
        for did in range(1,num_d_ids+1):
            all_new_nucseqs[did] = []

v_data = {}
j_data = {}
vj_pairs = {}

all_nucseqs = {}
germline_nucseq = {}

gene_counts = {}
rep_counts = {}
countrep_counts = {}
u_gene_counts = {}
u_rep_counts = {}
gene_frequencies_total = 0
all_cdr3s = []

v_gene_list = []
v_rep_list = []
v_countrep_list = []
j_gene_list = []
j_rep_list = []
j_countrep_list = []

seen = {}
ncols = 0
dats_names = []
counter=0
num_bad = 0
num_bad_d = 0
numlines = 0
for line in open( logfile,'r'):
    if not line.startswith('GENES '): continue
    numlines += 1
    if subsample_lines and numlines%subsample_lines: continue
    counter += 1
    if counter%100000==0:
        Log('numlines: %d num_uniqs: %d num_bad_D: %d num_bad: %d'\
            %(counter,len(seen.keys()),num_bad_d,num_bad))
    l = line.split()
    if chain == 'A':
        expected_line_length = 21
    else:
        expected_line_length = 16 + 10 * num_d_ids ## 36 (mouse) or 46 (human) for example
    if len(l) != expected_line_length:
        Log('bad line: len= {} line= {}'.format(len(l),line))
        continue
    if max_lines and counter>max_lines: break
    line_info = l[-1].split(':')
    assert len(line_info)==3
    tsv_file = line_info[1]
    cdr3_protseq = l[7]
    cdr3_new_nucseq = l[12]
    cdr3_protseq_masked = l[13]
    cdr3_nucseq = l[-2]
    if uniq_by_nucseq:
        uniqer = cdr3_nucseq
    elif uniq_by_nucseq_and_tsv_file:
        uniqer = (cdr3_nucseq,tsv_file)
    else:
        uniqer = cdr3_protseq + cdr3_new_nucseq + cdr3_protseq_masked
    if (not allow_repeats) and uniqer in seen: continue
    badseq = False
    for a in cdr3_protseq:
        if a not in amino_acids:
            badseq = True
    if badseq:
        print 'skip bad cdr3_protseq:',cdr3_protseq
        continue
    seen[ uniqer ] = seen.get(uniqer,0) + 1
    v_gene = l[1]
    v_rep = l[2]
    j_gene = l[4]
    j_rep = l[5]
    assert chain == all_genes[organism][v_gene].chain
    jno = 0 # no filtering for D/J compatibility
    if chain == 'B' and j_gene[2] == 'B':
        jno = int(j_rep[4])
        assert jno in [1,2]
    all_v_genes = l[-4].split(',')
    all_j_genes = l[-3].split(',')
    assert v_gene in all_v_genes
    assert j_gene in all_j_genes
    if chain == 'A':
        v_trim = -1*int(l[14])
        j_trim = -1*int(l[15])
        vj_insert = int(l[16])
        dats = ( v_trim,j_trim,vj_insert,v_trim+j_trim,vj_insert)

        dats = ( dats, )

        if not dats_names:
            dats_names = ( 'v_trim','j_trim','vj_insert',
                           'tot_trim',
                           'tot_insert' )
        assert len(dats[0]) == len(dats_names)
        all_new_nucseqs[ fake_did ].append( cdr3_new_nucseq )

    else:
        assert chain == 'B'
        assert l[14][0] == 'D' and l[24][0] == 'D'
        if l[14] == 'D0': continue ## bad if best d is 0
        line_dats = []
        #all_num_d_nucs = [ int(l[15+10*x]) for x in range(num_d_ids) ]
        all_num_d_nucs = []
        for r in range(num_d_ids):
            start = 14+10*r
            assert l[start][0] == 'D'
            d_id = int(l[start][1])
            if not d_id: continue ## trims not well defined...
            if force_good_did:
                if organism=='mouse' and jno == 1 and d_id == 2: continue ## hardly ever see this pairing
                if organism=='human' and jno == 1 and d_id >= 2: continue ## ditto
            num_d_nucs,v_trim,d0_trim,d1_trim,j_trim,vd_insert,dj_insert = map(int,l[start+1:start+8])
            v_trim,d0_trim,d1_trim,j_trim = [-1*x for x in (v_trim,d0_trim,d1_trim,j_trim)]
            all_num_d_nucs.append( num_d_nucs )
            line_dats.append( ( v_trim, d0_trim, d1_trim, j_trim, vd_insert, dj_insert,
                                v_trim+d0_trim+d1_trim+j_trim,
                                vd_insert+dj_insert,
                                d0_trim+d1_trim,
                                d_id ) )
            # if vd_insert + dj_insert>max_vj_insert:
            #     print 'too long:',vd_insert + dj_insert, max_vj_insert

        if len(line_dats) == 0:
            #Log('bad D id '+line)
            num_bad_d += 1
            continue
        elif dump_probs and len(all_num_d_nucs)>=2 and all_num_d_nucs[0] == all_num_d_nucs[1]:
            if len(all_num_d_nucs)>=3 and all_num_d_nucs[1] == all_num_d_nucs[2]:
                dats = ( line_dats[0], line_dats[1], line_dats[2] )
            else:
                dats = ( line_dats[0], line_dats[1] )
        else:
            dats = ( line_dats[0], )

            if not dats_names:
                dats_names = ( 'v_trim','d0_trim','d1_trim','j_trim','vd_insert','dj_insert',
                               'tot_trim',
                               'tot_insert',
                               'tot_d_trim',
                               'd_id')
            assert len(dats[0]) == len(dats_names)

        if all_num_d_nucs[0] >= min_d_nucs_for_imotifs:
            did = line_dats[0][-1]
            all_new_nucseqs[did].append( cdr3_new_nucseq )


    if v_rep not in v_data: v_data[v_rep] = []
    if j_rep not in j_data: j_data[j_rep] = []
    v_data[v_rep].append( dats )
    j_data[j_rep].append( dats )

    vj = (v_rep,j_rep)
    vj_pairs[vj] = vj_pairs.get(vj,0)+1

    all_cdr3s.append( cdr3_protseq )

    if gene_frequencies:

        gene_frequencies_total += 1

        all_v_reps = set( [ all_genes[organism][x].rep for x in all_v_genes ] )
        all_j_reps = set( [ all_genes[organism][x].rep for x in all_j_genes ] )

        all_v_countreps = set( [ util.get_mm1_rep_gene_for_counting(x,organism) for x in all_v_genes ] )
        all_j_countreps = set( [ util.get_mm1_rep_gene_for_counting(x,organism) for x in all_j_genes ] )

        for gene in all_v_genes:
            gene_counts[gene] = gene_counts.get(gene,0)+1
            if gene not in v_gene_list: v_gene_list.append( gene )
        for gene in all_j_genes:
            gene_counts[gene] = gene_counts.get(gene,0)+1
            if gene not in j_gene_list: j_gene_list.append( gene )

        for rep in all_v_reps:
            rep_counts[rep] = rep_counts.get(rep,0)+1
            if rep not in v_rep_list: v_rep_list.append( rep )

        for rep in all_j_reps:
            rep_counts[rep] = rep_counts.get(rep,0)+1
            if rep not in j_rep_list: j_rep_list.append( rep )

        for countrep in all_v_countreps:
            countrep_counts[countrep] = countrep_counts.get(countrep,0)+1
            if countrep not in v_countrep_list: v_countrep_list.append( countrep )

        for countrep in all_j_countreps:
            countrep_counts[countrep] = countrep_counts.get(countrep,0)+1
            if countrep not in j_countrep_list: j_countrep_list.append( countrep )

        if len(all_v_genes)==1: u_gene_counts[v_gene] = u_gene_counts.get(v_gene,0)+1
        if len(all_v_reps )==1:  u_rep_counts[v_rep ] =  u_rep_counts.get(v_rep ,0)+1
        if len(all_j_genes)==1: u_gene_counts[j_gene] = u_gene_counts.get(j_gene,0)+1
        if len(all_j_reps )==1:  u_rep_counts[j_rep ] =  u_rep_counts.get(j_rep ,0)+1
        #print len(all_v_genes), len(all_v_reps), len(all_j_genes), len(all_j_reps)


    if check_nucseqs:
        nucseq = l[-2] ## actually cdr3_nucseq

        for gene in all_v_genes:
            if gene not in germline_nucseq:
                germline_nucseq[gene] = tcr_sampler.get_v_cdr3_nucseq( organism, gene, paranoid = True )
        for gene in all_j_genes:
            if gene not in germline_nucseq:
                germline_nucseq[gene] = tcr_sampler.get_j_cdr3_nucseq( organism, gene, paranoid = True )

        v_nucseqs = set( [germline_nucseq[x] for x in all_v_genes] )
        j_nucseqs = set( [germline_nucseq[x] for x in all_j_genes] )

        if len(v_nucseqs) == 1:
            if v_gene not in all_nucseqs:
                all_nucseqs[v_gene] = []
            minlen = len(germline_nucseq[v_gene] )
            if len(nucseq) >= minlen:
                all_nucseqs[v_gene].append( nucseq[:minlen] )

        if len(j_nucseqs) == 1:
            if j_gene not in all_nucseqs:
                all_nucseqs[j_gene] = []
            minlen = len(germline_nucseq[j_gene] )
            if len(nucseq) >= minlen:
                all_nucseqs[j_gene].append( nucseq[-1*minlen:] )


if imotifs: ## look for motifs in the insertion sequences
    def shuffle_seq( seq ):
        iseq = list(seq.replace('+',''))
        random.shuffle(iseq)
        rseq = []
        iseq_pos=0
        for a in seq:
            if a== '+':
                rseq.append(a)
            else:
                rseq.append( iseq[ iseq_pos ] )
                iseq_pos +=1
        assert iseq_pos == len(iseq)
        return ''.join( rseq )


    from collections import Counter
    mlen = imotif_len
    for did in all_new_nucseqs:
        if chain == 'A':
            dseq = '-'
            dseq_maxlen = 1
        else:
            dseq_maxlen = max( len(x) for x in tcr_rearrangement_new.all_trbd_nucseq[organism].values() )
            dseq = tcr_rearrangement_new.all_trbd_nucseq[organism][did]
        word_counts = Counter()
        rword_counts = Counter()
        for real_seq in all_new_nucseqs[did]:
            rand_seq = shuffle_seq(real_seq)
            #print 'rand_seq:',real_seq,rand_seq
            for seq,counts in [(real_seq,word_counts), (rand_seq,rword_counts)]:
                L = len(seq)
                for pos in range(L-mlen+1):
                    word = seq[pos:pos+mlen]
                    if '+' in word: continue
                    counts[word] += 1

        for tag, counts in [('real',word_counts),('rand',rword_counts)]:
            print '{:{}s} {}'.format( dseq, dseq_maxlen, tag ),
            for word,count in counts.most_common(10):
                print word,count,
            print
    exit()




def get_style_and_color( counter ):
    """
    '-'	solid line style
    '--'	dashed line style
    '-.'	dash-dot line style
    ':'	dotted line style
    '.'	point marker
    ','	pixel marker
    'o'	circle marker
    'v'	triangle_down marker
    '^'	triangle_up marker
    '<'	triangle_left marker
    '>'	triangle_right marker
    '1'	tri_down marker
    '2'	tri_up marker
    '3'	tri_left marker
    '4'	tri_right marker
    's'	square marker
    'p'	pentagon marker
    '*'	star marker
    'h'	hexagon1 marker
    'H'	hexagon2 marker
    '+'	plus marker
    'x'	x marker
    'D'	diamond marker
    'd'	thin_diamond marker
    '|'	vline marker
    '_'	hline marker
    """
    colors = 'rgbcmyk'
    styles = ['-','--','-.',':','-o','.','o','v','<','>']
    c = colors[ counter%(len(colors)) ]
    s = styles[ (counter/(len(colors)))%(len(styles)) ]
    return s+c




if dump_probs:
    ## 01/11/17 adding v2 to reflect countrep counts and timestamp
    if not outfile:
        uniqtag = 'ar' if allow_repeats else 'ubn' if uniq_by_nucseq else 'ubntsv' if uniq_by_nucseq_and_tsv_file else \
                  'ubcdr3'
        outfile = '{}.N{}.U{}.max+{}-{}-{}.{}.{}{}v2_dump_probs'\
            .format( logfile.split('/')[-1],
                     counter,
                     len(seen.keys()),
                     max_vj_insert,
                     max_v_trim,
                     max_j_trim,
                     uniqtag,
                     'ssl{}.'.format(subsample_lines) if subsample_lines else '',
                     'fgdid.' if force_good_did else '' )
    print 'making',outfile
    out = open( outfile,'w')

    out.write('#CMD: {}\n'.format( ' '.join( sys.argv ) ) )

    if gene_frequencies:
        # if allow_repeats:
        #     total = float( sum( seen.values() ) )
        # else:
        #     total = float( len(seen.keys())) ## for the possibly non-unique counts
        total = float( gene_frequencies_total )

        u_v_gene_total = float( sum( [u_gene_counts.get(x,0) for x in v_gene_list ] ))
        u_j_gene_total = float( sum( [u_gene_counts.get(x,0) for x in j_gene_list ] ))
        u_v_rep_total  = float( sum( [ u_rep_counts.get(x,0) for x in v_rep_list ] ))
        u_j_rep_total  = float( sum( [ u_rep_counts.get(x,0) for x in j_rep_list ] ))
        for rep in v_rep_list:
            freq   = rep_counts.get(rep,0) / total
            u_freq = u_rep_counts.get(rep,0) / u_v_rep_total
            out.write( '%sV_REP_FREQ: %9.3f %9.3f %s\n'%( chain, 100*freq, 100*u_freq,rep))
        for rep in j_rep_list:
            freq   = rep_counts.get(rep,0) / total
            u_freq = u_rep_counts.get(rep,0) / u_j_rep_total
            out.write( '%sJ_REP_FREQ: %9.3f %9.3f %s\n'%( chain, 100*freq, 100*u_freq,rep))
        totalfreq=0.
        for countrep in v_countrep_list:
            freq   = countrep_counts.get(countrep,0) / total
            out.write( '%sV_COUNTREP_FREQ: %9.3f %s\n'%( chain, 100*freq, countrep))
            totalfreq+=freq
        out.write( '%sV_COUNTREP_TOTALFREQ: %9.3f\n'%( chain, 100*totalfreq))
        totalfreq=0.
        for countrep in j_countrep_list:
            freq   = countrep_counts.get(countrep,0) / total
            out.write( '%sJ_COUNTREP_FREQ: %9.3f %s\n'%( chain, 100*freq, countrep))
            totalfreq+=freq
        out.write( '%sJ_COUNTREP_TOTALFREQ: %9.3f\n'%( chain, 100*totalfreq))
        for gene in v_gene_list:
            freq   = gene_counts.get(gene,0) / total
            u_freq = u_gene_counts.get(gene,0) / u_v_gene_total
            out.write( '%sV_GENE_FREQ: %9.3f %9.3f %s\n'%( chain, 100*freq, 100*u_freq,gene))
        for gene in j_gene_list:
            freq   = gene_counts.get(gene,0) / total
            u_freq = u_gene_counts.get(gene,0) / u_j_gene_total
            out.write( '%sJ_GENE_FREQ: %9.3f %9.3f %s\n'%( chain, 100*freq, 100*u_freq,gene))

        if False: # don't have cdr3_properties file in the new repository yet...
            ## compute some CDR3 distributions
            for prop in cdr3_properties.cdr3_properties:
                vals = [ cdr3_properties.get_cdr3_fval( prop, x ) for x in all_cdr3s ]
                mn,sdev = get_mean_and_sdev( vals )
                median = get_median( vals )
                counts = {}
                for val in vals:
                    ival = int( floor( 0.5+val) )
                    counts[ival] = counts.get(ival,0) + 1

                ## show distribution
                total = float( len(vals ) )
                dist = ' '.join( [ '{}: {:9.6f}'.format(x,counts.get(x,0)/total) for x in range( min(counts.keys()),
                                                                                                 max(counts.keys())+1)])

                out.write( 'cdr3_{}_distribution: N: {} mean: {:9.6f} sdev: {:9.6f} median: {:9.6f} distribution: {}\n'\
                           .format( prop, len(vals), mn, sdev, median, dist ) )


    plt.figure(1,figsize=(14,14))
    plt.suptitle('{}\n{}'.format( ' '.join(sys.argv), getcwd() ) )
    assert max_vj_insert
    assert max_v_trim
    assert max_j_trim
    numfigs = 1
    if chain == 'B':
        numfigs = 2
        plt.figure(2,figsize=(23,14))
        plt.suptitle('{}\n{}'.format( ' '.join(sys.argv), getcwd() ))
        all_dats = {}
        all_weights = {}
        skip_count = {}
        for did in possible_d_ids:
            all_dats[did] = []
            all_weights[did] = []
            skip_count[did] = 0
        for rep in v_data:
            #print rep, len(v_data[rep])
            for rdats in v_data[rep]:
                assert len(rdats) in range(1,num_d_ids+1)
                weight = 1.0/len(rdats)
                for r in range(len(rdats)):
                    dats = rdats[r]
                    v_trim,d0_trim,d1_trim,j_trim,vd_insert,dj_insert,tot_trim,tot_insert,tot_d_trim,d_id = dats
                    if tot_insert <= max_vj_insert and v_trim<=max_v_trim and j_trim<=max_j_trim:
                        all_dats[d_id].append( ( v_trim, d0_trim, d1_trim, j_trim, tot_d_trim,
                                                 vd_insert, dj_insert, tot_insert ) )
                        all_weights[d_id].append( weight )
                    else:
                        skip_count[d_id] += 1
                    #    print 'skip:',tot_insert,v_trim,j_trim
        for did in possible_d_ids:
            print 'd_id: {} num_dats: {} skip_fraction: {:.3f}%'\
                .format( did, len(all_dats[did]), float( 100*skip_count[did] ) / ( skip_count[did] + len(all_dats[did])))

        counter=-1
        for d_id in sorted( all_dats.keys()):
            ## fit a 2d model of the d-trims
            plt.figure(2)
            plt.subplot(1,num_d_ids,d_id)
            counts = {}
            for d,wt in zip( all_dats[d_id], all_weights[d_id] ):
                dd = (d[1],d[2])
                counts[dd] = counts.get(dd,0)+wt
            xvals = sorted(counts.keys())
            total = float(sum(counts.values()))
            probs = [ (counts[x]/total) for x in xvals]
            label = 'PROB_{}_D{}_d01_trim'.format( chain, d_id )
            outline = '{} {}'.format( label, ' '.join( [ '%d,%d: %9.6f'%(x[0],x[1],y) for x,y in zip(xvals,probs) ] ))
            #print outline
            out.write(outline+'\n')


            dseq = tcr_rearrangement_new.all_trbd_nucseq[organism][d_id]
            L = len(dseq)+1
            A = np.zeros( ( L,L))
            for i in range(L):
                for j in range(L):
                    A[i][j] = float( counts.get((i,j),0))/total

            A = A.transpose()

            plt.imshow( A, origin = 'lower', interpolation='nearest' )#, cmap=plt.get_cmap('bwr'),vmin=vmin, vmax=vmax )
            plt.xticks( range(L), range(L) )
            plt.yticks( range(L), range(L) )
            plt.title('D{} trims'.format(d_id))


            #maxvals = [ max_v_trim, max_j_trim, max_vj_insert ]
            plt.figure(1)
            tags = [ 'v_trim', 'd0_trim','d1_trim', 'j_trim', 'tot_d_trim', 'vd_insert','dj_insert','tot_insert' ]
            color_counter=-1
            for ii,tag in enumerate(tags):
                if 'd0' in tag or 'd1' in tag: continue
                color_counter += 1
                color = 'rgbcmky'[color_counter]
                style = ['-','--',':','-.','-o','.','o','v','<','>'][ d_id-1 ]
                counter+=1
                counts = {}
                for d,wt in zip( all_dats[d_id], all_weights[d_id] ):
                    counts[d[ii]] = counts.get( d[ii],0) + wt
                xvals = sorted(counts.keys())
                total = float(sum(counts.values()))
                probs = [ (counts[x]/total) for x in xvals]
                label = 'PROB_{}_D{}_{}'.format( chain, d_id, tag )
                plt.plot( xvals,probs,style+color,label=label)
                outline = '{} {}'.format( label, ' '.join( [ '%d: %9.6f'%(x,y) for x,y in zip(xvals,probs) ] ) )
                #print outline
                out.write(outline+'\n')

        plt.legend()

    elif chain == 'A':
        assert chain == 'A'
        all_dats = []
        for rep in v_data:
            for (( v_trim,j_trim,vj_insert,tot_trim,tot_insert),) in v_data[rep]:
                if tot_insert <= max_vj_insert and v_trim<=max_v_trim and j_trim<=max_j_trim:
                    all_dats.append( ( v_trim, j_trim, vj_insert ) )
        maxvals = [ max_v_trim, max_j_trim, max_vj_insert ]
        tags = [ 'v_trim', 'j_trim', 'vj_insert' ]
        for ii in range(3):
            mx = maxvals[ii]
            xvals = range(mx+1)
            counts = dict( zip( xvals, [0]*(mx+1) ) )
            for d in all_dats:
                counts[d[ii]] += 1
            total = float(sum(counts.values()))
            probs = [ (counts[x]/total) for x in xvals]
            label = 'PROB_{}_{}'.format( chain, tags[ii] )
            plt.plot( xvals,probs,label=label)
            outline = '{} {}'.format(  label, ' '.join( [ '%9.6f'%x for x in probs ] ) )
            #print outline
            out.write(outline+'\n')

        plt.legend()

    out.close()

    for fig in range(1,numfigs+1):
        plt.figure(fig)
        pngfile = '{}.F{}.png'.format(outfile,fig)
        print 'making',pngfile
        plt.savefig(pngfile)

    if not make_png:
        plt.show()


    exit()


if check_nucseqs:
    for vj in 'VJ':
        for gene,nucseqs in all_nucseqs.iteritems():
            if gene[3] != vj: continue
            if not nucseqs: continue

            rep = all_genes[organism][gene].rep

            germline = germline_nucseq[gene]
            L = len(germline)
            pwm = {}
            for i in range(L):
                pwm[i] = dict(zip(bases,[0]*4))
            if vj=='V':
                start=L-1 ; direction = -1
            else:
                start = 0 ; direction = 1
            for nucseq in nucseqs:
                for i,a in enumerate(nucseq):
                    if a not in bases: continue
                    pwm[start+i*direction][a] += 1
            warn = False
            expected = []
            consensus = []
            for i in range(L):
                expected_base = germline[start+i*direction]
                total = float( sum( pwm[i].values() ) )
                if not total:
                    print 'huh?',i,L,len(nucseqs)
                    continue
                l = [ (pwm[i][x]/total,x) for x in bases]
                l.sort()
                l.reverse()
                expected.append( expected_base )
                consensus.append( l[0][1] )
                if i>=2:
                    top_base = l[0][1]
                    if top_base!=expected_base or l[1][0]>0.35:
                        print 'whoah %2d act: %s %5.1f exp: %s %5.1f %s %s %d'\
                            %( i,l[0][1],100*l[0][0],expected_base,100*pwm[i][expected_base]/total,
                               gene,rep,
                               len(nucseqs) )
                        warn = True
            if warn:
                print 'whoah: expected: {} consensus: {} {} {}'\
                    .format( ''.join( expected ), ''.join( consensus ), gene,len(nucseqs))
    exit()









v_repsl = [ ( len(v_data[x]),x) for x in v_data if len(v_data[x]) >= min_count ]
v_repsl.sort()
v_repsl.reverse()

j_repsl = [ ( len(j_data[x]),x) for x in j_data if len(j_data[x]) >= min_count ]
j_repsl.sort()
j_repsl.reverse()


v_reps = [ x[1] for x in v_repsl ]
j_reps = [ x[1] for x in j_repsl ]
#v_reps = sorted( [ x[1] for x in v_repsl ] )
#j_reps = sorted( [ x[1] for x in j_repsl ] )

big_total = sum( vj_pairs.values() )

if True:
    plt.figure(1,figsize=(14,14))
    A = np.zeros( ( len(v_reps), len(j_reps) ) )
    for ii,v in enumerate( v_reps ):
        p_v = float( len( v_data[v] ) ) / big_total
        for jj,j in enumerate( j_reps ):
            p_j = float( len( j_data[j] ) ) / big_total
            expected = p_v * p_j * big_total
            actual = vj_pairs.get((v,j),min( expected, 0.25 ) )## 0.25 pseudo count
            enrich = math.log( actual/expected )
            A[ii][jj] = enrich

    A = A.transpose()

    vmax = math.log(5)
    vmin = -1*vmax
    v_rep_names = [ '{} ({:.2f}%)'.format( x, (100.0*len(v_data[x]))/big_total ) for x in v_reps ]
    j_rep_names = [ '{} ({:.2f}%)'.format( x, (100.0*len(j_data[x]))/big_total ) for x in j_reps ]
    plt.imshow( A, origin = 'lower', interpolation='nearest', cmap=plt.get_cmap('bwr'),vmin=vmin, vmax=vmax )
    plt.xticks( range(len(v_reps)), v_rep_names, rotation='vertical' )
    plt.yticks( range(len(j_reps)), j_rep_names )

    plt.suptitle('{}\n{}'.format( ' '.join(sys.argv), getcwd() ) )
    pngfile = 'tmp.read_read_nextgen_matches.gene_correlations.{}.png'.format('_'.join(logfile.split('/')))
    print 'making',pngfile
    plt.savefig(pngfile)


if True:
    ncols = len(dats_names)

    plt.figure(2,figsize=(23,14))
    all_all_dats = []
    for counter,(v_rep_count,v_rep) in enumerate( v_repsl ):
        all_dats = v_data[v_rep]
        all_all_dats.extend( all_dats )
        if skip_raw_data: continue
        for i in range(ncols):
            plt.subplot(2,ncols,i+1)
            ## histogram
            vals = [x[0][i] for x in all_dats] ## each entry in all_dats is a tuple ( dats0, ) or (dats0,dats1) if tied
            counts = {}
            for x in vals:
                counts[x] = counts.get(x,0)+1
            total = sum(counts.values())
            xmin=0
            if chain == 'B' and 'tot' in dats_names[i]:
                xmax = default_xmax+10
                #legend_loc = 'upper left'
            elif chain == 'B' and dats_names[i] == 'd_id':
                xmin = 1
                xmax = num_d_ids
                legend_loc = 'lower center'
                legend_loc = 'upper center'
            else:
                xmax = default_xmax
                legend_loc = 'upper right'
            xvals = range(xmin,xmax+1)
            yvals = [ float( counts.get(x,0))/total for x in xvals]
            plt.plot( xvals, yvals, get_style_and_color( counter ), label = '{} {}'.format(v_rep,v_rep_count))
            if chain == 'B' and dats_names[i] == 'd_id':
                plt.ylim((0,1.5))
            plt.legend( loc=legend_loc, fontsize = 7 )



    for counter,(j_rep_count,j_rep) in enumerate( j_repsl ):
        all_dats = j_data[j_rep]
        for i in range(ncols):
            plt.subplot(2,ncols,i+ncols+1)
            ## histogram
            vals = [x[0][i] for x in all_dats]
            counts = {}
            for x in vals:
                counts[x] = counts.get(x,0)+1
            total = sum(counts.values())
            xmin=0
            if chain == 'B' and 'tot' in dats_names[i]:
                xmax = default_xmax+10
            elif chain == 'B' and dats_names[i] == 'd_id':
                xmin = 1
                xmax = num_d_ids
                legend_loc = 'upper center'
            else:
                xmax = default_xmax
                legend_loc = 'upper right'

            xvals = range(xmin,xmax+1)
            yvals = [ float( counts.get(x,0))/total for x in xvals]
            plt.plot( xvals, yvals, get_style_and_color( counter ), label = '{} {}'.format(j_rep,j_rep_count))
            if chain == 'B' and dats_names[i] == 'd_id':
                plt.ylim((0,1.5))
            plt.legend( loc=legend_loc, fontsize = 7 )
    plt.subplots_adjust( hspace=0.25, wspace = 0.22, left=0.02, right = 0.98, bottom=0.02, top=0.95 )

    for i in range(ncols):
        plt.subplot(2,ncols,i+1)
        plt.title(dats_names[i])
        plt.subplot(2,ncols,i+ncols+1)
        plt.title(dats_names[i])

    plt.suptitle('{}\n{}'.format( ' '.join(sys.argv), getcwd() ) )
    pngfile = 'tmp.read_read_nextgen_matches.trims_and_inserts_by_gene.{}.png'.format('_'.join(logfile.split('/')))
    print 'making',pngfile
    plt.savefig(pngfile)



if not make_png:
    plt.show()



    # ## show poisson fits
    # if show_poisson:
    #     for i in range(ncols):
    #         plt.subplot(2,ncols,i+1)
    #         ## histogram
    #         vals = [x[i] for x in all_all_dats]
    #         mean = sum(vals)/len(vals)
    #         counts = {}
    #         for x in vals:
    #             counts[x] = counts.get(x,0)+1
    #         total = sum(counts.values())
    #         if chain == 'B' and dats_names[i] == 'tot_trim':
    #             xmax = default_xmax+5
    #             legend_loc = 'upper left'
    #         elif chain == 'B' and dats_names[i] == 'd_id':
    #             xmax = 3
    #         else:
    #             xmax = default_xmax
    #             legend_loc = 'upper right'
    #         xvals = range(xmax)
    #         yvals = [ float( counts.get(x,0))/total for x in xvals]
    #         plt.plot( xvals, yvals, '-ok', label = 'all_data')

    #         rv = poisson(mean)
    #         yvals = rv.pmf(xvals)
    #         plt.plot( xvals, yvals, '-or', label = 'poisson' )

    #         plt.legend( loc=legend_loc, fontsize = 7 )
