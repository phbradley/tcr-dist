from basic import *
from amino_acids import amino_acids
from all_genes import all_genes
import re
import logo_tools
import svg_basic
import tcr_sampler
import util
import parse_tsv
from basic import *
import random


nbr_distance_default = pipeline_params[ 'distance_threshold_25' ]

with Parser(locals()) as p:
    p.str('clones_file').required()
    p.str('organism')
    p.str('outfile_prefix')
    p.float('min_chi_squared').default(200)
    p.float('max_expected_fraction_for_clustering').default(0.05)
    p.float('min_top_chi_squared').default(75) ## allow lower chi-seq if we aren't going to see any otherwise
    #p.float('distance_scale_factor').default(0.01)
    p.float('nbr_distance').default( nbr_distance_default ) ## single-chain nbr distance threshold
    p.float('motifs_clustering_threshold').default(0.3)
    p.flag('verbose')
    p.flag('constant_seed')
    p.flag('paper_figs')
    p.flag('paper_supp')
    p.flag('junction_bars')
    p.int('max_ng_lines')
    p.int('target_num')
    p.multiword('epitopes').cast(lambda x:x.split())
    p.multiword('ABs').cast(lambda x:x.split())

if constant_seed:
    random.seed(1)

if outfile_prefix is None:
    outfile_prefix = clones_file[:-4]

# if paper_figs and max_expected_fraction==0.05:
#     max_expected_fraction = 1.0 #dont filter

if paper_figs or paper_supp:
    junction_bars = True

junction_bars_color = { 'V':  'black',
                        'N1': 'red',
                        'N':  'red',
                        'N2': 'red',
                        'D':  'blue',
                        'J':  'gray' }

junction_bars_color = { 'V':  'silver',
                        'N1': 'red',
                        'N':  'red',
                        'N2': 'red',
                        'D':  'darkslategray', #'black',
                        'J':  'dimgray' }

junction_bars_color = { 'V':  '#D3D3D3',
                        'N1': 'red',
                        'N':  'red',
                        'N2': 'red',
                        'D':  '#696969', #'black',
                        'J':  '#A9A9A9' }

junction_bars_color = { 'V':  'silver',
                        'N1': 'red',
                        'N':  'red',
                        'N2': 'red',
                        'D':  'black',
                        'J':  'dimgray' }

junction_bars_order = { 'B': ['V','N1','D','N2','J'],
                        'A': ['V','N','J'] }


big = True

## this is the naming convention in run_cdr3_motifs.faster.nextgen.py
##
#motifs_files = glob('{}_cdr3_motifs_*.log'.format(clones_file[:-4] ) )
#assert motifs_files

#min_chi_squared = 200
max_overlap = 50


ypad = 30
gene_stack_width = 100
pwm_stack_width = 300
pwm_stack_height = 80
xmargin = 0 if paper_figs else 15
sep = 5 if paper_figs else 15
fontsize = pwm_stack_height/3

junction_bars_ypad = 2

max_column_relent_for_scaling = 3.0
min_prob_for_relent_for_scaling = 1e-3

num_nextgen_samples = 100
num_nextgen_samples = 25

### this is stolen from cdr3_motifs.faster.nextgen.py
groups = dict( zip( amino_acids, amino_acids ) )

groups['k'] = '[KR]'
groups['d'] = '[DE]'
groups['n'] = '[NQ]'
groups['s'] = '[ST]'
groups['f'] = '[FYWH]'
groups['a'] = '[AGSP]'
groups['v'] = '[VILM]'

begin = '^'
end = '$'
X = '[A-Z]'
dot = '.'

groups[begin] = begin
groups[end  ] = end
groups[dot] = X

def get_amino_acid_consensus_character( counts ):
    topcount,topaa = max( ( (y,x) for x,y in counts.iteritems() ) )
    frac = float( topcount)/ sum(counts.values())
    if frac >= 0.75:
        return topaa
    elif frac >= 0.4:
        return topaa.lower()
    else:
        return  '.'


def get_amino_acid_consensus( seqs ):
    numseqs = len(seqs)
    L = len(seqs[0])
    consensus = ''
    for i in range(L):
        counts = {}
        for seq in seqs:
            assert len(seq) == L
            counts[seq[i]] = counts.get(seq[i],0)+1
        consensus += get_amino_acid_consensus_character( counts )
    assert len(consensus) == L
    return consensus


#logfile = 'all_good_probs_exact_genes.nonuniqd.txt' ; l_offset = 13
#logfile2 = 'all_good_probs_exact_genes.nonuniqd.problines'

#motifs_files = glob('tmp.cdr3_motifs.faster.nextgen.*.groups.fixlenFalse.25samples.max100.redo.log')
#motifs_files = glob('tmp.cdr3_motifs.faster.nextgen.*.groups.fixlenFalse.25samples.max100.mc10.redo.ngsample.log')


#motifs_files.sort()


## new approach
all_tcr_infos = parse_tsv.parse_tsv_file( clones_file, ['epitope'], [], True )

if epitopes is None:
    epitopes = all_tcr_infos.keys()[:]
    epitopes.sort()


all_tcrs = {}
all_rep2label_rep = {}
all_rep2label_rep_color = {}

for epitope in epitopes:
    infos = all_tcr_infos[epitope]
    util.assign_label_reps_and_colors_based_on_most_common_genes_in_repertoire( infos, organism )

    all_tcrs[epitope] = []
    all_rep2label_rep[epitope] = {}
    all_rep2label_rep_color[epitope] = {}

    for l in infos:

        epitope = l['epitope']
        va_gene = l['va_gene']
        ja_gene = l['ja_gene']
        vb_gene = l['vb_gene']
        jb_gene = l['jb_gene']
        va = l['va_rep']
        ja = l['ja_rep']
        vb = l['vb_rep']
        jb = l['jb_rep']
        cdr3a = l['cdr3a']
        cdr3b = l['cdr3b']
        cdr3a_nucseq = l['cdr3a_nucseq']
        cdr3b_nucseq = l['cdr3b_nucseq']

        ## note-- we are using mm1 reps here, same as in find_cdr3_motifs.py
        ##
        va_rep = all_genes[organism][va].mm1_rep
        ja_rep = all_genes[organism][ja].rep
        vb_rep = all_genes[organism][vb].mm1_rep
        jb_rep = all_genes[organism][jb].rep

        a_junction_results = tcr_sampler.analyze_junction( organism, va_gene, ja_gene, cdr3a, cdr3a_nucseq,
                                                           return_cdr3_nucseq_src=True )
        b_junction_results = tcr_sampler.analyze_junction( organism, vb_gene, jb_gene, cdr3b, cdr3b_nucseq,
                                                           return_cdr3_nucseq_src=True )

        cdr3a_new_nucseq, cdr3a_protseq_masked, cdr3a_protseq_new_nucleotide_countstring,\
            a_trims,a_inserts,cdr3a_nucseq_src = a_junction_results
        cdr3b_new_nucseq, cdr3b_protseq_masked, cdr3b_protseq_new_nucleotide_countstring,\
            b_trims,b_inserts,cdr3b_nucseq_src = b_junction_results

        assert len(cdr3a_nucseq_src) == 3*len(cdr3a)
        assert len(cdr3b_nucseq_src) == 3*len(cdr3b)

        assert type(cdr3b_nucseq_src) == type([])## note

        if junction_bars: ## try to distinguish between N before D and N after D
            for i in range(len(cdr3b_nucseq_src)):
                if cdr3b_nucseq_src[i] == 'N':
                    if cdr3b_nucseq_src[:i].count('D')==0:
                        cdr3b_nucseq_src[i] = 'N1'
                    else:
                        cdr3b_nucseq_src[i] = 'N2'

        ## let's reconstruct where everything comes from in the cdr3a and cdr3b sequences
        all_tcrs[ epitope ].append( ( va, ja, vb, jb, cdr3a, cdr3b, va_rep, ja_rep, vb_rep, jb_rep,
                                      cdr3a_nucseq_src, cdr3b_nucseq_src, l ) )

        all_rep2label_rep[ epitope ][ va_rep ] = l['va_label_rep']
        all_rep2label_rep[ epitope ][ ja_rep ] = l['ja_label_rep']
        all_rep2label_rep[ epitope ][ vb_rep ] = l['vb_label_rep']
        all_rep2label_rep[ epitope ][ jb_rep ] = l['jb_label_rep']
        all_rep2label_rep_color[ epitope ][ va_rep ] = l['va_label_rep_color']
        all_rep2label_rep_color[ epitope ][ ja_rep ] = l['ja_label_rep_color']
        all_rep2label_rep_color[ epitope ][ vb_rep ] = l['vb_label_rep_color']
        all_rep2label_rep_color[ epitope ][ jb_rep ] = l['jb_label_rep_color']


ng_tcrs = { 'A':{}, 'B':{} }

## index these by the v_rep and the j_rep

for ab in 'AB':
    ng_logfile = '{}/new_nextgen_chains_{}_{}.tsv'.format( path_to_current_db_files(), organism, ab )
    if not exists(ng_logfile):
        Log('WARNING:: read_motifs.py: missing nextgen TCR chains file: {}'.format(ng_logfile))
        continue

    counter=0
    num_chains=0
    ab_chains = {}

    for line in open(ng_logfile,'r'):
        counter+=1
        if max_ng_lines and counter>max_ng_lines:break
        l = line[:-1].split('\t')
        if counter==1:
            assert l==['v_reps','j_reps','cdr3','cdr3_nucseq']
            continue
        if not counter%500000:Log(`counter`+' '+`num_chains`+' '+ng_logfile)
        v_reps = set( ( util.get_mm1_rep( x, organism ) for x in l[0].split(',') ) ) ## mm1 reps
        j_reps = l[1].split(',')
        cdr3,cdr3_nucseq = l[2:4]

        ## now add to the different places
        for v_rep in v_reps:
            for j_rep in j_reps:
                if v_rep not in ab_chains: ab_chains[v_rep] = {}
                if j_rep not in ab_chains[v_rep]: ab_chains[v_rep][j_rep] = []
                ab_chains[v_rep][j_rep].append( (cdr3, cdr3_nucseq ))

        num_chains += 1
    Log('read {} {}-chains from {}'.format(num_chains,ab,ng_logfile))
    ng_tcrs[ab] = ab_chains


total_y = ypad

cmds = []


svg_width = 0

motifid_counter = 0

motif_trees_info = {}

for epitope in epitopes:

    if paper_supp: ## separate fig for each epitope,chain
        cmds = []
        total_y = 0
        ab_cmds = { 'A':[], 'B':[] }
        ab_total_y = {}

    mfile = '{}_cdr3_motifs_{}.log'.format( clones_file[:-4], epitope )
    if not exists(mfile):
        Log('missing motif file: {} {}'.format(epitope,mfile))
        continue

    #epitope = mfile[:-4].split('_')[-1]
    tcrs = all_tcrs[epitope]
    num_tcrs = len(tcrs)

    ## hack so we can use a consistent rep coloring scheme relative to kpca, trees, gene plots etc
    rep2label_rep = all_rep2label_rep[ epitope ]
    rep2label_rep_color = all_rep2label_rep_color[ epitope ]

    ## load tcr distances
    all_distances = {}
    all_neighbors = {}
    for ab in 'AB':
        distfile = '{}_{}_{}.dist'.format(clones_file[:-4],ab,epitope)
        assert exists(distfile)
        N=0
        all_nbrs = []
        all_dists = []
        for line in open( distfile,'r'):
            l = line.split()
            clone_id = l[0]
            index = len(all_dists)
            assert tcrs[ index ][-1]['clone_id'] == clone_id
            dists = [ float(x) for x in l[1:] ]
            if not N:
                N = len(dists)
            else:
                assert N == len(dists)

            nbrs = []
            for ii,d in enumerate(dists):
                if d <= nbr_distance:
                    nbrs.append( ii )
            all_dists.append( dists )
            all_nbrs.append( nbrs )

        assert len(all_nbrs) == len(tcrs)

        all_distances[ab] = all_dists
        all_neighbors[ab] = all_nbrs


    ## for making a dendrogram of all the motifs
    motif_dists = {'A':[], 'B':[] }
    motif_matches = {'A':[], 'B':[] }
    motif_infos = {'A':[], 'B':[] }
    motif_positions = {'A':[], 'B':[] }


    def create_wtd_pwm_from_sequences( seqs, alphabet, target_reps, reps ):
        assert len(seqs) == len(reps)
        num_target_reps = len(target_reps)
        # repwts = {}
        # for ii in range(2):
        #     ii_target_reps = [x[ii] for x in target_reps]
        #     ii_reps = [x[ii] for x in reps]
        #     #ii_repwts = {}
        #     #total_wt = 0.
        #     for rep in ii_reps:
        #         count = ii_reps.count(rep)
        #         target_count = ii_target_reps.count(rep)
        #         assert target_count
        #         wt = target_count / float(count)
        #         #total_wt += wt
        #         #ii_repwts[rep] = wt
        #         repwts[rep] = wt
        #     print 'total_wt:',sum( (repwts[x] for x in ii_reps) ), ii,len(seqs),len(target_reps)
        #     # for rep in ii_repwts:
        #     #     repwts[rep] = float( ii_repwts[rep] ) / total_wt

        # for rep in sorted( repwts.keys()):
        #     print 'repwt:',rep,repwts[rep]

        ## now iteratively adjust the wts on rep-pairs

        for bigrepeat in range(10):
            reppair_wts = {}
            if bigrepeat==0:
                for rp in reps:
                    reppair_wts[rp] = 1.0 ## starting guess
            else:
                for rp in reps:
                    reppair_wts[rp] = 0.75 + 0.5 * random.random() ## starting guess

            prev_dev = 1e6
            for repeat in range(100):

                ## what's the deviation
                dev = 0.0
                for ii in range(2):
                    ii_target_reps = [x[ii] for x in target_reps]
                    ii_reps = [x[ii] for x in reps]

                    scale_factor = float( len(reps ) )/ len(target_reps)

                    counts = {}
                    for rp in reps:
                        counts[rp[ii]] = counts.get(rp[ii],0) + reppair_wts[rp]

                    for rep,count in counts.iteritems():
                        desired_count = scale_factor * ii_target_reps.count(rep)
                        dev += abs( desired_count - count )
                        fac = float(desired_count)/count
                        adjust = fac**(0.25)

                        #print 'desired_count:',desired_count,'count:',count,'fac:',fac,'adjust:',adjust,rep

                        for rp in reppair_wts:
                            if rp[ii] == rep:
                                reppair_wts[rp] *= adjust
                #print 'repeat:',repeat,'dev:',dev
                if abs(prev_dev-dev)<1e-3 and dev<1e-1:
                #if abs(dev)<1e-1:
                    break
                prev_dev = dev

            print 'final_dev:',bigrepeat,dev
            if dev<1e-1:
                break


        L = len(seqs[0])
        pwm = {}
        for i in range(L):
            pwm[i] = dict(zip(alphabet,[0.0]*len(alphabet)))

        for seq,rp in zip( seqs, reps ):
            assert len(seq) == L
            seqwt = reppair_wts[ rp ]
            #print seq, rp, seqwt
            for i,a in enumerate(seq):
                pwm[i][a] += seqwt

        for i in range(L):
            tot = sum( pwm[i].values() )
            for a in alphabet:
                pwm[i][a] /= tot
        return pwm


    def analyze_matches_using_ngseqs( matches, matched_tcrs, ab ):
        global tcrs

        ng_lenseqs = []
        ng_fwdseqs = []
        ng_revseqs = []

        ng_fwdseq_reps = []
        ng_lenseq_reps = []
        matched_reps = []

        seen = set() ## no repeats of ngseqs
        seen_samelen = set() ## no repeats of ngseqs

        for (mseq,nseq,positions,rpositions),ii in zip( matches, matched_tcrs ):
            tcr = tcrs[ii]
            if ab == 'A':
                my_cdr3,vrep,jrep = tcr[4:5]+tcr[6: 8]
            else:
                my_cdr3,vrep,jrep = tcr[5:6]+tcr[8:10]
            matched_reps.append( ( vrep, jrep ) )

            mylen = len(my_cdr3)
            if vrep in ng_tcrs[ab] and jrep in ng_tcrs[ab][vrep]:
                ngl = [ x for x in ng_tcrs[ab][vrep][jrep] if x not in seen ]
                if not ngl:
                    print 'empty ngl!'

                for ngseq in random.sample( ngl, min(num_nextgen_samples,len(ngl)) ):
                    seen.add(ngseq)
                    (cdr3,cdr3_nucseq) = ngseq
                    L = len(cdr3)
                    fseq = ''
                    rseq = ''
                    for pos in positions:
                        if pos>=L:
                            fseq += '-'
                        else:
                            fseq += cdr3[pos]
                    for pos in rpositions:
                        if pos>=L:
                            rseq += '-'
                        else:
                            rseq += cdr3[L-1-pos]
                    ng_fwdseqs.append(fseq)
                    ng_revseqs.append(rseq)
                    ng_fwdseq_reps.append( ( vrep, jrep ) )

                ## cdr3s with the same length
                ngl_samelen = [ x for x in ng_tcrs[ab][vrep][jrep] if len(x[0]) == mylen and x not in seen_samelen ]
                if not ngl_samelen:
                    print 'empty ngl_samelen!'
                for ngseq in random.sample( ngl_samelen, min(num_nextgen_samples,len(ngl_samelen))):
                    seen_samelen.add( ngseq )
                    cdr3 = ngseq[0]
                    ng_lenseqs.append( ''.join( [ cdr3[x] for x in positions ] ) )
                    ng_lenseq_reps.append( ( vrep, jrep ) )

        pwm = logo_tools.create_protein_pwm_from_sequences( [x[0] for x in matches ])

        npwm_alphabet = junction_bars_order[ab] if junction_bars else ['V','N','D','J']
        npwm = logo_tools.create_pwm_from_sequences( [x[1] for x in matches ], npwm_alphabet )

        #nbr_pwm = logo_tools.create_protein_pwm_from_sequences( [x[0] for x in nbr_matches ])
        #nbr_npwm = logo_tools.create_pwm_from_sequences( [x[1] for x in nbr_matches ], ['V','D','J','N'] )

        if ng_lenseqs:
            ng_lenpwm = create_wtd_pwm_from_sequences( ng_lenseqs, amino_acids+['-'], matched_reps, ng_lenseq_reps )
        ng_fwdpwm = create_wtd_pwm_from_sequences( ng_fwdseqs, amino_acids+['-'], matched_reps, ng_fwdseq_reps )
        ng_revpwm = create_wtd_pwm_from_sequences( ng_revseqs, amino_acids+['-'], matched_reps, ng_fwdseq_reps )

        N = len(pwm)
        fwdpwm = {}
        revpwm = {}
        for i in range(N):
            fwdpwm[i] = {}
            revpwm[i] = {}
            incrememnt = 1.0/len(matches)
            for pos in [x[2][i] for x in matches]:
                fwdpwm[i][`pos`] = fwdpwm[i].get(`pos`,0)+incrememnt
            for pos in [x[3][i] for x in matches]:
                revpwm[i][`pos`] = revpwm[i].get(`pos`,0)+incrememnt

        ## look at relative entropies between nbrpwm and the fwd and rev pwms
        ## not nbr anymore since this is a subroutine
        ##
        scale_by_relent = {}
        for i in range(N):
            relents=[]
            for control_pwm in [ ng_fwdpwm[i], ng_revpwm[i] ]:
                relent = 0.0
                for a,pa in pwm[i].iteritems():
                    if pa>= min_prob_for_relent_for_scaling:
                        qa = max(min_prob_for_relent_for_scaling, control_pwm.get(a,min_prob_for_relent_for_scaling))
                        relent += pa * math.log(pa/qa, 2.0 )
                relents.append( relent )
            scale_by_relent[i] = max(0.,min(1., min(relents)/max_column_relent_for_scaling) )
            print 'RE {:2d} {:5.2f} {:5.2f} {:5.2f} {} {} {}'.format( i, min(relents), relents[0], relents[1], ab, epitope, ''.join(showmotif) )

        return pwm, npwm, ng_lenpwm, ng_fwdpwm, ng_revpwm, fwdpwm, revpwm, scale_by_relent, \
            ng_fwdseq_reps, ng_lenseq_reps, len( ng_lenseqs ), len( ng_fwdseqs )


    for line in open( mfile,'r'):
        l = line.split()

        if l[0] != 'MOTIF': continue

        count, expect_random, expect_nextgen, chi_squared, nfixed, showmotif, num, othernum, overlap, ep, ab, nseqs, v_rep_counts, j_rep_counts = l[1:]
        count,nfixed,num,othernum,overlap,nseqs = map(int,[count,nfixed,num,othernum,overlap,nseqs] )
        expect_random, expect_nextgen, chi_squared = map(float, [ expect_random, expect_nextgen, chi_squared ] )
        showmotif = list(showmotif)
        assert ep==epitope

        if target_num != None and target_num != num: continue

        if ABs != None and ab not in ABs: continue

        expected_fraction = max(expect_random,expect_nextgen)/num_tcrs

        # if num>1 and chi_squared<min_chi_squared: continue
        # if num==1 and chi_squared<min_top_chi_squared: continue
        # if overlap>max_overlap: continue

        print line[:-1]

        motif = [groups[x] for x in showmotif]

        # v_reps = frozenset( [x.split(':')[0] for x in l[-2].split(',') ] )
        # j_reps = frozenset( [x.split(':')[0] for x in l[-1].split(',') ] )

        ## look for matches in the ab cdr3s
        prog = re.compile(''.join(motif))

        total=0
        matches = []
        nbr_matches = []

        matched_tcrs = []
        matched_tcrs_plus_nbrs = []

        for ii, tcr in enumerate( tcrs ):
            if ab == 'A':
                cdr3,cdr3_nucseq_src = tcr[4], tcr[10]
            else:
                cdr3,cdr3_nucseq_src = tcr[5], tcr[11]
            m = prog.search(cdr3)
            if m:
                mseq = cdr3[ m.start():m.end() ]
                nseq = cdr3_nucseq_src[ 3*m.start():3*m.end() ]
                positions = range(m.start(),m.end())
                rpositions = [len(cdr3)-1-x for x in positions]
                matches.append( (mseq,nseq,positions,rpositions) )
                matched_tcrs.append( ii )
                for nbr in all_neighbors[ab][ii]:
                    if nbr not in matched_tcrs_plus_nbrs:
                        nbr_tcr = tcrs[nbr]
                        if ab == 'A':
                            nbr_cdr3, nbr_cdr3_nucseq_src = nbr_tcr[4], nbr_tcr[10]
                        else:
                            nbr_cdr3, nbr_cdr3_nucseq_src = nbr_tcr[5], nbr_tcr[11]
                        if len(nbr_cdr3) == len(cdr3):
                            matched_tcrs_plus_nbrs.append( nbr )

                            nbr_mseq = nbr_cdr3[ m.start():m.end() ]
                            nbr_nseq = nbr_cdr3_nucseq_src[ 3*m.start():3*m.end() ]
                            nbr_matches.append( (nbr_mseq,nbr_nseq,positions,rpositions) )


            total += 1
        assert count == len(matches)
        assert total == nseqs


        nbr_consensus = get_amino_acid_consensus( [ x[0] for x in nbr_matches ] )

        ## write out some info to the tsv file
        motifid_counter += 1
        motifid = 'motif{}'.format(motifid_counter)
        matched_tcrs_ids = [ tcrs[x][-1]['clone_id'] for x in matched_tcrs ]
        matched_tcrs_plus_nbrs_ids = [ tcrs[x][-1]['clone_id'] for x in matched_tcrs_plus_nbrs ]

        tsv_outl = { 'id': motifid,
                     'epitope': epitope,
                     'chain': ab,
                     'epitope_num':num,
                     'chi_squared':chi_squared,
                     'overlap':overlap,
                     'num_matches': len(matched_tcrs_ids ),
                     'num_matches_with_nbrs': len(matched_tcrs_plus_nbrs_ids ),
                     'expected_fraction': expected_fraction,
                     'expected_num_matches': max(expect_random,expect_nextgen),
                     'showmotif':''.join(showmotif),
                     'matches':','.join( matched_tcrs_ids ),
                     'matches_with_nbrs':','.join( matched_tcrs_plus_nbrs_ids ),
                     'matches_with_nbrs_consensus':nbr_consensus }


        #out.write( make_tsv_line( tsv_outl, tsv_outfields )+'\n' )


        ## compute distances to previous motifs
        matched_tcrs_plus_nbrs_set = frozenset( matched_tcrs_plus_nbrs )
        num_matched_tcrs_plus_nbrs = len( matched_tcrs_plus_nbrs )

        my_tree_label = '{:7.1f} {:4d} {:20s}'.format(chi_squared, num_matched_tcrs_plus_nbrs, ''.join(showmotif) )

        dists = []
        for ii,prev_matches in enumerate( motif_matches[ab] ):
            intersection = 0
            for a in matched_tcrs_plus_nbrs_set:
                if a in prev_matches:
                    intersection += 1
            dist = 1.0 - float(intersection)/max( len(prev_matches), num_matched_tcrs_plus_nbrs )
            dists.append( dist )
            if dist<1e-3:
                print 'zero-dist:',epitope,ab,ii,len(motif_matches[ab]),dist, my_tree_label, motif_infos[ab][ii][0]



        ## hacking
        for ind,(mseq,nseq,positions,rpositions) in zip( matched_tcrs_plus_nbrs, nbr_matches ):
            tcr = tcrs[ind]
            cdr3 = tcr[4] if ab == 'A' else tcr[5]
            assert mseq == cdr3[ min(positions):max(positions)+1]


        if expected_fraction <= max_expected_fraction_for_clustering:
            ## save for dendrogram
            motif_matches[ab].append( matched_tcrs_plus_nbrs_set )
            motif_positions[ab].append( dict( (x, (y[0],y[2]) ) for x,y in zip( matched_tcrs_plus_nbrs, nbr_matches ) ))
            motif_dists[ab].append( dists )
            motif_infos[ab].append( ( my_tree_label, tsv_outl ) ) ## we still need to add clustering info to tsv_outl

        ## filter here for inclusion in the graphic
        if num>1 and chi_squared<min_chi_squared: continue
        if num==1 and chi_squared<min_top_chi_squared: continue
        if overlap>max_overlap: continue


        pwm, npwm, ng_lenpwm, ng_fwdpwm, ng_revpwm, fwdpwm, revpwm, scale_by_relent, ng_fwdseq_reps, ng_lenseq_reps, \
            num_ng_lenseqs, num_ng_fwdseqs \
            = analyze_matches_using_ngseqs( matches, matched_tcrs, ab )


        nbr_pwm, nbr_npwm, nbr_ng_lenpwm, nbr_ng_fwdpwm, nbr_ng_revpwm, nbr_fwdpwm, nbr_revpwm, nbr_scale_by_relent, \
            nbr_ng_fwdseq_reps, nbr_ng_lenseq_reps, nbr_num_ng_lenseqs, nbr_num_ng_fwdseqs \
            = analyze_matches_using_ngseqs( nbr_matches, matched_tcrs_plus_nbrs, ab )



        ## make a v-gene logo
        def get_counts_list_condensing_alleles( counts_string ):
            global rep2label_rep
            global rep2label_rep_color
            counts ={}
            for tag,count in [x.split(':') for x in counts_string.split(',') ]:
                rc = ( rep2label_rep[ tag ][4:], rep2label_rep_color[ tag ] )
                counts[rc] = counts.get(rc,0)+float(count)
            return [ (y,x[0],x[1]) for x,y in counts.iteritems() ]

        def get_counts_lists_from_tcr_indices( indices ):
            vcounts = {}
            jcounts = {}
            for ii in indices:
                tcr = tcrs[ii]
                if ab == 'A':
                    vrep,jrep = tcr[6: 8]
                else:
                    vrep,jrep = tcr[8:10]
                vcounts[vrep] = vcounts.get(vrep,0)+1
                jcounts[jrep] = jcounts.get(jrep,0)+1
            vstring = ','.join( ['{}:{}'.format(x,y) for x,y in vcounts.iteritems()] )
            jstring = ','.join( ['{}:{}'.format(x,y) for x,y in jcounts.iteritems()] )
            return get_counts_list_condensing_alleles(vstring), get_counts_list_condensing_alleles(jstring)

        def get_counts_lists_from_rep_lists( reps ):
            vcounts = {}
            jcounts = {}
            for vrep,jrep in reps:
                vcounts[vrep] = vcounts.get(vrep,0)+1
                jcounts[jrep] = jcounts.get(jrep,0)+1
            vstring = ','.join( ['{}:{}'.format(x,y) for x,y in vcounts.iteritems()] )
            jstring = ','.join( ['{}:{}'.format(x,y) for x,y in jcounts.iteritems()] )
            return get_counts_list_condensing_alleles(vstring), get_counts_list_condensing_alleles(jstring)


        #vl = get_counts_list_condensing_alleles( v_rep_counts ) # [ (float(x.split(':')[1]), x.split(':')[0][2:] ) for x in v_rep_counts.split(',') ]
        #jl = get_counts_list_condensing_alleles( j_rep_counts ) # [ (float(x.split(':')[1]), x.split(':')[0][2:] ) for x in j_rep_counts.split(',') ]

        vl, jl = get_counts_lists_from_tcr_indices( matched_tcrs )
        nbr_vl, nbr_jl = get_counts_lists_from_tcr_indices( matched_tcrs_plus_nbrs )
        if ng_lenseq_reps: ng_lenseq_vl, ng_lenseq_jl = get_counts_lists_from_rep_lists( ng_lenseq_reps )
        ng_fwdseq_vl, ng_fwdseq_jl = get_counts_lists_from_rep_lists( ng_fwdseq_reps )



        ########################################################################### VERBOSE VIEW


        if paper_figs:
            cmds = []
            total_y = 0


        if paper_supp:
            ## we want to create two images for each epitope, one with all the alpha motifs and one with all the
            ## beta motifs
            if num==1 and ab=='B': # start over with image
                ab_cmds['A'] = cmds[:]
                ab_total_y['A'] = total_y
                cmds = []
                total_y = 0


        if verbose:




            x0=xmargin; x1 = x0+gene_stack_width

            y0=total_y
            y1=y0+    pwm_stack_height ## pwm
            y2=y1+0.5*pwm_stack_height ## VNDNJ pwm
            y3=y2+    pwm_stack_height ## nbr-pwm
            y4=y3+0.5*pwm_stack_height ## nbr-NVDNJ pwm
            y5=y4+0.5*pwm_stack_height ## lenseq pwm
            y6=y5+0.5*pwm_stack_height ## fwd pwm
            y7=y6+0.5*pwm_stack_height ## ng-fwd pwm
            y8=y7+0.5*pwm_stack_height ## rev pwm
            y9=y8+0.5*pwm_stack_height ## ng-rev pwm
            y10=y9+   pwm_stack_height ## nbr-pwm scaled by relent

            y_last = y10

            p0,p1 = [ ( x0, y0 ), ( x1, y1 ) ]
            cmds.append( svg_basic.make_stack( p0, p1, vl ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            p0,p1 = [ ( x0, y2 ), ( x1, y3 ) ]
            cmds.append( svg_basic.make_stack( p0, p1, nbr_vl ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            if ng_lenseq_reps:
                p0,p1 = [ ( x0, y4 ), ( x1, y5 )]
                cmds.append( svg_basic.make_stack( p0, p1, ng_lenseq_vl ) )
                cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            p0,p1 = [ ( x0, y6 ), ( x1, y7 ) ]
            cmds.append( svg_basic.make_stack( p0, p1, ng_fwdseq_vl ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            x0 = x1+sep; x1 = x0+pwm_stack_width
            ## pwm
            p0,p1 = [ ( x0, y0 ), (x1,y1 ) ]
            cmds.append( svg_basic.protein_logo( p0, p1, pwm ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            ## VNDNJ-pwm
            p0,p1 = [ ( x0, y1 ), ( x1, y2 ) ]
            cmds.append( svg_basic.protein_logo( p0, p1, npwm ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            ## nbr pwm
            p0,p1 = [ ( x0, y2 ), ( x1, y3 ) ]
            cmds.append( svg_basic.protein_logo( p0, p1, nbr_pwm ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            ## nbr N-pwm
            p0,p1 = [ ( x0, y3 ), ( x1, y4 ) ]
            cmds.append( svg_basic.protein_logo( p0, p1, nbr_npwm ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )


            ## ng fixlen pwm
            if ng_lenseq_reps:
                p0,p1 = [ ( x0, y4 ), ( x1, y5 ) ]
                cmds.append( svg_basic.protein_logo( p0, p1, ng_lenpwm ) )
                cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            ## fwd position pwm
            p0,p1 = [ ( x0, y5 ), ( x1, y6 ) ]
            cmds.append( svg_basic.generic_logo( p0, p1, fwdpwm ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            ## ng fwd pwm
            p0,p1 = [ ( x0, y6 ), ( x1, y7 ) ]
            cmds.append( svg_basic.protein_logo( p0, p1, ng_fwdpwm ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            ## rev position pwm
            p0,p1 = [ ( x0, y7 ), ( x1, y8 ) ]
            cmds.append( svg_basic.generic_logo( p0, p1, revpwm ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            ## ng rev pwm
            p0,p1 = [ ( x0, y8 ), ( x1, y9 ) ]
            cmds.append( svg_basic.protein_logo( p0, p1, ng_revpwm ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            ## nbr pwm SCALED by relent
            p0,p1 = [ ( x0, y9 ), ( x1, y10 ) ]
            cmds.append( svg_basic.protein_logo( p0, p1, nbr_pwm, nbr_scale_by_relent ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )




            x0=x1+sep ; x1 = x0+gene_stack_width
            p0,p1 = [ ( x0, y0 ), ( x1, y1 ) ]
            cmds.append( svg_basic.make_stack( p0, p1, jl ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            p0,p1 = [ ( x0, y2 ), ( x1, y3 ) ]
            cmds.append( svg_basic.make_stack( p0, p1, nbr_jl ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            if ng_lenseq_reps:
                p0,p1 = [ ( x0, y4 ), ( x1, y5 ) ]
                cmds.append( svg_basic.make_stack( p0, p1, ng_lenseq_jl ) )
                cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            p0,p1 = [ ( x0, y6 ), ( x1, y7 ) ]
            cmds.append( svg_basic.make_stack( p0, p1, ng_fwdseq_jl ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

        else: ############################################################################ COMPACT VIEW


            x0=xmargin; x1 = x0+gene_stack_width

            y0=total_y
            y1=y0 ## pwm
            y2=y1 ## VNDNJ pwm
            y3=y2+    pwm_stack_height ## nbr-pwm
            y4=y3+0.5*pwm_stack_height ## nbr-NVDNJ pwm
            y5=y4 ## lenseq pwm
            y6=y5 ## fwd pwm
            y7=y6 ## ng-fwd pwm
            y8=y7 ## rev pwm
            y9=y8 ## ng-rev pwm
            y10=y9+   pwm_stack_height ## nbr-pwm scaled by relent
            y_last = y10


            # p0,p1 = [ ( x0, y0 ), ( x1, y1 ) ]
            # cmds.append( svg_basic.make_stack( p0, p1, vl ) )
            # cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            p0,p1 = [ ( x0, y2 ), ( x1, y3 ) ]
            cmds.append( svg_basic.make_stack( p0, p1, nbr_vl ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            # if ng_lenseqs:
            #     p0,p1 = [ ( x0, y4 ), ( x1, y5 )]
            #     cmds.append( svg_basic.make_stack( p0, p1, ng_lenseq_vl ) )
            #     cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            # p0,p1 = [ ( x0, y6 ), ( x1, y7 ) ]
            # cmds.append( svg_basic.make_stack( p0, p1, ng_fwdseq_vl ) )
            # cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            x0 = x1+sep; x1 = x0+pwm_stack_width
            ## pwm
            # p0,p1 = [ ( x0, y0 ), (x1,y1 ) ]
            # cmds.append( svg_basic.protein_logo( p0, p1, pwm ) )
            # cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            # ## VNDNJ-pwm
            # p0,p1 = [ ( x0, y1 ), ( x1, y2 ) ]
            # cmds.append( svg_basic.protein_logo( p0, p1, npwm ) )
            # cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            ## nbr pwm
            p0,p1 = [ ( x0, y2 ), ( x1, y3 ) ]
            cmds.append( svg_basic.protein_logo( p0, p1, nbr_pwm ) )
            #cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            ## nbr N-pwm
            if junction_bars:
                junction_pwm = nbr_npwm
                junction_bars_height = y4-y3 - 2*junction_bars_ypad
                ncols = len( nbr_npwm.keys() )
                junction_bar_width = ( x1-x0 )/float(ncols)
                for j in range(ncols):
                    lcol = [ ( junction_pwm[j][x],x) for x in junction_bars_order[ab] ]
                    y1shift = y3+junction_bars_ypad
                    ## largest at the top
                    for frac,a in lcol:
                        y1shift_next = y1shift + frac * junction_bars_height
                        color = junction_bars_color[ a ]
                        p0 = [ x0+ j   *junction_bar_width, y1shift]
                        p1 = [ x0+(j+1)*junction_bar_width, y1shift_next ]
                        cmds.append( svg_basic.rectangle( p0, p1, fill=color, stroke=color ) )
                        y1shift = y1shift_next
                    assert abs( y1shift+junction_bars_ypad-y4)<1e-3
            else:
                p0,p1 = [ ( x0, y3 ), ( x1, y4 ) ]
                cmds.append( svg_basic.protein_logo( p0, p1, nbr_npwm ) )
                cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )


            ## ng fixlen pwm
            # if ng_lenseqs:
            #     p0,p1 = [ ( x0, y4 ), ( x1, y5 ) ]
            #     cmds.append( svg_basic.protein_logo( p0, p1, ng_lenpwm ) )
            #     cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            # ## fwd position pwm
            # p0,p1 = [ ( x0, y5 ), ( x1, y6 ) ]
            # cmds.append( svg_basic.generic_logo( p0, p1, fwdpwm ) )
            # cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            # ## ng fwd pwm
            # p0,p1 = [ ( x0, y6 ), ( x1, y7 ) ]
            # cmds.append( svg_basic.protein_logo( p0, p1, ng_fwdpwm ) )
            # cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            # ## rev position pwm
            # p0,p1 = [ ( x0, y7 ), ( x1, y8 ) ]
            # cmds.append( svg_basic.generic_logo( p0, p1, revpwm ) )
            # cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            # ## ng rev pwm
            # p0,p1 = [ ( x0, y8 ), ( x1, y9 ) ]
            # cmds.append( svg_basic.protein_logo( p0, p1, ng_revpwm ) )
            # cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            ## nbr pwm SCALED by relent
            p0,p1 = [ ( x0, y9 ), ( x1, y10 ) ]
            cmds.append( svg_basic.protein_logo( p0, p1, nbr_pwm, scale_by_relent ) )
            ## box around all 3 of them
            p0,p1 = [ ( x0, y2 ), ( x1, y10 ) ]
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )




            x0=x1+sep ; x1 = x0+gene_stack_width
            # p0,p1 = [ ( x0, y0 ), ( x1, y1 ) ]
            # cmds.append( svg_basic.make_stack( p0, p1, jl ) )
            # cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            p0,p1 = [ ( x0, y2 ), ( x1, y3 ) ]
            cmds.append( svg_basic.make_stack( p0, p1, nbr_jl ) )
            cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            # if ng_lenseqs:
            #     p0,p1 = [ ( x0, y4 ), ( x1, y5 ) ]
            #     cmds.append( svg_basic.make_stack( p0, p1, ng_lenseq_jl ) )
            #     cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )

            # p0,p1 = [ ( x0, y6 ), ( x1, y7 ) ]
            # cmds.append( svg_basic.make_stack( p0, p1, ng_fwdseq_jl ) )
            # cmds.append( svg_basic.rectangle( p0, p1, "none", "black" ) )




        x0=x1+sep
        n=len(matched_tcrs)
        nplus=len(matched_tcrs_plus_nbrs)

        if not paper_figs and not paper_supp:
            font_family="Droid Sans Mono"
            expect_very_small = 0.001
            messages = [ '{} {} #tcrs= {}'.format(epitope,ab,len(tcrs)),
                         'chi-sq: {:.1f}'.format(chi_squared),
                         'motif: {}'.format(''.join(showmotif)),
                         'match-: {} {:.1f}%'.format(n,(100.0*n)/len(tcrs)),
                         'match+: {} {:.1f}%'.format(nplus,(100.0*nplus)/len(tcrs)),
                         'expect: {:.3f}'.format(max(expect_random,expect_nextgen)),
                         'enrich: {:.1f}'.format(float(n)/max([expect_random,expect_nextgen,expect_very_small])) ]
            if verbose:
                messages += [ '---',
                              '#naive_fixlen: {}'.format(num_ng_lenseqs),
                              '#naive: {} {:.1f}'.format(num_ng_fwdseqs,float(num_ng_fwdseqs)/n) ]

            for ii,text in enumerate(messages):
                p0 = ( x0, total_y + (ii+1)*fontsize ) ## lower left of text

                width = fontsize * len(text) * 0.6
                svg_width = max( svg_width, p0[0] + width + 50 )

                cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

        elif paper_supp:
            font_family="Droid Sans Mono"
            ## add some text
            ## epitope
            ## ab
            ## %matched
            ## see http://unicode.org/charts/PDF/U0370.pdf  for greek unicode
            ##
            greek_alpha = '&#x3b1;'
            greek_beta  = '&#x3b2;'
            greek_chi   = '&#x3c7;'
            alpha_beta = greek_alpha if ab == 'A' else greek_beta
            fontsize = 20
            p0 = [ xmargin + gene_stack_width/2.0 - 20, y3 + fontsize*0.9 ]
            text = 'V'+alpha_beta
            cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

            if False and ab == 'A':
                p0 = [ xmargin + 5.0, y_last - 5.0 ]
                text = epitope
                cmds.append( svg_basic.make_text( text, p0, 35.0, font_family=font_family ) )

            # p0 = [ xmargin + sep, y3 + 3*fontsize*0.9 ]
            # text = '{}, top'.format(epitope)
            # cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

            # p0 = [ xmargin + sep, y3 + 4*fontsize*0.9 ]
            # text = '{}-motif'.format(alpha_beta)
            # cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

            x0_for_J_gene_stack = xmargin+gene_stack_width+sep+pwm_stack_width+sep
            p0 = [ x0_for_J_gene_stack + gene_stack_width/2.0 - 20, y3 + fontsize*0.9 ]
            text = 'J'+alpha_beta
            cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

            xoff=2
            p0 = [ x0_for_J_gene_stack+xoff, y3 + 3*fontsize*0.9 ]
            text = 'Motif#{}'.format(num)
            cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

            p0 = [ x0_for_J_gene_stack+xoff, y3 + 4.5*fontsize*0.9 ]
            text = '{}2={:d}'.format(greek_chi,int(chi_squared))
            cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

            p0 = [ x0_for_J_gene_stack+xoff, y3 + 6.0*fontsize*0.9 ]
            text = 'f={:.1f}%'.format((100.0*nplus)/len(tcrs))
            cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family))



        else: #if paper_figs:
            font_family="Droid Sans Mono"
            ## add some text
            ## epitope
            ## ab
            ## %matched
            ## see http://unicode.org/charts/PDF/U0370.pdf  for greek unicode
            ##
            greek_alpha = '&#x3b1;'
            greek_beta  = '&#x3b2;'
            greek_chi   = '&#x3c7;'
            alpha_beta = greek_alpha if ab == 'A' else greek_beta
            fontsize = 20
            p0 = [ xmargin + gene_stack_width/2.0 - 20, y3 + fontsize*0.9 ]
            text = 'V'+alpha_beta
            cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

            if ab == 'A':
                p0 = [ xmargin + 5.0, y_last - 5.0 ]
                text = epitope
                cmds.append( svg_basic.make_text( text, p0, 35.0, font_family=font_family ) )

            # p0 = [ xmargin + sep, y3 + 3*fontsize*0.9 ]
            # text = '{}, top'.format(epitope)
            # cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

            # p0 = [ xmargin + sep, y3 + 4*fontsize*0.9 ]
            # text = '{}-motif'.format(alpha_beta)
            # cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

            x0_for_J_gene_stack = xmargin+gene_stack_width+sep+pwm_stack_width+sep
            p0 = [ x0_for_J_gene_stack + gene_stack_width/2.0 - 20, y3 + fontsize*0.9 ]
            text = 'J'+alpha_beta
            cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

            xoff=2
            p0 = [ x0_for_J_gene_stack+xoff, y3 + 3*fontsize*0.9 ]
            text = 'Motif#1'
            cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

            p0 = [ x0_for_J_gene_stack+xoff, y3 + 4.5*fontsize*0.9 ]
            text = '{}2={:d}'.format(greek_chi,int(chi_squared))
            cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family ) )

            p0 = [ x0_for_J_gene_stack+xoff, y3 + 6.0*fontsize*0.9 ]
            text = 'f={:.1f}%'.format((100.0*nplus)/len(tcrs))
            cmds.append( svg_basic.make_text( text, p0, fontsize, font_family=font_family))



            ## create our own little svg file with this guy
            svgfile = '{}_{}_{}_{}.svg'.format( outfile_prefix, epitope, ab, num )
            svg_basic.create_file( cmds, x1, y_last, svgfile, create_png=True )


        total_y = y_last + ypad

    if paper_supp:
        ## save all the recent set of cmds
        ab_cmds[ab] = cmds[:]
        ab_total_y[ab] = total_y

        ## create our own little svg file with this guy
        for ab in 'AB':
            if not ab_cmds[ab]:
                print 'no {} motifs for {}'.format(ab,epitope)
                continue
            svgfile = '{}_{}_{}.svg'.format( outfile_prefix, epitope, ab )
            svg_basic.create_file( ab_cmds[ab], x1, ab_total_y[ab], svgfile, create_png=True )


    motif_trees_info[ epitope ] = [ motif_dists, motif_infos, motif_positions, motif_matches ]

if paper_figs or paper_supp:
    exit() ## dont do the other junk ###################################################################################


if big:
    tag = 'ngbig'
else:
    tag = 'ngsmall'

prefix = '{}_motif_summary_{}'.format(outfile_prefix,tag)
svg_basic.create_file( cmds, svg_width, total_y+25, prefix+'.svg', create_png=True )

util.readme(prefix+'.png',"""This figure contains summary information on CDR3 motifs discovered by a simple
algorithm that basically looks for sequence patterns that occur more often in the epitope-specific datasets
than in a background set of sequences which is built by combining (1) random TCRs generated with a simple
model of the rearrangement process and (2) random subsets of next-gen unpaired data that is not selected for
epitope specificity. The challenge is to find motifs that are non-trivial in the sense that they do not just
come from germline sequence.
<br><br>
The main measure of significance I am using here is a simple "chi-squared" statistic which is basically <br>
<br>
 (observed-expected)*(observed-expected) / expected
<br><br>
Expected is judged based on two populations of background TCR sequences, one generated by a random model of the rearrangement process and one drawn from nextgen data. To be conservative, the value of expected used to calculate
significance is the maximum of the expected values based on these two background sets.
These background sets have the same V and J gene frequencies as the epitope specific set.
There is no correction for all the multiple testing going here, of course!
<br>
<br>
match- is the number of seed motif matches,<br>
match+ is the number of matched TCRs when I expand out to include nearby neighbors using the TCRdist distance measure
with a neighbor distance threshold of {:.2f}.
<br><br>
The match+ set are used to construct the sequence logos shown. The idea here is to account for the fact that the motif alphabet is not terribly expressive and motif discovery is not perfect, so we might be missing 'true' instances of the motif just using the initial seed motif (which is shown in "regular expression" syntax on the right):

<br><br>

^ = start of CDR3<br>
$ = end of CDR3<br>
. = any amino acid<br>
k = KR<br>
d = DE<br>
n = NQ<br>
s = ST<br>
f = FYWH<br>
a = AGSP<br>
v = VILM<br>
A = A<br>
C = C<br>
...<br>
Y = Y<br>
<br>
<br>
enrich is observed/expected
<br>
<br>
The color bars between the two logos indicate the inferred nucleotide source for the corresponding positions, as follows:<br>
<br>
light gray = V region<br>
black = D region<br>
dark gray = J region<br>
red = N insertions<br>
<br>
""".format( nbr_distance ))




### make some dendrograms
from scipy.cluster import hierarchy
from scipy.spatial import distance

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


## trick is that we don't know how big the figure should be
##
##


tree_infos = []

for epitope in sorted( motif_trees_info.keys() ):
    motif_dists, motif_infos, motif_positions, motif_matches = motif_trees_info[epitope]
    for chain in sorted( motif_infos.keys() ):
        if motif_infos[chain]:
            tree_infos.append( ( epitope, chain, len(motif_infos[chain] ) ) )

if not tree_infos:
    exit()


num_trees = len(tree_infos)
num_motifs = sum( (x[2] for x in tree_infos ) )

top_margin_inches = 0.5
bottom_margin_inches = 0.5
tree_vspacer_inches = 0.5

label_font_size = 8.0
single_motif_height_inches = label_font_size / 72.0 ## 72 points per inch, roughly

fig_height = top_margin_inches + bottom_margin_inches + (num_trees-1)*tree_vspacer_inches + num_motifs*single_motif_height_inches
fig_width = 8.0

print 'fig_height=',fig_height,'fig_width=',fig_width


left_margin_inches = 2.5
right_margin_inches = 0.5


fig = plt.figure(1,figsize=(fig_width,fig_height))

inches_from_top = top_margin_inches

for ( epitope, chain, num_motifs_ec ) in tree_infos:

    tcrs = all_tcrs[epitope]
    num_tcrs = len(tcrs)

    ## where do we start
    tree_height_inches = num_motifs_ec * single_motif_height_inches
    top = float( fig_height - inches_from_top ) / fig_height
    bottom = float( fig_height - ( inches_from_top + tree_height_inches ) ) / fig_height
    left = float( left_margin_inches ) / fig_width
    right = float( fig_width - right_margin_inches ) / fig_width
    rect = ( left, bottom, (right-left), (top-bottom) ) ## l,b,w,h
    ax = fig.add_axes( rect )

    #motif_dists, motif_infos, motif_positions, motif_matches = motif_trees_info[epitope]
    dists     = motif_trees_info[epitope][0][chain]
    infos     = motif_trees_info[epitope][1][chain]
    positions = motif_trees_info[epitope][2][chain]
    matches   = motif_trees_info[epitope][3][chain]

    assert len(infos) == num_motifs_ec
    assert len(dists) == num_motifs_ec
    assert len(positions) == num_motifs_ec

    D = np.zeros( (num_motifs_ec, num_motifs_ec ) )

    for i in range(num_motifs_ec):
        for j in range(num_motifs_ec):
            if j<i:
                continue
            elif j==i:
                D[i,i] = 0
            else:
                D[i,j] = dists[j][i]
                D[j,i] = dists[j][i]
                assert -1e-3 <= D[i,j] <= 1.001

    if num_motifs_ec>1:
        Z = hierarchy.average( distance.squareform(D,force='tovector') )

        hierarchy.dendrogram( Z, ax=ax, orientation='right',labels = [x[0] for x in infos],
                              leaf_font_size=label_font_size )

    inches_from_top += tree_height_inches + tree_vspacer_inches
    plt.title('{} {} {}'.format(epitope,chain,num_motifs_ec))


    ## try hacky clustering protocol with the distance matrix
    threshold = motifs_clustering_threshold
    min_cluster_size = 1
    N = num_motifs_ec

    all_nbrs = []
    for i in range(N):
        nbrs = []
        assert D[i,i] <= threshold
        for j in range(N):
            if D[i,j] <= threshold:
                nbrs.append( j )
        all_nbrs.append( nbrs )

    deleted = [False]*N

    centers = []
    all_members = []

    while True:
        clusterno = len(centers)

        best_nbr_count =0
        for i in range(N):
            if deleted[i]: continue
            nbr_count = 0
            for nbr in all_nbrs[i]:
                if not deleted[nbr]:
                    nbr_count+=1
            if nbr_count > best_nbr_count:
                best_nbr_count = nbr_count
                center = i

        if best_nbr_count < min_cluster_size:
            break

        centers.append( center )
        members = [center]

        print 'new-cluster:',epitope, chain, threshold, len(centers), best_nbr_count, infos[center][0]

        deleted[center] = True
        for nbr in all_nbrs[center]:
            if not deleted[nbr]:
                deleted[nbr] = True
                members.append( nbr )

        assert len(members) == best_nbr_count
        all_members.append( frozenset(members) )


    ## store member, center info in the infos/outl array
    assert len(centers) == len(all_members)
    for ii,(center, members) in enumerate( zip( centers, all_members ) ):
        for m in members:
            infos[m][1]['cluster_number'] = ii # 0-indexed
            infos[m][1]['is_cluster_center'] = 1 if m == center else 0

        ## work toward building a big pwm/consensus for this cluster of motifs
        center_positions = positions[center] ## dict from index to match_poslist
        center_info = infos[center][1]
        center_matches = matches[center]
        assert center_info['num_matches_with_nbrs'] == len(center_positions)
        assert center_info['num_matches_with_nbrs'] == len(center_matches)

        all_aa_counts = {}

        for m in members:
            m_positions = positions[m]
            m_matches = matches[m]
            assert infos[m][1]['num_matches_with_nbrs'] == len(m_positions)
            assert infos[m][1]['num_matches_with_nbrs'] == len(m_matches)
            expected = infos[m][1]['expected_fraction'] * num_tcrs
            my_weight = 1.0 if expected < 1 else 1.0/expected
            offsets = []
            for (ind,(mseq,posl)) in m_positions.iteritems():
                if ind in center_matches:
                    center_posl = center_positions[ind][1]
                    offsets.append( center_posl[0] - posl[0] ) ## add this to posl[0] to get center_posl[0]
            ## sanity check
            d = 1.0 - float( len( offsets ) )/(max(len(center_matches),len(m_matches)))
            assert abs( D[center,m] - d )<1e-3
            offsets.sort()
            ## OK, now we have a consensus offset
            offset = int(get_median( offsets ))
            for ind,(mseq,posl) in m_positions.iteritems():
                tcr = tcrs[ind]
                cdr3 = ( tcr[4] if chain == 'A' else tcr[5] )
                assert mseq == cdr3[min(posl):max(posl)+1]
                ## where would center_posl0 most likely match to?
                center_posl0 = posl[0] + offset
                for p in posl:
                    aa = cdr3[p]
                    center_ipos = p - center_posl0
                    if center_ipos not in all_aa_counts: all_aa_counts[center_ipos] = {}
                    all_aa_counts[center_ipos][aa] = all_aa_counts[center_ipos].get(aa,0)+my_weight


        ## now figure out which positions we should include in the logo...
        l = [ (sum( all_aa_counts[x].values() ), x ) for x in all_aa_counts ]
        l.sort()
        l.reverse()
        print 'all_aa_counts by position:',l
        max_total = l[0][0]
        posl = all_aa_counts.keys()[:]
        posl.sort()
        consensus = ''
        started = False
        finished = False
        min_total = 0.75 * max_total
        for pos in posl:
            total = sum( all_aa_counts[pos].values() )
            if total < min_total:
                if started:
                    finished = True
            else:
                assert not finished ## shouldnt go out and then come back in
                started = True
                consensus += get_amino_acid_consensus_character( all_aa_counts[pos] )

        ## now we've got the consensus
        for m in members:
            infos[m][1]['cluster_consensus'] = consensus



pngfile = '{}_motif_trees.png'.format(outfile_prefix)

print 'making:',pngfile
plt.savefig(pngfile)


## write out all the info
tsv_outfile = '{}_motifs.tsv'.format(outfile_prefix )
out = open(tsv_outfile,'w')

tsv_outfields = ['id','epitope','chain',
                 'epitope_num','chi_squared','overlap',
                 'num_matches','num_matches_with_nbrs',
                 'expected_fraction','expected_num_matches',
                 'showmotif','matches','matches_with_nbrs','matches_with_nbrs_consensus',
                 'cluster_number','is_cluster_center','cluster_consensus']

out.write('\t'.join( tsv_outfields )+'\n' )
for ( epitope, chain, num_motifs_ec ) in tree_infos:
    #motif_dists, motif_infos = motif_trees_info[epitope]
    dists = motif_trees_info[epitope][0][chain]
    infos = motif_trees_info[epitope][1][chain]

    for (label,tsv_outl) in infos:
        out.write( make_tsv_line( tsv_outl, tsv_outfields )+'\n' )
out.close()

