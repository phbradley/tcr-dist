## look at frequencies and enrichments of genes
##
## what makes this a little tricky is that gene assignments can be ambiguous: sequence reads
## are often too short to uniquely define the gene, especially in mouse V-alpha where there are
## highly sequence-similar genes.

from basic import *
#from util import get_rep, get_mm1_rep, get_mm1_rep_gene_for_counting
import util
from tcr_rearrangement_new import all_countrep_pseudoprobs
from all_genes import all_genes

with Parser(locals()) as p:
    #p.str('args').unspecified_default().multiple().required()
    p.str('clones_file').required()
    p.str('organism').required()


summary_tsvfile = clones_file[:-4] + '_JS_divergence.tsv'
summary_fields = ['epitope']+[x+'_jsd_normed' for x in segtypes_lowercase ]

probs_tsvfile = clones_file[:-4] + '_gene_probs.tsv'
probs_fields = ['epitope','gene','label_prob','jsd_prob','jsd_prob_enrich','pseudoprob','pseudoprob_enrich']

def get_freq_from_tuple_counts_and_bias( counts, bias ):
    reps = set()

    total = sum( counts.values() )
    for t in counts:
        for rep in t:
            reps.add(rep)
    nreps = len(reps)

    single_counts = dict( zip(reps, [0.]*nreps ) )

    for t,count in counts.iteritems():
        tfreqs = [ bias.get(x,0.) for x in t ]
        tfreqs_total = sum(tfreqs)
        if tfreqs_total==0:
            tfreqs = [ 1.0 ]*len(t)
            tfreqs_total = float(len(t))
        for rep,freq in zip(t,tfreqs):
            single_counts[ rep ] += count * freq / tfreqs_total

    rep_freq = {}
    for rep in reps:
        rep_freq[rep] = float( single_counts[rep] ) / total

    return rep_freq


def get_baseline_tuple_count_frequencies( counts ):
    reps = set()
    for t in counts:
        for rep in t:
            reps.add(rep)
    nreps = len(reps)

    bias = dict( zip( reps, [1.0/nreps]*nreps ) )

    return get_freq_from_tuple_counts_and_bias( counts, bias )


def js_divergence( P, Q ):
    total=0.
    for rep in frozenset( P.keys() + Q.keys() ):
        p = P.get(rep,0.0)
        q = Q.get(rep,0.0)
        m = 0.5*(p+q)
        if p: total += 0.5 * p * math.log( p / m, 2. )
        if q: total += 0.5 * q * math.log( q / m, 2. )
    return total


def shannon_entropy( P ):
    total=0.
    for rep,p in P.iteritems():
        if p:
            total -= p * math.log( p, 2.0 )
    return total


def get_jsd_normed( P, Q ):
    return  2.0 * js_divergence(P,Q) / ( shannon_entropy(P) + shannon_entropy(Q) )

# def get_relent_traditional( freq, base_freq ):
#     relent=0
#     for rep,prob in freq.iteritems():
#         if prob==0: continue
#         naive_prob = base_freq.get(rep,0)
#         if naive_prob==0:
#             print 'whoah:: rep={} prob={} naive_prob={}'.format(rep,prob,naive_prob)
#             continue
#         assert naive_prob>0
#         relent += prob * math.log(prob/naive_prob,2)
#     return relent


# def get_relent( freq, base_freq ):
#     return get_jsd_normed( freq, base_freq ) ## hacking!

# since the files don't have info on chain/region, set up a mapping (note that the tuples are of genes not alleles)
gene2segtype = {}
for id,g in all_genes[organism].iteritems():
    gene = id[:id.index('*')]
    segtype = g.region + g.chain
    if gene in gene2segtype:
        assert gene2segtype[gene] == segtype
    else:
        gene2segtype[gene] = segtype


## read the background gene-tuple counts
all_background_tuple_counts = {}

tuplecountsfile = path_to_current_db_files()+'/nextgen_tuple_counts_v2_{}_max10M.log'.format(organism)

if exists( tuplecountsfile ):
    Log('reading {}'.format( tuplecountsfile ))
    for line in open( tuplecountsfile,'r'):
        if line.startswith('TUPLE_COUNT'):
            l = line.split()
            count = int(l[1] )
            tup = tuple( sorted( l[2].split(',') ) )
            segtype = gene2segtype[ tup[0] ]
            #segtype = tup[0][3] + tup[0][2] ## go from 'TRAV*' to 'VA'
            assert segtype in segtypes_uppercase
            if segtype not in all_background_tuple_counts:
                all_background_tuple_counts[segtype] = {}
            logfile = l[-1]
            if logfile not in all_background_tuple_counts[segtype]:
                all_background_tuple_counts[segtype][logfile] = {}
            assert tup not in all_background_tuple_counts[segtype][logfile] ## should not be any repeats in the output
            all_background_tuple_counts[segtype][logfile][ tup ] = count
else:
    Log('WARNING: making up fake tuple counts since dbfile ({}) is missing'.format(tuplecountsfile))
    logfile = 'fake_logfile'
    for id,g in all_genes[organism].iteritems():
        segtype = g.region + g.chain
        if segtype in segtypes_uppercase:
            if segtype not in all_background_tuple_counts:
                all_background_tuple_counts[segtype] = {logfile:{}}
            tup = tuple( util.get_mm1_rep_gene_for_counting( id, organism ) )
            all_background_tuple_counts[segtype][logfile][ tup ] = 1 ## flat counts


## here we compute J-S divergences of the gene distributions to background
## the trick is that both our gene assigments and the bg gene assigments are ambiguous, and the degree of
## ambiguity depends on experiment-specific things like sequence read length. So to be conservative, when
## comparing our epitope-specific counts to the background counts we collapse ambiguity so as to minimize
## the divergence between the two distributions
##

## read the tcrs from the clones file
all_tcrs = parse_tsv_file( clones_file, ['epitope'], [], True )



print 'making:',summary_tsvfile
out_summary = open( summary_tsvfile, 'w')
out_summary.write( '\t'.join( summary_fields )+'\n' )

print 'making:',probs_tsvfile
out_probs = open( probs_tsvfile, 'w')
out_probs.write( '\t'.join( probs_fields )+'\n' )

for epitope,tcrs in all_tcrs.iteritems():
    ## this collapses the ambiguity of segment assignment by going with the most popular one in the repertoire
    ## only doing this here for completeness in the tsv output
    ## not used for JS-DIV calculation
    util.assign_label_reps_and_colors_based_on_most_common_genes_in_repertoire( tcrs, organism )

    out_summary_l = {'epitope':epitope}


    for segtype in segtypes_uppercase:
        region = segtype[0]
        chain = segtype[1]
        assert region in 'VJ' and chain in 'AB'

        tcr_counts = {}
        num_tcrs = len(tcrs)

        label_rep_counts = {}
        tcr_pseudoprob_counts = {}

        for tcr in tcrs:
            countreps = tuple( sorted( tcr['{}_countreps'.format(segtype.lower())].split(';') ) )
            tcr_counts[ countreps ] = tcr_counts.get( countreps, 0)+1
            for rep in countreps:
                tcr_pseudoprob_counts[rep] = tcr_pseudoprob_counts.get(rep,0)+1

            label_rep = tcr['{}_label_rep'.format(segtype.lower()) ]
            label_rep_counts[label_rep] = label_rep_counts.get(label_rep,0)+1

        best_jsd = 100.
        best_tcr_freq = {}
        best_bg_freq = {}

        for logfile in all_background_tuple_counts[ segtype ]:
            #print epitope,segtype,logfile

            bg_counts = all_background_tuple_counts[ segtype ][ logfile ]


            bg_freq  = get_baseline_tuple_count_frequencies( bg_counts )
            tcr_freq = get_baseline_tuple_count_frequencies( tcr_counts )

            jsd = get_jsd_normed( tcr_freq, bg_freq )
            #print 'jsd0:',jsd

            for niter in range(1,6):
                new_bg_freq  = get_freq_from_tuple_counts_and_bias(  bg_counts, tcr_freq )
                new_tcr_freq = get_freq_from_tuple_counts_and_bias( tcr_counts, new_bg_freq )

                tcr_freq = new_tcr_freq
                bg_freq = new_bg_freq

                jsd = get_jsd_normed( tcr_freq, bg_freq )

                #print 'jsd{}: {:9.6f}'.format(niter,jsd),

                if jsd<best_jsd:
                    #print '******'
                    best_logfile = logfile
                    best_jsd = jsd
                    best_bg_freq = bg_freq
                    best_tcr_freq = tcr_freq
                else:
                    pass
                    #print


        out_summary_l[ segtype.lower()+'_jsd_normed' ] = '{:9.3f}'.format( best_jsd )

        # print 'best_jsd: {:9.3f} {} {} {}'\
        #     .format( best_jsd, epitope, segtype, best_logfile )

        ## let's look at enrichments of the genes
        bg_freq, tcr_freq = best_bg_freq, best_tcr_freq

        l = [ (y,x) for x,y in tcr_freq.iteritems() ]
        l.sort()
        l.reverse()

        countreps = sorted( set( g.count_rep for g in all_genes[organism].values() \
                                 if g.chain == chain and g.region == region ) )

        num_tcrs = len(tcrs)

        for rep in countreps:
            label_prob = float( label_rep_counts.get(rep,0) )/num_tcrs
            jsd_prob = tcr_freq.get(rep,0.)
            bg_prob = bg_freq.get(rep,0.)
            jsd_prob_enrich = 1.0 if bg_prob==0. else jsd_prob / bg_prob
            pseudoprob = float( tcr_pseudoprob_counts.get(rep,0))/num_tcrs
            bg_pseudoprob = all_countrep_pseudoprobs[ organism ][ chain ][ region ][ rep ]
            pseudoprob_enrich = 1 if bg_pseudoprob==0. else pseudoprob / bg_pseudoprob
            outl = { 'epitope':epitope,
                     'gene':rep,
                     'label_prob':'{:.6f}'.format(label_prob),
                     'jsd_prob':'{:.6f}'.format(jsd_prob),
                     'jsd_prob_enrich':'{:9.3f}'.format(jsd_prob_enrich),
                     'pseudoprob':'{:.6f}'.format(pseudoprob),
                     'pseudoprob_enrich':'{:9.3f}'.format(pseudoprob_enrich),
            }

            out_probs.write( make_tsv_line( outl, probs_fields ) + '\n' )
    out_summary.write( make_tsv_line( out_summary_l, summary_fields ) + '\n' )


out_probs.close()
out_summary.close()





