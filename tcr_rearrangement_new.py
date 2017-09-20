from basic import *
import sys
from all_genes import all_genes
from operator import add
import util

## these guys are for external sharing, as well as
## get_alpha_trim_probs
## get_beta_trim_probs
##
all_trim_probs = {}
all_countrep_pseudoprobs = {} ## pseudoprobs because they do not sum to 1.0, due to ambiguity of gene assignments
all_trbd_nucseq = {}

verbose = ( __name__ == '__main__' )

organism_chains_with_missing_probs = []

for organism in all_genes:
    d_ids = sorted( ( id for id,g in all_genes[organism].iteritems() if g.region == 'D' and g.chain == 'B' ) )
    all_trbd_nucseq[ organism ] = dict( ( ( d_ids.index(id)+1, all_genes[organism][id].nucseq ) for id in d_ids ) )


for organism in all_genes:
    rep_freq_files = {}
    trim_prob_lines = {}
    for chain in 'AB':
        probs_files = glob('{}/probs_files_{}_{}/*'.format(path_to_current_db_files(),organism,chain))
        trim_prob_lines[chain] = []
        rep_freq_files[chain] = probs_files
        if not probs_files:
            Log('WARNING:: tcr_rearrangement_new: no probs files for {} {}'.format(organism,chain))
            organism_chains_with_missing_probs.append( ( organism, chain ) )
        else:
            for file in probs_files:
                trim_prob_lines[chain].extend( popen('grep "^PROB_{}_" {}'.format( chain, file )).readlines())



    trim_probs = {}
    if trim_prob_lines['A']:
        for line in trim_prob_lines['A']:
            l = line.split()
            if not l:continue
            assert l[0].startswith('PROB_A')
            tag = l[0][5:]
            vals = map(float,l[1:])
            if tag not in trim_probs:
                trim_probs[tag] = []
            trim_probs[tag].append( dict( zip( range(len(vals)), vals ) ) )

    if trim_prob_lines['B']:
        for line in trim_prob_lines['B']:
            l = line.split()
            if not l:continue
            assert l[0].startswith('PROB_B')
            tag = l[0][5:]
            if tag not in trim_probs:
                trim_probs[tag] = []
            assert len(l)%2==1
            num_vals = (len(l)-1)/2
            D = {}
            for i in range(num_vals):
                assert l[2*i+1][-1] == ':'
                key = l[2*i+1][:-1]
                if ',' in key:
                    key = tuple(map(int,key.split(',')))
                else:
                    key = int(key)
                #trim_probs[tag][key] = float(l[2*i+2])
                D[key] = float(l[2*i+2])
            trim_probs[tag].append( D )

    ## now average to get a single prob distn
    for tag in trim_probs.keys():
        Dlist = trim_probs[tag]
        ks = sorted( set( reduce( add, [D.keys() for D in Dlist] ) ) )
        newD = {}
        for k in ks:
            newD[k] = sum( ( D.get(k,0.) for D in Dlist ) )
        total = sum(newD.values())
        if verbose:
            print 'tag:',organism,tag,'len(Dlist):',len(Dlist),'total:',total
        for k in ks:
            newD[k] /= total
        trim_probs[tag] = newD


    if trim_prob_lines['B']:
        ## fake probability for total trimming of the D gene
        for did,nucseq in all_trbd_nucseq[organism].iteritems():
            trimtag = 'B_D{}_d01_trim'.format(did)
            prob_trim_all_but_1 = 0.0
            for d0_trim in range(len(nucseq)):
                d1_trim = (len(nucseq)-1)-d0_trim
                assert d0_trim + d1_trim == len(nucseq)-1
                prob_trim_all_but_1 += trim_probs[trimtag].get((d0_trim,d1_trim),0)
            prob_trim_all = 0.0
            for d0_trim in range(len(nucseq)+1):
                d1_trim = (len(nucseq))-d0_trim
                prob_trim_all += trim_probs[trimtag].get((d0_trim,d1_trim),0)
            assert prob_trim_all <1e-6
            #print 'old_prob_trim_all:',prob_trim_all,'prob_trim_all_but_1:',prob_trim_all_but_1,'D',did
            new_prob_trim_all = 0.75 * prob_trim_all_but_1
            for d0_trim in range(len(nucseq)+1):
                d1_trim = (len(nucseq))-d0_trim
                if d0_trim == 0:
                    #print 'new_prob_trim_all:',new_prob_trim_all
                    trim_probs[trimtag][ (d0_trim,d1_trim) ] = new_prob_trim_all ## concentrate all here
                else:
                    trim_probs[trimtag][ (d0_trim,d1_trim) ] = 0.0
            total = sum( trim_probs[trimtag].values())
            for k in trim_probs[trimtag]:
                trim_probs[trimtag][k] /= total
            assert abs( 1.0 - sum( trim_probs[trimtag].values()) ) < 1e-6



        beta_prob_tags_single = ['v_trim','j_trim','vd_insert','dj_insert']
        for tag in beta_prob_tags_single:
            tags = [ 'B_D{}_{}'.format(x,tag) for x in all_trbd_nucseq[organism] ]
            #tag1 = 'B_D1_{}'.format(tag)
            #tag2 = 'B_D2_{}'.format(tag)
            avgtag = 'B_{}'.format(tag)
            trim_probs[avgtag] = {}
            ks = sorted( set( reduce( add, [ trim_probs[x].keys() for x in tags ] ) ) )
            #print organism,tag,ks
            for k in ks:
                trim_probs[avgtag][k] = sum( ( trim_probs[x].get(k,0) for x in tags ) ) / float(len(tags))

    all_countrep_pseudoprobs[organism] = { 'A': { 'V':{}, 'J': {} },
                                           'B': { 'V':{}, 'J': {} } }
    for ab in rep_freq_files:
        files = rep_freq_files[ab]
        for vj in 'VJ':
            probs ={}
            if files:
                for file in files:
                    assert exists(file)
                    for line in popen('grep "^{}{}_COUNTREP_FREQ" {}'.format(ab,vj,file)):
                        l = line.split()
                        assert len(l) == 3
                        nonuniq_freq = float( l[1] ) / 100.0 ## now from 0 to 1
                        rep = l[2]
                        #assert rep[2:4] == ab+vj
                        if rep not in probs:probs[rep] = []
                        probs[rep].append( nonuniq_freq )
            else:
                #make up some fake probs
                #ids = [ x for x,g in all_genes[organism].iteritems() if g.chain == ab and g.region == vj ]
                reps = set( g.count_rep for g in all_genes[organism].values() if g.chain == ab and g.region == vj )
                flatprob = 1.0/len(reps)
                for rep in reps:
                    probs[rep] = [ flatprob ]

            avg_probs = {}
            for rep in probs:
                vals = probs[rep] + [0.0]*(len(files) - len(probs[rep]) )
                if len(vals) == 2:
                    avg_probs[rep] = sum( vals )/2.0
                else:
                    #assert len(vals) == 3 ## hack
                    avg_probs[rep] = get_median( vals)

            ## probs may have gone slightly below 1.0 due to combination of multiple datasets
            total = min( 1.0, sum( avg_probs.values() ) ) ##only increase probabilities...
            if verbose:
                print 'countrep_pseudoprobs total {:9.6f} actual_sum {:9.6f} {}{} {}'\
                    .format(total, sum(avg_probs.values()), vj, ab, organism )
            all_countreps = sorted( set( g.count_rep for g in all_genes[organism].values()
                                         if g.chain == ab and g.region == vj ) )
            #for rep in probs:
            for rep in all_countreps:
                if rep not in avg_probs and verbose:
                    print 'WARNING: rep with 0.0 pseudoprob: {} {}'.format(organism,rep)
                #countrep_probs[rep] = avg_probs[rep] / total
                prob = avg_probs.get(rep,0.)/total
                all_countrep_pseudoprobs[organism][ab][vj][rep] = prob
                if verbose: # __name__ == '__main__' and len(sys.argv) == 1:
                    print 'countrep_pseudoprobs: %12.6f %s %s'%(100.0*prob, organism, rep )

    ## normalize trim_probs
    for tag,probs in trim_probs.iteritems():
        if type(probs) == type({}):
            total = sum( probs.values())
            assert abs(1.0-total)<1e-2
            #print 'normalize trim_probs:',tag,total
            for k in probs:
                probs[k] = probs[k] / total
        else:
            assert False
            assert type(probs) == type([])
            total = sum( probs )
            assert abs(1.0-total)<1e-2
            #print 'normalize trim_probs:',tag,total
            for i in range(len(probs)):
                probs[i] = probs[i]/total

    all_trim_probs[organism] = trim_probs


def get_alpha_trim_probs( organism, v_trim, j_trim, vj_insert ):
    if (organism,'A') in organism_chains_with_missing_probs:
        return 1.0
    total_prob = 1.0
    for ( val, tag ) in zip( [v_trim, j_trim, vj_insert], ['A_v_trim','A_j_trim','A_vj_insert'] ):
        probs = all_trim_probs[organism][tag]
        if val >= len(probs):
            return 0.0
        total_prob *= probs[val]
    return total_prob

def get_beta_trim_probs( organism, d_id, v_trim, d0_trim, d1_trim, j_trim, vd_insert, dj_insert ): ## work in progress
    if (organism,'B') in organism_chains_with_missing_probs:
        return 1.0
    assert d_id in all_trbd_nucseq[organism]
    dd = (d0_trim,d1_trim)
    d_trim_tag = 'B_D{}_d01_trim'.format(d_id)
    total_prob = all_trim_probs[organism][d_trim_tag].get(dd,0)
    #total_prob = trim_probs[d_trim_tag][dd] ## what about full trims?? will get an error
    for ( val, tag ) in zip( [v_trim, j_trim, vd_insert, dj_insert], beta_prob_tags_single ):
        probs = all_trim_probs[organism]['B_'+tag] ## a dictionary for beta (a list for alpha)
        if val not in probs:
            return 0.0
        total_prob *= probs[val]
    return total_prob


def probs_data_exist( organism, chain ):
    assert chain in 'AB'
    return (organism,chain) not in organism_chains_with_missing_probs

