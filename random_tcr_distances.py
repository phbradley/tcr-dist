import sys
from basic import *
import tcr_distances
import cdr3s_human
import parse_tsv
from paths import path_to_db

with Parser(locals()) as p:
    p.str('organism').required()
    p.str('clones_file').required()
    p.int('nrandom').default(10000)
    p.multiword('nbrdist_percentiles').cast(lambda x: [int(val) for val in x.split()]).default("5 10 25")#-1 -3 -5 5 10 25")
    p.flag('constant_seed')
    p.str('distance_params')

if constant_seed: random.seed(1)

distance_params = tcr_distances.DistanceParams( config_string = distance_params )


## read the epitope-specific TCRs
all_tcrs = parse_tsv.parse_tsv_file( clones_file, ['epitope'], ['va_genes','vb_genes','cdr3a','cdr3b'] )

## read the random tcrs
random_tcrs = []
random_info = []

if organism == 'mouse':
    random_chains = {}
    for ab in 'AB':
        random_chains[ab] = []
        random_tcrs_file = '{}/new_nextgen_chains_{}_{}.tsv'.format(path_to_db,organism,ab)
        assert exists( random_tcrs_file )
        Log('reading {}'.format( random_tcrs_file ))
        header = []
        for line in open(random_tcrs_file,'r'):
            l = line[:-1].split('\t')
            if not header:
                header = l[:]
                assert header == ['v_reps','j_reps','cdr3','cdr3_nucseq']
            else:
                random_chains[ab].append( ( frozenset( l[0].split(',') ), l[2] ) ) ## (v_reps, cdr3)
        random.shuffle( random_chains[ab] )
        assert len(random_chains[ab])>=nrandom
    for a,b in zip( random_chains['A'][:nrandom], random_chains['B'][:nrandom] ):
        random_tcrs.append( ( a[0],b[0],a[1],b[1] ) ) #va_reps, vb_reps, cdr3a, cdr3b
        random_info.append( { 'va_reps': ';'.join( a[0] ),
                              'vb_reps': ';'.join( b[0] ),
                              'cdr3a':a[1],
                              'cdr3b':b[1] } )

else:
    ##
    human_paired_tcrs_file = path_to_db+'/randhuman_paired.tsv'
    assert exists( human_paired_tcrs_file )
    rtcrs = parse_tsv.parse_tsv_file( human_paired_tcrs_file, [], ['va_reps','vb_reps','cdr3a','cdr3b'] )
    for l in rtcrs:
        random_tcrs.append( ( frozenset( l[0].split(';') ), frozenset( l[1].split(';') ), l[2], l[3] ) )
        random_info.append( { 'va_reps': l[0], 'vb_reps': l[1], 'cdr3a': l[2], 'cdr3b': l[3] } )


print 'precomputing v-region distances'
rep_dists = tcr_distances.compute_all_v_region_distances( organism, distance_params )
print 'done precomputing v-region distances'


for epitope in all_tcrs:
    tcrs = []
    for l in all_tcrs[epitope]:
        ## ( va_reps, vb_reps, cdr3a, cdr3b )
        tcrs.append( ( frozenset( [cdr3s_human.all_loopseq_representative[ organism ][x] for x in l[0].split(';') ] ),
                       frozenset( [cdr3s_human.all_loopseq_representative[ organism ][x] for x in l[1].split(';') ] ),
                       l[2], l[3] ) )

    print 'num_tcrs:',epitope,len(tcrs)

    for chains in ['A','B','AB']:
        for ii, rtcr in enumerate( random_tcrs ):
            edists = []
            for etcr in tcrs:
                edists.append( tcr_distances.compute_distance( etcr, rtcr, chains, rep_dists, distance_params ) )

            edists.sort()
            for nbrdist_percentile in nbrdist_percentiles:
                nbrdist = tcr_distances.sort_and_compute_nbrdist_from_distances( edists, nbrdist_percentile, dont_sort=True )
                wtd_nbrdist = tcr_distances.sort_and_compute_weighted_nbrdist_from_distances( edists, nbrdist_percentile,
                                                                                      dont_sort=True)
                random_info[ii]['{}_{}_nbrdist{}'.format(epitope,chains,nbrdist_percentile)] = '{:.3f}'.format(nbrdist)
                random_info[ii]['{}_{}_wtd_nbrdist{}'.format(epitope,chains,nbrdist_percentile)] = '{:.3f}'.format(wtd_nbrdist)

            # for distance_sdev in distance_sdevs:
            #     sdev = float(distance_sdev) * len(chains)
            #     density = 0.
            #     for d in edists:
            #         density += math.exp( -1.0 * ( d/sdev)**2 )
            #     density /= len(edists)
            #     random_info[ii]['{}_{}_dens{}'.format(epitope,chains,distance_sdev)] = density

## write the output
outfields = ['va_reps','vb_reps','cdr3a','cdr3b' ]

for epitope in all_tcrs:
    for chains in ['A','B','AB']:
        for nbrdist_percentile in nbrdist_percentiles:
            outfields.append( '{}_{}_nbrdist{}'.format(epitope,chains,nbrdist_percentile) )
            outfields.append( '{}_{}_wtd_nbrdist{}'.format(epitope,chains,nbrdist_percentile) )
        # for distance_sdev in distance_sdevs:
        #     outfields.append( '{}_{}_dens{}'.format(epitope,chains,distance_sdev) )


outfile = '{}_random_nbrdists.tsv'.format(clones_file[:-4])
out =open( outfile,'w')
out.write('\t'.join( outfields )+'\n' )

for outl in random_info:
    out.write(make_tsv_line(outl,outfields)+'\n' )
out.close()
