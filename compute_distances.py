## this script will compute distances between all the tcrs and write out a distance matrix
## the order of rows and columns in the distance matrix will match that in the input file
## it will also compute the rank scores for each tcr with respect to all the epitopes present
##

from basic import *
from util import get_top_genes
from all_genes import all_genes
#from scipy import stats, optimize
import tcr_distances

#import matplotlib
#if make_png: matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#import numpy as np

with Parser(locals()) as p:
    p.str('clones_file') ## tsv
    p.str('distfile_prefix') ## will be white-space separated like used for trees code
    p.str('outfile') ## tsv
    p.str('organism')
    p.flag('intrasubject_nbrdists')
    p.multiword('clones_files').cast(lambda x:x.split())
    p.multiword('epitope_prefixes').cast(lambda x:x.split())
    p.multiword('nbrdist_percentiles').cast(lambda x: [int(val) for val in x.split()]).default("5 10 25")#"5 10 25 -1 -5 -10")
    p.str('distance_params')

#internal legacy hack
new_nbrdists = not intrasubject_nbrdists


distance_params = tcr_distances.DistanceParams( config_string = distance_params )


if not clones_files:
    assert clones_file
    clones_files = [clones_file]

if not epitope_prefixes:
    epitope_prefixes = ['']*len(clones_files)

assert len(epitope_prefixes) == len(clones_files)

print 'precomputing v-region distances'
rep_dists = tcr_distances.compute_all_v_region_distances( organism, distance_params )
print 'done precomputing v-region distances'

if not distfile_prefix:
    assert len(clones_files) == 1
    distfile_prefix = clones_files[0][:-4]

if not outfile:
    assert len(clones_files) == 1
    outfile = '{}_nbrdists.tsv'.format(clones_files[0][:-4] )



all_tcrs = []
all_info = []

all_infields = []

for clones_file, epitope_prefix in zip(clones_files,epitope_prefixes):
    infields = []
    for line in open( clones_file,'rU'):
        if not infields:
            if line[0] == '#':
                infields = line[1:-1].split('\t')
            else:
                infields = line[:-1].split('\t')
            all_infields.append( infields )
            continue
        assert infields

        l = parse_tsv_line( line[:-1], infields )

        assert 'epitope' in l
        l['epitope'] = epitope_prefix + l['epitope']

        va_genes = l['va_genes'].split(';')
        vb_genes = l['vb_genes'].split(';')

        va_reps = frozenset( ( all_genes[organism][x].rep for x in va_genes ))
        vb_reps = frozenset( ( all_genes[organism][x].rep for x in vb_genes ))

        all_info.append( l )

        tcr = [ va_reps, vb_reps, l['cdr3a'], l['cdr3b'] ]

        all_tcrs.append( tcr )


## we will add new fields here

outfields = []
for f in all_infields[0]:
    common = True
    for infields in all_infields[1:]:
        if f not in infields:
            common = False
            break
    if common:
        outfields.append( f )



for chains in ['A','B','AB']:
    epitopes = list( set( [x['epitope'] for x in all_info ] ) )
    epitopes.sort()

    subjects = list( set( [x['subject'] for x in all_info ] ) )

    if len(subjects)==1 and new_nbrdists:
        Log("""WARNING:: compute_distances.py: It looks like there's only a single subject...
        Now including intra-subject distances in computing NNdistance values.
        This just means that any classification results may be a little too rosy""")
        new_nbrdists = False

    ## for computing rank scores
    epitope_self_nbrdists = {}
    for wtd in ['','wtd_']:
        epitope_self_nbrdists[wtd] = {}
        for e in epitopes:
            epitope_self_nbrdists[wtd][e] = {}
            for p in nbrdist_percentiles:
                epitope_self_nbrdists[wtd][e][p] = []

    # epitope_self_densities = {}
    # for e in epitopes:
    #     epitope_self_densities[e] = {}
    #     for sdev in distance_sdevs:
    #         epitope_self_densities[e][sdev] = []

    all_clone_ids = [x['clone_id'] for x in all_info]

    ## make a distances file
    out = open('{}_{}.dist'.format(distfile_prefix,chains),'w')
    #out.write('\t'.join(['clone_id']+all_clone_ids)+'\n')

    e_out = {}
    for e in epitopes:
        e_out[e] = open('{}_{}_{}.dist'.format(distfile_prefix,chains,e),'w')

    total_tcrs = len(all_tcrs)

    for ii,t1 in enumerate(all_tcrs):
        if not ii%100: Log('computing {} distances {} {}'.format(chains,ii,total_tcrs))
        epitope_distances = {}
        for e in epitopes:
            epitope_distances[e] = []

        e1 = all_info[ii]['epitope']
        m1 = all_info[ii]['subject']
        clone_id1 = all_info[ii]['clone_id']

        dists=[]
        my_epitope_dists = []
        for jj,t2 in enumerate(all_tcrs):
            dist = tcr_distances.compute_distance( t1,t2,chains,rep_dists,distance_params)
            dists.append(dist)

            e2 = all_info[jj]['epitope']
            m2 = all_info[jj]['subject']
            if ii!=jj: ## dont include self distance in the nbrdists, rank scores
                if new_nbrdists:
                    if m1 != m2:
                        epitope_distances[e2].append(dist)
                else:
                    epitope_distances[e2].append(dist)
            if e2==e1:
                my_epitope_dists.append(dist)

        ## write out a new line in the distance file
        assert len(dists) == len(all_clone_ids)

        out.write(' '.join( [clone_id1]+[ '{:.3f}'.format(x) for x in dists] )+'\n' )

        e_out[e1].write(' '.join( [clone_id1]+[ '{:.3f}'.format(x) for x in my_epitope_dists] )+'\n' )

        ## get an nbrdist
        for e,edists in epitope_distances.iteritems():
            edists.sort()
            for nbrdist_percentile in nbrdist_percentiles:
                if edists:
                    nbrdist = tcr_distances.sort_and_compute_nbrdist_from_distances( edists, nbrdist_percentile, dont_sort=True )
                    wtd_nbrdist = tcr_distances.sort_and_compute_weighted_nbrdist_from_distances( edists, nbrdist_percentile, dont_sort=True)
                else:
                    nbrdist = 0.0
                    wtd_nbrdist = 0.0
                all_info[ii]['{}_{}_nbrdist{}'.format(e,chains,nbrdist_percentile)] = nbrdist
                all_info[ii]['{}_{}_wtd_nbrdist{}'.format(e,chains,nbrdist_percentile)] = wtd_nbrdist
                if e == e1:
                    epitope_self_nbrdists[''][e][nbrdist_percentile].append( nbrdist )
                    epitope_self_nbrdists['wtd_'][e][nbrdist_percentile].append( wtd_nbrdist )

            # ## now also compute a gaussian weighted density
            # for distance_sdev in distance_sdevs:
            #     sdev = float(distance_sdev) * len(chains)
            #     density = 0.
            #     for d in edists:
            #         density += math.exp( -1.0 * ( d/sdev)**2 )
            #     density /= len(edists)
            #     all_info[ii]['{}_{}_dens{}'.format(e,chains,distance_sdev)] = density
            #     if e == e1:
            #         epitope_self_densities[e][distance_sdev].append( density )

    out.close()
    for o in e_out.values(): o.close()


    ## record new fields in outfields
    for wtd in ['','wtd_']:
        for nbrdist_percentile in nbrdist_percentiles:
            for suffix in ['_{}_{}nbrdist{}'.format(chains,wtd,nbrdist_percentile),
                           '_{}_{}nbrdist{}rank'.format(chains,wtd,nbrdist_percentile) ]:
                for epitope in epitopes:
                    outfields.append( epitope+suffix)
    # for distance_sdev in distance_sdevs:
    #     for suffix in ['_{}_dens{}'.format(chains,distance_sdev),
    #                    '_{}_rdens{}'.format(chains,distance_sdev) ]:
    #         for epitope in epitopes:
    #             outfields.append( epitope+suffix)

    ## compute and store the rank scores
    for wtd in ['','wtd_']:
        for e in epitopes:
            for p in nbrdist_percentiles:
                epitope_self_nbrdists[wtd][e][p].sort()

    for ii in range(total_tcrs):
        for epitope in epitopes:
            for wtd in ['','wtd_']:
                for nbrdist_percentile in nbrdist_percentiles:
                    self_nbrdists = epitope_self_nbrdists[wtd][epitope][nbrdist_percentile]
                    assert self_nbrdists[0] <= self_nbrdists[-1] #confirm sorted
                    nbrdist = all_info[ii]['{}_{}_{}nbrdist{}'.format(epitope,chains,wtd,nbrdist_percentile)]
                    rank = tcr_distances.get_rank( nbrdist, self_nbrdists )
                    all_info[ii]['{}_{}_{}nbrdist{}rank'.format(epitope,chains,wtd,nbrdist_percentile)] \
                        = '{:.3f}'.format(rank)

    # ## now rank versions of densities
    # for e in epitopes:
    #     for sdev in distance_sdevs:
    #         epitope_self_densities[e][sdev].sort()

    # for ii in range(total_tcrs):
    #     for epitope in epitopes:
    #         for distance_sdev in distance_sdevs:
    #             self_densities = epitope_self_densities[epitope][distance_sdev]
    #             assert self_densities[0] <= self_densities[-1] #confirm sorted
    #             dens = all_info[ii]['{}_{}_dens{}'.format(epitope,chains,distance_sdev)]
    #             rank = tcr_distances.get_rank( dens, self_densities )
    #             all_info[ii]['{}_{}_rdens{}'.format(epitope,chains,distance_sdev)] = '{:.3f}'.format(rank)



out = open( outfile,'w')
out.write('\t'.join(outfields)+'\n')

for outl in all_info:
    out.write( make_tsv_line( outl, outfields )+'\n' )

out.close()


