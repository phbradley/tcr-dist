def get_top_genes( blast_hits_string ):
    hits = dict( [ ( x.split(':')[0], int( x.split(':')[1] ) ) for x in blast_hits_string.split(';') ] )
    top_score = max( hits.values() )
    return set( [ x for x,y in hits.iteritems() if y >= top_score ] )

def get_top_reps( blast_hits_string, organism ):
    hits = dict( [ ( x.split(':')[0], int( x.split(':')[1] ) ) for x in blast_hits_string.split(';') ] )
    top_score = max( hits.values() )
    # vj = hits.keys()[0][3]
    # if vj == 'V':
    #     rep_map = cdr3s_human.all_loopseq_representative[ organism ]
    # else:
    #     assert vj == 'J'
    #     rep_map = cdr3s_human.all_jseq_representative[ organism ]
    return set( [ all_genes[organism][x].rep for x,y in hits.iteritems() if y >= top_score ] )


def reps_from_genes( genes, organism, mm1=False, trim_allele=False ):
    ## if genes is a set we can't index into it
    # vj = [ x[3] for x in genes ][0]

    # if vj == 'V':
    #     if mm1:
    #         rep_map = cdr3s_human.all_loopseq_representative_mm1[ organism ]
    #     else:
    #         rep_map = cdr3s_human.all_loopseq_representative[ organism ]
    # else:
    #     assert vj == 'J'
    #     rep_map = cdr3s_human.all_jseq_representative[ organism ]

    # reps = set( [ rep_map[x] for x in genes ] )
    reps = set( ( all_genes[organism][x].mm1_rep for x in genes ) ) if mm1 else \
           set( ( all_genes[organism][x].rep     for x in genes ) )
    if trim_allele:
        reps = set( ( x[:x.index('*')] for x in reps ) )
    return reps