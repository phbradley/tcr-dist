from basic import *
import util
import logo_tools

with Parser(locals()) as p:
    #p.str('args').unspecified_default().multiple().required()
    p.str('organism')
    p.str('infile')
    p.str('outfile')
    #p.int('min_mice').default(2)
    p.float('none_score_for_averaging').shorthand('none_score')     # --float_arg 9.6
    #p.flag('plot')       # --flag_arg  (no argument passed)
    p.flag('verbose').shorthand('v')       # --flag_arg  (no argument passed)
    #p.flag('allow_stop_codons')       # --flag_arg  (no argument passed)
    #p.flag('allow_X')       # --flag_arg  (no argument passed)
    p.flag('clobber').shorthand('c')       # --flag_arg  (no argument passed)
    #p.flag('add_masked_seqs')       # --flag_arg  (no argument passed)
    #p.flag('filter')       # --flag_arg  (no argument passed)
    p.int('min_quality_for_singletons').default(20)
    #p.flag('find_exact_matches')       # --flag_arg  (no argument passed)
    #p.range('range_arg')     # --range_arg 1:2
    #p.multiword('multi_arg') # --multi_arg hello world
    #p.file('file_arg')       # --file_arg README.txt
    #p.directory('dir_arg')   # --dir_arg /tmp/
    #p.str('floatlist').cast(lambda x: [float(val) for val in x.split(',')])
    p.multiword('average_clone_scores').cast(lambda x:x.split())

if exists(outfile):assert clobber

out = open(outfile,'w')

infields = []
outfields = []

all_tcrs = {}

total_lines, total_clones, skipcount = (0,0,0)



#TRAV10N*01:111;TRAV10D*01:111;TRAV10*02:111;TRAV10*01:111;TRAV10D*02:105;TRAV10*05:105;TRAV10*04:105;TRAV10*03:98


for line in open(infile,'rU'):
    if not infields:
        if line[0] == '#':
            infields = line[1:-1].split('\t')
        else:
            infields = line[:-1].split('\t')

        outfields = infields[:]

        assert 'id' in outfields
        outfields[ outfields.index('id') ] = 'clone_id'

        outfields.extend( ['members','clone_size'] )
                           #'va_genes','ja_genes','vb_genes','jb_genes',
                           #'num_va_reps','num_ja_reps','num_vb_reps','num_jb_reps' ] )

        out.write('\t'.join( outfields )+'\n' )
        continue
    assert infields
    total_lines += 1


    l = parse_tsv_line( line[:-1], infields )

    epitope = l['epitope']
    mouse = l['subject']
    va_gene = l['va_gene']
    ja_gene = l['ja_gene']
    cdr3a_nucseq = l['cdr3a_nucseq']

    vb_gene = l['vb_gene']
    jb_gene = l['jb_gene']
    cdr3b_nucseq = l['cdr3b_nucseq']

    l['cdr3a_min_qual'] = min( [int(x) for x in l['cdr3a_quals'].split('.') ] )
    l['cdr3b_min_qual'] = min( [int(x) for x in l['cdr3b_quals'].split('.') ] )
    l['cdr3_min_qual'] = min( l['cdr3a_min_qual'], l['cdr3b_min_qual'] )

    genesets = []
    for ab in 'ab':
        for vj in 'vj':
            genesets.append( set( l[vj+ab+'_genes'].split(';')))

    em = (epitope,mouse)

    if em not in all_tcrs:
        all_tcrs[em] = {}

    tcrseq = (va_gene,ja_gene,vb_gene,jb_gene,cdr3a_nucseq,cdr3b_nucseq)

    if tcrseq not in all_tcrs[em]:
        all_tcrs[em][tcrseq] = []

    all_tcrs[em][tcrseq].append( [l, genesets] )


## 2nas:  wskmyr
## 4nas:  n


def count_mismatches( a,b):
    assert len(a) == len(b)
    mismatches =0
    for x,y in zip(a,b):
        if not logo_tools.nucleotide_symbols_match(x,y):
            mismatches += 1
    return mismatches

def get_common_genes( tcrs ):
    all_genesets = []
    first_genesets = tcrs[0][1]
    for ii in range(4):
        genes = []
        for g in first_genesets[ii]:
            allfound=True
            for (l,genesets) in tcrs:
                if g not in genesets[ii]:
                    allfound=False
            if allfound:
                genes.append( g )
        all_genesets.append( set(genes) )
    return all_genesets


for em in all_tcrs:

    nbrs = {}
    for t1 in all_tcrs[em]:
        nbrs[t1] = [t1]

    quals={}
    for t1 in all_tcrs[em]:
        tcrs1 = all_tcrs[em][t1]

        qa1 = max( [x[0][ 'cdr3a_min_qual' ] for x in tcrs1 ] )
        qb1 = max( [x[0][ 'cdr3b_min_qual' ] for x in tcrs1 ] )

        quals[t1] = (qa1,qb1)

        cdr3a_new_nucseqs = list( set( [ x[0]['cdr3a_new_nucseq'] for x in tcrs1 ] ) )
        cdr3b_new_nucseqs = list( set( [ x[0]['cdr3b_new_nucseq'] for x in tcrs1 ] ) )

        cdr3_nucseq_prob1  = float( tcrs1[0][0][ 'a_nucseq_prob'  ] ) * float( tcrs1[0][0][ 'b_nucseq_prob' ] )
        cdr3_protseq_prob1 = float( tcrs1[0][0][ 'a_protseq_prob' ] ) * float( tcrs1[0][0][ 'b_protseq_prob' ] )

        cdr3_nucseq_prob1  = int( math.log10( cdr3_nucseq_prob1  ) ) if cdr3_nucseq_prob1 >0 else -99
        cdr3_protseq_prob1 = int( math.log10( cdr3_protseq_prob1 ) ) if cdr3_protseq_prob1>0 else -99

        assert len(cdr3a_new_nucseqs)==1
        assert len(cdr3b_new_nucseqs)==1

        cdr3a_new_nucseq1 = cdr3a_new_nucseqs[0]
        cdr3b_new_nucseq1 = cdr3b_new_nucseqs[0]

        ## what genes are present for all of this guys hits?
        genesets1 = get_common_genes( tcrs1 )

        for ii in range(4):
            assert t1[ii] in genesets1[ii]

        ## look for nearby guys
        for t2 in all_tcrs[em]:
            if t2<=t1:continue

            tcrs2 = all_tcrs[em][t2]
            genesets2 = get_common_genes( tcrs2 )

            qa2 = max( [x[0][ 'cdr3a_min_qual' ] for x in tcrs2 ] )
            qb2 = max( [x[0][ 'cdr3b_min_qual' ] for x in tcrs2 ] )

            cdr3a_new_nucseq2 = list( set( [ x[0]['cdr3a_new_nucseq'] for x in tcrs2 ] ) )[0]
            cdr3b_new_nucseq2 = list( set( [ x[0]['cdr3b_new_nucseq'] for x in tcrs2 ] ) )[0]

            cdr3_nucseq_prob2  = float( tcrs2[0][0][ 'a_nucseq_prob'  ] ) * float( tcrs2[0][0][ 'b_nucseq_prob' ] )
            cdr3_protseq_prob2 = float( tcrs2[0][0][ 'a_protseq_prob' ] ) * float( tcrs2[0][0][ 'b_protseq_prob' ] )

            cdr3_nucseq_prob2  = int( math.log10( cdr3_nucseq_prob2  ) ) if cdr3_nucseq_prob2 >0 else -99
            cdr3_protseq_prob2 = int( math.log10( cdr3_protseq_prob2 ) ) if cdr3_protseq_prob2>0 else -99

            mismatches=0
            samelen = True
            for ii in [4,5]:
                if len(t1[ii] ) != len(t2[ii] ):
                    samelen = False
                    break
                mismatches += count_mismatches( t1[ii], t2[ii] )

            s1,s2 = ( cdr3a_new_nucseq1 + ' ' + cdr3b_new_nucseq1,
                      cdr3a_new_nucseq2 + ' ' + cdr3b_new_nucseq2 )
            new_mismatches = count_mismatches(s1,s2) if len(s1)==len(s2) else 9

            common_genes = []
            for ii,genes1 in enumerate( genesets1 ):
                genes=[]
                for g in genes1:
                    if g in genesets2[ii]:
                        genes.append( g )
                common_genes.append( genes )
            clones_have_common_genes = ( sum( [len(x)>0 for x in common_genes ] )==4 )


            if samelen and mismatches<3 and clones_have_common_genes:
                if verbose:
                    print 'close by1: {:2d} {:2d} {} {} {:2d} {} {} {} {}'\
                        .format( qa1, qb1, mismatches, new_mismatches, len(all_tcrs[em][t1]),
                                 cdr3_nucseq_prob1, cdr3_protseq_prob1, s1, ' '.join(t1[:4]) )
                    print 'close by2: {:2d} {:2d} {} {} {:2d} {} {} {} {}'\
                        .format( qa2, qb2, mismatches, new_mismatches, len(all_tcrs[em][t2]),
                                 cdr3_nucseq_prob2, cdr3_protseq_prob2, s2, ' '.join(t2[:4]) )

                ## plan: only merge if one gene is perfect
                if t1[4] == t2[4] or t1[5] == t2[5]:
                    mmgene = 'A' if t1[4] != t2[4] else 'B'
                    q1 = qa1 if mmgene=='A' else qb1
                    q2 = qa2 if mmgene=='A' else qb2
                    new_nucseq1 = cdr3a_new_nucseq1 if mmgene=='A' else cdr3b_new_nucseq1
                    new_nucseq2 = cdr3a_new_nucseq2 if mmgene=='A' else cdr3b_new_nucseq2

                    minq_size = len(all_tcrs[em][t1]) if q1<q2 else len(all_tcrs[em][t2])
                    min_size = min( len(all_tcrs[em][t1]), len(all_tcrs[em][t2]) )

                    if verbose:
                        print 'merge1: {:2d} {} {} {:2d} {} {} {}'\
                            .format( q1, mismatches, new_mismatches, len(all_tcrs[em][t1]),
                                     cdr3_nucseq_prob1, new_nucseq1, ' '.join(t1[:4]) )
                        print 'merge2: {:2d} {} {} {:2d} {} {} {}'\
                            .format( q2, mismatches, new_mismatches, len(all_tcrs[em][t2]),
                                     cdr3_nucseq_prob2, new_nucseq2, ' '.join(t2[:4]) )


                    do_merge = ( mismatches==0 or
                                 ( mismatches==1 and min(q1,q2) < 20 and minq_size ==1 and
                                   ( new_mismatches<=1 or min(cdr3_nucseq_prob1,cdr3_nucseq_prob2)<-15 ) ) )

                    if do_merge:
                        print 'domerge1: {:2d} {} {} {:2d} {} {} {} {} {}'\
                            .format( q1, mismatches, new_mismatches, len(all_tcrs[em][t1]),
                                     cdr3_nucseq_prob1, new_nucseq1, em[0], em[1], ' '.join(t1[:4]) )
                        print 'domerge2: {:2d} {} {} {:2d} {} {} {} {} {}'\
                            .format( q2, mismatches, new_mismatches, len(all_tcrs[em][t2]),
                                     cdr3_nucseq_prob2, new_nucseq2, em[0], em[1], ' '.join(t2[:4]) )

                        ## which should we take
                        nbrs[ t1 ].append( t2 )
                        nbrs[ t2 ].append( t1 )


    seen = []

    for t1 in all_tcrs[em]:
        if t1 in seen: continue

        ## get the big set-- single linkage...
        all_nbrs = [t1]
        while True:
            old = all_nbrs[:]
            new = all_nbrs[:]
            for t in old:
                for nbr in nbrs[t]:
                    if nbr not in new:
                        new.append(nbr)
            all_nbrs = new[:]
            if len(new) == len(old):
                break

        ## now which one should be the representative?
        ## now we are doing quality filtering here
        sizel = []
        clone_size = 0
        members = []
        member_tcrs = []
        for t in all_nbrs:
            size = len(all_tcrs[em][t])
            clone_size += size
            sizel.append( ( size, min(quals[t]) , t ) )
            members.extend( [x[0]['id'] for x in all_tcrs[em][t] ] )
            member_tcrs.extend( all_tcrs[em][t] )
            assert t not in seen

        sizel.sort()
        sizel.reverse()
        assert len(members) == clone_size and len(member_tcrs) == clone_size


        if len(sizel)>1:
            print 'sizel:',[(x[0],x[1]) for x in sizel]

        aq,bq = quals[t1]

        if clone_size==1 and ( aq < min_quality_for_singletons or bq<min_quality_for_singletons ):
            print 'skipping singleton because min_quality lower than ' + str(min_quality_for_singletons)+':',aq,bq,t1[:4]
            skipcount+=1
            continue

        trep = sizel[0][-1]
        if t1 != trep:
            if verbose: print 'nonrep:',aq,bq,t1[:4],'rep:',trep[:4]
            continue

        ## ok, we are taking this guy as the rep, so mark all members as seen
        for t in all_nbrs:
            assert t not in seen
            seen.append( t )


        ## write out a new line
        l = all_tcrs[em][t1]

        outl = dict( l[0][0] )## copy

        clone_id = outl['id']+'.clone'
        #members = ';'.join( [ x[0]['id'] for x in l ] )

        outl['clone_id'] = clone_id
        del outl['id']

        outl['members'] = ';'.join(members)
        outl['clone_size'] = clone_size

        if average_clone_scores:
            for tsvtag in average_clone_scores:
                scores = []
                for t in all_nbrs:
                    for t_l, t_genesets in all_tcrs[em][t]:
                        score = float( t_l[tsvtag] )
                        if none_score_for_averaging == None or abs(score-none_score_for_averaging)>1e-3:
                            scores.append( score )
                if none_score_for_averaging==None:
                    assert len(scores) == clone_size
                if scores:
                    outl[tsvtag] = sum(scores) / len(scores)
                else:
                    assert none_score_for_averaging != None
                    outl[tsvtag] = none_score_for_averaging

        genesets = get_common_genes( member_tcrs )

        for ii in range(4):
            if not genesets[ii]: ## whoah-- no overlap??
                counts = {}
                for (l,gsets) in member_tcrs:
                    for g in gsets[ii]:
                        counts[g] = counts.get(g,0)+1
                mx = max(counts.values())
                genesets[ii] = set( [ x for x,y in counts.iteritems() if y==mx ] )
                print 'empty common genes:',ii,'clone_size:',clone_size,'mx-genecount',mx,'newgeneset:',genesets[ii],em

        for genes,segtype in zip( genesets, segtypes_lowercase ):
            assert genes
            tag = segtype+'_genes'
            assert tag in outl # should already be there, now over-writing
            outl[tag] = ';'.join(sorted(genes))

            ## update reps
            reps = sorted( set( ( util.get_rep(x,organism) for x in genes ) ) )
            tag = segtype+'_reps'
            assert tag in outl # should already be there, now over-writing
            outl[tag] = ';'.join(reps)

            ## update countreps
            countreps = sorted( set( ( util.get_mm1_rep_gene_for_counting(x,organism) for x in genes ) ) )
            tag = segtype+'_countreps'
            assert tag in outl # should already be there, now over-writing
            outl[tag] = ';'.join(countreps)

        total_clones += 1

        out.write( make_tsv_line( outl, outfields, '-' )+'\n' )
        out.flush()

out.close()
print 'skipcount:',skipcount,'total_lines:',total_lines,'total_clones:',total_clones
