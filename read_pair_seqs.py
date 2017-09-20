from basic import *
import read_sanger_data
import sys
import logo_tools
import util

with Parser(locals()) as p:
    p.str('infile').required()
    p.str('outfile').required()
    p.str('organism').required()
    #p.str('Achain').default('A').described_as("Use 'G' for gamma delta") ### coming soon!
    #p.str('Bchain').default('B').described_as("Use 'D' for gamma delta") ### coming soon!
    p.flag('verbose')       # --flag_arg  (no argument passed)
    p.flag('clobber').shorthand('c')       # --flag_arg  (no argument passed)
    p.flag('dry_run')       # --flag_arg  (no argument passed)
    p.flag('save_results')       # --flag_arg  (no argument passed)
    p.flag('make_fake_alpha').described_as("Create a fictitious alpha chain sequence")
    p.flag('make_fake_beta').described_as("Create a fictitious beta chain sequence")
    p.flag('make_fake_quals').described_as("Create a fictitious quality string for each nucleotide sequence")
    p.flag('make_fake_ids').described_as("Create an id for each line based on index in the file")
    p.set_help_prefix("""
    This script reads a paired sequence file and tries to assign V/J genes and parse the CDR3 amino acid
    sequences of the TCRs in the file.

    The required fields in the input .tsv (tab-separated values) file are

    id -- a unique identifier for each line
    epitope -- the epitope for which the tcr is specific
    subject -- an identifier for the subject from whom the tcr was sampled
    a_nucseq -- the alpha chain sequence read
    b_nucseq -- the beta chain sequence read
    a_quals -- a '.'-separated list of the quality scores for a_nucseq; e.g. 19.25.36.40.38.20 (etc)
    b_quals -- like a_quals but for the b_nucseq sequence read

    If you don't have read quality info or it's a nuisance to add it to the file, you can use the option

        --make_fake_quals

    which will assign fictitious (high) quality scores to each read. Note that this will prevent any
    filtering of TCRs with low-quality CDR3 regions, which would ordinarily happen in find_clones.py

    If you don't have ids, you can use the option

        --make_fake_ids

    which will assign an id of the form tcr<N> to each tcr based on index in the file.

    """)

default_outfields = ['id','epitope','subject',
                     'va_gene','va_rep','va_mismatches','va_alignlen','va_evalue','va_bitscore_gap',
                     'ja_gene','ja_rep','ja_mismatches','ja_alignlen','ja_evalue','ja_bitscore_gap',
                     'vb_gene','vb_rep','vb_mismatches','vb_alignlen','vb_evalue','vb_bitscore_gap',
                     'jb_gene','jb_rep','jb_mismatches','jb_alignlen','jb_evalue','jb_bitscore_gap',
                     'cdr3a','cdr3a_nucseq','cdr3a_quals',
                     'cdr3b','cdr3b_nucseq','cdr3b_quals',
                     'va_genes','va_reps','va_countreps',
                     'ja_genes','ja_reps','ja_countreps',
                     'vb_genes','vb_reps','vb_countreps',
                     'jb_genes','jb_reps','jb_countreps',
                     'va_blast_hits','ja_blast_hits',
                     'vb_blast_hits','jb_blast_hits',
                     'status' ]


assert outfile.endswith('tsv') ## stupid

if exists(outfile): assert clobber


fake_nucseqs = {
    'human': { 'A':'CGTACATAAGNTATTTGTCTGTGATATACACATCAGAATCCTTACTTTGTGACACATTTGTTTGAGAATCAAAATCGGTGAATAGGCAGACAGACTTGTCACTGGATTTAGAGTCTCTCAGCTGGTACACGGCAGGGTCAGGCTTCTGGATATTGGGCTTCACCACCAGCTGAGTTCCATCTCCAAACATGAGTCTGGCATTGTTATTCCGGACAGCACAGAAGTACAAAGCGGAGTCGCTCACAAGGGCAGATGGTTTCTTCAGGTGGAAGGAGGTTTGACTCTAGTTAAATTCAGCTTCAAAGCCATAGCTGCCTTTAACCAGGTTATCCCCTGTGATGTATTTCAGAAGGAACTGGAGGCCTCGGTTGGGGTATTGAACATACCAAAAAAGATAAGGGTTTCCAGAGACTGAATATGTGCATAAAA',
               'B':'CCGGGGGTGAAACCTTGATTAGGTCCTCTACAACTGTGAGTCTGGTGCCTTGTCCAAAGAAAGCTTCGTTTACAATGCTGCTGGCACAGAAGTACACAGCTGAGTCCCTGGGTTCTGAGGGCTGGATCTTCAGAGTGGAAACAAN',
               'G':'',
               'D':'',
               },
    'mouse': {'A':'GGGAAAGCGTAGGGGTGCTGTCCTGAGACCGAGGATCTTTTAACTGGTACACAGCAGGTTCTGGGTTCTGGATGTCTGGCTTTATAATTAGCTTGGTCCCAGAGCCCCAGATCAACTGATAGTTGCTATCCATGTCAGTAGCACAGAAGTACACAGCGGTGTCTTCACACCGAGTAGCAGTGATGGACAAGGAGCTGCTCTGGCTGGAGGTGTCAAGGGTGGCTCTAAAAN',
              'B':'CACCGGGGGGGTCGGGGAAGAAGCCCCTNGGCCAGCACACGAGGGTAGCCTTTTGTTTGTTTGCAATCTCTGCTTTTGATGGCTCAAACAAGGAGACCTTGGGTGGAGTCACATTTCTCAGATCCTCCAGGACAGACAGCTTGGTTCCATGACCGAAAAATAATCTTTCGTTGGCCCCCATACTGCTGGCACAGAGAAAAACGGCCATCTCGTTCTTCTGGGCAGATGTCACAGTGAGAGAAAAAGATGACTTCTTCTCTCGAGACGCATCATAGCCTTCAGATAGATCGCCTTTTTGAAGATCGTTTTCAGTTATTGAATAGTAGATCAGTCTCAATCCTTTCCCTGAATCCTGTCGGTACCAGN',
              'G':'',
              'D':''
              }
    }


def find_largest_substring_contained_in( seq, bigseq ):
    best = [len(seq),len(seq)]
    for ntrim in range(len(seq)):
        for ctrim in range(len(seq)):
            if ntrim+ctrim>=len(seq):
                break
            if seq[ntrim:len(seq)-ctrim] in bigseq:
                if ntrim+ctrim < best[0]+best[1]:
                    best = [ntrim,ctrim]
                break # no sense increasing ctrim
    return best[0],best[1],seq[best[0]:len(seq)-best[1]]


def get_qualstring( cdr3seqtag, nucseq_in, quals_in ):
    # the tricky part here (new: 09/19/17) is that we might not be getting the entire CDR3 sequence from the read
    # We're now allowing a little to come from the gene
    # also nucseq and quals may be in the opposite orientation
    nucseq = nucseq_in[:].upper()
    quals = quals_in[:]

    if '-' not in cdr3seqtag: return ''
    if len( cdr3seqtag.split('-') )!= 2: return ''
    cdr3nucseq = cdr3seqtag.split('-')[1].upper()
    if not cdr3nucseq: return ''
    revnucseq = logo_tools.reverse_complement(nucseq)
    if cdr3nucseq not in nucseq:
        if cdr3nucseq in revnucseq:
            nucseq = logo_tools.reverse_complement(nucseq)
            quals.reverse()
        else:
            #doesn't match either
            ntrimf,ctrimf,substring_f = find_largest_substring_contained_in( cdr3nucseq, nucseq    ) #fwd
            ntrimr,ctrimr,substring_r = find_largest_substring_contained_in( cdr3nucseq, revnucseq ) #rev
            if ntrimr+ctrimr < ntrimf+ctrimf:
                # better in reverse direction
                nucseq = revnucseq
                quals.reverse()
                ntrimf,ctrimf,substring_f = ntrimr,ctrimr,substring_r
            assert substring_f in nucseq
            pos = nucseq.find(substring_f)
            good_quality = 60
            cdr3quals = [good_quality]*ntrimf + quals[pos:pos+len(substring_f)] + [good_quality]*ctrimf
            return '.'.join( [ '%d'%x for x in cdr3quals ] )


    if nucseq.count( cdr3nucseq ) == 0:
        # print cdr3nucseq
        # print nucseq
        # print nucseq_in
        return ''
    elif nucseq.count( cdr3nucseq ) == 1:
        pos = nucseq.index( cdr3nucseq)
    else:
        ## multiple occurrences of cdr3nucseq in nucseq, choose the one with the best quality
        max_min_qual = 0
        bestpos = 0
        offset = 0
        while nucseq.count( cdr3nucseq, offset ):
            pos = nucseq.find( cdr3nucseq, offset )
            min_qual = min( quals[pos:pos+len(cdr3nucseq)] )
            if verbose:
                print 'multiple matches:',pos,min_qual,max_min_qual
            if min_qual>max_min_qual:
                max_min_qual = min_qual
                bestpos = pos
            offset = pos+1
        pos = bestpos

    return '.'.join( [ '%d'%x for x in quals[pos:pos+len(cdr3nucseq)] ] )

num_lines = 0
successes = 0
infields = []

saved_results = {}


for line in open( infile,'rU'):
    if line[0] == '#' or not infields:
        assert not infields
        if line[0] == '#':
            infields = line[1:-1].split('\t')
        else:
            infields = line[:-1].split('\t')
        skip_infields = ['a_nucseq','b_nucseq','a_quals','b_quals'] ## dont put these guys into the outfile
        extra_outfields = []
        for field in infields:
            if field not in skip_infields and field not in default_outfields:
                extra_outfields.append( field )

        out = open(outfile,'w')
        out.write('#{}\n'.format('\t'.join(default_outfields+extra_outfields)))

        continue
    assert infields

    if num_lines and num_lines%100==0:
        print 'num_lines: {} successes: {}'.format(num_lines,successes)
        sys.stdout.flush()

    num_lines += 1

    l = parse_tsv_line( line[:-1], infields )
    #l = line[:-1].split('\t')
    #assert len(l) == len(infields)

    if make_fake_ids:
        id = 'tcr{:d}'.format(num_lines)
    else:
        id = l[ 'id' ]
    epitope = l[ 'epitope' ]
    mouse = l[ 'subject' ]

    if not epitope or not mouse or epitope=='-' or mouse=='-':
        print 'skipping line with bad epitope/subject',epitope,mouse
        continue

    ## nucseqs and quals
    if make_fake_alpha:
        aseq = fake_nucseqs[organism]['A']
        aquals = [60]*len(aseq)
    else:
        aseq = l[ 'a_nucseq' ]
        if not aseq or aseq == '-':
            print 'skipping line with bad a_nucseq:',aseq,line[:-1]
            continue
        if make_fake_quals:
            aquals = [60]*len(aseq)
        else:
            aquals = map( int, l[ 'a_quals' ].split('.') )

    if make_fake_beta:
        bseq = fake_nucseqs[organism]['B']
        bquals = [60]*len(bseq)
    else:
        bseq = l[ 'b_nucseq' ]
        if not bseq or bseq == '-':
            print 'skipping line with bad b_nucseq:',bseq,line[:-1]
            continue
        if make_fake_quals:
            bquals = [60]*len(bseq)
        else:
            bquals = map( int, l[ 'b_quals' ].split('.') )

    assert len(aseq) == len(aquals)
    assert len(bseq) == len(bquals)

    if dry_run: continue

    program = 'blastn'
    #program = 'blastx' ## hacking!!!!!!!!

    if save_results and (aseq,bseq) in saved_results:
        ab_genes,evalues,status,all_hits = saved_results[ (aseq,bseq)]
        #print 'recover saved'
    else:
        ab_genes,evalues,status,all_hits = read_sanger_data.parse_paired_dna_sequences( program, organism, aseq, bseq,
                                                                                        info=id, verbose=verbose,
                                                                                        extended_cdr3 = True,
                                                                                        return_all_good_hits = True,
                                                                                        nocleanup = verbose )
        if save_results:
            saved_results[ (aseq,bseq) ] = ab_genes,evalues,status,all_hits

    ahits,bhits = all_hits

    if not status:
        status = ['OK']
        successes += 1
    #evalues_info = ';'.join(['%s:%d:%.3g'%(x,evalues[x][1],evalues[x][0]) for x in evtags] )

    assert len(ab_genes)==2
    va_gene, va_rep, va_mm, ja_gene, ja_rep, ja_mm, cdr3a_plus  = ab_genes['A']
    vb_gene, vb_rep, vb_mm, jb_gene, jb_rep, jb_mm, cdr3b_plus  = ab_genes['B']

    cdr3a,cdr3a_nucseq,cdr3a_quals=['-','-','-']

    if '-' in cdr3a_plus:
        cdr3a,cdr3a_nucseq = cdr3a_plus.split('-')
        cdr3a_quals = get_qualstring( cdr3a_plus, aseq, aquals )

    cdr3b,cdr3b_nucseq,cdr3b_quals=['-','-','-']

    if '-' in cdr3b_plus:
        cdr3b,cdr3b_nucseq = cdr3b_plus.split('-')
        if verbose:
            r = cdr3b_nucseq.upper()
            s1 = bseq.upper()
            s2 = logo_tools.reverse_complement( s1 )
            print 'cdr3b_nucseq in b_nucseq:',r in s1, r in s2
        cdr3b_quals = get_qualstring( cdr3b_plus, bseq, bquals )

    ## the problem is that blast may return multiple hits with equal bitscore, so store also the full list of top hits
    if ahits and len(ahits)==2 and ahits[0] and ahits[1]:
        va_blast_hits = ';'.join( '{}:{}'.format(x[0],x[1]) for x in ahits[0] )
        ja_blast_hits = ';'.join( '{}:{}'.format(x[0],x[1]) for x in ahits[1] )
        va_genes = util.get_top_genes( va_blast_hits ) ## a set
        ja_genes = util.get_top_genes( ja_blast_hits )
        va_genestring = ';'.join( sorted( va_genes ) )
        ja_genestring = ';'.join( sorted( ja_genes ) )
        va_repstring  = ';'.join( sorted( util.get_top_reps( va_blast_hits, organism ) ) )
        ja_repstring  = ';'.join( sorted( util.get_top_reps( ja_blast_hits, organism ) ) )
        va_countrepstring = ';'.join( sorted( set( (util.get_mm1_rep_gene_for_counting(x,organism) for x in va_genes ))))
        ja_countrepstring = ';'.join( sorted( set( (util.get_mm1_rep_gene_for_counting(x,organism) for x in ja_genes ))))
    else:
        va_blast_hits = '-'
        ja_blast_hits = '-'
        va_genestring = '-'
        ja_genestring = '-'
        va_repstring  = '-'
        ja_repstring  = '-'
        va_countrepstring = '-'
        ja_countrepstring = '-'

    if bhits and len(bhits)==2 and bhits[0] and bhits[1]:
        vb_blast_hits = ';'.join( '{}:{}'.format(x[0],x[1]) for x in bhits[0] )
        jb_blast_hits = ';'.join( '{}:{}'.format(x[0],x[1]) for x in bhits[1] )
        vb_genes = util.get_top_genes( vb_blast_hits ) ## a set
        jb_genes = util.get_top_genes( jb_blast_hits )
        vb_genestring = ';'.join( sorted( vb_genes ) )
        jb_genestring = ';'.join( sorted( jb_genes ) )
        vb_repstring  = ';'.join( sorted( util.get_top_reps( vb_blast_hits, organism ) ) )
        jb_repstring  = ';'.join( sorted( util.get_top_reps( jb_blast_hits, organism ) ) )
        vb_countrepstring = ';'.join( sorted( set( (util.get_mm1_rep_gene_for_counting(x,organism) for x in vb_genes ))))
        jb_countrepstring = ';'.join( sorted( set( (util.get_mm1_rep_gene_for_counting(x,organism) for x in jb_genes ))))
    else:
        vb_blast_hits = '-'
        jb_blast_hits = '-'
        vb_genestring = '-'
        jb_genestring = '-'
        vb_repstring  = '-'
        jb_repstring  = '-'
        vb_countrepstring = '-'
        jb_countrepstring = '-'



    ## the order here has to match the order of default_outfields
    vals = [ id, epitope, mouse,
             va_gene, va_rep, str(va_mm[0]), str( va_mm[0]+va_mm[1] ),'%.3g'%(evalues['VA'][0]), str(evalues['VA'][1]),
             ja_gene, ja_rep, str(ja_mm[0]), str( ja_mm[0]+ja_mm[1] ),'%.3g'%(evalues['JA'][0]), str(evalues['JA'][1]),
             vb_gene, vb_rep, str(vb_mm[0]), str( vb_mm[0]+vb_mm[1] ),'%.3g'%(evalues['VB'][0]), str(evalues['VB'][1]),
             jb_gene, jb_rep, str(jb_mm[0]), str( jb_mm[0]+jb_mm[1] ),'%.3g'%(evalues['JB'][0]), str(evalues['JB'][1]),
             cdr3a,cdr3a_nucseq,cdr3a_quals,
             cdr3b,cdr3b_nucseq,cdr3b_quals,
             va_genestring, va_repstring, va_countrepstring,
             ja_genestring, ja_repstring, ja_countrepstring,
             vb_genestring, vb_repstring, vb_countrepstring,
             jb_genestring, jb_repstring, jb_countrepstring,
             va_blast_hits,ja_blast_hits,vb_blast_hits,jb_blast_hits,
             ';'.join(status) ]


    assert len(vals) == len(default_outfields)

    for tag in extra_outfields:
        vals.append( l[ tag ] )

    vals = [ (x if x else '-') for x in vals ] # remove empty strings, makes awk happier

    out.write('\t'.join( vals )+'\n' )



    # print 'blastn_pairs %s %s %s %d %d %s %s %d %d %s   %s %s %d %d %s %s %d %d %s  %s %s %d %s %s %s %s %s %s %s'\
    #     %( epi,
    #        va_gene, va_rep, va_mm[0], va_mm[0]+va_mm[1], ja_gene, ja_rep, ja_mm[0], ja_mm[0]+ja_mm[1], cdr3a,
    #        vb_gene, vb_rep, vb_mm[0], vb_mm[0]+vb_mm[1], jb_gene, jb_rep, jb_mm[0], jb_mm[0]+jb_mm[1], cdr3b,
    #        mouse, evalues_info,
    #        counter, ';'.join( status ),
    #        aseq, bseq,
    #        '.'.join( [`x` for x in aquals ] ),
    #        '.'.join( [`x` for x in bquals ] ),
    #        afile, bfile )
out.close()

print 'num_lines: {} successes: {}'.format(num_lines,successes)
