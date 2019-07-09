from basic import *
import parse_tsv
from collections import Counter


# yes, this is very silly, should just add pandas dependency
def parse_csv_file( csvfile ):
    header = None
    all_info = []
    for line in open(csvfile,'rU'):
        l = parse_tsv.safely_split_csv_line( line[:-1] )
        if header is None:
            header = l
        else:
            all_info.append( dict( zip( header,l ) ) )
    return header, all_info

def read_tcr_data( organism, contig_annotations_csvfile, consensus_annotations_csvfile,
                   include_gammadelta = False, allow_unknown_genes = False ):
    """ Parse tcr data, only taking 'productive' tcrs

    Returns:

    clonotype2tcrs, clonotype2barcodes

    """
    from all_genes import all_genes

    expected_gene_names = all_genes[organism].keys()

    #from cdr3s_human import all_align_fasta

    gene_suffix = '*01' # may not be used


    # read the contig annotations-- map from clonotypes to barcodes
    # barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id
    # AAAGATGGTCTTCTCG-1,True,AAAGATGGTCTTCTCG-1_contig_1,True,695,TRB,TRBV5-1*01,TRBD2*02,TRBJ2-3*01,TRBC2*01,True,True,CASSPLAGYAADTQYF,TGCGCCAGCAGCCCCCTAGCGGGATACGCAGCAGATACGCAGTATTTT,9427,9,clonotype14,clonotype14_consensus_1
    assert exists( contig_annotations_csvfile )

    _, lines = parse_csv_file(contig_annotations_csvfile)
    clonotype2barcodes = {}
    for l in lines:
        bc = l['barcode']
        clonotype = l['raw_clonotype_id']
        if clonotype =='None':
            if l['productive'] not in [ 'None','False' ]:
                assert l['productive'] == 'True'
                #print 'clonotype==None: unproductive?',l['productive']
            continue
        if clonotype not in clonotype2barcodes:
            clonotype2barcodes[clonotype] = []
        if bc in clonotype2barcodes[clonotype]:
            pass
            #print 'repeat barcode'
        else:
            clonotype2barcodes[clonotype].append( bc )


    ## now read details on the individual chains for each clonotype
    # ==> tcr/human/JCC176_TX2_TCR_consensus_annotations.csv <==
    # clonotype_id,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis
    # clonotype100,clonotype100_consensus_1,550,TRB,TRBV24-1*01,TRBD1*01,TRBJ2-7*01,TRBC2*01,True,True,CATSDPGQGGYEQYF,TGTGCCACCAGTGACCCCGGACAGGGAGGATACGAGCAGTACTTC,8957,9

    assert exists(consensus_annotations_csvfile)
    _, lines = parse_csv_file( consensus_annotations_csvfile )


    ## first get clonotypes with one alpha and one beta
    clonotype2tcrs = {}

    for l in lines:
        if l['productive'] == 'True':
            id = l['clonotype_id']
            if id not in clonotype2tcrs:
                # dictionaries mapping from tcr to umi-count
                clonotype2tcrs[id] = { 'A':Counter(), 'B':Counter() } #, 'G':[], 'D': [] }
                assert id in clonotype2barcodes

            ch = l['chain']
            if not ch.startswith('TR'):
                print 'skipline:', consensus_annotations_csvfile, ch, l['v_gene'], l['j_gene']
                continue
            ab = ch[2]
            if ab not in 'AB':
                print 'skipline:', consensus_annotations_csvfile, ch, l['v_gene'], l['j_gene']
                continue
            vg = l['v_gene']
            if '*' not in vg:
                vg += gene_suffix
            if 'DV' in vg and vg not in expected_gene_names:
                #print 'DV?',vg
                vg = vg[:vg.index('DV')]+'/'+vg[vg.index('DV'):]
            jg = l['j_gene']
            if '*' not in jg:
                jg += gene_suffix
            # if vg in tcr_gene_remap[organism]:
            #     vg = tcr_gene_remap[organism][vg]
            # if jg in tcr_gene_remap[organism]:
            #     jg = tcr_gene_remap[organism][jg]

            if vg not in expected_gene_names:
                print 'unrecognized V gene:', organism, vg
                if not allow_unknown_genes:
                    continue
            if vg not in expected_gene_names or jg not in expected_gene_names:
                print 'unrecognized J gene:', organism, jg
                if not allow_unknown_genes:
                    continue
            #assert vg in all_align_fasta[organism]
            #assert jg in all_align_fasta[organism]

            tcr_chain = ( vg, jg, l['cdr3'], l['cdr3_nt'].lower() )

            if tcr_chain not in clonotype2tcrs[id][ab]:
                clonotype2tcrs[id][ab][ tcr_chain ] = int( l['umis'] )
            else:
                print 'repeat?',id,ab,tcr_chain
        else:
            if l['productive'] not in [ 'None','False' ]:
                print 'unproductive?',l['productive']

    return clonotype2tcrs, clonotype2barcodes

def make_clones_file( organism, outfile, clonotype2tcrs, clonotype2barcodes ):
    ''' Make a clones file with information parsed from the 10X csv files

    organism is one of ['mouse','human']

    outfile is the name of the clones file to be created

    '''

    tmpfile = outfile+'.tmp' # a temporary intermediate file

    out = open(tmpfile,'w')
    outfields = 'clone_id subject clone_size va_gene ja_gene vb_gene jb_gene cdr3a cdr3a_nucseq cdr3b cdr3b_nucseq'\
        .split()
    extra_fields = 'alpha_umi beta_umi num_alphas num_betas'.split()
    outfields += extra_fields

    out = open(tmpfile,'w')
    out.write('\t'.join( outfields )+'\n' )

    for clonotype, tcrs in clonotype2tcrs.iteritems():
        if len(tcrs['A']) >= 1 and len(tcrs['B']) >= 1:
            atcrs = tcrs['A'].most_common()
            btcrs = tcrs['B'].most_common()
            if len(atcrs)>1:
                print 'multiple alphas, picking top umi:',' '.join( str(x) for _,x in atcrs )
            if len(btcrs)>1:
                print 'multiple  betas, picking top umi:',' '.join( str(x) for _,x in btcrs )
            atcr, atcr_umi = atcrs[0]
            btcr, btcr_umi = btcrs[0]
            outl = {}
            outl['clone_id'] = clonotype
            outl['subject'] = 'UNK_S'
            outl['clone_size'] = len(clonotype2barcodes[clonotype])
            outl['va_gene']      = atcr[0]
            outl['ja_gene']      = atcr[1]
            outl['cdr3a']        = atcr[2]
            outl['cdr3a_nucseq'] = atcr[3]
            outl['alpha_umi']    = str(atcr_umi)
            outl['num_alphas']   = str(len(atcrs))
            outl['vb_gene']      = btcr[0]
            outl['jb_gene']      = btcr[1]
            outl['cdr3b']        = btcr[2]
            outl['cdr3b_nucseq'] = btcr[3]
            outl['beta_umi']     = str(btcr_umi)
            outl['num_betas']    = str(len(btcrs))
            out.write( make_tsv_line(outl,outfields)+'\n' )
    out.close()


    cmd = 'python {}/file_converter.py --input_format clones --output_format clones --input_file {} --output_file {}  --organism {} --clobber --epitope UNK_E --check_genes --extra_fields {} '\
        .format( paths.path_to_scripts, tmpfile, outfile, organism, ' '.join(extra_fields) )
    print cmd
    system(cmd)


#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################

if __name__ == '__main__':

    with Parser(locals()) as p:
        p.str('filtered_contig_annotations_csvfile').shorthand('f').required()
        p.str('consensus_annotations_csvfile').shorthand('c').required()
        p.str('clones_file').shorthand('o').required()
        p.str('organism').required()
        p.flag('clobber')

    if exists(clones_file) and not clobber:
        print 'ERROR -- clones_file {} already exists; use --clobber to overwrite'
        exit(1)

    assert organism in ['human','mouse']

    clonotype2tcrs, clonotype2barcodes = read_tcr_data( organism, filtered_contig_annotations_csvfile,
                                                        consensus_annotations_csvfile )

    make_clones_file( organism, clones_file, clonotype2tcrs, clonotype2barcodes )


