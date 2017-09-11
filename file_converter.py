##
##
##



from basic import *
import util
import copy
import sys
#import matplotlib
#if make_png: matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#import numpy as np

with Parser(locals()) as p:
    #p.str('args').unspecified_default().multiple().required()
    p.str('organism').required()
    p.str('input_file').required()
    p.str('output_file').required()
    p.str('input_format').required().described_as('Should be one of: parsed_seqs -OR- clones -OR cdrblast')
    p.str('output_format').required().described_as('Should be one of: parsed_seqs -OR- clones')
    p.str('single_chain').described_as( 'Only generate info for a single chain; should be one of: alpha -OR- beta')
    p.str('epitope').described_as( 'Use this value for the epitope field')
    p.str('subject').described_as( 'Use this value for the subject field')
    p.flag('auto_ids').described_as('Auto-generate numbered TCR ids')
    p.str('id_base').described_as('If using --auto_ids, you can specify a base name for the IDs.')
    p.flag('clobber').shorthand('c')
    p.set_help_prefix("""

    This script is for converting old TCR-dist formats to the new format and for converted MIGEC output to an input format suitable for TCR-dist.

    Usage examples:

python ../file_converter.py  --input_format cdrblast --output_format clones --input_file MIGEC_d2_SH2_3.filtered.cdrblast.txt --output_file tmp.d2.tsv --organism human -c --auto_ids --epitope ring1 --subject frodo --single_chain beta

python ../file_converter.py --input_format parsed_seqs --output_format parsed_seqs --input_file JCC42_EpiMouse_AsEpi.tsv --output_file tmp.tsv --organism mouse -c

    """)

warnings = set()

prob_warning = '[WARNING] Setting a value of 1.0 for missing TCR probabilities'
clone_size_warning = '[WARNING] Setting a default clone_size of 1 for lines without clone_size info'
cdrblast_clone_size_warning = '[WARNING] Using the cdrblast "Count" field as the TCR clone_size value'
rep_warning = '[WARNING] Blast results do not seem to be specific to allele'

if exists( output_file ) and not clobber:
    print output_file,'already exists, use --clobber or remove it'
    exit(1)

## these are the fields we are going to include in the output file
##
required_fields_parsed_seqs_file = """
id epitope subject
va_gene va_genes va_rep va_reps va_countreps
ja_gene ja_genes ja_rep ja_reps ja_countreps
vb_gene vb_genes vb_rep vb_reps vb_countreps
jb_gene jb_genes jb_rep jb_reps jb_countreps
cdr3a cdr3a_nucseq cdr3a_quals
cdr3b cdr3b_nucseq cdr3b_quals
""".split()


required_fields_clones_file = """
clone_id epitope subject
clone_size
va_gene va_genes va_rep va_reps va_countreps
ja_gene ja_genes ja_rep ja_reps ja_countreps
vb_gene vb_genes vb_rep vb_reps vb_countreps
jb_gene jb_genes jb_rep jb_reps jb_countreps
cdr3a cdr3a_nucseq
cdr3b cdr3b_nucseq
a_protseq_prob
cdr3a_protseq_prob
va_rep_prob
ja_rep_prob
a_nucseq_prob
b_protseq_prob
cdr3b_protseq_prob
vb_rep_prob
jb_rep_prob
b_nucseq_prob
""".split()
## not doing 'members' field



## these are the fields we are going to include in the output file
## now stored in a dictionary by file-type
all_required_fields = {
    'parsed_seqs': required_fields_parsed_seqs_file,
    'clones': required_fields_clones_file,
    }

assert output_format in all_required_fields

outfields = all_required_fields[ output_format ]

if single_chain:
    assert single_chain in ['alpha','beta']
#    for pos in reversed( range( len(outfields ) ) ):
#        field = outfields[pos]
#        if   single_chain=='alpha' and ( field[:3] in ['vb_','jb_'] or field.startswith('cdr3b') ):
#            del outfields[pos]
#        elif single_chain== 'beta' and ( field[:3] in ['va_','ja_'] or field.startswith('cdr3a') ):
#            del outfields[pos]


listsep = ';' # for lists within fields


def remove_wonky_characters_from_the_end_of_line( inline ):
    outline = inline[:]
    while not ( outline[-1] == '\t' or outline[-1].split() ):
        outline = outline[:-1]
    return outline



# def exit_from_reconstruct_field_from_data( field, l ):
#     exit(1)

def reconstruct_field_from_data( field, l, organism ):
    try:
        if field in l:
            return l[ field ]

        if field == 'subject':
            if 'mouse' in l:
                return l['mouse']
        elif field.endswith('_prob'):
            if prob_warning not in warnings:
                print prob_warning
                warnings.add( prob_warning )
            return 1.0
        elif field == 'clone_size':
            if clone_size_warning not in warnings:
                print clone_size_warning
                warnings.add( clone_size_warning )
            return 1
        elif field[:3] in ['va_','ja_','vb_','jb_']: ## these are all pretty similar
            prefix = field[:3]
            tag = field[3:]
            #print 'prefix:',prefix,tag,l.keys()
            if tag == 'gene':
                # only one place to get this
                genes_field = prefix+'genes'
                if genes_field not in l:
                    return None
                genes = l[ genes_field ].split( listsep )
                #print 'genes:',genes
                return sorted( genes )[0]
            elif tag == 'genes':
                if prefix+'blast_hits' in l:
                    hits = l[ prefix+'blast_hits' ]
                    if hits == '-': #failed to find any hits
                        return None
                    else:
                        return listsep.join( sorted( util.get_top_genes( hits ) ) )
                elif prefix+'gene' in l: ## just take the one gene we know
                    return l[ prefix+'gene' ]
            elif tag == 'rep':
                if prefix+'gene' in l:
                    if "*" in l[prefix + 'gene']:
                        return util.get_rep( l[ prefix+'gene' ], organism )
                    else:
                        if rep_warning not in warnings:
                            print rep_warning
                            warnings.add( rep_warning )
                        return l[prefix + 'gene']
            elif tag == 'reps':
                ## we should already have hit 'genes' in the list of fields we are trying to fill !!!
                if "*" in l[prefix + 'gene']:
                    return listsep.join( sorted( ( util.get_rep(x,organism) for x in l[prefix+'genes'].split(listsep) ) ) )
                else:
                    if rep_warning not in warnings:
                        print rep_warning
                        warnings.add( rep_warning )
                    return l[prefix + 'gene']
            elif tag == 'countreps':
                if "*" in l[prefix + 'gene']:
                    return listsep.join( sorted( util.countreps_from_genes( l[prefix+'genes'].split(listsep), organism ) ) )
                else:
                    if rep_warning not in warnings:
                        print rep_warning
                        warnings.add( rep_warning )
                    return l[prefix + 'gene']
        elif field.endswith('_quals'):
            seqfield = field[:5]+'_nucseq'
            if seqfield not in l:
                return None
            return '.'.join( ['60']*len(l[seqfield] ) )

    except Exception as inst: ## this is not the best way to handle it...
        print 'Hit an exception trying to get field {} from line'.format(field),inst

    print 'Failed to reconstruct {} from the input fields: {}'\
        .format( field, ' '.join( sorted( l.keys() ) ) )
    return None



def map_field_names( l, input_format, output_format ):
    outl = copy.deepcopy( l )

    if input_format == 'cdrblast':
        ## the problem is that this line could be either alpha or beta; they can occur in the same file
        if 'V segments' in l.keys():
            v_genestring = l['V segments']
            j_genestring = l['J segments']
            cdr3aa = l['CDR3 amino acid sequence']
            cdr3nt = l['CDR3 nucleotide sequence']
            clonecount = l['Count']
        elif "cdr3nt" in l.keys():
            v_genestring = l['v']
            j_genestring = l['j']
            cdr3aa = l['cdr3aa']
            cdr3nt = l['cdr3nt']
            clonecount = l['count']
        else:
            print "Bad format detected."
            sys.exit(1)
        outl = {} # ditch the old info
        if 'TRA' in v_genestring or 'TRA' in j_genestring:
            outl['va_genes'] = listsep.join( v_genestring.split(',') )
            outl['ja_genes'] = listsep.join( j_genestring.split(',') )
            outl['cdr3a'] = cdr3aa
            outl['cdr3a_nucseq'] = cdr3nt
            if single_chain=='alpha':
                if organism == "mouse":
                    outl['vb_genes'] = "TRBV19*01"
                    outl['jb_genes'] = "TRBJ1-4*02"
                    outl['cdr3b'] = "CASSMGANERLFF"
                    outl['cdr3b_nucseq'] = "tgtgccagcagtatgggggccaacgaaagattatttttc"
                else:
                    print "error: need to add organism info in script"
                    sys.exit()
#            outl['cdr3a'] = l['CDR3 amino acid sequence']
#            outl['cdr3a_nucseq'] = l['CDR3 nucleotide sequence']
        elif 'TRB' in v_genestring or 'TRB' in j_genestring:
            outl['vb_genes'] = listsep.join( v_genestring.split(',') )
            outl['jb_genes'] = listsep.join( j_genestring.split(',') )
            outl['cdr3b'] = cdr3aa
            outl['cdr3b_nucseq'] = cdr3nt
            if single_chain == 'beta':
                if organism == "mouse":
                    outl['va_genes'] = "TRAV8D-1*02"
                    outl['ja_genes'] = "TRAJ33*01"
                    outl['cdr3a'] = "CATDMDSNYQLIW"
                    outl['cdr3a_nucseq'] = "tgtgctactgacatggatagcaactatcagttgatctgg"
            #outl['cdr3b'] = l['CDR3 amino acid sequence']
            #outl['cdr3b_nucseq'] = l['CDR3 nucleotide sequence']
        else:
            print 'cant determine whether this cdrblast line is alpha or beta:',l
            return None ## NOTE EARLY RETURN

        outl['clone_size'] = clonecount
        if cdrblast_clone_size_warning not in warnings:
            print cdrblast_clone_size_warning
            warnings.add( cdrblast_clone_size_warning )

    return outl




infields = []


print 'making:',output_file

out = open( output_file,'w')
out.write( '\t'.join( outfields ) + '\n' )

idcounter=0
for inline in open( input_file,'rU'):
    line = remove_wonky_characters_from_the_end_of_line( inline )
    if not infields:
        if line[0] == '#':
            line = line[1:]
        infields = line.split('\t')
        assert infields
        ## we may want to add some fields to outfields
    else:
        l = parse_tsv_line( line, infields )
        l = map_field_names( l, input_format, output_format )

        if l is None: continue

        if auto_ids:
            if not id_base:
                id_base = 'id'
            idcounter += 1
            if output_format == "clones":
                l['id'] = '{}{}{}'.format(id_base,idcounter,".clone")
            else:
                l['id'] = '{}{}'.format(id_base,idcounter)

        if epitope: ## danger if we are using this cmdline option as a variable...
            l['epitope'] = epitope
        if subject: ## somewhere else in the script...
            l['subject'] = subject

        badline = False
        for field in outfields:
            if field == "clone_id":
                l['clone_id'] = l['id']
            val = reconstruct_field_from_data( field, l, organism )
            if val is None:
                if single_chain:
                    if single_chain=='alpha' and ( field[:3] in ['vb_','jb_'] or field.startswith('cdr3b') ):
                        val = "-"
                        l[ field ] = val
                        break
                    elif single_chain== 'beta' and ( field[:3] in ['va_','ja_'] or field.startswith('cdr3a') ):
                        val = "-"
                        l[ field ] = val
                        break
                print 'ERROR bad line: field= {} failed, line= {}'.format( field, line[:-1] )
                badline = True
                break
            if field in ["cdr3a", "cdr3b"]:
                if len(l[field]) <= 5:
                    print 'bad line: length of ' + field + " is only " + str(len(l[field]))
                    badline = True
                    break
            l[ field ] = val

        if badline:
            continue
        outl = {}
        for field in outfields:
            if not field in l.keys():
                if "prob" in field:
                    outl[field] = 1.0
                    if prob_warning not in warnings:
                        print prob_warning
                        warnings.add( prob_warning )
                else:
                    outl[field] = "-"
                continue
            outl[ field ] = l[ field ]

        out.write( make_tsv_line( outl, outfields ) + '\n' )

out.close()
