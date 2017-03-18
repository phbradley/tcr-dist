# from basic import *
# import math
import amino_acids
# import random
import sys
#from operator import add

# weblogo_paths = False
# weblogo_exe = '/home/pbradley/download/weblogo/weblogo/seqlogo'

def get_alphabet( pwm ):
    alphabet = pwm[0].keys()[:]
    alphabet.sort()
    return alphabet

def check_pwm( pwm, tol= 0.001 ):
    L = len(pwm)
    alphabet = get_alphabet( pwm )
    for pos in range(L):
        for aa,val in pwm[pos].iteritems(): assert val > -1e-6
        total = sum( ( pwm[pos][aa] for aa in alphabet ) )
        assert abs( total - 1.0 ) - tol

# def boost_pwm_column( col, exponent ):
#     sum = 0.0
#     for a in col:
#         col[a] = ( col[a] ) ** exponent
#         sum += col[a]
#     for a in col:
#         col[a] /= sum
#     return

# def boost_pwm( pwm, exponent ):
#     for i in pwm:
#         boost_pwm_column( pwm[i], exponent )
#     return



def create_protein_pwm_from_sequences( seqs, pseudocounts = 0.0 ):
    return create_pwm_from_sequences( seqs, amino_acids.amino_acids, pseudocounts=pseudocounts )

def create_dna_pwm_from_sequences( seqs, pseudocounts=0.0 ):
    return create_pwm_from_sequences( seqs, list('acgt'), pseudocounts=pseudocounts )

def create_pwm_from_sequences( seqs, alphabet, pseudocounts=0.0 ):
    pwm = {}
    if len(seqs) == 0: return pwm
    L = len( seqs[0] )

    for pos in range(L):
        pwm[ pos ] = dict( zip( alphabet, [pseudocounts]*len(alphabet) ) )

    for s in seqs:
        assert len(s) == L
        for pos in range(L):
            if s[pos] not in alphabet:
                sys.stderr.write('logo_tools.create_pwm_from_sequences: skipping bad character %s\n'%(s[pos]))
                continue
            pwm[ pos ][ s[pos] ] += 1

    for pos in range(L):
        norm = 1.0 / sum( pwm[pos].values() )
        for a in alphabet: pwm[ pos ][ a ] *= norm

    check_pwm( pwm )

    return pwm


# def write_pwm_to_file( pwm, filename ):
#     check_pwm( pwm, tol = 0.01 )
#     out = open( filename, 'w' )
#     alphabet = get_alphabet( pwm )
#     out.write( ' '.join( alphabet )+'\n')
#     for pos in range(len(pwm)):
#         out.write( ' '.join( [ '%.6f'%(pwm[pos][x]) for x in alphabet ] ) +'\n' )
#     out.close()

# def read_pwm_from_file( filename ):
#     lines = map( string.split, open( filename, 'r' ).readlines() )
#     alphabet = lines[0]
#     pos = 0
#     pwm = {}
#     for line in lines[1:]:
#         pwm[pos] = {}
#         assert len(line) == len( alphabet )
#         total = 0.0
#         for i,a in enumerate( alphabet ):
#             pwm[pos][a] = float( line[i] )
#             total += pwm[pos][a]
#         for a in alphabet: pwm[pos][a] /= total

#         pos += 1
#     return pwm


# def weblogo_pwm_to_file( pwm, outfile, stretch = False, show_xaxis= True, number_fontsize = 8, show_yaxis = True ):
#     global weblogo_paths
#     if not weblogo_paths: ## only do this once, yes I know it's super-hacky; this code is not used by tcr stuff
#         #sys.path.append( '/home/pbradley/download/weblogo-3.3/' )
#         sys.path.append( '/home/pbradley/download/weblogo-3.0/' )
#         sys.path.append( '/home/pbradley/local/python_packages/lib/python2.6/site-packages/')
#         sys.path.append( '/home/pbradley/local/python_packages/lib64/python2.6/site-packages/' )
#         weblogo_paths = True

#     import numpy as np
#     import weblogolib
#     import corebio
#     first_index = 1
#     ## create a LogoData object
#     lcounts = []
#     L = len( pwm )
#     #alphabet = ['a','c','g','t'] ## has to match the order in unambiguous_dna_alphabet
#     pwm_alphabet = get_alphabet( pwm )

#     pwm_is_dna = True
#     pwm_is_protein = False
#     if len( pwm_alphabet ) == 4:
#         alphabet = corebio.seq.unambiguous_dna_alphabet
#         pwm_is_dna = True
#         pwm_is_protein = False
#     else:
#         assert len( pwm_alphabet ) == 20
#         alphabet = corebio.seq.unambiguous_protein_alphabet
#         pwm_is_dna = False
#         pwm_is_protein = True
#     for i in range(L):
#         col = []
#         for a in alphabet:
#             if a in pwm_alphabet: aa = a
#             else: aa = a.lower()
#             col.append( pwm[i][aa] )
#         lcounts.append( col )
#     counts = np.array( lcounts )

#     data = weblogolib.LogoData.from_counts( alphabet, counts )
#     #maxent = math.log( 20.0 )

#     if stretch:
#         for i in range(L):
#             data.entropy[i] = math.log( len(alphabet) ) #20.0 )

#     #for i in range(L):
#     #    data.entropy[i] = jsd_values[i]
#     #    assert jsd_values[i] < max_divergence + 1e-3

#     options = weblogolib.LogoOptions()
#     #options.title = "A Logo Title"
#     #import colorscheme
#     if pwm_is_dna:
#         options.color_scheme = weblogolib.classic
#     else:
#         options.color_scheme = weblogolib.std_color_schemes[ "chemistry" ]

#     #options.size = weblogolib.std_sizes[ "medium" ]
#     #options.size = weblogolib.std_sizes[ "large" ]
#     #options.size = weblogolib.std_sizes[ "huge" ]
#     sizescale = 8
#     options.size = weblogolib.LogoSize( stack_width = 5.4*sizescale,  stack_height = 5.4*sizescale*2 )
#     # options.size = weblogolib.LogoSize( stack_width = 5.4*sizescale,  stack_height = 5.4*sizescale*5 )
#     #options.creator_text = "",
#     options.show_fineprint = False
#     options.number_interval = 1
#     options.resolution = 305

#     options.show_xaxis = show_xaxis
#     options.show_yaxis = show_yaxis
#     options.logo_margin = 1 ## default is 2
#     options.first_index = first_index
#     options.number_fontsize = number_fontsize

#     options.stacks_per_line = 100

#     #options.show_yaxis = True
#     #options.rotate_numbers = True
#     #options.yaxis_label = 'JSD'
#     #options.unit_name = 'jsd'
#     #options.yaxis_scale = max_divergence
#     #options.yaxis_tic_interval = 0.05


#     format = weblogolib.LogoFormat(data, options)

#     fout = open(outfile, 'w')
#     outfiletype = outfile.split('.')[-1]
#     if outfiletype == 'eps':
#         weblogolib.eps_formatter( data, format, fout)
#     elif outfiletype == 'png':
#         weblogolib.png_formatter( data, format, fout)
#     elif outfiletype == 'jpeg':
#         weblogolib.jpeg_formatter( data, format, fout)
#     else:
#         print 'unrecognized image type!',outfiletype, outfile
#         exit()

#     fout.close()


# def make_logo_sequences_from_pwm( pwm, suggested_nlogoseqs ):
#     check_pwm( pwm )

#     nlogoseqs = suggested_nlogoseqs

#     big_counter = 0
#     while True:
#         big_counter += 1
#         sequences = [ '' ] * nlogoseqs
#         L = len( pwm )

#         any_pos_failed = False
#         for pos in range(L):

#             prob = pwm[ pos ]
#             ## now try to get the right number...
#             target = nlogoseqs
#             stepsize = 0.01
#             last_count = nlogoseqs
#             counter = 0
#             this_pos_failed = True
#             while True:
#                 count = 0
#                 for aa in prob:
#                     nseq = int( math.floor( 0.5 + target * prob[aa] ) )
#                     count += nseq
#                 if count  == nlogoseqs:
#                     this_pos_failed = False
#                     break
#                 elif count < nlogoseqs: target += stepsize
#                 elif count > nlogoseqs: target -= stepsize
#                 if ( count < nlogoseqs and last_count > nlogoseqs ) or ( count > nlogoseqs and last_count < nlogoseqs ):
#                     counter += 1
#                     if counter%50000==0:print 'count=',count,'target=',target,'pos=',pos
#                     stepsize = 0.1 / counter
#                     if counter > 5000: break
#                 last_count = count

#             if this_pos_failed:
#                 any_pos_failed = True
#                 break
#             count = 0
#             for aa in prob:
#                 nseq = int( math.floor( 0.5 + target * prob[aa] ) )
#                 for seq in range( count, count+nseq ):
#                     sequences[ seq ] += aa
#                 count += nseq
#             assert count == nlogoseqs
#         if any_pos_failed:
#             nlogoseqs += 1 ## try again with a different number of logo sequences
#             if not big_counter%10: print 'Try again: suggested_nlogoseqs=',suggested_nlogoseqs,'nlogoseqs=',nlogoseqs
#         else: break

#     return sequences


# def get_prob_from_aa_scores( aa_scores_bound, aa_scores_unbound, percentile, temperature, min_index, sdevs,
#                              average_decoy_scores = False, min_index_penalty = 0.0, verbose = False ):
#     avg_score = {}
#     for aa in sorted( aa_scores_bound.keys() ): ## note that aa might be a pair or triple of bases...

#         total_scores = aa_scores_bound[ aa ]
#         peptide_scores = aa_scores_unbound[ aa ]

#         total_scores.sort()
#         peptide_scores.sort()

#         n = len( total_scores )
#         assert n == len( peptide_scores )

#         score_index = ( n * percentile ) / 100

#         penalty = 0.0
#         if score_index < min_index:
#             assert percentile <= 50
#             new_score_index = min( min_index, n / 2 )
#             Log('aa= %s n= %d score_index < min_index: %d < %d ... using %d'\
#                     %(aa,n,score_index,min_index,new_score_index))
#             score_index = new_score_index
#             penalty = min_index_penalty

#         if average_decoy_scores:
#             avg_total_score   = reduce( add, total_scores  [ : score_index+1 ] ) / ( score_index+1 )
#             avg_peptide_score = reduce( add, peptide_scores[ : score_index+1 ] ) / ( score_index+1 )
#         else:
#             avg_total_score = total_scores[ score_index ]
#             avg_peptide_score = peptide_scores[ score_index ]


#         ## err estimate
#         errcount=0
#         bound_err = 0.0
#         unbound_err = 0.0
#         if score_index>0:
#             errcount+=1
#             bound_err   += abs(   total_scores[ score_index ] -   total_scores[ score_index-1 ] )
#             unbound_err += abs( peptide_scores[ score_index ] - peptide_scores[ score_index-1 ] )
#         if score_index < len(total_scores)-1:
#             errcount+=1
#             bound_err   += abs(   total_scores[ score_index ] -   total_scores[ score_index+1 ] )
#             unbound_err += abs( peptide_scores[ score_index ] - peptide_scores[ score_index+1 ] )

#         if errcount:
#             bound_err /= errcount
#             unbound_err /= errcount

#         avg_score[ aa ] = avg_total_score - avg_peptide_score + penalty

#         sdevs[aa] = bound_err + unbound_err

#         if verbose:
#             print '%s %9.3f %9.3f %9.3f %3d %4d bound_err: %5.2f unbound_err: %5.2f'\
#                 %( aa, avg_score[aa], avg_total_score, avg_peptide_score, score_index, n,
#                    bound_err,unbound_err )

#         ## for diagnostics, compute mean+sdev of scores
#         # P_remove = 20 ## remove 20 percent worst outliers -- same as in job21 on hyrax

#         # n_remove = (P_remove * len( total_scores ) ) / 100
#         # for i in range(n_remove):
#         #     del total_scores[-1]
#         #     del peptide_scores[-1]

#         # bound_mean = reduce( add, total_scores ) / len( total_scores )
#         # bound_sdev = math.sqrt( reduce( add, [ ( x - bound_mean ) ** 2 for x in total_scores ] ) /
#         #                         len( total_scores ) )

#         # unbound_mean = reduce( add, peptide_scores ) / len( peptide_scores )
#         # unbound_sdev = math.sqrt( reduce( add, [ ( x - unbound_mean ) ** 2 for x in peptide_scores ] ) /
#         #                           len( peptide_scores ) )

#         # sdevs[ aa ] = ( bound_sdev, unbound_sdev )


#     ## now get probabilities by boltzmann averaging
#     mn = min( avg_score.values() )
#     Z = reduce( add, [ math.exp( ( mn - avg_score[ aa ] ) / temperature ) for aa in avg_score ] )
#     prob = {}
#     for aa in avg_score:
#         prob[ aa ] = math.exp( ( mn - avg_score[ aa ] ) / temperature ) / Z

#     return prob


# def make_png_file_from_pwm( pwm, pngfile, stretch=False, x_axis_numbering = False, width=40, height=12 ):

#     suggested_nlogoseqs = 1000
#     sequences = make_logo_sequences_from_pwm( pwm, suggested_nlogoseqs )

#     seqfile = pngfile+'.tmpseqs'
#     out = open(seqfile,'w')
#     out.write( '\n'.join( sequences )+'\n' )
#     out.close()

#     alphabet = pwm[0].keys()
#     if len(alphabet) in [4,5]:
#         kstring = '-k 1'
#     elif len(alphabet) in [20,21]:
#         kstring = '-k 0'
#     else:
#         kstring = ''

#     cmd = '%s -f %s -c %s %s %s -F PNG -w %d -h %d > %s'\
#         %(weblogo_exe,seqfile,kstring,' -S '*stretch,' -n '*x_axis_numbering,width,height,pngfile)
#     #print cmd
#     system(cmd)

#     remove( seqfile )


# def kullback_leibler_divergence( p, q ):
#     assert len(p) == len(q)
#     div = 0.0
#     for a in p:
#         if p[a] != 0.0 and q[a] != 0.0:
#             div += p[a] * math.log( p[a] / q[a] )
#     return div

# def symmetric_kullback_leibler_divergence( p, q ): ## average of both directions
#     assert len(p) == len(q)
#     div = 0.0
#     for a in p:
#         if p[a] != 0.0 and q[a] != 0.0:
#             div += 0.5 * ( p[a] * math.log( p[a] / q[a] ) + q[a] * math.log( q[a] / p[a] ) )
#     return div

# def safe_symmetric_kullback_leibler_divergence( p, q ): ## sum of both directions
#     epsilon = 0.005
#     assert len(p) == len(q)
#     div = 0.0
#     for a in p:
#         pa = max( p[a], epsilon )
#         qa = max( q[a], epsilon )
#         div += pa * math.log( pa / qa ) + qa * math.log( qa / pa )
#     return div

# def jensen_shannon_divergence( p, q ):
#     assert len(p) == len(q)
#     pq = {}
#     for a in p: pq[a] = 0.5 * ( p[a] + q[a] )

#     return 0.5 * kullback_leibler_divergence( p, pq ) + 0.5 * kullback_leibler_divergence( q, pq )

# def blic_similarity( p, q ):
#     uniform_prob = 1.0 / len(p)
#     bg = {}
#     pq = {}
#     for a in p:
#         bg[a] = uniform_prob
#         pq[a] = 0.5 * ( p[a] + q[a] )
#     return -1.0 * jensen_shannon_divergence( p, q ) + jensen_shannon_divergence( pq, bg )

# def IC_of_column( col ):
#     alpha = col.keys()
#     ic = math.log(len(alpha),2)
#     for a in alpha:
#         if col[a]>1e-6:
#             ic += col[a] * math.log( col[a], 2.0 )
#     return ic

# def AAD_between_columns( ecol, pcol ):
#     alpha = ecol.keys()
#     assert len(pcol.keys()) == len(alpha)
#     aad = 0.0
#     for a in alpha:
#         aad += abs( ecol[a] - pcol[a] )
#     aad /= float(len(alpha))
#     return aad

# def rank_top_between_columns( ecol, pcol ):
#     elist = [(ecol[x],x) for x in ecol ]
#     plist = [(pcol[x],x) for x in pcol ]
#     (p,top_exp_aa) = max( elist )
#     plist.sort()
#     plist.reverse()
#     for i,(p,aa) in enumerate( plist ):
#         if aa == top_exp_aa:
#             return i+1 ## 1-indexed
#     assert False
#     return 0

# def AUC_between_columns( ecol, pcol ):
#     alpha = ecol.keys()
#     threshold = 0.1
#     positives = []
#     negatives = []
#     for a in alpha:
#         if ecol[a] >= threshold: positives.append( a )
#         else: negatives.append( a )
#     if not positives:
#         return -1
#     l = [ (pcol[x],x) for x in alpha ]
#     l.sort()
#     l.reverse()
#     TP=0
#     auc = 0.0
#     num_positives = len(positives)
#     num_negatives = len(negatives)
#     for (p,aa) in l:
#         if aa in positives:
#             TP += 1
#         else:
#             auc += float(TP) / ( num_positives * num_negatives )
#     return auc


# def random_pwm_column( alphabet ):
#     col = {}
#     total = 0.0
#     for a in alphabet:
#         col[a] = random.random()
#         total += col[a]

#     for a in alphabet: col[a] /= total

#     return col

# def blic_similarity_pvalue( prediction, experiment, n_iter ):
#     blic = blic_similarity( prediction, experiment )
#     nbetter = 0
#     alphabet = prediction.keys()
#     for i in range(n_iter ):
#         random_p = random_pwm_column( alphabet )
#         random_blic = blic_similarity( random_p, experiment )
#         if random_blic > blic: nbetter += 1
#     return float( nbetter ) / n_iter


# def jensen_shannon_pwm_divergence( p, q ):
#     L = len(p)
#     assert L == len(q)
#     div = 0.0
#     for pos in range( L ):
#         div += jensen_shannon_divergence( p[pos], q[pos] )
#     return div / L

# def symmetric_kullback_leibler_pwm_divergence( p, q ):
#     L = len(p)
#     assert L == len(q)
#     div = 0.0
#     for pos in range( L ):
#         div += symmetric_kullback_leibler_divergence( p[pos], q[pos] )
#     return div / L

# def blic_pwm_similarity( p, q ):
#     L = len(p)
#     assert L == len(q)
#     sim = 0.0
#     for pos in range( L ):
#         sim += blic_similarity( p[pos], q[pos] )
#     return sim / L

base_partner = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n',
                'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                'R': 'Y', 'Y': 'R',
                'S': 'S', 'W': 'W',
                'K': 'M', 'M': 'K',
                '.': '.' }

nucleotide_classes_lower_case = { 'a':'a',
                                  'c':'c',
                                  'g':'g',
                                  't':'t',
                                  'w':'at',
                                  's':'cg',
                                  'k':'gt',
                                  'm':'ac',
                                  'y':'ct',
                                  'r':'ag',
                                  'n':'acgt' }

nuc_match_lower_case = {}

for x in nucleotide_classes_lower_case:
    for y in nucleotide_classes_lower_case:
        nuc_match_lower_case[(x,y)] = False
        for x1 in nucleotide_classes_lower_case[x]:
            for y1 in nucleotide_classes_lower_case[y]:
                if x1==y1:
                    nuc_match_lower_case[(x,y)] = True
                    break

def nucleotide_symbols_match( a_in, b_in ):
    a = a_in.lower()
    b = b_in.lower()
    if a==b: return True
    return nuc_match_lower_case.get( (a,b), False )


# def preferred_base( col ):
#     return max( [ ( col[x], x ) for x in col ] )[1]


def reverse_complement( seq ):
    newseq = ''
    L = len(seq)
    for pos in range( L-1, -1, -1 ):
        newseq += base_partner[ seq[ pos ] ]
    assert len( newseq ) == L
    return newseq

# def reverse_complement_pwm( pwm ):
#     global base_partner
#     revpwm = {}
#     L = len( pwm )
#     for pos in range( L ):
#         revpos = L - 1 - pos
#         revpwm[ revpos ] = {}
#         for a in pwm[ pos ]:
#             revpwm[ revpos ][ base_partner[ a ] ] = pwm[ pos ][ a ]
#     return revpwm

# def trim_pwm( pwm, ntrim, ctrim ): ## now can handle ntrim or ctrim less than 0
#     newpwm = {}
#     L = len(pwm)

#     bgpwm = {}
#     for a in pwm[0]:
#         bgpwm[a] = 1.0 / len( pwm[0] )

#     for pos in range( L - ntrim - ctrim ):
#         if pos+ntrim < 0 or pos + ntrim >= L:
#             newpwm[ pos ] = deepcopy( bgpwm )
#         else:
#             newpwm[ pos ] = deepcopy( pwm[ pos + ntrim ] )

#     return newpwm

# def optimal_sequence( pwm ):
#     optseq = ''
#     for pos in range(len(pwm)):
#         optseq += max( [ ( pwm[pos][x], x) for x in pwm[pos] ] )[1]
#     return optseq



# """
# Gene:  Mig1-primary  Motif:  AA.GCGGGG  Enrichment Score:  0.496895755125614
# A:      0.227803086037971       0.304418113260319       0.44776757594722        0.566742600501985       0.430392652747271       0.337945502146704       0.710375457778704       0.766379107186479    0.0707446974457731      0.0062051023266637      0.00435938817004998     0.00148409972114262     0.0028021622810555      0.0522110731507016      0.00694422229983977     0.0808012505318237   0.45911214105643        0.306829780079599       0.268049829465637       0.345512634189367       0.174044615110447
# C:      0.188417132500032       0.162314633863411       0.145438889197767       0.0382188719267571      0.0496972441933571      0.0991614951822804      0.0665402611695172      0.0422863400965448   0.132213135735347       0.0420472789205912      0.779220922990459       0.00451917770380136     0.00275486372527724     0.00117433802545397     0.00137356133144941     0.0148304505793713   0.0375183121749853      0.327750484316737       0.225856837565124       0.11015988991179        0.0741805792935567
# G:      0.315278171588745       0.219630246033875       0.205050099448902       0.138191914259755       0.253399955191042       0.178927492767318       0.141406870446047       0.129483396934758    0.0584035327065899      0.947652419039066       0.0035197477195064      0.989896859524604       0.992662066214169       0.945895796256702       0.988431699990139       0.605540725689319    0.447558108311111       0.187091371555164       0.150269710331836       0.0620633886534441      0.401765283318147
# T:      0.268501609873252       0.313637006842396       0.201743435406111       0.256846613311503       0.26651014786833        0.383965509903697       0.0816774106057316      0.0618511557822183   0.73863863411229        0.00409519971367926     0.212899941119985       0.00409986305045171     0.0017809077794988      0.000718792567142971    0.00325051637857227     0.298827573199486    0.0558114384574731      0.1783283640485 0.355823622637403       0.482264087245399       0.350009522277849

# """
# def read_bulyk_pwm( filename ):
#     all_lines = map( string.split, open( filename, 'r' ).readlines() )
#     lines = []
#     for line in all_lines:
#         if line:
#             if len(line)==2 and line[0] == 'Probability':
#                 print 'read_bulyk_pwm: skipping previous lines:',line
#                 lines = []
#             elif len(line[0]) ==2 and line[0][1] == ':':
#                 lines.append(line)
#             elif 'Enrichment' not in line:
#                 id = line[0].split('_')[0]
#                 if not filename.count(id):
#                     print 'read_bulyk_pwm: bad line',line

#     if len(lines)>4:
#         print 'read_bulyk_pwm: skipping lines:',''.join([ ' '.join(x) for x in lines[5:] ] )
#         lines = lines[:4]
#     assert len(lines) == 4
#     pwm = {}
#     L = 0
#     for line in lines:
#         a = line[0][:-1].lower()
#         if L:
#             assert len(line) == L+1
#         else:
#             L = len(line) - 1
#             for pos in range(L): pwm[pos] = {}
#         for pos in range(L): pwm[pos][a] = float( line[pos+1] )
#     return pwm


# def relative_entropy_to_background( pwm ):
#     alphabet = get_alphabet( pwm )
#     background = 1.0 / len( alphabet )
#     relent = 0.0
#     for i in range( len(pwm )):
#         for a in alphabet:
#             if pwm[i][a] > 1e-6:
#                 relent += pwm[i][a] * math.log( pwm[i][a] / background )
#             else:
#                 msg = 'pwm[i][a] out of range: %9.3f'%pwm[i][a]
#                 print msg
#                 Log(msg)

#     return relent


# def matrix2string( pwm ):
#     alphabet = get_alphabet( pwm )
#     L = len(pwm)
#     s = 'L: %d alphabet: %s'%(L,''.join(alphabet))
#     for i in range(L):
#         s+= ' %d '%(i+1)+ ' '.join( [ '%.3f'%(pwm[i][x]) for x in alphabet ] ) ## note 1-indexing for output
#     return s


# def consensus_sequence( pwm ):
#     alphabet = get_alphabet( pwm )
#     assert len( alphabet ) == 4  ## the thresholds here only make sense for DNA
#     consensus = ''
#     for i in range( len( pwm )):
#         (p,a) = max( [ (pwm[i][x],x) for x in alphabet ] )
#         if p >0.75: consensus += a.upper()
#         elif p>0.5: consensus += a.lower()
#         else: consensus += 'n'
#     return consensus

# def create_pseudocounts_pwm( pwm_in, pwm_nseq, pseudocounts ):
#     pwm = {}
#     for i in range( len( pwm_in ) ):
#         pwm[i] = {}
#         total = 0.0
#         for a in pwm_in[i]:
#             pwm[i][a] = pwm_in[i][a] * pwm_nseq + pseudocounts
#             total += pwm[i][a]
#         for a in pwm_in[i]: pwm[i][a] /= total
#     return pwm

# def matrix_score( pwm, sequence ):
#     alphabet = get_alphabet( pwm )
#     L = len( pwm )
#     background = 1.0 / len( alphabet )
#     assert L == len(sequence)
#     score = 0.0
#     for i,base in enumerate( sequence):
#         score += math.log( pwm[i][ base ] / background )
#     return score


# def find_best_match( pwm_in, pwm_nseq, sequence, tag = '', pseudocounts = 0.25 ):
#     alphabet = get_alphabet( pwm_in )
#     L = len( pwm_in )
#     background = 1.0 / len( alphabet )
#     pwm = create_pseudocounts_pwm( pwm_in, pwm_nseq, pseudocounts )

#     relent = relative_entropy_to_background( pwm )

#     ## what range of offsets do we need
#     seqlen = len( sequence )
#     minl = min( L, seqlen )

#     min_overlap = min( minl-2, ( 4*minl )/5 )

#     best_score = -999
#     for seqstart in range( -seqlen+1, L-1 ):
#         ## whats the overlap for this offset
#         ## pwm goes from 0        --> L-1
#         ## seq goes from seqstart --> seqstart + seqlen - 1
#         overlap = 0
#         score = 0.0
#         for seqpos in range(seqlen):
#             a = sequence[ seqpos ]
#             pwmpos = seqpos + seqstart
#             if pwmpos >=0 and pwmpos < L:
#                 overlap += 1
#                 score += math.log( pwm[pwmpos][ a ] / background )
#         if overlap >= min_overlap:
#             if score > best_score:
#                 best_score = score
#                 best_seqstart = seqstart

#     consensus = consensus_sequence( pwm )

#     alcons = ''
#     alseq = ''
#     if best_seqstart <0: alcons = ' '*(-1*best_seqstart)
#     else: alseq += ' '*best_seqstart

#     print 'best_score: seqstart= %4d score= %9.3f relent= %9.3f %4d seqlen= %4d L= %4d %s  cons= %s seq= %s'\
#         %(best_seqstart, best_score, relent, pwm_nseq, seqlen, L, tag, consensus, sequence )
#     print 'best_align1 seqstart= %4d score= %9.3f relent= %9.3f %4d seqlen= %4d L= %4d %s  align= %s%s'\
#         %(best_seqstart, best_score, relent, pwm_nseq, seqlen, L, tag, alcons, consensus )
#     print 'best_align2 seqstart= %4d score= %9.3f relent= %9.3f %4d seqlen= %4d L= %4d %s  align= %s%s'\
#         %(best_seqstart, best_score, relent, pwm_nseq, seqlen, L, tag, alseq, sequence )


#     return best_seqstart, best_score, relent



