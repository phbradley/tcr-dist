import sys
from glob import glob
import logging
import os.path as op

import logo_tools
from genetic_code import genetic_code
import cdr3s_human
from paths import path_to_db, path_to_blast_executables

AB = 'AB'

bases_plus = logo_tools.nucleotide_classes_lower_case.keys()

for a in bases_plus:
    for b in bases_plus:
        for c in bases_plus:
            codon = a+b+c
            if codon in genetic_code: continue

            aas = []
            for a1 in logo_tools.nucleotide_classes_lower_case[a]:
                for b1 in logo_tools.nucleotide_classes_lower_case[b]:
                    for c1 in logo_tools.nucleotide_classes_lower_case[c]:
                        aas.append( genetic_code[ a1+b1+c1 ] )
            if min(aas) == max(aas):
                genetic_code[codon] = aas[0]
            else:
                genetic_code[codon] = 'X'

def get_translation( seq, frame ):
    assert frame[0] in '+-'
    if frame[0] == '-': seq = logo_tools.reverse_complement( seq )
    offset = abs( int( frame ))-1
    assert offset in range(3)
    seq = seq[offset:].lower()
    naa = len(seq)/3
    protseq = ''
    codons = []
    for i in range(naa):
        codon = seq[3*i:3*i+3]
        codons.append( codon )
        if '#' in codon:
            protseq += '#'
        else:
            protseq += genetic_code.get( codon, 'X' )
    return protseq,codons


all_offsets = {}
all_fasta = {}

for organism in ['mouse','human']:
    all_offsets[organism] = {}
    all_fasta[organism] = {}
    prot = 'protein'
    nuc = 'nucleotide'
    for ab in AB:
        all_offsets[organism][ab] ={}
        all_fasta[organism][ab] ={}
        for vj in 'VJ':
            all_offsets[organism][ab][vj] ={}
            all_fasta[organism][ab][vj] ={}
            myfasta = {}
            all_fastafile = {}
            for np in [prot,nuc]:
                myfasta[np] ={}
                fastafile = op.join(path_to_db, '/fasta/imgt_%s_TR_%s_sequences.fasta.TR%s%s.fasta'\
                            %( organism, np, ab, vj ))
                all_fastafile[np] = fastafile
                assert exists( fastafile )

                ## make sure the dbs are here
                dbfiles = glob('{}.{}*'.format(fastafile,np[0]))
                if not dbfiles:
                    cmd = '{} -p {} -i {}'.format( op.join(path_to_blast_executables, 'formatdb'), 'T' if np==prot else 'F',
                                                            fastafile )
                    logging.info(cmd)
                    system(cmd)
                    dbfiles = glob('{}.{}*'.format(fastafile,np[0]))
                    if not dbfiles:
                        logging.error('blast db creation failed!')
                        exit()

                id = ''
                for line in open( fastafile,'r'):
                    if line[0] == '>':
                        id = line[1:-1]
                        myfasta[np][id] = ''
                    else:
                        assert id
                        myfasta[np][id] += line[:-1]
            all_fasta[organism][ab][vj] = myfasta
            for id in myfasta[prot]:
                assert id in myfasta[nuc]
                pseq = myfasta[prot][id]
                nseq = myfasta[nuc][id]
                myframe = -1
                for i in range(3):
                    tseq = get_translation( nseq, '+{}'.format(i+1) )[0]
                    if pseq in tseq:
                        myframe = i + 3*tseq.index(pseq)
                assert myframe >= 0
                all_offsets[organism][ab][vj][id] = myframe

def reverse_q2hmap( qseq, hseq, hit ):
    assert hit.h_strand == -1
    q2hmap = {}
    for qpos,(hpos,na) in hit.q2hmap.iteritems():
        if hpos>=0:
            new_qpos = len(qseq)-qpos-1
            new_na = logo_tools.base_partner[ na ]
            assert new_na == hseq[ hpos ]
            q2hmap[new_qpos] = (hpos,new_na)
    return q2hmap

## should be just a single query sequence
def get_all_hits_with_evalues_and_scores( blastfile ):
    query = ''
    in_hits = False
    hits = []
    for line in open( blastfile,'r'):
        if line.startswith('Query='):
            assert not query
            query= line.split()[1]
        elif line.startswith('Sequences producing'):
            in_hits = True
        elif line.startswith('>'):
            assert in_hits
            break
        elif in_hits:
            l = line.split()
            if len(l) >= 3:
                evalue = l[-1]
                if evalue.startswith('e'): evalue = '1'+evalue
                bitscore = int( l[-2] )
                hitid = l[0]
                if '|' in hitid: hitid = hitid.split('|')[1]
                hits.append( ( hitid, bitscore, float(evalue) ) )
                #print 'new_hit:',hits[-1]
    return hits

