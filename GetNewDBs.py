from Bio import SeqIO
import sys
from subprocess import call
import os.path
from os.path import isfile, join
from os import listdir
from paths import path_to_db, path_to_scripts
from shutil import copyfile
import filecmp

defaultdb = path_to_db
newdbbase = path_to_scripts + "/tempdbs/"
newdbdir = path_to_scripts + "/tempdbs/fasta/"
if not os.path.isdir(newdbbase):
    os.makedirs(newdbbase)
if not os.path.isdir(newdbdir):
    os.makedirs(newdbdir)
if not os.path.isdir(newdbbase + "/probs_files/"):
    os.makedirs(newdbbase + "/probs_files")
    
##Ensure other necessary files are in the new directory.
dbfiles = ["human_rand_divs_new.txt", "mouse_rand_divs_new.txt", "new_nextgen_chains_human_A.tsv", "new_nextgen_chains_human_B.tsv", "new_nextgen_chains_mouse_A.tsv", "new_nextgen_chains_mouse_B.tsv", "nextgen_tuple_counts_v2_human_max10M.log", "nextgen_tuple_counts_v2_mouse_max10M.log", "randhuman_paired.tsv"]
for dbfile in dbfiles:
    if not os.path.isfile(defaultdb + "/" + dbfile):
        print "Error: " +dbfile + " missing from " + defaultdb
        sys.exit(1)
    if not os.path.isfile(newdbbase + "/" + dbfile):
        copyfile(defaultdb + "/" + dbfile, newdbbase + "/" + dbfile)
    if not filecmp.cmp(defaultdb + "/" + dbfile, newdbbase + "/" + dbfile):
        print "WARNING: " + defaultdb + "/" + dbfile + " and " + newdbbase + "/" + dbfile + " are not the same. Attempting to copy default version over. If you have not altered this file, you may want to re-run setup.py"
        copyfile(defaultdb + "/" + dbfile, newdbbase + "/" + dbfile)
probfiles = [pf for pf in listdir(defaultdb +"/probs_files/") if isfile(join(defaultdb +"/probs_files/", pf))]
for pf in probfiles:
    if not os.path.isfile(newdbbase + "/probs_files/" + pf):
         copyfile(defaultdb + "/probs_files/" + pf, newdbbase + "/probs_files/" + pf)
    if not filecmp.cmp(defaultdb + "/probs_files/" + pf, newdbbase + "/probs_files/" + pf):
         print "WARNING: " + defaultdb + "/probs_files/" + pf + " and " + newdbbase + "/probs_files/" + pf + " are not the same. Attempting to copy default version over. If you have not altered this file, you may want to re-run setup.py"
         copyfile(defaultdb + "/probs_files/" + pf, newdbbase + "/probs_files/" + pf)


filelist = ["IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP",
            "IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP",
            "IMGTGENEDB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP"]
nosplitfile = "IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP"
spplist = ["Homo_sapiens", "Mus_musculus"]
mainseglist = ["TRAJ", "TRAV", "TRBJ", "TRBV", "TRDV", "TRDJ", "TRGV", "TRGJ"]

for sp in spplist:
    if os.path.isfile(newdbdir + "imgt_" + sp + "_TR_protein_sequences_with_gaps.fasta"):
        call(['rm',newdbdir + "imgt_" + sp + "_TR_protein_sequences_with_gaps.fasta"])                    

for f in filelist:
    species = {}
    floc = "http://www.imgt.org/download/GENE-DB/" + f
    if os.path.isfile(newdbdir + f):
        print "Replacing old version of " + newdbdir + f
        call(["rm", newdbdir + f])
    call(["wget" , floc, "-P", newdbdir])
    if f == nosplitfile:
        seglist = ["."]
    else:
        seglist = mainseglist
    inf = open(newdbdir + f, "r")
    count = 0
    for record in SeqIO.parse(inf, "fasta"):
        count +=1
        ds = record.description.split("|")[2]
        ds = ds.split("/")[0].split("_")[0]
        ds = "_".join(ds.split())
        if ds not in spplist:
            continue
        if f != nosplitfile:
            record.name = record.name.split("|")[1]
            record.description = record.name
            record.id = record.name
        if ds in species.keys():
            for seg in seglist:
                if seg in record.name:
                    species[ds][seg].append(record)
                elif f == nosplitfile and "TR" in record.name.split("|")[1]:
                    species[ds][seg].append(record)
                else:
                    species[ds]["OTHERSPECIES"].append(record)
        else:
            species[ds] = {}
            species[ds]["OTHERSPECIES"] = []
            for seg in seglist:
                species[ds][seg] = []
                if seg in record.name:
                    species[ds][seg].append(record)
                else:
                    species[ds]["OTHERSPECIES"].append(record)
    inf.close()
    for specie, segs in species.iteritems():
        metarecs = []
        if specie == "Homo_sapiens":
            common = "human"
        elif specie == "Mus_musculus":
            common = "mouse"
        else:
             print "Specie " + specie + " not recognized."
             sys.exit()
        for seg, records in segs.iteritems():
            if f == nosplitfile:
                if seg in seglist or "TR" in seg:
                    fn = newdbdir + "imgt_" + common + "_TR_protein_sequences_with_gaps.fasta"
                    with open(fn, "a") as outfile:
                        SeqIO.write(records, outfile, "fasta")
                    continue
                else:
                    continue
            elif "-nt-WithoutGaps-F+ORF+allP" in f:
                base = "_TR_nucleotide_sequences.fasta." + seg + ".fasta"
            elif "-AA-WithoutGaps-F+ORF+inframeP" in f:
                base = "_TR_protein_sequences.fasta." + seg + ".fasta"
            #elif f == nosplitfile:
            #    base = "_TR_protein_sequences_with_gaps.fasta"
            else:
                print "PROBLEM--exiting"
                sys.exit(1)
            with open(newdbdir + "imgt_" + common + base, "w") as outfile:
                SeqIO.write(records, outfile, "fasta")
            if seg != "OTHERSPECIES":
                metarecs += records
        with open(newdbdir + "imgt_" + common + base.split(".fasta.")[0] + ".fasta", "w") as outmeta:
            SeqIO.write(metarecs, outmeta, "fasta")
