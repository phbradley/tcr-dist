import sys
import os.path as op

path_to_scripts = op.dirname(op.realpath(__file__))

##
## the directories db/ and external/ do not live in the github repository
## They will hopefully get created when you run the setup.py script.
##

path_to_db = op.join(path_to_scripts, 'db')
assert op.isdir( path_to_db )

## used for making nice sortable tables
path_to_tablesorter_files = op.join(path_to_scripts, 'external/tablesorter')
assert op.isdir( path_to_tablesorter_files)

## NCBI BLAST
path_to_blast_executables = op.join(path_to_scripts, 'external/blast-2.2.16/bin')
assert op.isdir( path_to_blast_executables )

## handy commandline parser
path_to_blargs = op.join(path_to_scripts, 'external/blargs')
assert op.isdir( path_to_blargs )

sys.path.append(path_to_blargs)

db_file = 'alphabeta_db.tsv' # db file corresponding to original publication
# db_file = gammadelta_db.tsv' # db file corresponding to gamma deltas

