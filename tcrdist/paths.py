import os
import sys

path_to_scripts = os.path.dirname(os.path.realpath(__file__))
assert not path_to_scripts.endswith('/')

##
## the directories db/ and external/ do not live in the github repository
## They will hopefully get created when you run the setup.py script.
##

path_to_db = path_to_scripts+'/db'
assert os.path.isdir( path_to_db )

## used for making nice sortable tables
path_to_tablesorter_files = path_to_scripts+'/external/tablesorter'
assert os.path.isdir( path_to_tablesorter_files)

## NCBI BLAST
path_to_blast_executables = path_to_scripts+'/external/blast-2.2.16/bin'
assert os.path.isdir( path_to_blast_executables )

## handy commandline parser
path_to_blargs = path_to_scripts+'/external/blargs'
assert os.path.isdir( path_to_blargs )

sys.path.append(path_to_blargs)

if __name__ == '__main__':
    for p in sys.path:
        print 'sys.path:',p

