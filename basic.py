## LAZY-- some basic imports and shared params
##
##
from glob import glob
from os import popen, system, chdir, remove, getcwd, mkdir
from os.path import exists, isdir, isfile
import math
from math import floor,sqrt
from sys import stderr,argv,exit
import random
import paths
from blargs import Parser
from parse_tsv import *

############################################################################################
############################################################################################
## some shared TCRdist analysis pipeline parameters

pipeline_params= {
    ## Use the update to probability calculation? (March 2017)
    'new_probs':True,

    ## A length-scale factor that lets us fiddle with the distance calculation
    ## and still get similar clustering and TCRdiv calculations
    ## The history is that with the original TCRdist we used nice round values of 25/50 for
    ## several distance thresholds. When we moved to the updated TCRdist, on average
    ## distances decreased by a factor of 1.355, so as a temporary hack we introduced
    ## this parameter so that clustering and diversity calculations can be run with
    ## default parameters and give similar results.
    #'distance_threshold_25': 25.0 / 1.355,
    'distance_threshold_25': 25.0 / 1.298, ## now moving to a gap penalty in the CDR3 region of 12

    ## the gene sets database file -- should move to make this commandline switchable
    ## this is the file that has the sequences and alignments of all the TR genes
    'db_file':'alphabeta_db.tsv' # db file corresponding to original publication
    #'db_file':'gammadelta_db.tsv' # db file corresponding to gamma deltas, uncomment to analyze them (and comment above)
}


## naming scheme for the gene segment types, occasionally useful for iterating
segtypes_uppercase = ['VA','JA','VB','JB']
segtypes_lowercase = ['va','ja','vb','jb']


############################################################################################
############################################################################################

def path_to_current_db_files():
    """Without the trailing /"""
    db_file = paths.path_to_db+'/'+pipeline_params['db_file']
    assert exists(db_file)
    db_files_dir = db_file+'_files'
    if not exists(db_files_dir):
        mkdir(db_files_dir)
    return db_files_dir




## you could modify this function if you have a different cmdline tool for converting svg to png
## like inkscape or cairosvg
##
def convert_svg_to_png( svgfile, pngfile, verbose=True, allow_missing=False, allow_failure=True ):
    if not isfile(svgfile):
        errmsg = 'Error: convert_svg_to_png: svgfile does not exist: {}'.format(svgfile)
        print errmsg
        Log( errmsg )
        if allow_missing:
            return
        else:
            exit()
    cmd = 'convert {} {}'.format( svgfile, pngfile )
    if verbose:
        print cmd
    system(cmd)

    if isfile( pngfile ):
        ## success
        return

    ## cmdline inkscape
    cmd = 'inkscape --export-png {} {}'.format( pngfile, svgfile )
    if verbose:
        print cmd
    system(cmd)

    if isfile( pngfile ):
        ## success
        return

    ## this is probably a long-shot, but in case inkscape is installed on mac
    inkscape_exe = '/Applications/Inkscape.app/Contents/Resources/bin/inkscape'
    if isfile( inkscape_exe ):
        from os.path import abspath
        svgfile_full = abspath( svgfile )
        pngfile_full = abspath( pngfile )

        cmd = '{} --export-png {} {}'.format( inkscape_exe, pngfile_full, svgfile_full )
        if verbose:
            print cmd
        system(cmd)

        if isfile( pngfile ):
            ## success
            return


    ## another possibility
    cmd = 'rsvg-convert {} -o {}'.format( svgfile, pngfile )
    if verbose:
        print cmd
    system(cmd)

    if isfile( pngfile ):
        ## success
        return


    ## this might also occur if the svgfile were empty...
    errmsg = 'Error: convert command failed: cmd="{}" -- is the "convert" cmdline tool (Imagemagick) installed?'\
             .format( cmd )
    print errmsg
    Log( errmsg )
    if not allow_failure:
        exit()

def get_mean_and_sdev( l ):
    N = len(l)
    assert N>0
    mean = float( sum( l ) ) / N
    sdev = 0.0
    for x in l: sdev += ( x- mean )**2

    sdev = sqrt( float(sdev)/N )
    return mean, sdev

def get_median(l_in):
    l = l_in[:] ##make copy
    l.sort()
    n = len(l)
    if n%2: return l[n/2]  ## n is odd
    else: return 0.5 * ( l[n/2] + l[n/2 - 1 ] ) ## n is even


def Log(s): ## silly legacy helper function
    stderr.write(s)
    if s and not s.endswith('\n'):
        stderr.write('\n')


