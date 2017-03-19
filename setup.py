##
## Run this script in the main repository directory by typing
##
## python setup.py
##
## To clean up and start again, just remove the db/ and external/ directories and re-run.
##
## If you are getting errors, clean up as described above and re-run, saving the stdout and stderr in logfiles
## and contact pbradley@fredhutch.org for help trouble-shooting.
##
##
from os import popen, system, chdir, mkdir
from os.path import exists, isdir, isfile
from sys import stderr,exit,platform

# I don't know how reliable this is:
mac_osx = ( platform.lower() == "darwin" )

if mac_osx:
    print 'Detected mac_osx operating system -- if not, hardcode mac_osx=False in setup.py'


def download_web_file( address ):
    newfile = address.split('/')[-1]

    if exists(newfile):
        print 'download_web_file: {} already exists, delete it to re-download'
        return

    ## try with wget
    cmd = 'wget '+address
    print cmd
    system(cmd)

    if not exists( newfile ):
        print 'wget failed, trying curl'
        cmd = 'curl -L {} -o {}'.format(address,newfile)
        print cmd
        system(cmd)

    if not exists( newfile ):
        print '[ERROR] unable to download (tried wget and curl) the link '+address



## check for python modules
try:
    import numpy
except:
    print '[ERROR] failed to import numpy'
    exit()

try:
    import scipy
except:
    print '[ERROR] failed to import scipy'
    exit()

try:
    import matplotlib
except:
    print '[ERROR] failed to import matplotlib'
    exit()

try:
    import sklearn
except:
    print """
=============================================================================
=============================================================================
[ERROR]
[ERROR] Failed to import the python module sklearn (scikit-learn)
[ERROR] Some analyses (kernelPCA plots, adjusted_mutual_information) will fail
[ERROR] Take a look at http://scikit-learn.org/stable/install.html
[ERROR]
=============================================================================
=============================================================================
"""
    #exit() ## not exiting since most stuff will probably still work...




## setup the
print 'Making the ./external/ directory'
external_dir = 'external/'
if not isdir( external_dir ):
    mkdir( external_dir )

chdir( external_dir )

## download blast
blastdir = './blast-2.2.16'

if not isdir( blastdir ):
    if mac_osx:
        address = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/2.2.16/blast-2.2.16-universal-macosx.tar.gz'
    else:
        address = 'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/2.2.16/blast-2.2.16-x64-linux.tar.gz'

    tarfile = address.split('/')[-1]
    if not exists( tarfile ):
        print 'Downloading a rather old BLAST tool'
        download_web_file( address )

        if not exists( tarfile ):
            print '[ERROR] download BLAST failed!'
            exit()

    cmd = 'tar -xzf '+tarfile
    print cmd
    system(cmd)


## download other db files
address = 'https://www.dropbox.com/s/6go61t51vpn2xiu/tcrdist_extras_v1.tgz'
tarfile = address.split('/')[-1]

if not exists( tarfile ):
    print 'Downloading database files'
    download_web_file( address )

    if not exists( tarfile ):
        print '[ERROR] download database files failed'
        exit()

## md5sum check
lines = popen('md5sum '+tarfile).readlines()
if lines and len(lines[0].split()) == 2:
    # rhino1 public_release$ md5sum tcrdist_extras_v1.tgz
    # 3f3d3ad768f1c68b02847696255d89ad  tcrdist_extras_v1.tgz
    checksum = lines[0].split()[0]
    expected_checksum = '3f3d3ad768f1c68b02847696255d89ad'
    if checksum == expected_checksum:
        print "md5sum checksum for tarfile matches expected..."
    else:
        print "[ERROR] OH NO! md5sum checksum for tarfile does not match: actual={} expected={}"\
            .format( checksum, expected_checksum )
else:
    print '[WARNING] md5sum command failed or gave unparseable output, unable to check the tarfile...'



download_dir = tarfile[:-4]
if not isdir( download_dir ):
    cmd = 'tar -xzf '+tarfile
    print cmd
    system(cmd)

    if not isdir( download_dir ):
        print '[ERROR] tar failed or the database download was corrupted!'
        exit()

cmd = 'mv {}/external/* .'.format(download_dir)
print cmd
system(cmd)

cmd = 'mv {}/db ../'.format(download_dir)
print cmd
system(cmd)

cmd = 'mv {}/datasets ../'.format(download_dir)
print cmd
system(cmd)

cmd = 'mv {}/testing_ref ../'.format(download_dir)
print cmd
system(cmd)

