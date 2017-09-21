##
## Run this script in the main repository directory by typing
##
## python setup.py
##
## If you are getting errors, re-run the script saving the output
## and contact pbradley@fredhutch.org for help trouble-shooting.
##
##
from os import popen, system, chdir, mkdir
from os.path import exists, isdir, isfile
from sys import stderr,exit,platform

msg = """
This script will download a set of compatible BLAST executables, parameter and
database files used by the tcr-dist pipeline, and some TCR datasets for
testing and analysis. It should be run in the main tcr-dist/ directory.

Altogether, it will end up taking about 500 Megabytes of space.

(To reduce this a bit you can delete the .tgz files in external/ after
it completes successfully.)

Do you want to proceed? [Y/n] """

ans = raw_input(msg)

if ans and ans not in 'Yy':
    print 'Setup aborted.'
    exit()


old_directories = ['db','external','datasets','testing_ref']
found_old_directory = False
for d in old_directories:
    if exists(d):
        found_old_directory = True

if found_old_directory:
    msg = """
It looks like you have some old directories from a previous setup.

I need to remove db/ external/ datasets/ and testing_ref/

Is that OK? [Y/n] """
    ans = raw_input(msg)

    if ans and ans not in 'Yy':
        print 'Setup aborted.'
        exit()
    for d in old_directories:
        if exists(d):
            cmd = 'rm -rf '+d
            print cmd
            system(cmd)

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
    exit(1)

try:
    import scipy
except:
    print '[ERROR] failed to import scipy'
    exit(1)

try:
    import matplotlib
except:
    print '[ERROR] failed to import matplotlib'
    exit(1)

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
            exit(1)

    cmd = 'tar -xzf '+tarfile
    print cmd
    system(cmd)


## download other db files
# switching to dropbox as the default since some users networks don't like the port 7007 address
address = 'https://www.dropbox.com/s/kivfp27gbz2m2st/tcrdist_extras_v2.tgz'
backup_address = 'http://xfiles.fhcrc.org:7007/bradley_p/pub/tcrdist_extras_v2.tgz'

tarfile = address.split('/')[-1]
assert tarfile == backup_address.split('/')[-1]

if not exists( tarfile ):
    print 'Downloading database files'
    download_web_file( address )

    if not exists( tarfile ):
        print '[ERROR] download database files failed, trying a backup location'
        download_web_file( backup_address )

        if not exists( tarfile ):
            print '[ERROR] download database files failed'
            exit(1)

## md5sum check
lines = popen('md5sum '+tarfile).readlines()
if lines and len(lines[0].split()) == 2:
    # rhino1 tcr-dist$ md5sum tcrdist_extras_v2.tgz
    # 2705f3a79152cd0382aa6c5d4a81ad0b  tcrdist_extras_v2.tgz
    checksum = lines[0].split()[0]
    expected_checksum = '2705f3a79152cd0382aa6c5d4a81ad0b'
    if checksum == expected_checksum:
        print "\n[SUCCESS] md5sum checksum for tarfile matches expected, phew!\n"
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
        exit(1)

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
