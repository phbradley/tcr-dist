####################################################################################################
TCRdist pipeline, version 0.0.1 beta
####################################################################################################

Thanks for your interest in TCRdist! This software is very much a
work in progress, so we welcome any/all feedback and will work to
add features and fix bugs as quickly as we can.

If you want to directly add to and modify the repository on github, feel free
to contact me (Phil); we would welcome the input.

Phil Bradley (pbradley@fredhutch.org)

See LICENSE.txt for details on the (MIT) license.


###########################
REQUIREMENTS
###########################

python:
				I have tested with 2.7.11

## PYTHON MODULES #########

I use anaconda to keep these updated on my system; I think pip will also work, although
it sounds like mixing and matching anaconda and pip can be problematic.

scipy:
				https://www.scipy.org/install.html
				tested with version 0.16.0
				
scikit-learn:
				aka "sklearn" for KernelPCA, adjusted_mutual_info_score
				http://scikit-learn.org/stable/install.html
				testing with version 0.17
				
matplotlib:
				tested with version 1.4.3
				
numpy:
				tested with version 1.10.1


## COMMAND LINE TOOLS #####

convert:
				from Imagemagick is used to convert svg files to png files
				if you have an alternative you can modify the function
				"convert_png_to_svg" in basic.py

wget:
				for downloading database and other files
				If you have something else that works similarly on your system,
				feel free to modify setup.py or contact me to add that as an option.
				
## NCBI Blast #############
The setup script will automatically download a compatible version into the ./external/ directory.


###########################
INSTALLATION
###########################

1) Go to the TCRdist/ directory (main code directory)

2) run the command:

	 python setup.py

3) cross your fingers.

There are some potentially useful comments at the top of setup.py

###########################
USAGE
###########################

For an overview of what the analysis is supposed to do, consult the publication
listed below in the "CITING" section.

The basic workflow starting from a sequence file would be to run

python run_basic_analysis.py --organism <organism> --pair_seqs_file <filename>

where <organism> is either mouse or human, and <filename> is the name of a .tsv
(tab-separated values) file with the following fields:



###########################
TESTING
###########################


###########################
CITING
###########################


###########################
THANKS
###########################

Part of this analysis uses parameters derived from nextgen data in publicly
released studies. We are grateful to the authors of those studies for making their
data available. See the README_db.txt file in ./db/ (after running setup.py)

The code uses the command line parsing toolkit "blargs". See the license and info in external/blargs/

The tables in the .html results can be sorted thanks to "tablesorter". See the license and info in
external/tablesorter/

Sequence parsing relies on the BLAST suite, see info in external/blast-2.2.16/

