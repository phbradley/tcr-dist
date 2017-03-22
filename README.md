# TCRdist pipeline, version 0.0.1 beta

Thanks for your interest in TCRdist! This software is very much a
work in progress, so we welcome any/all feedback and will work to
add features and fix bugs as quickly as we can.

If you would like to contribute to the project on github, feel free
to contact us.

Phil Bradley (pbradley@fredhutch.org)
Jeremy Chase Crawford (Jeremy.Crawford@stjude.org)

See LICENSE.txt for details on the (MIT) license.

---
# REQUIREMENTS

## System:

 - So far we have mainly tested on LINUX systems. I was also able to get things to work on my mac laptop after updating scipy, numpy, and sklearn and installing inkscape.

## Python version:
 - We have tested with 2.7.11

## Package dependencies

 - I use anaconda to keep these updated on my system; I think pip will also work, although it sounds like mixing and matching anaconda and pip can be problematic.

 - scipy:
   https://www.scipy.org/install.html
   tested with version 0.16.0
				
 - scikit-learn:
   aka "sklearn" for KernelPCA, adjusted_mutual_info_score
   http://scikit-learn.org/stable/install.html
   tested with version 0.17
				
 - matplotlib:
   tested with version 1.4.3
				
 - numpy:
   tested with version 1.10.1


## External command line tools

 - convert: (or rsvg-convert or inkscape)
   from Imagemagick is used to convert svg files to png files if you have an alternative you can modify the function "convert_svg_to_png" in basic.py

 - wget: (or curl)
   for downloading database and other files
   If you have something else that works similarly on your system, feel free to modify setup.py or contact me to add that as an option.
				
## NCBI Blast
 - The setup script will automatically download a compatible version into the ./external/ directory.

---
# INSTALLATION

1. Go to the tcr-dist/ directory (main code directory)

2. Run the command:
```
python setup.py
```

3) Cross your fingers.

There are some potentially useful comments at the top of setup.py

---
# USAGE

For an overview of what the analysis is supposed to do, consult the publication listed below in the "CITING" section.

The basic workflow starting from a sequence file would be to run:
```
python run_basic_analysis.py --organism <organism> --pair_seqs_file <filename>
```
where `<organism>` is either mouse or human, and `<filename>` is the name of a .tsv (tab-separated values) file with the following fields:
```
id	epitope	subject	a_nucseq	b_nucseq	a_quals	b_quals
```
`a_quals` and `b_quals` are '.' separated lists of the quality scores for the corresponding nucleotide sequences (a_nucseq and b_nucseq).

Try running:
```
python run_basic_analysis.py -h
```
for some help text.

You can also run individual steps in the analysis from the command line. Running a script with the -h option will print a very basic help message listing the command line arguments. We are currently working hard to flesh out the help messages -- our apologies for the lack of clarity. Some scripts depend on the output of previous steps (for example compute_distances.py generates distance matrices which are used by downstream programs). The source code for run_basic_analysis.py gives an example of an appropriate order for calling the scripts.

---
# TESTING

(After running setup.py)

In the tcr-dist/ directory you can run the command:
```
python create_and_run_testing_datasets.py
```
which will create a directory called testing/, make two small dataset files there, and run the analysis pipeline on them.

Running `setup.py` should have created a directory called `tcr-dist/testing_ref/` which should contain examples of the `*_web/index.html` output generated by the pipeline. So if you compare the results files that correspond between
```
testing/*web/index.html
```
and
```
testing_ref/*web/index.html
```
they should look pretty similar (at least all the deterministic components of
the analysis). You can also search for the text "missing" in those html
results to see if any of the results files are missing from the `testing/`
version but present in the `testing_ref/` version.

You can also re-run the analysis in the original paper by typing:
```
python rerun_paper_analysis.py
```
which may take a couple of hours. Add the `--from_pair_seqs` option to restart from nucleotide sequences. And/or the `--multicore` option to let the pipeline spawn multiple, independent processes.

---
# CITING

---
# THANKS

Part of this analysis uses parameters derived from nextgen data in publicly released studies. We are grateful to the authors of those studies for making their data available. See the `README_db.txt` file in `./db/` (after running setup.py)

The code uses the command line parsing toolkit `blargs`. See the license and info in `external/blargs/`

The tables in the .html results can be sorted thanks to "tablesorter". See the license and info in `external/tablesorter/`

Sequence parsing relies on the BLAST suite, see info in `external/blast-2.2.16/`
