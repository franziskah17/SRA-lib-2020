Documentation for "sralib.py"

This program uses .fastq files as an input and outputs two .txt files with tab-separated tables of NC contents and Di-NC contents with means and standard deviations (one table for all spots individually and one for the total library). It also outputs a boxplot and a histogram for NC means and stdev.

_____________________________________________________________________________
Technical requirements: 
(program written with these specifications, on Linux Mint 20)

-) Python 3.8

-) rpy2 package version 3.3.5
https://rpy2.github.io/doc/v3.0.x/html/index.html

-) R 3.6.3

-) R ggplot2 package version 3.3.2
https://cran.r-project.org/web/packages/ggplot2/index.html

-) SRA Toolkit (latest version, as of June 29, 2020 this is version 2.10.8 ). 
https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
______________________________________________________________________________
Installation Help: 

R was installed with 
$ sudo apt-get install r-base
(r-base-dev should be installed with it)

ggplot2 was installed by downloading ggplot2_3.3.2.tar.gz from CRAN (link above) and the following command
$ sudo R CMD INSTALL -l /usr/lib/R/site-library /home/Your/Path/ggplot2_3.3.2.tar.gz
this path was given, as this is where R searched for packages when using rpy2, but is not the default for installing. If you get an error message about packages not being installed, copy the path stated in the message into here.

rpy was installed with
$ pip install rpy2

Add the sralib folder to a path you can execute from (e.g. ~/.local/bin) or call program from current dir. Check valid paths with 
$ less .bashrc

make the python and bash scripts executable with 
$ chmod u+x sralib.py
$ chmod u+x sralib_pipeline.sh
_____________________________________________________________________________
HOW TO OBTAIN YOUR .fastq FILE from NCBI's SRA (short read archive):

The newest verison of SRA Toolkit is required. 


Search for your Project of Interest on https://www.ncbi.nlm.nih.gov/sra/
eg. a search like this:
(((shape[All Fields] AND map[All Fields]) OR (structure[All Fields] AND probing[All Fields])) AND biomol rna[Properties]) AND "Homo sapiens"[Organism]
Pick a project by ticking the box, then chose "send to">"file">"Format: Accession List" at the top.
Save the accession list. 

This contains accession numbers of all runs done in this project. Each run contains many spots (~1mio), as can be seen when opening the project on SRA. 


Obtain the .fastq files from each run by using the following parameters: 
$ fastq-dump
--split-spot (separates forward and backward reads into the same file)
--readids (creates different read IDs for fw and bw. Spot name will be "accession.spot.readid")
--dumpbase (Formats seq using base space, no colorspace)
--skip-technical (dumps only biological reads, no techincal reads)
--clip (apply left and right clips of tags)
--read-filter pass (should filter out pure N reads. program does doublecheck that itself too)
--outdir <path> (default is current wd)
SRR_ID_Input (no flag for input required)

The result should be one .fastq file per run (per entry in your accession list)
It contains fw and bw reads and quality scores. 
you can save them compressed, using the parameter --gzip, but will need to unzip it before use.
Nice explanation of parameters here: https://edwards.sdsu.edu/research/fastq-dump/

_____________________________________________________________________________
				RUNNING sralib.py
_____________________________________________________________________________

CURRENTLY ONLY RUNNING ONE RUN AT A TIME IS SUPPORTED BY SRA LIB!

If you want to run all runs of an SRA Project at the same time (e.g. the whole list of accession numbers you have downloaded from NCBI), please use the provided 
"sralib_pipeline.sh"

It will automatically download all .fastq files from your input accession list 
(so no pre-use of SRA Toolkit required)
and pipe them into the sralib.py program. 

To change parameters of either program, you can do so by editing sralib_pipeline.sh.
-------
DEFAULT of the bash script for fastq-dump is as above, 
for sralib the input is each .fastq from fastq-dump, -o without custom filename and -g (GRAPHICS ON!!).
-------


Help message from sralib.py:

$ ./myprogramplot.py -h
usage: myprogramplot.py [-h] -i INPUTFILE [-o [OUTFILENAME]] [-g] [-s SPOTLIST]

DESCRIPTION: This program takes .fastq input from the NCBI SRA and gives you general information about all spots (NC
contents, means, stdev)

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE, --input INPUTFILE
                        (REQUIRED) specify your input file (fastq)
  -o [OUTFILENAME], --output [OUTFILENAME]
                        specify outfilename with -o filename . Extension is automatical. If only -o is used, outfile
                        is created with INPUTFILE.txt as default name automatically. If argument is not used, output
                        is printed in stdout (commandline)
  -g, --graphics        add if you want graphics to be created. Default is False to reduce runtime.
  -s SPOTLIST, --spotlist SPOTLIST
                        If you want statistics for a subset of spots only, enter a .txt file with line separated
                        spot names.

First specify your input file! Then add additional arguments. A new folder for your project will be created in your
wd, based on the inputfile name. To download and run multiple runs of a Project at once, please use
sralib_pipeline.sh


