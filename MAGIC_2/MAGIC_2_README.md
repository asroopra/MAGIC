# MAGIC

The manuscript accompanying this script is:
MAGIC: A tool for predicting transcription factors and cofactors driving gene sets using ENCODE data.
Roopra, A. 2020
PLoS Computational Biology


Implementation of MAGIC_2.

Drag _magic_2.py into a terminal or type in the file path and hit return.

Alternatively, _magic_2 accepts command line inputs.  If none are provided, a file path to the lists file will be requested by the script upon running.
Command line options are:

-f	File path to lists file [REQUIRED].
-p:	Maximum acceptable adjusted p value for figure generation. Default=0.25  [OPTIONAL].
-t:	Minimum number of target genes for Factor to be called. Default=10 [OPTIONAL].
-m:	Path to matrix file. Default = /Matrices/1Kb_gene.pkl.gz [OPTIONAL].
-g:	Keep General Transcription Factors and RNApol in analysis (y/n).  Default = y [OPTIONAL].  
-z:	Keep zero ChiP values (y/n). Default = n [OPTIONAL]

NOTE on z option:  Keeping zero Chip signals will recreate MAGIC outputs very similar to MAGIC outlined in Roopra,2020.  Eliminating zero chip values will place more weight on Factors with higher ChIP values and de-emphasize those Factors that may regulate genes via moderate/lower binding affinities.  This tends to produce outputs more like - but not the same as - tests based on Fisher Exact Tests with gene lists and gene sets such as EnrichR.

A tab delimited text file is requested by MAGIC (lists file).  The first column contains a list of all genes expressed in the experiment (background list).  Any number of other columns are then added containing query lists.  For example, a query list may be all genes that go up under some criterion and another may be all genes that go down.  The first row is the header and must have unique names for each column.  These names will be used to generate file and folder names so do not use special characters, spaces, period etc in column names.

MAGIC analyzes each query list and generate a series of output files and sub-directories. A series of directories are generated named after each query list in the lists file.  Each directory contains sub-directories and files for the stand-alone analysis of that query list. A Query_List_Details.xlsx file contains statistical information for all non-triaged ENCODE experiments. Data reported in Query_List_Details.xls are:

Experiment: The ENCODE experiment.  format = cell-line_experiment-ID:Factor e.g K562_ENCFF558VPP:REST
D:					Kolmogorov Smirnov statistic
argD:				Argument (ChipValue) of D
p:					Kolmogorov Smirnov p value
Score:			-log(p)xD
TF:					Factor name 

A summary file contains the best scoring Experiment for each Factor.  Each Factor appears once in this file.  The format is:
TF:					Factor Name
Experiment: The best scoring ENCODE experiment for a Factor.  format = cell-line_experiment-ID:Factor e.g K562_ENCFF558VPP:REST
D:					Kolmogorov Smirnov statistic
argD:				Argument (ChipValue) of D
p:					Kolmogorov Smirnov p value
Score:			-log(p)xD
padj:				Benjamini-Hochberg corrected p value 

Query_list_Drivers.gmx is a tab delimited file in the GMX format utilized by GSEA.  The first column contains the background list of genes.  Subsequent columns contain all target genes of each Factor.  Target genes are defined as those genes in the query list whose ChIP signal is greater than the argD i.e. the ChIP at which there is the maximal difference between the population and query cumulative.

An Auxiliary_Files directory contains 2 sub directories:
'Distributions' contains html files showing the PDFs and CDFs for ChIP values in the query and master list for each gene in the ENCODE experiment
'Target_Data_' contains comma separated text files for each Factor with padj less than user defined or default cutoff.  Each file contains a list of target genes for that Factor and associated MAGIC Matrix ChIP value.

-----
To iterate through multiple lists files, for example lists files for each cell type in a single cell seq analysis:

import subprocess

# place following in loop to cycle through list file paths
	# Arguments to pass to MAGIC_2.py
	args = ['python', 'MAGIC_2.py', '-f', 'path to lists file', other arguments if needed]

	# Run the script with the specified arguments
	subprocess.run(args, check=True)

