# MAGIC

The manuscript accompanying this script is:
MAGIC: A tool for predicting transcription factors and cofactors driving gene sets using ENCODE data.
Roopra, A. 2020
PLoS Computational Biology



MAGIC requires Matrix files that can be downloaded from:
Download Matrices.zip and unpack.
Place unzipped folder in the same directory/folder as MAGIC.py prior to running MAGIC.py


Implementation of MAGIC.
A tab delimited text file is requested by MAGIC (lists file).  The first column contains a list of all genes expressed in the experiment (background list).  Any number of other columns are then added containing query lists.  For example, a query list may be all genes that go up under some criterion and another may be all genes that go down.  The first row is the header and must have unique names for each column.  MAGIC then requests which Matrix to use.

MAGIC analyzes each query list and generate a series of output files and sub-directories in the directory containing the original lists file. An ‘Accepted_Lists.txt’ file is generated which is the original input file filtered for genes in the Matrix.  A Platform_Matrix.txt file is the matrix file filtered for genes in the the background list.  A series of directories are generated named after each query list in the lists file.  Each directory contains sub-directories and files for the stand-alone analysis of that query list. A Query_List_Details.xls file contains statistical information for all non-triaged Factors. All Factors associated with Pcorr < 10% are highlighted in bold red. Data reported in Query_List_Details.xls are:

Factor			Name of Factor
Description		The ENCODE cell line, tissue or experiment description
Critical ChIP: 		argDsup i.e. ChIP value at dsup (This value is used to determine target genes in the list)
Obs Tail Mean: 	Average of the 95th percentile ChIP values (n values) in the query list.
Exp Tail Mean: 	Average of the top n ChIP values in the backgound.
Tail Enrichment: 	Ratio of the Obs and Exp Tail Means (r)
Raw P:  		Kolmogorov-Smirnov p value
Corrected P:		Benjamini Hochberg corrected p value (Pcorr)
Score: 			Score (-log(Pcorr)x r)

Where ENCODE contains multiple Chip-seq tracks for a Factor, the best scoring is also reported in a Summary.xls file with the same layout as above – each Factor appears once in this file.  A query_list_summary.pdf contains a bar graph of Factors and Scores with Pcorr <10%.

Query_list_Drivers.gmx is a tab delimited file in the GMX format utilized by GSEA.  The first column contains the background list of genes.  The second contains the query list.  Subsequent columns contain all target genes of each Factor.  Target genes are defined as those genes in the query list whose ChIP signal is greater than the Critical ChIP (argDsup) i.e. the ChIP at which there is the maximal difference between the population and query cumulative.

A sub-directory called ‘CDFs’ contains graphical displays of the analysis for all non-triaged Factors.  The naming format is ‘rank (integer)’_’factor (string)’.pdf (e.g. 1_NRSF.pdf; rank = 1, factor = NRSF) (S1 Fig) where ranking is determined by Score.  Two cumulative functions are displayed: the black curve is the fractional cumulative of all genes in the background list against ChIP values, red is the same for query genes.  A blue vertical line denotes the ChIP value at dsup i.e. argDsup.  Red ticks along the x-axis represent each gene in the query list and black ticks are all genes in the background.  Red ticks with circles (‘lollipops’)  are the n=0.05X best chiped genes.  Black lollipops are genes in the background list with the n highest ChIP values.

A second sub-directory called ‘Auxiliary_Files’ is populated with data behind the summary files.  The ‘query_list_raw_results.CSV’ file contains the same columns as the Query_list_Summary.xls file but has raw data for all factors including those that were triaged and not considered for further statistical analysis. It also contains the Kolmogorov-Smirnov D statistic for each factor.  D statistics with a negative sign denote D values for triaged factors; the negative sign is used by the algorithm for triage sorting.  

Query_list_Sub_Matrix.txt is the MAGIC Matrix filtered for genes in the query list.  

Triaged_Factors.txt is a list of factors that were not considered (triaged).

Triaged_Genes.txt contains all genes in the query that were not in the MAGIC Matrix and therefore eliminated from analysis.

A sub-directory named ‘Target_Data’ contains comma separated text files for each Factor with Pcorr < 10%.  Each file contains a list of target genes for that Factor and associated MAGIC Matrix ChIP value.

