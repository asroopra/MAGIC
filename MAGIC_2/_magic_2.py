#!/usr/bin/env python

import pandas as pd
import numpy as np
import scipy.stats as sc

import math

import plotly.express as px

import sys
from pathlib import Path

import file_handler as fh
from magic_calculations import get_d_arg, calculate_padj
from magic_graphics import histplot, summary_fig
from magic_results import write_target_data, write_gmx_file




######################################################################################################

def get_cmd_line_args(params):
	
	'''
	#default params.
			
	params = {
	"script_path": script_path,
	"ls_path": "",
	"matrix": Path(f"{script_path}/Matrices/1Kb_gene.pkl.gz"),
	"padj_cutoff": 0.1,
	"min_targets": 10,
	"keep_gtfs": "Y",
	"keep_zeros": "Y"
	
	}
	'''
	
	cmd = sys.argv
	script_path = cmd.pop(0)

	help_text = """
 -f:\tPath to lists file [REQUIRED].
 -p:\tMaximum acceptable adjusted p value for figure generation. Default=0.25  [OPTIONAL].
 -t:\tMinimum number of target genes for Factor to be called. Default=10 [OPTIONAL].
 -m:\tPath to matrix file. Default = /Matrices/1Kb_gene.pkl.gz [OPTIONAL].
 -g:\tKeep General Transcription Factors and RNApol in analysis (y/n).  Default = y [OPTIONAL].
 -z:\tKeep zero chip values (y/n). Default = y [OPTIONAL]  
	"""
	
	# cycle through cmd line, get keys and check vals
	bad_options = []
	
	while len(cmd)>0:
		_opt = cmd.pop(0).lower()
		
		if _opt == "-f":
			params["ls_path"] = cmd.pop(0)
			continue
		
		if _opt == "-p":
			params["padj_cutoff"] = float(cmd.pop(0))
			continue
					
		if _opt == "-t":
			params["min_targets"] = int(cmd.pop(0))			
			continue
		
		if _opt == "-g":
			params["keep_gtfs"] = cmd.pop(0).upper()
			continue
			
		if _opt == "-z":
			params["keep_zeros"] = cmd.pop(0).upper()		
			continue
			
		if _opt in ["-h","--h", "-help","--help"]:
			print (help_text)
			exit()

		bad_options.append(_opt)
	
	if len(bad_options) == 0:
		return params
	
	else:
		print ("\n FOLLOWING TERMS NOT FOUND.")
		print (bad_options)
		
		print ("\n ALLOWED OPTIONS:\n")
		print (help_text)
		exit()
			
	return params

######################################################################################################
def get_ls_path():

	print ('''
	This program will generate a series of output files contained in a folder.
	Make a tab delimited text file (*.txt), Comma Separated Values (*.csv) 
	or *.gmx file with the following format:\n\n
	All	List_1	List_2	....List_n
	A1BG	SNAP25	CHRGB	....SYP1
	A1CF	UNC13	CHRGA	....ESR1
	...	...	...	...
	\nThe first row contains UNIQUE labels for your gene lists.
	The first column contains all genes on your platform i.e the Background list.
	The 'tighter' this list, the more robust the statistics will be.
	Subsequent columns contain gene lists of interest.
	\nDrag the *.txt or *.csv file into this terminal or enter it's file path here.\n
	 ''')	

	while True:

		list_path = input().strip()

		p = Path(list_path)

		if p.exists:
			break

		print (" - ENTER VALID FILE PATH - ")
	
	return p
######################################################################################################
def get_gtf_option():
	
	while True:
		keep_gtfs = input("\n - Analyze GTFs? [y/n]:")
		if keep_gtfs[0].upper() in ["Y","N"]:
			break
		else:
			print (" - ENTER Y or N - ")
	
	return (keep_gtfs[0].upper())
	

######################################################################################################
def keep_zero_option():
	
	while True:
		keep_zero = input("\n - Keep Zeros? [y/n]:")
		if keep_zero[0].upper() in ["Y","N"]:
			break
		else:
			print (" - ENTER Y or N - ")
	
	return (keep_zero[0].upper())
	

######################################################################################################
def make_ls_df(ls_path):

	ls_path = str(ls_path)
	
	if ls_path.split(".")[-1].upper() in ["TXT","TSV"]:
		ls_df = pd.read_csv(ls_path, sep="\t")

	if ls_path.split(".")[-1].upper() == "CSV":
		ls_df = pd.read_csv(ls_path)

	if ls_path.split(".")[-1].upper() in ["XLS","XLSX"]:
		ls_df = pd.read_excel(ls_path)
	
	# make 1st column name 'All'
	cols = ls_df.columns.to_list()
	ls_df.columns = ["All"] + cols[1:]
	
	return ls_df
######################################################################################################	

script_path = sys.argv[0]

#default params.

gtfs = ["TBP","POL","GTF","TAF"]
	
params = {
"script_path": script_path,
"ls_path": "",
"matrix": Path(f"{script_path}/Matrices/1Kb_gene.pkl.gz"),
"keep_gtfs":"Y",
"keep_zeros":"N",
"padj_cutoff": 0.1,
"min_targets": 10
}

if __name__== "__main__":
	
	# get params from command line if options provided
	if len(sys.argv)>1:
		params = get_cmd_line_args(params)

	# get gtfs option and list file path from get_ls_path module if
	# no cmd line options provided
	if len(sys.argv)==1:
		params["ls_path"] = get_ls_path()
		params["keep_gtfs"] = get_gtf_option()
		params["keep_zeros"] = keep_zero_option()
		
	
# get working directory and params
_dir = Path(params["script_path"]).parent

mtx_path    = Path(f"{_dir}/Matrices/1Kb_gene.pkl.gz")
keep_gtfs   = params["keep_gtfs"]
padj_cutoff = params["padj_cutoff"]
min_targets = params["min_targets"]
ls_path     = params["ls_path"]
keep_zeros  = params["keep_zeros"]

# get matrix and filter out GTFs if required
mtx = pd.read_pickle(mtx_path)
mtx.iloc[:,1:] = mtx.iloc[:,1:].astype("float64") # do this for plotly compatibilty

if keep_gtfs == "N":
	keep_cols = ["GENE"]
	for e in mtx.columns[1:]:
		tf = e.split(":")[-1]
		if tf[:3] not in gtfs:
			keep_cols.append(e)

	mtx = mtx[keep_cols]


# format matrix genes to upper case
mtx_ls = [g.strip("\'") for g in mtx["GENE"].to_list()]
mtx["GENE"] = mtx_ls


# make results folder
MAGIC_output_path = fh.make_magic_folder(ls_path)


# make lists df
ls_df = make_ls_df(ls_path)

# get list names
lists = ls_df.columns[1:].to_list()

# make master list from lists file
master_ls = ls_df.iloc[:,0].str.upper().to_list()
if master_ls[0] =="!!":
	master_ls = mtx_ls

	
# trim matrix to master list and modify master list accordingly
mtx = mtx[mtx["GENE"].isin(master_ls)]
master_ls = mtx["GENE"].to_list()

# get experiments from matrix
expts = list(mtx.columns)[1:]	
	
list_counter = 1
for gene_ls in lists:

	# make results folder, subfolders and file paths for this gene list
	paths = fh.make_subfolders(output_folder = MAGIC_output_path,
	                                    analysis_name = gene_ls)
	'''
	paths["target_fol"]  folder for target data
	paths["distributions"] folder for distribution figs
	paths["summary"]     excel file
	paths["details"]     excel file
	paths["Factor_gmx"]  TSV file
	paths["summary_fig"] pdf file
	'''
	
	# get query genes for current list
	q_ls = ls_df[gene_ls].dropna().str.upper().to_list()

	# make results df to populate later.  rows = expts in matrix, cols = stats
	res_df = pd.DataFrame(index=expts, columns = ["D","argD","p","Score"])

	expt_counter = 1

	# make dict to hold expt name and lists of non-zero chip vals for master and query list
	# Will be used to generate CDF figures.
	chip_vals = dict()# key = expt, val = [master df, query df]
	
	# make dict to hold target genes and chip values as df for each experiment in matrix
	# will be used to generate Target files (genes and chip values for positive expts)
	# key = expt, val = df of genes and chip vals for current expt
	target_data = dict()

	for e in expts:

		sys.stdout.write(f" - List {list_counter}/{len(lists)}. Expt {expt_counter} of {len(expts)}           \r")
		sys.stdout.flush()
		
		# print (f" - List {list_counter}/{len(lists)}. Expt {expt_counter} of {len(expts)}.  {e}")
		# make df of genes and chip vals for current expt in matrix
		df = mtx[["GENE",e]]

		# and eliminate zero binders if required
		if keep_zeros == "N":
			df = df[df[e]>0]
		
		# get master chip values
		m_vals = df[e].to_list()
		
		# make df with just query genes
		qdf = df[df["GENE"].isin(q_ls)]
		
		# ... and get query chip vals
		q_vals = qdf[e].to_list()
		
		# store master and query df in dict. key = expt
		chip_vals[e] = [df, qdf]

		# perform 1 tailed kolmogorov-Smirnov test
		try:
			D,p = sc.ks_2samp(m_vals, q_vals, alternative="greater", method="asymp")
		except:
			continue

		# get arg_d and arg_d_index (arg_d = chip value, arg_d_index is posn in array)
		arg_d, arg_d_index = get_d_arg(m_vals=m_vals, q_vals=q_vals)
		
		# calculate scores
		score = -1*math.log(p)*D	
		
		# get target genes
		# trim df to genes with chip>arg_d
		high_val_df = df[df[e]>arg_d]

		# sort df by chip value ranks (0=poorest, 1=best)
		high_val_df = high_val_df.sort_values(by=e, ascending=True)
		
		# trim for genes in query list
		target_df = high_val_df[high_val_df["GENE"].isin(q_ls)]
		
		# only keep TF if more than min_targets target genes.
		if len(target_df) < min_targets:
			continue
		
		# reverse order so df goes from best to worst chip value
		target_df = target_df.sort_values(by=e, ascending=False)

		# add target df to targets dict.  key = expt, val = target df
		target_data[e] = target_df
	
		res_df.loc[e] = [D,arg_d,p,score]
	
		expt_counter +=1

	# sort by p value 
	res_df = res_df.sort_values(by="p")

	# get TFs from expt names and add to results df
	expts = list(res_df.index) # order of expts changed from mtx due to sorting by p
	res_df["TF"] =[e.split(":")[-1] for e in expts]
	
	# save res_df
	res_df.to_excel(paths["details"])
	
	# add expt (currently index in form expt:TF) as column
	res_df["Experiment"] = expts

	# get best expt for each TF for summary file
	summary_df = res_df.groupby("TF", as_index=False).first().sort_values(by="p")
	summary_df = summary_df.dropna(axis=0, how="any")
	
	# reorder columns
	summary_df = summary_df[["TF","Experiment","D","argD","p","Score"]]
	
	# calculate padj
	summary_df = calculate_padj(summary_df)
	
	# save summary file
	summary_df.to_excel(paths["summary"], index=False)
	
	# draw figs and generate target data is any TFs pass padj_cutoff
	if summary_df["padj"].min() <=padj_cutoff:	
	
		print ("\n - Writing target data.")
		# write target data (gene and chip value for each TF)
		write_target_data(summary_df = summary_df,
										 target_data = target_data,
										 padj_cutoff = padj_cutoff,
										 folder = paths["target_fol"])
									 
		print (" - Writing gmx file.")
		write_gmx_file(summary_df = summary_df,
									 master_ls = master_ls,
									 target_data = target_data,
									 padj_cutoff = padj_cutoff,
									 f_path = paths["Factor_gmx"])
	
		print (" - Generating summary figure.\n")
											 
		print (" - Generating figures.\n")
		summary_fig(summary_df=summary_df,
		            padj_cutoff=padj_cutoff,
		            fig_path=paths["summary_figure"])
		            
		histplot(summary_df=summary_df,
											chip_vals=chip_vals,
											padj_cutoff=padj_cutoff,
											fig_folder=paths["distributions"])
		
	list_counter += 1
	
