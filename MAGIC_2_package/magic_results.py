#!/usr/bin/env python

import numpy as np

import pandas as pd
from pathlib import Path

def write_target_data(summary_df, target_data, padj_cutoff, folder):
	
	# target_data = dict of dfs.
	# key =expt, val = target_df with cols [["GENE","ChIP"]]
	# genes already sorted by chip value

	# make df with TFs<=padj_cutoff
	df = summary_df[summary_df["padj"]<=padj_cutoff]
	
	# list of expts
	expts = df["Experiment"].to_list()
	
	# for each expt, get factor, make file name and save as df
	for e in expts:
		f_name = e.split(":")[1] + ".csv"
		_path = str(Path(folder/f_name))
		
		# target_data[e] is a data frame
		_df = target_data[e].copy()
		_df.columns = ["Gene","ChIP"]
		
		_df.to_csv(_path, index = False)

####################################################################################
def 	write_gmx_file(summary_df, master_ls, target_data, padj_cutoff, f_path):
	# target data = dict. key = expt, val = gene/chip val df
	
	_df = summary_df[summary_df["padj"]<=padj_cutoff]
	
	_arr = [["All"] + master_ls]
	
	l = len(_arr[0])
	
	experiments = _df["Experiment"].to_list()
	
	for expt in experiments:
		tf = expt.split(":")[1]
		genes = target_data[expt]["GENE"].to_list()
		
		ln = [tf]+genes + [""]*(l - len([tf]+genes))
		
		_arr.append(ln)
	
	_arr = np.array(_arr)
	#gmx = pd.DataFrame(column = _arr[:,0])
	
	gmx =pd.DataFrame(_arr[:,1:])
	gmx = gmx.T
	
	gmx.columns = [_arr[:,0]]
	
	gmx.to_csv(f_path, sep = "\t", index=False)
	