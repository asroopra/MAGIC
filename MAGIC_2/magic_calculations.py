#!/usr/bin/env python

import numpy as np
import pandas as pd

def calculate_padj(df):

	df = df.sort_values(by="p")
	df["rank"] = np.arange(1,len(df)+1)
	df["padj"] = df["p"]*len(df)/df["rank"]
	df["padj"] = np.where(df["padj"]>=1., 1., df["padj"])
	
	df.drop(labels="rank", axis=1, inplace=True)
	
	df = df.sort_values(by="padj")
	
	return df
############################################################################

def get_d_arg(m_vals, q_vals):

	# mHist = master list histogram
	mHistogram  = np.histogram(m_vals,
                               bins = len(set(m_vals)),
                               density = True)
	mHist = mHistogram[0]
	norm_mHist = mHist/np.max(mHist)
	mBins = mHistogram[1]                       
	
	# qHist = query histogram   			
	qHistogram = np.histogram(q_vals,
                              bins    = len(set(m_vals)),
                              density = True,
                              range   = (min(m_vals), max(m_vals)))  		          
	qHist = qHistogram[0]
	norm_qHist = qHist/np.max(qHist)
	qBins = qHistogram[1]


	m_cumsum = np.cumsum(mHist)
	norm_m_cumsum = m_cumsum/m_cumsum[-1]
	
	q_cumsum = np.cumsum(qHist)
	norm_q_cumsum = q_cumsum/q_cumsum[-1]
	
	d = norm_m_cumsum - norm_q_cumsum
	
	# get argument of d as posn
	arg_d_pos = np.argmax(d)+1
	# convert to chip value
	arg_d = qBins[arg_d_pos]

	
	return arg_d, arg_d_pos