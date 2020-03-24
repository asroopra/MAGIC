#!/usr/bin/env python3

#####################################################################################################
# Version notes
# 1-tailed KS test using Hodges eqn
# Uses 5 percentile enrichment
# Uses c= uq/ub for coefficient
# lollipop markers for 5th percentile genes
# uses 'Target' rather than 'Driver' terminology
# saves Factor data for all non-triaged Factors
# uses .mtx files for matrix

# need to:
# fix zero errors for lists with no hits
# replace windows, osx file handler modules with speicifed divider
# make magicPathFinder_NEW windows compatible

#####################################################################################################


import os
import sys
from   sys import platform
import time
from   operator import itemgetter
import math
import xlwt
import random
import numpy as np
import matplotlib.pyplot as plt
from   matplotlib import gridspec as gridspec

#####################################################################################################
def test_platform():

	if platform =='darwin':
		pl = 'darwin'
		
	else:
		pl = 'win'

	return pl
#####################################################################################################

def Universal_FileHandler(operating_system, inputPath):
	
	if operating_system == "darwin":
		separator = "/"
		
	if operating_system == "win":
		separator = "\\"
		
		
	timestamp = str(int(time.time()))[-4:]

	home = inputPath[:inPath.rfind(separator)+1]

	folderPaths         = dict()
	folderPaths['home'] = home
	
	# open file and get list names (header elements from 1 to end)
	if inputPath[-3:].upper() == 'TXT' or inputPath[-3:].upper() == 'GMX':
		listNames  = open(inputPath).read().splitlines()[0].split('\t')[1:]
		
	if inputPath[-3:].upper() == 'CSV':
		listNames  = open(inputPath).read().splitlines()[0].split(',')[1:]
		
	# check all names are unique to avoid 'file already exists error'
	names = set(listNames)
	
	if len(listNames) > len(names):
		print ('''		
***********************************************
***       List names must be unique:        ***
***********************************************
		''')
		
		for name in names:
			c = listNames.count(name)
			if c > 1:
				print ('"%s" appears %d times.' %(name,c))

		print ('''
***********************************************
***  Please check list names and try again  ***
***********************************************\n\n
		''')
		
		sys.exit()

	# make subfolder paths
	for name in listNames:
		
		folderPaths[name]                          = folderPaths['home'] + name + '_' + timestamp + separator
		folderPaths[name + '_CDF_Folder']          = folderPaths[name] + 'CDFs' + separator
		folderPaths[name + '_NonSigCdfFolder']     = folderPaths[name + '_CDF_Folder'] + 'FDR_above_10_percent' + separator
		
		folderPaths[name + '_auxFolder']           = folderPaths[name] + 'Auxiliary_Files' + separator
		folderPaths[name + '_targetDataFolder']    = folderPaths[name + '_auxFolder'] + 'Target_Data' + separator
		folderPaths[name + '_NonSigFactorFolder']  = folderPaths[name + '_targetDataFolder'] + 'FDR_above_10_percent' + separator


		os.mkdir(folderPaths[name])
		os.mkdir(folderPaths[name + '_auxFolder'])
		os.mkdir(folderPaths[name + '_targetDataFolder'])
		os.mkdir(folderPaths[name + '_NonSigFactorFolder'])
		os.mkdir(folderPaths[name + '_CDF_Folder'])
		os.mkdir(folderPaths[name + '_NonSigCdfFolder'])

	return folderPaths
#####################################################################################################

def magicPathFinder(operating_system):

	# find root folder where *.py exists.
	# generate path to Matrices folder which should be in root folder
	# ifMatrices folder not in root folder, make one and ask for .mtxs matrix file
	
	fileNotFoundMessage = '''
****************************************************************************
NO MATRIX FILES FOUND.  PLEASE DRAG MAGIC MATRICES INTO \"MATRICES\" FOLDER.
****************************************************************************\n'''
	
	appPath = sys.argv[0] # list of files in cmd line starting with program itself.
	
	if operating_system  == 'darwin':
		separator = "/"

	if operating_system  == 'win':
		separator = "\\"	
	
	dirPath = appPath[:appPath.rfind(separator)+1]
		
	# check magic folder exists.  Make one if not
	magicFolderPath = dirPath + 'Matrices'
	
	if os.path.isdir(magicFolderPath) == False:
		os.mkdir(magicFolderPath)
		
		print (fileNotFoundMessage)
		
		sys.exit()
	
	# if folder exists, get list of folder contents
	content = os.listdir(magicFolderPath)		
	
	# keep only .mtx files.  ignore windows-generated system files that start with '._'
	matrices = [c for c in content if c[-3:].upper() =='MTX' and c[:2] != '._']
	
	# if more than 1 matrix file, get files with .mtx extension
	if len(matrices) > 1:	
		print ('\n- Enter the number corresponding to the matrix to be used:')

		m_count = 0
		for m in matrices:
			matrix_name = m[:-4]
			print ('[%s]...%s' %(str(m_count + 1), matrix_name))
			m_count += 1
		print (' \n')
		
		# get matrix to use
		while True:
			try:
				to_use = int(input()) - 1
			except:
				print (' ENTER A NUMBER BETWEEN 1 AND %s' %(len(matrices)))
				continue
				
			if to_use <0 or to_use > len(matrices) - 1:
				print (' ENTER A NUMBER BETWEEN 1 AND %s' %(len(matrices)))
				continue
				
			else:
				break

		magicMatrixPath = magicFolderPath + separator +  matrices[to_use]
	
	if len(matrices) == 1:
		magicMatrixPath = magicFolderPath + separator +  matrices[0]		
		
	if len(matrices) == 0:
		print(fileNotFoundMessage)
		sys.exit()
	
	return magicMatrixPath

#####################################################################################################
def magicExtractor(path):

	print ('- Extracting Matrix Data.')
	data   = open(path).read().splitlines()
	matrix = [ln.split('\t') for ln in data]
	
	factors = matrix[0][1:]
	symbols = [row[0].strip("\'") for row in matrix[1:]]
	
	return factors, symbols

#####################################################################################################
def listsMaker(inputPath, matrixGenes, directories):

	print ('- Organizing Lists.')
	
	suffix = inputPath[-3:].upper()
	
	if suffix == 'TXT' or suffix == 'GMX':
		inputData  = [l.split('\t') for l in open(inputPath).read().splitlines()]
		
	if suffix == 'CSV':
		inputData  = [l.split(',')  for l in open(inputPath).read().splitlines()]
		
	listNames  = inputData[0][1:] # first column is all genes in analysis

	masterList = [row[0]
	              .upper()
	              .replace('\'','')
	              .replace('"','')
	              .replace(' ','')
	               for row in inputData if row[0] != '']

	if '!!' in masterList: platformGenes = list(matrixGenes)
			
	else:
		platformGenes = list(set(masterList) & set(matrixGenes))

	raw_lists = [[row[column]
	              .upper()
	              .replace('\'','')
	              .replace('"','')
	              .replace(' ','')
	               for row in inputData[1:]] for column in range(1,len(listNames)+1)]
	              
	lists     = [list(set(l) & set(platformGenes)) for l in raw_lists]
	discarded = [list(set(l) - set(platformGenes)) for l in raw_lists]

	# remove blanks and 'na' and 'NA' from lists
	for ls in lists:
		if  ''  in ls: ls.remove('')
		if 'NA' in ls: ls.remove('NA')

	for d in discarded:
		if  ''  in d: d.remove('')
		if 'NA' in d: d.remove('NA')
		
	# combine lists, make square and rotate for saving to file
	l         = len(platformGenes)
	sqru_lists= [row + ['']*(l-len(row)) for row in lists]

	combinedLs=[platformGenes] + sqru_lists
	rotated   =[[row[column] for row in combinedLs] for column in range(len(platformGenes))]

	# make accepted lists file path and open file
	acceptListsPath = directories['home'] + 'Accepted_lists.txt'
	f_acceptedLists = open(acceptListsPath, 'w')
	
	# write header
	f_acceptedLists.write('Platform_Genes\t' + '\t'.join(h for h in listNames) + '\n')
	
	# write lists
	[f_acceptedLists.write('\t'.join('"%s"'%(i) for i in row) + '\n') for row in rotated]
	f_acceptedLists.close()

# Make discarded list file paths for each input list and save in relevant Aux folder
	for name in listNames:
		discardPath = directories[name + '_auxFolder'] + 'Triaged_Genes.txt'

		f_discard   = open(discardPath,'w')
		[f_discard.write('\n'.join(g for g in discarded[listNames.index(name)]))]
		f_discard.close()

	return listNames, lists, platformGenes

#####################################################################################################
def subMatrixGenerator(mtx_identifier, mtx, queryList):
	# mtx_identifier = 'Path' or 'Matrix' depending on whether mtx is a path or the actual matrix
	#queryList:   accepted background list
	#mtx = matrix
	#pltfrm = platform
	
	# if path provided, generate mtx else use mtx provided.
	
	if mtx_identifier.upper() == 'PATH':
		m_block = open(mtx).read().splitlines()
		mtx = [ln.split('\t') for ln in m_block]

	pltfrm  = [mtx[0]] # header

	[pltfrm.append([row[0].strip("\'")] + list(map(float,row[1:]))) for row in mtx if row[0].strip("\'") in queryList]
	
	return pltfrm
	
#####################################################################################################	
def writeSubMatrix(sub_mtx, path):
	# open file
	f_file = open(path,'w')
	
	# write header
	f_file.write('\t'.join(i for i in sub_mtx[0]) + '\n')
	
	# write matrix
	[f_file.write('\'%s\'\t'%(ln[0]) + '\t'.join(str(e) for e in ln[1:]) + '\n') for ln in sub_mtx[1:]]
	f_file.close()


#####################################################################################################
def Hodges_approximation(n1, n2, d):
	# took this from scipy.stats.ks_2samp v1.3.0 source module
	
	# Product n1*n2 is large.  Use Smirnov's asymptoptic formula
	# we have around 10K background genes
	
	# n1 = size of background
	# n2 = size of query
	# d = D statistics
	
	# Use Hodges' suggested approximation Eqn 5.3
	# np.exp(x) give e^x
	
	m = max(n1, n2)
	n = min(n1, n2)
	z = np.sqrt(m*n/(m+n)) * d

	prob = np.exp(-2 * z**2 - 2 * z * (m + 2*n)/np.sqrt(m*n*(m+n))/3.0)

	return prob    
            
#####################################################################################################
def Kolmogorov_Smirnov_Test(qB, mB):
	# qB = query data block (columns = TFs, rows = genes)
	# mB = matrix data block
	
	# return dict of:
	# KS statistics (U,p)
	# matrix and query histograms for each factor
	# arg of D (bin number)
	# ChIP value at D
	
	statistics = dict() # key = stats type, vals = lists of stats

	# get number of TFs
	nTFs = len(qB[0])
	
	q = np.array(qB)
	m = np.array(mB)
	
	# collect KS statistics for each TF as a list of lists [[D1,p1] , [D2,p2]...]

	
	# Determine if query distribution is to the right of matrix distribution
	#   -calculate distributions of chip signals for query list per factor
	#   -generate cumulatives
	#   -get differences, calculate maximum (max_d) and minimum difference (min_d)
	#   -if query cumulative to right of platform cumulative (|Dmax|>|Dmin|), assign a
	#    negative value to D statistic (make it invalid)
	
	d_args  = []
	d_chips = []
	mHistograms = []
	qHistograms = []
	ks_ls = []
	for z in range(nTFs):
		# mHist = platform histogram. mHist[0] = values, mHist[1] = bin edges (rhs)
		# i want 1 bin per value to get highest resolution cdf.  bins = len(set(m[:,z]))
                           
		mHistogram  = np.histogram(m[:,z],
		                               bins = len(set(m[:,z])),
		                               density = True)
		mHist = mHistogram[0]
		mBins = mHistogram[1]                       
		
		# qHist = query histogram   
                          
		qHistogram = np.histogram(q[:,z],
		                            bins = len(set(m[:,z])),
		                            density = True,
		                            range   = (min(m[:,z]), max(m[:,z])))  		          
		qHist = qHistogram[0]
		qBins = qHistogram[1]
		
		mHistograms.append(mHistogram) 
		qHistograms.append(qHistogram)
		
		# calculate normalized cumulatives of sample numbers

		m_cumsum = np.cumsum(mHist)
		norm_m_cumsum = m_cumsum/m_cumsum[-1]
		
		q_cumsum = np.cumsum(qHist)
		norm_q_cumsum = q_cumsum/q_cumsum[-1]
		
		# get difference in cumulatives and maxima and minima
		d     = norm_m_cumsum - norm_q_cumsum
		max_d = max(d)
		min_d = min(d)

		# get RHS most index of max_d and min_d
		max_d_pos = np.where(d == max_d)[0][-1]
		min_d_pos = np.where(d == min_d)[0][-1]

		# get chip values at max or min position
		max_d_chip = qBins[max_d_pos]
		min_d_chip = qBins[min_d_pos]
		
		# calculate KS p using Hodge equation.
		p = Hodges_approximation(len(m[:,z]), len(q[:,z]), max_d)
		ks_ls.append([max_d,p/2.]) # divide p by 2 for 1 tailed result

		# assign sign to D statistic
		#   (+ve if q_cumsum to right of m_cumsum AND more than 5 of q are targets else -ve)
			
		# sort query list and get value of 5th best chip
		q_sorted = np.sort(q[:,z])
		q_best = q_sorted[-5] 
		
		# filter for Factors with right shifted curves and more than 5 targets
		#      assign min_d to traiged factors
		if abs(min_d)> max_d and max_d_chip > q_best:

			ks_ls[z][0] = min_d	
			d_args.append(min_d_pos)
			d_chips.append(min_d_chip)
			
		else:
			d_args.append(max_d_pos)
			d_chips.append(max_d_chip)
			
	# return dict of : KS outputs as list
	#                  d args
	#                  chip at d
	
	statistics['KS']       = ks_ls
	statistics['d_args']   = d_args
	statistics['d_chips']  = d_chips
	statistics['plt_hist'] = mHistograms
	statistics['q_hist']   = qHistograms
	
	return statistics
#####################################################################################################

def getCoefficients(qB, mB):
	# qB = query data block
	# mB = matrix data block
	
	# make arrays where each column is sorted (independently of other columns)
	
	q = np.sort(qB, axis=0)
	m = np.sort(mB, axis=0)
	
	# get PERCentile of values in q and that number of top values in m
	# qMax = list of n values in 5th percentile of query block
	# mMax = list of  top n values in matrix block
	
	perc = 0.05
	qMax = q[int(len(q)*(1-perc)):]   
	mMax = m[len(m) - len(qMax):]
	
	# calculate means over the top values
	qU = np.mean(qMax, axis = 0)
	mU = np.mean(mMax, axis = 0)
	
	# replace zeros in mU to prevent div by zero error
	mU = [u if u > 0 else 0.000001 for u in mU  ]
	
	c = qU/mU
	
	return qU, mU, c, qMax, mMax
	
#####################################################################################################
	
def benjamini_Hoschberg_corrector(data_array): # nested array with last element of row being raw p value
	
	# to account for some samples with same p value, need to find number of ranks
	pvals = list(set([row[-1] for row in data_array])) # set of all p values
	pvals.sort()
	
	e_pvals = (e for e in enumerate(pvals)) # index is rank
	
	#make dictionary: key = p value, element = rank
	ranks = dict()
	for v in e_pvals:
		ranks[v[1]] = v[0]
	
	# sort by p value
	data_array.sort(key = itemgetter(-1))
	
	for row in data_array:
		adjusted_p = row[-1]*len(data_array)/(ranks[row[-1]] + 1)
		
		if adjusted_p > 1: adjusted_p = 1
		
		row.append(adjusted_p)
		
	return data_array
	
#####################################################################################################
def format_results(arr, res_type):
	# arr = results array with column[0] being Desription:Factor ordered by Score
	# produce array with column[0] = TF, column[1] = description and best factor
	# res_type = Summary or Details
	
	# generate array with column[0] = TF, column[1] = descrtiption

	block = [[row[0].split(":")[1], row[0].split(":")[0]] + row[1:] for row in arr]

	
	#collect rows with best scoring condition per factor
	trimmed_block     = []
	collected_factors = []
	
	for r in block:

		if r[0] not in collected_factors:
			collected_factors.append(r[0])
			
			trimmed_block.append(r)
	
	if res_type == "Summary":
		return trimmed_block
	
	if res_type == "Details":
		return block

#####################################################################################################

def writeExcelFile(path, header, res_type, _array):
	# header = list of elements to go in header
	#_array = results array
	# res_type = tab label
	# sort by Score (last element in list)
	
	_arr =[r for r in _array if len(r)>7]
	
	_arr.sort(key = itemgetter(-1), reverse = True)
	
	bk = xlwt.Workbook()
	
	header_Style = xlwt.easyxf('font: name Calibri, color-index black, height 240, bold on, italic on')
	FDR_10_Style = xlwt.easyxf('font: name Calibri, height 240, color-index red, bold on')
	FDR_NS_Style = xlwt.easyxf('font: name Calibri, height 240, color-index black')
	
	sheet = bk.add_sheet('res_type')
	
	# cycle though header and write to file.
	[sheet.write(0, c, header[c], header_Style) for c in range(len(header))]
	
	#write key
	sheet.write(0, len(header) + 0, 'FDR<10%', FDR_10_Style)
	sheet.write(0, len(header) + 1, 'FDR>10%', FDR_NS_Style)
	
	# write results to file
	for ln in _arr:
		
		# get rid of accession in annotation for summary
		if res_type == "Overview":
			ln[1] = ln[1].split("_")[0]  

		row   = _arr.index(ln)+1
		raw_p = ln[6]
		fdr   = ln[7]
		
		if fdr < 0.1:
			style = FDR_10_Style
		
		if fdr >= 0.1:
			style = FDR_NS_Style
		
		[sheet.write(row, col, ln[col], style) for col in range(len(ln))]
		
	bk.save(path)

#####################################################################################################

def writeResultsCSVFile(path, header_row, _array):
	#_array = results array
	# sort by Score (last element in list)
	_array.sort(key = itemgetter(-1), reverse = True)
	
	data = [[row[0].split(":")[1], row[0].split(":")[0]] + row[1:] for row in _array]
	
	f_csv = open(path,'w')
	
	# write header row
	f_csv.write(','.join(q for q in header_row) + '\n')
	
	# write data
	for row in data:
		if len(row) == 10:	f_csv.write(','.join(str(i) for i in row) + '\n')
		else:	            f_csv.write(','.join(str(i) for i in row) + ',NA,NA\n')
	
	f_csv.close()

#####################################################################################################

def makeTargetLists(subMtx,listName, allGenes, res, path):
	# subMtx    = sub Matrix for list
	# listName  = name of list
	# AllGenes  = all platform genes
	# res_arr   = results array
	# path      = target list path

	TFs       = subMtx[0]
	chips     = subMtx[1:]
	listGenes = list(np.array(subMtx)[1:,0])
	
	# make targets array and add 1st row (all expressed genes)
	# ultimately, array will be rotated 90 degrees so this becomes 1st column
	targets = [['All','na'] + allGenes]
	
	# add second row (genes in gene list)
	targets.append([listName,'na'] + listGenes) # 2d array of drivers for each TF

	for ln in res:
		
		f = ln[1] + ":" + ln[0] # condition/description : factor
		tf= ln[0] # transcription factor
		d = ln[1] # description
		c = ln[2] # critical chip value
		a = ln[7] # raw p
		b = ln[8] # FDR

		if a >= 0.05: continue
		if b >= 0.1: continue
		
		# get column for current TF
		col = TFs.index(f)
		
		# order matrix by chip signal in column
		chips.sort(key = itemgetter(col), reverse = True)
		
		# make list of targets (genes with ChIP >= mean ChIP for list)
		t = [l[0] for l in chips if float(l[col]) >= c]
		
		targets.append([tf,d] + t)
	
	# make targets array square
	maxL = len(targets[0])
	sqr  = [row + ['']*(len(targets[0])-len(row)) for row in targets]

	# rotate for saving to file
	rot  =[[row[column] for row in sqr] for column in range(maxL)]

	#save gene lists as *.gmx file
	f_targets = open(path,'w')
	
	# write TF header row and 'na' row
	[f_targets.write('\t'.join(i for i in row) + '\n') for row in rot[:1]]
	
	# write gene lists with quotes protection
	[f_targets.write('\t'.join('"%s"'%(g) for g in row) + '\n') for row in rot[1:]]

	f_targets.close()
#####################################################################################################

def saveTargetData(subMtx, res, fol):
	# subMtx = sub Matrix for list
	# res    = results array
	# fol    = path to Aux/Target_Gene_Data folder
	
	TFs       = subMtx[0]
	chips     = subMtx[1:]
	
	for ln in res:
		
		f = ln[1] + ":" + ln[0] # factor
		c = ln[2] # critical chip value
		a = ln[7] # raw p value
		b = ln[8] # FDR
		
		# get column for current TF
		col = TFs.index(f)
		
		# order matrix by chip signal in column
		chips.sort(key = itemgetter(col), reverse = True)
		
		# make list of targets and chip signals (genes with ChIP >= mean ChIP for list)
		gene_Val = [[l[0], l[col]] for l in chips if float(l[col]) >= c]
		
		if b < 0.1:	path = fol + '%s.CSV'%(ln[0])
		if b >=0.1:	path = fol + 'FDR_above_10_percent/' + '%s.CSV'%(ln[0])
		

		f_path = open(path,'w')
		f_path.write('Gene,ChIP\n')
		
		[f_path.write(','.join(str(i) for i in ln) + '\n') for ln in gene_Val]
		
		f_path.close()

#####################################################################################################
def save_Triaged_Factors(res, list_name, path):
	# res = all results array
	# list_name = name of list
	# path = file path
	
	triaged = []
	for row in res:
		if row[2] < 0: triaged.append(row[0])
	
	f_path = open(path,'w')
	for t in triaged:
		f_path.write(t + '\n')
	
	f_path.close()
	
#####################################################################################################
def drawCDFs(factors, platform_histograms, query_histograms, res_array, ls_name, query_blk, matrix_blk, query_tails, matrix_tails, dirs):
	# factors   = TFs in Matrix header
	# res_array = results array
	# query_blk = query datablock
	# query_tails  = 5th percentile query values block
	# matrix_tails = 5th percentile matrix values block
	
	max_diffs     = []
	rank          = 0
	draw_count    = 0
	non_sig_count = 0
	
	# get all FDRs to see how mant curves to draw (all sigs + 20 more)
	fdrs     = [float(r[8]) for r in res_array]
	num_sigs = len(np.where(np.array(fdrs) < 0.1)[0])
	
	CDFs_to_draw  = num_sigs + 20
	
	for row in res_array:
		rank += 1
		
		alpha  = row[7]
		fdr    = row[8]
		tf     = row[1] + ":" + row[0]
		crt_chp= row[2]
		
		if fdr < 0.1:
			filePath = dirs[ls_name + '_CDF_Folder'] + str(rank) + '_' + row[0] + '.pdf'
			draw_count += 1
			
		else:
			# limit number of non sig CDFs to 20
			non_sig_count += 1
			if non_sig_count >20:
				continue
								
			filePath = dirs[ls_name + '_NonSigCdfFolder'] + row[0] + '.pdf'
			draw_count += 1
			
			
		column = factors.index(tf)
		q_chips= np.array(query_blk)[:,column]
		p_hist = platform_histograms[column]
		q_hist = query_histograms[column]
		
		m_chips= np.array(matrix_blk)[:,column]
		
		q_tail_chips = np.array(query_tails)[:,column]
		m_tail_chips = np.array(matrix_tails)[:,column]
		
		# get matrix (expected) cumulative
		expCumulative = np.cumsum(p_hist[0])/sum(p_hist[0])
		exp_bins      = q_hist[1]

		# get query (observed) cumulative 
		obsCumulative = np.cumsum(q_hist[0])/sum(q_hist[0])
		obs_bins      = p_hist[1]

		sys.stdout.write(' Drawing %s [%s/%s CDFs]                \r' %(row[0], str(draw_count), str(CDFs_to_draw)))
		sys.stdout.flush()
		
		plt.figure(figsize = (6,4))
		
		# draw CDFs with grid
		plt.grid(linestyle = 'dashed')

		# set axis limits
		plt.xlim(left = exp_bins[-1]/-250, right = exp_bins[-1]+(exp_bins[-1]/250))
		plt.ylim(bottom = -0.1, top = 1)
		
		# draw CDF curves
		plt.plot(exp_bins[1:], expCumulative,
		         color = 'black',
		         linewidth = 2,
		         label = 'All Genes')
		         
		plt.plot(obs_bins[1:], obsCumulative,
		         color = 'red',
		         linewidth = 2,
		         label = ls_name)

		# draw line at maximal divergence
		plt.plot([crt_chp, crt_chp], [-0.1,1],
		         color = 'blue',
		         linewidth = 1 )
		
		# fill gene tick bar
		plt.fill_between([0, exp_bins[-1]], [0,0], [-0.1,-0.1],
		                facecolor = 'yellow',
		                alpha = 0.25)
		
		# draw horizontal line at y=0
		plt.axhline(0, color = 'black',
		            linewidth = 0.5)
		
		# draw query gene ticks
		plt.plot(q_chips, [-0.05 for q in q_chips],
		                   '|',
		                   color = 'red',
		                   markersize = 20)
		                   
		# draw query tail markers
		plt.plot(q_tail_chips,[0 for i in q_tail_chips],
		                      'o',
		                      color = 'red',
		                      markersize = 4)
		                   
		# draw matrix gene ticks
		plt.plot(m_chips, [-0.075 for m in m_chips],
		                   '|',
		                   color = 'black',
		                   markersize = 10)	

		# draw matrix tail markers
		plt.plot(m_tail_chips, [-0.05 for j in m_tail_chips],
		                   'o',
		                   color = 'black',
		                   markersize = 4)	

		# write axis labels
		plt.ylabel('Cumulative Gene Fraction',
		           fontname = 'Arial',
		           weight   = 'bold')
		           
		plt.xlabel('ChIP Signal',
		           fontname = 'Arial',
		           weight   = 'bold')
		
		# legend
		plt.legend(fontsize = 'small')
		
		# title
		plt.title(row[0])
		
		plt.savefig(filePath)
		plt.close()

		
#####################################################################################################
def drawScoreSummary(ls_name, res, dirs):

	path = dirs[ls_name] + ls_name + '_Score Summary.pdf'
	
	# get tf and score for TFs with FDR<0.10 to plot
	tfs    = [row[0]  for row in res if row[-2] < 0.10]
	scores = [row[-1] for row in res if row[-2] < 0.10]
	
	if len(tfs) > 0:
		
		n = range(len(tfs), 0, -1)
	
		if len(tfs) < 30:
			lineWidth = 10
			fontSize  = 11
			gs = gridspec.GridSpec(40,1, left = 0.2)            # set num of subplots in fig to 40
		
		if len(tfs) >=30:
			lineWidth = 10 * 35/len(tfs)            # set width to decrease with increasing num of TFs
			fontSize  = 9 * 35/len(tfs)             # set font to decrease with increasing num of TFs
			gs = gridspec.GridSpec(len(tfs) + 4,1, left = 0.2 )  # set num of subplots in fig to 4 more than TFs
		
		ax1 = plt.subplot(gs[:len(tfs)*2,0])        # define sublot height as 2x more that TFs to plot
		                                            # by combining sublots that I set above
		                                          
		ax1.hlines(n, [0], scores,		        
				      lw=lineWidth,
				      color = 'red')
		                  
		ax1.set_yticks(n)
		ax1.set_yticklabels(tfs, fontsize = fontSize)
	
		ax1.set_ylim(bottom = 0.5, top = len(tfs) + 0.5)
		ax1.set_xlim(left = 0)
	
		ax1.tick_params(axis  = 'y',
						left  = False,
						right = False)
	
		ax1.tick_params(axis  = 'x',
						top   = False)
					
		ax1.set_xlabel('Score')
		ax1.set_title ('Enriched Factors')
		
	
		plt.savefig(path)
		plt.close()
	
	
#####################################################################################################

OpSys = test_platform()

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
\nPlace the *.txt or *.csv file into a folder where you wish the results to appear.
Then drag the file into this terminal or enter it's file path here.\n
 ''')
 
if OpSys == 'darwin':
	while True:
		inPath = input().replace('\"','').strip().replace('\\','')
		if os.path.isfile(inPath) != False: break
		
if OpSys != 'darwin':
	while True:
		inPath = input().replace('\"','').strip()
		if os.path.isfile(inPath) != False: break
		
folders = Universal_FileHandler(OpSys, inPath)

magicPath = magicPathFinder(OpSys)


TFs, mtxGenes  = magicExtractor(magicPath)

geneLists   = listsMaker(inPath, mtxGenes, folders)
lsNames     = geneLists[0]
lsGenes     = geneLists[1]
pltfrmGenes = geneLists[2]


# generate Platform Matrix
print ('- Generating Platform Matrix.')

pltfrmMtxPath = folders['home'] + 'Platform_Matrix.txt'

pltfrmMtx = subMatrixGenerator('path',magicPath, pltfrmGenes)

pltfrmBlk = [list(map(float, r[1:])) for r in pltfrmMtx[1:]]

for name in lsNames:

	print ('\n- Generating Sub Matrix For Genes In List: %s.' %(name))
	# Generate Sub Matrices for each list
	geneList   = lsGenes[lsNames.index(name)]
	subMtxPath = folders[name + '_auxFolder'] + name + '_Sub_Matrix.txt'
	subMatrix  = subMatrixGenerator('platform matrix',pltfrmMtx, geneList)
	
	# Generate data block
	dataBlock = [list(map(float,row[1:])) for row in subMatrix[1:]]
	
	# dataBlock is vals for query genes per tf
	# pltfrmBlk = matrix values block
	
	# perform KS tests and get chip vals at D
	
	raw_stats = Kolmogorov_Smirnov_Test(dataBlock, pltfrmBlk)
	ks        = raw_stats['KS'] 
	d_list    = [m[0] for m in ks]   # list of D statistics
	p_list    = [m[1] for m in ks]   # 1 tailed p values
	
	crit_Chp  = raw_stats['d_chips'] # critical chip values (chip at max d)
	
	platform_histograms = raw_stats['plt_hist']
	query_histograms    = raw_stats['q_hist']
	
	# calculate ratio of Obs and Exp max chips
	coeff_data   = getCoefficients(dataBlock, pltfrmBlk)
	
	obs_tail_mean = coeff_data[0]
	exp_tail_mean = coeff_data[1]
	coefficients  = coeff_data[2]
	obs_tails     = coeff_data[3]
	exp_tails     = coeff_data[4]
	
	# make a results block:
	resultsBlock = [TFs, crit_Chp, d_list, obs_tail_mean, exp_tail_mean, coefficients, p_list]
	
	# rotate resultsBlock
	# call this 'all_results' because all subsequent list additions point to this list
	# and it gets modified.  In the end, it has scores etc for relevant rows!
	all_results   = [[row[column] for row in resultsBlock] for column in range(len(TFs))]
	
	# filter for TFs with d>0
	filtered_results = [row for row in all_results if row[2] > 0]
	
	results = benjamini_Hoschberg_corrector(filtered_results)
	
	# append Score to results.  if corrected p=0, corrected p = 1e-16
	for r in results:
		ratio = r[5]
		pBH   = r[7]
		
		if pBH > 0: score = -1*ratio*math.log10(pBH)
		else:       score = ratio*16
		
		# if pBH = 1, score becomes -0 and need to be set at 0.
		if score == 0: score = 0 
		
		r.append(score)
	
	results.sort(key = itemgetter(-1), reverse = True)
	
	top_results = format_results(results, "Summary")
	
	xls_header_row = [
		  'Factor',
		  'Description',
		  'Critical ChIP',
		  'Obs Tail Mean',
		  'Exp Tail Mean',
		  'Tail Enrichment',
		  'Raw p',
		  'Corrected p',
		  'Score'
			  ]
	
	txt_header_row = [
		  'Factor',
		  'Description',
		  'Critical ChIP',
		  'D statistic',
		  'Obs Tail Mean',
		  'Exp Tail Mean',
		  'Tail Enrichment',
		  'Raw p',
		  'Corrected p',
		  'Score'
			  ]
	
	# remove D statistics (row[2]) from top_results for overview xls file
	overview = [row[:3] + row[4:] for row in top_results]
	
	# write overview to xls file
	file_path = folders[name] + '%s_Summary.xls' %(name)
	writeExcelFile(file_path, xls_header_row, "Overview", overview)
	
	# remove D statistics (row[2]) for from results for detailed xls file
	details = [row[:3] + row[4:] for row in all_results]
	
	detailed_results = format_results(details, "Details")
	
	# write details to xls file
	file_path = folders[name] + '%s_Details.xls' %(name)
	writeExcelFile(file_path, xls_header_row, "Details", detailed_results)


	# write all results to txt file in Auxiliary folder for distribution graphs
	txtResPath = folders[name + '_auxFolder'] + name + '_raw_results.CSV'
	writeResultsCSVFile( txtResPath, txt_header_row, all_results)
	
	#Draw summary graph
	drawScoreSummary(name, top_results, folders)	

	print ('-\tWriting Target Data')
	
	# write Target summary
	target_path = folders[name] + '%s_Targets.gmx' %(name)
	makeTargetLists(subMatrix, name, pltfrmGenes, top_results, target_path)

	targetDataFolderPath = folders[name + '_targetDataFolder']
	saveTargetData(subMatrix, top_results, targetDataFolderPath)
	
	# save list matrix
	print ('-\tWriting Sub Matrix For %s:' %(name))
	writeSubMatrix(subMatrix, subMtxPath)

	#draw CDFs
	print ('-\tSaving CDF Plots For %s:' %(name))
	drawCDFs(TFs, platform_histograms, query_histograms, top_results, name, dataBlock, pltfrmBlk, obs_tails, exp_tails, folders)
	print ('                                                                ')

	# save triaged Factors list
	f_path = folders[name + '_auxFolder'] + 'Triaged_Factors.txt'
	save_Triaged_Factors(all_results, name, f_path)

# save platform matrix
print ('- Writing Platform Matrix.')
writeSubMatrix(pltfrmMtx, pltfrmMtxPath)