import os,sys
import cooler
import numpy as np
import pandas as pds
import scipy.sparse as sps
#from multiprocessing import Pool
#import argparse

def getBins(coolfile):
	binsInfo = {}
	chroms = coolfile.chroms()["name"][:]
	for chrom in chroms:
		idxarray = np.where(coolfile.bins()["chrom"][:] == chrom)
		chromstart = idxarray[0][0]
		chromend = idxarray[0][-1]
		binsInfo[chrom] = [chromstart,chromend]
	#print binsInfo
	return binsInfo

def getSumPerChrom(binsInfo,chroms,chrom1,coolfile,outputname):
	binstart1,binend1 = binsInfo[chrom1]
	chromRowsums = np.zeros((binend1-binstart1+1,1))
	for j in range(len(chroms)):
		chrom2 = chroms[j]
		print(chrom1,chrom2)
		binstart2,binend2 = binsInfo[chrom2]
		matrix = coolfile.matrix(balance=False)[binstart1:(binend1+1), binstart2:(binend2+1)]
		rowSum = np.nansum(matrix, axis=1)
		rowSum = np.reshape(rowSum,(binend1-binstart1+1,1))
		chromRowsums += rowSum
	writeGenomeRowSums(coolfile,chromRowsums,chrom1,outputname,name)

def writeGenomeRowSums(coolfile,transRowSum,chrom,outputname,name):
	print("Writing rowsums of %s!"%chrom)
	bins = coolfile.bins()[:]
	allbins = bins[bins['chrom']== chrom].ix[:,0:3]
	genomeRowSum = np.reshape(transRowSum,(len(allbins),1))
	dfrowsums = pds.DataFrame(genomeRowSum,columns=[name],index=allbins.index)
	dfrowsums.to_csv(outputname,sep="\t")

def calRowsums(params):
	coolfile,binsInfo,targetchrom,chunksize,name,outputname = params
	s,e = binsInfo[targetchrom]
	length = e - s + 1
	RowSums = np.zeros((length,1))
	steps = int(length/chunksize)
	#print(steps)
	rowsums = np.array([])
	if steps > 0:
		for j in range(0,steps-1):
			#print j
			start = s + j*chunksize
			end = s + (j+1)*chunksize
			#print start,end
			rowSum = getSumPerChunk(coolfile,start,end)
			rowsums = np.append(rowsums,rowSum)
		start = s + chunksize*(steps-1)
		if start > e:
			pass
		else:
			end = e
			rowSum = getSumPerChunk(coolfile,start,end+1)
			rowsums = np.append(rowsums,rowSum)
		rowsums = np.reshape(rowsums,(len(rowsums),1))
	else:
		rowsums = getSumPerChunk(coolfile,s,e+1)
		rowsums = np.reshape(rowsums,(len(rowsums),1))
	genomeRowSum = RowSums + rowsums
	#print(np.shape(genomeRowSum))
	writeGenomeRowSums(coolfile,genomeRowSum,targetchrom,outputname,name)
	outInfo = targetchrom + '\t' + outputname +'\n'
	return outInfo

def getallChromsRowSums(coolpath,name,outputdir,resolution):
	coolfile = cooler.Cooler(coolpath)
	binsInfo = getBins(coolfile)
	chroms = binsInfo.keys()
	rowSumFilesInfo = {}
	chroms = list(binsInfo.keys())
	excludechroms = ['chrY','chrM','Y','M']
	for excludechrom in excludechroms:
		if excludechrom in chroms:
			chroms.remove(excludechrom)
	print(chroms)
	for i,targetchrom in enumerate(chroms):
		outputname = os.path.join(outputdir,name + '_%s_%skb_GenomeRowSums.txt'%(targetchrom,str(resolution)))
		rowSumFilesInfo[targetchrom] = outputname
		getSumPerChrom(binsInfo,chroms,targetchrom,coolfile,outputname)

	#print results
	return rowSumFilesInfo
