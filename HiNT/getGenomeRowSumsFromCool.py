import os,sys
import cooler
import numpy as np
import pandas as pds
import scipy.sparse as sps
from multiprocessing import Pool
import argparse

def getBins(coolfile):
	binsInfo = {}
	chroms = coolfile.chroms()["name"][:]
	for chrom in chroms:
		idxarray = np.where(coolfile.bins()["chrom"][:] == chrom)
		chromstart = idxarray[0][0]
		chromend = idxarray[0][-1]
		binsInfo[chrom] = [chromstart,chromend]
	#print binsInfo
	return chroms,binsInfo

def getSumPerChunk(coolfile,start,end):
	length = end - start
	rowSums = np.zeros((length,1))
	#print np.shape(rowSums)
	binNum = np.shape(coolfile.bins())[0]
	size = 100
	steps = int(binNum/size)
	for s in range(steps):
		start2 = s*size
		end2 = (s+1)*size
		matrix = coolfile.matrix(balance=False)[start:end, start2:end2]
		rowSum = np.sum(matrix, axis=1)
		rowSum = np.reshape(rowSum, (len(rowSum),1))
		rowSums += rowSum
	if steps*size < binNum:
		lastmatrix = coolfile.matrix(balance=False)[start:end, steps*size:]
		rowSum = np.sum(lastmatrix, axis=1)
		rowSum = np.reshape(rowSum, (len(rowSum),1))
		rowSums += rowSum
	else:
		pass

	return rowSums

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

def getallChromsRowSums(coolpath,name,outputdir,resolution,threads):
	chunksize=int(1000/int(resolution))
	coolfile = cooler.Cooler(coolpath)
	chroms,binsInfo = getBins(coolfile)
	rowSumFilesInfo = {}

	allparamsInfo = []
	for i,targetchrom in enumerate(chroms):
		outputname = os.path.join(outputdir,name + '_%s_%skb_GenomeRowSums.txt'%(targetchrom,str(resolution)))
		allparamsInfo.append([coolfile,binsInfo,targetchrom,chunksize,name,outputname])
		rowSumFilesInfo[targetchrom] = outputname
	results = []
	p = Pool(threads)
	result = p.map_async(calRowsums, allparamsInfo, callback=results.append)
	p.close()
	p.join()
	#print results

	return rowSumFilesInfo
