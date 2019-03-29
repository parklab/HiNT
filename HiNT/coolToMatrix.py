import os,sys
import numpy as np
import cooler

def getBins(coolfile):
	binsInfo = {}
	chroms = coolfile.chroms()["name"][:]
	for chrom in chroms:
		idxarray = np.where(coolfile.bins()["chrom"][:] == chrom)
		chromstart = idxarray[0][0]
		chromend = idxarray[0][-1]
		binsInfo[chrom] = [chromstart,chromend]
		#print chrom, (chromend - chromstart + 1)
	return binsInfo

def dumpMatrix(resolution,coolfile,binsInfo,chrom1,chrom2,outdir,name):

	chrom1start,chrom1end = binsInfo[chrom1]
	chrom2start,chrom2end = binsInfo[chrom2]
	matrixName = '_'.join([name, str(resolution)+'kb',chrom1,chrom2,"InterMap_matrix.txt"])
	matrixfile = os.path.join(outdir,matrixName)
	matrix = coolfile.matrix(balance=True)[chrom1start:(chrom1end + 1), chrom2start:(chrom2end + 1)]
	np.savetxt(matrixfile,matrix,fmt='%.5f',delimiter='\t')

	return matrixfile

def coolToMatrix(matrixFile,resolution,outdir,name):
	MatrixInfo = {}
	coolfile = cooler.Cooler(matrixFile)
	binsInfo = getBins(coolfile)
	rankedChroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
	outdir = os.path.join(outdir,"InterMap_matrix")
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	for i in range(0,len(rankedChroms)-2):
		for j in range((i+1), len(rankedChroms)-1):
			chrom1 = rankedChroms[i]
			chrom2 = rankedChroms[j]
			#print chrom1,chrom2
			matrixfile = dumpMatrix(resolution,coolfile,binsInfo,chrom1,chrom2,outdir,name)
			MatrixInfo[chrom1+'_'+chrom2] = matrixfile
	return MatrixInfo
