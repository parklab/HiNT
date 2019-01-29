import os,sys
import numpy as np
from scipy.sparse import coo_matrix
from HiNT.straw import *

def get_chromInfo(chromlf):
	chroms = []
	infos = {}
	inf = open(chromlf)
	for line in inf:
		line = line.strip().split('\t')
		infos[line[0]] = int(line[1])
		chroms.append(line[0])
	return chroms,infos

def dumpMatrix(chrom1, chrom2, resolution, hicfile, chromInfo, outputMatrixFile):
    chrom1length = chromInfo[chrom1]
    chrom2length = chromInfo[chrom2]
    binsize = resolution*1000
    binnumber1 = np.divide(chrom1length,binsize) + 1
    binnumber2 = np.divide(chrom2length,binsize) + 1
    chr1 = chrom1.lstrip('chr')
    chr2 = chrom2.lstrip('chr')
    result = straw('KR', hicfile, str(chr1), str(chr2), 'BP', binsize)
    row = np.divide(result[0],binsize)
    col = np.divide(result[1],binsize)
    data = result[2]
    #print max(row),max(col),binnumber1,binnumber2
    res = coo_matrix((data, (row, col)), shape=(binnumber1, binnumber2)).toarray()
    np.savetxt(outputMatrixFile, res, fmt='%.5f', delimiter='\t')

def hicToMatrix(hicfile, resolution, chromlf, outputdir,name):
    MatrixInfo = {}
    chroms,chromInfo = get_chromInfo(chromlf)
    for i in range(len(chroms)-3):
        chrom1length = chromInfo[chroms[i]]
        for j in range(i+1,len(chroms)-2):
            chrom1 = chroms[i]
            chrom2 = chroms[j]
            outmatrixdir = os.path.join(outputdir,"InterMap_matrix")
            if not os.path.isdir(outmatrixdir):
                os.mkdir(outmatrixdir)
            outputname = os.path.join(outmatrixdir,name + '_%skb_%s_%s_InterMap_matrix.txt'%(resolution,chrom1,chrom2))
            dumpMatrix(chrom1, chrom2, resolution, hicfile, chromInfo,outputname)
            MatrixInfo[chrom1+'_'+chrom2] = outputname
    return MatrixInfo
