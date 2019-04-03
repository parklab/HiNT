import numpy as np
from HiNT.straw import *
import os,sys

def get_chromInfo(chromlf):
	chroms = []
	infos = {}
	inf = open(chromlf)
	for line in inf:
		line = line.strip().split('\t')
		infos[line[0]] = int(line[1])
		chroms.append(line[0])
	return chroms,infos

def getSumPerChrom(i, j, hicfile, binsize, chroms, chromInfo, sumInfo):
    chrom1 = chroms[i] #the primary chromosome
    chrom2 = chroms[j] #the supplementary chromosome
    chr1 = chrom1.lstrip('chr')
    chr2 = chrom2.lstrip('chr')
    result = straw('NONE', hicfile, str(chr1), str(chr2), 'BP', binsize)
    #x is the coordinate for i, and y is the coordinate for j
    if i < j:
        xs = np.divide(result[0],binsize)
        ys = np.divide(result[1],binsize)
    else:
        xs = np.divide(result[1],binsize)
        ys = np.divide(result[0],binsize)
    values = np.array(result[2])
    chrom1length = chromInfo[chrom1]

    for n in sumInfo:
        idx = np.where(xs == n)
        nvalues = values[idx]
        nsum = np.sum(nvalues)
        sumInfo[n] += nsum
    return sumInfo

def writeGenomeRowSums(sumInfo,outputname,name,baseIdx):
    outf = open(outputname,'w')
    header=['',name]
    outf.write('\t'.join(header) + '\n')
    binsIdx = list(sumInfo.keys())
    binsIdx.sort()
    for idx in binsIdx:
        newidx = idx + baseIdx
        res = [str(newidx), str(sumInfo[idx])]
        outf.write('\t'.join(res) + '\n')
    outf.close()
    lastIdx = binsIdx[-1] + baseIdx
    return lastIdx

def getGenomeRowSums(resolution, hicfile, chromlf, outputdir,name):
    rowSumFilesInfo = {}
    chroms,chromInfo = get_chromInfo(chromlf)
    binsize = resolution * 1000
    baseIdx = 0
    for i in range(len(chroms)-2):
        sumInfo = {}
        chrom1length = chromInfo[chroms[i]]
        binnumber = int(chrom1length/binsize) + 1
        for n in range(binnumber):
            sumInfo[n] = 0
        for j in range(len(chroms)-2):
            sumInfo = getSumPerChrom(i, j, hicfile, binsize, chroms, chromInfo, sumInfo)
        outputname = os.path.join(outputdir,name + '_%s_%skb_GenomeRowSums.txt'%(chroms[i],str(resolution)))
        writeGenomeRowSums(sumInfo,outputname,name,baseIdx)
        rowSumFilesInfo[chroms[i]] = outputname
    return rowSumFilesInfo
