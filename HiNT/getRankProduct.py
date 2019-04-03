import os,sys
import numpy as np
from scipy.stats import rankdata
from operator import itemgetter

def readBackgroundMatrix(backgroundMatrixDir):
    print(backgroundMatrixDir)
    backgroundMatrixInfo = {}
    matrices = os.listdir(backgroundMatrixDir)
    for matrix in matrices:
        nameInfo = matrix.split('_')
        chrompair = nameInfo[2] + '_' + nameInfo[3]
        matrixf = os.path.join(backgroundMatrixDir,matrix)
        backgroundMatrixInfo[chrompair] = matrixf

    return backgroundMatrixInfo

def gini(x):
    # (Warning: This is a concise implementation, but it is O(n**2)
    # in time and memory, where n = len(x).  *Don't* pass in huge
    # samples!)

    # Mean absolute difference
    mad = np.nanmean(np.abs(np.subtract.outer(x, x)))
    # Relative mean absolute difference
    rmad = mad/np.nanmean(x)
    # Gini coefficient
    g = 0.5 * rmad
    return g

def getGini(mat1,mat2):
    matrix1 = np.genfromtxt(mat1,delimiter="\t")
    matrix2 = np.genfromtxt(mat2,delimiter="\t")
    matrix1[np.isfinite(matrix1)==0] = 0
    matrix2[np.isfinite(matrix2)==0] = 0
    rowsum1 = np.sum(matrix1,axis=1)
    rowsum2 = np.sum(matrix2,axis=1)
    colsum1 = np.sum(matrix1,axis=0)
    colsum2 = np.sum(matrix2,axis=0)
    ridx1 = np.where(rowsum1==0)
    cidx1 = np.where(colsum1==0)
    ridx2 = np.where(rowsum2==0)
    cidx2 = np.where(colsum2==0)
    ridx = np.union1d(ridx1[0], ridx2[0])
    cidx = np.union1d(cidx1[0], cidx2[0])

    temp1 = np.delete(matrix1,ridx,0)
    temp2 = np.delete(matrix2,ridx,0)
    selectedData1 = np.delete(temp1,cidx,1)
    selectedData2 = np.delete(temp2,cidx,1)

    average1 = np.mean(selectedData1)
    average2 = np.mean(selectedData2)
    tm1 = np.divide(selectedData1,average1)
    tm2 = np.divide(selectedData2,average2)
    division = np.divide(tm1,tm2)
    giniIndex = gini(np.asarray(division).reshape(-1))
    maximum = np.nanmax(np.asarray(division).reshape(-1))

    return giniIndex,maximum

def getRankProduct(matrix1MbInfo,background1MbInfo,outdir,name):
    rpout = os.path.join(outdir,name + '_chrompairs_rankProduct.txt')
    outf = open(rpout,'w')
    ginis = []
    maximums = []
    chrompairs = []
    for chrompair in matrix1MbInfo:
        #print chrompair
        matrix1 = matrix1MbInfo[chrompair]
        matrix2 = background1MbInfo[chrompair]
        giniIndex,maximum = getGini(matrix1,matrix2)
        chrompairs.append(chrompair)
        ginis.append(giniIndex)
        maximums.append(maximum)
    rankgini = len(ginis) - rankdata(ginis)
    rankmaximum = len(maximums) - rankdata(maximums)
    #print rankgini,rankmaximum
    rps = (np.divide(rankgini,len(ginis)*1.0))*(np.divide(rankmaximum,len(maximums)*1.0))
    result = np.stack((chrompairs,ginis,maximums,rps),axis=-1)
    sortedResult = sorted(result, key=itemgetter(-1))
    outf.write('\t'.join(['ChromPair',"GiniIndex","Maximum","RankProduct"]) + '\n')
    for res in sortedResult:
        chrompair, gini, maximum, rp = res
        newres = [chrompair, str(gini), str(maximum), str(rp)]
        outf.write('\t'.join(newres) + '\n')
    outf.close()
   
    return rpout

def getRPinfo(rpout):
    rpInfo = {}
    inf = open(rpout)
    for line in inf:
        if line.startswith('Chrompair'):
            pass
        else:
            line = line.strip().split('\t')
            chrompair, gini, maximum, rp = line
            rpInfo[chrompair] = rp
    return rpInfo
