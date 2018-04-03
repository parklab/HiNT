import os,sys
from HiNT.corelib import *
import numpy as np
import math

def getDivisionMatrix(mat1,mat2,fname):
	observedMatrix = np.genfromtxt(mat1,delimiter="\t")
	expectedMatrix = np.genfromtxt(mat2,delimiter="\t")
	expectedFinitevalues = expectedMatrix[np.isfinite(expectedMatrix)]
	observedFinitevalues = observedMatrix[np.isfinite(observedMatrix)]
	expected = expectedMatrix/np.mean(expectedFinitevalues)
	observed = observedMatrix/np.mean(observedFinitevalues)
	expected[np.where(expected==0)]=1
	with np.errstate(all='ignore'): divisionMatrix = observed/expected
	np.savetxt(fname,divisionMatrix,fmt="%.3f",delimiter='\t')
	return

def getAllRoughBreakpoints(matrix100kbInfo,background100kbInfo,rpInfo,outdir,name,cutoff,RscriptBPcaller):
	DivisionMatrixInfo = {}
	outputsubdir = os.path.join(outdir,"AdjustedMatrixFiles100kb")
	if not os.path.isdir(outputsubdir):
		os.mkdir(outputsubdir)

	for chrompair in matrix100kbInfo:
		if chrompair not in background100kbInfo:
			pass
		else:
			if float(rpInfo[chrompair]) <= float(cutoff):
				mat1 = matrix100kbInfo[chrompair]
				mat2 = background100kbInfo[chrompair]
				fname = os.path.join(outputsubdir,name + '_' + chrompair + '_DivisionMatrix.txt')
				'''
				getDivisionMatrix(mat1,mat2,fname)
				'''
				DivisionMatrixInfo[chrompair] = fname
			else:
				pass
	#print DivisionMatrixInfo
	bpoutputfile = os.path.join(outdir,name + '_roughBP_100kb.txt')
	command = "Rscript %s %s %s"%(RscriptBPcaller,outputsubdir,bpoutputfile)
	print command
	'''
	run_cmd(command)
	'''
	return bpoutputfile,DivisionMatrixInfo

################################################################################
#############     below are doing the breakpoints filtering   ##################

def getchromsize(chromlengthf,resolution):
	chromSizeInfo = {}
	inf = open(chromlengthf)
	for line in inf:
		line = line.strip().split('\t')
		chrom,size = line
		chromSizeInfo[chrom] = [int(size),(int(size)/int(resolution)+1)]
	return chromSizeInfo

def readBPs(breakpointsf,chromSizeInfo):
	bpInfo = {}
	inf = open(breakpointsf)
	for line in inf:
		if line.startswith('chrom'):
			pass
		else:
			line = line.strip().split('\t')
			chrom,bprow,maxrow,bpcol,maxcol = line
			chroms = chrom.split('_')
			chr1,chr2 = chroms
			size1 = chromSizeInfo[chr1][0]
			size2 = chromSizeInfo[chr2][0]
			if size1 > size2:
				bpInfo[chrom] = [(chr1,bpcol,maxcol),(chr2,bprow,maxrow)]
			else:
				bpInfo[chrom] = [(chr1,bprow,maxrow),(chr2,bpcol,maxcol)]
	return bpInfo

def mergeBPs(bps,matrix,axis,windowsize=5):
	bps = [int(i)-1 for i in bps] #python index starts with 0
	bps.sort()
	mergedbps = [str(bps[0])]
	for i in range(1,len(bps)):
		a = int(bps[i])
		b = int(mergedbps[-1])
		#print a,b
		#a is always larger than b
		if (a - b) <= windowsize:
			if axis == "row":
				suma = np.nansum(matrix[a,:])
				sumar = np.nansum(matrix[a+1,:])
				sumb = np.nansum(matrix[b,:])
				sumbr = np.nansum(matrix[b+1,:])
			elif axis == 'column':
				suma = np.nansum(matrix[:,a])
				sumar = np.nansum(matrix[:,a+1])
				sumb = np.nansum(matrix[:,b])
				sumbr = np.nansum(matrix[:,b+1])
			else:
				print "give the correct axis, row or column"
			if suma > sumar and sumb < sumbr: # __|--|__
				fca = suma/sumar
				fcb = sumbr/sumb
				if suma > sumbr and fca > fcb:
					print "merge",np.nansum(matrix[a-1,:]), suma, sumar, np.nansum(matrix[b-1,:]),sumb, sumbr, fca, fcb
					mergedbps.append(str(a))
					mergedbps.remove(str(b))
				elif suma < sumbr and fca > fcb:
					mergedbps.append(str(a))
				else:
					continue
		else:
			mergedbps.append(str(a))
	return mergedbps

def relativefold(a,b):
	rfold = (float(a) - float(b))/float(b)
	return rfold

def filtering(bps,matrixfile,chromSizeInfo,windowsize):
	matrix = np.loadtxt(matrixfile)
	rowsum = np.nansum(matrix,axis=1)
	colsum = np.nansum(matrix,axis=0)
	chrom1Info,chrom2Info = bps
	chr1,bp1,max1 = chrom1Info
	chr2,bp2,max2 = chrom2Info
	bp1 = bp1.split(',')
	bp2 = bp2.split(',')
	bp1s = bp1
	bp2s = bp2
	bp1s = list(set(bp1s))
	bp2s = list(set(bp2s))
	#print bp1s,bp2s
	mergedbp1s = mergeBPs(bp1s,matrix,'row')
	#print mergedbp1s
	if str(int(max1)-1) not in mergedbp1s:
		mergedbp1s = mergedbp1s + [str(int(max1)-1)]
	mergedbp2s = mergeBPs(bp2s,matrix,'column')
	#print mergedbp2s

	if str(int(max2)-1) not in mergedbp2s:
		mergedbp2s = mergedbp2s + [str(int(max2)-1)]
	#print mergedbp1s,mergedbp2s

	validBPs = []
	binsize1 = chromSizeInfo[chr1]
	binsize2 = chromSizeInfo[chr2]
	cutoff = np.nanpercentile(matrix,99)
	if cutoff == 0:
		pass
	else:
		for bp1 in mergedbp1s:
			for bp2 in mergedbp2s:
				#print bp1,bp2
				if rowsum[int(bp1)] > rowsum[int(bp1)+1]:
					x1 = int(bp1)-windowsize+1
					x2 = int(bp1)+1
				else:
					x1 = int(bp1)-windowsize
					x2 = int(bp1)
				if colsum[int(bp2)] > colsum[int(bp2)+1]:
					y1 = int(bp2)-windowsize+1
					y2 = int(bp2)+1
				else:
					y1 = int(bp2)-windowsize
					y2 = int(bp2)
				if x1 < 0:
					x1 = 0
				if y1 < 0:
					y1 = 0
				#print x1,x2,y1,y2
				r1 = matrix[x1:x2,y1:y2]
				mr1 = np.nanmean(r1)
				if rowsum[int(bp1)] > rowsum[int(bp1)+1]:
					x1 = int(bp1)+1
					x2 = int(bp1)+1+windowsize
				else:
					x1 = int(bp1)
					x2 = int(bp1)+windowsize
				if colsum[int(bp2)] > colsum[int(bp2)+1]:
					y1 = int(bp2)+1
					y2 = int(bp2)+1+windowsize
				else:
					y1 = int(bp2)
					y2 = int(bp2)+windowsize
				if x2 >= binsize1:
					x2 = binsize1
				if y2 >= binsize2:
					y1 = binsize2
				#print x1,x2,y1,y2
				r2 = matrix[x1:x2,y1:y2]
				mr2 = np.nanmean(r2)
				if rowsum[int(bp1)] > rowsum[int(bp1)+1]:
					x1 = int(bp1)-windowsize+1
					x2 = int(bp1)+1
				else:
					x1 = int(bp1)-windowsize
					x2 = int(bp1)
				if colsum[int(bp2)] > colsum[int(bp2)+1]:
					y1 = int(bp2)+1
					y2 = int(bp2)+1+windowsize
				else:
					y1 = int(bp2)
					y2 = int(bp2)+windowsize
				if x1 < 0:
					x1 = 0
				if y2 > binsize2:
					y1 = binsize2
				#print x1,x2,y1,y2
				r3 = matrix[x1:x2,y1:y2]
				mr3 = np.nanmean(r3)
				if rowsum[int(bp1)] > rowsum[int(bp1)+1]:
					x1 = int(bp1)+1
					x2 = int(bp1)+1+windowsize
				else:
					x1 = int(bp1)
					x2 = int(bp1)+windowsize
				if colsum[int(bp2)] > colsum[int(bp2)+1]:
					y1 = int(bp2)-windowsize+1
					y2 = int(bp2)+1
				else:
					y1 = int(bp2)-windowsize
					y2 = int(bp2)
				if x2 >= binsize1:
					x2 = binsize1
				if y1 < 0:
					y1 = 0
				#print x1,x2,y1,y2
				r4 = matrix[x1:x2,y1:y2]
				mr4 = np.nanmean(r4)
				###   mr1 | mr3
				###   ----|-----
				###   mr4 | mr2
				#print mr1,mr2,mr3,mr4
				count = 0
				mrs = [mr1,mr2,mr3,mr4]
				for mr in mrs:
					if np.isnan(mr):
						pass
					else:
						if mr > cutoff:
							count += 1
				#print count
				if count == 0:
					pass
				elif count >= 3:
					pass
				elif count == 2:
					mrs.sort()
					if (mr1 > cutoff and mr2 > cutoff) or (mr3> cutoff and mr4 > cutoff):
						rfold1 = relativefold(mrs[-2],mrs[-3])
						if rfold1 > 5:
							#print 'TRUE',bp1,bp2
							if rowsum[int(bp1)] > rowsum[int(bp1)+1]:
								correctbp1 = bp1
							else:
								correctbp1 = str(int(bp1)+1)
							if colsum[int(bp2)] > colsum[int(bp2)+1]:
								correctbp2 = bp2
							else:
								correctbp2 = str(int(bp2)+1)
							validBPs.append([chr1,str(int(correctbp1)+1),chr2,str(int(correctbp2)+1)])
					else:
						pass
				else:
					#print 'count: ', count
					sortedmrs = sorted(mrs, key = lambda x : float('-inf') if math.isnan(x) else x)
					#print sortedmrs
					rfold1 = 1
					rfold2 = 1
					if not math.isnan(sortedmrs[-2]) and sortedmrs[-2] != 0:
						rfold1 = relativefold(sortedmrs[-1],sortedmrs[-2])
						rfold2 = 1
						#print 'rfold1: ', rfold1
					else:
						rfold2 = relativefold(sortedmrs[-1],cutoff)
						#print 'rfold2: ', rfold2
						rfold1 = 1
					if rfold1 > 5 or rfold2 > 2:
						#print 'TRUE',bp1,bp2
						if rowsum[int(bp1)] > rowsum[int(bp1)+1]:
							correctbp1 = bp1
						else:
							correctbp1 = str(int(bp1)+1)
						if colsum[int(bp2)] > colsum[int(bp2)+1]:
							correctbp2 = bp2
						else:
							correctbp2 = str(int(bp2)+1)
						validBPs.append([chr1,str(int(correctbp1)+1),chr2,str(int(correctbp2)+1)])

	return validBPs

def getValidRoughBP(chromlengthf,matrix100kbInfo,background100kbInfo,rpInfo,outdir,name,cutoff,RscriptBPcaller):
	bpoutputfile,DivisionMatrixInfo = getAllRoughBreakpoints(matrix100kbInfo,background100kbInfo,rpInfo,outdir,name,cutoff,RscriptBPcaller)
	chromSizeInfo = getchromsize(chromlengthf,resolution = 100000)
	bpInfo = readBPs(bpoutputfile,chromSizeInfo)
	windowsize = 10
	validbpoutputfile = os.path.join(outdir,name + '_roughBP_100kb_filtered.txt')
	outf = open(validbpoutputfile,'w')
	for chrompair in bpInfo:
		if chrompair in DivisionMatrixInfo:
			print chrompair
			bps = bpInfo[chrompair]
			matrixfile = DivisionMatrixInfo[chrompair]
			validBPs = filtering(bps,matrixfile,chromSizeInfo,windowsize)
			for validbp in validBPs:
				res = [chrompair] + validbp
				outf.write('\t'.join(res) + '\n')
		else:
			print chrompair
	outf.close()

	return validbpoutputfile
