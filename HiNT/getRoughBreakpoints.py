import os,sys
from HiNT.corelib import *
from multiprocessing import Pool
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

def runBPcaller(params):
	#print params
	RscriptBPcaller,outputsubdir,matrixfile,tempbpoutputfile = params
	#print RscriptBPcaller,outputsubdir,matrixfile,tempbpoutputfile
	command = "Rscript %s %s %s %s"%(RscriptBPcaller,outputsubdir,matrixfile,tempbpoutputfile)
	print(command)
	run_cmd(command)
	chrompair_bps = open(tempbpoutputfile).readlines()
	os.remove(tempbpoutputfile)
	return chrompair_bps

def getAllRoughBreakpoints(matrix100kbInfo,background100kbInfo,rpInfo,outdir,name,cutoff,RscriptBPcaller,threads):
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
				getDivisionMatrix(mat1,mat2,fname)
				DivisionMatrixInfo[chrompair] = fname
			else:
				pass
	#print DivisionMatrixInfo
	bpoutputfile = os.path.join(outdir,name + '_roughBP_100kb.txt')
	matrixfiles = os.listdir(outputsubdir)
	allparamsInfo = []
	for matrixfile in matrixfiles:
		matrixfile_prefix = matrixfile.replace('_DivisionMatrix.txt','')
		tempbpoutputfile = os.path.join(outdir,matrixfile_prefix + '_roughBP_100kb.txt')
		allparamsInfo.append([RscriptBPcaller,outputsubdir,matrixfile,tempbpoutputfile])
	#print allparamsInfo[0]
	results = []
	p = Pool(threads)
	result = p.map_async(runBPcaller, allparamsInfo, callback=results.append)
	p.close()
	p.join()

	#print results
	with open(bpoutputfile,'w') as outf:
		for i,line in enumerate(results[0]):
			if i == 0:
				outf.write(line[0])
				outf.write(line[1])
			else:
				outf.write(line[1])
	return bpoutputfile,DivisionMatrixInfo

##################################################################################################
##############    scripts below this line are doing the breakpoints filtering   ##################

def getchromsize(chromlengthf,resolution):
	chromSizeInfo = {}
	inf = open(chromlengthf)
	for line in inf:
		line = line.strip().split('\t')
		chrom,size = line
		chromSizeInfo[chrom] = [int(size),int(int(size)/int(resolution)+1)]
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

def mergeValidBPs(validBPs,matrixfile):
	matrix = np.loadtxt(matrixfile)
	rowsum = np.nansum(matrix,axis=1)
	colsum = np.nansum(matrix,axis=0)

	tempBPs = []
	for validbp in validBPs:
		if validbp not in tempBPs:
			tempBPs.append(validbp)
		else:
			pass
	bppairs = [(i[1],i[3]) for i in tempBPs]
	chrompair = list(set([(i[0],i[2]) for i in tempBPs]))
	if len(chrompair) != 1:
		print(chrompair)

	#print tempBPs,chrompair

	bppairs.sort()
	bp1s = [int(i[0]) for i in bppairs]
	bp2s = [int(i[1]) for i in bppairs]

	differences1 = ['NA']+[(bp1s[i+1]-bp1s[i]) for i in range(len(bp1s)-1)]
	differences2 = ['NA']+[(bp2s[i+1]-bp2s[i]) for i in range(len(bp2s)-1)]

	#print differences1,differences2
	tempMerged = [i for i in bppairs]
	#print tempMerged
	for j in range(1,len(differences1)):
		if differences1[j] == 0 and differences2[j] == 1:
			#print colsum[bp2s[j-1]],colsum[bp2s[j]]
			if colsum[bp2s[j-1]] > colsum[bp2s[j]]:
				tempMerged.remove(bppairs[j])
			else:
				tempMerged.remove(bppairs[j-1])
		elif differences1[j] == 1 and differences2[j] == 0:
			#print rowsum[bp2s[j-1]] > rowsum[bp2s[j]]
			if rowsum[bp2s[j-1]] > rowsum[bp2s[j]]:
				tempMerged.remove(bppairs[j])
			else:
				tempMerged.remove(bppairs[j-1])
		else:
			pass

	#print tempMerged

	mergedValidBPs = []
	for bp in tempMerged:
		allinfo = [chrompair[0][0],str(bp[0]),chrompair[0][1],str(bp[1])]
		mergedValidBPs.append(allinfo)

	return mergedValidBPs

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
				print("give the correct axis, row or column")
			if suma > sumar and sumb < sumbr: # __|--|__
				fca = suma/sumar
				fcb = sumbr/sumb
				if suma > sumbr and fca > fcb:
					print("merge",np.nansum(matrix[a-1,:]), suma, sumar, np.nansum(matrix[b-1,:]),sumb, sumbr, fca, fcb)
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
	if float(b) != 0:
		rfold = (float(a) - float(b))/float(b)
	else:
		rfold = (float(a) - float(b))/(float(b)+1)
	return rfold

def filtering(bps,matrixfile,chromSizeInfo,windowsize):
	matrix = np.loadtxt(matrixfile)
	rowsum = np.nansum(matrix,axis=1)
	#print np.shape(rowsum)
	colsum = np.nansum(matrix,axis=0)
	#print np.shape(colsum)
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

	print(mergedbp1s,mergedbp2s)

	validBPs = []
	binsize1 = chromSizeInfo[chr1][1]
	binsize2 = chromSizeInfo[chr2][1]

	cutoff = np.nanpercentile(matrix,99)
	#print "cutoff:", cutoff

	for bp1 in mergedbp1s:
		for bp2 in mergedbp2s:
			sm1 = np.nanmedian(matrix[int(bp1)-5:int(bp1),int(bp2)+1:int(bp2)+6])
			sm2 = np.nanmedian(matrix[int(bp1)+1:int(bp1)+6,int(bp2)+1:int(bp2)+6])
			sm3 = np.nanmedian(matrix[int(bp1)-5:int(bp1),int(bp2)-5:int(bp2)])
			sm4 = np.nanmedian(matrix[int(bp1)+1:int(bp1)+6,int(bp2)-5:int(bp2)])
			#print sm1,sm2,sm3,sm4
			sms = np.sort(np.asarray([sm1,sm2,sm3,sm4]))
			translocationType = "unbalanced"
			if sms[2] > min(np.nanpercentile(matrix,80),2) and sms[1] < max(np.nanpercentile(matrix,60),1.5):
				translocationType = "balanced"
			#print translocationType

			if translocationType == "NULL":
				continue
			else:

				if translocationType == "balanced":
					if sm3 > sm4: #bp1 belongs to the left
					#if rowsum[int(bp1)] > rowsum[int(bp1)+1]:
						x1 = int(bp1)-windowsize+1
						x2 = int(bp1)+1
					else:
						x1 = int(bp1)-windowsize
						x2 = int(bp1)
				if translocationType == "unbalanced":
					if np.nanmean(rowsum[int(bp1)-5:int(bp1)]) > np.nanmean(rowsum[int(bp1)+1:int(bp1)+6]):
						x1 = int(bp1)-windowsize+1
						x2 = int(bp1)+1
					else:
						x1 = int(bp1)-windowsize
						x2 = int(bp1)

				if translocationType == "balanced":
					if sm3 > sm1:
						#if colsum[int(bp2)] > colsum[int(bp2)+1]:
						y1 = int(bp2)-windowsize+1
						y2 = int(bp2)+1
					else:
						y1 = int(bp2)-windowsize
						y2 = int(bp2)

				if translocationType == "unbalanced":
					if np.nanmean(colsum[int(bp2)-5:int(bp2)]) > np.nanmean(colsum[int(bp2)+1:int(bp2)+6]):
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

				if translocationType == "balanced":
					if sm1 > sm2:
						x1 = int(bp1)+1
						x2 = int(bp1)+1+windowsize
					else:
						x1 = int(bp1)
						x2 = int(bp1)+windowsize
				if translocationType == "unbalanced":
					if np.nanmean(rowsum[int(bp1)-5:int(bp1)]) > np.nanmean(rowsum[int(bp1)+1:int(bp1)+6]):
						x1 = int(bp1)+1
						x2 = int(bp1)+1+windowsize
					else:
						x1 = int(bp1)
						x2 = int(bp1)+windowsize
				if translocationType == "balanced":
					if sm4 > sm2:
						y1 = int(bp2)+1
						y2 = int(bp2)+1+windowsize
					else:
						y1 = int(bp2)
						y2 = int(bp2)+windowsize
				if translocationType == "unbalanced":
					if np.nanmean(colsum[int(bp2)-5:int(bp2)]) > np.nanmean(colsum[int(bp2)+1:int(bp2)+6]):
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

				if translocationType == "balanced":
					if sm1 > sm2:
						x1 = int(bp1)-windowsize+1
						x2 = int(bp1)+1
					else:
						x1 = int(bp1)-windowsize
						x2 = int(bp1)
				if translocationType == "unbalanced":
					if np.nanmean(rowsum[int(bp1)-5:int(bp1)]) > np.nanmean(rowsum[int(bp1)+1:int(bp1)+6]):
						x1 = int(bp1)-windowsize+1
						x2 = int(bp1)+1
					else:
						x1 = int(bp1)-windowsize
						x2 = int(bp1)
				if translocationType == "balanced":
					if sm3 > sm1:
						y1 = int(bp2)+1
						y2 = int(bp2)+1+windowsize
					else:
						y1 = int(bp2)
						y2 = int(bp2)+windowsize
				if translocationType == "unbalanced":
					if np.nanmean(colsum[int(bp2)-5:int(bp2)]) > np.nanmean(colsum[int(bp2)+1:int(bp2)+6]):
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
				mr3 = np.nanpercentile(r3,80,interpolation="lower")
				mr3 = np.nanmean(r3)

				if translocationType == "balanced":
					if sm3 > sm4:
						x1 = int(bp1)+1
						x2 = int(bp1)+1+windowsize
					else:
						x1 = int(bp1)
						x2 = int(bp1)+windowsize
				if translocationType == "unbalanced":
					if np.nanmean(rowsum[int(bp1)-5:int(bp1)]) > np.nanmean(rowsum[int(bp1)+1:int(bp1)+6]):
						x1 = int(bp1)+1
						x2 = int(bp1)+1+windowsize
					else:
						x1 = int(bp1)
						x2 = int(bp1)+windowsize
				if translocationType == "balanced":
					if sm4 > sm2:
						y1 = int(bp2)-windowsize+1
						y2 = int(bp2)+1
					else:
						y1 = int(bp2)-windowsize
						y2 = int(bp2)
				if translocationType == "unbalanced":
					if np.nanmean(colsum[int(bp2)-5:int(bp2)]) > np.nanmean(colsum[int(bp2)+1:int(bp2)+6]):
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
				#print(np.shape(r4))
				mr4 = np.nanpercentile(r4,80,interpolation="lower")
				mr4 = np.nanmean(r4)
				###   mr1 | mr3
				###   ----|-----
				###   mr4 | mr2

				#print mr1,mr2,mr3,mr4,cutoff
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
				elif count == 2: #balanced translocation, need a submatrix to calculate the rowsums and colsums
					mrs.sort()
					if (mr1 > cutoff and mr2 > cutoff) or (mr3> cutoff and mr4 > cutoff):
						rfold1 = relativefold(mrs[-2],mrs[-3])
						#print rfold1
						if rfold1 > 5:
							#print 'TRUE',bp1,bp2
							submatrix = matrix[:,0:int(bp2)]
							subrowsum = np.nansum(submatrix,1)
							submatrix = matrix[0:int(bp1),:]
							subcolsum = np.nansum(submatrix,0)

							if subrowsum[int(bp1)] > subrowsum[int(bp1)+1]:
								correctbp1 = bp1
							else:
								correctbp1 = str(int(bp1)+1)
							if subcolsum[int(bp2)] > subcolsum[int(bp2)+1]:
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

def getValidRoughBP(chromlengthf,matrix100kbInfo,background100kbInfo,rpInfo,outdir,name,cutoff,RscriptBPcaller,threads):
	bpoutputfile,DivisionMatrixInfo = getAllRoughBreakpoints(matrix100kbInfo,background100kbInfo,rpInfo,outdir,name,cutoff,RscriptBPcaller,threads)
	chromSizeInfo = getchromsize(chromlengthf,resolution = 100000)
	bpInfo = readBPs(bpoutputfile,chromSizeInfo)
	windowsize = 10
	validbpoutputfile = os.path.join(outdir,name + '_roughBP_100kb_filtered.txt')
	outf = open(validbpoutputfile,'w')
	#print DivisionMatrixInfo.keys()
	for chrompair in bpInfo:
		if chrompair in DivisionMatrixInfo:
			#print chrompair
			bps = bpInfo[chrompair]
			matrixfile = DivisionMatrixInfo[chrompair]
			validBPs = filtering(bps,matrixfile,chromSizeInfo,windowsize)
			if len(validBPs) == 1:
				mergedValidBPs = validBPs
			else:
				mergedValidBPs = mergeValidBPs(validBPs,matrixfile)
			for validbp in mergedValidBPs:
				res = [chrompair] + validbp
				outf.write('\t'.join(res) + '\n')
		else:
			print("This chromosomal pair (%s) is not valid"%chrompair)
	outf.close()

	return validbpoutputfile

