import os,sys
import numpy as np
import math

matrixdir = sys.argv[1]
bpfile = sys.argv[2]
chromlengthf = sys.argv[3]
outbpfile = sys.argv[4]
resolution = sys.argv[5]

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

def getMatrixInfo(matrixdir):
	matrixInfo = {}
	matrixfile = os.listdir(matrixdir)
	for mf in matrixfile:
		mfInfo = mf.split('_')
		chrompair = '_'.join(mfInfo[0:2])
		#chrompair = '_'.join(mfInfo[2:4])
		matrixInfo[chrompair] = os.path.join(matrixdir,mf)
	return matrixInfo

def mergeBPs(bps,matrix,axis,windowsize=5):
	bps = [int(i)-1 for i in bps] #python index starts with 0
	bps.sort()
	#print bps

	mergedbps = [str(bps[0])]
	for i in range(1,len(bps)):
		a = int(bps[i])
		b = int(mergedbps[-1])
		print a,b
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
			'''
			if suma > sumar:
				sumaForCompare = suma
			else:
				sumaForCompare = sumar
			if sumb > sumbr:
				sumbForCompare = sumb
			else:
				sumbForCompare = sumbr
			if sumaForCompare >= sumbForCompare:
				mergedbps.append(str(a))
				mergedbps.remove(str(b))
			else:
				continue
			'''
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
	#print mergedbps
	return mergedbps

def relativefold(a,b):
	rfold = (float(a) - float(b))/float(b)
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
	#bp1s = [str(int(bp)-1) for bp in bp1] # all the bps add 1 when output in the R script
	#bp2s = [str(int(bp)-1) for bp in bp2] # all the bps add 1 when output in the R script
	bp1s = bp1
	bp2s = bp2
	bp1s = list(set(bp1s))
	bp2s = list(set(bp2s))

	print bp1s,bp2s
	mergedbp1s = mergeBPs(bp1s,matrix,'row')
	print mergedbp1s
	if str(int(max1)-1) not in mergedbp1s:
		mergedbp1s = mergedbp1s + [str(int(max1)-1)]

	mergedbp2s = mergeBPs(bp2s,matrix,'column')
	print mergedbp2s
	if str(int(max2)-1) not in mergedbp2s:
		mergedbp2s = mergedbp2s + [str(int(max2)-1)]

	print mergedbp1s,mergedbp2s

	validBPs = []
	binsize1 = chromSizeInfo[chr1]
	binsize2 = chromSizeInfo[chr2]

	cutoff = np.nanpercentile(matrix,99)
	print cutoff

	for bp1 in mergedbp1s:
		for bp2 in mergedbp2s:
			print bp1,bp2
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
			print x1,x2,y1,y2
			r1 = matrix[x1:x2,y1:y2]
			#print(np.shape(r1))
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
			print x1,x2,y1,y2
			r2 = matrix[x1:x2,y1:y2]
			#print(np.shape(r2))
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
			print x1,x2,y1,y2
			r3 = matrix[x1:x2,y1:y2]
			#print(np.shape(r3))
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
			print x1,x2,y1,y2
			r4 = matrix[x1:x2,y1:y2]
			#print(np.shape(r4))
			mr4 = np.nanmean(r4)

			###   mr1 | mr3
			###   ----|-----
			###   mr4 | mr2

			print mr1,mr2,mr3,mr4

			count = 0
			mrs = [mr1,mr2,mr3,mr4]
			for mr in mrs:
				if np.isnan(mr):
					pass
				else:
					if mr > cutoff:
						count += 1
			print count

			if count == 0:
				pass
			elif count >= 3:
				pass
				#mrs.sort()
				#
				#rfold1 = relativefold(mrs[-1],mrs[-2])
				#rfold2 = relativefold(mrs[-2],mrs[-3])
				#print rfold1,rfold2
				#if rfold1 >= 0.8:
				#	print 'TRUE',bp1,bp2
				#	validBPs.append([chr1,str(int(bp1)+1),chr2,str(int(bp2)+1)])
				#elif rfold2 >= 0.8 and rfold1 <= 0.5:
				#	print 'TRUE',bp1,bp2
				#	validBPs.append([chr1,str(int(bp1)+1),chr2,str(int(bp2)+1)])
				#else:
				#	pass

			elif count == 2:
				mrs.sort()
				if (mr1 > cutoff and mr2 > cutoff) or (mr3> cutoff and mr4 > cutoff):
					rfold1 = relativefold(mrs[-2],mrs[-3])
					if rfold1 > 5:
						print 'TRUE',bp1,bp2
						if rowsum[int(bp1)] > rowsum[int(bp1)+1]:
							correctbp1 = bp1
						else:
							correctbp1 = str(int(bp1)+1)
						if colsum[int(bp2)] > colsum[int(bp2)+1]:
							correctbp2 = bp2
						else:
							correctbp2 = str(int(bp2)+1)

						validBPs.append([chr1,str(int(correctbp1)+1),chr2,str(int(correctbp2)+1)])

				#elif relativefold(mrs[-1],mrs[-2]) >= 1:
				#	print 'TRUE',bp1,bp2
				#	validBPs.append([chr1,str(int(bp1)+1),chr2,str(int(bp2)+1)])
				else:
					pass
			else:
				print 'count: ', count
				sortedmrs = sorted(mrs, key = lambda x : float('-inf') if math.isnan(x) else x)
				print sortedmrs
				rfold1 = 1
				rfold2 = 1
				if not math.isnan(sortedmrs[-2]) and sortedmrs[-2] != 0:
					rfold1 = relativefold(sortedmrs[-1],sortedmrs[-2])
					rfold2 = 1
					print 'rfold1: ', rfold1
				else:
					rfold2 = relativefold(sortedmrs[-1],cutoff)
					print 'rfold2: ', rfold2
					rfold1 = 1

				if rfold1 > 5 or rfold2 > 2:
					print 'TRUE',bp1,bp2
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

def main():
	chromSizeInfo = getchromsize(chromlengthf,resolution)
	bpInfo = readBPs(bpfile,chromSizeInfo)
	matrixInfo = getMatrixInfo(matrixdir)
	outf = open(outbpfile,'w')

	if int(resolution) == 100000:
		windowsize = 10
	if int(resolution) == 1000000:
		windowsize = 5

	#print bpInfo
	#print matrixInfo
	for chrompair in bpInfo:
		if chrompair in matrixInfo:
		#if chrompair in matrixInfo and chrompair == 'chr3_chr10':
			bps = bpInfo[chrompair]
			print chrompair,bps
			matrixfile = matrixInfo[chrompair]
			validBPs = filtering(bps,matrixfile,chromSizeInfo,windowsize)
			#print validBPs
			for validbp in validBPs:
				res = [chrompair] + validbp
				outf.write('\t'.join(res) + '\n')
		else:
			print chrompair
		#break
	outf.close()

if __name__ == '__main__':
	main()
