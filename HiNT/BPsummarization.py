import os,sys
def readchemiricBPs(hicbpf):
	inf = open(hicbpf)
	Info = {}
	for line in inf:
		if line.startswith('chromosome1'):
			pass
		else:
			line = line.strip().split('\t')
			chr1,sr1,bp1,chr2,sr2,bp2,r1,r2,rt = line
			chrompair = chr1 + '_' + chr2
			if chrompair not in Info:
				Info[chrompair]=[[sr1,sr2,bp1,bp2,r1,r2,rt]]
			else:
				Info[chrompair].append([sr1,sr2,bp1,bp2,r1,r2,rt])
	# Info = {chrompair1:[[sr1,sr2,bp1,bp2,r1,r2,rt],[...],...],chrompair2:[[...],[...]]}
	return Info

def mergeCloseBPs(Info):
	#merge breakpoints that located within a 20bp windown centered from the compared one
	mergedInfo = {}
	for chrompair in Info:
		print(chrompair)
		bps = Info[chrompair]
		bps_sorted = sorted(bps, key=lambda x:(int(x[2]),int(x[3])))
		merged_bps = []
		merged_bps.append(bps_sorted[0])
		for i in range(1,len(bps_sorted)):
			sr11,sr21,bp11,bp21,r11,r21,rt1 = merged_bps[-1]
			sr12,sr22,bp12,bp22,r12,r22,rt2 = bps_sorted[i]
			if abs(int(bp12) - int(bp11) <= 10) and abs(int(bp22) - int(bp21) <= 10):
				if int(rt2) > int(rt1):
					merged_bps.remove(merged_bps[-1])
					merged_bps.append(bps_sorted[i])
				else:
					pass
			else:
				merged_bps.append(bps_sorted[i])
		sorted_mergedbps = sorted(merged_bps, key=lambda x:(int(x[-1])))
		mergedInfo[chrompair] = sorted_mergedbps

	return mergedInfo

def readchimericInfo(mergedInfo):
	chimericInfo = {}
	for k in mergedInfo:
		chr1,chr2 = k.split('_')
		values = mergedInfo[k]
		for value in values:
			searchregion1,searchregion2,bp1,bp2,supportread1,supportread2,supportreadtotal = value
			res = [chr1,searchregion1,bp1,chr2,searchregion2,bp2,supportread1,supportread2,supportreadtotal]
		try:
			chimericInfo[k].append([chr1,searchregion1,bp1,chr2,searchregion2,bp2,supportread1,supportread2,supportreadtotal])
		except:
			chimericInfo[k] = [[chr1,searchregion1,bp1,chr2,searchregion2,bp2,supportread1,supportread2,supportreadtotal]]

	return chimericInfo

def readGiniIndexf(giniIndexf):
	giniInfo = {}
	inf = open(giniIndexf)
	for line in inf:
		if line.startswith('GiniIndex'):
			pass
		else:
			line = line.strip().split('\t')
			chrompair, gini, maximum, rankproduct = line
			giniInfo[chrompair] = rankproduct
	return giniInfo

def readbpregionf(bpregionf):
	bpregionInfo = {}
	inf = open(bpregionf)
	for line in inf:
		line = line.strip().split('\t')
		chrompair,chr1,bp1,chr2,bp2 = line
		try:
			bpregionInfo[chrompair].append([chr1,bp1,chr2,bp2])
		except:
			bpregionInfo[chrompair] = [[chr1,bp1,chr2,bp2]]
	return bpregionInfo

def integration(bpregionInfo,chimericInfo,giniInfo,outbpf,resolution):
	outf = open(outbpf,'w')

	allres = []
	for chrompair in bpregionInfo:
		record = []
		regionbps = bpregionInfo[chrompair]
		rankproduct = giniInfo[chrompair]
		if chrompair not in chimericInfo:
			for regionbp in regionbps:
				chr1,bp1,chr2,bp2 = regionbp
				region1 = [str((int(bp1)-1)*resolution),str((int(bp1)+1)*resolution)]
				region2 = [str((int(bp2)-1)*resolution), str((int(bp2)+1)*resolution)]
				res = [chr1,'_'.join(region1),'-',chr2,'_'.join(region2),'-','0','0','0']
				allres.append((res + [str("%.2e" % float(rankproduct))]))
		else:
			chimericbps = chimericInfo[chrompair]
			for regionbp in regionbps:
				chr1,bp1,chr2,bp2 = regionbp
				regionbp1 = int(bp1)*resolution
				regionbp2 = int(bp2)*resolution
				n = 0
				for chimericbp in chimericbps:
					search_interval1 = chimericbp[1].split('_')
					search_interval2 = chimericbp[4].split('_')
					search_interval1 = [int(i) for i in search_interval1]
					search_interval2 = [int(i) for i in search_interval2]
					if (regionbp1 > search_interval1[0] and regionbp1 < search_interval1[1]) and (regionbp2 > search_interval2[0] and regionbp2 < search_interval2[1]):
						if chimericbp not in record:
							res = chimericbp
							allres.append((res + [str("%.2e" % float(rankproduct))]))
							n += 1
				if n ==0:
					region1 = [str((int(bp1)-1)*resolution),str((int(bp1)+1)*resolution)]
					region2 = [str((int(bp2)-1)*resolution), str((int(bp2)+1)*resolution)]
					res = [chr1,'_'.join(region1),'-',chr2,'_'.join(region2),'-','0','0','0']
					allres.append((res + [str("%.2e" % float(rankproduct))]))
	allres_sorted = sorted(allres, key=lambda x: (-float(x[-1]),float(x[-2])), reverse=True)

	title = ['chrom1','region1','breakpoint1','chrom2','region2','breakpoint2','supportRead1','supportRead2',"supportReadTotal","RankProduct(p-value)"]
	outf.write('\t'.join(title) + '\n')

	written = []
	for ar in allres_sorted:
		if ar not in written:
			outf.write('\t'.join(ar) + '\n')
			written.append(ar)
		else:
			pass
	outf.close()

def getSummarizedBP(validChimericbpf,outdir,name,giniIndexf,BPregionf):
	resolution = 100000
	summarizedChimericBPf = os.path.join(outdir,name + '_summarizedBP_fromChimeras.txt')
	chimericBPInfo = readchemiricBPs(validChimericbpf)
	mergedInfo = mergeCloseBPs(chimericBPInfo)
	FormatedChimericInfo = readchimericInfo(mergedInfo)
	giniInfo = readGiniIndexf(giniIndexf)
	bpregionInfo = readbpregionf(BPregionf)
	IntegratedBPf = os.path.join(outdir,name + '_Translocation_IntegratedBP.txt')
	integration(bpregionInfo,FormatedChimericInfo,giniInfo,IntegratedBPf,resolution)
	return IntegratedBPf
