import os,sys
bpregionf = sys.argv[1]
bpchimericf = sys.argv[2]
giniIndexf = sys.argv[3]
outbpf = sys.argv[4]
resolution = int(sys.argv[5])

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

def readchimericf(bpchimericf):
	chimericInfo = {}
	inf = open(bpchimericf)
	for line in inf:
		line = line.strip().split('\t')
		chr1,searchregion1,bp1,chr2,searchregion2,bp2,supportread1,supportread2,supportreadtotal = line
		chrompair = chr1 + '_' + chr2
		try:
			chimericInfo[chrompair].append([chr1,searchregion1,bp1,chr2,searchregion2,bp2,supportread1,supportread2,supportreadtotal])
		except:
			chimericInfo[chrompair] = [[chr1,searchregion1,bp1,chr2,searchregion2,bp2,supportread1,supportread2,supportreadtotal]]

	return chimericInfo

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
				#print str((int(bp1)-1)*resolution),str((int(bp1)+1)*resolution)
				region1 = [str((int(bp1)-1)*resolution),str((int(bp1)+1)*resolution)]
				region2 = [str((int(bp2)-1)*resolution), str((int(bp2)+1)*resolution)]
				res = [chr1,'_'.join(region1),'-',chr2,'_'.join(region2),'-','0','0','0']
				#res = '\t'.join(res) + '\n'
				allres.append((res + [str("%.2e" % float(rankproduct))]))
				#outf.write('\t'.join(res) + '\n')
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
							#res = '\t'.join(chimericbp) + '\n'
							allres.append((res + [str("%.2e" % float(rankproduct))]))
							#outf.write('\t'.join(chimericbp) + '\n')
							#record.append(chimericbp)
							n += 1

				if n ==0:
					region1 = [str((int(bp1)-1)*resolution),str((int(bp1)+1)*resolution)]
					region2 = [str((int(bp2)-1)*resolution), str((int(bp2)+1)*resolution)]
					res = [chr1,'_'.join(region1),'-',chr2,'_'.join(region2),'-','0','0','0']
					#res = '\t'.join(res) + '\n'
					allres.append((res + [str("%.2e" % float(rankproduct))]))
					#outf.write('\t'.join(res) + '\n')
	allres_sorted = sorted(allres, key=lambda x: (-float(x[-1]),float(x[-2])), reverse=True)

	title = ['chrom1','region1','breakpoint1','chrom2','region2','breakpoint2','supportRead1','supportRead2',"supportReadTotal","RankProduct(p-value)"]
	outf.write('\t'.join(title) + '\n')

	written = []
	for ar in allres_sorted:
		if ar not in written:
			#outf.write('\t'.join(ar[:-1]) + '\n')
			outf.write('\t'.join(ar) + '\n')
			written.append(ar)
		else:
			pass
	outf.close()

def main():
	bpregionInfo = readbpregionf(bpregionf)
	chimericInfo = readchimericf(bpchimericf)
	giniInfo = readGiniIndexf(giniIndexf)
	integration(bpregionInfo,chimericInfo,giniInfo,outbpf,resolution)

if __name__ == '__main__':
	main()
