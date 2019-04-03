#We would like to regress rowsums with mappablity, GC content and the number of fragment is this region
import os,sys
import numpy as np
from sklearn import linear_model
import pandas as pd
import tabix
from HiNT.corelib import *
from pkg_resources import resource_filename

def get_chromInfo(chromlf):
	chroms = []
	infos = {}
	inf = open(chromlf)
	for line in inf:
		line = line.strip().split('\t')
		infos[line[0]] = line[1]
		chroms.append(line[0])
	return chroms,infos

def getGCpercent(gcfile,chrom,resolution):
	GCdata = pd.read_csv(gcfile, header = 0,sep="\t",index_col=0)
	GCper = GCdata.loc[chrom:chrom,'5_pct_gc']
	end = GCdata.loc[chrom:chrom,'3_usercol']
	length = end[-1:][0]
	steps = int(length/resolution)
	GC = []
	GC += [0]*(steps+1)
	for i in range(steps):
		binstart = i*(int(resolution/1000))
		binend = (i+1)*(int(resolution/1000))
		averageGC = np.average(GCper[binstart:binend])
		GC[i] = averageGC
	lastbinStart = steps*(int(resolution/1000))
	lastbinEnd = len(GCper)
	GC[-1] = np.average(GCper[lastbinStart:lastbinEnd])
	#print GCper
	GC = np.asarray(GC)

	return GC

def getmappability(mappablity_track, matrixchrom, matrixstart, matrixend, resolution):
	tb = tabix.open(mappablity_track)
	records = tb.query(matrixchrom, int(matrixstart), int(matrixend))
	Allstart = int(matrixstart)
	Allend = int(matrixend)
	values = []
	binvalues = []
	values += [0]*(Allstart%resolution)
	for r in records:
		chrom,start,end,value = r
		length = int(end) - int(start)
		if int(start) == Allstart:
			temp = [float(value)]*length
		else:
			temp = [0]*(int(start)-Allstart) + [float(value)]*length
		values += temp
		Allstart = int(end)
	if Allstart == Allend:
		pass
	else:
		values += [0]*(Allend-Allstart)
	values += [0]*(resolution-Allend%resolution)
	tempcount = 0
	for i,v in enumerate(values):
		if i%resolution == (resolution-1):
			binvalues.append(float(tempcount)/float(resolution))
			tempcount = 0
		else:
			tempcount += v
	return binvalues

def getRestrictionSitesInfo(restrictionSites):
	inf = open(restrictionSites)
	sitesInfo = {}
	for line in inf:
		line = line.strip().split()
		chrom = line[0]
		sites = line[1:]
		sitesInfo[chrom] = sites
	return sitesInfo

def getFragmentsNumber(sitesInfo, matrixchrom, resolution,chromlength):
	sites = sitesInfo[matrixchrom]
	bincounts = [0.0] * (int(int(chromlength)/int(resolution)) + 1)
	for site in sites:
		idx = int(int(site) / int(resolution))
		if idx >= len(bincounts):
			pass
		else:
			bincounts[idx] += 1.0
	return bincounts

def Regression(rowsums, mappability,gc,fragmentCounts,outputfilename, filesName):
	matrix = rowsums.values
	label = filesName[0]
	label = label.lstrip('#').strip()
	outRegressionfile = outputfilename + '%s.txt'%label

	#print outRegressionfile
	rowsum = np.asarray(matrix[:,0])
	nrowdim = len(rowsum)
	GC = np.asarray(gc[:nrowdim])
	mapscore = np.asarray(mappability[:nrowdim])
	fragmentscore = np.asarray(fragmentCounts[:nrowdim])
	#print len(rowsum),len(GC),len(mapscore),len(fragmentscore)

	d = {'rowsum':rowsum, 'gcper':GC, 'mapscore': mapscore,'cutSitesNumber':fragmentscore}
	mydata = pd.DataFrame(d,columns=['cutSitesNumber','gcper','mapscore','rowsum'])
	mydata.to_csv(outRegressionfile,sep="\t",header=False)

	return outRegressionfile

def getnonzeros(outputfilename,nonzerooutputfilename):
	inf = open(outputfilename)
	outf = open(nonzerooutputfilename,'w')
	for line in inf:
		if line.startswith(' '):
			pass
		else:
			newline = line.strip().split('\t')
			if len(newline) != 5:
				pass
				#print line
			else:
				id,cn,gc,ms,rs = newline
				if float(cn) == float(gc) == float(ms) == 0:
					pass
				else:
					outf.write('\t'.join(newline) + '\n')
	outf.close()

def mergeAllchroms(regressionFileInfo,chroms,outfile,headerfile):
	command = "cat %s "%headerfile
	for chrom in chroms:
		if chrom in regressionFileInfo:
			datafile = regressionFileInfo[chrom]
			command += datafile
			command += ' '
		else:
			pass
	command += '> %s'%outfile
	print(command)
	run_cmd(command)

def prepareData(name,outdir,referencedir,chromlf,rowSumFilesInfo,binsize,hg19_1k_GCPercent,mappablity_track,restrictionSites):
	regressionFileInfo = {}
	chroms,chromInfo = get_chromInfo(chromlf)
	newchroms = []
	suboutdir = os.path.join(outdir,name+'_dataForRegression')
	if os.path.isdir(suboutdir):
		outputname = os.path.join(suboutdir,"regressionInfo")
	else:
		os.mkdir(suboutdir)
		outputname = os.path.join(suboutdir,"regressionInfo")
	sitesInfo = getRestrictionSitesInfo(restrictionSites)
	for chrom in chromInfo:
		if chrom not in rowSumFilesInfo or chrom == 'chrY' or chrom == 'chrM':
			pass
		else:
			#print chrom
			chromLength = chromInfo[chrom]
			file = rowSumFilesInfo[chrom]
			rowsums = pd.read_csv(file, header = 0,sep="\t",index_col=0)
			filesName = list(rowsums)
			#print filesName
			resolution = int(binsize)*1000
			mapscore = getmappability(mappablity_track, chrom, 1, chromLength, resolution)
			fragmentCounts = getFragmentsNumber(sitesInfo, chrom, resolution, chromLength)
			gc = getGCpercent(hg19_1k_GCPercent,chrom,resolution)
			outputfilename = outputname + '_%s_'%chrom
			outRegressionfile = Regression(rowsums, mapscore, gc, fragmentCounts, outputfilename, filesName)
			nonzerooutputfilename = outRegressionfile.replace('.txt', '_nonzero.txt')
			getnonzeros(outRegressionfile, nonzerooutputfilename)
			regressionFileInfo[chrom] = nonzerooutputfilename
			newchroms.append(chrom)
	outfile = os.path.join(outputname + '_nonZeros_allchroms.txt')
	headerfile = os.path.join(referencedir, 'regressionDataHeader.txt')
	mergeAllchroms(regressionFileInfo,chroms,outfile,headerfile)
	regressionData = outfile
	return newchroms,regressionData,regressionFileInfo
