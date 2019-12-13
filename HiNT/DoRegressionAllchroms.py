import os,sys
from HiNT.corelib import *
import numpy as np
import math

def write_Rscript(Rscript,datafile,outname,outpdfname,regressionLog):
	Rscript += 'library(mgcv)\n'
	Rscript += 'inf <- "%s"\n'%datafile
	Rscript += 'sink("%s",append=T)\n'%regressionLog
	Rscript += 'inf\n'
	Rscript += 'sink()\n'
	Rscript += 'data <- read.table(inf,sep="\t",header=T)\n'
	Rscript += 'data <- data.frame(data)\n'
	Rscript += 'lm1 <- gam(rowsum ~ s(gcper) + s(mapscore) + s(cutSitesNumber),data=data, family=poisson(link=log))\n'
	Rscript += 'sink("%s",append=T)\n'%regressionLog
	Rscript += 'summary(lm1)\n'
	Rscript += 'sink()\n'
	Rscript += 'rs <- residuals.gam(lm1,type="working")\n'
	Rscript += 'write.table(rs,file="%s",quote=F,row.names=F,col.names=F)\n'%outname
	Rscript += 'pdf("%s")\n'%outpdfname
	Rscript += 'par(mfrow=c(2,2))\n'
	Rscript += 'plot(lm1)\n'
	Rscript += 'dev.off()\n'
	return Rscript

def DoRegression(name,outdir,regressionData):

	Rfile = os.path.join(outdir,name+'_dataForRegression',name+'_GAMregression.R')
	regressionLog = os.path.join(outdir,name+'_dataForRegression',name+'_GAMregression.log')
	Rscript = ''
	filep = regressionData
	filen = regressionData.rstrip('.txt')
	outname = filen+'_residuals.r'
	outpdfname = os.path.join(outdir, name + '_GAMPoisson.pdf')
	Rscript = write_Rscript(Rscript,filep,outname,outpdfname,regressionLog)
	outf = open(Rfile,'w')
	outf.write(Rscript)
	outf.close()
	command = "Rscript %s"%Rfile
	print(command)
	run_cmd(command)
	residualf = outname
	return residualf

def prepareIterativeRegression(regressionData,residualf,regression2Data):
	outf = open(regression2Data,'w')
	title = ["binIndex","cutSitesNumber","gcper","mapscore","rowsum","observed","expected","copyRatio"]
	outf.write('\t'.join(title)+'\n')
	regressionf = regressionData
	inbinf = open(regressionf)
	allbins = []
	allCutsizes = []
	allGCPer = []
	allMapScore = []
	allRowSum = []
	for line in inbinf:
		line = line.strip().split('\t')
		if len(line) != 5:
			pass
		else:
			bins,cutSites,gcPer,mapScore,rowsum = line
			allbins.append(int(bins))
			allCutsizes.append(cutSites)
			allGCPer.append(gcPer)
			allMapScore.append(mapScore)
			allRowSum.append(rowsum)
	matrix = np.loadtxt(residualf,delimiter='\t')
	mint = np.min(matrix)
	normtreat = matrix - mint
	expect = [0 - mint]*len(normtreat)
	normFactor = sum(expect)/sum(normtreat)
	if np.shape(normtreat)[0] != len(allbins):
		print(np.shape(normtreat)[0],len(allbins))
		print("size of the residuals and bins are not match!")
		sys.exit()
	else:
		for i,rowsum in enumerate(normtreat):
			binID = allbins[i]
			cutsites = allCutsizes[i]
			gcper = allGCPer[i]
			mapscore = allMapScore[i]
			regressionrowsum = allRowSum[i]
			obs = rowsum
			expected = expect[i]
			copyRatio = math.log((obs + 0.0000001)/(expected+0.0000001)*normFactor,2)
			row = [str(binID),cutsites,gcper,mapscore,regressionrowsum,str(obs),str(expected),str(copyRatio)]
			#print(row)
			outf.write('\t'.join(row) + '\n')
	outf.close()

def write_IterativeGAMRscript(Rscript,regression2Data,outname,outpdfname,regressionLog):
	Rscript += 'library(mgcv)\n'
	Rscript += 'library(gtools)\n'
	Rscript += 'inf <- "%s"\n'%regression2Data
	Rscript += 'sink("%s",append=T)\n'%regressionLog
	Rscript += 'inf\n'
	Rscript += 'sink()\n'
	Rscript += 'allRes <- read.table(inf,sep="\t",header=T)\n'
	Rscript += 'ValidData <- allRes[(allRes[,ncol(allRes)]<=0.3)&(allRes[,ncol(allRes)]>=-0.3),]\n'
	Rscript += 'data <- data.frame(ValidData)\n'
	Rscript += 'lm1 <- gam(rowsum ~ s(gcper) + s(mapscore) + s(cutSitesNumber),data=data, family=poisson(link=log))\n'
	Rscript += 'sink("%s",append=T)\n'%regressionLog
	Rscript += 'summary(lm1)\n'
	Rscript += 'sink()\n'
	Rscript += 'predictValues <- predict(lm1,allRes,type="response")\n'
	Rscript += 'newrs <- allRes$rowsum - predictValues\n'
	Rscript += 'write.table(newrs,file="%s",quote=F,row.names=F,col.names=F)\n'%outname
	Rscript += 'pdf("%s")\n'%outpdfname
	Rscript += 'par(mfrow=c(2,2))\n'
	Rscript += 'plot(lm1)\n'
	Rscript += 'dev.off()\n'
	return Rscript

def DoIterativeRegression(name,outdir,regressionData,residualf):
	Rfile = os.path.join(outdir,name+'_dataForRegression',name+'_GAMIterativeRegression.R')
	regressionLog = os.path.join(outdir,name+'_dataForRegression',name+'_GAMIterativeRegression.log')
	Rscript = ''
	regression2Data = os.path.join(outdir,name+'_dataForRegression','IterativeRegressionInfo_allChroms.txt')
	prepareIterativeRegression(regressionData,residualf,regression2Data)
	filep = regression2Data
	filen = regression2Data.rstrip('.txt')
	outname = filen+'_Residuals.r'
	outpdfname = os.path.join(outdir, name + '_GAMIterativeRegressionPoisson.pdf')
	Rscript = write_Rscript(Rscript,regression2Data,outname,outpdfname,regressionLog)
	outf = open(Rfile,'w')
	outf.write(Rscript)
	outf.close()
	command = "Rscript %s"%Rfile
	print(command)
	run_cmd(command)
	residualf = outname
	return residualf

def sepResidualsByChrom(residualf,chroms,regressionFileInfo):
	residualFileInfo = {}
	allresiduals = np.loadtxt(residualf)
	for chromfile in regressionFileInfo:
		infile = regressionFileInfo[chromfile]
		infileDir = os.path.dirname(infile)
		outfilename = os.path.basename(infile).replace('.txt','_residuals.r')
		residualDir = os.path.join(infileDir,"residualsPerChrom")
		if not os.path.isdir(residualDir):
			os.mkdir(residualDir)
		outfile = os.path.join(residualDir,outfilename)
		residualFileInfo[chromfile] = [infile,outfile]
	start = 0
	for chrom in chroms:
		if chrom not in residualFileInfo:
			pass
		else:
			infname,outfname = residualFileInfo[chrom]
			#print infname,outfname
			inf = np.loadtxt(infname)
			nrow,ncol = np.shape(inf)
			end = start + nrow
			chromresidual = allresiduals[start:end]
			np.savetxt(outfname, chromresidual,fmt='%.10f')
			start += nrow
	return residualFileInfo

def calculateResiduals(name,outdir,regressionData,regressionFileInfo,chroms, doIterative=False):
	AllchromResidualf = DoRegression(name,outdir,regressionData)
	residualChromFilesInfo = sepResidualsByChrom(AllchromResidualf,chroms,regressionFileInfo)
	if doIterative == True:
		AllchromResidualf = DoIterativeRegression(name,outdir,regressionData,AllchromResidualf)
		residualChromFilesInfo = sepResidualsByChrom(AllchromResidualf,chroms,regressionFileInfo)
	return residualChromFilesInfo
