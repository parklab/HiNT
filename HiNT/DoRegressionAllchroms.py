import os,sys
from HiNT.corelib import *
import numpy as np

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

def calculateResiduals(name,outdir,regressionData,regressionFileInfo,chroms):
	AllchromResidualf = DoRegression(name,outdir,regressionData)
	residualChromFilesInfo = sepResidualsByChrom(AllchromResidualf,chroms,regressionFileInfo)
	return residualChromFilesInfo
