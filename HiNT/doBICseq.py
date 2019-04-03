import os,sys
import numpy as np
from HiNT.corelib import *

def makebinFiles(residualf,regressionf,outbinfile,binsize):
	outf = open(outbinfile,'w')
	title = ['start','end','obs','expected','var']
	outf.write('\t'.join(title)+'\n')
	inbinf = open(regressionf)
	allbins = []
	for line in inbinf:
		line = line.strip().split('\t')
		if len(line) != 5:
			pass
		else:
			bins = int(line[0])
			allbins.append(bins)
	matrix = np.loadtxt(residualf,delimiter='\t')
	mint = np.min(matrix)
	normtreat = matrix - mint
	if np.shape(normtreat)[0] != len(allbins):
		print(np.shape(normtreat)[0],len(allbins))
		print("size of the residuals and bins are not match!")
		sys.exit()
	else:
		for i,rowsum in enumerate(normtreat):
			binID = allbins[i]
			start = binID*binsize
			end = (binID + 1)*binsize - 1
			obs = rowsum
			expected = 0 - mint
			var = obs - expected
			row = [str(start),str(end),str(obs),str(expected),str(var)]
			outf.write('\t'.join(row) + '\n')
	outf.close()

def BICseqPrepare(name,outdir,resolution,resiudalChromFilesInfo,chroms):
	subDir = os.path.join(outdir,"segmentation")
	if not os.path.isdir(subDir):
		os.mkdir(subDir)
	outputconfigfile = os.path.join(subDir,name + '_GAMregression_allchroms_%skb.cfg'%str(resolution))
	outf2 = open(outputconfigfile,'w')
	title2 = ['chrom','sample']
	outf2.write('\t'.join(title2) + '\n')
	outbinDir = os.path.join(subDir,'b%s'%str(resolution*1000))
	if not os.path.isdir(outbinDir):
		os.mkdir(outbinDir)
	for chrom in chroms:
		formatedchrom = chrom.lstrip('chr')
		if len(formatedchrom) == 1 and formatedchrom.isdigit():
			formatedchrom = '0' + formatedchrom
		if chrom not in resiudalChromFilesInfo:
			pass
		else:
			regressionf,residualf = resiudalChromFilesInfo[chrom]
			outbinfile = os.path.join(outbinDir,"%s_GAM_%skb.chrm_%s.b%s.bin"%(name,str(resolution),formatedchrom,str(resolution*1000)))
			makebinFiles(residualf,regressionf,outbinfile,resolution*1000)
			outf2.write(chrom+ '\t' + outbinfile + '\n')
	return subDir,outbinDir,outputconfigfile

def run_BICseq(name,outdir,resolution,resiudalChromFilesInfo,chroms,BICseqpath,CNVplotPath,genome):
	#print resiudalChromFilesInfo
	subDir,outbinDir,outputconfigfile = BICseqPrepare(name,outdir,resolution,resiudalChromFilesInfo,chroms)
	bicseqOut = os.path.join(subDir,name + '_CNV_segments.txt')
	tmp = os.path.join(subDir,'tmp')
	if not os.path.isdir(tmp):
		os.mkdir(tmp)
	outfig = os.path.join(subDir,name + '_%skb_CNV_segments.png'%str(resolution))
	command = 'perl %s %s %s --tmp %s --title %s_%skb_CNV --bootstrap --fig %s --lambda 3'%(BICseqpath,outputconfigfile,bicseqOut,tmp,name,str(resolution),outfig)
	print(command)
	run_cmd(command)
	#command2 = 'Rscript %s --lamda 3 --title %s_%skb_GAMPoissonResiduals -i %s -o %s --bin_size %s --chrm'%(CNVplotPath,name,str(resolution),outbinDir,bicseqOut,str(int(resolution)*1000))
	if genome.startswith('hg'):
		mouse = "0"
	if genome.startswith('mm'):
		mouse = "1"
	command2 = 'Rscript %s %s %s %s 3 %s %s_%skb_GAMPoissonResiduals --chrm'%(CNVplotPath, outbinDir, bicseqOut, str(int(resolution)*1000), mouse, name,str(resolution))
	print(command2)
	run_cmd(command2)
	outfig2 = os.path.join(subDir,name + '_CNV_segments.l2r.pdf')
	return bicseqOut,outfig,outfig2
