#!/home/sw229/miniconda2/bin/env python

#########################
####     Su Wang      ###
####    1-11-2018     ###
####   version 1.0    ###
#########################

"""Script Description:

HiNT - HiC for copy Number variations and Translocations detection

Part1:
Hi-C data preprocessing (fastq or bam format).
1, align to the genome via bwa-mem if the data format is fastq (skip if the data format is bam)
2, create valid read pairs (.pairs format) and chimeric read pairs (.pairsam format) with pairsamtools
3, create contact matrix with cooler or juicer tools in .cool or .hic format

Part2:
HiNT for CNV detection
1, extract 50kb (default) unnormalized contact matrices (.cool/.hic) or with input as a dense or sparse matrix format
2, calculate 1D coverage, GC content, mappability, and the number of restriction cutting sites in each 50kb bin (as default)
3, get residuals from GAM with Poisson linked function
4, segmentation with BIC-seq segmentation algorithm
5, visualization of CNVs

Part3:
HiNT for translocations detection
1, extract normalized Hi-C contact matrices in 100kb and 1Mb resolution with .cool or .hic format, or with the input as dense or sparse matrix fotmat
2, calculate the adjusted contact matrices with the in-house background contact matrices
3, calculate the rank product for each 1Mb resolution inter-chromosomal contact matrice to define the translocated chromosomal pairs
4, calculate the breakpoints based on contact matrices in 100kb resolution using a 1D coverage profile
5, Integrate with chimeric read pairs, refine the breakpoints to 1bp resolution


This code is free software; you can redistribute it and/or modify it.

@version: $1.0$
@author: Su Wang
@contact: wangsu0623@gmail.com
"""

import sys, os, time
import re
from pkg_resources import resource_filename
from HiNT.ArgsValidator import *
from HiNT.corelib import *

def prerun(args):
    opts = opt_validate_hintpre(args)
    from HiNT.HiCprocessingPip import *
    dataInfo = {}
    if opts.inputformat == "fastq":
        Info("Align to %s via BWA mem"%opts.genome)
        dataInfo = runBWA(opts,dataInfo)
    else:
        Info("Input data is in bam format, skip alignment step")
        dataInfo['bam'] = opts.hicdata
    dataInfo = runpairsamtools(dataInfo,opts)
    dataInfo = sortpairsam(dataInfo,opts)
    dataInfo = pairsIndex(dataInfo,opts)
    if opts.outputformat == 'cooler':
        dataInfo = runcooler(dataInfo,opts)
    elif opts.outputformat == 'juicer':
        dataInfo = runjuicer(dataInfo,opts)
    return

def cnvrun(args):
    opts = opt_validate_hintcnv(args)
    from HiNT.prepare_regression import *
    from HiNT.DoRegressionAllchroms import *
    from HiNT.doBICseq import *
    chromlf = resource_filename('HiNT', 'references/%s.len'%opts.genome)
    if opts.format == 'cooler':
        from HiNT.getGenomeRowSumsFromCool import *
        rowSumFilesInfo = getallChromsRowSums(opts.matrixfile,opts.name,opts.outdir,opts.resolution) #Calculate rowsums
    if opts.format == 'juicer':
        from HiNT.getGenomeRowSumsFromHiC import *
        rowSumFilesInfo = getGenomeRowSums(opts.resolution, opts.matrixfile, chromlf, opts.outdir,opts.name)
    #if opts.format == 'sparse':
    #    from HiNT.getGenomeRowSumsFromSparse import *

    binsize = opts.resolution
    GCPercent_1kb = resource_filename('HiNT', 'references/%s_1k_GCPercent.bed'%opts.genome)
    mappablity_track = resource_filename('HiNT', 'references/%s_wgEncodeCrgMapabilityAlign50mer.bdg.gz'%opts.genome)
    mappablity_trackIndex = resource_filename('HiNT', 'references/%s_wgEncodeCrgMapabilityAlign50mer.bdg.gz.tbi'%opts.genome)
    restrictionSites = resource_filename('HiNT', 'references/%s_%s_enzymeSites.txt'%(opts.genome,opts.enzyme))
    chroms,regressionFileAllchroms,regressionChromFilesInfo = prepareData(opts.name,opts.outdir,chromlf,rowSumFilesInfo,binsize,GCPercent_1kb,mappablity_track,restrictionSites) #Prepare the other data Information for regression
    print chroms
    resiudalChromFilesInfo = calculateResiduals(opts.name,opts.outdir,regressionFileAllchroms,regressionChromFilesInfo,chroms) #form fuction DoregressionAllchroms, get residuals per chromosome
    #Segmentation step
    BICseqPath = os.path.join(opts.bicseq, 'BICseq-seg.pl')
    CNVplotPath = resource_filename('HiNT', 'externalScripts/plot.l2r.ms_zoom.R')
    run_BICseq(opts.name,opts.outdir,opts.resolution,resiudalChromFilesInfo,chroms,BICseqPath,CNVplotPath)
    Info("Done! Find your CNV results from %s. ;)"%validBPregionOutf)
    return

def translrun(args):
    opts = opt_validate_hinttransl(args)
    chromlengthf = resource_filename('HiNT', 'references/%s.len'%opts.genome)
    from HiNT.getRankProduct import *
    from HiNT.getRoughBreakpoints import *
    if opts.format == 'cooler':
        from HiNT.coolToMatrix import *
        matrixfile1Mb,matrixfile100kb = opts.matrixfile
        matrix1MbInfo = coolToMatrix(matrixfile1Mb,1000,opts.outdir,opts.name)
        matrix100kbInfo = coolToMatrix(matrixfile100kb,100,opts.outdir,opts.name)
    if opts.format == 'juicer':
        from HiNT.hicToMatrix import *
        matrix1MbInfo = hicToMatrix(opts.matrixfile, 1000, chromlengthf, opts.outdir, opts.name)
        matrix100kbInfo = hicToMatrix(opts.matrixfile, 100, chromlengthf, opts.outdir, opts.name)
    backgroundMatrix1MbDir = os.path.join(opts.backgroundInterChromMatrixDir,'1Mb')
    backgroundMatrix100kbDir = os.path.join(opts.backgroundInterChromMatrixDir,'100kb')
    background1MbInfo = readBackgroundMatrix(backgroundMatrix1MbDir)
    background100kbInfo = readBackgroundMatrix(backgroundMatrix100kbDir)
    rpoutfile = getRankProduct(matrix1MbInfo,background1MbInfo,opts.outdir,opts.name)
    rpInfo = getRPinfo(rpoutfile)
    RscriptBPcallerPath = resource_filename('HiNT', 'externalScripts/getBreakPoints2steps.R')
    validBPregionOutf = getValidRoughBP(chromlengthf,matrix100kbInfo,background100kbInfo,rpInfo,opts.outdir,opts.name,opts.cutoff,RscriptBPcallerPath)
    if not opts.chimeric:
        Info("Done! Find your translocation breakpoints file from %s. ;)"%validBPregionOutf)
    else:
        from HiNT.getBPfromChimericReads import *
        from HiNT.BPsummarization import *
        restrictionSites = resource_filename('HiNT', 'references/%s_%s_enzymeSites.txt'%(opts.genome,opts.enzyme)) #This may need to be modified according to Dhawal's script
        BP_fromChimerasfile = getBPfromChimeras(opts.chimeric,restrictionSites,opts.enzyme,opts.outdir,opts.name,chromlengthf,validBPregionOutf,opts.pairixpath)
        IntegratedBPf = getSummarizedBP(BP_fromChimerasfile,opts.outdir,opts.name,rpoutfile,validBPregionOutf)
        Info("Done! Find your translocation summarized breakpoints file from %s. ;)"%IntegratedBPf)
    return True
