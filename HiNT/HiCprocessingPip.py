import os,sys
from HiNT.corelib import *
from pkg_resources import resource_filename
#samtools path, pairixPath,bgzip,

def runBWA(opts,dataInfo):
    fq1,fq2 = opts.hicdata
    outbam = os.path.join(opts.outdir, opts.name + '.bam')
    command = '%s mem -SP5M -t %s %s %s %s | %s view -Shb - > %s'%(opts.bwapath,str(opts.threads),opts.bwaIndex,fq1,fq2,opts.samtoolspath,outbam)
    Info("Alignment with BWA-mem")
    print(command)
    run_cmd(command)
    dataInfo['fastq'] = opts.hicdata
    dataInfo['bam'] = outbam

    return dataInfo

def runpairsamtools(dataInfo,opts):
    inbam = dataInfo['bam']
    chromsizef = os.path.join(opts.referencedir, '%s.len'%opts.genome)
    outpairsam = os.path.join(opts.outdir,opts.name + '.pairsam.gz')
    command = '%s view -h %s | %s parse -c %s -o %s --nproc-in %s --nproc-out %s --assembly %s'%(opts.samtoolspath,inbam,opts.pairsampath,chromsizef,outpairsam, str(opts.threads), str(opts.threads), opts.genome)
    print(command)
    run_cmd(command)
    dataInfo['pairsam']=outpairsam

    return dataInfo

def sortpairsam(dataInfo,opts):
    chromsizef = os.path.join(opts.referencedir, '%s.len'%opts.genome)
    pairsamf = dataInfo['pairsam']
    sortedpairsam = pairsamf.rstrip('.pairsam.gz') + '.sorted.pairsam.gz'
    tmpdir = os.path.join(opts.outdir,'tmp')
    if not os.path.isdir(tmpdir):
        os.mkdir(tmpdir)
    command1 = "%s sort -o %s --memory 20G --compress-program gzip --nproc-in %s --nproc-out %s --tmpdir %s %s "%(opts.pairsampath,sortedpairsam,str(opts.threads), str(opts.threads),tmpdir,pairsamf)
    print(command1)
    run_cmd(command1)

    validpairsam = pairsamf.rstrip('.pairsam.gz') + '_valid.sorted.pairsam.gz'
    command2 = "%s select --nproc-in %s --nproc-out %s --chrom-subset %s -o %s '(pair_type == %s) or (pair_type == %s) or (pair_type == %s)' %s"%(opts.pairsampath,str(opts.threads),str(opts.threads),chromsizef,validpairsam, '"UU"','"UR"','"RU"',sortedpairsam)
    print(command2)
    run_cmd(command2)

    chimericpairsam = pairsamf.rstrip('.pairsam.gz') + '_chimeric.sorted.pairsam.gz'
    command3 = "%s select --nproc-in %s --nproc-out %s -o %s '(pair_type == %s) or (pair_type == %s) or (pair_type == %s)' %s"%(opts.pairsampath,str(opts.threads),str(opts.threads),chimericpairsam, '"NR"','"CC"','"MR"',sortedpairsam)
    print(command3)
    run_cmd(command3)

    validdeduppairsam = pairsamf.rstrip('.pairsam.gz') + '_valid.sorted.deduped.pairsam.gz'
    command4 = '%s dedup --nproc-in %s --nproc-out %s --output %s %s'%(opts.pairsampath,str(opts.threads),str(opts.threads),validdeduppairsam,validpairsam)
    print(command4)
    run_cmd(command4)

    validpairs = pairsamf.rstrip('.pairsam.gz') + '_merged_valid.pairs.gz' #output pairs format
    command5 = '%s split --nproc-in %s --nproc-out %s --output-pairs %s %s'%(opts.pairsampath,str(opts.threads),str(opts.threads),validpairs,validdeduppairsam)
    print(command5)
    run_cmd(command5)

    os.remove(pairsamf)
    os.remove(validpairsam)
    os.remove(validdeduppairsam)

    dataInfo['validPairs'] = validpairs
    dataInfo['chimericPairsam'] = chimericpairsam

    return dataInfo

def pairsIndex(dataInfo,opts):
    validpairs = dataInfo["validPairs"]
    #command1 = "bgzip -f %s"%(validpai^irs)
    #print(command1)
    #run_cmd(command1)

    command2 = "pairix -f %s"%(validpairs)
    print(command2)
    run_cmd(command2)
    dataInfo['bgzippedValidPairs'] = validpairs

    return dataInfo

def runcooler(dataInfo,opts):
    pairsfile = dataInfo['bgzippedValidPairs']
    chromsizef = os.path.join(opts.referencedir, '%s.len'%opts.genome)
    if not opts.resolution:
        resolutions = [50,100,1000]
    else:
        resolutions = [50,100,1000] + [opts.resolution]

    for i in range(len(resolutions)):
        binsizefile = os.path.join(opts.outdir,"bins.%skb.bed"%str(resolutions[i]))
        command1 = '%s makebins %s %s > %s'%(opts.coolerpath,chromsizef,resolutions[i]*1000,binsizefile)
        print(command1)
        run_cmd(command1)
        outcool = os.path.join(opts.outdir,opts.name+'_%skb.cool'%str(resolutions[i]))
        command2 = '%s cload pairix %s %s %s'%(opts.coolerpath,binsizefile,pairsfile,outcool)
        print(command2)
        run_cmd(command2)
        command3 = '%s balance %s'%(opts.coolerpath, outcool)
        print(command3)
        run_cmd(command3)
        dataInfo['%skb_cool'%(str(resolutions[i]))] = outcool

    return dataInfo

def runjuicer(dataInfo,opts):
    #Juicer Tools Version 1.7.6, java is required to be installed in advance
    pairsfile = dataInfo['bgzippedValidPairs']
    outhic = os.path.join(opts.outdir,opts.name+'.hic')
    command = 'java -jar %s pre -n %s %s %s'%(opts.juicerpath,pairsfile,outhic,opts.genome)
    print(command)
    run_cmd(command)
    command2 = 'java -jar %s addNorm -d -w 5000 -F %s'%(opts.juicerpath,outhic)
    print(command2)
    run_cmd(command2)
    dataInfo['hic'] = outhic

    return dataInfo
