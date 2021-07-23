# HiNT 
[![DOI](https://zenodo.org/badge/240967519.svg)](https://zenodo.org/badge/latestdoi/240967519)
## A computational method for detecting copy number variations and translocations from Hi-C data

## Summary
**HiNT** (**Hi**-C for copy **N**umber variation and **T**ranslocation detection), a computational method to detect CNVs and Translocations from Hi-C data. HiNT has three main components: **HiNT-PRE**, **HiNT-CNV**, and **HiNT-TL**. HiNT-PRE preprocesses Hi-C data and computes the contact matrix, which stores contact frequencies between any two genomic loci; both HiNT-CNV and HiNT-TL starts with HI-C contact matrix, predicts copy number segments, and inter-chromosomal translocations, respectively 

#### Overview of HiNT workflow: 
<img src="https://github.com/suwangbio/HiNT/blob/master/images/HiNT_workflow.png" width="800">

## Installation

### Dependencies
R and R packages

1. [R >= 3.4](https://www.r-project.org/)
2. [mgcv](https://cran.r-project.org/web/packages/mgcv/index.html), [strucchange](https://cran.r-project.org/web/packages/strucchange/index.html), [doParallel](https://www.rdocumentation.org/packages/parallel/versions/3.4.1), [Cairo](https://cran.r-project.org/web/packages/Cairo/index.html), [foreach](https://cran.r-project.org/web/packages/optparse/index.html)

Python and Python packages

1. [python >= 3.5](https://www.python.org/)
2. [pyparix >= 0.3.0](https://github.com/4dn-dcic/pairix#pypairix), [cooler >= 0.7.4](https://github.com/mirnylab/cooler), [pairtools >= 0.2.2](https://github.com/mirnylab/pairtools), [numpy](https://www.scipy.org/install.html), [scipy](https://www.scipy.org/install.html), [pandas](https://pandas.pydata.org/), [sklearn](https://scikit-learn.org/stable/install.html), [multiprocessing](https://pypi.org/project/multiprocess/)

Java and related tools (Optional: required when want to process Hi-C data with juicer tools)

1. [Java (version >= 1.7)](https://www.java.com/en/download/)
2. [Juicer tools (1.8.9 is recommended)](https://github.com/aidenlab/juicer/wiki/Download) 

Perl

1. [Perl (version >= 5)](https://www.perl.org/)


Other dependencies

1. [samtools](http://www.htslib.org/doc/) (1.3.1+)
2. [BIC-seq2](http://www.math.pku.edu.cn/teachers/xirb/downloads/software/BICseq2/BICseq2/BICseq2-seg_v0.7.3.tar.gz) (0.7.3) ! This is optional: if you don't want to run HiNT-CNV, you don't need this package. [Download BICseq2, unzip it, and give the path of BICseq2-seg_v0.7.3 (/path/to/BICseq2-seg_v0.7.3)].
3. [bwa](https://sourceforge.net/projects/bio-bwa/files/) (0.7.16+) ! This is optional: required only when your input is fastq
4. [tabix](https://sourceforge.net/projects/samtools/files/tabix/) (0.2.6)

### Install HiNT

* Method1: Install using conda (highly recommended)

	``` $ conda install -c su hint```

	or

	``` $ conda install hint```
	
* Method2: Install from PyPI using pip.

	``` $ pip install HiNT-Packages```

* Method3: Install manually 
  1. Install HiNT dependencies
  2. Download HiNT ```git clone https://github.com/parklab/HiNT.git```
  3. Go to HiNT directory, install it by ```$ python setup.py install ```
  
  *** Type ```$ hint``` to test if HiNT successfully installed

* Method 4: Run HiNT in a Docker container (highly recommended)

	``` $ docker pull suwangbio/hint```

	``` $ docker run suwangbio/hint hint```

	See details of the usage on HiNT page at [docker hub](https://hub.docker.com/r/suwangbio/hint) 

### Download reference files used in HiNT [HERE](http://compbio.med.harvard.edu/hint/)

1. Download HiNT references [HERE](http://compbio.med.harvard.edu/hint/refData/). Only hg19, hg38 and mm10 are available currently. Unzip it ```$ unzip hg19.zip ```
2. Download HiNT background matrices [HERE](http://compbio.med.harvard.edu/hint/backgroundMatrices/). Only hg19, hg38 and mm10 are available currently. Unzip it ```$ unzip hg19.zip ```
3. Download BWA index files [HERE](http://compbio.med.harvard.edu/hint/BWAIndex/). Only hg19, hg38 and mm10 are available currently. Unzip it ```$ unzip hg19.zip ```

## Quick Start

* Download the test datasets from [HERE](https://www.dropbox.com/sh/z1rceh8ddnsdtj7/AAC0VuDu48eh_RtzKHipztkLa?dl=0)

### HiNT-PRE
HiNT pre: Preprocessing Hi-C data. HiNT pre does alignment, contact matrix creation and normalization in one command line.

```$ hint pre -d /path/to/hic_1.fastq.gz,/path/to/hic_2.fastq.gz -i /path/to/bwaIndex/hg19/hg19.fa --refdir /path/to/refData/hg19 --informat fastq --outformat cooler -g hg19 -n test -o /path/to/outputdir --pairtoolspath /path/to/pairtools --samtoolspath /path/to/samtools --coolerpath /path/to/cooler```

```$ hint pre -d /path/to/test.bam --refdir /path/to/refData/hg19 --informat bam --outformat juicer -g hg19 -n test -o /path/to/outputdir --pairtoolspath /path/to/pairtools --samtoolspath /path/to/samtools --juicerpath /path/to/juicer_tools.1.8.9_jcuda.0.8.jar```

use ```$ which samtools ``` ```$ which pairtools ``` ```$ which cooler ``` to get the absolute path of these tools, and ```/path/to/juicer_tools.1.8.9_jcuda.0.8.jar``` should be the path where you store this file

see details and more options

```$ hint pre -h ```

### HiNT-CNV
HiNT cnv: prediction of copy number information, as well as segmentation from Hi-C.

```$ hint cnv -m contactMatrix.cool -f cooler --refdir /path/to/refDir/hg19 -r 50 -g hg19 -n test -o /path/to/outputDir --bicseq /path/to/BICseq2-seg_v0.7.3 -e MboI```

```$ hint cnv -m /path/to/4DNFIS6HAUPP.mcool::/resolutions/50000 -f cooler --refdir /path/to/refDir/hg38 -r 50 -g hg38 -n HepG2 --bicseq /path/to/BICseq2-seg_v0.7.3 -e DpnII --maptrack 36mer```

```$ hint cnv -m /path/to/4DNFICSTCJQZ.hic -f juicer --refdir /path/to/refDir/hg38 -r 50 -g hg38 -n HepG2 --bicseq /path/to/BICseq2-seg_v0.7.3 -e DpnII```

```$ hint cnv -m /path/to/4DNFICSTCJQZ.hic -f juicer --refdir /path/to/refDir/hg38 -r 50 -g hg38 -n HepG2 --bicseq /path/to/BICseq2-seg_v0.7.3 -e DpnII --doiter```

```/path/to/BICseq2-seg_v0.7.3``` should be the path where you store this package

see details and more options

```$ hint cnv -h ```

### HiNT-TL
HiNT tl: interchromosomal translocations and breakpoints detection from
Hi-C inter-chromosomal interaction matrices.

```$ hint tl -m /path/to/data_1Mb.cool,/path/to/data_100kb.cool --chimeric /path/to/test_chimeric.sorted.pairsam.gz --refdir /path/to/refDir/hg19 --backdir /path/to/backgroundMatrices/hg19 --ppath /path/to/pairix -f cooler -g hg19 -n test -o /path/to/outputDir```

```$ hint tl -m /path/to/4DNFIS6HAUPP.mcool::/resolutions/1000000,/path/to/4DNFIS6HAUPP.mcool::/resolutions/100000 -f cooler --refdir /path/to/refDir/hg38 --backdir /path/to/backgroundMatrices/hg38 -g hg38 -n 4DNFICSTCJQZ -c 0.05 --ppath /path/to/pairix -p 12```

```$ hint tl -m /path/to/4DNFICSTCJQZ.hic -f juicer --refdir /path/to/refData/hg38 --backdir /path/to/backgroundMatrices/hg38 -g hg38 -n 4DNFICSTCJQZ -c 0.05 --ppath /path/to/pairix -p 12 -o HiNTtransl_juicerOUTPUT```

use ```$ which pairix ``` to get the absolute path of pairix

see details and more options

```$ hint tl -h ```

## Output of HiNT
### HiNT-PRE output
In the HiNT-PRE output directory, you will find

1. ```jobname.bam``` aligned lossless file in bam format
2. ```jobname_merged_valid.pairs.gz``` reads pairs in pair format
3. ```jobname_chimeric.sorted.pairsam.gz``` ambiguous chimeric read pairs used for breakpoint detection in [pairsam](https://github.com/mirnylab/pairtools) format
4. ```jobname_valid.sorted.deduped.pairsam.gz``` valid read pairs used for Hi-C contact matrix creation in [pairsam](https://github.com/mirnylab/pairtools) format
5. ```jobname.mcool``` Hi-C contact matrix in [cool](https://github.com/mirnylab/cooler) format
6. ```jobname.hic``` Hi-C contact matrix in [hic](https://github.com/aidenlab/juicer) format

### HiNT-CNV output
In the HiNT-CNV output directory, you will find

1. ```jobname_GAMPoisson.pdf``` the GAM regression result
2. ```segmentation/jobname_bicsq_allchroms.txt``` CNV segments with log2 copy ratio and p-values in txt file
3. ```segmentation/jobname_resolution_CNV_segments.png``` figure to visualize CNV segments
4. ```segmentation/jobname_bicseq_allchroms.l2r.pdf``` figure to visualize log2 copy ration in each bin (bin size = resolution you set)
5. ```segmentation/other_files``` intermediate files used to run BIC-seq
6. ```jonname_dataForRegression/*``` data used for regression as well as residuals after removing Hi-C biases

### HiNT-TL output
In the HiNT-TL output directory, you will find

1. ```jobname_Translocation_IntegratedBP.txt``` the final integrated translocation breakpoint
2. ```jobname_chrompairs_rankProduct.txt``` rank product predicted potential translocated chromosome pairs
3. ```otherFolders``` intermediate files used to identify the translocation breakpoints
