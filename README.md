# HiNT
### Hi-C for copy Number variations and Translocations detection

## Summary
**HiNT** (**Hi**-C for copy **N**umber variation and **T**ranslocation detection), a computational method to detect CNVs and Translocations from Hi-C data. HiNT has three main components: **HiNT-PRE**, **HiNT-CNV**, and **HiNT-TL**. HiNT-PRE preprocesses Hi-C data and computes the contact matrix, which stores contact frequencies between any two genomic loci; both HiNT-CNV and HiNT-TL starts with HI-C contact matrix, predicts copy number segments, and inter-chromosomal translocations, respectively 

#### overview of HiNT workflow: 

<img src="https://github.com/suwangbio/HiNT/blob/master/images/HiNT_workflow.png" width="600">



## Installation

### Dependencies
R and R packages

1. [R >= 3.4](https://www.r-project.org/)
2. [strucchange](https://cran.r-project.org/web/packages/strucchange/index.html) 
3. [parallel](https://www.rdocumentation.org/packages/parallel/versions/3.4.1)
4. [Cairo](https://cran.r-project.org/web/packages/Cairo/index.html)
5. [optparse](https://cran.r-project.org/web/packages/optparse/index.html)


Python and Python packages

1. [python >= 2.7](https://www.python.org/)
2. [pyparix](https://github.com/4dn-dcic/pairix#pypairix)
3. [cooler](https://github.com/mirnylab/cooler)
4. [pairtools](https://github.com/mirnylab/pairtools)
5. [numpy](https://www.scipy.org/install.html), [scipy](https://www.scipy.org/install.html), [pandas](https://pandas.pydata.org/), [sklearn](https://scikit-learn.org/stable/install.html), [multiprocessing](https://pypi.org/project/multiprocess/)

### Install
```$ python setup.py install ```


## Quick Start
### HiNT-PRE
HiNT pre: Preprocessing Hi-C data. HiNT pre does alignment, contact matrix creation and normalization in one command line.

```$ hint pre -d /path/to/hic_1.fastq.gz,/path/to/hic_2.fastq.gz -i /path/to/bwaIndex --informat fastq --outformat cooler -g hg19 -n test -o /path/to/outputdir --pairsampath /path/to/pairsamtools```

see details and more options

```$ hint pre -h ```

### HiNT-CNV
HiNT cnv: prediction of copy number information, as well as segmentation from Hi-C.

```$ hint cnv -m contactMatrix.mcool -f cooler -r 50 -g hg19 -n test -o /path/to/outputDir```

see details and more options

```$ hint cnv -h ```

### HiNT-TL
HiNT transl: interchromosomal translocations and breakpoints detection from
Hi-C inter-chromosomal interaction matrices.

```hint transl -m /path/to/data_1Mb.cool,/path/to/data_100kb.cool -c chimericReads.pairsam -f cooler -g hg19 -n test -o /path/to/outputDir```

see details and more options

```$ hint transl -h ```

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
3. ```segmentation/jobname_resolution_CNV_l03.png``` figure to visualize CNV segments
4. ```segmentation/jobname_bicseq_allchroms.l2r.pdf``` figure to visualize log2 copy ration in each bin (bin size = resolution you set)
5. ```segmentation/other_files``` intermediate files used to run BIC-seq
6. ```jonname_dataForRegression/*``` data used for regression as well as residuals after removing Hi-C biases

### HiNT-TL output
In the HiNT-TL output directory, you will find

1. ```jobname_Translocation_IntegratedBP.txt``` the final integrated translocation breakpoint
2. ```jobname_chrompairs_rankProduct.txt``` rank product predicted potential translocated chromosome pairs
3. ```otherFolders``` intermediate files used to identify the translocation breakpoints
