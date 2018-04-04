# HiNT
### Hi-C for copy Number variations and Translocations detection

## Summary
HiNT is a computational approach, which can be applied to detect CNVs and inter-chromosomal translocations with breakpoints at single base pair resolution, within Hi-C data. HiNT has two main stages, stage1 does Hi-C data preprocessing, create Hi-C contact matrix, and stage2 detects CNVs or translocations or both in parallel.

#### overview of HiNT workflow: 

<img src="https://github.com/suwangbio/HiNT/blob/master/images/HiNToverview.png" width="600">



## Installation

### Dependencies
R and R packages

1. [R >= 3.4] (https://www.r-project.org/)
2. [strucchange] (https://cran.r-project.org/web/packages/strucchange/index.html) 
3. [parallel] (https://www.rdocumentation.org/packages/parallel/versions/3.4.1)
4. [Cairo] (https://cran.r-project.org/web/packages/Cairo/index.html)
5. [optparse] (https://cran.r-project.org/web/packages/optparse/index.html)


Python and Python packages

1. [python >= 2.7] (https://www.python.org/)

### Install
```$ python setup.py install ```
