Args <- commandArgs()
#print(Args)
divisionMatrixDir = Args[6]
matrix = Args[7]
outputfile = Args[8]

#start_time <- Sys.time()

library('strucchange')
#library(pastecs)
library(foreach)
library(doParallel)

quantileFunc <- function(x,q) {
  x <- sort(x)
  n <- length(x)
  ln = floor(n*q)
  res = x[ln]
  return(res)
 }

svtype <- function(matrixf,chrompair){
	traltype="unbalancedTRAL"
	mat <- read.table(matrixf)
	rowsum = rowSums(mat,na.rm=T)
	colsum = colSums(mat,na.rm=T)
	maxRIndex = which(rowsum == max(rowsum,na.rm=T))
	highRIndex = which(rowsum == quantileFunc(rowsum,0.999))
	maxCIndex = which(colsum == max(colsum,na.rm=T))
	highCIndex = which(colsum == quantileFunc(colsum,0.999))

	if(abs(maxRIndex-highRIndex)<=length(rowsum)*0.005 && (median(rowsum[(highRIndex-5):(highRIndex-1)],na.rm=T) > quantileFunc(rowsum,0.9)) && (median(rowsum[(maxRIndex+1):(maxRIndex+5)],na.rm=T) > quantileFunc(rowsum,0.9)) && abs(maxCIndex-highCIndex)<=length(colsum)*0.005 && (median(colsum[(maxCIndex-5):(maxCIndex-1)],na.rm=T) > quantileFunc(colsum,0.9)) && (median(colsum[(maxCIndex+1):(maxCIndex+5)],na.rm=T) > quantileFunc(colsum,0.9))){
		traltype = "balancedTRAL"
		#print("Probably a balanced translocation")
	}
	return(traltype)
}

filtering <- function(breakpoints, array){
	bps = c()
	cutoff <- quantile(array[array>0],0.9)
	bcutoff <- median(array[array>0])

	#cutoff <- sort(array)[ceiling(0.9*length(sort(array)))]
	for (bp in breakpoints){
		if(bp==length(array)|(bp==1)|bp==length(array)-1|(bp==2)){bps <- bps}
		else if((array[bp]<=bcutoff & array[bp+1]>cutoff)|(array[bp]>cutoff & array[bp+1]<=bcutoff)){bps <- c(bps,bp)}
		#else if((array[bp]<=bcutoff & array[bp+2]>cutoff)|(array[bp]>cutoff & array[bp+2]<=bcutoff)){bps <- c(bps,bp)}
		else if((array[bp-1]>=cutoff & array[bp]<bcutoff)|(array[bp-1]<bcutoff & array[bp]>=cutoff)){bps <- c(bps,bp)}
		#else if((array[bp-2]>=cutoff & array[bp]<bcutoff)|(array[bp-2]<bcutoff & array[bp]>=cutoff)){bps <- c(bps,bp)}

		else if((array[bp-1]>=cutoff & array[bp+1]<bcutoff)|(array[bp-1]<bcutoff & array[bp+1]>=cutoff)){bps <- c(bps,bp)}

		else if((min(array[bp],array[bp+1])!=0) & (abs(array[bp+1]-array[bp])/(max(array[bp+1],array[bp])+1)>0.3)){bps <- c(bps,bp)}
		#else if((min(array[bp],array[bp+2])!=0) & (abs(array[bp+2]-array[bp])/(max(array[bp+2],array[bp])+1)>0.3)){bps <- c(bps,bp)}
		else if((min(array[bp],array[bp-1])!=0) & (abs(array[bp-1]-array[bp])/(max(array[bp-1],array[bp])+1)>0.3)){bps <- c(bps,bp)}
		#else if((min(array[bp],array[bp-2])!=0) & (abs(array[bp-2]-array[bp])/(max(array[bp-2],array[bp])+1)>0.3)){bps <- c(bps,bp)}
		else {bps <- bps}
	}

	if(is.null(bps)){bps <- c(which(array==max(array)))}
	else if(length(bps)==1){bps <- bps}
	else{
		ids <- c()
		for(i in 1:(length(bps)-1)){
			if((bps[i+1]-bps[i])==1){if(array[bps[i]]<array[bps[i+1]]){ids<-c(ids,i)}else{ids<-c(ids,i+1)}}
			else{ids <- ids}
		}
		if(is.null(ids)){bps <- bps}
		else{bps <- bps[-(ids)]}
	}
	return(bps)
}

searchMatric <- function(matrixf,chrompair,threads){
	mat <- read.table(matrixf)
	h = dim(mat)[1]
	w = dim(mat)[2]
	# To check the breakpoints
	if(h>w){
		mat <- t(mat)
	}else{
		mat <- mat
	}
	#print(dim(mat))
	rowsum = rowSums(mat,na.rm=T)
	colsum = colSums(mat,na.rm=T)

	registerDoParallel(cores=strtoi(threads))
	breakpoints_row <- breakpoints(rowsum ~ 1, h=5, breaks=10, hpc="foreach")$breakpoints
	breakpoints_col <- breakpoints(colsum ~ 1, h=5, breaks=10, hpc="foreach")$breakpoints

	max_row <- which(rowsum==max(rowsum))
	max_col <- which(colsum==max(colsum))
	validRowBP <- filtering(breakpoints_row,rowsum)
	validColBP <- filtering(breakpoints_col,colsum)
	bprow <- paste(validRowBP,collapse = ',')
	bpcol <- paste(validColBP,collapse = ',')
	#bprow <- paste(breakpoints_row, collapse = ',')
	#bpcol <- paste(breakpoints_col, collapse = ',')
	res <- cbind(chrompair,bprow,max_row,bpcol,max_col)
	return(res)
}

searchBalancedMatric <- function(matrixf,chrompair,threads){
	mat <- read.table(matrixf)
	h = dim(mat)[1]
	w = dim(mat)[2]
	# To check the breakpoints
	if(h>w){
		mat <- t(mat)
	}else{
		mat <- mat
	}
	#print(dim(mat))
	rowsum = rowSums(mat,na.rm=T)
	colsum = colSums(mat,na.rm=T)

	maxRIndex = which(rowsum == max(rowsum,na.rm=T))
	if(sum(is.finite(rowsum[1:maxRIndex]))>sum(is.finite(rowsum[maxRIndex:length(rowsum)]))){
		submatrix <- mat[1:(maxRIndex-1),]
	}else{
		submatrix <- mat[(maxRIndex+1):length(rowsum),]
	}
	subcolSum <- colSums(submatrix,na.rm=T)

	maxCIndex = which(colsum == max(colsum,na.rm=T))
	if(sum(is.finite(colsum[1:maxCIndex]))>sum(is.finite(colsum[maxCIndex:length(colsum)]))){
		submatrix <- mat[,1:(maxCIndex-1)]
	}else{
		submatrix <- mat[,(maxCIndex+1):length(colsum)]
	}
	subrowSum <- rowSums(submatrix,na.rm=T)

	registerDoParallel(cores=strtoi(threads))
	breakpoints_row <- breakpoints(subrowSum ~ 1, h=5, breaks=10, hpc="foreach")$breakpoints
	breakpoints_col <- breakpoints(subcolSum ~ 1, h=5, breaks=10, hpc="foreach")$breakpoints

	max_row <- maxRIndex
	max_col <- maxCIndex

	#validRowBP <- filtering(breakpoints_row,rowsum)
	#validColBP <- filtering(breakpoints_col,colsum)
	#bprow <- paste(validRowBP,collapse = ',')
	#bpcol <- paste(validColBP,collapse = ',')

	bprow <- paste(breakpoints_row, collapse = ',')
	bpcol <- paste(breakpoints_col, collapse = ',')

	res <- cbind(chrompair,bprow,max_row,bpcol,max_col)
	return(res)
}

#print(divisionMatrixDir)
#print(matrixfiles)
#matrixfiles = c("chr9_chr13_DivisionMatrix.txt")
threads = 8

matrixf <- file.path(divisionMatrixDir,matrix)
chrompair = paste(strsplit(matrix,'_')[[1]][2],strsplit(matrix,'_')[[1]][3],sep="_")
print(chrompair)
traltype <- svtype(matrixf,chrompair)
print(traltype)
#res <- searchMatric(matrixf,chrompair)
if(traltype == "unbalancedTRAL"){
	res <- searchMatric(matrixf,chrompair,threads)
}else{
	res <- searchBalancedMatric(matrixf,chrompair,threads)
}

write.table(res,file=outputfile,sep="\t",quote=F,row.names=F)

#end_time <- Sys.time()
#print(start_time)
#print(end_time)
