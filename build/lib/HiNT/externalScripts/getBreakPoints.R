Args <- commandArgs()
print(Args)
divisionMatrixDir = Args[6]
outputfile = Args[7]

library('strucchange')
library(pastecs)

allres <- c()

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

searchMatric <- function(matrixf,chrompair){
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
	breakpoints_row <- breakpoints(rowsum ~ 1, h=5, breaks=10)$breakpoints
	breakpoints_col <- breakpoints(colsum ~ 1, h=5, breaks=10)$breakpoints
	#breakpoints_row <- paste(breakpoints_row, collapse = ',')
	#breakpoints_col <- paste(breakpoints_col, collapse = ',')
	max_row <- which(rowsum==max(rowsum))
	max_col <- which(colsum==max(colsum))
	validRowBP <- filtering(breakpoints_row,rowsum)
	validColBP <- filtering(breakpoints_col,colsum)
	bprow <- paste(validRowBP,collapse = ',')
	bpcol <- paste(validColBP,collapse = ',')
	res <- cbind(chrompair,bprow,max_row,bpcol,max_col)
	return(res)
}

#print(divisionMatrixDir)
matrixfiles = list.files(divisionMatrixDir)
#print(matrixfiles)
#matrixfiles = c("test_chr17_chr20_DivisionMatrix.txt","test_chr3_chr10_DivisionMatrix.txt","test_chr9_chr15_DivisionMatrix.txt")
for (matrix in matrixfiles){
	matrixf <- file.path(divisionMatrixDir,matrix)
	chrompair = paste(strsplit(matrix,'_')[[1]][2],strsplit(matrix,'_')[[1]][3],sep="_")
	print(chrompair)
	res <- searchMatric(matrixf,chrompair)
	allres <- rbind(allres,res)
}

#print(allres)
#print(outputfile)
write.table(allres,file=outputfile,sep="\t",quote=F,row.names=F)
