
#### merge the neigboring bins whose log2 copy ratio are very close

bin.merge.chr = function(segs,min_diff=0.1,adjust=0)
        {
          copy.diff = diff(segs$log2.copyRatio)
          indx = which.min(abs(copy.diff))
          seg.tmp = segs
          min.diff.tmp = c()
          k = 0
          while(nrow(seg.tmp)>1 && abs(copy.diff[indx])<=min_diff){
                min.diff.tmp = c(min.diff.tmp ,copy.diff[indx])
                seg.tmp$end[indx] = seg.tmp$end[indx+1]
		seg.tmp$binNum[indx] = seg.tmp$binNum[indx] + seg.tmp$binNum[indx+1]
                seg.tmp$observed[indx] = seg.tmp$observed[indx] + seg.tmp$observed[indx+1]
		seg.tmp$expected[indx] = seg.tmp$expected[indx] + seg.tmp$expected[indx+1]
		if(ncol(seg.tmp)>7) seg.tmp[indx,-c(1:7)] = seg.tmp[indx,-c(1:7)] + seg.tmp[indx+1,-c(1:7)]
                seg.tmp = seg.tmp[-(indx+1),]
		seg.tmp$log2.copyRatio = log2((seg.tmp$observed/seg.tmp$expected)+1e-10) - adjust
		
                copy.diff = diff(seg.tmp$log2.copyRatio)
                indx = which.min(abs(copy.diff))
                k = k+1
                }
          return(seg.tmp)
        }


bin.merge = function(segs,min_diff=0.1){
	chrom.code = aggregate(1:nrow(segs),by=list(segs$chrom),min) ### make sure the chromosomes are ordered as loaded
	rk = order(chrom.code[,2])
	chrom.code = chrom.code[rk,]

	for(j in c(5:ncol(segs))){
	        segs[,j] = as.numeric(segs[,j])
		}
        cutoff = 0.3
        indx = which(abs(segs$log2.copyRatio)<=cutoff)
        adjust = 0
        if(sum(segs$observed[indx])>=1000){
               adjust = log2(sum(segs$observed[indx])/sum(segs$expected[indx]))
               }
        segs$log2.copyRatio = log2(segs$observed/segs$expected+1e-10) - adjust




	chroms = as.character(segs$chrom)
	segs.tmp = NULL
	for(k in c(1:nrow(chrom.code))){
		chr = as.character(chrom.code[k,1])
		indx = which(chroms==chr)
		segs.chr = segs[indx,]
		segs.chr.tmp = bin.merge.chr(segs.chr,min_diff,adjust=adjust)
		segs.tmp = rbind(segs.tmp,segs.chr.tmp)
		}

	return(segs.tmp)
	}



args <- commandArgs(TRUE)
if(length(args)!=3)
        { print("Usage: \'R --slave --args <InputData> <RemoveFDs(Y/N)> <output> < reportOneSample.R\'")
          q(status = 1)
        }

infile1 = args[[1]];
removeFD = args[[2]];
outfile = args[[3]];

if(removeFD!="Y" && removeFD!="N"){
	print("The second argument must be 'Y' or 'N'");
	q(status=1)
	}

x1 = read.table(infile1,header=T)


num_Sample = (ncol(x1) - 4)/2;
if(num_Sample<1 || num_Sample!=floor(num_Sample)){
	print("reportOneSample.R: The input file must have 4+2n columns");
	q(status = 1)
	}

if(num_Sample==1){
	#log2copyratio = log2((tumor[,2]/tumor_expect[,2])/(normal[,2]/normal_expect[,2])) - log2(weight.tumor/weight.normal)
	log2copyratio = log2((x1$observed+1e-10)/(x1$expected + 1e-10)) - log2(sum(x1$observed)/sum(x1$expected)+1e-10)
	segs = data.frame(x1,log2.copyRatio=log2copyratio)

	if(removeFD=="Y"){
		meanReadPerbp = sum(segs$expected)/sum(as.numeric(segs$end-segs$start+1))
		ReadPerbp = segs$expected/(as.numeric(segs$end-segs$start+1))
		ind = which(ReadPerbp < meanReadPerbp/5)
		if(length(ind)>0) segs = segs[-ind,]
		}

	segs = bin.merge(segs)
	write.table(segs,file=outfile,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
	}else{
	ind.obs = 4 + c(1:num_Sample)*2-1
	ind.exp = 4 + c(1:num_Sample)*2
	cnt.obs = rowSums(x1[,ind.obs])
	cnt.exp = rowSums(x1[,ind.exp])

	log2copyratio = log2((cnt.obs+1e-10)/(cnt.exp+1e-10)) - log2(sum(cnt.obs)/sum(cnt.exp)+1e-10)
	segs = data.frame(x1[,1:4],observed=cnt.obs,expected=cnt.exp,log2.copyRatio=log2copyratio,x1[,-c(1:4)])
	if(removeFD=="Y"){
		meanReadPerbp = sum(segs$expected)/sum(as.numeric(segs$end-segs$start+1))
		ReadPerbp = segs$expected/(as.numeric(segs$end-segs$start+1))
		ind = which(ReadPerbp < meanReadPerbp/5)
		if(length(ind)>0) segs = segs[-ind,]
		}
	segs = bin.merge(segs)
	write.table(segs,file=outfile,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
	}

