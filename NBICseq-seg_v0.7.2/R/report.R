#### merge the neigboring bins whose log2 copy ratio are very close

bin.merge.chr = function(segs,min_diff=0.1,adjust=0)
        {
          copy.diff = diff(segs$log2.TumorExpectRatio)
          indx = which.min(abs(copy.diff))
          seg.tmp = segs
          min.diff.tmp = c()
          k = 0
          while(nrow(seg.tmp)>1 && abs(copy.diff[indx])<=min_diff){
                min.diff.tmp = c(min.diff.tmp ,copy.diff[indx])
                seg.tmp$end[indx] = seg.tmp$end[indx+1]
		seg.tmp$binNum[indx] = seg.tmp$binNum[indx] + seg.tmp$binNum[indx+1]
                seg.tmp$tumor[indx] = seg.tmp$tumor[indx] + seg.tmp$tumor[indx+1]
		seg.tmp$tumor_expect[indx] = seg.tmp$tumor_expect[indx] + seg.tmp$tumor_expect[indx+1]
		seg.tmp$normal[indx] = seg.tmp$normal[indx] + seg.tmp$normal[indx+1]
		seg.tmp$normal_expect[indx] = seg.tmp$normal_expect[indx] + seg.tmp$normal_expect[indx+1]
                seg.tmp = seg.tmp[-(indx+1),]
		seg.tmp$log2.TumorExpectRatio = log2((seg.tmp$tumor/seg.tmp$tumor_expect)+1e-10) - adjust
		
                copy.diff = diff(seg.tmp$log2.TumorExpectRatio)
                indx = which.min(abs(copy.diff))
                k = k+1
                }
          return(seg.tmp)
        }


bin.merge = function(segs,min_diff=0.1){
	chrom.code = aggregate(1:nrow(segs),by=list(segs$chrom),min) ### make sure the chromosomes are ordered as loaded
	rk = order(chrom.code[,2])
	chrom.code = chrom.code[rk,]

        segs$tumor = as.numeric(segs$tumor)
        segs$tumor_expect = as.numeric(segs$tumor_expect)
        segs$normal = as.numeric(segs$normal)
        segs$normal_expect = as.numeric(segs$normal_expect)
        adjust = log2(sum(segs$tumor)/sum(segs$tumor_expect)+1e-10)
        segs$log2.TumorExpectRatio = log2((segs$tumor/segs$tumor_expect)+1e-10) - adjust


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
        { print("Usage: \'R --slave --args <InputData1> <InputData2> <output> < report.R\'")
          q(status = 1)
        }

infile1 = args[[1]];
infile2 = args[[2]];
outfile = args[[3]];


x1 = read.table(infile1,header=T)
x2 = read.table(infile2,header=T)



ind = rep(c(1:nrow(x2)),times=x2$binNum)

tumor = aggregate(x1$tumor,by=list(ind),sum)
tumor_expect = aggregate(x1$tumor_expected,by=list(ind),sum)


normal = aggregate(x1$normal,by=list(ind),sum)
normal_expect = aggregate(x1$normal_expected,by=list(ind),sum)

binNum = aggregate(x1$binNum,by=list(ind),sum)

start = aggregate(x1$start,by=list(ind),min)
end = aggregate(x1$end,by=list(ind),max)
chroms = aggregate(x1$chrom,by=list(ind),FUN=function(z){z[1]})

if(sum(start[,2]-x2$start)!=0||sum(end[,2]-x2$end)!=0){
	stop(paste("Error",infile1, infile2, "are not compatible"))
	}


weight.tumor = sum(as.numeric(tumor[,2]))/sum(as.numeric(tumor_expect[,2]))
weight.normal = sum(as.numeric(normal[,2]))/sum(as.numeric(normal_expect[,2])) 

log2copyratio = log2((tumor[,2]/tumor_expect[,2])/(normal[,2]/normal_expect[,2] + 1e-10)+1e-10) - log2(weight.tumor/weight.normal + 1e-10)
log2TumorExpectRatio = log2((tumor[,2]/tumor_expect[,2])/(normal[,2]/normal_expect[,2]) +1e-10) - log2(weight.tumor + 1e-10)
segs = data.frame(chrom=chroms[,2],start=start[,2],end=end[,2],binNum=binNum[,2],tumor=tumor[,2],tumor_expect = tumor_expect[,2],normal=normal[,2],normal_expect=normal_expect[,2],log2.copyRatio=log2copyratio, log2.TumorExpectRatio = log2TumorExpectRatio)
segs = bin.merge(segs)
segs$log2.copyRatio = log2((segs$tumor/segs$tumor_expect)/(segs$normal/segs$normal_expect + 1e-10)+1e-10) - log2(weight.tumor/weight.normal + 1e-10)


#meanReadPerbp_tumor = sum(segs$tumor_expect)/sum(as.numeric(segs$end-segs$start+1))
#ReadPerbp_tumor = segs$tumor_expect/(as.numeric(segs$end-segs$start+1))
#ind = which(ReadPerbp_tumor < meanReadPerbp_tumor/5)
#if(length(ind)>0) segs = segs[-ind,] ### remove regions whose mappability is not good; CNV calls in those regions are less confident and I remove them


write.table(segs,file=outfile,col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")


