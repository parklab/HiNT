#### plot functions
highlightChrom <- function(adjustments, min, max){
  for(index in 1:(length(adjustments)-1)){
    if(index %% 2 == 1){
      polygon(c(adjustments[index], adjustments[index + 1],
                adjustments[index + 1], adjustments[index]),
              c(min, min, max, max), col = "gray", border = "white")
    }
  }
  return(invisible())
}



markChrom <- function(adjustments, min,chromNames,height=10){
  chromLocs <- NULL
  ind = NULL
  for(i in 1:(length(adjustments) - 1)){
    if(i %% 2 == 1){
      chromLocs <- c(chromLocs, mean(c(adjustments[i], adjustments[i + 1])))
      ind = c(ind,i)
    }
  }
  if(length(chromLocs)!=length(chromNames[ind])){stop("Incorrect number of chromsome names")}
  text(chromLocs, rep(min - 0.125, length(chromLocs)), chromNames[ind], cex = 1.2)
  return(invisible())
}


markboundary <- function(adjustments,min,max){
    if(length(adjustments)<3) return;
    for(i in c(2:(length(adjustments-1)))){
        lines(rep(adjustments[i],2),c(min,max))
        }
}




args <- commandArgs(TRUE)
if(length(args)!=3 && length(args)!=2)
	{ print("       \'R --slave --args <SegFile> <FigName> < report.R\'")
	  print("or")
	  print("       \'R --slave --args <SegFile> <FigName> <description> < report.R\'")
          q(status = 1)
        }

segfile = args[[1]];
segs = read.table(segfile,header=T)


description = ""
figname = args[[2]]
if(length(args)==3) description = args[[3]]


png(filename=figname,width=1200,height=500)
chrom.code = aggregate(1:nrow(segs),by=list(segs$chrom),min) ### make sure the chromosomes are ordered as loaded
rk = order(chrom.code[,2])
chrom.length = aggregate(as.numeric(segs$end),by=list(segs$chrom),max)
chrom.length = chrom.length[rk,]
num.segs.chrom = aggregate(segs$end,by=list(segs$chrom),length)
num.segs.chrom = num.segs.chrom[rk,]
chrom.boundary = c(0,cumsum(chrom.length[,2]))
plot(0, 0, type = "n", main = description, xlab = "Chromosome",
ylab = " Log2 ratio", ylim = c(-5, 5), axes = FALSE,
xlim = c(0, max(chrom.boundary) + 10))

highlightChrom(chrom.boundary,-5,5)
lines(c(min(chrom.boundary), max(chrom.boundary)), rep(0, 2), lwd = 1, col = "blue")
lines(c(min(chrom.boundary), max(chrom.boundary)), rep(-1, 2), lwd = 1, col = "blue")
lines(c(min(chrom.boundary), max(chrom.boundary)), rep(1, 2), lwd = 1, col = "blue")
axis(2)
box()

adjust = rep(chrom.boundary[-length(chrom.boundary)],num.segs.chrom[,2])
loc.start = segs$start + adjust
loc.end = segs$end + adjust
for(k in c(1:length(loc.start))){
	if(!is.null(segs$pvalue)){ 
		if(segs$pvalue[k]<0.01 && abs(segs$log2.copyRatio[k])>0.2) col="red" else col = "green"
		}else{
		if(abs(segs$log2.copyRatio[k])>0.2)  col="red" else col = "green"
		}
	lines(c(loc.start[k],loc.end[k]),rep(segs$log2.copyRatio[k],2),col=col,lwd=2)
	}
markChrom(chrom.boundary,-5,chrom.length[,1])
dev.off()


