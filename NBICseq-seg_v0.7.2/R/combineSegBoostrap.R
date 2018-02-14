
args <- commandArgs(TRUE)
if(length(args)!=2)
        { print("Usage: \'R --slave --args <SegFile> <BootStrapFile> < report.R\'")
          q(status = 1)
        }

segfile = args[[1]];
bootfile = args[[2]];

segs = read.table(segfile,header=T)
bootdata = read.table(bootfile,header=T)

if(nrow(segs)!=nrow(bootdata)||sum(segs$start!=bootdata$start)!=0 || sum(segs$end!=bootdata$end)!=0){
	print("Error <SegFile> and <BootStrapFile are not compatible>");
	q(status=1)
	}

#segs$zscore = bootdata$zscore
segs$pvalue = bootdata$pvalue
segs$pvalue.TumorVsExpected = bootdata$pvalue.TumorVsExpected

write.table(segs,file=segfile,row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")


