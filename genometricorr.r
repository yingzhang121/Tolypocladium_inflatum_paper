library(GenometriCorr)

source("tools.r")
s="31671 31975 567_84 714_70 8044 824_70"
#sp=scan(textConnection(s),what='')
#sp=c("8044", "31671", "31975")
seqtable=read.table("seqinfo")

for (s in sp) {
    gcf <- paste( s, "genecluster", sep=".")
    rpf <- paste( s, "min.out", sep=".")
    gc <- read.table(gcf, header=F)
    rp <- read.table(rpf, header=F, fill=T)
    subs = seqtable[ seqtable$V4==s, ]
    seqinfo <- Seqinfo( as.vector(subs[,1]), subs[,2], subs[,3], subs[,4] )
    xgr = GRanges(seqnames=gc$V1, ranges=IRanges(start=gc$V2, end=gc$V3), values=gc$V4, seqinfo=seqinfo)
    ygr = GRanges(seqnames=rp$V1, ranges=IRanges(start=rp$V2, end=rp$V3), values=rp$V5, seqinfo=seqinfo)
    res <- analysis_gc( xgr, ygr )       
    fnm <- paste(s,"result", sep=".")
    sink(fnm)
    print(res)
    sink()
}

sessionInfo()
