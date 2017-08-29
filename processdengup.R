rm(list=ls())
setwd('D:\\cell\\data\\processed')
load("D:/cell/data/processed/raw.deng.read.rdata")
load("D:/cell/data/processed/raw.RPKM.deng.read.rdata")
bio.group <- c(rep(1,35), rep(2,19))
batch.info <- c(rep(1,27),rep(2,8),rep(1,10),rep(2,9))
total.count <- apply(count.data[,3:56],2,sum)
lengthm <- matrix(nrow=22958, ncol=0)
for(i in 1:54){
  leng <- count.data[,i+2]/(fpkm.data[,i+2]*total.count[i])
  lengthm <- cbind(lengthm, leng)
}

exon.length <- apply(lengthm,1,FUN= function(x) median(x,na.rm=T))
exon.length <- exon.length*10^8
count.data$gene.length <- exon.length
rm(list=c('i', 'leng', 'fpkm.data', 'lengthm', 'exon.length'))
##########################################################################################################
mean.data <- apply(count.data , 1, FUN = function(x) mean(as.integer(x[3:56])))
threshold <- quantile(mean.data, 0.95)
index.processed <- mean.data < threshold
processed.data <- count.data[index.processed,]
##########################################################################################################
amount <- apply(processed.data,1, FUN=function(x) sum(as.integer(x[3:56])!=0))
index.p <- amount > 10
processed.data <- processed.data[index.p,]
rm(list=c('amount', 'index.p', 'index.processed', 'mean.data', 'threshold'))
###############################
index.pr <- processed.data$gene.length != Inf & processed.data$gene.length < 5000
processed.data <- processed.data[index.pr,]
rm(index.pr)
###############################
save(list=ls(), file='deng.data.rdata')
