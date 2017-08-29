rm(list=ls())
setwd('D:\\cell\\data\\processed')
sixteenrun88 <- dir('D:\\cell\\data\\deng\\possion data\\16cellrun0088')
eightrun88 <- dir('D:\\cell\\data\\deng\\possion data\\8cellrun0088')
sixteenrun193 <- dir('D:\\cell\\data\\deng\\possion data\\16cellrun00193')
eightrun193 <- dir('D:\\cell\\data\\deng\\possion data\\8cellrun00193')
example <- read.csv('D:\\cell\\data\\deng\\possion data\\16cellrun0088\\GSM1112490_16cell_1-10_expression.txt',sep='\t')

fpkm.data <- data.frame(X.Gene_symbol = example$X.Gene_symbol, Reseq_IDs = example$Refseq_IDs)
i <- 1
for(filename in sixteenrun88){
  raw.data <- read.csv(paste0('D:\\cell\\data\\deng\\possion data\\16cellrun0088\\', filename), sep='\t')
  fpkm.data = cbind(fpkm.data, raw.data$RPKM)
}
for(filename in sixteenrun193){
  raw.data <- read.csv(paste0('D:\\cell\\data\\deng\\possion data\\16cellrun00193\\', filename), sep='\t')
  fpkm.data = cbind(fpkm.data, raw.data$RPKM)
}
for(filename in eightrun88){
  raw.data <- read.csv(paste0('D:\\cell\\data\\deng\\possion data\\8cellrun0088\\', filename), sep='\t')
  fpkm.data = cbind(fpkm.data, raw.data$RPKM)
}
for(filename in eightrun193){
  raw.data <- read.csv(paste0('D:\\cell\\data\\deng\\possion data\\8cellrun00193\\', filename), sep='\t')
  fpkm.data = cbind(fpkm.data, raw.data$RPKM)
}

#############################################################################################################

colnames(fpkm.data)[3:29] = sapply(1:27, FUN = function(x) paste0('16cell-0088-', x))
colnames(fpkm.data)[30:37] = sapply(1:8, FUN = function(x) paste0('16cell-00193-', x))
colnames(fpkm.data)[38:47] = sapply(1:10, FUN = function(x) paste0('8cell-0088-', x))
colnames(fpkm.data)[48:56] = sapply(1:9, FUN = function(x) paste0('8cell-00193-', x))

save(fpkm.data, file='raw.RPKM.read.rdata')

fpkm.data <- data.frame(X.Gene_symbol = example$X.Gene_symbol, Reseq_IDs = example$Refseq_IDs)
i <- 1
for(filename in sixteenrun88){
  raw.data <- read.csv(paste0('D:\\cell\\data\\deng\\possion data\\16cellrun0088\\', filename), sep='\t')
  fpkm.data = cbind(fpkm.data, raw.data$reads)
}
for(filename in sixteenrun193){
  raw.data <- read.csv(paste0('D:\\cell\\data\\deng\\possion data\\16cellrun00193\\', filename), sep='\t')
  fpkm.data = cbind(fpkm.data, raw.data$reads)
}
for(filename in eightrun88){
  raw.data <- read.csv(paste0('D:\\cell\\data\\deng\\possion data\\8cellrun0088\\', filename), sep='\t')
  fpkm.data = cbind(fpkm.data, raw.data$reads)
}
for(filename in eightrun193){
  raw.data <- read.csv(paste0('D:\\cell\\data\\deng\\possion data\\8cellrun00193\\', filename), sep='\t')
  fpkm.data = cbind(fpkm.data, raw.data$reads)
}


colnames(fpkm.data)[3:29] = sapply(1:27, FUN = function(x) paste0('16cell-0088-', x))
colnames(fpkm.data)[30:37] = sapply(1:8, FUN = function(x) paste0('16cell-00193-', x))
colnames(fpkm.data)[38:47] = sapply(1:10, FUN = function(x) paste0('8cell-0088-', x))
colnames(fpkm.data)[48:56] = sapply(1:9, FUN = function(x) paste0('8cell-00193-', x))

save(fpkm.data, file='raw.read.rdata')

####################################################################################################
load('D:\\cell\\data\\processed\\raw.deng.read.rdata')

processed.data <- fpkm.data[,-43]
processed.data <- processed.data[apply(processed.data, 1, FUN = function(x) any(as.integer(x[3:55]) != 0)),]
mean.data <- apply(processed.data , 1, FUN = function(x) mean(as.integer(x[3:55])))
threshold <- quantile(mean.data, 0.95)
index.processed <- mean.data < threshold
processed.data <- processed.data[index.processed,]
save(processed.data, file='processed.deng.read.rdata')

####################################################################################################

load("D:/cell/data/processed/processed.read.rdata")
amount <- apply(processed.data,1, FUN=function(x) sum(as.integer(x[3:56])!=0))
index.p <- amount > 40

###################################################################################################
rm(list=ls())
setwd('D:\\cell\\data\\processed')
load("D:/cell/data/processed/raw.deng.read.rdata")
load("D:/cell/data/processed/raw.RPKM.deng.read.rdata")
total.count <- apply(count.data[,3:56],2, sum)
