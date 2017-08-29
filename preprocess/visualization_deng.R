rm(list=ls())
setwd('D:\\cell\\data\\visualization')
load("D:/cell/data/processed/processed.deng.read.rdata")


###########################################################
##pca######################################################
###########################################################

expression <- t(as.matrix(processed.data[, 3:55]))
expression <- log(expression+1)
X <- t(scale(expression))
singular <- svd(X)
pc1d <- singular$u[,1]
pc2d <- singular$u[,2]
pc1 <- expression%*%pc1d
pc2 <- expression%*%pc2d
pc1.16.88 <- pc1[1:27]
pc1.16.193 <- pc1[28:35]
pc1.8.88 <- pc1[36:44]
pc1.8.193 <- pc1[45:53]
pc2.16.88 <- pc2[1:27]
pc2.16.193 <- pc2[28:35]
pc2.8.88 <- pc2[36:44]
pc2.8.193 <- pc2[45:53]
library(ggplot2)
p <- ggplot()
p <- p + geom_point(data=data.frame(pc1.16.88, pc2.16.88), aes(x=pc1.16.88, y=pc2.16.88, color="16cellrun0088"), size=3)
p <- p + geom_point(data=data.frame(pc1.16.193, pc2.16.193), aes(x=pc1.16.193, y=pc2.16.193, color="16cellrun00193"), size=3)
p <- p + geom_point(data=data.frame(pc1.8.88, pc2.8.88), aes(x=pc1.8.88, y=pc2.8.88, color="8cellrun0088"), size=3)
p <- p + geom_point(data=data.frame(pc1.8.193, pc2.8.193), aes(x=pc1.8.193, y=pc2.8.193, color="8cellrun00193"), size=3)
p <- p + ggtitle('PC1 and PC2 for 16cell and 8cell in two batches')


########################################################################################################################################

fenmu <- nrow(processed.data)
detected.rate <- apply(processed.data[,3:55], 2, f <- function(x)sum(x!=0)/fenmu)
cell <- c(rep('16cell', 35), rep('8cell', 18))
batch <- c(rep('run0088', 27), rep('run00193', 8), rep('run0088', 9), rep('run00193', 9))
detec.frame <- data.frame(cell,batch,detected.rate)
p<-ggplot(data=detec.frame, aes(x=cell,y=detected.rate))+geom_boxplot(aes(fill=batch))

########################################################################################################################################

pc1.frame <- data.frame(cell,batch,pc1)
p<-ggplot(data=pc1.frame, aes(x=cell,y=pc1))+geom_boxplot(aes(fill=batch))

########################################################################################################################################

p <- ggplot()
p <- p + geom_point(data=data.frame(detected.rate[1:27], pc1.16.88), aes(x=detected.rate[1:27], y=pc1.16.88, color='16cellrun0088'))
p <- p + geom_point(data=data.frame(detected.rate[28:35], pc1.16.193), aes(x=detected.rate[28:35], y=pc1.16.193, color='16cellrun00193'))
p <- p + geom_point(data=data.frame(detected.rate[36:44], pc1.8.88), aes(x=detected.rate[36:44], y=pc1.8.88, color='8cellrun0088'))
p <- p + geom_point(data=data.frame(detected.rate[45:53], pc1.8.193), aes(x=detected.rate[45:53], y=pc1.8.193, color='8cellrun00193'))
