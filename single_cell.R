#Read the data, we use deng's data as a test of our method.
data_all <- as.matrix(new.processed.data)
Y <- as.matrix(data_all[,c(3:56)])
G_num <- length(Y[,1])
C_num <- length(Y[1,])
Y <- matrix(as.numeric(Y), G_num, C_num)
gene_length <- as.numeric(as.vector(data_all[,57])) * 10

tech_para <- list(iternum = 100, error = 1e-5, mhnum = 300, jump_rate = 0.15,
                  jump = 0.03, burnin = 50)
source('fitDABB_function.R')
#fit the model:
#Y is the G*C read count matrix, total.count is a vector of each cell's total read counts,
#batch.info is a vector, the index of the batch each cell belongs to, bio.group is a 
#vector, the index of the biological group each cell belongs to, gene_length is a vector,
#length of each RNA, and tech_para is a group of parameters that are used to fit the model.
#Among these parameters, jump_rate and jump are used to control the acceptance rate of 
#Metropolis-Hasting alogrithm. A recommanded rate ranges from 0.2 to 0.45.
results <- fitDABB(Y, total.count, batch.info, bio.group, gene_length, tech_para)
results$coef
#Qualtiy Control, alternative hypothesis can be chosen as 'left', 'right' and 'two side'.
DABB_QC(results, alternative = 'right')
#refit after deleting outliers
results <- fitDABB(Y[,-c(22, 33, 35)], total.count[-c(22, 33, 35)], batch.info[-c(22, 33, 35)], bio.group[-c(22, 33, 35)], 
                   gene_length, tech_para)

#Differential Expression
pvl <- DABB_DE(Y[,-c(22, 33, 35)], total.count[-c(22, 33, 35)], bio.group[-c(22, 33, 35)], gene_len = gene_length, results, 
               sample_num = 200)
sort(pvl$p.value, index.return = T)
sort(p.adjust(pvl$p.value))

#Visualization, method can be chosen as 'PCA' and 'ISOmap'
library(ggplot2)
library(vegan)

visual <- DABB_visualize(Y, total.count, gene_length, results$pweight,
                         results$bsample, method = 'PCA', k = 10)
##visual <- visual$points
pc1 <- visual[,1]
pc2 <- visual[,2]
pc1.16.88 <- pc1[1:27]
pc1.16.193 <- pc1[28:35]
pc1.8.88 <- pc1[36:45]
pc1.8.193 <- pc1[46:54]
pc2.16.88 <- pc2[1:27]
pc2.16.193 <- pc2[28:35]
pc2.8.88 <- pc2[36:45]
pc2.8.193 <- pc2[46:54]

p <- ggplot()
p <- p + geom_point(data=data.frame(pc1.16.88, pc2.16.88), aes(x=pc1.16.88, y=pc2.16.88, color="16cell", shape='run0088'), size=3)
p <- p + geom_point(data=data.frame(pc1.16.193, pc2.16.193), aes(x=pc1.16.193, y=pc2.16.193, color="16cell",shape='run00193'), size=3)
p <- p + geom_point(data=data.frame(pc1.8.88, pc2.8.88), aes(x=pc1.8.88, y=pc2.8.88, color="8cell",shape='run0088'), size=3)
p <- p + geom_point(data=data.frame(pc1.8.193, pc2.8.193), aes(x=pc1.8.193, y=pc2.8.193, color="8cell",shape='run00193'), size=3)
p <- p + ggtitle('PC1 and PC2 for 16cell and 8cell in two batches')

####################################################################
#Visualization of the unfitted data.
G_num <- length(Y[,1])
C_num <- length(Y[1,])
sample_num <- length(results$bsample[1,])
visual <- DABB_visualize(Y, total.count, gene_length, matrix(1, G_num, C_num),
                         b_sample_mat = matrix(0, C_num, sample_num), method = 'PCA', k = 10)
##visual <- visual$points
pc1 <- visual[,1]
pc2 <- visual[,2]
pc1.16.88 <- pc1[1:27]
pc1.16.193 <- pc1[28:35]
pc1.8.88 <- pc1[36:45]
pc1.8.193 <- pc1[46:54]
pc2.16.88 <- pc2[1:27]
pc2.16.193 <- pc2[28:35]
pc2.8.88 <- pc2[36:45]
pc2.8.193 <- pc2[46:54]

p <- ggplot()
p <- p + geom_point(data=data.frame(pc1.16.88, pc2.16.88), aes(x=pc1.16.88, y=pc2.16.88, color="16cell", shape='run0088'), size=3)
p <- p + geom_point(data=data.frame(pc1.16.193, pc2.16.193), aes(x=pc1.16.193, y=pc2.16.193, color="16cell",shape='run00193'), size=3)
p <- p + geom_point(data=data.frame(pc1.8.88, pc2.8.88), aes(x=pc1.8.88, y=pc2.8.88, color="8cell",shape='run0088'), size=3)
p <- p + geom_point(data=data.frame(pc1.8.193, pc2.8.193), aes(x=pc1.8.193, y=pc2.8.193, color="8cell",shape='run00193'), size=3)
p <- p + ggtitle('PC1 and PC2 for 16cell and 8cell in two batches')

