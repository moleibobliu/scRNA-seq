## install package
source("https://bioconductor.org/biocLite.R")
biocLite("scater")
library(scater, quietly = TRUE)

## example from the package (do not run)
data("sc_example_counts")
data("sc_example_cell_info")
summary(sc_example_cell_info)
pd <- new("AnnotatedDataFrame", data = sc_example_cell_info)
rownames(pd) <- pd$Cell
example_sceset <- newSCESet(countData = sc_example_counts, phenoData = pd)

keep_feature <- rowSums(exprs(example_sceset) > 0) > 0
example_sceset <- example_sceset[keep_feature,]

example_sceset <- calculateQCMetrics(example_sceset, feature_controls = 1:40)
scater_gui(example_sceset)

## Zeisel data
# create secset data
data <- rbind.data.frame(data.CA1Pyr1, data.CA1Pyr2)
data1 <- t(apply(data, c(1,2), as.numeric))
subtype <- c(rep('Pyr1',380), rep('Pyr2',447))
metadata <- data.frame(subtype)
row.names(metadata) <- row.names(data)
pd <- new("AnnotatedDataFrame", data = metadata)
sceset <- newSCESet(countData = data1, phenoData = pd)
# do quality control
qc_sceset <- calculateQCMetrics(sceset)
scater_gui(qc_sceset)
# no cell is dected as outlier

## note:
# 1, input countdata should be of p*n form
# 2, transform data into numeric type