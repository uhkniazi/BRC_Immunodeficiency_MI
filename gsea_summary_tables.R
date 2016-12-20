# Name: gsea_summary_tables.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 20/12/2016
# Desc: merge the gsea tables into one


lFiles = list.files('Results/', pattern='*mSigDb_c2.xls', full.names = T, ignore.case = T)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = c('D_UN_v_D_S', 'H_S_v_D_S', 'H_UN_v_D_UN', 'H_UN_v_H_S') 
names(ldfData.up) = sn

sn = c('D_UN_v_D_S', 'H_S_v_D_S', 'H_UN_v_D_UN', 'H_UN_v_H_S') 
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
cn = sn
colnames(mMerged.up) = paste(cn, 'up', sep='-')
colnames(mMerged.down) = paste(cn, 'down', sep='-')

mMerged.c2 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c2)
o = c(1, 5, 2, 6, 3, 7, 4, 8)
mMerged.c2 = mMerged.c2[,o]

# remove na sections
dim(mMerged.c2)
mMerged.c2 = na.omit(mMerged.c2)
dim(mMerged.c2)
head(mMerged.c2)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c2.bin = getBinaryMatrix(mMerged.c2)

## group this matrix into combinations
mMerged.c2.bin.grp = mMerged.c2.bin
set.seed(123)
dm = dist(mMerged.c2.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c2.bin.grp = cbind(mMerged.c2.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c2.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c2.bin)
dfMerged.c2 = data.frame(round(mMerged.c2, 3), sig.pvals, groups, DB='C2')
str(dfMerged.c2)
head(dfMerged.c2)
tail(dfMerged.c2)

################### repeat for c5 database

lFiles = list.files('Results/', pattern='*mSigDb_c5.xls', full.names = T, ignore.case = T)

# load the files
ldfData = lapply(lFiles, function(x) read.csv(x, row.names=1))
names(ldfData) = lFiles
sapply(ldfData, nrow)
# put everything in one order by row names
rn = rownames(ldfData[[1]])
head(rn)

ldfData = lapply(ldfData, function(df){
  df = df[rn,]
})

# get the upregulated/downregulated in order
ldfData.up = ldfData[grepl(lFiles, pattern = 'upregulated')]
ldfData.down = ldfData[grepl(lFiles, pattern = 'downregulated')]

## set the names for each contrast
sn = c('D_UN_v_D_S', 'H_S_v_D_S', 'H_UN_v_D_UN', 'H_UN_v_H_S') 
names(ldfData.up) = sn

sn = c('D_UN_v_D_S', 'H_S_v_D_S', 'H_UN_v_D_UN', 'H_UN_v_H_S') 
names(ldfData.down) = sn

## create a table/matrix of p-values
mMerged.up = sapply(ldfData.up, function(x) x$p.val)
rownames(mMerged.up) = rownames(ldfData.up[[1]])
## create similar table for downregulated
mMerged.down = sapply(ldfData.down, function(x) x$p.val)
rownames(mMerged.down) = rownames(ldfData.down[[1]])

# sanity check
identical(rownames(mMerged.up), rownames(mMerged.down))
cn = sn
colnames(mMerged.up) = paste(cn, 'up', sep='-')
colnames(mMerged.down) = paste(cn, 'down', sep='-')

mMerged.c5 = cbind(mMerged.up, mMerged.down)
# reorder the columns
colnames(mMerged.c5)
o = c(1, 5, 2, 6, 3, 7, 4, 8)
mMerged.c5 = mMerged.c5[,o]

# remove na sections
dim(mMerged.c5)
mMerged.c5 = na.omit(mMerged.c5)
dim(mMerged.c5)
head(mMerged.c5)

### create a binary matrix based on cutoffs
getBinaryMatrix = function(mat, cutoff=0.01){
  mat2 = apply(mat, 2, function(x) round(x, 3) <= cutoff)
}

mMerged.c5.bin = getBinaryMatrix(mMerged.c5)

## group this matrix into combinations
mMerged.c5.bin.grp = mMerged.c5.bin
set.seed(123)
dm = dist(mMerged.c5.bin.grp, method='binary')
hc = hclust(dm)
plot(hc, labels=F)
# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mMerged.c5.bin.grp = cbind(mMerged.c5.bin.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mMerged.c5.bin.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

## map these names to the cp
groups = cp
sig.pvals = rowSums(mMerged.c5.bin)
dfMerged.c5 = data.frame(round(mMerged.c5, 3), sig.pvals, groups, DB='C5')
str(dfMerged.c5)
head(dfMerged.c5)
tail(dfMerged.c5)

# sanity check
identical(colnames(dfMerged.c2), colnames(dfMerged.c5))

write.csv(dfMerged.c2, file='Results/gsea_msigdb_c2_merged.xls')
write.csv(dfMerged.c5, file='Results/gsea_msigdb_c5_merged.xls')

## merge together into one dataframe
# remove the group that has no significant p-values
table(dfMerged.c2$groups)
dfMerged.c2.sub = dfMerged.c2[dfMerged.c2$groups != 9,]

table(dfMerged.c5$groups)
dfMerged.c5.sub = dfMerged.c5[dfMerged.c5$groups != 5,]

dfMerged = rbind(dfMerged.c2.sub, dfMerged.c5.sub)
dfMerged = droplevels.data.frame(dfMerged)

write.csv(dfMerged, file='Results/gsea_msigdb_significant_c2_c5_merged.csv')

### heatmaps
df = dfMerged
head(df)
mMat = as.matrix(df[,c(1:8)])
head(mMat)
mMat = -10*log10(mMat+1e-16)
g1 = df[,'groups']
g1 = factor(as.character(g1))
levels(g1)
g2 = df[,'DB']
g2 = factor(as.character(g2))
levels(g2)

ann = data.frame(DB=g2, Group=g1 )
range(mMat)
# threshold the values for plotting
mMat[mMat < 15] = 0 
mMat[mMat > 100] = 100

library(NMF)
library(RColorBrewer)
aheatmap(mMat, annRow = ann, scale = 'none', Rowv = order(g2:g1), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))

pdf('Results/gsea_msigdb_significant_merged.pdf')
aheatmap(mMat, annRow = ann, scale = 'none', Rowv = order(g2:g1), Colv=NA, cexRow=5, cexCol = 0.6, #labCol=c('C2vC1', 'K1vC1', 'K2vC2', 'K2vK1'), 
         col=c('white', brewer.pal(9, 'YlOrRd')))
dev.off(dev.cur())