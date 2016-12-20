# Name: import.R
# Auth: umar.niazi@kcl.ac.uk
# Date: 19/12/2016
# Desc: import and examine data with some eda

# 4 contrasts

## HEALTHY_UNSTIMUTATION and HEALTHY_STIMULATION
dfData.1 = read.csv('Data_external/HEALTHY_UNSTIMUTATION_v_HEALTHY_STIMULATION.diff', sep='\t', header=T)
dim(dfData.1)

table(dfData.1$significant)

## HEALTHY_UNSTIMUTATION and DISEASED_UNSTIMULATION
dfData.2 = read.csv('Data_external/HEALTHY_UNSTIMUTATION_v_DISEASED_UNSTIMULATION.diff', sep='\t', header=T)
dim(dfData.2)

table(dfData.2$significant)

## HEALTHY_STIMULATION and DISEASED_STIMULATION
dfData.3 = read.csv('Data_external/HEALTHY_STIMULATION_v_DISEASED_STIMULATION.diff', sep='\t', header=T)
dim(dfData.3)

table(dfData.3$significant)

## DISEASED_UNSTIMULATION and DISEASED_STIMULATION
dfData.4 = read.csv('Data_external/DISEASED_UNSTIMULATION_v_DISEASED_STIMULATION.diff', sep='\t', header=T)
dim(dfData.4)

table(dfData.4$significant)

lData = list('H_US_vs_H_S'=dfData.1, 'H_US_vs_D_US'=dfData.2, 'H_S_vs_D_S'=dfData.3, 'D_US_vs_D_S'=dfData.4)

par(mfrow=c(2,2))
sapply(lData, function(df) hist(df$p_value))
sapply(lData, function(df) hist(df$q_value))
sapply(lData, function(df) hist(df$log2.fold_change))
sapply(lData, function(df) {
  df = df[is.finite(df$log2.fold_change) & df$status == 'OK',]
  p.val = -1 * log10(df$p_value)
  fc = df$log2.fold_change
  fc[fc > 5] = 5
  fc[fc < -5] = -5
  col = rep('lightgrey', times=length(p.val))
  c = which(df$significant == 'yes')
  col[c] = 'red'
  plot(fc, jitter(p.val), pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col)
})

## perform gsea analysis
library(gage)
library(org.Hs.eg.db)
oMsigGS = readList('~/Data/MetaData/msigdb_c2.all.v5.2.entrez.gmt')
## or load a different c5 database
oMsigGS = readList('~/Data/MetaData/msigdb_c5.all.v5.2.entrez.gmt')

dfContrast = lData$D_US_vs_D_S
# get enterez names of the genes
dfEnterez = AnnotationDbi::select(org.Hs.eg.db, keys = as.character(dfContrast$test_id), columns = c('ENTREZID'), keytype = 'ENSEMBL')
# remove duplicated ids
table(duplicated(dfEnterez$ENSEMBL))
dfEnterez = dfEnterez[!duplicated(dfEnterez$ENSEMBL),]
## about half the ensembl trainscripts do not map to the entrez ids
## so drop those from analysis
i = is.na(dfEnterez$ENTREZID)
table(i)
dfEnterez = dfEnterez[!i,]
# match these ids with the ids in the contrast table
i = match(as.character(dfContrast$test_id), dfEnterez$ENSEMBL)
dfContrast$ENTREZID = as.character(dfEnterez$ENTREZID[i])
# remove NA rows 
i = is.na(dfContrast$ENTREZID)
table(i)
dfContrast = dfContrast[!i,]
# how many ids are duplicated?
table(duplicated(dfContrast$ENTREZID))
# drop the genes where no tests are performed
dfContrast = dfContrast[dfContrast$status == 'OK',]
# only 38 or so genes are now duplicated
table(duplicated(dfContrast$ENTREZID))
# drop these genes
dfContrast = dfContrast[!duplicated(dfContrast$ENTREZID),]
dim(dfContrast)
# drop genes with infinite fold change
table(is.finite(dfContrast$log2.fold_change.))
dfContrast = dfContrast[is.finite(dfContrast$log2.fold_change.),]
dim(dfContrast)
rownames(dfContrast) = dfContrast$ENTREZID
# get the fold changes to perform gsea
head(dfContrast)
iContFc = dfContrast$log2.fold_change.
names(iContFc) = rownames(dfContrast)
head(iContFc)
hist(iContFc)
# threshold the contrasts
iContFc[iContFc > 5] = 5
iContFc[iContFc < -5] = -5
hist(iContFc)

oGage = gage(iContFc, oMsigGS)

dfGreater = data.frame(oGage$greater)
str(dfGreater)
i = which(dfGreater$p.val < 0.01)
rownames(dfGreater[i,])

dfLess = data.frame(oGage$less)
str(dfLess)
i = which(dfLess$p.val < 0.01)
rownames(dfLess[i,])

# write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file='Results/Keloid:2vsControl:2_upregulated_pathways_mSigDb_c2_curated.xls')
# write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file='Results/Keloid:2vsControl:2_downregulated_pathways_mSigDb_c2_curated.xls')

write.csv(dfGreater[,c('p.val', 'q.val', 'set.size')], file='Results/DISEASED_UNSTIMULATION_v_DISEASED_STIMULATION_upregulated_msigdb_c5.xls')
write.csv(dfLess[,c('p.val', 'q.val', 'set.size')], file='Results/DISEASED_UNSTIMULATION_v_DISEASED_STIMULATION_downregulated_msigdb_c5.xls')


