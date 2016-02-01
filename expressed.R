# Name: expressed.R
# Auth: u.niazi@imperial.ac.uk
# Date: 1/02/2016
# Desc: list of expressed genes in each group


library(DESeq2)
library(org.Hs.eg.db)

## data loading and setting
dfDat = read.csv('Data_external/Counts/toc.csv', header=T)
rownames(dfDat) = dfDat$emtrez_id
dfDat = dfDat[,-1]

cn = colnames(dfDat)
fGroups = gsub('(\\w)\\d+', '\\1', cn)
fGroups = factor(fGroups, levels=c('C', 'H', 'D'))

mDat = as.matrix(dfDat)
colnames(mDat) = fGroups
mDat = mDat# + 1
sf = estimateSizeFactorsForMatrix(mDat)
mCounts = sweep(mDat, 2, sf, '/')
mCounts = na.omit(mCounts)

mMeans = matrix(NA, nrow = nrow(mCounts), ncol=3, dimnames = list(rownames(mCounts), levels(fGroups)))

for (i in seq_along(levels(fGroups))){
  mMeans[,i] = rowMeans(mCounts[,colnames(mCounts) %in% levels(fGroups)[i]])
}

# get a vector for each group
C = mMeans[,1]
head(C)
C = C[C >= 1]
C = log(C)
bins = hist(C, probability = T)
summary(C)
dn = dnorm(bins$mids, median(C), mad(C))
lines(bins$mids, dn, col=2, type='b')
# cut data into quantiles
cut.pts = cut(C, breaks = quantile(C, 0:10/10), labels = 1:10)
boxplot(C ~ cut.pts, main='Quantiles vs Expression - Group C', xlab='quantile', ylab='log expression')
df = select(org.Hs.eg.db, names(C), columns = c('GENENAME', 'SYMBOL'), keytype = 'ENTREZID')
df$group = cut.pts
write.csv(df, 'Results/expression_quantile_C.csv')

C = mMeans[,2]
head(C)
C = C[C >= 1]
C = log(C)
bins = hist(C, probability = T, main='H')
summary(C)
dn = dnorm(bins$mids, median(C), mad(C))
lines(bins$mids, dn, col=2, type='b')
# cut data into quantiles
cut.pts = cut(C, breaks = quantile(C, 0:10/10), labels = 1:10)
boxplot(C ~ cut.pts, main='Quantiles vs Expression - Group H', xlab='quantile', ylab='log expression')
df = select(org.Hs.eg.db, names(C), columns = c('GENENAME', 'SYMBOL'), keytype = 'ENTREZID')
df$group = cut.pts
write.csv(df, 'Results/expression_quantile_H.csv')

colnames(mMeans)
C = mMeans[,3]
head(C)
C = C[C >= 1]
C = log(C)
bins = hist(C, probability = T, main='D')
summary(C)
dn = dnorm(bins$mids, median(C), mad(C))
lines(bins$mids, dn, col=2, type='b')
# cut data into quantiles
cut.pts = cut(C, breaks = quantile(C, 0:10/10), labels = 1:10)
boxplot(C ~ cut.pts, main='Quantiles vs Expression - Group D', xlab='quantile', ylab='log expression')
df = select(org.Hs.eg.db, names(C), columns = c('GENENAME', 'SYMBOL'), keytype = 'ENTREZID')
df$group = cut.pts
write.csv(df, 'Results/expression_quantile_D.csv')