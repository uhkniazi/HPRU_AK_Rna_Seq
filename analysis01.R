# Name: analysis01.R
# Auth: u.niazi@imperial.ac.uk
# Date: 12/01/2016
# Desc: DE analysis for count matrix data using RNA-Seq


library(DESeq2)


## data loading and setting
dfDat = read.csv('Data_external/Counts/toc.csv', header=T)
rownames(dfDat) = dfDat$emtrez_id
dfDat = dfDat[,-1]

cn = colnames(dfDat)
fGroups = gsub('(\\w)\\d+', '\\1', cn)
fGroups = factor(fGroups, levels=c('C', 'H', 'D'))

mDat = as.matrix(dfDat)
dfDesign = data.frame(condition=fGroups, row.names = colnames(mDat))

## DE analysis
# call deseq2 constructor
oDseq = DESeqDataSetFromMatrix(mDat, dfDesign, design = ~ condition)
oDseq = DESeq(oDseq)
plotDispEsts(oDseq)

# get the results for each comparison
# where all three comparisons are performed
oRes.D.vs.C = results(oDseq, contrast = c('condition', 'D', 'C'))
oRes.H.vs.C = results(oDseq, contrast = c('condition', 'H', 'C'))
oRes.D.vs.H = results(oDseq, contrast = c('condition', 'D', 'H'))

plotMA(oRes.H.vs.C, main='H vs C')
plotMA(oRes.D.vs.C, main='D vs C')
plotMA(oRes.D.vs.H, main='D vs H')

# get results with significant p-values
dfD.vs.C = as.data.frame(oRes.D.vs.C[which(oRes.D.vs.C$padj < 0.1),])
dfH.vs.C = as.data.frame(oRes.H.vs.C[which(oRes.H.vs.C$padj < 0.1),])
dfD.vs.H = as.data.frame(oRes.D.vs.H[which(oRes.D.vs.H$padj < 0.1),])

nrow(dfD.vs.C)
nrow(dfD.vs.H)
nrow(dfH.vs.C)


## choose the comparison for plotting
library(org.Hs.eg.db)
# add annotation to the data set after selecting comparison
res = as.data.frame(oRes.D.vs.C)

rn = rownames(res)
df = select(org.Hs.eg.db, as.character(rn), c('SYMBOL'), 'ENTREZID')
df = df[!duplicated(df$ENTREZID),]
rownames(df) = df$ENTREZID
dfPlot = res
dfPlot = cbind(dfPlot[rn,], df[rn,])
dfPlot = na.omit(dfPlot)

## write csv file
write.csv(dfPlot, file='Results/DEAnalysis_D.vs.C.csv')

dfGenes = data.frame(P.Value=dfPlot$pvalue, logFC=dfPlot$log2FoldChange, adj.P.Val = dfPlot$padj, SYMBOL=dfPlot$SYMBOL)

f_plotVolcano = function(dfGenes, main, p.adj.cut = 0.1, fc.lim = c(-3, 3)){
  p.val = -1 * log10(dfGenes$P.Value)
  fc = dfGenes$logFC
  # cutoff for p.value y.axis
  y.cut = -1 * log10(0.01)
  col = rep('lightgrey', times=length(p.val))
  c = which(dfGenes$adj.P.Val < p.adj.cut)
  col[c] = 'red'
  plot(fc, p.val, pch=20, xlab='Fold Change', ylab='-log10 P.Value', col=col, main=main, xlim=fc.lim)
  abline(v = 0, col='grey', lty=2)
  abline(h = y.cut, col='red', lty=2)
  # second cutoff for adjusted p-values
  y.cut = quantile(p.val[c], probs=0.95)
  abline(h = y.cut, col='red')
  # identify these genes
  g = which(p.val > y.cut)
  lab = dfGenes[g, 'SYMBOL']
  text(dfGenes$logFC[g], y = p.val[g], labels = lab, pos=2, cex=0.6)
}

f_plotVolcano(dfGenes, 'D vs C')

## repeat for the second comparison
res = as.data.frame(oRes.D.vs.H)

rn = rownames(res)
df = select(org.Hs.eg.db, as.character(rn), c('SYMBOL'), 'ENTREZID')
df = df[!duplicated(df$ENTREZID),]
rownames(df) = df$ENTREZID
dfPlot = res
dfPlot = cbind(dfPlot[rn,], df[rn,])
dfPlot = na.omit(dfPlot)

## write csv file
write.csv(dfPlot, file='Results/DEAnalysis_D.vs.H.csv')

dfGenes = data.frame(P.Value=dfPlot$pvalue, logFC=dfPlot$log2FoldChange, adj.P.Val = dfPlot$padj, SYMBOL=dfPlot$SYMBOL)
f_plotVolcano(dfGenes, 'D vs H')

## group the genes by expression profile i.e. DE or not DE
cvCommonGenes = unique(c(rownames(dfD.vs.H), rownames(dfD.vs.C)))
mCommonGenes = matrix(NA, nrow=length(cvCommonGenes), ncol=2)
mCommonGenes[,1] = cvCommonGenes %in% rownames(dfD.vs.H)
mCommonGenes[,2] = cvCommonGenes %in% rownames(dfD.vs.C)
rownames(mCommonGenes) = cvCommonGenes
colnames(mCommonGenes) = c('D.vs.H', 'D.vs.C')

#### analysis by grouping genes
# create groups in the data based on 2^2-1 combinations
mCommonGenes.grp = mCommonGenes
set.seed(123)
dm = dist(mCommonGenes.grp, method='binary')
hc = hclust(dm)

# cut the tree at the bottom to create groups
cp = cutree(hc, h = 0.2)
# sanity checks
table(cp)
length(cp)
length(unique(cp))

mCommonGenes.grp = cbind(mCommonGenes.grp, cp)

### print and observe this table and select the groups you are interested in
temp = mCommonGenes.grp
temp = (temp[!duplicated(cp),])
temp2 = cbind(temp, table(cp))
rownames(temp2) = NULL
print(temp2)

# write csv file with gene names in each group
rn = rownames(mCommonGenes.grp[mCommonGenes.grp[,'cp'] == '1',])
length(rn)
head(mCommonGenes.grp[rn,])

df = select(org.Hs.eg.db, as.character(rn), c('SYMBOL'), 'ENTREZID')
df = df[!duplicated(df$ENTREZID),]
rownames(df) = df$ENTREZID
write.csv(df, 'Results/DEAnalysis_Genes_Group1.csv')

# repeat for other 2 groups
rn = rownames(mCommonGenes.grp[mCommonGenes.grp[,'cp'] == '2',])
length(rn)
head(mCommonGenes.grp[rn,])

df = select(org.Hs.eg.db, as.character(rn), c('SYMBOL'), 'ENTREZID')
df = df[!duplicated(df$ENTREZID),]
rownames(df) = df$ENTREZID
write.csv(df, 'Results/DEAnalysis_Genes_Group2.csv')

rn = rownames(mCommonGenes.grp[mCommonGenes.grp[,'cp'] == '3',])
length(rn)
head(mCommonGenes.grp[rn,])

df = select(org.Hs.eg.db, as.character(rn), c('SYMBOL'), 'ENTREZID')
df = df[!duplicated(df$ENTREZID),]
rownames(df) = df$ENTREZID
write.csv(df, 'Results/DEAnalysis_Genes_Group3.csv')


fSamples = fGroups
## get the count matrix
mCounts = counts(oDseq, normalized=T)
mCounts = na.omit(log(mCounts))
f = is.finite(rowSums(mCounts))
mCounts = mCounts[f,]
fGroups = fSamples
# data quality check
plot(density(rowMeans(mCounts)))

### create graph and clusters
library(org.Hs.eg.db)
library(downloader)
source('../CGraphClust/CGraphClust.R')

## change the significant clusters function
# get the significant clusters matrix and the p.values
setGeneric('getSignificantClusters', def = function(obj, mCounts, fGroups, ...) standardGeneric('getSignificantClusters'))
setMethod('getSignificantClusters', signature='CGraphClust', definition = function(obj, mCounts, fGroups, ...){
  #   # stabalize the data before performing DE
  #   if (bStabalize){
  #     mCounts = t(apply(mCounts, 1, function(x) f_ivStabilizeData(x, fGroups)))
  #     colnames(mCounts) = fGroups
  #   }  
  # get the marginal of each cluster
  mCent = getClusterMarginal(obj, mCounts, bScaled = F)
  # check which cluster shows significant p-values
  #p.vals = na.omit(apply(mCent, 1, function(x) pairwise.t.test(x, fGroups, p.adjust.method = 'BH')$p.value))
  #fSig = apply(p.vals, 2, function(x) any(x < 0.01))
  p.val = apply(mCent, 1, function(x) anova(lm(x ~ fGroups))$Pr[1])
  #p.val = apply(mCent, 1, function(x) oneway.test(x ~ fGroups)$p.value)
  p.val = p.adjust(p.val, method = 'BH')
  fSig = p.val < 0.01
  mCent = mCent[fSig,]
  p.val = p.val[fSig]
  # reorder the matrix based on range of mean
  rSort = apply(mCent, 1, function(x){ m = tapply(x, fGroups, mean); r = range(m); diff(r)}) 
  mCent = mCent[order(rSort, decreasing = T),]
  p.val = p.val[order(rSort, decreasing = T)]
  lRet = list(clusters=mCent, p.val=p.val)  
  return(lRet)
})


# plotting parameters
p.old = par()

# try different combinations of graphs
rn = rownames(mCommonGenes.grp[mCommonGenes.grp[,'cp'] == '1',])
length(rn)
# or 
rn = rownames(mCommonGenes.grp)

# select significant genes and prepare data for graphing
mCounts = mCounts[rownames(mCounts) %in% rn,]
colnames(mCounts) = fGroups
mCounts = mCounts[,order(fGroups)]
fGroups = fGroups[order(fGroups)]
mCounts = t(mCounts)
mCounts.bk = mCounts

dfMap = AnnotationDbi::select(org.Hs.eg.db, colnames(mCounts), 'UNIPROT', 'ENTREZID')
dfMap = na.omit(dfMap)

### load the uniprot2reactome mapping obtained from
# http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt
# get reactome data
url = 'http://www.reactome.org/download/current/UniProt2Reactome_All_Levels.txt'
dir.create('Data_external', showWarnings = F)
csReactomeFile = 'Data_external/UniProt2Reactome_All_Levels.txt'
# download the reactome file if it doesnt exist
if (!file.exists(csReactomeFile)) download(url, csReactomeFile)
dfReactome = read.csv(csReactomeFile, header = F, stringsAsFactors=F, sep='\t')
x = gsub('\\w+-\\w+-(\\d+)', replacement = '\\1', x = dfReactome$V2, perl = T)
dfReactome$V2 = x

## map reactome ids to uniprot ids
dfReactome.sub = dfReactome[dfReactome$V1 %in% dfMap$UNIPROT,]
# get the matching positions for uniprot ids in the reactome table
i = match(dfReactome.sub$V1, dfMap$UNIPROT)
dfReactome.sub$ENTREZID = dfMap$ENTREZID[i]
dfGraph = dfReactome.sub[,c('ENTREZID', 'V2')]
dfGraph = na.omit(dfGraph)

n = unique(dfGraph$ENTREZID)
mCounts = mCounts[,n]
print(paste('Total number of genes with Reactome terms', length(n)))
levels(fGroups)

# create a correlation matrix to decide cor cutoff
mCor = cor(mCounts)

# check distribution 
hist(sample(mCor, 1000, replace = F), prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
axis(1, at = seq(-1, 1, by=0.1), las=2)

# stabalize the data and check correlation again
mCounts.bk = mCounts
# stabalize the data
mCounts.st = apply(mCounts, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(mCounts.st) = fGroups

# create a correlation matrix
mCor = cor(mCounts.st)
# check distribution 
hist(sample(mCor, 1000, replace = F), prob=T, main='Correlation of genes', xlab='', family='Arial', breaks=20, xaxt='n')
axis(1, at = seq(-1, 1, by=0.1), las=2)

# use the unstabalized version
# create the graph cluster object
# using absolute correlation vs actual values lead to different clusters
oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.65, bSuppressPlots = F)

## general graph structure
set.seed(1)
plot.final.graph(oGr)

## community structure
## overview of how the commuinties look like
# plot the main communities in 2 different ways
ig = getFinalGraph(oGr)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = T, iSize = 10)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, mark.groups=NULL, edge.color='lightgrey')
set.seed(1)
ig = getFinalGraph(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 10)
plot(getCommunity(oGr), ig, vertex.label=NA, layout=layout_with_fr, 
     vertex.frame.color=NA, edge.color='darkgrey')


## centrality diagnostics
## centrality parameters should not be correlated significantly and the location of the central
## genes can be visualized
# look at the graph centrality properties
set.seed(1)
ig = plot.centrality.graph(oGr)
# plot the genes or vertex sizes by fold change
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F, iSize = 10)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='darkgrey')
par(p.old)

## the diagnostic plots show the distribution of the centrality parameters
# these diagnostics plots should be looked at in combination with the centrality graphs
plot.centrality.diagnostics(oGr)

# get the centrality parameters
mCent = mPrintCentralitySummary(oGr)

## top vertices based on centrality scores
## get a table of top vertices 
dfTopGenes.cent = dfGetTopVertices(oGr, iQuantile = 0.90)
rownames(dfTopGenes.cent) = dfTopGenes.cent$VertexID
# assign metadata annotation to these genes and clusters
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
df = f_dfGetGeneAnnotation(as.character(dfTopGenes.cent$VertexID))
dfTopGenes.cent = cbind(dfTopGenes.cent[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
dfCluster = dfCluster[as.character(dfTopGenes.cent$VertexID),]
dfTopGenes.cent = cbind(dfTopGenes.cent, Cluster=dfCluster$cluster)

dir.create('Results', showWarnings = F)
write.csv(dfTopGenes.cent, file='Results/Top_Centrality_Genes.csv')

## if we want to look at the expression profiles of the top genes
# plot a heatmap of these top genes
library(NMF)
m1 = mCounts[,as.character(dfTopGenes.cent$VertexID)]
m1 = scale(m1)
m1 = t(m1)
# threshhold the values
m1[m1 < -3] = -3
m1[m1 > 3] = 3
rownames(m1) = as.character(dfTopGenes.cent$SYMBOL)
# draw the heatmap  color='-RdBu:50'
aheatmap(m1, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = TRUE, 
         annColors=NA, Colv=NA)


## in addition to heatmaps the graphs can be plotted
# plot a graph of these top genes
# plot for each contrast i.e. base line vs other level
lev = levels(fGroups)[-1]
m = mCounts
m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = induced_subgraph(getFinalGraph(oGr), vids = as.character(dfTopGenes.cent$VertexID))
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=30)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, vertex.label.cex=0.2, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey', 
       main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}


### Looking at the largest clique can be informative in the graph
# plot the graph with location of the clique highlighted
set.seed(1)
ig = plot.graph.clique(oGr)
ig = f_igCalculateVertexSizesAndColors(ig, t(mCounts), fGroups, bColor = F)
par(mar=c(1,1,1,1)+0.1)
set.seed(1)
plot(ig, vertex.label=NA, layout=layout_with_fr, vertex.frame.color=NA, edge.color='lightgrey')

# plot the largest clique at each grouping contrast
lev = levels(fGroups)[-1]
m = mCounts
#m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
#rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = induced_subgraph(getFinalGraph(oGr), vids = unlist(getLargestCliques(oGr)))
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=80)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, layout=layout_with_fr, main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}


## instead of looking at individual genes we can look at clusters
## we can look at the problem from the other direction and look at clusters instead of genes
# some sample plots
# mean expression of groups in every cluster
par(p.old)
plot.mean.expressions(oGr, t(mCounts), fGroups, legend.pos = 'bottomleft', main='Total Change in Each Cluster', cex.axis=0.7)
# only significant clusters
par(mar=c(7, 3, 2, 2)+0.1)
plot.significant.expressions(oGr, t(mCounts), fGroups, main='Significant Clusters', lwd=1, bStabalize = T, cex.axis=0.7)
# principal component plots
pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = T)
par(mar=c(4,2,4,2))
biplot(pr.out, cex=0.8, cex.axis=0.8, arrow.len = 0)
# plot summary heatmaps
# marginal expression level in each cluster
plot.heatmap.significant.clusters(oGr, t(mCounts), fGroups, bStabalize = F)
# plot variance of cluster
m = getSignificantClusters(oGr, t(mCounts), fGroups)$clusters
#m = getClusterMarginal(oGr, t(mCounts))
# plot.cluster.variance(oGr, m[c('1280218', '1280215'),], fGroups, log = F)

csClust = rownames(m)
length(csClust)

pdf('Temp/clusters.var.pdf')
i = 1
temp = t(as.matrix(m[csClust[i],]))
rownames(temp) = csClust[i]
plot.cluster.variance(oGr, temp, fGroups, log=FALSE); i = i+1
par(mfrow=c(2,2))
boxplot.cluster.variance(oGr, m, fGroups, log=T, iDrawCount = length(csClust))
dev.off(dev.cur())
# cluster names
i = which(dfReactome.sub$V2 %in% csClust)
dfCluster.name = dfReactome.sub[i,c('V2', 'V4')]
dfCluster.name = dfCluster.name[!duplicated(dfCluster.name$V2),]
rownames(dfCluster.name) = NULL
dfCluster.name


#### plot a graph of clusters 
#m = getSignificantClusters(oGr, t(mCounts), fGroups, bStabalize = T)
dfCluster = getClusterMapping(oGr)
colnames(dfCluster) = c('gene', 'cluster')
rownames(dfCluster) = dfCluster$gene
# how many genes in each cluster
sort(table(dfCluster$cluster))
#csClust = rownames(m$clusters)
csClust = as.character(unique(dfCluster$cluster))

# graph
lev = levels(fGroups)[-1]
m = mCounts
#m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
#rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1)
for(i in 1:length(lev)){
  ig = getClusterSubgraph(oGr, csClust)
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=10)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, vertex.label.cex=0.14, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey',
       main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}

df = f_dfGetGeneAnnotation(as.character(dfCluster$gene))
dfCluster = cbind(dfCluster[as.character(df$ENTREZID),], SYMBOL=df$SYMBOL, GENENAME=df$GENENAME)
write.csv(dfCluster, file='Results/Clusters.csv')

##### Various plots for one cluster of choice
csClust = '211859'

lev = levels(fGroups)[-1]
m = mCounts
#m = apply(m, 2, function(x) f_ivStabilizeData(x, fGroups))
#rownames(m) = rownames(mCounts)
par(mar=c(1,1,1,1)+0.1, mfrow=c(1,2))
for(i in 1:length(lev)){
  ig = getClusterSubgraph(oGr, csClust)
  fG = factor(fGroups, levels= c(levels(fGroups)[1], lev[-i], lev[i]) )
  ig = f_igCalculateVertexSizesAndColors(ig, t(m), fG, bColor = T, iSize=60)
  n = V(ig)$name
  lab = f_dfGetGeneAnnotation(n)
  V(ig)$label = as.character(lab$SYMBOL)
  set.seed(1)
  plot(ig, vertex.label.cex=0.7, layout=layout_with_fr, vertex.frame.color='darkgrey', edge.color='lightgrey',
       main=paste(lev[i], 'vs', levels(fGroups)[1]))
  legend('topright', legend = c('Underexpressed', 'Overexpressed'), fill = c('lightblue', 'pink'))
}

par(p.old)
# heatmap of the genes
ig.sub = getClusterSubgraph(oGr, csClustLabel = csClust)
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
mC = t(mCounts)
mC = mC[n$ENTREZID,]
rownames(mC) = n$SYMBOL
mC = t(scale(t(mC)))
# threshhold the values
mC[mC < -3] = -3
mC[mC > +3] = +3
# draw the heatmap
hc = hclust(dist(mC))
aheatmap(mC, color=c('blue', 'black', 'red'), breaks=0, scale='none', Rowv = hc, annRow=NA, 
         annColors=NA, Colv=NA)


# if we want to plot variance of one gene at a time
n = f_dfGetGeneAnnotation(V(ig.sub)$name)
mC = t(mCounts)
mC = mC[n$ENTREZID,]
rownames(mC) = n$SYMBOL
rn = rownames(mC)
length(rn)
i = 1

par(p.old)
par(mfrow=c(2,2))
boxplot.cluster.variance(oGr, (mC), fGroups, iDrawCount = length(rn))
temp = t(as.matrix(mC[rn[i],]))
rownames(temp) = rn[i]
plot.cluster.variance(oGr, temp, fGroups, log=FALSE); i = i+1

