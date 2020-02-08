###################################################
#### Differential Expression Analysis by limma ####
###################################################

#source("http://www.bioconductor.org/biocLite.R")
#biocLite()
#biocLite('statmod')
#biocLite('edgeR')
#biocLite('limma')

library(statmod)
library(limma)
library(edgeR)


######### Comparison Among Different Sample Types ##########
setwd('~/Documents/Research/Cancer_Genomics/Data/Bladder/TCGA-BLCA-TEST/Gene_Expression_Quantification/Analysis/')
readCounts <- read.table('limma_input_RNAseq_counts_data_sample_type_TCGA_BLCA.txt', header = T, stringsAsFactors = F)
dim(readCounts)
readCounts[1:5,1:5]

readCounts <- readCounts[13:nrow(readCounts),]

readCounts <- data.frame(lapply(readCounts, as.numeric), row.names = rownames(readCounts))
readCounts[1:6,1:6]
dim(readCounts)


sampleType = as.matrix(methy['sample_type',])
sampleType
sampleType = as.factor(sampleType)

### summary of sample type
table(sampleType)


design <- model.matrix(~0+sampleType)
for (i in 1:length(colnames(design))) {
  colnames(design)[i] = sub('sampleType', '',colnames(design)[i])
}
head(design)


#if (length(levels(sampleType)) == 2) {
#  comp <- paste(levels(sampleType)[1], levels(sampleType)[2], sep = '-')
#  print (comp)


contrast.matrix <- makeContrasts(PrimaryTumor-SolidTissueNormal, levels=design) # do it manually
contrast.matrix

### TMM normalization on the raw counts by edgeR
dge = DGEList(counts = readCounts, group = sampleType)
dge$samples


### filtering and normalization
#keep <- rowSums(dge$counts) > 378
keep <- rowSums(cpm(dge) > 1) >= 0.5*length(sampleType) ## not sure how much is better
summary(keep)

dge <- dge[keep,,keep.lib.sizes = FALSE]
dge$samples
#write.table(rownames(dge), file='Genes_with_detectable_expression_value.txt', sep='\t', quote=FALSE)

dge = calcNormFactors(dge)
dge$samples

### voom transformation
v <- voom(dge, design, plot = TRUE)
v$E[1:5,1:5]

### differential expression analysis
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
results_type <- decideTests(fit2, adjust.method = 'BH', p.value = 0.01, lfc = 1, method='global')
rownames(results_type)
head(results_type)
summary(results_type)
vennDiagram(results_type, include=c('up', 'down'))

DEG_N2P <- topTable(fit2, coef=1, n = Inf, sort.by = 'P', p=0.01, lfc = 1)
DEG_N2P <- na.omit(DEG_N2P)

write.table(DEG_N2P, file = 'limma_001_2X_miRNAs_sample_type_TCGA_BLCA.txt', quote = FALSE, sep = '\t')

### volcano plot
geneList <- topTable(fit2, coef=1, n=Inf, adjust.method = 'BH', sort.by='P')
head(geneList$adj.P.Val)
geneList$threshold <- as.factor(abs(geneList$logFC)>1 & geneList$adj.P.Val<0.01)
head(geneList$threshold)
volcano <- ggplot(data=geneList, aes(x=geneList$logFC, y = -log10(geneList$adj.P.Val)))
#volcano <- ggplot(data=geneList, aes(x=-log10(geneList$adj.P.Val), y = geneList$logFC))
volcano+geom_point(aes(color=geneList$threshold), alpha=1, size=0.2) + xlab("log2 Fold Change") + ylab("-log10 adj.p-value") +
  scale_colour_manual(breaks = geneList$threshold, values = c('black', 'red')) + theme(legend.position="none") +
  geom_vline(xintercept = c(-1,1), color='darkgreen', linetype=3) + geom_hline(yintercept = 2, color='darkgreen',linetype=3)#+ xlim(0, 20)





##### Comparison Among Different Stages #####
setwd('~/Documents/Research/Cancer_Genomics/Data/Bladder/TCGA-BLCA-TEST/Gene_Expression_Quantification/Analysis/')
readCounts <- read.table('limma_input_RNAseq_counts_data_tumor_stage_TCGA_BLCA.txt', header = T)
dim(readCounts)
readCounts[1:5,1:5]

sampleType <- c()
for (columnName in colnames(readCounts)) {
  sampleType <- c(sampleType, strsplit(columnName, '.', fixed = TRUE)[[1]][1])
}
sampleType <- as.factor(sampleType)
table(sampleType)


design <- model.matrix(~0+sampleType)
head(design)

for (i in 1:length(colnames(design))) {
  colnames(design)[i] = sub('sampleType', '',colnames(design)[i])
}
head(design)

contrast.matrix <- makeContrasts(II2N=stageii-normal, III2N=stageiii-normal,
                                 IV2N=stageiv-normal, III2II=stageiii-stageii,
                                 IV2III=stageiv-stageiii,levels=design) # do it manually
contrast.matrix

### TMM normalization on the raw counts by edgeR
dge = DGEList(counts = readCounts, group = sampleType)
dge$samples


### filtering and normalization
#keep <- rowSums(dge$counts) > 378
keep <- rowSums(cpm(dge) > 1) >= 0.5*length(sampleType) ## not sure how much is better
summary(keep)

dge <- dge[keep,,keep.lib.sizes = FALSE]
dge$samples

dge = calcNormFactors(dge)
dge$samples

### voom transformation
v <- voom(dge, design, plot = TRUE)
v

### differential expression analysis
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results_stage <- decideTests(fit2, adjust.method = 'BH', p.value = 0.1, lfc = log(1.2,2), 
                             method='separate')

head(results_stage)
summary(results_stage)

write.table(results_stage, file='limma_001_2X_miRNAs_stage_summary_TCGA_BLCA.txt',quote = FALSE, sep = '\t')
### Stage vs Normal
vennDiagram(results_stage[,1:3], names=c('StageII-Normal','StageIII-Normal','StageIV-Normal'),
            include=c('up', 'down'),counts.col=c('red','blue'),
            circle.col=c('pink','lightblue','orange'),cex=1)
### Stage vs Stage
vennDiagram(results_stage[,4:5], names=c('StageIII-StageII','StageIV-StageIII'),
            include=c('up', 'down'),counts.col=c('red','blue'),
            circle.col=c('pink','lightblue'),cex=1)

?vennDiagram

pairs <- c('II2N', 'III2N', 'IV2N', 'III2II', 'IV2III')
pairs
for (pair in pairs) {
  print (pair)
  fileName <- paste('limma_001_2X_miRNAs', pair, 'TCGA_BLCA.txt', sep='_')
  #pairName = strsplit(pair,'2',fix=TRUE)
  DEG_stage <- topTable(fit2, coef=pair, n = Inf, sort.by = 'P')#, p=0.01, lfc = 1)
  DEG_stage <- na.omit(DEG_stage)
  
  write.table(DEG_stage, file = fileName, quote = FALSE, sep = '\t')
}


### combination stage and normal venn diagram
head(results_type)
colnames(results_type)
resultCombinatinon = cbind(results_type,results_stage[,4:5])
colnames(resultCombinatinon) = c(colnames(results_type), colnames(results_stage[,4:5]))
vennDiagram(resultCombinatinon, names=c('Tumor-Normal','StageIII-StageII','StageIV-StageIII'),
            include=c('up', 'down'),counts.col=c('red','blue'),
            circle.col=c('pink','lightblue','green3','orange'), cex=1)



### combine stage III and IV
newDesign <- cbind(design[,1], design[,2], design[,3]+design[,4])
colnames(newDesign) <- c('normal','stageii','stageiiiandiv')
head(newDesign)

newContrastMatrix <- makeContrasts(II2N=stageii-normal, IIIandIV2N=stageiiiandiv-normal,
                                 IIIandIV2II=stageiiiandiv-stageii,levels=newDesign)


### voom transformation
v <- voom(dge, newDesign, plot = TRUE)
v

### differential expression analysis
fit <- lmFit(v, newDesign)
fit2 <- contrasts.fit(fit, newContrastMatrix)
fit2 <- eBayes(fit2)

results_stage <- decideTests(fit2, adjust.method = 'BH', p.value = 0.01, lfc = 1, 
                             method='separate')
head(results_stage)
summary(results_stage)

#write.table(results_stage, file='limma_001_2X_miRNAs_stage_summary_TCGA_BLCA.txt',quote = FALSE, sep = '\t')
vennDiagram(results_stage, names=c('StageII-Normal','StageIIIandIV-Normal','StageIIIandIV-StageII'),
            include=c('up', 'down'),counts.col=c('red','blue'),
            circle.col=c('pink','lightblue','orange'),cex=1)






###### miRNA and mRNA correlation
cor( ,method=pearson) # to get correlation coefficient
cor.test(,method=pearson) # significance test

x = matrix(c(1,2,3,4,5,3,1,2,3,4), nrow=5, ncol=2)
y = matrix(c(6,5,3,1,2,3,4,5,3,2), nrow=5, ncol=2)
install.packages('Hmisc')
library(Hmisc)
rcorr(x,y)


cor(x=),y=as.matrix(c(6,5,3,1,2),c(3,4,5,3,2))
y=
y

n=20
df <- data.frame(tyrosine=runif(n), urea=runif(n), glucose=runif(n), inosine=runif(n))
df
df[,-1]
COR <- cor(as.matrix(df[,1]), as.matrix(df[,-1]))
COR

### gene ontology analysis
go.fisher <- goana(fit, universe = NULL, species = 'Hs')
topGo(go.fisher, sort = 'up', n =30)
topGo(go.fisher, sort = 'down', n =30)

class(fit)


### heatmap
library(gplots)
overlaps <- read.table('overlapping_gene_list_for_heatmap.txt', header = F)
overlaps <- c(as.matrix(overlaps[,1]))

logcpm <- cpm(dge, log = TRUE) # prior.count=2???
#logcpm <- logcpm[rownames(geneList) %in% rownames(dge),]
#logcpm <- logcpm[match(rownames(geneList), rownames(logcpm)),]
#logcpm <- logcpm[overlaps %in% rownames(dge),]
logcpm <- logcpm[match(overlaps, rownames(logcpm)),]



sampleNames <- c()
for (i in 1:length(colnames(dge))){
  sampleNames <- c(sampleNames, paste("A",i,sep = '-'))
}
colnames(logcpm) <- sampleNames

heatmap.2(as.matrix(logcpm[1:1000,1:30]), col=redgreen(75), trace='none', cexCol=0.9, labRow=NA)

### manually check cpm differences in each sample
o <- order(DEG_voom$adj.P.Val)
head(DEG_voom)
cpm(dge)[o[1:10],]
write.table(cpm(dge)[o[1:10],], file = 'check_cpm_in_each_sample.txt', quote = FALSE, sep = '\t')


### to give more weight to fold-changes in the ranking
fit <- treat(fit, lfc=log2(1.2))
xxx <- topTreat(fit, coef=ncol(design), n = Inf)

write.table(xxx, file = 'test.txt', quote = FALSE, sep = '\t')


1/(2^-0.52)

