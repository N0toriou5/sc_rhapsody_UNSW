###### Code written by Daniele Mercatelli, University of Bologna 2022
##### Create working environment -----
setwd("D:/Projects/Orazio/")
library(apeglm)
library(babelgene)
library(cowplot)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
library(enrichR)
library(ggforce)
library(glmGamPoi)
library(grid)
library(gridExtra)
library(harmony)
library(magrittr)
library(MAST)
library(patchwork)
library(scales)
library(Seurat)
library(scran)
library(stringr)
library(xlsx)
library(writexl)

#d3c
# follow HTO labelling as in # https://www.sciencedirect.com/science/article/pii/S2666166720302203
counts <- read.table("data/D3C2/Combined_CONTROL_RSEC_MolsPerCell.csv", skip = 7,
                     sep = ",", header = TRUE, row.names = 1)
d3c <- CreateSeuratObject(counts = t(counts), project = "D3C")
d3c_sampletags <- d3c@assays$RNA@counts@Dimnames[[2]]
# load SampleTags info (CD45+ cells and multiplets)
sampletags <- read.table("data/D3C2/_2_CONTROL_Sample_Tag_Calls.csv",
                         sep = ",", header = TRUE, row.names = 1)
d3c@meta.data$tags <- sampletags$Sample_Tag
table(d3c$tags)

# old tepa d7
counts <- read.table("data/D7T/Combined_D7_TEPA_RSEC_MolsPerCell.csv", skip = 7,
                     sep = ",", header = TRUE, row.names = 1)
d7t <- CreateSeuratObject(counts = t(counts),project = "D7T")
d7t_sampletags <- d7t@assays$RNA@counts@Dimnames[[2]]
# load SampleTags info (CD45+ cells and multiplets)
sampletags <- read.table("data/D7T/D7_TEPA_Sample_Tag_Calls.csv",
                         sep = ",", header = TRUE, row.names = 1)
d7t@meta.data$tags <- sampletags$Sample_Tag
table(d7t$tags)
# Undetermined 
# 11722
# Multiplet SampleTag10_mm SampleTag11_mm   Undetermined 
# 40           1744            799           6922 


### Sample Tags have not been used for D3C2: Can't use tags to discrimate 
### tumor vs. immune cells

### reload
#d3c
counts <- read.table("data/D3C2/Combined_CONTROL_RSEC_MolsPerCell.csv", skip = 7, sep = ",", header = TRUE, row.names = 1)
d3c <- CreateSeuratObject(counts = t(counts), project="D3C2", min.cells=3, min.features = 500)
#d7t
# old tepa d7
counts <- read.table("data/D7T/Combined_D7_TEPA_RSEC_MolsPerCell.csv", skip = 7, sep = ",", header = TRUE, row.names = 1)
d7t <- CreateSeuratObject(counts = t(counts), project = "D7T", min.cells=3, min.features = 500)
#### Change folder for git 
setwd("git/sc_rhapsody_UNSW")
seuset <- merge(d3c, y = d7t, add.cell.ids = c("D3C2", "D7T"))
save(seuset, file = "results/000_rawcounts_tumor.rda")
head(colnames(seuset))
table(seuset$orig.ident)
# D3C2   D7T 
# 11662  9412

### Library size
counts <- as.matrix(seuset@assays$RNA@counts)
counts_per_cell <- colSums(counts)
cat("counts per cell: ", counts_per_cell[1:5], "\n") ## counts for first 5 cells
summary(counts_per_cell)
counts_per_gene <- rowSums(counts)
cat("counts per gene: ", counts_per_gene[1:5], "\n")  ## counts for first 5 genes
summary(counts_per_gene)
genes_per_cell <- colSums(counts > 0) # count gene only if it has non-zero reads mapped.
cat("counts for non-zero genes: ", genes_per_cell[1:5])  ## counts for first 5 genes
hist(log10(counts_per_cell+1), main = 'counts per cell')
Mycn_exp <- counts["Mycn", ]
length(Mycn_exp) #21074
length(Mycn_exp[Mycn_exp>0]) #15548
tum <- names(Mycn_exp[Mycn_exp>0])
seuset@meta.data$rna_Mycn <- ifelse(colnames(seuset) %in% tum, 1, 0)
data <- subset(seuset, subset = rna_Mycn >= 1)

### QC filtering ----
# Add percent mito data in seurat object
genenames <- rownames(data@assays$RNA@counts)
mitogenes <- grep("^mt", genenames, value = TRUE)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt.")
# Visualize QC metrics as a violin plot
png("plots/000_pre_filter_QC_Tumor.png", w = 4000, h = 2000, res = 300)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


### Visualize QC metrics as histograms: Pre-filtering

png("plots/000_QC_histograms_tumor.png", w = 2000, h = 3000, res = 300)
par(mfrow = c(3, 1))
hist(data@meta.data$nCount_RNA, main = "Number of UMIs / Cell", xlab = "Counts")
abline(v = 50000, col = "red3", lty = 3)
hist(data@meta.data$nFeature_RNA, main = "Number of Genes / Cell", xlab = "Counts")
abline(v = 2500, col = "red3", lty = 3)
abline(v = 7500, col = "red3", lty = 3)
hist(data@meta.data$percent.mt, main = "Mitochondrial Reads", xlab = "Percent")
abline(v = 20, col = "red3", lty = 3)
dev.off()

### Visualize QC metrics as histograms: Post-filtering
data <- subset(data, subset = nFeature_RNA > 2500 & nFeature_RNA < 7500 & percent.mt < 20 & nCount_RNA < 50000)
table(data@meta.data$orig.ident)
# D3C2  D7T 
# 6965 7684 

png("plots/000_post_filter_QC.png",w=4000,h=2000,res=300)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# seurat standard normalization
# Normalize data with custom Seurat pipeline
data <- NormalizeData(data)
# Seurat Method
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000) #Inf
# Data Scaling
data <- ScaleData(data)
data <- RunPCA(data, pc.genes = data@var.genes, npcs = 20, verbose = FALSE)

p1 <- DimPlot(object = data, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p1

png("plots/000_PCA_tumor.png", h = 1500, w = 1500, res = 300)
DimPlot(object = data, reduction = "pca", pt.size = .1, group.by = "orig.ident")
dev.off()

#### No batch correction required!

# Clustering on number of significant PCs
# Determine percent of variation associated with each PC
pct <- data@reductions$pca@stdev / sum(data@reductions$pca@stdev) * 100
# Calculate cumulative percents for each PC
cum <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cum > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),  decreasing = T)[1] + 1 # last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2) # change to any other number
pcs #17

cell.num <- table(data@meta.data$orig.ident)
#Data dimensionality reduction and cluster-forming analysis
data <- FindNeighbors(object = data, dims = 1:17, reduction = 'pca')
#data <- FindClusters(object = data, resolution = 0.5, algorithm = 3,group.singletons = FALSE,n.start = 25)
data <- FindClusters(object = data, resolution = 0.5) # original plot is 0.7
data <- RunTSNE(data, dims = 1:17,tsne.method = "Rtsne", reduction = "pca")
data <- RunUMAP(data, dims = 1:17, reduction = "pca")
DimPlot(data, reduction = "tsne", label = T, pt.size = 1)
DimPlot(data, reduction = "umap", label = T, pt.size = 1, group.by = 'orig.ident')
png("plots/000_umap_origin_tumor.png", w = 2000, h = 2000, res = 300)
DimPlot(data, group.by = "orig.ident", reduction = "umap", pt.size = 1.2) +
  ggtitle("Original Sample")
dev.off()
png("plots/000_umap_all_tumor.png", w = 2000, h = 2000, res = 300)
DimPlot(object = data,pt.size = 1.2, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE) +
  ggtitle(paste(as.character(nrow(data@meta.data)), "cells (TEPA and Control)")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Find cells in each cluster

all.markers <- FindAllMarkers(object = data, only.pos = FALSE, min.pct = 0.25, min.diff.pct = 0.25)
save(all.markers, file = "results/000_all.markers_tumor.rda")
print(all.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC) %>% data.frame)

### how many markers per cluster
all.markers %>% group_by(cluster) %>% summarize(n_p05 = sum(p_val_adj<0.05),
                                                npe10 = sum(p_val_adj<1e-10),
                                                npe100 = sum(p_val_adj<1e-100),
                                                np0 = sum(p_val_adj==0))
### write file with markers having adj p < 0.05
top05 <- all.markers %>% group_by(cluster) %>% filter(p_val_adj<0.05)
### create a single file
marker.excel.pages <- list('clust0' = top05 %>% filter(cluster==0),
                           'clust1' = top05 %>% filter(cluster==1),
                           'clust2' = top05 %>% filter(cluster==2),
                           'clust3' = top05 %>% filter(cluster==3),
                           'clust4' = top05 %>% filter(cluster==4),
                           'clust5' = top05 %>% filter(cluster==5),
                           'clust6' = top05 %>% filter(cluster==6),
                           'clust7' = top05 %>% filter(cluster==7),
                           'clust8' = top05 %>% filter(cluster==8),
                           'clust9' = top05 %>% filter(cluster==9),
                           'clust10' = top05 %>% filter(cluster==10))
write_xlsx(marker.excel.pages, "results/000_topgenes_cluster_tumor.xlsx")

# exclude immune cells
data <- subset(data, ident = c(8,9,10), invert = TRUE)
DimPlot(data)
png("plots/000_umap_origin_tumor_no_immune.png", w = 2000, h = 2000, res = 300)
DimPlot(data, group.by = "orig.ident", reduction = "umap", pt.size = 1.2) +
  ggtitle("Original Sample")
dev.off()

png("plots/000_umap_clusters_tumor_no_immune.png", w = 2000, h = 2000, res = 300)
DimPlot(data, group.by = "seurat_clusters", reduction = "umap", pt.size = 1.2) +
  ggtitle("Original Sample")
dev.off()

png("plots/000_PCA_tumor_only.png", h = 1500, w = 1500, res = 300)
DimPlot(data, reduction = "pca", group.by = "orig.ident")
dev.off()
png("plots/000_features.png", w = 3000, h = 1500, res = 300)
FeaturePlot(data, features = c("Trbc2","Cd68"))
dev.off()

#### Taking out all immune cells left, the sign of a batch effect comes out!

save(data, file="results/000_clustering_tumor.rda")

##### The DE Tumor side
### Label treatment and ctrl
data$Treat <- ifelse(str_detect(colnames(data), "C")==TRUE, "CTRL", "TEPA")
table(Idents(data), data@meta.data$Treat)
rawcounts <- data@assays$RNA@counts
annotation <- matrix(nrow = ncol(rawcounts), ncol = 1, dimnames = list(colnames(rawcounts), "Treat"))
annotation[, "Treat"] <- ifelse(str_detect(colnames(rawcounts), "C")==TRUE, "Ctrl", "TEPA")
# DESeq2 block (filter out poorly expressed genes)
dds <- DESeqDataSetFromMatrix(countData = rawcounts, colData = annotation, design = ~Treat)
# low count filter - at least 10 with count of 5 or more
# https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
# example: keep <- rowSums(counts(sim) >= 5) >= 10
dds <- dds[rowSums(counts(dds)>=5)>=10,]
dds$Treat <- relevel(dds$Treat, ref = "Ctrl")
dds <- estimateSizeFactors(dds, type = "poscounts")

scr <- computeSumFactors(dds)
# use scran's sum factors:
sizeFactors(dds) <- sizeFactors(scr)
# Since our 'full' model only has one factor (Treat), the 'reduced' model (removing that factor) leaves us with nothing in our design formula.
# DESeq2 cannot fit a model with nothing in the design formula, and so in the scenario where you have no additional covariates the intercept
# is modeled using the syntax ~ 1.
dea <- DESeq(dds, parallel = FALSE, test = "LRT", useT = TRUE, minmu = 1e-6,
             minReplicatesForReplace = Inf, reduced = ~1, fitType = "glmGamPoi")

resultsNames(dea) #TreatTEPA"
res <- results(dea)
res <- as.data.frame(results(dea,name="TreatTEPA"))
# Volcano
res <- res[, -3]
res <- na.omit(res)
res <- res[res$baseMean>=0.1,] #min = 0.0115, 1st 0.6425
res <- res[-grep("Rpl|Rps", rownames(res)),]
res <-res[order(res$log2FoldChange),]
dn <- rownames(res)[1:25]
up <- rownames(res)[(nrow(res)-24):nrow(res)]
labels <- c(up,dn)

png("plots/000_NB_Tumor.png", w=2500,h=2500, res=300)
#EnhancedVolcano(res,x="log2FoldChange",y="padj",lab=rownames(res))
gp<-EnhancedVolcano(res, subtitle = "",
                    lab = rownames(res),
                    selectLab = labels,
                    x = 'log2FoldChange',
                    y = 'padj',
                    xlim = c(-5, 5),
                    #ylim = c(0,100),
                    title = "Whole Tumor, TEPA vs. CTRL ",
                    pCutoff = 0.05, #0.05 cutoff
                    FCcutoff = 0.25, # 2-fold change
                    labFace = "bold",
                    labSize = 3,
                    col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                    colAlpha = 4/5,
                    legendLabSize = 14,
                    legendIconSize = 4.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.3, colConnectors = 'gray51', maxoverlapsConnectors = Inf,
                    caption = paste0('Upregulated = ', nrow(res[res$log2FoldChange>0.25&res$padj<=0.05,]), ' genes',"\n",'Downregulated = ',
                                     nrow(res[res$log2FoldChange< -0.25&res$padj<=0.05,]), ' genes'))+
  theme(plot.title = element_text(hjust = 0.5))
print(gp)
dev.off()
res <- res[order(res$padj),]
write.xlsx(res, file = "results/000_DESEQ_whole_tumor.xlsx", row.names = T)
save(res, file = "results/000_DESEQ_whole_tumor.rda")

# GO

dbs <- listEnrichrDbs()
dbs <- c("WikiPathways_2019_Mouse", "KEGG_2019_Mouse", "HDSigDB_Mouse_2021")

de <- res
up <- rownames(de[de$log2FoldChange>0.25&de$padj<=0.05,])
dn <- rownames(de[de$log2FoldChange< -0.25&de$padj<=0.05,])

eup <- enrichr(up, dbs)
edn <- enrichr(dn, dbs)

# Wikipathways
up <- eup$WikiPathways_2019_Mouse
down <- edn$WikiPathways_2019_Mouse

up$type<-"up"
down$type<-"down"

up<-up[c(1:10),]
up<-up[order(up$Combined.Score),]
down<-down[c(1:10),]
down$Combined.Score<- (-1)*down$Combined.Score
down<-down[order(down$Combined.Score),]
gos<-rbind(down,up)
gos$Term<-gsub("WP(.*)","",gos$Term)
gos$Term<-paste0(gos$Term,ifelse(gos$P.value<=0.05,"*",""))
if (sum(duplicated(gos$Term))==0){
  #Diverging Barchart
  gos$Term<-factor(gos$Term,levels = gos$Term)
  png("plots/000_GO_NB_sc.png",w=2500,h=1500,res=300)
  gp<-ggplot(gos,aes(x=Term,y=Combined.Score,label=Combined.Score))+
    geom_bar(stat='identity',aes(fill=type),width=.5,position='dodge')+
    scale_fill_manual(name="Expression",
                      labels=c("Down regulated","Up regulated"),
                      values=c("down"="lightblue","up"="#f8766d"))+
    labs(subtitle="",
         title=paste0("Enriched in 7 days scTH-MYCN (","whole_tumor",")"))+
    coord_flip()+
    theme_bw()+ylab("EnrichR Combined Score")+
    theme(plot.title = element_text(hjust = 0.5))
  print(gp)
  dev.off()
} else {
  # Wikipathways
  up<-eup$WikiPathways_2019_Mouse
  down<-edn$WikiPathways_2019_Mouse
  
  up$type<-"up"
  down$type<-"down"
  
  up<-up[c(1:5),]
  up<-up[order(up$Combined.Score),]
  down<-down[c(1:5),]
  down$Combined.Score<- (-1)*down$Combined.Score
  down<-down[order(down$Combined.Score),]
  gos<-rbind(down,up)
  gos$Term<-gsub("WP(.*)","",gos$Term)
  gos$Term<-paste0(gos$Term,ifelse(gos$P.value<=0.05,"*",""))
  gos$Term<-factor(gos$Term,levels = gos$Term)
  #Diverging Barchart
  png(paste0("plots/000_GO_DESEQ_","whole_tumor","_sc.png"),w=2500,h=1500,res=300)
  gp<-ggplot(gos,aes(x=Term,y=Combined.Score,label=Combined.Score))+
    geom_bar(stat='identity',aes(fill=type),width=.5,position='dodge')+
    scale_fill_manual(name="Expression",
                      labels=c("Down regulated","Up regulated"),
                      values=c("down"="lightblue","up"="#f8766d"))+
    labs(subtitle="",
         title=paste0("Enriched in 7 days scTH-MYCN (","whole_tumor",")"))+
    coord_flip()+
    theme_bw()+ylab("EnrichR Combined Score")+
    theme(plot.title = element_text(hjust = 0.5))
  print(gp)
  dev.off()
}

### What if we remove the batch effect?
# split the dataset into a list of two seurat objects (TEPA and CTRL)
tumor.list <- SplitObject(data, split.by = "Treat")
# normalize and identify variable features for each dataset independently
tumor.list <- lapply(X = tumor.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = tumor.list)

# Perform integration
gene.anchors <- FindIntegrationAnchors(object.list = tumor.list, anchor.features = features)

data.combined <- IntegrateData(anchorset = gene.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(data.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)
DimPlot(data.combined, group.by = 'orig.ident')
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindClusters(data.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(data.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(data.combined, reduction = "umap", label = TRUE, repel = TRUE)

png("plots/000_dumplot_integrated.png", w = 4000, h = 1500, res = 300)
print(p1 + p2)
dev.off()

clusters<-levels(Idents(data.combined))
#clusters<-clusters[-c(11,14)]
for (cluster in clusters){
  message(paste0("Doing ", cluster))
  set<-subset(x = data.combined, idents = cluster)
  rawcounts<-set@assays$RNA@counts
  annotation<-matrix(nrow=ncol(rawcounts),ncol=1, dimnames=list(colnames(rawcounts),"Treat"))
  annotation[,"Treat"]<-ifelse(str_detect(colnames(rawcounts),"C")==TRUE,"Ctrl","TEPA")
  
  # DESeq2 block (filter out poorly expressed genes)
  dds<-DESeqDataSetFromMatrix(countData=rawcounts,colData=annotation,design=~Treat)
  # low count filter - at least 10 with count of 5 or more
  # https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
  # example: keep <- rowSums(counts(sim) >= 5) >= 10
  dds<-dds[rowSums(counts(dds)>=5)>=10,]
  dds$Treat<-relevel(dds$Treat,ref="Ctrl")
  dds <- estimateSizeFactors(dds, type="poscounts")
  
  scr <- computeSumFactors(dds)
  # use scran's sum factors:
  sizeFactors(dds) <- sizeFactors(scr)
  # Since our 'full' model only has one factor (Treat), the 'reduced' model (removing that factor) leaves us with nothing in our design formula.
  # DESeq2 cannot fit a model with nothing in the design formula, and so in the scenario where you have no additional covariates the intercept
  # is modeled using the syntax ~ 1.
  dea<-DESeq(dds,parallel=FALSE,test="LRT",useT=TRUE, minmu=1e-6,minReplicatesForReplace=Inf,reduced=~1,fitType = "glmGamPoi")
  
  resultsNames(dea) #TreatTEPA"
  res <- results(dea)
  res<-as.data.frame(results(dea,name="TreatTEPA"))
  # res[1:5,]
  # res <- results(dea)
  # #res <- lfcShrink(dea,coef = "TreatTEPA",
  # #                res=res)
  # res<-as.data.frame(res)
  # Volcano
  res<-res[,-3]
  res<-na.omit(res)
  res <- res[res$baseMean>=0.5,]
  res<-res[-grep("Rpl|Rps",rownames(res)),]
  res<-res[order(res$log2FoldChange),]
  dn<-rownames(res)[1:25]
  up<-rownames(res)[(nrow(res)-24):nrow(res)]
  labels<-c(up,dn)
  # Volcano Plots
  
  png(paste0("plots/006_DESEQ_",cluster,".png"), w=2500,h=2500, res=300)
  #EnhancedVolcano(res,x="log2FoldChange",y="padj",lab=rownames(res))
  gp<-EnhancedVolcano(res, subtitle = "",
                      lab = rownames(res),
                      selectLab = labels,
                      x = 'log2FoldChange',
                      y = 'padj',
                      xlim = c(-5, 5),
                      #ylim = c(0,100),
                      title = paste0(cluster,', TEPA vs. CTRL '),
                      pCutoff = 0.05, #0.05 cutoff
                      FCcutoff = 0.25, # 2-fold change
                      labFace = "bold",
                      labSize = 3,
                      col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                      colAlpha = 4/5,
                      legendLabSize = 14,
                      legendIconSize = 4.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.3,colConnectors = 'gray51',maxoverlapsConnectors = Inf,
                      caption = paste0('Upregulated = ', nrow(res[res$log2FoldChange>0.25&res$padj<=0.05,]), ' genes',"\n",'Downregulated = ',
                                       nrow(res[res$log2FoldChange< -0.25&res$padj<=0.05,]), ' genes'))+
    theme(plot.title = element_text(hjust = 0.5))
  print(gp)
  dev.off()
  res<-res[order(res$padj),]
  write.xlsx(res,file=paste0("results/006_DESEQ_",cluster,".xlsx"),row.names=T)
  save(res,file=paste0("results/006_DESEQ_",cluster,".rda"))
  # GO
  
  dbs<-listEnrichrDbs()
  dbs<-c("WikiPathways_2019_Mouse","KEGG_2019_Mouse","HDSigDB_Mouse_2021")
  
  de<-res
  up<-rownames(de[de$log2FoldChange>0.25&de$padj<=0.05,])
  dn<-rownames(de[de$log2FoldChange< -0.25&de$padj<=0.05,])
  
  eup<-enrichr(up,dbs)
  edn<-enrichr(dn,dbs)
  
  # Wikipathways
  up<-eup$WikiPathways_2019_Mouse
  down<-edn$WikiPathways_2019_Mouse
  
  up$type<-"up"
  down$type<-"down"
  
  up<-up[c(1:10),]
  up<-up[order(up$Combined.Score),]
  down<-down[c(1:10),]
  down$Combined.Score<- (-1)*down$Combined.Score
  down<-down[order(down$Combined.Score),]
  gos<-rbind(down,up)
  gos$Term<-gsub("WP(.*)","",gos$Term)
  gos$Term<-paste0(gos$Term,ifelse(gos$P.value<=0.05,"*",""))
  if (sum(duplicated(gos$Term))==0){
    #Diverging Barchart
    gos$Term<-factor(gos$Term,levels = gos$Term)
    png(paste0("plots/006_GO_DESEQ_",cluster,"_sc.png"),w=2500,h=1500,res=300)
    gp<-ggplot(gos,aes(x=Term,y=Combined.Score,label=Combined.Score))+
      geom_bar(stat='identity',aes(fill=type),width=.5,position='dodge')+
      scale_fill_manual(name="Expression",
                        labels=c("Down regulated","Up regulated"),
                        values=c("down"="lightblue","up"="#f8766d"))+
      labs(subtitle="",
           title=paste0("Enriched in 7 days scTH-MYCN (",cluster,")"))+
      coord_flip()+
      theme_bw()+ylab("EnrichR Combined Score")+
      theme(plot.title = element_text(hjust = 0.5))
    print(gp)
    dev.off()
  } else {
    # Wikipathways
    n <- which(duplicated(gos$Term))
    n <- n[1]
    n <- as.integer((n - 2)/2)
    up<-eup$WikiPathways_2019_Mouse
    down<-edn$WikiPathways_2019_Mouse
    
    up$type<-"up"
    down$type<-"down"
    
    up<-up[c(1:2),]
    up<-up[order(up$Combined.Score),]
    down<-down[c(1:2),]
    down$Combined.Score<- (-1)*down$Combined.Score
    down<-down[order(down$Combined.Score),]
    gos<-rbind(down,up)
    gos$Term<-gsub("WP(.*)","",gos$Term)
    gos$Term<-paste0(gos$Term,ifelse(gos$P.value<=0.05,"*",""))
    gos$Term<-factor(gos$Term,levels = gos$Term)
    #Diverging Barchart
    png(paste0("plots/006_GO_DESEQ_",cluster,"_sc.png"),w=2500,h=1500,res=300)
    gp<-ggplot(gos,aes(x=Term,y=Combined.Score,label=Combined.Score))+
      geom_bar(stat='identity',aes(fill=type),width=.5,position='dodge')+
      scale_fill_manual(name="Expression",
                        labels=c("Down regulated","Up regulated"),
                        values=c("down"="lightblue","up"="#f8766d"))+
      labs(subtitle="",
           title=paste0("Enriched in 7 days scTH-MYCN (",cluster,")"))+
      coord_flip()+
      theme_bw()+ylab("EnrichR Combined Score")+
      theme(plot.title = element_text(hjust = 0.5))
    print(gp)
    dev.off()
  }
}


# Test for DE features using the MAST package
Idents(data) <- data@meta.data$Treat
res2 <- FindMarkers(data, ident.1 = "TEPA", ident.2 = "CTRL", test.use = "MAST")
dim(res)
dim(res2)

# GO
res2 <- res2[-grep("Rpl|Rps", rownames(res2)),]
res2<-res2[order(res2$avg_log2FC),]
dn<-rownames(res2)[1:25]
up<-rownames(res2)[(nrow(res2)-24):nrow(res2)]
labels<-c(up,dn)
# Volcano Plots

png("plots/000_MAST_volcano.png", w=2500,h=2500, res=300)
#EnhancedVolcano(res,x="log2FoldChange",y="padj",lab=rownames(res))
gp<-EnhancedVolcano(res2, subtitle = "",
                    lab = rownames(res2),
                    selectLab = labels,
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    xlim = c(-5, 5),
                    #ylim = c(0,100),
                    title = "MAST, TEPA vs. Ctrl",
                    pCutoff = 0.05, #0.05 cutoff
                    FCcutoff = 0.25, # 2-fold change
                    labFace = "bold",
                    labSize = 3,
                    col = c('lightgrey', 'pink', 'lightblue', 'salmon'),
                    colAlpha = 4/5,
                    legendLabSize = 14,
                    legendIconSize = 4.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.3,colConnectors = 'gray51',maxoverlapsConnectors = Inf,
                    caption = paste0('Upregulated = ', nrow(res2[res2$avg_log2FC>0.25&res2$p_val_adj<=0.05,]), ' genes',"\n",'Downregulated = ',
                                     nrow(res2[res2$avg_log2FC< -0.25&res2$p_val_adj<=0.05,]), ' genes'))+
  theme(plot.title = element_text(hjust = 0.5))
print(gp)
dev.off()

de <- res2
up <- rownames(de[de$avg_log2FC>0.25&de$p_val_adj<=0.05,])
dn <- rownames(de[de$avg_log2FC< -0.25&de$p_val_adj<=0.05,])

eup <- enrichr(up, dbs)
edn <- enrichr(dn, dbs)

# Wikipathways
up <- eup$WikiPathways_2019_Mouse
down <- edn$WikiPathways_2019_Mouse

up$type<-"up"
down$type<-"down"

up<-up[c(1:10),]
up<-up[order(up$Combined.Score),]
down<-down[c(1:10),]
down$Combined.Score<- (-1)*down$Combined.Score
down<-down[order(down$Combined.Score),]
gos<-rbind(down,up)
gos$Term<-gsub("WP(.*)","",gos$Term)
gos$Term<-paste0(gos$Term,ifelse(gos$P.value<=0.05,"*",""))
if (sum(duplicated(gos$Term))==0){
  #Diverging Barchart
  gos$Term<-factor(gos$Term,levels = gos$Term)
  png(paste0("plots/000_GO_DESEQ_","whole_tumor","_sc.png"),w=2500,h=1500,res=300)
  gp<-ggplot(gos,aes(x=Term,y=Combined.Score,label=Combined.Score))+
    geom_bar(stat='identity',aes(fill=type),width=.5,position='dodge')+
    scale_fill_manual(name="Expression",
                      labels=c("Down regulated","Up regulated"),
                      values=c("down"="lightblue","up"="#f8766d"))+
    labs(subtitle="",
         title=paste0("Enriched in 7 days scTH-MYCN (","whole_tumor",")"))+
    coord_flip()+
    theme_bw()+ylab("EnrichR Combined Score")+
    theme(plot.title = element_text(hjust = 0.5))
  print(gp)
  dev.off()
} else {
  # Wikipathways
  up<-eup$WikiPathways_2019_Mouse
  down<-edn$WikiPathways_2019_Mouse
  
  up$type<-"up"
  down$type<-"down"
  
  up<-up[c(1:1),]
  up<-up[order(up$Combined.Score),]
  down<-down[c(1:1),]
  down$Combined.Score<- (-1)*down$Combined.Score
  down<-down[order(down$Combined.Score),]
  gos<-rbind(down,up)
  gos$Term<-gsub("WP(.*)","",gos$Term)
  gos$Term<-paste0(gos$Term,ifelse(gos$P.value<=0.05,"*",""))
  gos$Term<-factor(gos$Term,levels = gos$Term)
  #Diverging Barchart
  png(paste0("plots/000_GO_DESEQ_","whole_tumor","_sc.png"),w=2500,h=1500,res=300)
  gp<-ggplot(gos,aes(x=Term,y=Combined.Score,label=Combined.Score))+
    geom_bar(stat='identity',aes(fill=type),width=.5,position='dodge')+
    scale_fill_manual(name="Expression",
                      labels=c("Down regulated","Up regulated"),
                      values=c("down"="lightblue","up"="#f8766d"))+
    labs(subtitle="",
         title=paste0("Enriched in 7 days scTH-MYCN (","whole_tumor",")"))+
    coord_flip()+
    theme_bw()+ylab("EnrichR Combined Score")+
    theme(plot.title = element_text(hjust = 0.5))
  print(gp)
  dev.off()
}

### why not using the slalom thing?
##### use slalom and see what's next
library(slalom)
exprs_matrix <- as.matrix(data@assays$RNA@data)
dim(exprs_matrix) #13903 cols
table(data@meta.data$Treat)
# CTRL TEPA 
# 6442 7461 
tepa <- exprs_matrix[, colnames(data)[data$Treat=="TEPA"]]
sce <- SingleCellExperiment::SingleCellExperiment(assays =  list(logcounts =  tepa))
gmtfile <- "results/mh.all.v2022.1.Mm.symbols.gmt"
genesets <- GSEABase::getGmt(gmtfile)
# Generate a f-scLVM model
model <- newSlalomModel(sce, genesets)
# 50 annotated factors retained;  393 annotated factors dropped.
# 3920 genes retained for analysis.
# Initialize it
model <- initSlalom(model)
# Train it
model <- trainSlalom(model, nIterations =  10000) # train the model until it converges
save(model, file =  "results/000_slalom.rda")