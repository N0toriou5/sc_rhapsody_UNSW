# complete sc analysis on 7 days
#### Code written by Daniele Mercatelli, University of Bologna
### Github: https://github.com/N0toriou5/sc_rhapsody_UNSW.git

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
library(patchwork)
library(scales)
library(Seurat)
library(scran)
library(stringr)
library(xlsx)
library(writexl)

##-- Data loading from CSV files derived from Seven Bridges into Seurat
# https://github.com/rbpatt2019/chooseR/blob/master/examples/1_seurat_pipeline.R for clustering automation
# follow HTO labelling as in # https://www.sciencedirect.com/science/article/pii/S2666166720302203
# follow HTO labelling as in # https://www.sciencedirect.com/science/article/pii/S2666166720302203
counts <- read.table("data/D7CM/Combined_D7_CONTROLMYE_RSEC_MolsPerCell.csv", skip = 7, sep = ",", header = TRUE, row.names = 1)
d7cm <- CreateSeuratObject(counts = t(counts), project="D7CM", min.cells = 3, min.features = 500)
#7Tm
counts <- read.table("data/D7TM/Combined_D7_TEPAMYE_RSEC_MolsPerCell.csv", skip = 7, sep = ",", header = TRUE, row.names = 1)
d7tm <- CreateSeuratObject(counts = t(counts), project = "D7TM",min.cells = 3, min.features = 500)
#d3c
counts <- read.table("data/D3C2/Combined_CONTROL_RSEC_MolsPerCell.csv", skip = 7, sep = ",", header = TRUE, row.names = 1)
d3c <- CreateSeuratObject(counts = t(counts), project="D3C2", min.cells=3, min.features = 500)
#d7t
# old tepa d7
counts <- read.table("data/D7T/Combined_D7_TEPA_RSEC_MolsPerCell.csv", skip = 7, sep = ",", header = TRUE, row.names = 1)
d7t <- CreateSeuratObject(counts = t(counts), project = "D7T", min.cells=3, min.features = 500)

#### Change folder for git 
setwd("git/sc_rhapsody_UNSW")
dir.create("plots")
dir.create("results")
# create list object to be merged in a large Seurat object
data.list <- list(D7CM=d7cm,D7TM=d7tm,D7T=d7t)
# merge the objects
combined <- merge(
  x = d3c,
  y = data.list, # can add a list like list(pbmc1k, pbmc5k, pbmc10k)
  add.cell.ids = c("D3C", "D7CM", "D7TM", "D7T")
)

seuset <- combined
save(seuset, file = "results/000_rawcounts.rda")

head(colnames(seuset))
table(seuset$orig.ident)
#  D3C2  D7CM   D7T  D7TM 
# 11662  7377  9412 11174

### QC filtering ----
# Add percent mito data in seurat object
genenames<-rownames(seuset@assays$RNA@counts)
mitogenes<-grep("^mt",genenames,value=TRUE)
seuset[["percent.mt"]] <- PercentageFeatureSet(seuset, pattern = "^mt.")
table(seuset@meta.data$orig.ident)
#  D3C2  D7CM   D7T  D7TM 
# 11662  7377  9412 11174

# Visualize QC metrics as a violin plot
png("plots/000_pre_filter_QC_MYE.png", w = 4000, h = 2000, res = 300)
VlnPlot(seuset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

#### Cell coming from different batches look very different, we would need to check for batch effects
plot1 <- FeatureScatter(seuset, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  geom_hline(yintercept=20,linetype="dashed") +
  geom_vline(xintercept = 50000,linetype="dashed") + 
  xlab("UMIs")+ylab("% of mitochondrial genes")
plot2 <- FeatureScatter(seuset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept=500,linetype="dashed") +
  geom_vline(xintercept = 50000,linetype="dashed") + 
  geom_hline(yintercept=7500,linetype="dashed")+xlab("UMIs") +
  ylab("Genes / Cell")

png("plots/000_featureScatter.png", w = 4000, h = 2000, res = 300)
plot1 + plot2
dev.off()

### Visualize QC metrics as histograms: Pre-filtering

png("plots/000_QC_histograms_MYE.png", w = 2000, h = 3000, res = 300)
par(mfrow = c(3, 1))
hist(seuset@meta.data$nCount_RNA, main = "Number of UMIs / Cell", xlab = "Counts")
abline(v = 50000, col = "red3", lty = 3)
hist(seuset@meta.data$nFeature_RNA, main = "Number of Genes / Cell", xlab = "Counts")
abline(v = 500, col = "red3", lty = 3)
abline(v = 7500, col = "red3", lty = 3)
hist(seuset@meta.data$percent.mt, main = "Mitochondrial Reads", xlab = "Percent")
abline(v = 20, col = "red3", lty = 3)
dev.off()

### Visualize QC metrics as histograms: Post-filtering
data <- subset(seuset, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 20 & nCount_RNA < 50000)
table(data@meta.data$orig.ident)
# D3C2  D7CM   D7T  D7TM 
# 11396  7282  9303 11040

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

png("plots/000_PCA_ident.png", h = 1500, w = 1500, res = 300)
DimPlot(object = data, reduction = "pca", pt.size = .1, group.by = "orig.ident")
dev.off()

# Batch correction
data <- RunHarmony(data, "orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(data, 'harmony')
harmony_embeddings[1:5, 1:5]
p1 <- DimPlot(object = data, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
png("plots/000_PCA_batchcorr.png", h = 1500, w = 1500, res = 300)
p1
dev.off()

# Clustering on number of significant PCs
# Determine percent of variation associated with each PC
pct <- data@reductions$harmony@stdev / sum(data@reductions$harmony@stdev) * 100
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
data <- FindNeighbors(object = data, dims = 1:17, reduction = 'harmony')
#data <- FindClusters(object = data, resolution = 0.5, algorithm = 3,group.singletons = FALSE,n.start = 25)
data <- FindClusters(object = data, resolution = 1.2, reduction= "harmony") # original plot is 0.7
data <- RunTSNE(data, dims = 1:17,tsne.method = "Rtsne", reduction = "harmony")
data <- RunUMAP(data, dims = 1:17, reduction = "harmony")
DimPlot(data, reduction = "tsne", label = T)
png("plots/000_tsne_origin.png", w = 2000, h = 2000, res = 300)
DimPlot(data, group.by = "orig.ident", reduction = "tsne", pt.size = 1.2) +
  ggtitle("Original Sample")
dev.off()
png("plots/000_tsne_all.png", w = 2000, h = 2000, res = 300)
DimPlot(object = data,pt.size = 1.2, reduction = 'tsne', group.by = 'seurat_clusters', label = TRUE) +
  ggtitle(paste(as.character(nrow(data@meta.data)), "cells (TEPA and Control)")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
png("plots/000_ident_tsne_aggregate.png", w = 2000, h = 2000, res = 300)
#DimPlot(object = data,pt.size=1.2,reduction='tsne',group.by='seurat_clusters',split.by = "orig.ident")
DimPlot(object = data, pt.size = .1, reduction = 'tsne', group.by = 'orig.ident')
dev.off()
png("plots/000_umap_all.png", w = 2000, h = 2000, res = 300)
DimPlot(object = data,pt.size = 1.2, reduction = 'umap',group.by = 'seurat_clusters', label = TRUE) +
  ggtitle(paste(as.character(nrow(data@meta.data)), "cells (TEPA and Control)")) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Find cells in each cluster

all.markers <- FindAllMarkers(object = data, only.pos = FALSE, min.pct = 0.25, min.diff.pct = 0.25)
save(all.markers, file = "results/000_all.markers.rda")
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
                           'clust10' = top05 %>% filter(cluster==10),
                           'clust11' = top05 %>% filter(cluster==11),
                           "clust12"= top05 %>% filter(cluster==12),
                           "clust13"= top05 %>% filter(cluster==13),
                           "clust14"= top05 %>% filter(cluster==14),
                           "clust15"= top05 %>% filter(cluster==15),
                           'clust16' = top05 %>% filter(cluster==16),
                           'clust17' = top05 %>% filter(cluster==17),
                           'clust18' = top05 %>% filter(cluster==18),
                           'clust19' = top05 %>% filter(cluster==19),
                           'clust20' = top05 %>% filter(cluster==20),
                           'clust21' = top05 %>% filter(cluster==21),
                           'clust22' = top05 %>% filter(cluster==22),
                           'clust23' = top05 %>% filter(cluster==23),
                           'clust24' = top05 %>% filter(cluster==24),
                           'clust25' = top05 %>% filter(cluster==25),
                           "clust26"= top05 %>% filter(cluster==26),
                           "clust27"= top05 %>% filter(cluster==27),
                           "clust28"= top05 %>% filter(cluster==28),
                           "clust29"= top05 %>% filter(cluster==29),
                           "clust30"= top05 %>% filter(cluster==30),
                           "clust31"= top05 %>% filter(cluster==31))

top05 <- top05[order(top05$avg_log2FC,decreasing = T),]
write_xlsx(marker.excel.pages, "results/000_topgenes_cluster.xlsx")

png("plots/000_Mycn_staining.png",w=2500,h=2500,res=300)
FeaturePlot(data,pt.size = 1, features = "Mycn",cols = c("grey","red"),
            max.cutoff = "q90", reduction = "tsne")
dev.off()
save(data, file="results/000_clustering.rda")

# Stage 1: unveil tumor heterogeneity
tumor <- subset(data, idents = c(2, 4, 6, 7, 8, 11, 15, 20, 25))
tumor <- subset(tumor, subset = orig.ident == c("D3C2","D7T"))

# Visualize QC metrics as a violin plot for tumor
png("plots/000_post_filter_QC_Tumor.png", w = 4000, h = 2000, res = 300)
VlnPlot(tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

png("plots/000_PCA_ident_tumor.png", h = 1500, w = 1500, res = 300)
DimPlot(object = tumor, reduction = "pca", pt.size = .1, group.by = "orig.ident")
dev.off()

png("plots/000_PCA_ident_tumor_origIdent.png", h = 1500, w = 1500, res = 300)
DimPlot(object = tumor, reduction = "pca", pt.size = .1)
dev.off()

png("plots/000_TSNE_ident_tumor.png", h = 1500, w = 3000, res = 300)
DimPlot(object = tumor, reduction = "tsne", pt.size = 1.2, split.by = "orig.ident")
dev.off()

tumor <- subset(tumor, subset = nFeature_RNA > 2000 & nFeature_RNA < 7500 & percent.mt < 20 & nCount_RNA <= 35000
                & nCount_RNA >= 10000)

png("plots/000_post_filter_QC_Tumor2.png", w = 4000, h = 2000, res = 300)
VlnPlot(tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# ##### Restart clustering
# # Seurat Method
# tumor <- FindVariableFeatures(tumor, selection.method = "vst", nfeatures = 2000) #Inf
# # Data Scaling
# tumor <- ScaleData(tumor)
# tumor <- RunPCA(tumor, pc.genes = data@var.genes, npcs = 20, verbose = FALSE)
# 
# p1 <- DimPlot(object = tumor, reduction = "pca", pt.size = .1, group.by = "orig.ident")
# p1
# 
# png("plots/000_PCA_ident_tumor_post.png", h = 1500, w = 1500, res = 300)
# DimPlot(object = tumor, reduction = "pca", pt.size = .1, group.by = "orig.ident")
# dev.off()
# 
# # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# # segregate this list into markers of G2/M phase and markers of S phase
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# s.genes <- orthologs(s.genes, species = "mouse")[,5]
# g2m.genes <- orthologs(g2m.genes, species = "mouse")[,5]
# 
# # Assign cell cycle score
# tumor <- CellCycleScoring(tumor, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# 
# # view cell cycle scores and phase assignments
# head(tumor[[]])
# tumor <- RunPCA(tumor, features = c(s.genes, g2m.genes))
# png("plots/000_PCA_ident_tumor_cellcycle.png", h = 1500, w = 3000, res = 300)
# DimPlot(tumor, reduction = "tsne", split.by = "orig.ident")
# dev.off()

##### The DE Tumor side
### Label treatment and ctrl
tumor$Treat <- ifelse(str_detect(colnames(tumor), "C")==TRUE, "CTRL", "TEPA")
table(Idents(tumor), tumor@meta.data$Treat)
rawcounts <- tumor@assays$RNA@counts
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
res <- res[res$baseMean>=0.5,]
res <- res[-grep("Rpl|Rps", rownames(res)),]
res <-res[order(res$log2FoldChange),]
dn <- rownames(res)[1:25]
up <- rownames(res)[(nrow(res)-24):nrow(res)]
labels <- c(up,dn)

png(paste0("plots/000_DESEQ_","whole_tumor",".png"), w=2500,h=2500, res=300)
#EnhancedVolcano(res,x="log2FoldChange",y="padj",lab=rownames(res))
gp<-EnhancedVolcano(res, subtitle = "",
                    lab = rownames(res),
                    selectLab = labels,
                    x = 'log2FoldChange',
                    y = 'padj',
                    xlim = c(-5, 5),
                    #ylim = c(0,100),
                    title = paste0("Whole Tumor",', TEPA vs. CTRL '),
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
res <- res[order(res$padj),]
write.xlsx(res, file = paste0("results/000_DESEQ_","whole_tumor",".xlsx"), row.names = T)
save(res, file = paste0("results/000_DESEQ_","whole_tumor",".rda"))

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



clusters<-levels(Idents(tumor))
#clusters<-clusters[-c(11,14)]
for (cluster in clusters){
  message(paste0("Doing ", cluster))
  set<-subset(x = tumor, idents = cluster)
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

##### use slalom and see what's next
library(slalom)
exprs_matrix <- as.matrix(tumor@assays$RNA@data)
sce <- SingleCellExperiment::SingleCellExperiment(assays =  list(logcounts =  exprs_matrix))
gmtfile <- "results/mh.all.v2022.1.Mm.symbols.gmt"
genesets <- GSEABase::getGmt(gmtfile)
# Generate a f-scLVM model
model <- newSlalomModel(sce, genesets)
# 50 annotated factors retained;  393 annotated factors dropped.
# 4016  genes retained for analysis.
# Initialize it
model <- initSlalom(model)
# Train it
model <- trainSlalom(model, nIterations =  10000) # train the model until it converges
save(model, file =  "results/000_slalom.rda")

#### Stage 2: investigate the immune clusters
immune <- subset(data, idents = c(2, 4, 6, 7, 8, 11, 15, 20, 25), invert = TRUE)

# Visualize QC metrics as a violin plot for tumor
png("plots/000_post_filter_QC_Immune.png", w = 4000, h = 2000, res = 300)
VlnPlot(immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

png("plots/000_PCA_ident_Immune.png", h = 1500, w = 1500, res = 300)
DimPlot(object = immune, reduction = "pca", pt.size = .1, group.by = "orig.ident")
dev.off()

png("plots/000_tsne_ident_Immune_labelled.png", h = 1500, w = 1500, res = 300)
DimPlot(object = immune, reduction = "tsne", pt.size = 1, label = TRUE)
dev.off()

png("plots/000_PCA_ident_immune_origIdent.png", h = 1500, w = 1500, res = 300)
DimPlot(object = immune, reduction = "pca", pt.size = .1)
dev.off()

png("plots/000_TSNE_ident_immune1.png", h = 1500, w = 4500, res = 300)
DimPlot(object = immune, reduction = "tsne", pt.size = 1.2, split.by = "orig.ident")
dev.off()

####### Show markers!
png("plots/000_cd4.png", h = 1500, w = 1500, res = 300)
FeaturePlot(immune, features = "Cd4", reduction = "tsne")
dev.off()

immune <- subset(immune, subset = nFeature_RNA > 1500 & nFeature_RNA < 4000 & percent.mt < 20 & nCount_RNA <= 20000
                & nCount_RNA >= 5000)

png("plots/000_post_filter_QC_Immune2.png", w = 4000, h = 2000, res = 300)
VlnPlot(immune, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

idx <- which(table(immune@meta.data$seurat_clusters)<=20)
groups <- names(idx)
immune <- subset(data, idents = c(0,1,2,4,5,6,7,8,11,15,16,20,22,25,29,30,31), invert = TRUE)

png("plots/000_tsne_real.png", h = 1500, w = 1500, res = 300)
DimPlot(object = immune, reduction = "tsne", pt.size = 1, label = TRUE)
dev.off()

markers <- c("Itgam",
             "Cxcr4",
             "Ly6g",
             "Itgax",
             "Il3ra",
             "Fcgr1",
             "Adgre1",
             "Csf1r",
             "Ly6c1",
             "Itga2b",
             "Fcer1a",
             "Ifitm1",
             "Cd3e",
             "Il7r",
             "Gzma",
             "Cd4",
             "Cd8a",
             "Cd79a",
             "Cd19",
             "Ms4a1",
             "Cd3d",
             "Cd8b1",
             "Trbc1",
             "Trbc2",
             "Trdc",
             "Ncam1",
             "Icam1")

png("plots/000_Dotplot_immune_clustering.png",h=2000,w=2500,res=300)
DotPlot(
  object = immune, features = markers
) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")+
  coord_flip() +xlab("")+ylab("")+ 
  scale_y_discrete(limits = levels(immune))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.off()
png("plots/000_Bcells.png", h = 1500, w = 2500, res = 300)
g <- c("Ms4a1", "Cd19")
FeaturePlot(immune, features = g, reduction = "tsne")
dev.off()

immune.markers <- FindAllMarkers(object = immune, only.pos = FALSE, min.pct = 0.25, min.diff.pct = 0.25)

### how many markers per cluster
all.markers %>% group_by(cluster) %>% summarize(n_p05 = sum(p_val_adj<0.05),
                                                npe10 = sum(p_val_adj<1e-10),
                                                npe100 = sum(p_val_adj<1e-100),
                                                np0 = sum(p_val_adj==0))
### write file with markers having adj p < 0.05
top05 <- all.markers %>% group_by(cluster) %>% filter(p_val_adj<0.05)

### create a single file
marker.excel.pages <- list('clust3' = top05 %>% filter(cluster==3),
                           'clust9' = top05 %>% filter(cluster==9),
                           'clust10' = top05 %>% filter(cluster==10),
                           'clust12' = top05 %>% filter(cluster==12),
                           'clust13' = top05 %>% filter(cluster==13),
                           'clust14' = top05 %>% filter(cluster==14),
                           'clust17' = top05 %>% filter(cluster==17),
                           'clust18' = top05 %>% filter(cluster==18),
                           'clust19' = top05 %>% filter(cluster==19),
                           'clust21' = top05 %>% filter(cluster==21),
                           'clust23' = top05 %>% filter(cluster==23),
                           'clust24' = top05 %>% filter(cluster==24),
                           "clust26"= top05 %>% filter(cluster==26),
                           "clust27"= top05 %>% filter(cluster==27),
                           "clust28"= top05 %>% filter(cluster==28))


write_xlsx(marker.excel.pages, "results/000_topgenes_cluster_immune.xlsx")

png("plots/000_Monocytes.png", h = 1500, w = 2500, res = 300)
g <- c("S100a9", "S100a8", "Cd68", "Cd14")
FeaturePlot(immune, features = g, reduction = "tsne")
dev.off()

png("plots/000_Tcells.png", h = 1500, w = 2500, res = 300)
g <- c("Cd4", "Cd8a", "Trdc", "Trbc2")
FeaturePlot(immune, features = g, reduction = "tsne")
dev.off()

library(Seurat)
library(SeuratObject)
library(SeuratDisk) #reference-based mapping, remotes::install_github("mojaveazure/seurat-disk")
library(stringr)
SaveH5Seurat(immune, filename = "results/Orazio.h5Seurat")
Convert("results/Orazio.h5Seurat", dest = "results/Orazio.h5ad")
source("D:/Archive/geneids.R")
convertMouseGeneList(rownames(immune))

library(dplyr)

mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")

convert_mouse_to_human <- function(gene_list){
  
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="mouse, laboratory"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="human"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  
  return (output)
}

save(mouse_human_genes, file = "results/convertTable.rda")
convert_mouse_to_human(rownames(immune)[1:5000])

# Plot markers of interest on a cluster basis
png("plots/000_ident_tsne_aggregate_clustering.png", w = 2000, h = 2000, res = 300)
DimPlot(object = data, pt.size = .1, reduction = 'tsne', group.by = 'seurat_clusters', label = TRUE)
dev.off()

markers <- c("Mycn",
             "Itgam",
             "Cxcr4",
             "Ly6g",
             "Itgax",
             "Il3ra",
             "Fcgr1",
             "Adgre1",
             "Csf1r",
             "Ly6c1",
             "Itga2b",
             "Fcer1a",
             "Ifitm1",
             "Cd3e",
             "Il7r",
             "Gzma",
             "Cd4",
             "Cd8a",
             "Cd79a",
             "Cd19",
             "Ms4a1",
             "Cd3d",
             "Cd8b1",
             "Trbc1",
             "Trbc2",
             "Trdc",
             #"Trgc1",
             #"Trgc2",
             "Ncam1",
             "Icam1")
png("plots/000_violin_clustering.png", w = 6000, h = 3000, res = 300)
VlnPlot(data, features = markers[1:6])
dev.off()
png("plots/000_violin_clustering2.png", w = 6000, h = 3000, res = 300)
VlnPlot(data, features = markers[7:12])
dev.off()
png("plots/000_violin_clustering3.png", w = 6000, h = 3000, res = 300)
VlnPlot(data, features = markers[13:18])
dev.off()
png("plots/000_violin_clustering4.png", w = 6000, h = 3000, res = 300)
VlnPlot(data, features = markers[19:24])
dev.off()
png("plots/000_violin_clustering5.png", w = 6000, h = 3000, res = 300)
VlnPlot(data, features = markers[25:28])
dev.off()


png("plots/000_Dotplot_clustering.png",h=2000,w=2500,res=300)
DotPlot(
  object = data, features = markers
) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")+
  coord_flip() +xlab("")+ylab("")+ 
  scale_y_discrete(limits = levels(data))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.off()

# Assign cell types according to Jourdin's schema
######### New custom clustering
levels(data$seurat_clusters) #32
data <- subset(data, idents = 31, invert = TRUE)
data <- subset(data, idents = 30, invert = TRUE)
named.clusters <- RenameIdents(data,
                               `0` = "CD4+ T cells",
                               `1` = "Monocytes",
                               `2` = "Tumor",
                               `3` = "Macrophages",
                               `4` = "Tumor",
                               `5` = "Monocytes",
                               `6` = "Tumor",
                               `7` = "Tumor",
                               `8` = "Tumor",
                               `9` = "Macrophages",
                               `10` = "Macrophages",
                               `11` = "Tumor",
                               `12` = "CD8+ T cells",
                               `13` = "DCs",
                               `14` = "B cells",
                               `15` = "Tumor",
                               `16` = "CD8+ T cells",
                               `17` = "Neutrophils",
                               `18` = "CD8+ T cells",
                               `19` = "Macrophages",
                               `20` = "Tumor",
                               `21` = "CD4+ T cells",
                               `22` = "Monocytes",
                               `23` = "DCs",
                               `24` = "DCs",
                               `25` = "Tumor",
                               `26` = "CD8+ T cells",
                               `27` = "B cells",
                               `28` = "Macrophages",
                               `29` = "Basophils")


png("plots/006_named_clustering_tsne_colored.png",w=2500,h=2000,res=300)
DimPlot(named.clusters, reduction = "tsne", label=T,pt.size=1,repel = T,
        cols = c('Tumor'= '#AFABAB', 'Neutrophils' = '#8DA9D7',
                 'DCs' = '#AC9DCB', 'Macrophages' = '#9ACA3C',
                 'Monocytes' = '#D85443', 'Basophils' = "#F388AE",
                 'CD4+ T cells' = '#EDC81A', 'CD8+ T cells' = '#FAA41A',
                 'B cells' = '#CDB18B'))
dev.off()

png("plots/006b_named_clustering_tsne_colored.png",w=2500,h=2000,res=300)
DimPlot(named.clusters, reduction = "tsne", label=F,pt.size=1,repel = T,
        cols = c('Tumor'= '#AFABAB', 'Neutrophils' = '#8DA9D7',
                 'DCs' = '#AC9DCB', 'Macrophages' = '#9ACA3C',
                 'Monocytes' = '#D85443', 'Basophils' = "#F388AE",
                 'CD4+ T cells' = '#EDC81A', 'CD8+ T cells' = '#FAA41A',
                 'B cells' = '#CDB18B'))
dev.off()


png("plots/006_named_clustering_umap_colored.png",w=2500,h=2000,res=300)
DimPlot(named.clusters, reduction = "umap", label=T,pt.size=.75, repel = T,
        cols = c('Tumor'= '#AFABAB', 'Neutrophils' = '#8DA9D7',
                 'DCs' = '#AC9DCB', 'Macrophages' = '#9ACA3C',
                 'Monocytes' = '#D85443', 'Basophils' = "#F388AE",
                 'CD4+ T cells' = '#EDC81A', 'CD8+ T cells' = '#FAA41A',
                 'B cells' = '#CDB18B'))
dev.off()

png("plots/006b_named_clustering_umap_colored.png",w=2500,h=2000,res=300)
DimPlot(named.clusters, reduction = "umap", label=F,pt.size=.75, repel = T,
        cols = c('Tumor'= '#AFABAB', 'Neutrophils' = '#8DA9D7',
                 'DCs' = '#AC9DCB', 'Macrophages' = '#9ACA3C',
                 'Monocytes' = '#D85443', 'Basophils' = "#F388AE",
                 'CD4+ T cells' = '#EDC81A', 'CD8+ T cells' = '#FAA41A',
                 'B cells' = '#CDB18B'))
dev.off()

### Label treatment and ctrl
named.clusters$Treat<-ifelse(str_detect(colnames(named.clusters),"C")==TRUE,"CTRL","TEPA")
table(Idents(named.clusters),named.clusters@meta.data$Treat)

data<-named.clusters
png("plots/000_CvsT.png", w = 2500, h = 2500, res = 300)
DimPlot(object = data, pt.size = 1, reduction = 'tsne',
        group.by = 'Treat', label = F, cols = c("grey", "red")) 
dev.off()

png("plots/000b_CvsT.png", w = 5000, h = 2500, res = 300)
DimPlot(object = data, pt.size = 1, reduction = 'tsne',
        group.by = 'Treat', split.by = 'Treat', label = F,
        cols = c("grey", "red")) 
dev.off()

### DotPlot
genes<-c("Mycn",
         "Itgam",
         "Cxcr4",
         "Ly6g",
         "Itgax",
         "Il3ra",
         "Fcgr1",
         "Adgre1",
         "Csf1r",
         "Ly6c1",
         "Itga2b",
         "Fcer1a",
         "Ifitm1",
         "Cd3e",
         "Il7r",
         "Gzma",
         "Cd4",
         "Cd8a",
         "Cd79a",
         "Cd19",
         "Ms4a1")
genes <- rev(genes)

png("plots/000_Dotplot.png",h=2000,w=2500,res=300)
DotPlot(
  object = data, features = genes
) + scale_colour_gradient2(low = "blue", mid = "white", high = "red")+
  coord_flip() +xlab("")+ylab("")+ 
  scale_y_discrete(limits = c('Tumor', 'Neutrophils','DCs', 'Macrophages','Monocytes',
                              'Basophils','CD4+ T cells', 'CD8+ T cells','B cells'))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
dev.off()

### Stacked barplot
subset <- subset(data, idents = 'Tumor', invert = TRUE)
my_table<-table(Idents(subset),subset@meta.data$Treat)

# relative freq table
totals<-colSums(my_table)
freq_tab_CTRL<-setNames((my_table[,1]/totals[1])*100,nm=rownames(my_table))
freq_tab_TEPA<-setNames((my_table[,2]/totals[2])*100,nm=rownames(my_table))
my_table2<-my_table
my_table2[,1]<-freq_tab_CTRL
my_table2[,2]<-freq_tab_TEPA
# Create vector of default ggplot2 colors

color_list <- c('#EDC81A','#D85443','#9ACA3C','#FAA41A','#AC9DCB','#CDB18B','#8DA9D7','#F388AE')


# Absolute frequency barplot
png("plots/000_stacked_barplot.png",w=3500,h=2500,res=300)
par(mar=c(5,4,4,20),xpd=T)
barplot(my_table2, main = "Relative frequency (%)",
        col = color_list,border = "white" )

legend("topright", 
       legend = rownames(my_table2), 
       fill = color_list,inset = c(- 0.4, 0),col=1:2)

dev.off()

### Differential expression analysis
# DESeq2
clusters<-levels(Idents(data))
#clusters<-clusters[-c(11,14)]
for (cluster in clusters){
  message(paste0("Doing ", cluster))
  set<-subset(x = data, idents = cluster)
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
  }
  else {
    # Wikipathways
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
