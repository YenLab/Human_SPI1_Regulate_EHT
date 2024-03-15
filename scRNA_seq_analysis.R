library(clusterProfiler)
library(ggplot2)
library(ggsignif)
library(ggthemes)
library(httr)
library(jsonlite)
library(magrittr)
library(magrittr)
library(monocle)
library(org.Hs.eg.db)
library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
library(SingleCellExperiment)
library(tidyverse)


setwd("./Input_scRNAseq")

############################# ----------------------- #############################
############################# In vivo scRNA-seq UMAP  #############################
############################# ----------------------- #############################

## Load cell meta data provided by Calvanese et al.
## Data source: https://github.com/mikkolalab/Human-HSC-Ontogeny - Fig3_25
meta <- read.csv("seurat_metadata.csv",row.names = 1)

## Load single-cell matrices outputted by Cell Ranger
## Rename cell barcode to match the cell meta data
## AGM4.5: XXXXX-3; AGM5: XXXXX-1; AGM5.5: XXXXXX-2
invivo <- cbind(Read10X_h5("AGM4.5.h5") %>% set_colnames(str_replace(colnames(.),"-1","-3")),
           Read10X_h5("AGM5.h5"),
           Read10X_h5("AGM5.5.h5") %>% set_colnames(str_replace(colnames(.),"-1","-2"))) %>% 
  .[,rownames(meta)] %>% 
  CreateSeuratObject()

## Assign cell information to seurat object
invivo$cluster <- meta$seurat_clusters # identical(colnames(invivo),rownames(meta)): TRUE
invivo$celltype <- meta$celltype
Idents(invivo) <- factor(invivo$celltype, levels = c("vEC","AEC","HEC","HSC","non-HSC"))

## Run scRNA-seq analysis pipline following instruction from Fig3_25_seurat_analysis.R
invivo[["percent.mt"]] <- PercentageFeatureSet(invivo, pattern = "^MT-")
invivo <- invivo %>% 
  subset(., subset = nFeature_RNA > 100 & percent.mt < 10) %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>% 
  RunUMAP(dims = 1:23)

## UMAP plot
invivo[['umap']]@cell.embeddings <- -invivo[['umap']]@cell.embeddings

color <- c("#20b2aa","#33a02c","#c51b7d","#e31a1c","darkgray") %>%
  setNames(c("vEC","AEC","HEC","HSC","non-HSC"))

DimPlot(invivo, label = T, label.size = 5, cols = color, pt.size = 2)+NoLegend()

## Marker expression
marker <- c("NRP2", "NR2F2", "IL33", "GJA4", "GJA5", "CXCR4", "SOX17",
            "TMEM100", "EDN1", "LTBP4", "ALDH1A1","MECOM","MLLT3","HOXA9",
            "MYCN", "RUNX1", "CD44", "KCNK17", "MYB", "SPINK2", "SPN", "PTPRC")

marker_subset <- c("NR2F2", "CDH5","GJA5", "RUNX1","PTPRC", "SPINK2", "HLF")

DotPlot(invivo, features = marker, scale = F, cols = c("lightgrey","#bd0026"))+
  theme(axis.title = element_blank())

DotPlot(invivo, features = marker_subset, scale = F, cols = c("lightgrey","#bd0026"))+
  theme(axis.title = element_blank())

saveRDS(invivo,"invivo.rds")

############################# ----------------------- #############################
############################# In vitro scRNA-seq UMAP #############################
############################# ----------------------- #############################

## Load single-cell matrices deposited on GEO by Shen et al.
## Cell number: D0-11934; D2-37586; D4-9122; D6-5379

invitro <- cbind(Read10X_h5("GSM4338092_D0_H1_Nor_10X.h5"),
                 Read10X_h5("GSM4338093_D2_Nor_10X.h5"),
                 Read10X_h5("GSM4338094_D4_Nor_10X.h5"),
                 Read10X_h5("GSM4338095_D6_Nor_10X.h5")) %>% CreateSeuratObject()

## Assign culturing days to seurat object
invitro$day <- c(rep("D0",11934),rep("D2",37586),rep("D4",9122),rep("D6",5379))

## Run scRNA-seq analysis pipline following instruction from Seurat:
## 'https://satijalab.org/seurat/articles/sctransform_vignette'
invitro[["percent.mt"]] <- PercentageFeatureSet(invitro, pattern = "^MT-")
invitro <- subset(invitro, subset = nFeature_RNA > 2000 & percent.mt < 10)
invitro <- invitro %>% 
  SCTransform(do.correct.umi = F) %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>% 
  FindNeighbors() %>% 
  FindClusters(resolution = 0.2)

## Rename seurat clusters as cell types
invitro <- RenameIdents(invitro,
                        "0"="MPC","4"="MPC","1"="hPSC","8"="hPSC","3"="D4_Mes",
                        "5"="D6_Mes","2"="D4_EC","6"="D6_EC","7"="HPC")

## UMAP plot
DimPlot(invitro, label = T, label.size = 5)+NoLegend()


## Marker expression
marker <- c('NANOG', "SOX2", "DNMT3B", "ZFP42", "TBXT", "MIXL1", "MESP1",
            "TBX6", "CDX1", "PDGFRA", "PDGFRB", "COL1A1", "CD34", "CDH5", "PECAM1",
            "FLT1", "RUNX1", "GATA2", "SPN", "SPI1", "GFI1B", "ITGA2B")

marker_subset <- c("SOX2","TBXT", "PDGFRA", "CDH5", "CD34", "RUNX1", "SPN")

DotPlot(invitro, features = marker, scale = F, cols = c("lightgrey","#bd0026"))+
  theme(axis.title = element_blank())

DotPlot(invitro, features = marker_subset, scale = F, cols = c("lightgrey","#bd0026"))+
  theme(axis.title = element_blank())


saveRDS(invitro,"invitro.rds")
############################# ---------------- #############################
############################# In vivo monocle2 #############################
############################# ---------------- #############################

## Load pre-built Seurat object
invivo <- readRDS("invivo.rds") %>% 
  SetIdent(value = .$celltype) %>% 
  `DefaultAssay<-`(value = "RNA")

invivo <- invivo[,Idents(invivo) %in% c("AEC","HEC","HSC")]

cds <- as.CellDataSet(invivo)

## Size factors are used to correct for differences in library sizes between cells
## Library size refers to the total number of reads obtained from a cell
cds <- estimateSizeFactors(cds)

## Estimates the dispersion values using a maximum likelihood estimation (MLE) approach,
## which models the relationship between the mean and variance of gene expression.
cds <- estimateDispersions(cds)

# Filter on genes
cds <- detectGenes(cds, min_expr = 0.5)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

# Calculate DE genes by clusters
diff_test_res <- differentialGeneTest(cds[expressed_genes,], cores=2,
                                      fullModelFormulaStr = '~celltype')

## Select the ordering genes
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))

## Reduce dimension and order cells
cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)
cds <- reduceDimension(cds)
cds <- orderCells(cds)

# Plot trajectory
color <- c("#33a02c","#c51b7d","#e31a1c") %>%
  setNames(c("AEC","HEC","HSC"))

plot_cell_trajectory(cds, color_by = "celltype",
                     show_branch_points=FALSE, cell_size=3)+
  theme(legend.position = c(.5,.3))+
  scale_color_manual(values = color)

plot_cell_trajectory(cds, color_by = "Pseudotime",
                     show_branch_points=FALSE, cell_size=3)+
  theme(legend.position = c(.5,.3))+
  scale_colour_gradientn(colours = c("#440154FF","#443A83FF","#31688EFF", "#21908CFF",
                                     "#35B779FF", "#8FD744FF", "#FDE725FF"))

## Plot gene
plot_cell_trajectory(cds, markers = "RUNX1",use_color_gradient = T,
                     show_branch_points=FALSE, cell_size=3)+
  theme(legend.position = c(.5,.3))+
  scale_color_gradient2(low = "#4292c6",mid = "white", midpoint = 0,high = "#bd0026")

save(file="monocle2_invivo.Rdata", cds)

## Pseudo Time
df_invivo <- pData(cds)
df_invivo <- cbind(df_invivo,t(invivo[["RNA"]]@data))
save(file = "invivo_pseudotime.Rdata", df_invivo)

############################# ----------------- #############################
############################# In vitro monocle2 #############################
############################# ----------------- #############################

## Load pre-built Seurat object
invitro <- readRDS("invitro.rds") %>% `DefaultAssay<-`(value = "RNA")
invitro <- invitro[,Idents(invitro) %in% c("D4_EC","D6_EC","HPC")]

invitro$idents <- Idents(invitro)
cds <- as.CellDataSet(invitro)

## Size factors are used to correct for differences in library sizes between cells
## Library size refers to the total number of reads obtained from a cell
cds <- estimateSizeFactors(cds)

## Estimates the dispersion values using a maximum likelihood estimation (MLE) approach,
## which models the relationship between the mean and variance of gene expression.
cds <- estimateDispersions(cds)

## Filter on genes
cds <- detectGenes(cds, min_expr = 0.5)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

## Calculate DE genes by clusters
diff_test_res <- differentialGeneTest(cds[expressed_genes,], cores=2,
                                      fullModelFormulaStr = '~idents')

## Select the ordering genes
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))

## Reduce dimension and order cells
cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)
cds <- reduceDimension(cds)
cds <- orderCells(cds)

## Specify root cells (input for root_state in orderCells)
## check state using pData(cds) after you run orderCells(cds)
coord <- reducedDimS(cds) %>% t()
root_cell <- rownames(coord)[coord[,1]>=5 & coord[,2]>=9]
cds$State <- as.vector(cds$State)
cds$State[colnames(cds) %in% root_cell] <- "root"

cds <- orderCells(cds,root_state = "root")

## Plot trajectory
color <- c("#33a02c","#807dba","#e31a1c") %>% setNames(c("D4_EC","D6_EC","HPC"))

plot_cell_trajectory(cds, color_by = "idents",
                     show_branch_points=FALSE, cell_size=2)+
  theme(legend.position = c(.5,.3))+
  scale_color_manual(values = color)

plot_cell_trajectory(cds, color_by = "Pseudotime",
                     show_branch_points=FALSE, cell_size=3)+
  theme(legend.position = c(.5,.3))+
  scale_colour_gradientn(colours = c("#440154FF","#443A83FF","#31688EFF","#21908CFF",
                                     "#35B779FF", "#8FD744FF", "#FDE725FF"))

save(file="monocle2_invitro.Rdata", cds)

## Pseudo Time and gene expression Table

df_invitro <- pData(cds) %>% .[.$idents!="D6_EC" & .$State %in% c("root",1,3),]

invitro <- invitro[,rownames(df_invitro)]
invitro <- invitro %>% NormalizeData()

df_invitro <- cbind(df_invitro,t(invitro[["RNA"]]@data))
save(file = "invitro_pseudotime.Rdata", df_invitro)

############################# ---------------------------- #############################
############################# Differential Expressed Genes #############################
############################# ---------------------------- #############################

## Load pre-built Seurat object
invitro <- readRDS("invitro.rds") %>% `DefaultAssay<-`(value = "RNA") %>% NormalizeData()
invitro <- invitro[,Idents(invitro) %in% c("D4_EC","D6_EC","HPC")]

invivo <- readRDS("invivo.rds") %>% SetIdent(value = .$celltype) %>% `DefaultAssay<-`(value = "RNA")
invivo <- invivo[,Idents(invivo) %in% c("AEC","HEC","HSC")]

## Set RUNX1+ EC at D4 as potential HECs
D4_RUNX1 <- subset(invitro, RUNX1>0 & Idents(invitro)=="D4_EC") %>% Cells()
invitro <- SetIdent(invitro, cells = D4_RUNX1, value = "RUNX1+D4_EC")

## In vitro DEGs set
DEG_invitro <- FindMarkers(invitro, ident.1 = "HPC", ident.2 = "RUNX1+D4_EC")
DEG_invitro <- DEG_invitro[order(DEG_invitro$avg_log2FC,decreasing = T),]
DEG_invitro <- DEG_invitro[DEG_invitro$p_val_adj<=0.0001,]

write.csv(DEG_invitro,"DEG_invitro.csv")

## In vivo DEGs set
DEG_invivo <- FindMarkers(invivo, ident.1 = "HSC", ident.2 = "HEC")
DEG_invivo <- DEG_invivo[order(DEG_invivo$avg_log2FC,decreasing = T),]
DEG_invivo <- DEG_invivo[DEG_invivo$p_val_adj<=0.0001,]

write.csv(DEG_invivo,"DEG_invivo.csv")

############################# ------------- #############################
############################# GO enrichment #############################
############################# ------------- #############################

DEG_invitro <- read.csv("DEG_invitro.csv",row.names = 1)
DEG_invivo <- read.csv("DEG_invivo.csv",row.names = 1)

## GO enrichment based on DEGs
GO_invitro_up <- enrichGO(gene= rownames(DEG_invitro[DEG_invitro$avg_log2FC>0,]),ont = "BP",
                  keyType = "SYMBOL",OrgDb=org.Hs.eg.db,readable= F)
GO_invitro_down <- enrichGO(gene= rownames(DEG_invitro[DEG_invitro$avg_log2FC<0,]),ont = "BP",
                    keyType = "SYMBOL",OrgDb=org.Hs.eg.db,readable= F)

GO_invivo_up <- enrichGO(gene= rownames(DEG_invivo[DEG_invivo$avg_log2FC>0,]),ont = "BP",
                          keyType = "SYMBOL",OrgDb=org.Hs.eg.db,readable= F)
GO_invivo_down <- enrichGO(gene= rownames(DEG_invivo[DEG_invivo$avg_log2FC<0,]),ont = "BP",
                            keyType = "SYMBOL",OrgDb=org.Hs.eg.db,readable= F)

## Visualization
df <- GO_invivo_up@result

ggplot(df[1:5,],aes(x=reorder(Description, -log10(p.adjust)), y=-log10(p.adjust))) +
  geom_segment(aes(xend=Description, yend=0)) +
  geom_point(size = 3, color="#bd0026")+ # "#08519c": down "#bd0026" : up
  coord_flip()+
  theme_bw()+
  xlab("")

############################# ------------------------ #############################
############################# Gene Expression Patterns #############################
############################# ------------------------ #############################

load("invitro_pseudotime.Rdata")
load("invivo_pseudotime.Rdata")

endo <- c("SOX17", "SOX7", "GATA3", "ETV2", "ERG", "ETS1")
hema <- c("SPI1", "RUNX1", "MYB", "GFI1B", "GATA2", "GFI1", "LYL1")

## Visualization
list1 <- lapply(endo,function(key) {
  tmp <- ggplot(df_invitro, aes_string("Pseudotime", key))+
    geom_smooth(method = "loess", show.legend = T,
                color = "#33a02c", size = 1.5, se = F)+
    theme_base()+
    ylab(key)+
    theme(axis.title.x = element_blank(),
          # axis.text.y = element_text(angle = 90),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  tmp})

list2 <- lapply(endo,function(key) {
  tmp <- ggplot(df_invivo, aes_string("Pseudotime", key))+
    geom_smooth(method = "loess", show.legend = T,
                color = "#33a02c", size = 1.5, se = F)+
    theme_base()+
    ylab(key)+
    theme(axis.title.x = element_blank(),
          # axis.text.y = element_text(angle = 90),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  tmp})

EC_pattern_invitro <- patchwork::wrap_plots(list1, ncol = 1)
EC_pattern_invivo <- patchwork::wrap_plots(list2, ncol = 1)

EC_pattern_invitro | EC_pattern_invivo

############################# ------------------------------------- #############################
############################# ChEA3 up-stream regulators enrichment #############################
############################# ------------------------------------- #############################

## Load pre-built DEGs set
DEG_invitro <- read.csv("DEG_invitro.csv",row.names = 1)
DEG_invivo <- read.csv("DEG_invivo.csv",row.names = 1)

## Set up ChEA3 R API
url = "https://maayanlab.cloud/chea3/api/enrich/"
encode <- "json"

payload_invitro <- list(query_name = "DEG_up_invitro",
                        gene_set = rownames(DEG_invitro[DEG_invitro$avg_log2FC>0,]))
payload_invivo <- list(query_name = "DEG_up_invivo",
                       gene_set = rownames(DEG_invivo[DEG_invivo$avg_log2FC>0,]))

## Grep results from ChEA3
ChEA3_invivo <- fromJSON(httr::content(POST(url = url,
                                            body = payload_invivo,
                                            encode = encode),
                                       as="text"))$`Integrated--meanRank`

ChEA3_invitro <- fromJSON(httr::content(POST(url = url,
                                             body = payload_invitro,
                                             encode = encode),
                                        as="text"))$`Integrated--meanRank`

## Summarize overlapping gene number
ChEA3_invitro$inter_gene <- ChEA3_invitro$Overlapping_Genes %>% strsplit("\\,") %>% sapply(length)
ChEA3_invivo$inter_gene <- ChEA3_invivo$Overlapping_Genes %>% strsplit("\\,") %>% sapply(length)

## Visualization

p1 <- ggplot(ChEA3_invivo[1:15,],
             aes(x=reorder(TF,-as.numeric(Rank)),
                 y=1/as.numeric(Score),
                 fill = inter_gene,text = TF))+
  geom_bar(stat = "identity")+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        legend.position = c(.8,.4),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_gradient2(low = "lightgrey",high = "#bd0026")+
  xlab("ChEA3 Rank")+
  ylab("1/Integrated Score")+
  coord_flip()+
  ggtitle("invivo")

p2 <- ggplot(ChEA3_invitro[1:15,],
             aes(x=reorder(TF,-as.numeric(Rank)),
                 y=1/as.numeric(Score),
                 fill = inter_gene))+
  geom_bar(stat = "identity")+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.title = element_blank(),
        legend.position = c(.8,.4),
        plot.title = element_text(hjust = 0.5))+
  scale_fill_gradient2(low = "lightgrey",high = "#bd0026")+
  xlab("ChEA3 Rank")+
  ylab("1/Integrated Score")+
  coord_flip()+
  ggtitle("invitro")

p1 | p2
############################# --------------------------------------- #############################
############################# Potential heterogeneity of in vitro HPC #############################
############################# --------------------------------------- #############################

## Load pre-built Seurat Object
invitro <- readRDS("invitro.rds") %>% `DefaultAssay<-`(value = "RNA") %>% NormalizeData()
HPC_invitro <- invitro[,Idents(invitro)=="HPC"] %>% DietSeurat(assays = "RNA")

SPI1_LYL1 <- Cells(subset(HPC_invitro, SPI1>0 & KLF1==0))
SPI1_KLF1 <- Cells(subset(HPC_invitro, SPI1>0 & KLF1>0))
other <- setdiff(Cells(HPC_invitro),c(SPI1_LYL1,SPI1_KLF1))

HPC_invitro <- SetIdent(HPC_invitro, cells = SPI1_LYL1, value = 'SPI1_LYL1')
HPC_invitro <- SetIdent(HPC_invitro, cells = SPI1_KLF1, value = 'SPI1_KLF1')
HPC_invitro <- SetIdent(HPC_invitro, cells = other, value = 'other')

Idents(HPC_invitro) <- factor(Idents(HPC_invitro), levels = c("SPI1_KLF1", "SPI1_LYL1", "other"))

## Calculate module score using lineage-sepecific gene sets
# GO:0030099 myeloid cell differentiation
# GO:0002320 lymphoid progenitor cell differentiation
# GO:0030218 erythrocyte differentiation

HPC_invitro <- AddModuleScore(HPC_invitro, name = "mye_score",
                              features = list(readLines("GO_0030099.txt")))
HPC_invitro <- AddModuleScore(HPC_invitro, name = "lym_score",
                              features = list(readLines("GO_0002320.txt")))
HPC_invitro <- AddModuleScore(HPC_invitro, name = "ery_score",
                              features = list(readLines("GO_0030218.txt")))
df <- data.frame(celltype = Idents(HPC_invitro),
                 mye_score = HPC_invitro$mye_score1,
                 lym_score = HPC_invitro$lym_score1,
                 ery_score = HPC_invitro$ery_score1)

plot_list <- lapply(c("mye_score","lym_score","ery_score"),function(key){
  ggplot(df, aes(x = as.factor(celltype), y = .data[[key]], fill = celltype))+
    geom_violin()+
    geom_signif(comparison = list(c("SPI1_KLF1", "SPI1_LYL1")))+
    theme_classic()+
    NoLegend()+
    theme(axis.title = element_blank())
})

patchwork::wrap_plots(plot_list, ncol = 1)


## Unsupervised clustering of HPC

HPC_invitro <- HPC_invitro %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:50) %>% 
  FindNeighbors()

HPC_invitro <- FindClusters(HPC_invitro, resolution = seq(0.1,1,0.1))

Idents(HPC_invitro) <- HPC_invitro$RNA_snn_res.0.2

## Feature genes expression
HPC_invitro <- ScaleData(HPC_invitro,features = rownames(HPC_invitro))
marker <- FindAllMarkers(HPC_invitro) %>% split(.$cluster)
marker <- lapply(marker,function(x){
  x <- x[order(x$avg_log2FC,decreasing = T),]
  x <- x[x$p_val_adj<=0.0001,]
  x
})

feature_gene <- sapply(marker,function(x){x$gene[1:20]}) %>% as.vector()

p1 <- DoHeatmap(HPC_invitro, features = feature_gene)+
  theme(legend.position = "bottom")+
  guides(color = "none")

## GO terms enriched by Feature genes
marker_GO <- lapply(marker,function(x){
  tmp <- x$gene[x$avg_log2FC>0]
  tmp_GO <- enrichGO(tmp[1:20], OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL", readable= F)
  tmp_GO
}) %>% set_names(names(marker))

# Visualization
col <- c("#F8766D","#00BA38","#619CFF")

p_list <- lapply(1:3,function(x){
  ggplot(marker_GO[[x]]@result[1:5,],aes(x=reorder(Description,-log10(pvalue)),y=-log10(pvalue)))+
    geom_bar(stat = "identity",color = "black",fill = col[x])+
    geom_text(aes(label = Description),position = position_fill(vjust = 2))+
    coord_flip()+
    theme_classic()+
    ggtitle("")+
    theme(panel.background = element_rect(fill = "#F0F0F0"),axis.line.y = element_blank(),
          legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_blank(),
          axis.ticks.y = element_blank())})

patchwork::wrap_plots(p1, patchwork::wrap_plots(p_list,ncol = 1),ncol = 2)


## Marker genes expression
marker <- c("BIRC5", "TPX2", "CENPE", "MKI67", "CDK1","TMEM14C", "KLF1", "GATA1","TAL1",
            "LYZ", "MEIS1", "HES1", "LCP1", "PLCG2", "HHEX", "RBPJ", "ZBTB1", "SOX4")

DotPlot(HPC_invitro, features = marker)+
  coord_flip()+
  scale_color_gradient2(low = "#6baed6",
                        mid = "white",
                        midpoint = 0,
                        high = "#bd0026")+
  theme(axis.title = element_blank(),
        axis.text.x.bottom = element_text(angle = 45))+
  scale_x_discrete(limits = rev(marker))

## Add cell cycle core to sub-clusters
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

HPC_invitro <- CellCycleScoring(HPC_invitro, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

DimPlot(HPC_invitro, group.by = "Phase",label = T,label.size = 5)+
  NoLegend()+
  theme(plot.title = element_blank())

## Rename identity based on results above
HPC_invitro <- RenameIdents(HPC_invitro, "0"="E/M-HPC", "2"="L/M-HPC", "1"="P-HPC")

DimPlot(HPC_invitro, pt.size = 3,label = T,label.size = 6)+
  NoLegend()+
  theme(plot.title = element_blank())

saveRDS(HPC_invitro,"HPC_invitro.rds")


############################# -------------------------------------- #############################
############################# Potential heterogeneity of in vivo HSC #############################
############################# -------------------------------------- #############################

## Load pre-built Seurat Object

invivo <- readRDS("invivo.rds") %>% `DefaultAssay<-`(value = "RNA")
HSC_invivo <- invivo[,invivo$celltype=="HSC"] %>% DietSeurat(assays = "RNA")

## Split HSC into sub-clusters based on monocle2 trajectory
load("monocle2_invivo.Rdata")
cds <- cds[,colnames(HSC_invivo)]
HSC_invivo$state <- cds$State
HSC_invivo$celltype_new <- HSC_invivo$celltype
HSC_invivo$celltype_new[HSC_invivo$celltype_new=="HSC" & HSC_invivo$state!=5] <- "HSC_1"
HSC_invivo$celltype_new[HSC_invivo$celltype_new=="HSC" & HSC_invivo$state==5] <- "HSC_2"

Idents(HSC_invivo) <- HSC_invivo$celltype_new
## Calculate module score using lineage-sepecific gene sets
# GO:0030099 myeloid cell differentiation
# GO:0030097 lymphocyte differentiation

HSC_invivo <- AddModuleScore(HSC_invivo, name = "mye_score",
                              features = list(readLines("GO_0030099.txt")))
HSC_invivo <- AddModuleScore(HSC_invivo, name = "lym_score",
                              features = list(readLines("GO_0030097.txt")))

df <- data.frame(celltype = Idents(HSC_invivo),
                 mye_score = HSC_invivo$mye_score1,
                 lym_score = HSC_invivo$lym_score1)

plot_list <- lapply(c("mye_score","lym_score"),function(key){
  ggplot(df, aes(x = as.factor(celltype), y = .data[[key]], fill = celltype)) +
    geom_violin()+
    geom_signif(comparison = list(c("HSC_1", "HSC_2")))+
    theme_classic()+
    NoLegend()+
    theme(axis.title = element_blank())
})

patchwork::wrap_plots(plot_list, ncol = 1)


## Unsupervised clustering of HSC

HSC_invivo <- HSC_invivo %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:50) %>% 
  FindNeighbors()

for(i in seq(0.1,1,0.1)){HSC_invivo <- FindClusters(HSC_invivo, resolution = i)}

Idents(HSC_invivo) <- HSC_invivo$RNA_snn_res.0.4

## Feature genes expression
HSC_invivo <- ScaleData(HSC_invivo,features = rownames(HSC_invivo))

HSC_0 <- FindMarkers(HSC_invivo,ident.1 = "0")
HSC_0 <- HSC_0[order(HSC_0$avg_log2FC,decreasing = T),]
HSC_0 <- HSC_0[HSC_0$p_val_adj<=0.0001,]
HSC_0_GO <- enrichGO(rownames(HSC_0)[1:20],
                     OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL", readable= F)

HSC_1 <- FindMarkers(HSC_invivo,ident.1 = "1",logfc.threshold = 0.2)
HSC_1 <- HSC_1[order(HSC_1$avg_log2FC,decreasing = T),]
HSC_1 <- HSC_1[HSC_1$p_val_adj<=0.05,]
HSC_1_GO <- enrichGO(rownames(HSC_1)[HSC_1$avg_log2FC>0],
                     OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL", readable= F)

HSC_2 <- FindMarkers(HSC_invivo,ident.1 = "2")
HSC_2 <- HSC_2[order(HSC_2$avg_log2FC,decreasing = T),]
HSC_2 <- HSC_2[HSC_2$p_val_adj<=0.0001,]
HSC_2_GO <- enrichGO(rownames(HSC_2)[1:20],
                     OrgDb = org.Hs.eg.db, ont = "BP", keyType = "SYMBOL", readable= F)

feature_gene <- c(rownames(HSC_0)[1:20],
                  rownames(HSC_1)[HSC_1$avg_log2FC>0],
                  rownames(HSC_2)[1:20]) %>% unique()

p1 <- DoHeatmap(HSC_invivo, features = feature_gene)+
  theme(legend.position = "bottom")+
  guides(color = "none")

## GO terms enriched by Feature genes
marker_GO <- list(HSC_0_GO,HSC_1_GO,HSC_2_GO)

# Visualization
col <- c("#F8766D","#00BA38","#619CFF")
p_list <- lapply(1:3,function(x){
  ggplot(marker_GO[[x]]@result[1:5,],aes(x=reorder(Description,-log10(pvalue)),y=-log10(pvalue)))+
    geom_bar(stat = "identity",color = "black",fill = col[x])+
    geom_text(aes(label = Description),position = position_fill(vjust = 2))+
    coord_flip()+
    theme_classic()+
    ggtitle("")+
    theme(panel.background = element_rect(fill = "#F0F0F0"),axis.line.y = element_blank(),
          legend.title = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_blank(),
          axis.ticks.y = element_blank())})

wrap_plots(p1, wrap_plots(p_list,ncol = 1),ncol = 2)

## Marker genes expression
marker <- c("BIRC5", "TPX2", "CENPE", "MKI67", "CDK1", "TAL1","CEBPB", "CSF1R",
            "CEBPA", "IKZF1", "WWP1", "ZBTB1", "RBPJ", "LYL1", "PLCG2", "RAC2")

DotPlot(HSC_invivo, features = marker)+
  coord_flip()+
  scale_color_gradient2(low = "#6baed6",
                        mid = "white",
                        midpoint = 0,
                        high = "#bd0026")+
  theme(axis.title = element_blank(),
        axis.text.x.bottom = element_text(angle = 45))+
  scale_x_discrete(limits = rev(marker))

## Add cell cycle core to sub-clusters
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

HSC_invivo <- CellCycleScoring(HSC_invivo, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

DimPlot(HSC_invivo, group.by = "Phase",label = T,label.size = 5)+
  NoLegend()+
  theme(plot.title = element_blank())

## Rename identity based on results above
HSC_invivo <- RenameIdents(HSC_invivo, "1"="M-HSC", "2"="L-HSC", "0"="P-HPC")

DimPlot(HSC_invivo, pt.size = 3,label = T,label.size = 6)+
  NoLegend()+
  theme(plot.title = element_blank())

saveRDS(HPC_invitro,"HSC_invivo.rds")