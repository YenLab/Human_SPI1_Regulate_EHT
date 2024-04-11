package <- c("airway",
             "BSgenome.Hsapiens.UCSC.hg38",
             "ChIPseeker", "circlize", "cluster", "clusterProfiler", "ComplexHeatmap",
             "DESeq2", "DiffBind", "DOSE", "dplyr",
             "edgeR", "EnhancedVolcano",
             "factoextra", "FindIT2",
             "genefilter", "GenomicFeatures", "ggimage", "ggplot2", "ggplotify", "ggrepel", "ggthemes", "ggupset",
             "JASPAR2020",
             "magrittr", "Mfuzz", "motifmatchr",
             "org.Hs.eg.db",
             "parallel", "pheatmap",
             "RColorBrewer", "reshape2",
             "tidyverse")
invisible(lapply(package, library, character.only = T))

setwd("./Input_RNA&ATAC-seq")

############################# -------------------------- #############################
#############################           RNA-seq          #############################
############################# -------------------------- #############################

#### Quality Control & differential expression analysis

PSC1 <- read.table("PSC_rep1number.txt", sep = "\t", header = F, col.names = c("gene","PSC1"))
PSC2 <- read.table("PSC_rep2number.txt", sep = "\t", header = F, col.names = c("gene","PSC2"))
VME1 <- read.table("VME_rep1number.txt", sep = "\t", header = F, col.names = c("gene","VME1"))
VME2 <- read.table("VME_rep2number.txt", sep = "\t", header = F, col.names = c("gene","VME2"))
EPC1 <- read.table("EPC_rep1number.txt", sep = "\t", header = F, col.names = c("gene","EPC1"))
EPC2 <- read.table("EPC_rep2number.txt", sep = "\t", header = F, col.names = c("gene","EPC2"))
HPC1 <- read.table("HPC_rep1number.txt", sep = "\t", header = F, col.names = c("gene","HPC1"))
HPC2 <- read.table("HPC_rep2number.txt", sep = "\t", header = F, col.names = c("gene","HPC2"))
raw_count1 <- merge(merge(PSC1, PSC2, by="gene"), merge(VME1, VME2, by="gene"))
raw_count2 <- merge(merge(EPC1, EPC2, by="gene"), merge(HPC1, HPC2, by="gene"))
raw_count <- merge(raw_count1, raw_count2, by = "gene")
raw_count <- raw_count[-1:-5,]
rownames(raw_count)<-raw_count[,1]
raw_count<-raw_count[,-1]
condition <- factor(c(rep("PSC",2),rep("VME",2),rep("EPC",2),rep("HPC",2)), 
                    levels = c("PSC","VME","EPC","HPC"))
colData <- data.frame(row.names=colnames(raw_count), condition)

dds <- DESeqDataSetFromMatrix(raw_count, colData, design= ~ condition)
dds <- dds[ rowSums(counts(dds))>10, ]
dds <- DESeq(dds)
norcount <- counts(dds,normalized=T)

res1 = results(dds, contrast = c("condition", "VME", "PSC"))
res2 = results(dds, contrast = c("condition", "EPC", "VME"))
res3 = results(dds, contrast = c("condition", "HPC", "EPC"))
table(res1$padj<0.05)
diff_gene1 <- subset(res1, padj < 0.05)
diff_gene2 <- subset(res2, padj < 0.05)
diff_gene3 <- subset(res3, padj < 0.05)
gene1 <- rownames(diff_gene1)
gene2 <- rownames(diff_gene2)
gene3 <- rownames(diff_gene3)

#### comparison between up-/down-regulated genes among cell types (Fig S2.B)
mydata <- c(1794,1742,1708,1905,2485,2774)
Stage <- c("PSC to VME","PSC to VME","VME to EPC","VME to EPC","EPC to HPC","EPC to HPC")
regulate <- c("up","down","up","down","up","down")
updown <- data.frame(Stage,regulate,mydata)
updown$Stage = factor(updown$Stage, levels = c("PSC to VME","VME to EPC","EPC to HPC"))
updown$regulate = factor(updown$regulate, levels = c("up","down"))
ggplot(updown,aes(Stage,mydata,fill=regulate))+geom_bar(stat="identity",position="dodge")+
  ggtitle("Number of differentially expressed genes")+
  theme_wsj()+
  scale_fill_wsj()


allgene <- union(gene1,gene2)
allgene <- union(allgene,gene3)

norcount <- as.data.frame(norcount)
norcount$PSC <- (norcount$PSC1+norcount$PSC2)/2
norcount$VME <- (norcount$VME1+norcount$VME2)/2
norcount$EPC <- (norcount$EPC1+norcount$EPC2)/2
norcount$HPC <- (norcount$HPC1+norcount$HPC2)/2
mergecount <- norcount[,9:12]
mergecount <- mergecount[allgene,]

#### Select TF from DEGs
TFlist <- read.csv("humanTFlist.csv",header=F)
TF <- bitr(TFlist$V1, fromType = "SYMBOL", toType = "ENSEMBL",OrgDb = org.Hs.eg.db)
TFoverlap <- intersect(TF$ENSEMBL,rownames(mergecount))

dynamicTFcount <- mergecount[TFoverlap,]
dynamicTFcount$ENSEMBL <- rownames(dynamicTFcount)
dynamicTFcount <- merge(dynamicTFcount,TF,by="ENSEMBL")
rownames(dynamicTFcount) <- dynamicTFcount$SYMBOL
dynamicTFcount <- dynamicTFcount[,-1]
dynamicTFcount <- dynamicTFcount[,-5]

quTF <- (dynamicTFcount - rowMeans(dynamicTFcount))/rowSds(dynamicTFcount)

#### Clustering of RNA-seq DE-TF using K-means (Fig S2.C)
fviz_nbclust(quTF, kmeans, method="wss")+geom_vline(xintercept = 4, linetype=2)
gap_stat <- clusGap(quTF, FUN = kmeans, nstart = 25, K.max = 10, B = 500)
# fviz_gap_stat(gap_stat)+geom_vline(xintercept = 7, linetype=2)
km <- kmeans(quTF,4,nstart = 25)
# fviz_cluster(object = km,
#              quTF,
#              ellipse.type = "euclid",
#              star.plot=T,
#              repel = T,
#              geom=c("point","text"),
#              palette='jco', main="",
#              ggtheme = theme_minimal())+
#   theme(axis.title = element_blank())
qucluster <- cbind(quTF,cluster=km$cluster)
qucluster <- data.frame(qucluster)
qucluster <- qucluster[order(qucluster$cluster,decreasing = F),]

#### visualization (Fig 2.A)
TFheatmap <- qucluster[,-5]
TFheatmap <- as.matrix(TFheatmap)
pheatmap(TFheatmap, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(250),
         main = "dynamic TFs", 
         fontsize = 25, 
         border_color = "NA", 
         cluster_cols = F,
         cluster_rows = F,
         gaps_row = c(100,271,442),
         show_rownames = F)

#### GO terms enrichment (Fig S2.E)
b <- rownames(qucluster[which(qucluster$cluster=="4"),])
C1GO <- enrichGO(gene       = b,
                 OrgDb      = org.Hs.eg.db,
                 keyType    = 'SYMBOL',
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
C1GO <- clusterProfiler::simplify(C1GO, cutoff = 0.7, by = "p.adjust", select_fun = min)
barplot(C1GO, showCategory = 10, title="C1 GO ", font.size=20)


############################# -------------------------- #############################
#############################          ATAC-seq          #############################
############################# -------------------------- #############################

#### Construct dynamic peak matrix
DBA_object <- dba(sampleSheet = "sample_sheet.csv")
Rep_consensus <- dba.peakset(DBA_object, consensus = DBA_TREATMENT, minOverlap = 2)
dba.plotVenn(DBA_object,DBA_object$masks$four)
Rep_consensus <- dba(Rep_consensus, mask=Rep_consensus$masks$Consensus,minOverlap = 1)
consensus_peaks <- dba.peakset(Rep_consensus, bRetrieve=TRUE)

PSC_consensus <- dba(Rep_consensus,mask = Rep_consensus$masks$one, minOverlap = 2)
PSC_peaks <- dba.peakset(PSC_consensus, bRetrieve=TRUE)
VME_consensus <- dba(Rep_consensus,mask = Rep_consensus$masks$two, minOverlap = 2)
VME_peaks <- dba.peakset(VME_consensus, bRetrieve=TRUE)
EPC_consensus <- dba(Rep_consensus,mask = Rep_consensus$masks$three, minOverlap = 2)
EPC_peaks <- dba.peakset(EPC_consensus, bRetrieve=TRUE)
HPC_consensus <- dba(Rep_consensus,mask = Rep_consensus$masks$four, minOverlap = 2)
HPC_peaks <- dba.peakset(HPC_consensus, bRetrieve=TRUE)

consensus_count <- dba.count(DBA_object,bUseSummarizeOverlaps = TRUE,peaks=consensus_peaks)
diffpeak <- dba.contrast(consensus_count, categories=DBA_TREATMENT,minMembers = 2)
diffpeak <- dba.analyze(diffpeak, method=DBA_ALL_METHODS) ####FDR<=0.5
dba.show(diffpeak,bContrasts = T)
dba.plotVenn(diffpeak,contrast=1,method=DBA_ALL_METHODS)
dba.plotVenn(diffpeak,contrast=4,method=DBA_ALL_METHODS)
dba.plotVenn(diffpeak,contrast=6,method=DBA_ALL_METHODS)
dba.plotMA(diffpeak,contrast=1,method=DBA_ALL_METHODS)
comp1deseq <- dba.report(diffpeak,method=DBA_DESEQ2,contrast = 1)
PSCVMEdiffpeak <- as.data.frame(comp1deseq)
comp2deseq <- dba.report(diffpeak,method=DBA_DESEQ2,contrast = 4)
VMEEPCdiffpeak <- as.data.frame(comp2deseq)
comp3deseq <- dba.report(diffpeak,method=DBA_DESEQ2,contrast = 6)
EPCHPCdiffpeak <- as.data.frame(comp3deseq)
norcountsmatrix <- as.data.frame(dba.peakset(diffpeak, bRetrieve=TRUE))
dynamicpeak <- union(rownames(PSCVMEdiffpeak),rownames(VMEEPCdiffpeak))
dynamicpeak <- union(dynamicpeak,rownames(EPCHPCdiffpeak))
dynamicpeakmatrix <- norcountsmatrix[dynamicpeak,]

#### Summarize diffpeak number in Fig.S2G
table(PSCVMEdiffpeak$Fold>0)
table(VMEEPCdiffpeak$Fold>0)
table(EPCHPCdiffpeak$Fold>0)

df <- data.frame(a=c(table(PSCVMEdiffpeak$Fold>0),
                     table(VMEEPCdiffpeak$Fold>0),
                     table(EPCHPCdiffpeak$Fold>0)) %>% as.numeric(),
                 b=paste0("G",rep(1:3,each = 2)),
                 c=rep(c("C1","C2"),3))

ggplot(df,aes(b,a,group=c,fill = c))+
  geom_bar(stat = "identity",position = 'dodge')+
  theme_classic()

peakset <- dynamicpeakmatrix[,c(1:3)] %>% 
  mutate(feature_id=paste0("Peak",1:nrow(.)),
         score=0)
peakset$seqnames <- paste0("chr",peakset$seqnames)

dynamicpeakmatrix <- dynamicpeakmatrix[,-c(1:5)] %>% 
  set_colnames(paste0(rep(c("PSC","VME","EPC","HPC"),each = 2),c(1,2))) %>% 
  set_rownames(paste0("Peak",1:nrow(.)))

write.csv(peakset,"ATACdynamic.csv")
write.csv(dynamicpeakmatrix,"dynamicpeakmatrix.csv")


#### Calculate TF activity (Fig.2D-E)
ATAC_norcount <- read.csv("dynamicpeakmatrix.csv",row.names = 1)
type <- factor(c(rep("PSC",2),rep("VME",2),rep("EPC",2),rep("HPC",2)), 
               levels = c("PSC","VME","EPC","HPC"))
ATAC_colData <- data.frame(row.names=colnames(ATAC_norcount), type)
ATAC_normCount_merge <- integrate_replicates(ATAC_norcount, ATAC_colData)
ATACdynamic <- readPeakFile("ATACdynamic.csv",header=T, sep=",")

getJasparMotifs <- function(species = "Homo sapiens", 
                            collection = "CORE", ...) {
  opts <- list()
  opts["species"] <- species
  opts["collection"] <- collection
  opts[["all_versions"]] <- TRUE
  opts <- c(opts, list(...))
  out <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020, opts)
  if (!isTRUE(all.equal(TFBSTools::name(out), names(out)))) 
    names(out) <- paste(names(out), TFBSTools::name(out), sep = "_")
  return(out)
}
motifs_JASPAR2020 <- getJasparMotifs()

## filter top20000 peak by sd for motif matching
apply(ATAC_normCount_merge, 1, sd) %>% 
  sort(decreasing = TRUE) %>% 
  names() %>% 
  "["(1:20000) -> peak_name
test <- as.data.frame(ATACdynamic)
test <- filter(test,feature_id %in% peak_name)

# write.csv(test,"top20000peak.csv")

ATACdynamictop <- readPeakFile("top20000peak.csv",header=T, sep=",")
motif_ix <- matchMotifs(motifs_JASPAR2020, 
                        ATACdynamictop, 
                        genome = "BSgenome.Hsapiens.UCSC.hg38", 
                        out = "positions", 
                        bg = "genome")
motif_GR <- unlist(motif_ix)
motif_GR$TF_id <- gsub("(MA[0-9]{4}\\.[0-9]_)", "", names(motif_GR))
set.seed(20220811)

## tf influence trend between samples
findIT_MARA(input_feature_id = peak_name,
            peak_GR = ATACdynamictop,
            peakScoreMt = ATAC_norcount[peak_name,],
            TF_GR_database = motif_GR) -> result_MARA

MARA_mt <- as.matrix(result_MARA[, -1])
rownames(MARA_mt) <- result_MARA$TF_id

integrate_replicates(mt = MARA_mt,
                     colData = ATAC_colData,
                     type = "rank_zscore") -> MARA_mt_merge
MARA_mt_merge <- as.data.frame(MARA_mt_merge)
MARA_mt_merge$a <- MARA_mt_merge$HPC - MARA_mt_merge$EPC
MARA_mt_merge$a <- abs(MARA_mt_merge$a)
MARA_mt_merge <- MARA_mt_merge[order(MARA_mt_merge$a,decreasing = T),]
MARA_mt_merge <- MARA_mt_merge[,1:4]
MARA_mt_merge <- as.matrix(MARA_mt_merge)

test <- (MARA_mt_merge - rowMeans(MARA_mt_merge))/rowSds(MARA_mt_merge)
pheatmap(test[1:50,], 
         color = colorRampPalette(c("#00473C", "white", "#5C3619"))(250),
         main = "TF activity", 
         fontsize = 13.5, 
         border_color = "NA", 
         cluster_cols = F,
         cluster_rows = T,angle_col = "0",
         show_rownames = T,
)

#### Clustering of ATAC-seq dynamic peak using K-means (Fig.2B-C)
Peakscale <- (ATAC_normCount_merge - rowMeans(ATAC_normCount_merge))/rowSds(ATAC_normCount_merge)
fviz_nbclust(Peakscale, kmeans, method="wss")+geom_vline(xintercept = 4, linetype=2)
km <- kmeans(Peakscale,4,nstart = 25)
cluster <- cbind(Peakscale,cluster=km$cluster)
cluster <- data.frame(cluster)
cluster <- cluster[order(cluster$cluster,decreasing = F),]
table(cluster$cluster) ###see Fig.2B


## Calculate correlation between gene expression and TF activity (Fig.S2H)
allpeakavtivity <- (MARA_mt_merge - rowMeans(MARA_mt_merge))/rowSds(MARA_mt_merge)
dynamicTF <- read.csv("E:/OneDrive/share_with_kyq/分析重复/dynamicTFcount.csv")
rownames(dynamicTF) <- dynamicTF[,1]
dynamicTF <- dynamicTF[,-1]
tfoverlap <- intersect(rownames(dynamicTF),rownames(allpeakavtivity))
allpeakavtivityoverlap <- allpeakavtivity[tfoverlap,]
dynamicTFoverlap <- dynamicTF[tfoverlap,]
allpeakavtivityoverlap <- as.matrix(allpeakavtivityoverlap)
dynamicTFoverlap <- as.matrix(dynamicTFoverlap)
dynamicTFoverlap <- (dynamicTFoverlap - rowMeans(dynamicTFoverlap))/rowSds(dynamicTFoverlap)
cor_results <- list()
for (gene in rownames(allpeakavtivityoverlap)) {
  gene_data_mat1 <- allpeakavtivityoverlap[gene,]
  gene_data_mat2 <- dynamicTFoverlap[gene,]
  correlation <- cor(gene_data_mat1, gene_data_mat2)
  cor_results[[gene]] <- correlation
}
# write.csv(cor_results,"Cor_TFAct_GeneExp.csv")

a <- read.csv("Cor_TFAct_GeneExp.csv",row.names = 1)
a <- a %>% t() %>% as.data.frame()
a <- a[order(a[,1],decreasing = T),,drop=F] %>% 
  mutate(rank=1:nrow(.),gene_id=rownames(.)) %>% 
  set_colnames(c("cor","rank","gene_id"))
b <- filter(a,gene_id %in% c("SPI1","RUNX1","POU5F1","TBXT","ETS1"))
ggplot(a,aes(x=rank,y=cor))+
  geom_point(color="red4",size=2)+
  theme_few()+
  labs(x="Rank",y="Correlation")+
  geom_text_repel(data = b, mapping = aes(label=gene_id))


#### peak annotation (Fig.S2I)
txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.104.gtf.gz",format="gtf")
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

## bed files output from peak calling based on merged diffbind consensus peak
files = list(PSC="E:/OneDrive/share_with_kyq/分析重复/peak bed/PSCintersect.bed",
             VME="E:/OneDrive/share_with_kyq/分析重复/peak bed/VMEintersect.bed",
             EPC="E:/OneDrive/share_with_kyq/分析重复/peak bed/EPCintersect.bed",
             HPC="E:/OneDrive/share_with_kyq/分析重复/peak bed/HPCintersect.bed")
options(ChIPseeker.ignore_1st_exon = T)
options(ChIPseeker.ignore_1st_intron = T)
options(ChIPseeker.ignore_downstream = T)
options(ChIPseeker.ignore_promoter_subcategory = T)
options(ChIPseeker.ignore_utr_subcategory = T)
peakAnnoList <- lapply(files, 
                       annotatePeak,
                       TxDb=txdb,
                       tssRegion=c(-3000, 3000),annoDb="org.Hs.eg.db")
plotAnnoBar(peakAnnoList)+coord_flip()


#### visualization of PECA2 output
a <- read.table("HPC_network.txt",sep = "\t",header = T)[,1:4]
a <- a[a$TF=="SPI1",]

ggplot()+
  geom_bar(data = a,aes(x = Score,y = reorder(TG,Score)),stat='identity')+
  geom_text()







