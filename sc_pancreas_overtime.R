load("~/Desktop/A1_outs/Sample_expr_FastMNN_Pancreas.Robj")

library(Seurat)
library(ggplot2)
library(dplyr)

Idents(Sample_expr_FastMNN)='RNA_snn_res.0.5'
clMono<-WhichCells(Sample_expr_FastMNN,idents = 9)
names<-rep("Classical Monocytes",length(clMono))
names(names)=clMono
ClMono<-names
rm(names)

Macro<-WhichCells(Sample_expr_FastMNN,idents = c(0,10,14))
names<-rep("Macrophages",length(Macro))
names(names)=Macro
Macro<-names
rm(names)

mDC<-WhichCells(Sample_expr_FastMNN,idents = 16)
names<-rep("mDCs",length(mDC))
names(names)=mDC
mDC<-names
rm(names)

cDC<-WhichCells(Sample_expr_FastMNN,idents = 12)
names<-rep("cDCs",length(cDC))
names(names)=cDC
cDC<-names
rm(names)

pDC<-WhichCells(Sample_expr_FastMNN,idents = 20)
names<-rep("pDCs",length(pDC))
names(names)=pDC
pDC<-names
rm(names)

ImmNeu<-WhichCells(Sample_expr_FastMNN,idents = 22)
names<-rep("Immature Neutrophiles",length(ImmNeu))
names(names)=ImmNeu
ImmNeu<-names
rm(names)


MatNeu<-WhichCells(Sample_expr_FastMNN,idents = 3)
names<-rep("Mature Neutrophiles",length(MatNeu))
names(names)=MatNeu
MatNeu<-names
rm(names)
NK<-WhichCells(Sample_expr_FastMNN,idents = 17)
names<-rep("NK",length(NK))
names(names)=NK
NK<-names
rm(names)
Tcells<-WhichCells(Sample_expr_FastMNN,idents = c(7,5,19))
names<-rep("T cells",length(Tcells))
names(names)=Tcells
Tcells<-names
rm(names)
Bcells<-WhichCells(Sample_expr_FastMNN,idents = c(4,23))
names<-rep("B cells",length(Bcells))
names(names)=Bcells
Bcells<-names
rm(names)
ImmB<-WhichCells(Sample_expr_FastMNN,idents = 18)
names<-rep("Immature B cells",length(ImmB))
names(names)=ImmB
ImmB<-names
rm(names)
Acinar<-WhichCells(Sample_expr_FastMNN,idents = 13)
names<-rep("Acinar cells",length(Acinar))
names(names)=Acinar
Acinar<-names
rm(names)
End<-WhichCells(Sample_expr_FastMNN,idents = 11)
names<-rep("Endothelial cells",length(End))
names(names)=End
End<-names
rm(names)
Epith_Tumor<-WhichCells(Sample_expr_FastMNN,idents = c(1,2,6,15))
names<-rep("Epith_Tumor",length(Epith_Tumor))
names(names)=Epith_Tumor
Epith_Tumor<-names
rm(names)


Idents(Sample_expr_FastMNN)='RNA_snn_res.0.5'
CAF<-WhichCells(Sample_expr_FastMNN,idents = c(21,8))
Idents(Sample_expr_FastMNN)='RNA_snn_res.0.6'
CAF1<-WhichCells(Sample_expr_FastMNN,idents = 15)
CAFFin<-unique(c(CAF,CAF1,CAF2))
names<-rep("CAF",length(CAFFin))
names(names)=CAFFin
CAFFin<-names
rm(names)

Idents(Sample_expr_FastMNN)='RNA_snn_res.0.6'
Cancer<-WhichCells(Sample_expr_FastMNN,idents = c(0,9,5,7))
names<-rep("Cancer",length(Cancer))
names(names)=Cancer
Cancer<-names
rm(names)

Idents(Sample_expr_FastMNN)='RNA_snn_res.0.6'
Ductal<-WhichCells(Sample_expr_FastMNN,idents = 18)
names<-rep("Ductal",length(Ductal))
names(names)=Ductal
Ductal<-names
rm(names)



annotation<-c(ClMono,Macro,mDC,cDC,pDC,ImmNeu,MatNeu,NK,Tcells,Bcells,ImmB,Acinar,End,CAFFin,Cancer,Ductal)
Sample_expr_FastMNN<-AddMetaData(Sample_expr_FastMNN,metadata = annotation,col.name = 'Annotation')
Idents(Sample_expr_FastMNN)='Annotation'
CAF2<-setdiff(rownames(Sample_expr_FastMNN@meta.data),names(annotation))
#redo and call all CAF

samplebyorig<-SplitObject(Sample_expr_FastMNN,split.by = 'orig.ident')

td10<-samplebyorig$Tumor_d10@assays$RNA@counts
td20<-samplebyorig$Tumor_d20@assays$RNA@counts
td30<-samplebyorig$Tumor_d30@assays$RNA@counts



header<-unique(sampleannot[,2])

sampleannot2<-data.frame(Cell=rownames(samplebyorig$Tumor_d30@meta.data),Type=samplebyorig$Tumor_d30@meta.data$Annotation,orig=samplebyorig$Tumor_d30@meta.data$orig.ident)

sampleannot1<-sampleannot2%>%
  filter(Type==c('Acinar cells','Endothelial cells','CAF','Epith_Tumor'))

sampleannot1[,2]=gsub('Epith_Tumor','malignant',sampleannot1[,2])

sampleannot1<-sampleannot1[,1:2]
names(sampleannot1)=NULL
sampleannot1<-as.data.frame(sampleannot1)

write.table(sampleannot1,row.names = F,col.names = F,file = 'sampleannot_d30.txt',sep='\t')
header<-c("Acinar cells","Endothelial cells","CAF")

infercnvSample_check<-CreateInfercnvObject(raw_counts_matrix = Sample_expr_FastMNN@assays$RNA@counts,
                                     gene_order_file = '/Users/leonardi.carlo/Desktop/A1_outs/infercnv/gene_ordering_infercnv.txt',
                                     annotations_file = '/Users/leonardi.carlo/Desktop/A1_outs/infercnv/sampleannot.txt',
                                     ref_group_names = header)

infercnvSample_check = infercnv::run(infercnvSample_check,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="/Users/leonardi.carlo/Desktop/A1_outs/infercnv/InferCNV_allTumors",
                             cluster_by_groups=TRUE,
                             denoise=TRUE,
                             HMM=TRUE)

library(NGCHM)
ngchm(infercnv_obj          = infercnvSample_check,
      out_dir              = "/Users/leonardi.carlo/Desktop/A1_outs/infercnv/InferCNV_allTumors",
      path_to_shaidyMapGen = '/Users/leonardi.carlo/Downloads/ShaidyMapGen.jar')

#infercnv sample-specific

infercnvSample_d10<-CreateInfercnvObject(raw_counts_matrix = samplebyorig$Tumor_d10@assays$RNA@counts,
                                           gene_order_file = '/Users/leonardi.carlo/Desktop/A1_outs/infercnv/gene_ordering_infercnv.txt',
                                           annotations_file = '/Users/leonardi.carlo/Desktop/A1_outs/infercnv/sampleannot_d10.txt',
                                           ref_group_names = header)

infercnvSample_d10 = infercnv::run(infercnvSample_d10,
                                     cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                     out_dir="/Users/leonardi.carlo/Desktop/A1_outs/infercnv/InferCNV_d10",
                                     cluster_by_groups=TRUE,
                                     denoise=TRUE,
                                     HMM=TRUE)

infercnvSample_d20<-CreateInfercnvObject(raw_counts_matrix = samplebyorig$Tumor_d20@assays$RNA@counts,
                                         gene_order_file = '/Users/leonardi.carlo/Desktop/A1_outs/infercnv/gene_ordering_infercnv.txt',
                                         annotations_file = '/Users/leonardi.carlo/Desktop/A1_outs/infercnv/sampleannot_d20.txt',
                                         ref_group_names = header)
infercnvSample_d20 = infercnv::run(infercnvSample_d20,
                                   cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                   out_dir="/Users/leonardi.carlo/Desktop/A1_outs/infercnv/InferCNV_d20",
                                   cluster_by_groups=TRUE,
                                   denoise=TRUE,
                                   HMM=TRUE)

infercnvSample_d30<-CreateInfercnvObject(raw_counts_matrix = samplebyorig$Tumor_d30@assays$RNA@counts,
                                         gene_order_file = '/Users/leonardi.carlo/Desktop/A1_outs/infercnv/gene_ordering_infercnv.txt',
                                         annotations_file = '/Users/leonardi.carlo/Desktop/A1_outs/infercnv/sampleannot_d30.txt',
                                         ref_group_names = header)

infercnvSample_d30 = infercnv::run(infercnvSample_d30,
                                   cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                   out_dir="/Users/leonardi.carlo/Desktop/A1_outs/infercnv/InferCNV_d30",
                                   cluster_by_groups=TRUE,
                                   denoise=TRUE,
                                   HMM=TRUE)

library(phylogram)
library(dendextend)
dendod10<-read.dendrogram('/Users/leonardi.carlo/Desktop/A1_outs/infercnv/InferCNV_d10/infercnv.preliminary.observations_dendrogram.txt')
dendod20<-read.dendrogram('/Users/leonardi.carlo/Desktop/A1_outs/infercnv/InferCNV_d20/infercnv.preliminary.observations_dendrogram.txt')
dendod30<-read.dendrogram('/Users/leonardi.carlo/Desktop/A1_outs/infercnv/InferCNV_d30/infercnv.preliminary.observations_dendrogram.txt')

cutd30<-cutree(dendod30,h = 17)

DimPlot(Sample_expr_FastMNN,cells.highlight = names(cutd30[cutd30==3]))
names(cutd30[cutd30==3])

cutd10<-cutree(dendod10, h = 6)

DimPlot(Sample_expr_FastMNN,cells.highlight = names(cutd10[cutd10==2]))


cutd20<-cutree(dendod20, h=8)

DimPlot(Sample_expr_FastMNN,cells.highlight = c(names(cutd20[cutd20==4]),
                                                names(cutd10[cutd10==2]),
                                                names(cutd30[cutd30==3]))
)


#analysis per sample

samplebyorig$Tumor_d30<-NormalizeData(samplebyorig$Tumor_d30)
samplebyorig$Tumor_d30<-RunPCA(samplebyorig$Tumor_d30, assay = "RNA", verbose = FALSE)
samplebyorig$Tumor_d30 <- FindNeighbors(samplebyorig$Tumor_d30, reduction = "pca", dims = 1:30)
samplebyorig$Tumor_d30 <- FindClusters(samplebyorig$Tumor_d30, verbose = FALSE,resolution = c(0.3,0.4,0.5,0.6,0.7,0.8))
samplebyorig$Tumor_d30 <- RunUMAP(samplebyorig$Tumor_d30, reduction = "pca", dims = 1:30)
DimPlot(samplebyorig$Tumor_d30,group.by = 'Annotation')

samplebyorig$Healthy<-NormalizeData(samplebyorig$Healthy)
samplebyorig$Healthy<-RunPCA(samplebyorig$Healthy, assay = "RNA", verbose = FALSE)
samplebyorig$Healthy <- FindNeighbors(samplebyorig$Healthy, reduction = "pca", dims = 1:30)
samplebyorig$Healthy <- FindClusters(samplebyorig$Healthy, verbose = FALSE,resolution = c(0.3,0.4,0.5,0.6,0.7,0.8))
samplebyorig$Healthy <- RunUMAP(samplebyorig$Healthy, reduction = "pca", dims = 1:30)
DimPlot(samplebyorig$Healthy,group.by = 'Annotation')

annotTAM<-read.csv('/Users/leonardi.carlo/Downloads/annotated_clusters_subset.csv')

ggarrange(DimPlot(samplebyorig$Tumor_d30,group.by = 'Annotation',label=T)+NoLegend(),
          DimPlot(samplebyorig$Tumor_d30,group.by = 'RNA_snn_res.0.6',label=T)+NoLegend(),
          DimPlot(samplebyorig$Tumor_d30,group.by = 'RNA_snn_res.0.7',label=T)+NoLegend()
)

annot1<-samplebyorig$Tumor_d30@meta.data$Annotation
names(annot1)<-rownames(samplebyorig$Tumor_d30@meta.data)
Idents(samplebyorig$Tumor_d30)='RNA_snn_res.0.6'

MatNeu<-WhichCells(samplebyorig$Tumor_d30,idents = 5)
names<-rep("Mature Neutrophiles",length(MatNeu))
names(names)=MatNeu
MatNeu<-names
rm(names)

annot1[names(MatNeu)]="Mature Neutrophiles"

Idents(samplebyorig$Tumor_d30)='RNA_snn_res.0.6'
cancer<-WhichCells(samplebyorig$Tumor_d30,idents = c(1,0,3))
names<-rep("cancer",length(cancer))
names(names)=cancer
cancer<-names
rm(names)

annot1[names(cancer)]="Cancer"

samplebyorig$Tumor_d30<-AddMetaData(samplebyorig$Tumor_d30,metadata = annot1,col.name = 'Annotation')
tumord30<-samplebyorig$Tumor_d30

tumord30 <- FindClusters(tumord30, verbose = FALSE,resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5))


#trial casper
library(dplyr)

a<-rownames(td10)
b<-rownames(td20)
c<-rownames(td30)

gtf <- rtracklayer::import('/Users/leonardi.carlo/Desktop/genes.gtf')
gtf_df1=as.data.frame(gtf)
gtf_df1<-gtf_df1%>%
  filter(type == 'gene')%>%
  select(gene_id,gene_name)
gtf_df1$gene_name<-make.unique(gtf_df1$gene_name)
gtf_df1<-gtf_df1[gtf_df1$gene_name%in%a,]

rownames(gtf_df1)=gtf_df1$gene_name
gtf_df1<-gtf_df1[a,]
rownames(gtf_df1)=NULL

rownames(td10)=gtf_df1$gene_id

gtf_df1=as.data.frame(gtf)
gtf_df1<-gtf_df1%>%
  filter_(gtf_df1$type=='gene')%>%
  select(gene_id,gene_name)
gtf_df1$gene_name<-make.unique(gtf_df1$gene_name)
gtf_df1<-gtf_df1[gtf_df1$gene_name%in%b,]

rownames(gtf_df1)=gtf_df1$gene_name
gtf_df1<-gtf_df1[b,]
rownames(gtf_df1)=NULL

rownames(td20)=gtf_df1$gene_id

gtf_df1=as.data.frame(gtf)
gtf_df1<-gtf_df1%>%
  filter(gtf_df1$type=='gene')%>%
  select(gene_id,gene_name)
gtf_df1$gene_name<-make.unique(gtf_df1$gene_name)
gtf_df1<-gtf_df1[gtf_df1$gene_name%in%c,]

rownames(gtf_df1)=gtf_df1$gene_name
gtf_df1<-gtf_df1[c,]
rownames(gtf_df1)=NULL

rownames(td30)=gtf_df1$gene_id

cytoband <- read.delim("/Users/leonardi.carlo/Downloads/cytoBand.txt", header=F)
cytoband <- data.frame(V1=gsub("chr", "", cytoband[,1]), V2=cytoband[,2], V3=cytoband[,3], V4=substring(cytoband$V4, 1, 1), stringsAsFactors=F)
start <- do.call(rbind, lapply(split(cytoband$V2, paste0(cytoband$V1, cytoband$V4)), min))
end <- do.call(rbind, lapply(split(cytoband$V3, paste0(cytoband$V1, cytoband$V4)), max))
cytoband <- data.frame(V1=gsub("p", "", gsub("q", "", rownames(start))), V2=start, V3=end, V4=rownames(start), stringsAsFactors=F)
cytoband <- cytoband [as.vector(unlist(sapply(c(1:22, "X"), function(x) which(cytoband$V1 %in% x)))), ]
cytoband$V4[grep("q", cytoband$V4)] <- "q"
cytoband$V4[grep("p", cytoband$V4)] <- "p"
rownames(cytoband) <- NULL

library(CaSpER)
generateAnnotation1 <- function (id_type = "ensembl_gene_id", genes)
{
  id_type = "ensembl_gene_id"
  mart <- useDataset("mmusculus_gene_ensembl", useEnsembl(biomart = "ensembl"))

  G_list <- getBM(filters = "ensembl_gene_id", attributes = c(id_type,
                                                              "mgi_symbol", "ensembl_gene_id", "chromosome_name", "start_position",
                                                              "end_position", "band"), values = genes, mart = mart)
  common <- intersect(genes, G_list[, id_type])
  ord <- match(common, G_list[, id_type])
  annotation <- G_list[ord, ]
  annotation <- annotation[order(annotation$start_position),
  ]
  annotation$cytoband <- paste0(annotation$chromosome_name,
                                substr((annotation$band), 0, 1))
  chr.list <- as.character(1:19)
  idx <- unlist(unique(as.vector(sapply(chr.list, function(x) as.vector(unlist(which(as.character(annotation$chromosome_name) ==
                                                                                       x)))))))
  annotation <- as.data.frame(annotation[idx, ])
  colnames(annotation)[c(1:6, 8)] <- c("Gene", "GeneSymbol",
                                       "EntrezID", "Chr", "start", "end", "cytoband")
  annotation$isCentromer <- rep("no", nrow(annotation))


  annotation$Position <- (as.numeric(annotation$start) + as.numeric(annotation$end))/2
  annotation$new_positions <- as.vector(unlist(lapply(lapply(split(annotation$cytoband,
                                                                   annotation$cytoband), length)[unique(annotation$cytoband)],
                                                      function(x) 1:x)))
  return(annotation)

}
annotation <- generateAnnotation1(id_type="ensembl_gene_id", genes=rownames(td10))
annotation1 <- generateAnnotation1(id_type="ensembl_gene_id", genes=rownames(td20))
annotation2 <- generateAnnotation1(id_type="ensembl_gene_id", genes=rownames(td30))

loh <- readBAFExtractOutput(path="/Users/leonardi.carlo/Desktop/bamfiles_tumor_TC/TD10/BAFs/",sequencing.type = 'single-cell')
names(loh) <- gsub(".snp", "", names(loh))
loh1 <- readBAFExtractOutput(path="/Users/leonardi.carlo/Desktop/bamfiles_tumor_TC/TD20/BAFs/",sequencing.type = 'single-cell')
names(loh1) <- gsub(".snp", "", names(loh1))
loh2 <- readBAFExtractOutput(path="/Users/leonardi.carlo/Desktop/bamfiles_tumor_TC/TD30/BAFs/",sequencing.type = 'single-cell')
names(loh2) <- gsub(".snp", "", names(loh2))


loh.name.mapping<-data.frame(loh.name=rep('TumorD10',dim(td10)[2]),sample.name=colnames(td10))
loh.name.mapping1<-data.frame(loh.name=rep('TumorD20',dim(td20)[2]),sample.name=colnames(td20))
loh.name.mapping2<-data.frame(loh.name=rep('TumorD30',dim(td30)[2]),sample.name=colnames(td30))


control.sample.ids<-WhichCells(samplebyorig$Tumor_d10,idents = c('Macrophages','Mature Neutrophiles',"Classical Monocytes","T cells","mDCs","B cells","NK","Endothelial cells",
                                                              "Acinar cells","cDCs","Immature B cells","CAF","pDCs","Immature Neutrophiles"))
control.sample.ids1<-WhichCells(samplebyorig$Tumor_d20,idents = c('Macrophages','Mature Neutrophiles',"Classical Monocytes","T cells","mDCs","B cells","NK","Endothelial cells",
                                                                 "Acinar cells","cDCs","Immature B cells","CAF","pDCs","Immature Neutrophiles"))
control.sample.ids2<-WhichCells(samplebyorig$Tumor_d30,idents = c('Macrophages','Mature Neutrophiles',"Classical Monocytes","T cells","mDCs","B cells","NK","Endothelial cells",
                                                                 "Acinar cells","cDCs","Immature B cells","CAF","pDCs","Immature Neutrophiles"))

object <- CreateCasperObject(raw.data=as.matrix(td10),loh.name.mapping=loh.name.mapping, sequencing.type="single-cell",
                             cnv.scale=3, loh.scale=3, matrix.type="raw", expr.cutoff=0.1,
                             annotation=annotation, method="iterative", loh=loh, filter="median",
                             control.sample.ids=control.sample.ids, cytoband=cytoband,genomeVersion='mm10')

object1 <- CreateCasperObject(raw.data=as.matrix(td20),loh.name.mapping=loh.name.mapping1, sequencing.type="single-cell",
                             cnv.scale=3, loh.scale=3, matrix.type="raw", expr.cutoff=0.1,
                             annotation=annotation1, method="iterative", loh=loh1, filter="median",
                             control.sample.ids=control.sample.ids1, cytoband=cytoband,genomeVersion='mm10')

object2 <- CreateCasperObject(raw.data=as.matrix(td30),loh.name.mapping=loh.name.mapping2, sequencing.type="single-cell",
                             cnv.scale=3, loh.scale=3, matrix.type="raw", expr.cutoff=0.1,
                             annotation=annotation2, method="iterative", loh=loh2, filter="median",
                             control.sample.ids=control.sample.ids2, cytoband=cytoband,genomeVersion='mm10')

final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")

final.objects1 <- runCaSpER(object1, removeCentromere=T, cytoband=cytoband, method="iterative")

final.objects2 <- runCaSpER(object2, removeCentromere=T, cytoband=cytoband, method="iterative")

#after that has been run,: large-scale CNV summarisation, segment-based CNV and then gene-based CNV summarisation (per cell visualisation)


extractLargeScaleEventsMouse <- function(final.objects, thr = 0.5) {

  mergeScales <- mergeScalesAndGenerateFinalEventSummary(final.objects)
  mergeScalesAmp <- mergeScales$mergeScalesAmp
  mergeScalesDel <- mergeScales$mergeScalesDel

  finalChrMat <- matrix(0, ncol =19, nrow = length(rownames(mergeScales$mergeScalesAmp)))
  colnames(finalChrMat) <- 1:19
  rownames(finalChrMat) <- rownames(mergeScales$mergeScalesAmp)

  finalChrMat[(mergeScalesAmp/length(final.objects)) >= thr] <- 1
  finalChrMat[(mergeScalesDel/length(final.objects)) >= thr] <- (-1)

  return(finalChrMat)
}

finalChrMat_thr5 <- extractLargeScaleEventsMouse (final.objects1, thr=0.5)
finalChrMat_thr75 <- extractLargeScaleEventsMouse (final.objects1, thr=0.75)

gamma <- 6
all.segments <- do.call(rbind, lapply(final.objects1, function(x) x@segments))
segment.summary <- extractSegmentSummary (final.objects1)
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loh <- segment.summary$all.summary.loh
loss.final <- loss[loss$count>gamma, ]
gain.final <- gain[gain$count>gamma, ]
loh.final <- loh[loh$count>gamma, ]

all.summary<- rbind(loss.final, gain.final)
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))),
                IRanges(all.summary$Start, all.summary$End))
ann.gr <- makeGRangesFromDataFrame(final.objects1[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "GeneSymbol")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects1[[1]]@annotation.filt[,2])
all.samples <- unique(as.character(final.objects1[[1]]@segments$ID))
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

plotHeatmapMM <- function (object, fileName, cnv.scale = 3, cluster_cols = F,
                           cluster_rows = T, show_rownames = T, only_soi = F)

{
  assignInNamespace(x = "draw_matrix", value = draw_matrix2,
                    ns = asNamespace("pheatmap"))
  assignInNamespace(x = "draw_colnames", value = "draw_colnames_45",
                    ns = asNamespace("pheatmap"))
  breaks <- seq(-2, 2, by = 0.2)
  color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks))
  idx <- cumsum(table(object@annotation.filt$Chr)[as.character(1:19)])
  xlabel <- rep("", length(rownames(object@data)))
  half <- round(table(object@annotation.filt$Chr)[as.character(1:19)]/2)[-1]
  xpos <- c(half[1], (idx[-19] + half))
  xlabel[xpos] <- 1:19
  data <- log2(object@control.normalized.noiseRemoved[[cnv.scale]])
  p<- pheatmap(t(data), cluster_cols = F, cluster_rows = T, gaps_col = idx,
               color = color, breaks = breaks, labels_col = xlabel,
               show_rownames = T, filename = fileName)
  return(p)
}

plotHeatmap10x <- function (object, fileName, cnv.scale = 3, cluster_cols = F,
                            cluster_rows = T, show_rownames = T, only_soi = T)
{
  assignInNamespace(x = "draw_matrix", value = draw_matrix2,
                    ns = asNamespace("pheatmap"))
  assignInNamespace(x = "draw_colnames", value = "draw_colnames_45",
                    ns = asNamespace("pheatmap"))
  data <- object@control.normalized.noiseRemoved[[cnv.scale]]
  x.center <- mean(data)
  quantiles = quantile(data[data != x.center], c(0.01, 0.99))
  delta = max(abs(c(x.center - quantiles[1], quantiles[2] -
                      x.center)))
  low_threshold = x.center - delta
  high_threshold = x.center + delta
  x.range = c(low_threshold, high_threshold)
  data[data < low_threshold] <- low_threshold
  data[data > high_threshold] <- high_threshold
  breaks <- seq(x.range[1], x.range[2], length = 16)
  color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks))
  idx <- cumsum(table(object@annotation.filt$Chr)[as.character(1:19)])
  xlabel <- rep("", length(rownames(object@data)))
  half <- round(table(object@annotation.filt$Chr)[as.character(1:19)]/2)[-1]
  xpos <- c(half[1], (idx[-19] + half))
  xlabel[xpos] <- 1:19
  if (only_soi)
    data <- data[, !(colnames(data) %in% object@control.sample.ids)]
  pheatmap(t(data), cluster_cols = F, cluster_rows = T, gaps_col = idx,
           color = color, breaks = breaks, labels_col = xlabel,
           show_rownames = T, filename = "heatmap.png")
}

obj <- final.objects1[[9]]
plotHeatmap(object=obj, fileName="heatmap.png",cnv.scale= 3, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = T)




