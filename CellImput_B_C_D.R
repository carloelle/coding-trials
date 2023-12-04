#annotation based on DestVi proportions

#retrieve DestVi proportions

proportionsA1<-read.csv('/Users/leonardi.carlo/Desktop/A1_outs/ProportionsDestVi/A1proportions.csv')
proportionsB1<-read.csv('/Users/leonardi.carlo/Desktop/A1_outs/ProportionsDestVi/B1proportions.csv')
proportionsC1<-read.csv('/Users/leonardi.carlo/Desktop/A1_outs/ProportionsDestVi/C1proportions.csv')
proportionsD1<-read.csv('/Users/leonardi.carlo/Desktop/A1_outs/ProportionsDestVi/D1proportions.csv')

rownames(proportionsA1)=proportionsA1$X
proportionsA1<-proportionsA1[,2:11]
colnames(proportionsA1)=gsub('$','_DestVI',colnames(proportionsA1))

rownames(proportionsB1)=proportionsB1$X
proportionsB1<-proportionsB1[,2:11]
colnames(proportionsB1)=gsub('$','_DestVI',colnames(proportionsB1))

rownames(proportionsC1)=proportionsC1$X
proportionsC1<-proportionsC1[,2:11]
colnames(proportionsC1)=gsub('$','_DestVI',colnames(proportionsC1))

rownames(proportionsD1)=proportionsD1$X
proportionsD1<-proportionsD1[,2:11]
colnames(proportionsD1)=gsub('$','_DestVI',colnames(proportionsD1))

SpatialA1FilteredFil<-AddMetaData(SpatialA1FilteredFil,metadata = proportionsA1,col.name = colnames(proportionsA1))
SpatialB1FilteredFil<-AddMetaData(SpatialB1FilteredFil,metadata = proportionsB1,col.name = colnames(proportionsB1))
SpatialC1FilteredFil<-AddMetaData(SpatialC1FilteredFil,metadata = proportionsC1,col.name = colnames(proportionsC1))
SpatialD1FilteredFil<-AddMetaData(SpatialD1FilteredFil,metadata = proportionsD1,col.name = colnames(proportionsD1))

tumor<-SpatialA1FilteredFil@meta.data$Tumor_DestVI+SpatialA1FilteredFil@meta.data$Acinar_DestVI+SpatialA1FilteredFil@meta.data$Ductal_DestVI
stromal<-SpatialA1FilteredFil@meta.data$CAF_DestVI+SpatialA1FilteredFil@meta.data$Endothelial_DestVI

onehot<-data.frame(Cell=rownames(SpatialA1FilteredFil@meta.data),Tumor=tumor,Stromal=stromal)
rownames(onehot)=onehot$Cell

residual<-onehot[intersect(onehot$Cell[onehot$Stromal>0.2],onehot$Cell[onehot$Tumor>0.5]),]
library(ggpubr)

ggarrange(
  SpatialDimPlot(SpatialA1FilteredFil,cells.highlight = onehot$Cell[onehot$Tumor>0.5]),
  SpatialDimPlot(SpatialA1FilteredFil,cells.highlight = onehot$Cell[onehot$Stromal>0.2]),
  SpatialDimPlot(SpatialA1FilteredFil,cells.highlight = residual$Cell[residual$Tumor>0.55]),
  SpatialDimPlot(SpatialA1FilteredFil,cells.highlight = intersect(stromal1,tumor1))
)

residual<-onehot[intersect(onehot$Cell[onehot$Stromal>0.2],onehot$Cell[onehot$Tumor>0.5]),]

stromal1<-setdiff(onehot$Cell[onehot$Stromal>0.2],residual$Cell[residual$Tumor>0.55])
tumor1<-setdiff(onehot$Cell[onehot$Tumor>0.5],residual$Cell[residual$Stromal>0.3])

tumor1<-setdiff(tumor1,intersect(stromal1,tumor1))


ggarrange(
  SpatialDimPlot(SpatialA1FilteredFil,cells.highlight = tumor1),
  SpatialDimPlot(SpatialA1FilteredFil,cells.highlight = stromal1)
)

included<-c(stromal1,tumor1)
excluded<-setdiff(rownames(SpatialA1FilteredFil@meta.data),included)

ggarrange(
  SpatialDimPlot(SpatialA1FilteredFil,cells.highlight = tumor1),
  SpatialDimPlot(SpatialA1FilteredFil,cells.highlight = stromal1),
  SpatialDimPlot(SpatialA1FilteredFil,cells.highlight = stilltumor),
  SpatialDimPlot(SpatialA1FilteredFil,cells.highlight = stillstromal)
)

#look at histograms

stilltumor<-excluded[onehot[excluded,'Tumor']>0.3]
stillstromal<-setdiff(excluded,stilltumor)

tumor1<-c(tumor1,stilltumor)
stromal1<-c(stromal1,stillstromal)

#final
ggarrange(SpatialDimPlot(SpatialA1FilteredFil,cells.highlight = tumor1)+ggtitle('Tumor area'),
          SpatialDimPlot(SpatialA1FilteredFil,cells.highlight = stromal1)+ggtitle('Stromal area'))

tum<-rep('Tumor',length(tumor1))
names(tum)=tumor1
str<-rep('Stromal',length(stromal1))
names(str)=stromal1

A1Annot<-c(tum,str)
A1Annot<-A1Annot[rownames(SpatialA1FilteredFil@meta.data)]
SpatialA1FilteredFil<-AddMetaData(SpatialA1FilteredFil,col.name = 'TumorOrStromal_DestVI',metadata = A1Annot)

tumor<-SpatialB1FilteredFil@meta.data$Tumor_DestVI+SpatialB1FilteredFil@meta.data$Acinar_DestVI+SpatialB1FilteredFil@meta.data$Ductal_DestVI
stromal<-SpatialB1FilteredFil@meta.data$CAF_DestVI+SpatialB1FilteredFil@meta.data$Endothelial_DestVI
onehot<-data.frame(Cell=rownames(SpatialB1FilteredFil@meta.data),Tumor=tumor,Stromal=stromal)
rownames(onehot)=onehot$Cell

ggarrange(
  SpatialDimPlot(SpatialB1FilteredFil,cells.highlight = onehot$Cell[onehot$Tumor>0.4]),
  SpatialDimPlot(SpatialB1FilteredFil,cells.highlight = onehot$Cell[onehot$Tumor>0.5]),
  SpatialDimPlot(SpatialB1FilteredFil,cells.highlight = onehot$Cell[onehot$Stromal>0.2]),
  SpatialDimPlot(SpatialB1FilteredFil,cells.highlight = onehot$Cell[onehot$Stromal>0.3])
)

#chose 0.5 0.3


residual<-onehot[intersect(onehot$Cell[onehot$Stromal>0.3],onehot$Cell[onehot$Tumor>0.5]),]

stromal1<-setdiff(onehot$Cell[onehot$Stromal>0.3],residual$Cell[residual$Tumor>0.6])
tumor1<-setdiff(onehot$Cell[onehot$Tumor>0.5],residual$Cell[residual$Stromal>0.36])

ggarrange(
  SpatialDimPlot(SpatialB1FilteredFil,cells.highlight = stromal1),
  SpatialDimPlot(SpatialB1FilteredFil,cells.highlight = tumor1),
  SpatialDimPlot(SpatialB1FilteredFil,cells.highlight = intersect(stromal1,tumor1))
)

stromal1<-setdiff(stromal1,intersect(stromal1,tumor1))

ggarrange(
  SpatialDimPlot(SpatialB1FilteredFil,cells.highlight = tumor1),
  SpatialDimPlot(SpatialB1FilteredFil,cells.highlight = intersect(tumor1,stromal1))
)

included<-c(stromal1,tumor1)
excluded<-setdiff(rownames(SpatialB1FilteredFil@meta.data),included)

stilltumor<-excluded[onehot[excluded,'Tumor']>0.5]
stillstromal<-excluded[onehot[excluded,'Stromal']>0.3]

#stilltumor and stillstromal are the same cells, check if they're more stromal or tumoral and sort them
#tumor is 0.6 in all, go with tumoral

tumor1<-c(tumor1,stilltumor)
stromal1<-c(stromal1)

stromal1<-setdiff(stromal1,intersect(tumor1,stromal1))

ggarrange(
  SpatialDimPlot(SpatialB1FilteredFil,cells.highlight = tumor1),
  SpatialDimPlot(SpatialB1FilteredFil,cells.highlight = stromal1)
)

excluded<-setdiff(rownames(SpatialB1FilteredFil@meta.data),c(tumor1,stromal1))

excluded_str<-excluded[onehot[excluded,'Stromal']>0.2]
excluded_tum<-excluded[onehot[excluded,'Tumor']>0.3]

excluded_str<-setdiff(excluded_str,intersect(excluded_str,excluded_tum))

stromal1<-c(stromal1,excluded_str)
tumor1<-c(tumor1,excluded_tum)

stromal1<-c(stromal1,setdiff(rownames(SpatialB1FilteredFil@meta.data),c(tumor1,stromal1)))

tum<-rep('Tumor',length(tumor1))
names(tum)=tumor1
str<-rep('Stromal',length(stromal1))
names(str)=stromal1


B1Annot<-c(tum,str)
B1Annot<-B1Annot[rownames(SpatialB1FilteredFil@meta.data)]
SpatialB1FilteredFil<-AddMetaData(SpatialB1FilteredFil,col.name = 'TumorOrStromal_DestVI',metadata = B1Annot)

---

tumor<-SpatialC1FilteredFil@meta.data$Tumor_DestVI+SpatialC1FilteredFil@meta.data$Acinar_DestVI+SpatialC1FilteredFil@meta.data$Ductal_DestVI
stromal<-SpatialC1FilteredFil@meta.data$CAF_DestVI+SpatialC1FilteredFil@meta.data$Endothelial_DestVI
onehot<-data.frame(Cell=rownames(SpatialC1FilteredFil@meta.data),Tumor=tumor,Stromal=stromal)
rownames(onehot)=onehot$Cell

ggarrange(
  SpatialDimPlot(SpatialC1FilteredFil,cells.highlight = onehot$Cell[onehot$Tumor>0.5]),
  SpatialDimPlot(SpatialC1FilteredFil,cells.highlight = onehot$Cell[onehot$Stromal>0.4])
)

residual<-onehot[intersect(onehot$Cell[onehot$Stromal>0.4],onehot$Cell[onehot$Tumor>0.5]),]

stromal1<-setdiff(onehot$Cell[onehot$Stromal>0.4],residual$Cell)
tumor1<-unique(c(onehot$Cell[onehot$Tumor>0.5],residual$Cell))

ggarrange(
  SpatialDimPlot(SpatialC1FilteredFil,cells.highlight = tumor1),
  SpatialDimPlot(SpatialC1FilteredFil,cells.highlight = stromal1)
)

included<-c(stromal1,tumor1)
excluded<-setdiff(rownames(SpatialC1FilteredFil@meta.data),included)

stilltumor<-excluded[onehot[excluded,'Tumor']>0.35]
stillstromal<-setdiff(excluded,stilltumor)

tumor1<-c(tumor1,stilltumor)
stromal1<-c(stromal1,stillstromal)

tum<-rep('Tumor',length(tumor1))
names(tum)=tumor1
str<-rep('Stromal',length(stromal1))
names(str)=stromal1

C1Annot<-c(tum,str)
C1Annot<-C1Annot[rownames(SpatialC1FilteredFil@meta.data)]
SpatialC1FilteredFil<-AddMetaData(SpatialC1FilteredFil,col.name = 'TumorOrStromal_DestVI',metadata = C1Annot)


tumor<-SpatialD1FilteredFil@meta.data$Tumor_DestVI+SpatialD1FilteredFil@meta.data$Acinar_DestVI+SpatialD1FilteredFil@meta.data$Ductal_DestVI
stromal<-SpatialD1FilteredFil@meta.data$CAF_DestVI+SpatialD1FilteredFil@meta.data$Endothelial_DestVI
onehot<-data.frame(Cell=rownames(SpatialD1FilteredFil@meta.data),Tumor=tumor,Stromal=stromal)
rownames(onehot)=onehot$Cell


ggarrange(
  SpatialDimPlot(SpatialD1FilteredFil,cells.highlight = onehot$Cell[onehot$Tumor>0.6]),
  SpatialDimPlot(SpatialD1FilteredFil,cells.highlight = onehot$Cell[onehot$Stromal>0.4])
)

residual<-onehot[setdiff(rownames(SpatialD1FilteredFil@meta.data),c(onehot$Cell[onehot$Stromal>0.4],onehot$Cell[onehot$Tumor>0.6])),]

stilltumor<-residual[residual$Tumor>0.4,]

tumor1<-c(stilltumor$Cell,onehot$Cell[onehot$Tumor>0.6])
stromal1<-c(onehot$Cell[onehot$Stromal>0.4])

ggarrange(
  SpatialDimPlot(SpatialD1FilteredFil,cells.highlight = stromal1),
  SpatialDimPlot(SpatialD1FilteredFil,cells.highlight = tumor1)
)

included<-c(tumor1,stromal1)
excluded<-setdiff(rownames(SpatialD1FilteredFil@meta.data),included)

stillstromal<-excluded[onehot[excluded,'Stromal']>0.25]
stilltumor<-setdiff(excluded,stillstromal)

tumor1<-c(tumor1,stilltumor)
stromal1<-c(stromal1,stillstromal)

tum<-rep('Tumor',length(tumor1))
names(tum)=tumor1
str<-rep('Stromal',length(stromal1))
names(str)=stromal1

D1Annot<-c(tum,str)
D1Annot<-D1Annot[rownames(SpatialD1FilteredFil@meta.data)]
SpatialD1FilteredFil<-AddMetaData(SpatialD1FilteredFil,col.name = 'TumorOrStromal_DestVI',metadata = D1Annot)

ggarrange(
  SpatialDimPlot(SpatialA1FilteredFil,group.by = 'TumorOrStromal_DestVI'),
  SpatialDimPlot(SpatialB1FilteredFil,group.by = 'TumorOrStromal_DestVI'),
  SpatialDimPlot(SpatialC1FilteredFil,group.by = 'TumorOrStromal_DestVI'),
  SpatialDimPlot(SpatialD1FilteredFil,group.by = 'TumorOrStromal_DestVI')
)


pdf("A1_stromalOrtumor_NoCountour.pdf")
SpatialDimPlot(SpatialA1FilteredFil,group.by = 'TumorOrStromal_DestVI',pt.size.factor = 1.8,image.alpha = 0,stroke=0)+guides(fill=guide_legend('Annotation'))
dev.off()


pdf("B1_stromalOrtumor_NoCountour.pdf")
SpatialDimPlot(SpatialB1FilteredFil,group.by = 'TumorOrStromal_DestVI',pt.size.factor = 1.8,image.alpha = 0,stroke=0)+guides(fill=guide_legend('Annotation'))
dev.off()

pdf("C1_stromalOrtumor_NoCountour.pdf")
SpatialDimPlot(SpatialC1FilteredFil,group.by = 'TumorOrStromal_DestVI',pt.size.factor = 1.8,image.alpha = 0,stroke=0)+guides(fill=guide_legend('Annotation'))
dev.off()

pdf("D1_stromalOrtumo_NoCountour.pdf")
SpatialDimPlot(SpatialD1FilteredFil,group.by = 'TumorOrStromal_DestVI',pt.size.factor = 1.8,image.alpha = 0,stroke=0)+guides(fill=guide_legend('Annotation'))
dev.off()


pdf('SpatialLocations_MonoMacro_A1_NoContour.pdf')
SpatialFeaturePlot(SpatialA1FilteredFil,features = 'Mono.Macrophages_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#FFCC99',high = 'red',midpoint=0.2,limits=c(0.1,NA),na.value = "#E0E0E0",name='% of Mono/Macrophages')
dev.off()

pdf('SpatialLocations_MonoMacro_B1_NoContour.pdf')
SpatialFeaturePlot(SpatialB1FilteredFil,features = 'Mono.Macrophages_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#FFCC99',high = 'red',midpoint=0.2,limits=c(0.1,NA),na.value = "#E0E0E0",name='% of Mono/Macrophages')
dev.off()


pdf('SpatialLocations_CAF_A1_NoContour_red.pdf')
SpatialFeaturePlot(SpatialA1FilteredFil,features = 'CAF_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#FFCC99',high = 'red',midpoint=0.6,limits=c(0.4,NA),na.value = "#E0E0E0",name='% of CAF')
dev.off()

pdf('SpatialLocations_CAF_B1_NoContour_red.pdf')
SpatialFeaturePlot(SpatialB1FilteredFil,features = 'CAF_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#FFCC99',high = 'red',midpoint=0.2,limits=c(0.25,NA),na.value = "#E0E0E0",name='% of CAF')
dev.off()


pdf('SpatialLocations_Tumor_A1_AmetistaMelanzana_Lim02.pdf')
SpatialFeaturePlot(SpatialA1FilteredFil,features = 'Tumor_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#9966CC',high = '#990066',midpoint=0.6,limits=c(0.2,NA),na.value = "#E0E0E0",name='% of Tumor')
dev.off()
pdf('SpatialLocations_Tumor_A1_AmetistaIndacodark_Lim02.pdf')
SpatialFeaturePlot(SpatialA1FilteredFil,features = 'Tumor_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#9966CC',high = '#310062',midpoint=0.6,limits=c(0.2,NA),na.value = "#E0E0E0",name='% of Tumor')
dev.off()
pdf('SpatialLocations_Tumor_B1_AmetistaMelanzana_Lim02.pdf')
SpatialFeaturePlot(SpatialB1FilteredFil,features = 'Tumor_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#9966CC',high = '#990066',midpoint=0.6,limits=c(0.2,NA),na.value = "#E0E0E0",name='% of Tumor')
dev.off()
pdf('SpatialLocations_Tumor_B1_AmetistaIndacodark_Lim02.pdf')
SpatialFeaturePlot(SpatialB1FilteredFil,features = 'Tumor_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#9966CC',high = '#310062',midpoint=0.6,limits=c(0.2,NA),na.value = "#E0E0E0",name='% of Tumor')
dev.off()





pdf('SpatialLocations_CAF_A1_AmetistaMelanzana.pdf')
SpatialFeaturePlot(SpatialA1FilteredFil,features = 'Tumor_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#9966CC',high = '#990066',midpoint=0.8,limits=c(0.6,NA),na.value = "#E0E0E0",name='% of Tumor')
dev.off()
pdf('SpatialLocations_CAF_A1_AmetistaIndacodark.pdf')
SpatialFeaturePlot(SpatialA1FilteredFil,features = 'Tumor_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#9966CC',high = '#310062',midpoint=0.8,limits=c(0.6,NA),na.value = "#E0E0E0",name='% of Tumor')
dev.off()

pdf('SpatialLocations_CAF_B1_AmetistaMelanzana.pdf')
SpatialFeaturePlot(SpatialB1FilteredFil,features = 'Tumor_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#9966CC',high = '#990066',midpoint=0.7,limits=c(0.5,NA),na.value = "#E0E0E0",name='% of Tumor')
dev.off()
pdf('SpatialLocations_CAF_B1_AmetistaIndacodark.pdf')
SpatialFeaturePlot(SpatialB1FilteredFil,features = 'Tumor_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#9966CC',high = '#310062',midpoint=0.7,limits=c(0.5,NA),na.value = "#E0E0E0",name='% of Tumor')
dev.off()



pdf('SpatialLocations_CAF_A1_NoContour_OlivaMirto.pdf')
SpatialFeaturePlot(SpatialA1FilteredFil,features = 'CAF_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#808000',high = '#21421E',midpoint=0.6,limits=c(0.3,NA),na.value = "#E0E0E0",name='% of CAF')
dev.off()
pdf('SpatialLocations_CAF_A1_NoContour_OlivastroMirto.pdf')
SpatialFeaturePlot(SpatialA1FilteredFil,features = 'CAF_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#6B8E23',high = '#21421E',midpoint=0.6,limits=c(0.3,NA),na.value = "#E0E0E0",name='% of CAF')
dev.off()
pdf('SpatialLocations_CAF_B1_NoContour_OlivaMirto.pdf')
SpatialFeaturePlot(SpatialB1FilteredFil,features = 'CAF_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#808000',high = '#21421E',midpoint=0.55,limits=c(0.25,NA),na.value = "#E0E0E0",name='% of CAF')
dev.off()
pdf('SpatialLocations_CAF_B1_NoContour_OlivastroMirto.pdf')
SpatialFeaturePlot(SpatialB1FilteredFil,features = 'CAF_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#6B8E23',high = '#21421E',midpoint=0.55,limits=c(0.25,NA),na.value = "#E0E0E0",name='% of CAF')



pdf('SpatialLocations_MonoMacro_A1_NoContour_BlufaenzaZaffiro.pdf')
SpatialFeaturePlot(SpatialA1FilteredFil,features = 'Mono.Macrophages_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#004F92',high = '#1D1E33',midpoint=0.4,limits=c(0.1,NA),na.value = "#E0E0E0",name='% of Mono/Macrophages')
dev.off()
pdf('SpatialLocations_MonoMacro_A1_NoContour_BlufaenzaNotteperlato.pdf')
SpatialFeaturePlot(SpatialA1FilteredFil,features = 'Mono.Macrophages_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#004F92',high = '#102C54',midpoint=0.4,limits=c(0.1,NA),na.value = "#E0E0E0",name='% of Mono/Macrophages')
dev.off()
pdf('SpatialLocations_MonoMacro_B1_NoContour_BlufaenzaZaffiro.pdf')
SpatialFeaturePlot(SpatialB1FilteredFil,features = 'Mono.Macrophages_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#004F92',high = '#1D1E33',midpoint=0.35,limits=c(0.1,NA),na.value = "#E0E0E0",name='% of Mono/Macrophages')
dev.off()
pdf('SpatialLocations_MonoMacro_B1_NoContour_BlufaenzaNotteperlato.pdf')
SpatialFeaturePlot(SpatialB1FilteredFil,features = 'Mono.Macrophages_DestVI',image.alpha = 0,stroke = 0)+scale_fill_gradient2(low="#E0E0E0", mid = '#004F92',high = '#102C54',midpoint=0.35,limits=c(0.1,NA),na.value = "#E0E0E0",name='% of Mono/Macrophages')
dev.off()


pdf('VlnPlot_Endothelial_A1.pdf')
VlnPlot(SpatialA1FilteredFil,group.by = 'TumorOrStromal_DestVI',features = 'Endothelial_DestVI')+ggtitle('% of Endothelial cells')+xlab('Annotation')
dev.off()

pdf('VlnPlot_Endothelial_B1.pdf')
VlnPlot(SpatialB1FilteredFil,group.by = 'TumorOrStromal_DestVI',features = 'Endothelial_DestVI')+ggtitle('% of Endothelial cells')+xlab('Annotation')
dev.off()

pdf('VlnPlot_CAF_B1.pdf')
VlnPlot(SpatialB1FilteredFil,group.by = 'TumorOrStromal_DestVI',features = 'CAF_DestVI')+ggtitle('% of CAF')+xlab('Annotation')
dev.off()

pdf('VlnPlot_CAF_A1.pdf')
VlnPlot(SpatialA1FilteredFil,group.by = 'TumorOrStromal_DestVI',features = 'CAF_DestVI')+ggtitle('% of CAF')+xlab('Annotation')
dev.off()


pdf('VlnPlot_DCs_B1.pdf')
VlnPlot(SpatialB1FilteredFil,group.by = 'TumorOrStromal_DestVI',features = 'DCs_DestVI')+ggtitle('% of DCs')+xlab('Annotation')
dev.off()

pdf('VlnPlot_DCs_A1.pdf')
VlnPlot(SpatialA1FilteredFil,group.by = 'TumorOrStromal_DestVI',features = 'DCs_DestVI')+ggtitle('% of DCs')+xlab('Annotation')
dev.off()


pdf('VlnPlot_Granulo_B1.pdf')
VlnPlot(SpatialB1FilteredFil,group.by = 'TumorOrStromal_DestVI',features = 'Granulocytes_DestVI')+ggtitle('% of Granulocytes')+xlab('Annotation')
dev.off()

pdf('VlnPlot_Granulo_A1.pdf')
VlnPlot(SpatialA1FilteredFil,group.by = 'TumorOrStromal_DestVI',features = 'Granulocytes_DestVI')+ggtitle('% of Granulocytes')+xlab('Annotation')
dev.off()


SpatialA1FilteredFil<-AddMetaData(SpatialA1FilteredFil,metadata = c(SpatialA1FilteredFil@meta.data$B_cells_DestVI+SpatialA1FilteredFil@meta.data$T_NK_DestVI),
                                                        col.name = 'Lymphoid_DestVI')
SpatialB1FilteredFil<-AddMetaData(SpatialB1FilteredFil,metadata = c(SpatialB1FilteredFil@meta.data$B_cells_DestVI+SpatialB1FilteredFil@meta.data$T_NK_DestVI),
                                  col.name = 'Lymphoid_DestVI')

pdf('VlnPlot_Lymphoid_A1.pdf')
VlnPlot(SpatialA1FilteredFil,group.by = 'TumorOrStromal_DestVI',features = 'Lymphoid_DestVI')+ggtitle('% of Lymphoid')+xlab('Annotation')
dev.off()

pdf('VlnPlot_Lymphoid_B1.pdf')
VlnPlot(SpatialB1FilteredFil,group.by = 'TumorOrStromal_DestVI',features = 'Lymphoid_DestVI')+ggtitle('% of Lymphoid')+xlab('Annotation')
dev.off()


pdf('VlnPlot_MonoMacro_A1.pdf')
VlnPlot(SpatialA1FilteredFil,group.by = 'TumorOrStromal_DestVI',features = 'Mono.Macrophages_DestVI')+ggtitle('% of Macrophages')+xlab('Annotation')
dev.off()

pdf('VlnPlot_MonoMacro_B1.pdf')
VlnPlot(SpatialB1FilteredFil,group.by = 'TumorOrStromal_DestVI',features = 'Mono.Macrophages_DestVI')+ggtitle('% of Macrophages')+xlab('Annotation')
dev.off()



Idents(SpatialA1FilteredFil)='TumorOrStromal_DestVI'

SpatialA1FilteredFil_str<-subset(x=SpatialA1FilteredFil,idents='Stromal')

A1str_Endo<-median(SpatialA1FilteredFil_str@meta.data$Endothelial_DestVI)
A1str_CAF<-median(SpatialA1FilteredFil_str@meta.data$CAF_DestVI)
A1str_Lymphoid<-median(SpatialA1FilteredFil_str@meta.data$B_cells_DestVI+SpatialA1FilteredFil_str@meta.data$T_NK_DestVI)
A1str_DC<-median(SpatialA1FilteredFil_str@meta.data$DCs_DestVI)
A1str_Granulo<-median(SpatialA1FilteredFil_str@meta.data$Granulocytes_DestVI)

SpatialA1FilteredFil_tum<-subset(x=SpatialA1FilteredFil,idents='Tumor')

A1tum_Endo<-median(SpatialA1FilteredFil_tum@meta.data$Endothelial_DestVI)
A1tum_CAF<-median(SpatialA1FilteredFil_tum@meta.data$CAF_DestVI)
A1tum_Lymphoid<-median(SpatialA1FilteredFil_tum@meta.data$B_cells_DestVI+SpatialA1FilteredFil_tum@meta.data$T_NK_DestVI)
A1tum_DC<-median(SpatialA1FilteredFil_tum@meta.data$DCs_DestVI)
A1tum_Granulo<-median(SpatialA1FilteredFil_tum@meta.data$Granulocytes_DestVI)


cbind(c(A1str_Endo,
        A1str_CAF,
        A1str_Lymphoid,
        A1str_DC,
        A1str_Granulo),c(A1tum_Endo,
                         A1tum_CAF,
                         A1tum_Lymphoid,
                         A1tum_DC,
                         A1tum_Granulo))


Idents(SpatialB1FilteredFil)='TumorOrStromal_DestVI'

SpatialB1FilteredFil_str<-subset(x=SpatialB1FilteredFil,idents='Stromal')
SpatialB1FilteredFil_tum<-subset(x=SpatialB1FilteredFil,idents='Tumor')

tumord30_ALRA<-SeuratWrappers::RunALRA(tumord30)

tumord30_ALRA<-NormalizeData(tumord30_ALRA,assay='alra')
tumord30_ALRA<-ScaleData(tumord30_ALRA)
tumord30_ALRA<-FindVariableFeatures(tumord30_ALRA)
tumord30_ALRA<-RunPCA(tumord30_ALRA)
tumord30_ALRA<-FindNeighbors(tumord30_ALRA)
tumord30_ALRA<-RunUMAP(tumord30_ALRA,dims = 1:20)


test_meta<-data.frame(Cell=gsub('-1','.1',rownames(tumord30@meta.data)),cell_type=gsub(' ','.',tumord30@meta.data$Annotation))

celltypes<-gsub(' ','.',sort(unique(tumord30@meta.data$Annotation)))
location<-c('Tumor','Stroma','Stroma','Tumor','Stroma','Stroma','Tumor','Stroma','Stroma','Stroma','Stroma','Stroma','Stroma','Stroma','Stroma','Stroma')

test_microenvironment=data.frame(cell_type=celltypes,microenviroment=location)



write.table(tumord30@assays$RNA@data,file = 'test_counts.txt',quote = F,row.names = T,sep = "\t")
write.table(tumord30@assays$SCT@data,file = 'test_countsSCT.txt',quote = F,row.names = T,sep = "\t")

write.table(test_microenvironment,file = 'test_microenviroments.txt',quote = F,row.names = F,sep = "\t")
write.table(test_meta,file = 'test_meta.txt',quote = F,row.names = F,sep = "\t")

mousegenes<-rownames(tumord30)

convertMouseGeneList <- function(x){
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

convertHumanGeneList <- function(x){
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

humanort<-vector()
for(i in 1:length(mousegenes)){
  b<-convertMouseGeneList(mousegenes[i])
  humanort<-c(humanort,b)}

filteredmouse<-vector()
for(i in 1:length(humanort)){
  b<-convertHumanGeneList(humanort[i])
  filteredmouse<-c(filteredmouse,b)}


write.table(humanort,file = 'genesHuman.txt',quote = F,row.names = F,sep = "\t")

SpatialA1FilteredFil_TS<-SplitObject(SpatialA1FilteredFil,split.by = 'TumorOrStromal')
SpatialB1FilteredFil_TS<-SplitObject(SpatialB1FilteredFil,split.by = 'TumorOrStromal')
SpatialC1FilteredFil_TS<-SplitObject(SpatialC1FilteredFil,split.by = 'TumorOrStromal')
SpatialD1FilteredFil_TS<-SplitObject(SpatialD1FilteredFil,split.by = 'TumorOrStromal')

SpatialA1FilteredFil_TS$Tumor<-SCTransform(SpatialA1FilteredFil_TS$Tumor, assay = "Spatial", verbose = FALSE)
SpatialA1FilteredFil_TS$Tumor<-RunPCA(SpatialA1FilteredFil_TS$Tumor, assay = "SCT", verbose = FALSE)
SpatialA1FilteredFil_TS$Tumor<-FindNeighbors(SpatialA1FilteredFil_TS$Tumor, reduction = "pca", dims = 1:30)
SpatialA1FilteredFil_TS$Tumor<-FindClusters(SpatialA1FilteredFil_TS$Tumor, verbose = FALSE,resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))
SpatialA1FilteredFil_TS$Tumor<-RunUMAP(SpatialA1FilteredFil_TS$Tumor, reduction = "pca", dims = 1:30)

SpatialA1FilteredFil_TS$Stromal<-SCTransform(SpatialA1FilteredFil_TS$Stromal, assay = "Spatial", verbose = FALSE)
SpatialA1FilteredFil_TS$Stromal<-RunPCA(SpatialA1FilteredFil_TS$Stromal, assay = "SCT", verbose = FALSE)
SpatialA1FilteredFil_TS$Stromal<-FindNeighbors(SpatialA1FilteredFil_TS$Stromal, reduction = "pca", dims = 1:30)
SpatialA1FilteredFil_TS$Stromal<-FindClusters(SpatialA1FilteredFil_TS$Stromal, verbose = FALSE,resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))
SpatialA1FilteredFil_TS$Stromal<-RunUMAP(SpatialA1FilteredFil_TS$Stromal, reduction = "pca", dims = 1:30)

SpatialB1FilteredFil_TS$Tumor<-SCTransform(SpatialB1FilteredFil_TS$Tumor, assay = "Spatial", verbose = FALSE)
SpatialB1FilteredFil_TS$Tumor<-RunPCA(SpatialB1FilteredFil_TS$Tumor, assay = "SCT", verbose = FALSE)
SpatialB1FilteredFil_TS$Tumor<-FindNeighbors(SpatialB1FilteredFil_TS$Tumor, reduction = "pca", dims = 1:30)
SpatialB1FilteredFil_TS$Tumor<-FindClusters(SpatialB1FilteredFil_TS$Tumor, verbose = FALSE,resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))
SpatialB1FilteredFil_TS$Tumor<-RunUMAP(SpatialB1FilteredFil_TS$Tumor, reduction = "pca", dims = 1:30)

SpatialB1FilteredFil_TS$Stromal<-SCTransform(SpatialB1FilteredFil_TS$Stromal, assay = "Spatial", verbose = FALSE)
SpatialB1FilteredFil_TS$Stromal<-RunPCA(SpatialB1FilteredFil_TS$Stromal, assay = "SCT", verbose = FALSE)
SpatialB1FilteredFil_TS$Stromal<-FindNeighbors(SpatialB1FilteredFil_TS$Stromal, reduction = "pca", dims = 1:30)
SpatialB1FilteredFil_TS$Stromal<-FindClusters(SpatialB1FilteredFil_TS$Stromal, verbose = FALSE,resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))
SpatialB1FilteredFil_TS$Stromal<-RunUMAP(SpatialB1FilteredFil_TS$Stromal, reduction = "pca", dims = 1:30)

SpatialC1FilteredFil_TS$Tumor<-SCTransform(SpatialC1FilteredFil_TS$Tumor, assay = "Spatial", verbose = FALSE)
SpatialC1FilteredFil_TS$Tumor<-RunPCA(SpatialC1FilteredFil_TS$Tumor, assay = "SCT", verbose = FALSE)
SpatialC1FilteredFil_TS$Tumor<-FindNeighbors(SpatialC1FilteredFil_TS$Tumor, reduction = "pca", dims = 1:30)
SpatialC1FilteredFil_TS$Tumor<-FindClusters(SpatialC1FilteredFil_TS$Tumor, verbose = FALSE,resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))
SpatialC1FilteredFil_TS$Tumor<-RunUMAP(SpatialC1FilteredFil_TS$Tumor, reduction = "pca", dims = 1:30)

SpatialC1FilteredFil_TS$Stromal<-SCTransform(SpatialC1FilteredFil_TS$Stromal, assay = "Spatial", verbose = FALSE)
SpatialC1FilteredFil_TS$Stromal<-RunPCA(SpatialC1FilteredFil_TS$Stromal, assay = "SCT", verbose = FALSE)
SpatialC1FilteredFil_TS$Stromal<-FindNeighbors(SpatialC1FilteredFil_TS$Stromal, reduction = "pca", dims = 1:30)
SpatialC1FilteredFil_TS$Stromal<-FindClusters(SpatialC1FilteredFil_TS$Stromal, verbose = FALSE,resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))
SpatialC1FilteredFil_TS$Stromal<-RunUMAP(SpatialC1FilteredFil_TS$Stromal, reduction = "pca", dims = 1:30)

SpatialD1FilteredFil_TS$Tumor<-SCTransform(SpatialD1FilteredFil_TS$Tumor, assay = "Spatial", verbose = FALSE)
SpatialD1FilteredFil_TS$Tumor<-RunPCA(SpatialD1FilteredFil_TS$Tumor, assay = "SCT", verbose = FALSE)
SpatialD1FilteredFil_TS$Tumor<-FindNeighbors(SpatialD1FilteredFil_TS$Tumor, reduction = "pca", dims = 1:30)
SpatialD1FilteredFil_TS$Tumor<-FindClusters(SpatialD1FilteredFil_TS$Tumor, verbose = FALSE,resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))
SpatialD1FilteredFil_TS$Tumor<-RunUMAP(SpatialD1FilteredFil_TS$Tumor, reduction = "pca", dims = 1:30)

SpatialD1FilteredFil_TS$Stromal<-SCTransform(SpatialD1FilteredFil_TS$Stromal, assay = "Spatial", verbose = FALSE)
SpatialD1FilteredFil_TS$Stromal<-RunPCA(SpatialD1FilteredFil_TS$Stromal, assay = "SCT", verbose = FALSE)
SpatialD1FilteredFil_TS$Stromal<-FindNeighbors(SpatialD1FilteredFil_TS$Stromal, reduction = "pca", dims = 1:30)
SpatialD1FilteredFil_TS$Stromal<-FindClusters(SpatialD1FilteredFil_TS$Stromal, verbose = FALSE,resolution = c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,2))
SpatialD1FilteredFil_TS$Stromal<-RunUMAP(SpatialD1FilteredFil_TS$Stromal, reduction = "pca", dims = 1:30)

SpatialA1FilteredFil_stromalMG<-SpatialA1FilteredFil_TS$Stromal[,SpatialA1FilteredFil_TS$Stromal@meta.data$Macrophages_C2L>0.1]

SpatialA1FilteredFil_stromalMG<-SpatialA1FilteredFil_stromalMG[,SpatialA1FilteredFil_stromalMG@meta.data$TAMIl1b_hyper>1]

a<-SpatialA1FilteredFil_stromalMG@meta.data$TAMIl1b_hyper


Sample.expr_data.info<-SpatialA1FilteredFil_stromalMG@meta.data

dend<-as.dendrogram(hclust(dist(Sample.expr_data.info[,c(39:54)])))
order<-labels(dend)
a<-Sample.expr_data.info[,c(39:54)]
a<-a[order,]
a<-a %>%
  mutate(ind = factor(row_number()))%>%
  gather(variable,value,-ind)
a$variable<-factor(a$variable,levels = c('Macrophages_C2L',
                                         'CAF_C2L',
                                         "Endothelial.cells_C2L",
                                         "Cancer_C2L",
                                         "Ductal_C2L",
                                         "Acinar.cells_C2L",
                                         "T.cells_C2L",
                                         "NK_C2L",
                                         "B.cells_C2L",
                                         "Immature.B.cells_C2L",
                                         "pDCs_C2L",
                                         "Mature.Neutrophiles_C2L",
                                         "Immature.Neutrophiles_C2L",
                                         "mDCs_C2L",
                                         "cDCs_C2L",
                                         "Classical.Monocytes_C2L"))


p1<-ggplot(dend,horiz=T)
p2<-ggplot(a,aes(x=ind,y=value,fill=variable))+geom_bar(position='fill',stat='identity')+scale_y_continuous(labels = scales::percent_format())+coord_flip()+scale_fill_manual(values = colors)
plot_grid(p1, p2, align = "h")



## violinplot with enrichment score for different signatures


Il1bInA1<-rownames(SpatialA1FilteredFil)[rownames(SpatialA1FilteredFil)%in%TAMIl1b]
SpatialA1FilteredFil<-AddModuleScore(SpatialA1FilteredFil,features = Il1bInA1)

load("~/Desktop/A1_outs/CellTypeAnnotD30_Dall.Robj")

setwd("/Users/leonardi.carlo/Desktop/Spatial/Giotto")
hyperA1<-read.csv('hypergeom_A1.csv')
hyperA1<-hyperA1[,2:9]
rownames(hyperA1)<-hyperA1$cell_ID
hyperA1=hyperA1[rownames(SpatialA1FilteredFil@meta.data),]
hyperA1<-hyperA1[,2:8]

setwd("/Users/leonardi.carlo/Desktop/Spatial/Giotto")
PAGE_A1<-read.csv('PAGE_A1.csv')
PAGE_A1<-PAGE_A1[,2:9]
rownames(PAGE_A1)<-PAGE_A1$cell_ID
PAGE_A1=PAGE_A1[rownames(SpatialA1FilteredFil@meta.data),]
PAGE_A1<-PAGE_A1[,2:8]

setwd("/Users/leonardi.carlo/Desktop/Spatial/Giotto")
hyperB1<-read.csv('hypergeom_B1.csv')
hyperB1<-hyperB1[,2:9]
rownames(hyperB1)<-hyperB1$cell_ID
hyperB1=hyperB1[rownames(SpatialB1FilteredFil@meta.data),]
hyperB1<-hyperB1[,2:8]

setwd("/Users/leonardi.carlo/Desktop/Spatial/Giotto")
PAGE_B1<-read.csv('PAGE_B1.csv')
PAGE_B1<-PAGE_B1[,2:9]
rownames(PAGE_B1)<-PAGE_B1$cell_ID
PAGE_B1=PAGE_B1[rownames(SpatialB1FilteredFil@meta.data),]
PAGE_B1<-PAGE_B1[,2:8]

setwd("/Users/leonardi.carlo/Desktop/Spatial/Giotto")
hyperC1<-read.csv('hypergeom_C1.csv')
hyperC1<-hyperC1[,2:9]
rownames(hyperC1)<-hyperC1$cell_ID
hyperC1=hyperC1[rownames(SpatialC1FilteredFil@meta.data),]
hyperC1<-hyperC1[,2:8]

setwd("/Users/leonardi.carlo/Desktop/Spatial/Giotto")
PAGE_C1<-read.csv('PAGE_C1.csv')
PAGE_C1<-PAGE_C1[,2:9]
rownames(PAGE_C1)<-PAGE_C1$cell_ID
PAGE_C1=PAGE_C1[rownames(SpatialC1FilteredFil@meta.data),]
PAGE_C1<-PAGE_C1[,2:8]

setwd("/Users/leonardi.carlo/Desktop/Spatial/Giotto")
hyperD1<-read.csv('hypergeom_D1.csv')
hyperD1<-hyperD1[,2:9]
rownames(hyperD1)<-hyperD1$cell_ID
hyperD1=hyperD1[rownames(SpatialD1FilteredFil@meta.data),]
hyperD1<-hyperD1[,2:8]

setwd("/Users/leonardi.carlo/Desktop/Spatial/Giotto")
PAGE_D1<-read.csv('PAGE_D1.csv')
PAGE_D1<-PAGE_D1[,2:9]
rownames(PAGE_D1)<-PAGE_D1$cell_ID
PAGE_D1=PAGE_D1[rownames(SpatialD1FilteredFil@meta.data),]
PAGE_D1<-PAGE_D1[,2:8]


SpatialA1FilteredFil<-AddMetaData(SpatialA1FilteredFil,metadata = hyperA1,col.name = gsub('$','_hyper',colnames(hyperA1)))
SpatialB1FilteredFil<-AddMetaData(SpatialB1FilteredFil,metadata = hyperB1,col.name = gsub('$','_hyper',colnames(hyperB1)))
SpatialC1FilteredFil<-AddMetaData(SpatialC1FilteredFil,metadata = hyperC1,col.name = gsub('$','_hyper',colnames(hyperC1)))
SpatialD1FilteredFil<-AddMetaData(SpatialD1FilteredFil,metadata = hyperD1,col.name = gsub('$','_hyper',colnames(hyperD1)))

SpatialA1FilteredFil<-AddMetaData(SpatialA1FilteredFil,metadata = PAGE_A1,col.name = gsub('$','_PAGE',colnames(PAGE_A1)))
SpatialB1FilteredFil<-AddMetaData(SpatialB1FilteredFil,metadata = PAGE_B1,col.name = gsub('$','_PAGE',colnames(PAGE_B1)))
SpatialC1FilteredFil<-AddMetaData(SpatialC1FilteredFil,metadata = PAGE_C1,col.name = gsub('$','_PAGE',colnames(PAGE_C1)))
SpatialD1FilteredFil<-AddMetaData(SpatialD1FilteredFil,metadata = PAGE_D1,col.name = gsub('$','_PAGE',colnames(PAGE_D1)))

meta<-SpatialA1FilteredFil@meta.data
meta0<-SpatialB1FilteredFil@meta.data
meta1<-SpatialC1FilteredFil@meta.data
meta2<-SpatialD1FilteredFil@meta.data

metaGlobal=merge(x=SpatialA1FilteredFil,y = c(SpatialB1FilteredFil,SpatialC1FilteredFil,SpatialD1FilteredFil),add.cell.ids = c('A1','B1','C1','D1'))

class<-colnames(meta)[39:54]
metaG<-metaGlobal@meta.data

metaG_Str<-metaG%>%
  filter(TumorOrStromal=='Stromal')%>%
  nest_by(orig.ident)

metaG_Tum<-metaG%>%
  filter(TumorOrStromal=='Tumor')%>%
  nest_by(orig.ident)

test_Macro<-list()
test_Endo<-list()
test_CAF<-list()
for(i in c(1:4)){
  tot<-c(metaG_Tum[[2]][[i]]$Macrophages_C2L,metaG_Str[[2]][[i]]$Macrophages_C2L)
  tum<-rep('Tumor',length(metaG_Tum[[2]][[i]]$Macrophages_C2L))
  str<-rep('Stromal',length(metaG_Str[[2]][[i]]$Macrophages_C2L))
  data<-data.frame(value=tot,group=c(tum,str))

  test_Macro[[i]]<-bartlett.test(value~group,data)

  tot<-c(metaG_Tum[[2]][[i]]$Endothelial.cells_C2L,metaG_Str[[2]][[i]]$Endothelial.cells_C2L)
  tum<-rep('Tumor',length(metaG_Tum[[2]][[i]]$Endothelial.cells_C2L))
  str<-rep('Stromal',length(metaG_Str[[2]][[i]]$Endothelial.cells_C2L))
  data<-data.frame(value=tot,group=c(tum,str))

  test_Endo[[i]]<-bartlett.test(value~group,data)

  tot<-c(metaG_Tum[[2]][[i]]$CAF_C2L,metaG_Str[[2]][[i]]$CAF_C2L)
  tum<-rep('Tumor',length(metaG_Tum[[2]][[i]]$CAF_C2L))
  str<-rep('Stromal',length(metaG_Str[[2]][[i]]$CAF_C2L))
  data<-data.frame(value=tot,group=c(tum,str))

  test_CAF[[i]]<-bartlett.test(value~group,data)
}

test_Myeloid<-list()
test_Linfoid<-list()
for(i in c(1:4)){
  myeloid_Tum<-rowSums(metaG_Tum[[2]][[i]][,c('Classical.Monocytes_C2L','Immature.Neutrophiles_C2L','Mature.Neutrophiles_C2L','mDCs_C2L','cDCs_C2L')])
  myeloid_Str<-rowSums(metaG_Str[[2]][[i]][,c('Classical.Monocytes_C2L','Immature.Neutrophiles_C2L','Mature.Neutrophiles_C2L','mDCs_C2L','cDCs_C2L')])
  tot<-c(myeloid_Tum,myeloid_Str)
  tum<-rep('Tumor',length(myeloid_Tum))
  str<-rep('Stromal',length(myeloid_Str))
  data<-data.frame(value=tot,group=c(tum,str))

  test_Myeloid[[i]]<-bartlett.test(value~group,data)

  linfoid_Tum<-rowSums(metaG_Tum[[2]][[i]][,c('B.cells_C2L','Immature.B.cells_C2L','T.cells_C2L','NK_C2L','pDCs_C2L')])
  linfoid_Str<-rowSums(metaG_Str[[2]][[i]][,c('B.cells_C2L','Immature.B.cells_C2L','T.cells_C2L','NK_C2L','pDCs_C2L')])
  tot<-c(linfoid_Tum,linfoid_Str)
  tum<-rep('Tumor',length(myeloid_Tum))
  str<-rep('Stromal',length(myeloid_Str))
  data<-data.frame(value=tot,group=c(tum,str))

  test_Linfoid[[i]]<-bartlett.test(value~group,data)}



metaG<-data.frame(metaG,MyeloidNoMG=rowSums(metaG[,c('Classical.Monocytes_C2L','Immature.Neutrophiles_C2L','Mature.Neutrophiles_C2L','mDCs_C2L','cDCs_C2L')]))
metaG<-data.frame(metaG,Lymphoid=rowSums(metaG[,c('B.cells_C2L','Immature.B.cells_C2L','T.cells_C2L','NK_C2L','pDCs_C2L')]))


grouped_ggbetweenstats(data = metaG,
                       x=TumorOrStromal,
                       y= Macrophages_C2L,
                       grouping.var = orig.ident,
                       point.path=F,
                       type="n",
                       output='subtitle')

grouped_ggbetweenstats(data = metaG,
                       x=TumorOrStromal,
                       y= Lymphoid,
                       xlab='Tissue Area',
                       ylab='% Lymphoid',
                       grouping.var = orig.ident,
                       point.path=F,
                       type="n",
                       results.subtitle=F,
                       ggtheme=ggplot2::theme_classic())

metaG_A1<-metaG%>%
  filter(orig.ident=='A1')
metaG_B1<-metaG%>%
  filter(orig.ident=='B1')
metaG_C1<-metaG%>%
  filter(orig.ident=='C1')
metaG_D1<-metaG%>%
  filter(orig.ident=='D1')


A1<-ggbetweenstats(data = metaG_A1,
                   x=TumorOrStromal,
                   y= Cancer_C2L,
                   xlab='Tissue annotation',
                   ylab='% Cancer',
                   grouping.var = orig.ident,
                   point.path=F,
                   type="n",
                   results.subtitle=F,
                   ggtheme=ggplot2::theme_classic(),
                   title='A1',
                   plot.type = 'violin')
A1<-A1+scale_color_manual(values=c('#5D97EA','#A179E2'))+theme(axis.text.x=element_text(face='bold'))+labs(face="bold")


B1<-ggbetweenstats(data = metaG_B1,
                   x=TumorOrStromal,
                   y= Cancer_C2L,
                   xlab='Tissue annotation',
                   ylab='% Cancer',
                   grouping.var = orig.ident,
                   point.path=F,
                   type="n",
                   results.subtitle=F,
                   ggtheme=ggplot2::theme_classic(),
                   title='B1',
                   plot.type = 'violin')
B1<-B1+scale_color_manual(values=c('#5D97EA','#A179E2'))+theme(axis.text.x=element_text(face='bold'))+labs(face="bold")


C1<-ggbetweenstats(data = metaG_C1,
                   x=TumorOrStromal,
                   y= Cancer_C2L,
                   xlab='Tissue annotation',
                   ylab='% Cancer',
                   grouping.var = orig.ident,
                   point.path=F,
                   type="n",
                   results.subtitle=F,
                   ggtheme=ggplot2::theme_classic(),
                   title='C1',
                   plot.type = 'violin')
C1<-C1+scale_color_manual(values=c('#5D97EA','#A179E2'))+theme(axis.text.x=element_text(face='bold'))+labs(face="bold")

D1<-ggbetweenstats(data = metaG_D1,
                   x=TumorOrStromal,
                   y= Cancer_C2L,
                   xlab='Tissue annotation',
                   ylab='% Cancer',
                   grouping.var = orig.ident,
                   point.path=F,
                   type="n",
                   results.subtitle=F,
                   ggtheme=ggplot2::theme_classic(),
                   title='D1',
                   plot.type = 'violin')
D1<-D1+scale_color_manual(values=c('#5D97EA','#A179E2'))+theme(axis.text.x=element_text(face='bold'))+labs(face="bold")


Lymph<-combine_plots(list(A1,B1,C1,D1))

Myelo<-combine_plots(list(A1,B1,C1,D1))

Macro<-combine_plots(list(A1,B1,C1,D1))

CAF<-combine_plots(list(A1,B1,C1,D1))

Endo<-combine_plots(list(A1,B1,C1,D1))


for(i in c(1:length(class))){
  grouped_ggwithinstats(data = metaG,
                        x=TumorOrStromal,
                        y= !!class[i],
                        grouping.var = orig.ident,
                        point.path=F,
                        type="n",
                        results.subtitle=T)
  dev.copy2pdf(file=paste(class[i],'_stat.pdf',sep=''), width = 10, height = 10)
}


p1<-ggplot(meta,aes(x=TumorOrStromal,y=mDCs_C2L))+
  geom_violin(aes(fill=TumorOrStromal),trim = TRUE,scale = "width", show.legend = FALSE,alpha=0.5)+
  geom_jitter(size=1,aes(alpha=0.6))+ggtitle('A1')+theme_classic()+theme(legend.position = "none")

p2<-ggplot(meta0,aes(x=TumorOrStromal,y=mDCs_C2L))+
  geom_violin(aes(fill=TumorOrStromal),trim = TRUE,scale = "width", show.legend = FALSE,alpha=0.5)+
  geom_jitter(size=1,aes(alpha=0.6))+ggtitle('B1')+theme_classic()+theme(legend.position = "none")

p3<-ggplot(meta1,aes(x=TumorOrStromal,y=mDCs_C2L))+
  geom_violin(aes(fill=TumorOrStromal),trim = TRUE,scale = "width", show.legend = FALSE,alpha=0.5)+
  geom_jitter(size=1,aes(alpha=0.6))+ggtitle('C1')+theme_classic()+theme(legend.position = "none")

p4<-ggplot(meta2,aes(x=TumorOrStromal,y=mDCs_C2L))+
  geom_violin(aes(fill=TumorOrStromal),trim = TRUE,scale = "width", show.legend = FALSE,alpha=0.5)+
  geom_jitter(size=1,aes(alpha=0.6))+ggtitle('D1')+theme_classic()+theme(legend.position = "none")
  
ggarrange(p1,p2,p3,p4,ncol = 4,nrow=1)



metaT<-meta%>%
  filter(meta$TumorOrStromal=='Tumor')

metaS<-meta%>%
  filter(meta$TumorOrStromal=='Stromal')



ggplot(meta,aes(x=TAMIl1b_hyper,fill=TumorOrStromal,alpha=0.5))+
  geom_freqpoly(bins = 100,position = 'identity',aes(col=TumorOrStromal))+
  theme_classic()

ggplot(meta,aes(x=TAMIl1b_hyper,fill=TumorOrStromal,alpha=0.3))+
  geom_histogram(bins = 100,position = 'identity')+
  theme_classic()


normS<-metaS$TAMIl1b_hyper/dim(metaS)[1]
normT<-metaT$TAMIl1b_hyper/dim(metaT)[1]
names(normS)=rownames(metaS)
names(normT)=rownames(metaT)
norm<-c(normS,normT)

norm<-norm[rownames(meta)]


p1<-ggplot(meta_01,aes(x=TumorOrStromal,y=Macrophages_C2L))+
  geom_violin(aes(fill=TumorOrStromal),trim = TRUE,scale = "width", show.legend = FALSE,alpha=0.5)+
  geom_jitter(size=1,aes(color=TAMIl1b_hyper,alpha=0.6))+
  scale_colour_gradient(low="#ADD8E6",high="red")+theme_classic()

p2<-ggplot(meta_02,aes(x=TumorOrStromal,y=Macrophages_C2L))+
  geom_violin(aes(fill=TumorOrStromal),trim = TRUE,scale = "width", show.legend = FALSE,alpha=0.5)+
  geom_jitter(size=1,aes(color=TAMIl1b_hyper,alpha=0.6))+
  scale_colour_gradient(low="#ADD8E6",high="red")+theme_classic()

p3<-ggplot(meta_02up,aes(x=TumorOrStromal,y=Macrophages_C2L))+
  geom_violin(aes(fill=TumorOrStromal),trim = TRUE,scale = "width", show.legend = FALSE,alpha=0.5)+
  geom_jitter(size=1,aes(color=TAMIl1b_hyper,alpha=0.6))+
  scale_colour_gradient(low="#ADD8E6",high="red")+theme_classic()


meta_01S<-meta%>%
  filter(Macrophages_C2L<=0.1)%>%
  filter(TumorOrStromal=='Stromal')
meta_01T<-meta%>%
  filter(Macrophages_C2L<=0.1)%>%
  filter(TumorOrStromal=='Tumor')

meta_02S<-meta%>%
  filter(Macrophages_C2L>0.1)%>%
  filter(TumorOrStromal=='Stromal')
meta_02T<-meta%>%
  filter(Macrophages_C2L>0.1)%>%
  filter(TumorOrStromal=='Tumor')

#dend<-as.dendrogram(hclust(dist(meta_01S[,c(39:54)])))
#order<-labels(dend)
a<-meta_02S[1:20,c(39:54)]
order<-meta_02S$TAMIl1b_hyper
names(order)=rownames(meta_02S)
order<-order[order(order,decreasing = T)]
a<-a[names(order),]
meta_02S$Cells<-rownames(meta_02S[names(order),])
meta_02S<-meta_02S[names(order),]

a1<-meta_02S[1:20,c(39:54)]
order<-meta_02S$TAMSFolr2_hyper
names(order)=rownames(meta_02S)
order<-order[order(order,decreasing = T)]
a1<-a1[names(order),]
meta_02S$Cells<-rownames(meta_02S[names(order),])
meta_02S<-meta_02S[names(order),]

b<-meta_02T[,c(39:54)]
order<-meta_02T$TAMIl1b_hyper
names(order)=rownames(meta_02T)
order<-order[order(order,decreasing = T)]
b<-b[names(order),]
meta_02T<-meta_02T[names(order),]
meta_02T$Cells<-rownames(meta_02T)

#a<-a[order,]
a<-a %>%
  mutate(ind = factor(row_number()))%>%
  gather(variable,value,-ind)
a$variable<-factor(a$variable,levels = c('Macrophages_C2L',
                                         'CAF_C2L',
                                         "Endothelial.cells_C2L",
                                         "Cancer_C2L",
                                         "Ductal_C2L",
                                         "Acinar.cells_C2L",
                                         "T.cells_C2L",
                                         "NK_C2L",
                                         "B.cells_C2L",
                                         "Immature.B.cells_C2L",
                                         "pDCs_C2L",
                                         "Mature.Neutrophiles_C2L",
                                         "Immature.Neutrophiles_C2L",
                                         "mDCs_C2L",
                                         "cDCs_C2L",
                                         "Classical.Monocytes_C2L"))

a1<-a1 %>%
  mutate(ind = factor(row_number()))%>%
  gather(variable,value,-ind)
a1$variable<-factor(a1$variable,levels = c('Macrophages_C2L',
                                           'CAF_C2L',
                                           "Endothelial.cells_C2L",
                                           "Cancer_C2L",
                                           "Ductal_C2L",
                                           "Acinar.cells_C2L",
                                           "T.cells_C2L",
                                           "NK_C2L",
                                           "B.cells_C2L",
                                           "Immature.B.cells_C2L",
                                           "pDCs_C2L",
                                           "Mature.Neutrophiles_C2L",
                                           "Immature.Neutrophiles_C2L",
                                           "mDCs_C2L",
                                           "cDCs_C2L",
                                           "Classical.Monocytes_C2L"))

b<-b %>%
  mutate(ind = factor(row_number()))%>%
  gather(variable,value,-ind)
b$variable<-factor(b$variable,levels = c('Macrophages_C2L',
                                         'CAF_C2L',
                                         "Endothelial.cells_C2L",
                                         "Cancer_C2L",
                                         "Ductal_C2L",
                                         "Acinar.cells_C2L",
                                         "T.cells_C2L",
                                         "NK_C2L",
                                         "B.cells_C2L",
                                         "Immature.B.cells_C2L",
                                         "pDCs_C2L",
                                         "Mature.Neutrophiles_C2L",
                                         "Immature.Neutrophiles_C2L",
                                         "mDCs_C2L",
                                         "cDCs_C2L",
                                         "Classical.Monocytes_C2L"))

#p1<-ggplot(dend,horiz=T)
p4<-ggplot(a,aes(x=ind,y=value,fill=variable))+
  geom_bar(position='fill',stat='identity',show.legend = T)+
  scale_y_continuous(labels = scales::percent_format())+
  coord_flip()+
  scale_fill_manual(values = colors)+
  ggtitle('Macrophages > 0.1%, Stroma - Top 20 Il1b enrich.')

p4_1<-ggplot(a1,aes(x=ind,y=value,fill=variable))+
  geom_bar(position='fill',stat='identity',show.legend = T)+
  scale_y_continuous(labels = scales::percent_format())+
  coord_flip()+
  scale_fill_manual(values = colors)+
  ggtitle('Macrophages > 0.1%, Stroma - Top 20 Folr2 enrich.')

p5<-ggplot(b,aes(x=ind,y=value,fill=variable))+
  geom_bar(position='fill',stat='identity')+
  scale_y_continuous(labels = scales::percent_format())+
  coord_flip()+
  scale_fill_manual(values = colors)+
  ggtitle('Macrophages > 0.1%, Tumor')


p6<-ggplot(meta_02S[1:20,], aes(x=seq(TAMIl1b_hyper), y=-TAMIl1b_hyper)) + geom_point()+coord_flip()+theme_classic()


stats<-data.frame(Celltype=colnames(meta_02S[1:20,c(39:54)]),MeanTop20Folr2=colMeans(meta_02S[1:20,c(39:54)]),MedianTop20Folr2=colMedians(as.matrix(meta_02S[1:20,c(39:54)])))
stats<-data.frame(stats,MeanTop20Il1b=colMeans(meta_02S[1:20,c(39:54)]),MedianTop20Il1b=colMedians(as.matrix(meta_02S[1:20,c(39:54)])))

WriteXLS(stats,ExcelFileName = 'stats_top20.xls')































