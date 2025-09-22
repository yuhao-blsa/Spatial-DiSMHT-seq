library(Seurat)
library(ggplot2)

# Ke LoH + Adult
color.ls5a = c("#7CAE00","#ffdf00","#AD0000","#7F7F7F","#C09B00","#ff4500","#00BAE0","#A3A500","#000080","#90EE90",
               "#EA8331","#00B0F6","#FF0000","#704241","#FFFF00","#C77CFF","#daa520","#adff2f","#8D2282","#E76BF3",
               "#39B600","#9590FF","#FF6A98","#FA62DB","#4169E1","#FF9999","#B22222","#FF62BC","#FFD032","#A2FF00",
               "#006400","#c0c0c0","#004A00","#FFB676","#87CEFA","#B34E00","#63E463","#6A5ACD")
names(color.ls5a) = 0:37
color.ls5a2 = c("#7CAE00","#ffdf00","#AD0000","#6A5ACD","#C09B00","#ff4500","#00BAE0","#A3A500","#000080","#90EE90",
                "#EA8331","#00B0F6","#FF0000","#704241","#FFFF00","#C77CFF","#daa520","#adff2f","#8D2282","#E76BF3",
                "#39B600","#9590FF","#FF6A98","#FA62DB","#4169E1","#FF9999","#B22222","#FF62BC","#FFD032","#A2FF00",
                "#006400","#c0c0c0","#004A00","#FFB676","#87CEFA","#B34E00","#63E463","#6A5ACD")
names(color.ls5a2) = c(0:37)
# load data
GW6_total = readRDS("Plot/1_SCT_transform&Merged_data.add_C17.filter.RDS") # "Plot/1_SCT_transform&Merged_data.add_C17.RDS"
Sample.ls = c("GW6_1","GW6_2")
size = 0.4;nchannel = 192

# 0 Spatial map svg
for (Samplei in Sample.ls) {
  svg(paste0("Plot_svg/1_spatial_map.",Samplei,".svg"))
  p1 = ggplot(GW6_total@meta.data[GW6_total$orig.ident==Samplei,],aes(x=nB,y=nA,color=C18)) +
    geom_point(shape=15,size=0.6,show.legend=F) + theme_void() +
    scale_y_continuous(trans = "reverse",breaks = seq(0,192,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,192,5),name = paste("nB",Samplei)) +
    expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() + scale_color_manual(values = color.ls5a) +
    theme(panel.background = element_rect(fill = "black"))
  print(p1)
  dev.off()
}

## 2 Run in shell ##
# Extract coordinate
# cd /Plot_svg/
# for si in GW6_1 GW6_2
# do
# cat 1_spatial_map.${si}.ai.svg | grep "path class" | sed 's/^.*M\([0-9]*\.\?[0-9]*,[0-9]*\.\?[0-9]*\).*$/\1/g' | sed 's/,/ /g' > 1_spatial_map.${si}.ai.coord
# cat 2_spatial_map.${si}.ai.filter.svg | grep "path class" | sed 's/^.*M\([0-9]*\.\?[0-9]*,[0-9]*\.\?[0-9]*\).*$/\1/g' | sed 's/,/ /g' > 2_spatial_map.${si}.ai.filter.coord
# done

## 3 Save in Seurat object ##
final_cell.ls = c()
for (si in Sample.ls) {
  cord.df = read.table(paste0("Plot_svg/1_spatial_map.",si,".ai.coord"))
  cord_crop.df = read.table(paste0("Plot_svg/2_spatial_map.",si,".ai.filter.coord"))
  
  cord.df$cellid = row.names(GW6_total@meta.data)[GW6_total$orig.ident==si]
  cord.df$isfinal = ifelse(paste0(cord.df$V1,"_",cord.df$V2) %in% 
                             paste0(cord_crop.df$V1,"_",cord_crop.df$V2),1,0)
  final_cell.ls = c(final_cell.ls,cord.df$cellid[cord.df$isfinal==1])
}
GW6_total$isfinal = 0
GW6_total@meta.data[final_cell.ls,"isfinal"] = 1

saveRDS(GW6_total,file = "Plot/1_SCT_transform&Merged_data.add_C17.filter.RDS")
write.csv(GW6_total@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","channel","stage","coord","nA","nB",
                                 "nCount_SCT","nFeature_SCT","seurat_clusters","C18","isfinal")],file = "Plot/1_SCT.cell_metadata.csv")

## Validate result #
GW6_total = GW6_total[,GW6_total$isfinal == 1]
pdf("Plot/9_cluster_UMAP.add_C17.filter2.pdf")
DimPlot(GW6_total, reduction = "umap", label = T,shuffle = T,group.by = "C18") + scale_color_manual(values = color.ls5a2)
DimPlot(GW6_total, reduction = "umap", label = T,shuffle = T,group.by = "orig.ident")
dev.off()
pdf("Plot/9_spatial_map.cluster.add_C17.filter.pdf")
for (Samplei in Sample.ls) {
  p1 = ggplot(GW6_total@meta.data[GW6_total$orig.ident==Samplei,],aes(x=nB,y=nA,color=C18)) +
    geom_point(shape=19,size=size,show.legend=F) + theme_void() +
    scale_y_continuous(trans = "reverse",breaks = seq(0,192,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,192,5),name = paste("nB",Samplei)) +
    expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() + scale_color_manual(values = color.ls5a) +
    theme(panel.background = element_rect(fill = "black"))
  print(p1)
}
dev.off()


# c("DNMT1","DNMT3A","DNMT3B","DNMT3L","UHRF1","UHRF2","MBD1","MBD2","MBD3","MBD4","MECP2",
#   "TET1","TET2","TET3","TDG","SMUG1","APEX1","APEX2","IDH1","IDH2","OGDH","SDHA","SDHB","SDHC","SDHD",
#   "GADD45A","KDM1A","MIR29A","MIR29B1","MIR29B2","MIR29C","HOTAIR","USP7","PHF20","PHF20L1")
GW6_total@assays$SCT@data = log1p(GW6_total@assays$SCT@data)

Idents(GW6_total) = factor(GW6_total$C18,levels = c(13,6,9,15,16,14,5,1,11,2,7,10,12,4,17,0,3,8))
GW6_total$tmp = factor(GW6_total$C18,levels = c(13,6,9,15,16,14,5,1,11,2,7,10,12,4,17,0,3,8))
levels(GW6_total$tmp) = c("C13"="craniofacial","C6"="forebrain","C9"="forebrain","C15"="midbrain","C16"="midbrain",
                          "C14"="midbrain","C5"="hindbrain","C1"="spinalcord","C11"="spinalcord","C2"="spinalcord",
                          "C7"="mesenchyme","C10"="forelimb","C12"="forelimb","C4"="heart","C17"="heart",
                          "C0"="hepatoblast","C3"="hindlimb","C8"="hindlimb")
Idents(GW6_total) = GW6_total$tmp
pdf("Plot/10_DNMT_expression.vlnplot.pdf",width = 10)
VlnPlot(GW6_total, features = c("DNMT1","DNMT3A","DNMT3B","DNMT3L","UHRF1","UHRF2","MBD1","MBD2","MBD3","MBD4","MECP2",
                                "TET1","TET2","TET3","TDG","SMUG1","APEX1","APEX2","IDH1","IDH2","OGDH","SDHA","SDHB","SDHC","SDHD",
                                "GADD45A","KDM1A","MIR29A","MIR29B1","MIR29B2","MIR29C","HOTAIR","USP7","PHF20","PHF20L1"),
        layer = "data",stack = T,flip = T,log = F,adjust = 1) + NoLegend()
dev.off()
pdf("Plot/10_DNMT_expression.dotplot.pdf",height = 6,width = 10)
DotPlot(GW6_total,assay = "SCT", features = c("DNMT1","DNMT3A","DNMT3B","DNMT3L","UHRF1","UHRF2","MBD1","MBD2","MBD3","MBD4","MECP2",
                                              "TET1","TET2","TET3","TDG","SMUG1","APEX1","APEX2","IDH1","IDH2","OGDH","SDHA","SDHB","SDHC","SDHD",
                                              "GADD45A","KDM1A","MIR29A","MIR29B1","MIR29B2","MIR29C","HOTAIR","USP7","PHF20","PHF20L1"), dot.scale = 6,scale = T) + RotatedAxis()
dev.off()
pdf("Plot/10a_DNMT_expression.dotplot.pdf",height = 3.5,width = 5.5)
DotPlot(GW6_total,assay = "SCT", features = c("DNMT1","DNMT3A","DNMT3B","DNMT3L","UHRF1","UHRF2","TET1","TET2","TET3"), dot.scale = 6,scale = T) + RotatedAxis() +
  scale_color_gradientn(colours=c("#F3F1E4FF","#F2F2DEFF","#EEF2DAFF","#EAF2D6FF","#E4F2D3FF","#DEF1D0FF","#D6F0CDFF","#CDEECAFF","#C4EBC8FF","#BAE9C6FF","#AFE5C4FF","#A3E2C3FF",
                                  "#97DDC1FF","#8BD8C0FF","#7DD2BFFF","#70CCBDFF","#62C5BCFF","#54BEBAFF","#45B6B8FF","#37ADB6FF","#28A4B3FF","#189AAFFF","#068FABFF","#0084A7FF",
                                  "#0078A2FF","#046C9CFF","#115F96FF","#1C518FFF","#254289FF","#2D3184FF"))
dev.off()

# library(paletteer)
for (genei in c("DNMT1","DNMT3A","DNMT3B","UHRF1","UHRF2","TET1","TET2","TET3","OGDH")){
  pdf(paste0("Plot/11_spatial_map.",genei,".SCT_data.min0.8.pdf"),width = 7,height = 7)
  # 1
  p1 = FeaturePlot(GW6_total, features = genei,slot = "data",order = T,max.cutoff = "q99")
  print(p1)
  # 2
  GW6_total$focus = GW6_total@assays$SCT@data[genei,]
  maxi = sort(GW6_total$focus,decreasing = T)[10]
  mini = 0.8
  GW6_total$focus = pmin(GW6_total$focus,maxi)
  GW6_total$focus = pmax(GW6_total$focus,mini)
  for (Samplei in Sample.ls) {
    p1 = ggplot(GW6_total@meta.data[GW6_total$orig.ident==Samplei,],aes(x=nB,y=nA,color=focus)) +
      geom_point(shape=19,size=0.4,show.legend=T) + theme_void() +
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",Samplei)) +
      expand_limits(y=c(0,96),x=c(0,96)) + coord_fixed() + 
      scale_color_gradientn(colours=c("#313695FF","#4575B4FF","#74ADD1FF","#ABD9E9FF","#E0F3F8FF","#FFFFBFFF","#FEE090FF","#FDAE61FF","#F46D43FF","#D73027FF","#A50026FF")) + 
      theme(panel.background = element_rect(fill = "black"))
    print(p1)
  }
  dev.off()
}


Idents(GW6_total) = GW6_total$C18
GW6_total <- JoinLayers(GW6_total,assay = "RNA")
GW6_total <- NormalizeData(GW6_total,assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000)
de_markers <- FindAllMarkers(GW6_total, assay = "RNA", only.pos = T,min.pct = 0.1,logfc.threshold = 0.5,min.diff.pct = 0.1)
write.csv(de_markers,row.names = F,file = "Plot/6_de_genes.total.csv")
de_markers <- FindAllMarkers(GW6_total, assay = "RNA", only.pos = T,min.pct = 0.1,logfc.threshold = 0.5,min.diff.pct = 0)
write.csv(de_markers,row.names = F,file = "Plot/6_de_genes.total_soft.csv")


Idents(GW6_total) = factor(GW6_total$C18,levels = c(13,6,9,15,16,14,5,1,11,2,7,10,12,4,17,0,3,8))
pdf("Plot/8_dot_plot.all_cluster_tmp.pdf",height = 6,width = 18)
DotPlot(GW6_total,assay = "RNA", features = c("ALX1","LHX2","NOL4","NR2E1","DLX6.AS1","SOX2.OT","RMST","MAP2","MAP1B","SYT1","GPM6A","NRCAM","CTNNA2","MAPK10",
                                              "EMX1","EMX2","SP8","ZIC2","OTX2","TCF7L2","LEF1","LMX1A",  "SHH","CMTM8","SLIT2","RFX4","RNF220","DCC","NTN1","SLIT1","SEMA5B","PLEKHA5"  ,"CRABP1","ELAVL4",
                                              "NOVA1","NTRK3","NEB","MYH3","PTCH1","SULF1",
                                             "PRRX1","DLX5","TBX5","ACTC1","MYL7","HBA1","HBA2","AFP","ALB","TBX4"), dot.scale = 6,scale = T) + RotatedAxis() + 
  scale_color_gradientn(colours=c("#F3F1E4FF","#F2F2DEFF","#EEF2DAFF","#EAF2D6FF","#E4F2D3FF","#DEF1D0FF","#D6F0CDFF","#CDEECAFF","#C4EBC8FF","#BAE9C6FF","#AFE5C4FF","#A3E2C3FF",
                                  "#97DDC1FF","#8BD8C0FF","#7DD2BFFF","#70CCBDFF","#62C5BCFF","#54BEBAFF","#45B6B8FF","#37ADB6FF","#28A4B3FF","#189AAFFF","#068FABFF","#0084A7FF",
                                  "#0078A2FF","#046C9CFF","#115F96FF","#1C518FFF","#254289FF","#2D3184FF"))
dev.off()
pdf("Plot/8_dot_plot.all_cluster_final.pdf",height = 6,width = 8)
DotPlot(GW6_total,assay = "RNA", features = unique(c("ALX1","LHX2","NOL4","EMX2","OTX2","LEF1","LMX1A","NTN1","NTRK3","MYH3","NPR3","DCN","PTCH1","DLX5","TBX5","ACTC1","MYH6","ALB","TBX4","COL1A1","ADAMTS20")), dot.scale = 6,scale = T,scale.by = "size") + RotatedAxis() + 
  scale_color_gradientn(colours=c("#F3F1E4FF","#F2F2DEFF","#EEF2DAFF","#EAF2D6FF","#E4F2D3FF","#DEF1D0FF","#D6F0CDFF","#CDEECAFF","#C4EBC8FF","#BAE9C6FF","#AFE5C4FF","#A3E2C3FF",
                                  "#97DDC1FF","#8BD8C0FF","#7DD2BFFF","#70CCBDFF","#62C5BCFF","#54BEBAFF","#45B6B8FF","#37ADB6FF","#28A4B3FF","#189AAFFF","#068FABFF","#0084A7FF",
                                  "#0078A2FF","#046C9CFF","#115F96FF","#1C518FFF","#254289FF","#2D3184FF"))
DotPlot(GW6_total,assay = "SCT", features = c("ALX1","LHX2","NOL4","EMX2","OTX2","LEF1","LMX1A","NTN1","NTRK3","MYH3","NPR3","DCN","PTCH1","DLX5","TBX5","ACTC1","MYH6","ALB","TBX4","COL1A1","ADAMTS20"), dot.scale = 6,scale = T,scale.by = "size") + RotatedAxis() + 
  scale_color_gradientn(colours=c("#F3F1E4FF","#F2F2DEFF","#EEF2DAFF","#EAF2D6FF","#E4F2D3FF","#DEF1D0FF","#D6F0CDFF","#CDEECAFF","#C4EBC8FF","#BAE9C6FF","#AFE5C4FF","#A3E2C3FF",
                                  "#97DDC1FF","#8BD8C0FF","#7DD2BFFF","#70CCBDFF","#62C5BCFF","#54BEBAFF","#45B6B8FF","#37ADB6FF","#28A4B3FF","#189AAFFF","#068FABFF","#0084A7FF",
                                  "#0078A2FF","#046C9CFF","#115F96FF","#1C518FFF","#254289FF","#2D3184FF"))
dev.off()

## export cell type average expression ##
GW6_total <- JoinLayers(GW6_total,assay = "RNA")
agg_ct <- as.matrix(AggregateExpression(GW6_total, slot = "counts", assays = "RNA", group.by = c("C18"), return.seurat = F)[[1]])
save(list = c("agg_ct"),file = "Plot/12_agg_ct.count.Rdata")
library(preprocessCore)
agg_ct.quantile = normalize.quantiles(agg_ct,copy=TRUE,keep.names = T)
save(list = c("agg_ct.quantile"),file = "Plot/12_agg_ct.quantile.Rdata")

## export cell type average expression - merged cell types ##
GW6_total@meta.data$merge_ct = c("C13"="craniofacial","C6"="forebrain","C9"="forebrain","C15"="midbrain","C16"="midbrain",
                                 "C14"="midbrain","C5"="hindbrain","C1"="spinalcord","C11"="spinalcord","C2"="spinalcord",
                                 "C7"="mesenchyme","C10"="forelimb","C12"="forelimb","C4"="heart","C17"="heart",
                                 "C0"="hepatoblast","C3"="hindlimb","C8"="hindlimb")[paste0("C",GW6_total$C18)]
GW6_total <- JoinLayers(GW6_total,assay = "RNA")
agg_ct <- as.matrix(AggregateExpression(GW6_total, slot = "counts", assays = "RNA", group.by = c("merge_ct"), return.seurat = F)[[1]])
library(preprocessCore)
agg_ct.quantile = normalize.quantiles(agg_ct,copy=TRUE,keep.names = T)
save(list = c("agg_ct.quantile"),file = "Plot/12_agg_ct.quantile.merged_ct.Rdata")

