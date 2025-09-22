options(future.globals.maxSize = 5*1024^3)

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(Matrix)


LoadPixelInTissue <- function(i,norm=T){
  print(Sample.ls[i])
  count.m = read.table(data_path.ls[i])
  count.m = as.matrix(count.m)
  count.m = count.m[,apply(count.m, 2, sum)>50]
  tmp = CreateSeuratObject(counts = t(count.m), project = Sample.ls[i], assay = "RNA")
  tmp$channel = channel.ls[i]
  tmp$stage = Stage.ls[i]
  tmp$coord = paste0(Sample.ls[i],"_",row.names(tmp@meta.data))
  tmp$nA = as.integer(str_split_fixed(tmp$coord,"_|x",4)[,3])
  tmp$nB = as.integer(str_split_fixed(tmp$coord,"_|x",4)[,4])
  
  tmp.subset = tmp[,tmp$nCount_RNA > Filter.m[i]]
  if(norm) tmp.subset <- SCTransform(tmp.subset, assay = "RNA", verbose = T,variable.features.n = 3000)
  return(tmp.subset)
}
# 1 Merge data
data_path.ls = c("RNA_exp4_250407_HEM2_stdata.symbol.tsv",
                 "RNA_exp4_250407_HEM3_stdata.symbol.tsv")
channel.ls = c("20um","20um")
Stage.ls = c("GW6","GW6")
Sample.ls = c("GW6_1","GW6_2")
Filter.m = c(1500,1500)
nchannel = 192
size = 0.4 # 1.2

# Run 1 by 1
for (i in seq(Sample.ls)) {
  assign(Sample.ls[i],LoadPixelInTissue(i,norm=F))
}
pdf("Plot/0_read_distribution.filtered.1500.pdf",width = 7,height = 7) ## Final ##
for (i in seq(Sample.ls)) {
  tmp.df = get(Sample.ls[i])@meta.data
  print(ggplot(tmp.df[tmp.df$nCount_RNA>1500,],aes(x=nB,y=nA,color=nCount_RNA)) +
          geom_point(shape=19,size = size,show.legend=T) + theme_void() +
          scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
          scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",Sample.ls[i])) +
          expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() + 
          scale_color_gradientn(colors=c("#313695FF","#4575B4FF","#74ADD1FF","#ABD9E9FF","#E0F3F8FF","#FFFFBFFF","#FEE090FF","#FDAE61FF","#F46D43FF","#D73027FF","#A50026FF"),
                                limits = c(min(tmp.df$nCount_RNA), quantile(tmp.df$nCount_RNA,0.99))) +
          theme(panel.background = element_rect(fill = "black")))
}
dev.off()
for (i in seq(Sample.ls)) {
  assign(Sample.ls[i],LoadPixelInTissue(i,norm=T))
}

GW6_total <- merge(get(Sample.ls[1]), y = unlist(lapply(2:length(Sample.ls),function(i) get(Sample.ls[i]))),
                    add.cell.ids = Sample.ls, project = "GW6_total", merge.data = T)
rm(list = c(Sample.ls,"i","tmp.df","LoadPixelInTissue","Stage.ls","channel.ls","data_path.ls","Filter.m"))


# 2 Clustering
# Seurat Issue 2814/5183
set.seed(2023)
VariableFeatures(GW6_total[["SCT"]]) <- rownames(GW6_total[["SCT"]]@scale.data) # 3733

GW6_total <- RunPCA(GW6_total, assay = "SCT", verbose = T,npcs = 50)
saveRDS(GW6_total,file = "Plot/1_SCT_transform&Merged_data.RDS")

k1 = 30;k2 = 0.3
GW6_total <- FindNeighbors(GW6_total, reduction = "pca", dims = 1:k1)
GW6_total <- FindClusters(GW6_total, verbose = T, resolution = k2)
GW6_total <- RunUMAP(GW6_total, reduction = "pca", dims = 1:k1, return.model = TRUE)

# Ke LoH + Adult
color.ls5a = c("#7CAE00","#ffdf00","#AD0000","#7F7F7F","#C09B00","#ff4500","#00BAE0","#A3A500","#000080","#90EE90",
               "#EA8331","#00B0F6","#FF0000","#704241","#FFFF00","#C77CFF","#daa520","#adff2f","#8D2282","#E76BF3",
               "#39B600","#9590FF","#FF6A98","#FA62DB","#4169E1","#FF9999","#B22222","#FF62BC","#FFD032","#A2FF00",
               "#006400","#c0c0c0","#004A00","#FFB676","#87CEFA","#B34E00","#63E463","#6A5ACD")
names(color.ls5a) = 0:37
pdf(paste0("Plot/2_cluster_UMAP.PCA",k1,"_res",k2,".pdf"))
DimPlot(GW6_total, reduction = "umap", label = T,shuffle = T) + scale_color_manual(values = color.ls5a)
DimPlot(GW6_total, reduction = "umap", label = T,shuffle = T,group.by = "orig.ident")
dev.off()
pdf(paste0("Plot/2a_cluster_UMAP.PCA",k1,"_res",k2,".sub.pdf"))
color.ls5b = color.ls5a
color.ls5b["3"] = "red"
for (i in levels(GW6_total$seurat_clusters)) {
  GW6_total$focus = GW6_total$seurat_clusters
  GW6_total$focus[GW6_total$focus != i] = NA
  print(DimPlot(GW6_total, reduction = "umap",group.by = "focus", label = T,shuffle = T) + scale_color_manual(values = color.ls5b))
}
dev.off()
pdf("Plot/2_cluster_UMAP.split.pdf",height = 5,width = 8)
DimPlot(GW6_total, reduction = "umap",group.by = "seurat_clusters",label = T, split.by = "orig.ident",label.size = 3) + scale_color_manual(values = color.ls5a)
dev.off()

# 3 Spatial map
Sample.ls = c("GW6_1","GW6_2")
pdf("Plot/3_spatial_map.cluster.pdf")
for (Samplei in Sample.ls) {
  p1 = ggplot(GW6_total@meta.data[GW6_total$orig.ident==Samplei,],aes(x=nB,y=nA,color=seurat_clusters)) +
    geom_point(shape=19,size=size,show.legend=F) + theme_void() +
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",Samplei)) +
    expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() + scale_color_manual(values = color.ls5a) +
    theme(panel.background = element_rect(fill = "black"))
  print(p1)
}
dev.off()
for (Ci in levels(GW6_total$seurat_clusters)) {
  pdf(paste0("Plot/3_splited_spatial_map.cluster",Ci,".pdf"))
  GW6_total$focus = GW6_total$seurat_clusters
  GW6_total$focus[GW6_total$focus != Ci] = NA
  for (Samplei in Sample.ls) {
    if(length(GW6_total$orig.ident==Samplei & GW6_total$seurat_clusters==Ci)==0) next
    p1 = ggplot(GW6_total@meta.data[GW6_total$orig.ident==Samplei,],aes(x=nB,y=nA,color=focus)) +
      geom_point(shape=19,size=size,show.legend=F) + theme_void() +
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",Samplei,Ci)) +
      expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() + scale_color_manual(values = color.ls5b)+
      theme(panel.background = element_rect(fill = "black"))
    print(p1)
  }
  dev.off()
}

# 4 Known markers
GW6_total$mean_depth_sample = tapply(GW6_total$nCount_SCT,GW6_total$orig.ident,mean)[GW6_total$orig.ident]
GW6_total@assays[["SCT"]]@data = Matrix(t(t(as.matrix(GW6_total@assays[["SCT"]]@counts))/GW6_total$mean_depth_sample*10000),sparse=T)
saveRDS(GW6_total,file = "Plot/1_SCT_transform&Merged_data.RDS")
write.csv(GW6_total@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","channel","stage","coord","nA","nB",
                                "nCount_SCT","nFeature_SCT","seurat_clusters")],file = "Plot/1_SCT.cell_metadata.csv")
GW6_total = readRDS("Plot/1_SCT_transform&Merged_data.RDS")
Sample.ls = c("GW6_1","GW6_2")
size = 0.4;nchannel = 192

pdf("Plot/4_Vln_UMI_data.pdf",width = 10,height = 5)
VlnPlot(GW6_total, features = "nCount_RNA")
VlnPlot(GW6_total, features = "nCount_RNA",log = T)
VlnPlot(GW6_total, features = "nCount_RNA",y.max = 25000)
dev.off()
# 5 Known markers
system("mkdir Plot/Known_marker")
refgene.df = read.table("../RefMarker/human_early_organogenesis_markers.total.tsv",sep = "\t",stringsAsFactors = F)
gene.ls = row.names(GW6_total@assays$SCT@counts)[apply(GW6_total@assays$SCT@counts,1,sum)>0]
for(i in 1:dim(refgene.df)[1]){
  genei = strsplit(refgene.df$V2[i],", |,")[[1]]
  genei = intersect(genei,gene.ls)
  if(length(genei)==0) next()
  pdf(paste0("Plot/Known_marker/",refgene.df$V1[i],"_",refgene.df$V4[i],".",refgene.df$V3[i],".pdf"),width = 10,height = 10)
  for (j in 1:ceiling(length(genei)/4)) {
    genej = genei[(4*j-3):min(4*j,length(genei))]
    if(length(genej)<4) genej = c(genej,c("SMOC2","NRXN3","IGFBP5")[1:(4-length(genej))])
    print(FeaturePlot(GW6_total, features = genej,pt.size = 3,slot = "data",order = T,max.cutoff = "q99",raster = T,raster.dpi = c(1024, 1024)))
  }
  for (j in 1:ceiling(length(genei)/10)) {
    genej = genei[(10*j-9):min(10*j,length(genei))]
    if(length(genej)==1) genej = c(genej,"SMOC2","NRXN3","IGFBP5")
    # genej = c(genej,rep(tail(genej,n=1),10-length(genej)))
    print(VlnPlot(GW6_total, features = genej,slot = "data",stack = T,flip = T))
  }
  dev.off()
}

# 6 DE genes
GW6_total <- JoinLayers(GW6_total,assay = "RNA")
GW6_total <- NormalizeData(GW6_total,assay = "RNA",normalization.method = "LogNormalize", scale.factor = 10000)
de_markers <- FindAllMarkers(GW6_total, assay = "RNA", only.pos = T,min.pct = 0.1,logfc.threshold = 0.5,min.diff.pct = 0.1)
de_markers %>% filter(!grepl("^ENSG[0-9]+$", gene) & !grepl("^LINC[0-9]+$", gene)) %>% group_by(cluster) %>% slice_max(n=6,order_by = avg_log2FC) -> top6
de_markers %>% filter(!grepl("^ENSG[0-9]+$", gene) & !grepl("^LINC[0-9]+$", gene)) %>% group_by(cluster) %>% slice_head(n=6) -> top6
write.csv(de_markers,row.names = F,file = "Plot/6_de_genes.csv")
# de_markers = read.csv("Plot/6_de_genes.csv",stringsAsFactors=F)
# de_markers_include_down <- FindAllMarkers(GW6_total, assay = "RNA", only.pos = F,min.pct = 0.1,logfc.threshold = 0.5)
# de_markers_include_down = de_markers_include_down[de_markers_include_down$avg_log2FC < 0,]
# de_markers_include_down = de_markers_include_down[(de_markers_include_down$pct.1 - de_markers_include_down$pct.2) < -0.1,]
# write.csv(de_markers_include_down,row.names = F,file = "Plot/6_de_genes_include_down.csv")

pdf("Plot/6_findmarker_UMAP.data.top_logp.pdf",height = 140,width = 15)
FeaturePlot(GW6_total, features = top6$gene,pt.size = 3,slot = "data",ncol = 3,order = T,max.cutoff = "q99",raster = T,raster.dpi = c(1024, 1024),coord.fixed = T)
dev.off()
pdf("Plot/6_findmarker_dot_plot.data.pdf",height = 5,width = 24)
# DotPlot use data matrix as input
# Mouse use RNA assay
DotPlot(GW6_total, features = unique(top6$gene)) + RotatedAxis()
DotPlot(GW6_total,assay = "RNA", features = unique(top6$gene)) + RotatedAxis()
dev.off()

# "Meis1","Col3a1","Pdgfra" # stromal signature
# Progenitor: "EYA1","SIX1","RET","GFRA1"
# Podocyte: "NPHS2","NPHS1","PTPRO"
# Stroma: "NTN1","PTN","SMOC2","NRXN3","IGFBP5"
# ligand/receptor: "EFNA5","EPHA3","EPHA4","EPHA7","WNT5A","FZD4","RYK","CORIN"
# library(paletteer)
for (genei in c("NEB","MYH3")){ # C4 marker
  pdf(paste0("Plot/7_spatial_map.",genei,".SCT_data.pdf"),width = 7,height = 7)
  # 1
  p1 = FeaturePlot(GW6_total, features = genei,slot = "data",order = T,max.cutoff = "q99")
  print(p1)
  # 2
  GW6_total$focus = GW6_total@assays$SCT@data[genei,]
  maxi = sort(GW6_total$focus,decreasing = T)[100]
  mini = 0
  GW6_total$focus = pmin(GW6_total$focus,maxi)
  GW6_total$focus = pmax(GW6_total$focus,mini)
  for (Samplei in Sample.ls) {
    p1 = ggplot(GW6_total@meta.data[GW6_total$orig.ident==Samplei,],aes(x=nB,y=nA,color=focus)) +
      geom_point(shape=19,size=size,show.legend=T) + theme_void() +
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",Samplei)) +
      expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() + 
      scale_color_gradientn(colours=c("#313695FF","#4575B4FF","#74ADD1FF","#ABD9E9FF","#E0F3F8FF","#FFFFBFFF","#FEE090FF","#FDAE61FF","#F46D43FF","#D73027FF","#A50026FF")) + # paletteer_c("grDevices::Blues 3", 30,direction = -1)
      theme(panel.background = element_rect(fill = "black"))
    print(p1)
  }
  dev.off()
}




