# Origin: 0_spatialMethy/4_merge_E9.5_RNA/4_joint_RNA_mC_analysis.signac.zscore0_as_NA.dmr_deg_feature.R

library(Seurat)
library(ggplot2)
library(Signac)
library(stringr)
library(dplyr)
library(parallel)

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

# 1 load RNA
GW6_total = readRDS("../6_merge_GW6_RNA/Plot/1_SCT_transform&Merged_data.add_C17.filter.RDS")
GW6_total = GW6_total[,GW6_total$isfinal==1]
Sample.ls = c("GW6_1","GW6_2")

# 2 load DMR mlevel
pixel_DMR.df = read.csv("6_analysis/HEM23/DMR_hmC/5_Total.pixel_level.DMR_hypo_hyper.csv")
pixel_DMR.df = pixel_DMR.df[pixel_DMR.df$coord %in% GW6_total$coord,]
pixel_DMR.df = pixel_DMR.df[apply(pixel_DMR.df[,5:7],1,sum)>10,] # 3

# 6 create matrix
total.df = pixel_DMR.df
total.m = matrix(0,nrow = length(unique(total.df$group)),ncol = length(unique(total.df$coord)))
row.names(total.m) = sort(unique(total.df$group))
colnames(total.m) = sort(unique(total.df$coord))
# mclapply not work
total.df$h_zscore = unsplit(
  tapply(total.df$h_level,total.df$group,function(x) as.vector(scale(x))),
  total.df$group
)
tmp = lapply(1:dim(total.df)[1],function(i){
  total.m[total.df$group[i],total.df$coord[i]] <<- total.df$h_zscore[i]
})

## filter! ##
cell.ls = names(table(total.df$coord))[table(total.df$coord) >= 10] # 3
total.m = total.m[,cell.ls]

# 6 create ATAC assay and add it to the object
# row.names(total.m) = gsub("_",".",row.names(total.m))
# row.names(total.m) = gsub("\\.([^.]+)\\.([^.]+)$", "_\\1_\\2",row.names(total.m)) # chrUn_JH584304v1_0_114452 -> chrUn.JH584304v1_0_114452
row.names(total.m) = gsub("_","",row.names(total.m))
row.names(total.m) = paste0(row.names(total.m),"_100_200")
GW6_total[["ATAC"]] <- CreateChromatinAssay(
  counts = total.m,
  sep = c("_", "_"),
  fragments = NULL,
  annotation = NULL
)

# 7 UMAP visualization
DefaultAssay(GW6_total) <- "ATAC"
# GW6_total <- FindTopFeatures(GW6_total, min.cutoff = 5) # not suitable for methylation
VariableFeatures(GW6_total) = row.names(total.m)
# GW6_total <- RunTFIDF(GW6_total, verbose = T) # normalize ATAC$data # not suitable for methylation # 设计用于“词频统计”，适合离散型词;TF-IDF不能体现"值的大小和方向性"
# GW6_total@assays$ATAC$data
GW6_total <- RunSVD(GW6_total,n=length(VariableFeatures(GW6_total))-1)
GW6_total <- RunUMAP(GW6_total, reduction = 'lsi', dims = 1:17, reduction.name = 'umap.meth')
pdf("6_analysis/HEM23/DMR_hmC/6_DMR_methy_level.signc_minCover10_UMAP17.pdf",width = 7,height = 7)
DimPlot(GW6_total, reduction = 'umap.meth', label = TRUE,pt.size = 0.3) + NoLegend() + ggtitle("CpG UMAP") + scale_color_manual(values = color.ls5a2)+ coord_fixed()
dev.off()





# build a joint neighbor graph using both assays
GW6_total <- FindMultiModalNeighbors(
  object = GW6_total,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)
# build a joint UMAP visualization
GW6_total <- RunUMAP(
  object = GW6_total,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)
DimPlot(GW6_total, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend() + scale_color_manual(values = color.ls5a)




