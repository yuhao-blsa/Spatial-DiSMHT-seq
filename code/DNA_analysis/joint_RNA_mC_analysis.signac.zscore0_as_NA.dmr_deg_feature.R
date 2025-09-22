# conda activate monocle3
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
## cover pixel
# hyperDMR_0  hyperDMR_1 hyperDMR_10 hyperDMR_11 hyperDMR_12 hyperDMR_13  hyperDMR_2  hyperDMR_3  hyperDMR_4  hyperDMR_5 
# 3997        2412        1955        1864        1270        1014        2521        2433        2578        2366 
# hyperDMR_6  hyperDMR_7  hyperDMR_8  hyperDMR_9   hypoDMR_0   hypoDMR_1  hypoDMR_10  hypoDMR_11  hypoDMR_12  hypoDMR_13 
# 1935        1942        2070        2159        4263        2532        3178        2167        1793        4845 
# hypoDMR_2   hypoDMR_3   hypoDMR_4   hypoDMR_5   hypoDMR_6   hypoDMR_7   hypoDMR_8   hypoDMR_9 
# 3264        3169        3529        2759        3097        2334        2358        2600
## cover pixel
# downDEG_0  downDEG_1 downDEG_10 downDEG_11 downDEG_12 downDEG_13  downDEG_2  downDEG_3  downDEG_4  downDEG_5  downDEG_6 
# 695       1422       1994        338        480       2625        271        165        148        741       1276 
# downDEG_7  downDEG_8  downDEG_9    upDEG_0    upDEG_1   upDEG_10   upDEG_11   upDEG_12   upDEG_13    upDEG_2    upDEG_3 
# 1365        566        253        422        966       2023        992       1293       1155        648        397 
# upDEG_4    upDEG_5    upDEG_6    upDEG_7    upDEG_8    upDEG_9 
# 1978        562        397        534       1217        680

# 2 load DMR mlevel
pixel_DMR.df = read.csv("6_analysis/HEM23/DMR/5_Total.pixel_level.DMR_hypo_hyper.csv")
pixel_DMR.df = pixel_DMR.df[pixel_DMR.df$coord %in% GW6_total$coord,] # 7501 pixels
pixel_DMR.df = pixel_DMR.df[apply(pixel_DMR.df[,5:7],1,sum)>10,]

# 6 create matrix
total.df = pixel_DMR.df
total.m = matrix(0,nrow = length(unique(total.df$group)),ncol = length(unique(total.df$coord)))
row.names(total.m) = sort(unique(total.df$group))
colnames(total.m) = sort(unique(total.df$coord))
# mclapply not work
total.df$m_zscore = unsplit(
  tapply(total.df$m_level,total.df$group,function(x) as.vector(scale(x))),
  total.df$group
)
tmp = lapply(1:dim(total.df)[1],function(i){
  total.m[total.df$group[i],total.df$coord[i]] <<- total.df$m_zscore[i]
})

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
GW6_total <- RunUMAP(GW6_total, reduction = 'lsi', dims = 1:20, reduction.name = 'umap.meth')
pdf("6_analysis/HEM23/DMR/6_DMR_methy_level.signc_minCover11_UMAP20_2.pdf",width = 10,height = 10)
DimPlot(GW6_total, reduction = 'umap.meth', label = TRUE,pt.size = 1.5) + NoLegend() + ggtitle("CpG UMAP") + scale_color_manual(values = color.ls5a2)+ coord_fixed()
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




