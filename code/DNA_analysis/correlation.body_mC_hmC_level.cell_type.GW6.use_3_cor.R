library(ggplot2)
library(parallel)
library(paletteer)
library(ComplexHeatmap)
library(dplyr)

# 1 input
ct_gene_mC.df = read.csv("6_analysis/HEM23/gene_level/1_all_cell_type.promoter.modif_level.csv",row.names = 1)
ct_gene_mC.df$cluster = gsub("C(\\d+)_.*","C\\1",ct_gene_mC.df$id) # c(13,6,9,15,16,14,5,1,11,2,7,10,12,4,17,0,3,8)
ct_gene_mC.df$gene = gsub("C\\d+_(.*)","\\1",ct_gene_mC.df$id)
ct_gene_mC.df$merge = c("C13"="craniofacial","C6"="forebrain","C9"="forebrain","C15"="midbrain","C16"="midbrain",
                        "C14"="midbrain","C5"="hindbrain","C1"="spinalcord","C11"="spinalcord","C2"="spinalcord",
                        "C7"="mesenchyme","C10"="forelimb","C12"="forelimb","C4"="heart","C17"="heart",
                        "C0"="hepatoblast","C3"="hindlimb","C8"="hindlimb")[ct_gene_mC.df$cluster]
merge.df = ct_gene_mC.df %>%
  group_by(gene,merge) %>%
  summarise(
    V4 = sum(V4, na.rm = TRUE),
    V5 = sum(V5, na.rm = TRUE),
    V6 = sum(V6, na.rm = TRUE),
    .groups = "drop"
  )

ct_gene_mC.cov.m = matrix(-1,nrow = length(unique(merge.df$gene)),ncol = length(unique(merge.df$merge)),dimnames = list(sort(unique(merge.df$gene)),
                                                                                                                        c("craniofacial","forebrain","midbrain","hindbrain","spinalcord","mesenchyme","forelimb","heart","hepatoblast","hindlimb")))
ct_gene_mC.level.m = matrix(-1,nrow = length(unique(merge.df$gene)),ncol = length(unique(merge.df$merge)),dimnames = list(sort(unique(merge.df$gene)),
                                                                                                                          c("craniofacial","forebrain","midbrain","hindbrain","spinalcord","mesenchyme","forelimb","heart","hepatoblast","hindlimb")))
ct_gene_hC.level.m = matrix(-1,nrow = length(unique(merge.df$gene)),ncol = length(unique(merge.df$merge)),dimnames = list(sort(unique(merge.df$gene)),
                                                                                                                          c("craniofacial","forebrain","midbrain","hindbrain","spinalcord","mesenchyme","forelimb","heart","hepatoblast","hindlimb")))

for (i in 1:dim(merge.df)[1]) {
  ct_gene_mC.cov.m[merge.df$gene[i],merge.df$merge[i]] = sum(merge.df[i,3:5])
  ct_gene_mC.level.m[merge.df$gene[i],merge.df$merge[i]] = as.numeric(merge.df[i,4]/sum(merge.df[i,3:5]))
  ct_gene_hC.level.m[merge.df$gene[i],merge.df$merge[i]] = as.numeric(merge.df[i,5]/sum(merge.df[i,3:5]))
}
table(apply(ct_gene_mC.cov.m,1,function(x) sum(x>0)))
# 1    2    3    4    5    6    7    8    9   10 
# 953 1718 3037 4895 6716 7841 7126 4719 2049  654 
# sort(apply(ct_gene_mC.cov.m,1,function(x) sum(x>0)),decreasing = T)[1:7422] # cover >= 8 cell type

ct_gene_mC.level.m = ct_gene_mC.level.m[apply(ct_gene_mC.cov.m,1,function(x) sum(x>0))>=8,]
ct_gene_hC.level.m = ct_gene_hC.level.m[apply(ct_gene_mC.cov.m,1,function(x) sum(x>0))>=8,]
ct_gene_mC.cov.m = ct_gene_mC.cov.m[apply(ct_gene_mC.cov.m,1,function(x) sum(x>0))>=8,]

save(list = c("ct_gene_mC.level.m","ct_gene_hC.level.m","ct_gene_mC.cov.m"),file = "6_analysis/HEM23/gene_mC_hmC/1_all_cell_type.genebody.modif_level.filter_coverage.Rdata")
load("6_analysis/HEM23/gene_mC_hmC/1_all_cell_type.genebody.modif_level.filter_coverage.Rdata")

## correlation
load("../6_merge_GW6_RNA/Plot/12_agg_ct.quantile.merged_ct.Rdata")
agg_ct.quantile = agg_ct.quantile[,colnames(ct_gene_mC.level.m)]
gene.ls = intersect(row.names(agg_ct.quantile),row.names(ct_gene_mC.level.m)) # 4057
agg_ct.quantile = agg_ct.quantile[gene.ls,]
agg_ct.quantile = t(scale(t(agg_ct.quantile)))
ct_gene_mC.level.m = ct_gene_mC.level.m[gene.ls,]
ct_gene_hC.level.m = ct_gene_hC.level.m[gene.ls,]
ct_gene_mC.cov.m = log1p(ct_gene_mC.cov.m[gene.ls,]+1)

## correlation
cor.ls1 = unlist(lapply(1:dim(agg_ct.quantile)[1],function(i){
  a=agg_ct.quantile[i,]
  b=ct_gene_mC.level.m[i,]
  a=a[b > -1]
  b=b[b > -1]
  return(cor(a,b))
}))
names(cor.ls1) = row.names(agg_ct.quantile)
cor.ls2 = unlist(lapply(1:dim(agg_ct.quantile)[1],function(i){
  a=ct_gene_hC.level.m[i,]
  b=ct_gene_mC.level.m[i,]
  a=a[b > -1]
  b=b[b > -1]
  return(cor(a,b))
}))
names(cor.ls2) = row.names(agg_ct.quantile)
cor.ls3 = unlist(lapply(1:dim(agg_ct.quantile)[1],function(i){
  a=agg_ct.quantile[i,]
  b=ct_gene_hC.level.m[i,]
  a=a[b > -1]
  b=b[b > -1]
  return(cor(a,b))
}))
names(cor.ls3) = row.names(agg_ct.quantile)
cor.ls1 = cor.ls1[!(is.na(cor.ls1) | is.na(cor.ls2) | is.na(cor.ls3))]
cor.ls2 = cor.ls2[names(cor.ls1)] # 3966
cor.ls3 = cor.ls3[names(cor.ls1)] # 3966

# statistic
df <- tibble(cor.ls1, cor.ls2, cor.ls3)
df2 <- df %>%
  mutate(
    cat_x = case_when(
      cor.ls1 > 0.4 ~ ">0.5",
      cor.ls1 < -0.4 ~ "<-0.5",
      TRUE ~ "other"),
    cat_y = case_when(
      cor.ls2> 0.4 ~ ">0.5",
      cor.ls2 < -0.4 ~ "<-0.5",
      TRUE ~ "other"),
    cat_z = case_when(
      cor.ls3> 0.4 ~ ">0.5",
      cor.ls3 < -0.4 ~ "<-0.5",
      TRUE ~ "other"))
df2 %>% count(cat_x, cat_y, cat_z)
# cat_x cat_y cat_z     n
# 1 <-0.5 <-0.5 <-0.5     1
# 2 <-0.5 <-0.5 >0.5     42
# 3 <-0.5 <-0.5 other    44
# 4 <-0.5 >0.5  <-0.5    55
# 5 <-0.5 >0.5  >0.5      1
# 6 <-0.5 >0.5  other   103
# 7 <-0.5 other <-0.5    61
# 8 <-0.5 other >0.5     68
# 9 <-0.5 other other   304
# 10 >0.5  <-0.5 <-0.5    24
# 11 >0.5  <-0.5 other    37
# 12 >0.5  >0.5  <-0.5     1
# 13 >0.5  >0.5  >0.5     59
# 14 >0.5  >0.5  other    52
# 15 >0.5  other <-0.5    42

agg_ct.quantile = agg_ct.quantile[c(names(which(cor.ls1 < -0.4 & cor.ls2 < -0.4 & cor.ls3 > 0.4)), # 42
                                    names(which(cor.ls1 < -0.4 & cor.ls2 > 0.4 & cor.ls3 < -0.4)), # 55
                                    names(which(cor.ls1 > 0.4 & cor.ls2 < -0.4 & cor.ls3 < -0.4)), # 24
                                    names(which(cor.ls1 > 0.4 & cor.ls2 > 0.4 & cor.ls3 > 0.4))),] # 59
ct_gene_mC.level.m = ct_gene_mC.level.m[row.names(agg_ct.quantile),]
ct_gene_hC.level.m = ct_gene_hC.level.m[row.names(agg_ct.quantile),]
ct_gene_mC.cov.m = ct_gene_mC.cov.m[row.names(agg_ct.quantile),]

ct_gene_mC.level.m[ct_gene_mC.level.m<0] = NA
ct_gene_hC.level.m[ct_gene_hC.level.m<0] = NA

# reorder
high_expr_bin.ls = unlist(lapply(1:nrow(agg_ct.quantile),function(i){
  tmp = agg_ct.quantile[i,]
  return(mean(which(tmp >= sort(tmp,decreasing = T)[1])))
}))
high_expr_distance.ls = unlist(lapply(1:nrow(agg_ct.quantile),function(i){
  tmp = agg_ct.quantile[i,]
  range = max(which(tmp >= sort(tmp,decreasing = T)[2])) - min(which(tmp >= sort(tmp,decreasing = T)[2]))
  return(range)
}))
names(high_expr_bin.ls) = row.names(agg_ct.quantile)
high_expr_bin.ls = high_expr_bin.ls[high_expr_distance.ls<5]
high_expr_bin.ls = sort(high_expr_bin.ls)
final_order = c(intersect(names(high_expr_bin.ls),names(which(cor.ls1 > 0.4 & cor.ls2 > 0.4 & cor.ls3 > 0.4))), # 47
                intersect(names(high_expr_bin.ls),names(which(cor.ls1 > 0.4 & cor.ls2 < -0.4 & cor.ls3 < -0.4))), # 19
                intersect(names(high_expr_bin.ls),names(which(cor.ls1 < -0.4 & cor.ls2 > 0.4 & cor.ls3 < -0.4))), # 39
                intersect(names(high_expr_bin.ls),names(which(cor.ls1 < -0.4 & cor.ls2 < -0.4 & cor.ls3 > 0.4)))) # 29
                
                

agg_ct.quantile = agg_ct.quantile[final_order,]
ct_gene_mC.level.m = ct_gene_mC.level.m[final_order,]
ct_gene_hC.level.m = ct_gene_hC.level.m[final_order,]
ct_gene_mC.cov.m = ct_gene_mC.cov.m[final_order,]

ct_gene_mC.level.norm = t(scale(t(ct_gene_mC.level.m),scale = F))
ct_gene_hC.level.norm = t(scale(t(ct_gene_hC.level.m),scale = F))
ct_gene_mC.level.scale = t(scale(t(ct_gene_mC.level.m),scale = T))
ct_gene_hC.level.scale = t(scale(t(ct_gene_hC.level.m),scale = T))

library(ComplexHeatmap)
# 构建热图列表
split_vec <- c(rep("group1", 47), rep("group2", 19), rep("group3", 39), rep("group4", 29))
ht_list <- Heatmap(agg_ct.quantile,cluster_columns = F,cluster_rows = F,col=paletteer_c("ggthemes::Red-Blue-White Diverging",10,direction = -1),
                   row_names_gp = gpar(fontsize = 2),row_split  = split_vec,column_title = "RNA") +
  Heatmap(ct_gene_mC.level.m,cluster_columns = F,cluster_rows = F,col=paletteer_c("grDevices::Viridis", 10),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5mC") +
  Heatmap(ct_gene_hC.level.m,cluster_columns = F,cluster_rows = F,col=paletteer_c("grDevices::Viridis", 10),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5hmC")
pdf("6_analysis/HEM23/gene_mC_hmC/7_RNA_and_5mC.heatmap.top1_cell_type_order.pdf",height = 5,width = 7)
draw(ht_list,heatmap_legend_side = "right",annotation_legend_side = "right")
dev.off()

library(circlize)
pdf("6_analysis/HEM23/gene_mC_hmC/7_RNA_and_5mC.heatmap.top1_cell_type_order.norm.pdf",height = 5,width = 7)
ht_list2 <- Heatmap(agg_ct.quantile,cluster_columns = F,cluster_rows = F,col=paletteer_c("ggthemes::Red-Blue-White Diverging",10,direction = -1),
                    row_names_gp = gpar(fontsize = 2),row_split  = split_vec,column_title = "RNA") +
  Heatmap(ct_gene_mC.level.norm,cluster_columns = F,cluster_rows = F,col=paletteer_c("grDevices::Viridis", 10),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5mC") +
  Heatmap(ct_gene_hC.level.norm,cluster_columns = F,cluster_rows = F,col=paletteer_c("grDevices::Viridis", 10),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5hmC")
draw(ht_list2,heatmap_legend_side = "right",annotation_legend_side = "right")
ht_list2 <- Heatmap(agg_ct.quantile,cluster_columns = F,cluster_rows = F,col=paletteer_c("ggthemes::Red-Blue-White Diverging",10,direction = -1),
                    row_names_gp = gpar(fontsize = 2),row_split  = split_vec,column_title = "RNA") +
  Heatmap(ct_gene_mC.level.norm,cluster_columns = F,cluster_rows = F,
          col=colorRamp2(breaks = seq(-0.2,0.2,length.out=10),colors = paletteer_c("grDevices::Viridis", 10)),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5mC") +
  Heatmap(ct_gene_hC.level.norm,cluster_columns = F,cluster_rows = F,
          col=colorRamp2(breaks = seq(-0.05,0.05,length.out=10),colors = paletteer_c("grDevices::Viridis", 10)),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5hmC")
draw(ht_list2,heatmap_legend_side = "right",annotation_legend_side = "right")
ht_list2 <- Heatmap(agg_ct.quantile,cluster_columns = F,cluster_rows = F,col=paletteer_c("ggthemes::Red-Blue-White Diverging",10,direction = -1),
                    row_names_gp = gpar(fontsize = 2),row_split  = split_vec,column_title = "RNA") +
  Heatmap(ct_gene_mC.level.norm,cluster_columns = F,cluster_rows = F,
          col=colorRamp2(breaks = seq(-0.15,0.15,length.out=10),colors = paletteer_c("grDevices::Viridis", 10)),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5mC") +
  Heatmap(ct_gene_hC.level.norm,cluster_columns = F,cluster_rows = F,
          col=colorRamp2(breaks = seq(-0.04,0.04,length.out=10),colors = paletteer_c("grDevices::Viridis", 10)),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5hmC")
draw(ht_list2,heatmap_legend_side = "right",annotation_legend_side = "right")
ht_list2 <- Heatmap(agg_ct.quantile,cluster_columns = F,cluster_rows = F,col=paletteer_c("ggthemes::Red-Blue-White Diverging",10,direction = -1),
                    row_names_gp = gpar(fontsize = 2),row_split  = split_vec,column_title = "RNA") +
  Heatmap(ct_gene_mC.level.norm,cluster_columns = F,cluster_rows = F,
          col=colorRamp2(breaks = seq(-0.1,0.1,length.out=10),colors = paletteer_c("grDevices::Viridis", 10)),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5mC") +
  Heatmap(ct_gene_hC.level.norm,cluster_columns = F,cluster_rows = F,
          col=colorRamp2(breaks = seq(-0.02,0.02,length.out=10),colors = paletteer_c("grDevices::Viridis", 10)),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5hmC")
draw(ht_list2,heatmap_legend_side = "right",annotation_legend_side = "right")
dev.off()

pdf("6_analysis/HEM23/gene_mC_hmC/7_RNA_and_5mC.heatmap.top1_cell_type_order.scale.pdf",height = 5,width = 7)
ht_list2 <- Heatmap(agg_ct.quantile,cluster_columns = F,cluster_rows = F,col=paletteer_c("ggthemes::Red-Blue-White Diverging",10,direction = -1),
                    row_names_gp = gpar(fontsize = 2),row_split  = split_vec,column_title = "RNA") +
  Heatmap(ct_gene_mC.level.scale,cluster_columns = F,cluster_rows = F,col=paletteer_c("grDevices::Viridis", 10),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5mC") +
  Heatmap(ct_gene_hC.level.scale,cluster_columns = F,cluster_rows = F,col=paletteer_c("grDevices::Viridis", 10),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5hmC")
draw(ht_list2,heatmap_legend_side = "right",annotation_legend_side = "right")
ht_list2 <- Heatmap(agg_ct.quantile,cluster_columns = F,cluster_rows = F,col=paletteer_c("ggthemes::Red-Blue-White Diverging",10,direction = -1),
                    row_names_gp = gpar(fontsize = 2),row_split  = split_vec,column_title = "RNA") +
  Heatmap(ct_gene_mC.level.scale,cluster_columns = F,cluster_rows = F,
          col=colorRamp2(breaks = seq(-2,2,length.out=10),colors = paletteer_c("grDevices::Viridis", 10)),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5mC") +
  Heatmap(ct_gene_hC.level.scale,cluster_columns = F,cluster_rows = F,
          col=colorRamp2(breaks = seq(-2,2,length.out=10),colors = paletteer_c("grDevices::Viridis", 10)),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5hmC")
draw(ht_list2,heatmap_legend_side = "right",annotation_legend_side = "right")
ht_list2 <- Heatmap(agg_ct.quantile,cluster_columns = F,cluster_rows = F,col=paletteer_c("ggthemes::Red-Blue-White Diverging",10,direction = -1),
                    row_names_gp = gpar(fontsize = 2),row_split  = split_vec,column_title = "RNA") +
  Heatmap(ct_gene_mC.level.scale,cluster_columns = F,cluster_rows = F,
          col=colorRamp2(breaks = seq(-1,1,length.out=10),colors = paletteer_c("grDevices::Viridis", 10)),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5mC") +
  Heatmap(ct_gene_hC.level.scale,cluster_columns = F,cluster_rows = F,
          col=colorRamp2(breaks = seq(-1,1,length.out=10),colors = paletteer_c("grDevices::Viridis", 10)),row_names_gp = gpar(fontsize = 2),
          row_split  = split_vec,column_title = "5hmC")
draw(ht_list2,heatmap_legend_side = "right",annotation_legend_side = "right")
dev.off()

# Heatmap(agg_ct.quantile,cluster_columns = F,cluster_rows = F,col=paletteer_c("ggthemes::Red-Blue-White Diverging",10,direction = -1),row_names_gp = gpar(fontsize = 7))
# Heatmap(ct_gene_mC.cov.m,cluster_columns = F,cluster_rows = F,col=paletteer_c("ggthemes::Red-Blue-White Diverging",10,direction = -1),row_names_gp = gpar(fontsize = 7))
# Heatmap(ct_gene_mC.level.m,cluster_columns = F,cluster_rows = F,col=paletteer_c("grDevices::Viridis", 10),row_names_gp = gpar(fontsize = 7))
# Heatmap(t(scale(t(ct_gene_mC.level.m),scale = F)),cluster_columns = F,cluster_rows = F,col=paletteer_c("grDevices::Viridis", 10),row_names_gp = gpar(fontsize = 7))
# tmp=t(scale(t(ct_gene_mC.level.m),scale = F))
# tmp[tmp > 0.4] = 0.4
# tmp[tmp < -0.4] = -0.4
# Heatmap(tmp,cluster_columns = F,cluster_rows = F,col=paletteer_c("grDevices::Viridis", 10),row_names_gp = gpar(fontsize = 7))
# tmp[tmp > 0.2] = 0.2
# tmp[tmp < -0.2] = -0.2
# Heatmap(tmp,cluster_columns = F,cluster_rows = F,col=paletteer_c("grDevices::Viridis", 10),row_names_gp = gpar(fontsize = 7))

