library(ggplot2)
library(parallel)
library(dplyr)

samplei = "HEM23"
bin_size = "2kb"
filter_CpG = 5
ct_modif.df = read.csv(paste0("6_analysis/",samplei,"/gene_level/2_all_cell_type.",bin_size,"_bin.modif_level.csv"),row.names = 1)
ct_modif.df$seurat_clusters = gsub("^C(\\d+)_.*","\\1",ct_modif.df$id)
ct_modif.df$windows = gsub("^C\\d+_(.*)","\\1",ct_modif.df$id)

ct_modif.df = ct_modif.df[,c("id","seurat_clusters","windows","V4","V5","V6")]
ct_modif.df = ct_modif.df[order(ct_modif.df$windows),]

ct_modif.df$total = apply(ct_modif.df[,c("V4","V5","V6")],1,sum)
ct_modif.df = ct_modif.df[ct_modif.df$total > 10,]
ct_modif.df$h_level = ct_modif.df$V6/ct_modif.df$total

high_hmc.df = ct_modif.df %>%
  group_by(seurat_clusters) %>%
  filter(h_level >= quantile(h_level, 0.9)) %>%
  ungroup()
colnames(high_hmc.df) = c("id","seurat_clusters","windows","un","mC","hmC","total","h_level")
system(paste0("mkdir 6_analysis/",samplei,"/DMR_hmC/"))
write.csv(high_hmc.df,paste0("6_analysis/",samplei,"/DMR_hmC/3_Total.",bin_size,"_bin.high_hmC_bins.top10percent.csv"),row.names = F)

# tapply(high_hmc.df$h_level, high_hmc.df$seurat_clusters, min)
# table(high_hmc.df$seurat_clusters)
high_hmc.df = read.csv(paste0("6_analysis/",samplei,"/DMR_hmC/3_Total.",bin_size,"_bin.high_hmC_bins.top10percent.csv"))

# Solve chr1_MU069434v1_random_4000_6000
split_region <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  len <- length(parts)
  chr <- paste(parts[1:(len - 2)], collapse = "_")
  start <- parts[len - 1]
  end <- parts[len]
  return(c(chr, start, end))
}
bed.df = as.data.frame(do.call(rbind, lapply(high_hmc.df$windows,split_region)))
bed.df$id = high_hmc.df$id
bed.df$seurat_clusters = high_hmc.df$seurat_clusters
bed.df$strand = "*"
# bed.df$diff_mlevel = round(high_hmc.df$diff_mlevel,2)
write.table(bed.df,paste0("6_analysis/",samplei,"/DMR_hmC/3_Total.",bin_size,"_bin.high_hmC.bed"),row.names = F,col.names = F,quote = F,sep = "\t")

bed.df$total = high_hmc.df$total
bed.df$total_back = high_hmc.df$total_back
write.table(bed.df,paste0("6_analysis/",samplei,"/DMR_hmC/3_Total.",bin_size,"_bin.high_hmC_add.bed"),row.names = F,col.names = F,quote = F,sep = "\t")
# awk '{print > "4_C"$5".2kb_bin.high_hmC.bed"}' 3_Total.2kb_bin.high_hmC.bed

# conda activate Deep_env
# cd /media/yuhao/0EFCCEC5FCCEA5F5/0_spatialMethy/6_sample4_250407_5mC/
# findMotifsGenome.pl 6_analysis/HEM23/DMR_hmC/4_C5.2kb_bin.high_hmC.bed hg38 6_analysis/HEM23/DMR_hmC/MotifOutput_C5/ -size given -mask -useNewBg -p 20 -nomotif

## homer output
result.df = read.csv("6_analysis/HEM23/DMR_hmC/MotifOutput_C5/knownResults.txt",sep = "\t",header = T)
result.df$gene = gsub("^([^(/]+).*", "\\1",result.df$Motif.Name)


