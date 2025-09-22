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
global_modif.df <- ct_modif.df %>%
  group_by(windows) %>%
  summarise(
    V1 = first(windows),
    V2 = sum(V4),
    V3 = sum(V5),
    V4 = sum(V6)
  )

system(paste0("mkdir -p 6_analysis/",samplei,"/DMR"))
for (cti in sort(unique(ct_modif.df$seurat_clusters))) {
  print(cti)
  # 1)
  cti_result.df = ct_modif.df[ct_modif.df$seurat_clusters==cti,c("seurat_clusters","windows","V4","V5","V6")] # mark1
  cti_result.df$total = apply(cti_result.df[,c("V4","V5","V6")],1,sum)
  colnames(cti_result.df) = c("seurat_clusters","windows","un","mC","hmC","total")
  # 2)
  cti_result.df = cbind(cti_result.df, global_modif.df[match(cti_result.df$windows,global_modif.df$windows),c("V2","V3","V4")])
  cti_result.df$total_back = apply(cti_result.df[,c("V2","V3","V4")],1,sum)
  colnames(cti_result.df) = c("seurat_clusters","windows","un","mC","hmC","total","un_back","mC_back","hmC_back","total_back")
  # 3)
  cti_result.df$un_back = cti_result.df$un_back - cti_result.df$un
  cti_result.df$mC_back = cti_result.df$mC_back - cti_result.df$mC
  cti_result.df$hmC_back = cti_result.df$hmC_back - cti_result.df$hmC
  cti_result.df$total_back = cti_result.df$total_back - cti_result.df$total
  # 4)
  cti_result.df = cti_result.df[cti_result.df$total>=filter_CpG & cti_result.df$total_back>=filter_CpG,] # mark2
  base_bin = dim(cti_result.df)[1]
  cti_result.df$diff_mlevel = cti_result.df$mC/cti_result.df$total - cti_result.df$mC_back/cti_result.df$total_back
  cti_result.df = cti_result.df[abs(cti_result.df$diff_mlevel)>=0.2,] # mark3
  # 5)
  cti_result.df$unmC <- cti_result.df$total - cti_result.df$mC
  cti_result.df$unmC_back <- cti_result.df$total_back - cti_result.df$mC_back
  cti_result.df$p_value <- unlist(mclapply(1:dim(cti_result.df)[1], function(i) {
    mat1 <- matrix(c(cti_result.df[i,"mC"], cti_result.df[i,"unmC"], cti_result.df[i,"mC_back"], cti_result.df[i,"unmC_back"]),
                   nrow = 2, byrow = TRUE)
    fisher.test(mat1)$p.value
  },mc.cores = 60))
  cti_result.df$adjustp <- p.adjust(cti_result.df$p_value, method = "BH")
  cti_result.df = cti_result.df[order(cti_result.df$adjustp),]
  write.csv(cti_result.df,paste0("6_analysis/",samplei,"/DMR/1_C",cti,".",bin_size,"_bin.DMR.csv"),row.names = F)
}

# table(cti_result.df$adjustp < 0.1)
# table(cti_result.df$adjustp < 0.05)
# table(cti_result.df$adjustp < 0.01)
# 
# round(sum(cti_result.df$adjustp < 0.1)/base_bin*100,2)
# round(sum(cti_result.df$adjustp < 0.05)/base_bin*100,2)
# round(sum(cti_result.df$adjustp < 0.01)/base_bin*100,2)

## summary
total_result.df = data.frame()
for (cti in sort(unique(ct_modif.df$seurat_clusters))) {
  cti_result.df = read.csv(paste0("6_analysis/",samplei,"/DMR/1_C",cti,".",bin_size,"_bin.DMR.csv"))
  cti_result.df = cti_result.df[cti_result.df$adjustp < 0.05,]
  total_result.df = rbind(total_result.df,cti_result.df)
}
write.csv(total_result.df,paste0("6_analysis/",samplei,"/DMR/2_Total.",bin_size,"_bin.DMR.csv"),row.names = F)
total_result.df = read.csv(paste0("6_analysis/",samplei,"/DMR/2_Total.",bin_size,"_bin.DMR.csv"))

# Solve chr1_MU069434v1_random_4000_6000
split_region <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  len <- length(parts)
  chr <- paste(parts[1:(len - 2)], collapse = "_")
  start <- parts[len - 1]
  end <- parts[len]
  return(c(chr, start, end))
}
bed.df = as.data.frame(do.call(rbind, lapply(total_result.df$windows,split_region)))
bed.df$seurat_clusters = total_result.df$seurat_clusters
bed.df$diff_mlevel = round(total_result.df$diff_mlevel,2)
write.table(bed.df,paste0("6_analysis/",samplei,"/DMR/2_Total.",bin_size,"_bin.DMR.bed"),row.names = F,col.names = F,quote = F,sep = "\t")

bed.df$total = total_result.df$total
bed.df$total_back = total_result.df$total_back
write.table(bed.df,paste0("6_analysis/",samplei,"/DMR/2_Total.",bin_size,"_bin.DMR_add.bed"),row.names = F,col.names = F,quote = F,sep = "\t")



# awk '{print > "3_C"$4".2kb_bin.DMR.bed"}' 2_Total.2kb_bin.DMR.bed
# awk -v OFS="\t" '{print $1, $2, $3 > "3_C"$4".2kb_bin.DMR.col3.bed"}' 2_Total.2kb_bin.DMR.bed
# awk -v OFS="\t" '$5 > 0{print $1, $2, $3 > "4_C"$4".2kb_bin.DMR.hyper.bed"} $5 < 0{print $1, $2, $3 > "4_C"$4".2kb_bin.DMR.hypo.bed"}' 2_Total.2kb_bin.DMR.bed


