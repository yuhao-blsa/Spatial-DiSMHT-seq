library(ggplot2)
library(parallel)

samplei = "HEM23"
# 1 filter input
CpG.df = read.table("5_Plot/CpG_site.single.HEM23.bed",stringsAsFactors = F)
CpG.df$V7 = gsub("_(\\d+)_","_\\1x",CpG.df$V7) # i701_10_10 -> i701_10x10
CpG.df$V7 = gsub("i702","GW6_1",CpG.df$V7) # i701_10x10 -> testis_1_10x10
CpG.df$V7 = gsub("i703","GW6_2",CpG.df$V7) # i701_10x10 -> testis_1_10x10

# cell_metadata = read.csv("../1a_RNA/Plot_monocle3/1_SCT.gene_metadata.csv",stringsAsFactors = F,row.names = 1)
# cell_metadata = read.csv("../6_merge_GW6_RNA/Plot/1_SCT.cell_metadata.csv",stringsAsFactors = F,row.names = 1)
# cell_metadata$seurat_clusters = cell_metadata$C18
cell_metadata = read.csv("../6_merge_GW6_RNA/Plot_soft_merge/5_meta_data.ai_filter.add_grey.csv",stringsAsFactors = F)
cell_metadata$seurat_clusters = cell_metadata$C18_add_pred
# cell_metadata$seurat_clusters[cell_metadata$seurat_clusters %in% c(13,6,9,15,16,14,5,1,11,2)] = 1
# cell_metadata$seurat_clusters[cell_metadata$seurat_clusters %in% c(10,12,3,8,4,17)] = 2
# cell_metadata$seurat_clusters[cell_metadata$seurat_clusters %in% c(0,7)] = 3
# cell_metadata = cell_metadata[!is.na(cell_metadata$seurat_clusters),]
# CpG.df = CpG.df[CpG.df$V7 %in% cell_metadata$coord,]
system(paste0("mkdir -p 6_analysis/",samplei,"/CpG_level/track"))
system(paste0("mkdir -p 6_analysis/",samplei,"/gene_level"))
system(paste0("mkdir -p 6_analysis/",samplei,"/global_level"))
write.csv(CpG.df,file = paste0("6_analysis/",samplei,"/CpG_level/0_all_pixel.CpG_site.modif_level.csv"))

# get CpG-promoter link
site.info = CpG.df[,1:3]
site.info$id = paste0(site.info$V1,"_",site.info$V2,"_",site.info$V3)
site.info = site.info[!duplicated(site.info),]
# write.table(site.info[,1:3],file = paste0("6_analysis/",samplei,"/CpG_level/CpG_site.bed"),row.names = F,col.names = F,quote = F,sep = "\t")
system(paste0("cd 6_analysis/",samplei,"/CpG_level;bedtools intersect -a CpG_site.bed -b promoter_ud2kb.bed -wa -wb > CpG_site_and_promoter.bed"))
system(paste0("cd 6_analysis/",samplei,"/CpG_level;bedtools intersect -a CpG_site.bed -b genome_2kb_bins.bed -wa -wb > CpG_site_and_2kb_bin.bed"))
system(paste0("cd 6_analysis/",samplei,"/CpG_level;bedtools intersect -a CpG_site.bed -b Transposon_elements/total_focused_famiy.satellite.bed -wa -wb > CpG_site_and_repeat_family.bed"))
system(paste0("cd 6_analysis/",samplei,"/CpG_level;bedtools intersect -a CpG_site.bed -b total_CGI_telomere_centromere.bed -wa -wb > CpG_site_and_other_family.bed"))
system(paste0("cd 6_analysis/",samplei,"/CpG_level;bedtools intersect -a CpG_site.bed -b genebody.bed -wa -wb > CpG_site_and_gene_body.bed"))
promoter_info = read.table(paste0("6_analysis/",samplei,"/CpG_level/CpG_site_and_promoter.bed"))
windows_info2 = read.table(paste0("6_analysis/",samplei,"/CpG_level/CpG_site_and_2kb_bin.bed"))
windows_info2$V8 = paste0(windows_info2$V4,"_",windows_info2$V5,"_",windows_info2$V6)
repeat_info = read.table(paste0("6_analysis/",samplei,"/CpG_level/CpG_site_and_repeat_family.bed"))
gap_info = read.table(paste0("6_analysis/",samplei,"/CpG_level/CpG_site_and_other_family.bed"))
body_info = read.table(paste0("6_analysis/",samplei,"/CpG_level/CpG_site_and_gene_body.bed"))
#     V1        V2        V3   V4        V5        V6 V7      V8
# 1 chr1 149977365 149977366 chr1 149973714 149977714  - Ptgs2os
# 2 chr1 149977365 149977366 chr1 149973782 149977782  +   Ptgs2

# 2 merge function
# input format:
#     V1       V2       V3 V4 V5 V6             V7
# 1 chr1 95100517 95100518  0  1  0 testis_1_10x10
# 2 chr1 95100531 95100532  0  0  1 testis_1_10x10
library(dplyr)
merge_modif_data <- function(CpG.df,cell_metadata=NULL,promoter_info=NULL){
  CpG.df$site = paste0(CpG.df$V1,"_",CpG.df$V2,"_",CpG.df$V3)
  if(!is.null(cell_metadata) & is.null(promoter_info)){
    cell_metadata = cell_metadata[!is.na(cell_metadata$seurat_clusters),]
    CpG.df = CpG.df[CpG.df$V7 %in% cell_metadata$coord,]
    CpG.df$ct = cell_metadata$seurat_clusters[match(CpG.df$V7,cell_metadata$coord)]
    CpG.df$id = paste0("C",CpG.df$ct,"_",CpG.df$site)
  }
  if(is.null(cell_metadata) & !is.null(promoter_info)){
    promoter_info$site = paste0(promoter_info$V1,"_",promoter_info$V2,"_",promoter_info$V3)
    CpG.df$gene = promoter_info$V8[match(CpG.df$site,promoter_info$site)]
    CpG.df = CpG.df[!is.na(CpG.df$gene),]
    CpG.df$id = paste0(CpG.df$V7,"_",CpG.df$gene)
  }
  if(!is.null(cell_metadata) & !is.null(promoter_info)){
    cell_metadata = cell_metadata[!is.na(cell_metadata$seurat_clusters),]
    CpG.df = CpG.df[CpG.df$V7 %in% cell_metadata$coord,]
    CpG.df$ct = cell_metadata$seurat_clusters[match(CpG.df$V7,cell_metadata$coord)]
    promoter_info$site = paste0(promoter_info$V1,"_",promoter_info$V2,"_",promoter_info$V3)
    CpG.df$gene = promoter_info$V8[match(CpG.df$site,promoter_info$site)]
    CpG.df = CpG.df[!is.na(CpG.df$gene),]
    CpG.df$id = paste0("C",CpG.df$ct,"_",CpG.df$gene)
  }
  output = CpG.df[,c("V1","V2","V3","V4","V5","V6","id")]
  ## Slow 1 6min
  system.time(output <- output %>%
                group_by(id) %>%
                summarise(
                  V1 = first(V1),
                  V2 = first(V2),
                  V3 = first(V3),
                  V4 = sum(V4),
                  V5 = sum(V5),
                  V6 = sum(V6)
                ))
  ## Slow 2 +oomin
  # id.ls = sort(unique(output$id))
  # tmp = mclapply(id.ls,function(i){
  #   sum.df = output[which(output$id==i)[1],]
  #   sum.df[1,4:6] = apply(output[which(output$id==i),4:6],2,sum)
  #   return(sum.df)
  # },mc.cores = 64)
  # output = Reduce(rbind,tmp)
  ## Slow 3 83min
  # output = output[order(output$id),]
  # template = 1
  # dup_count = 1
  # for (i in 2:dim(output)[1]) {
  #   if(output$id[i]==output$id[i-1]){
  #     output[template,4:6] = output[template,4:6] + output[i,4:6]
  #     dup_count = dup_count+1
  #     print(i)
  #   }
  #   else{
  #     template = i
  #   }
  # }
  # output = output[!duplicated(output$id),]
  return(output)
}
write.csv(merge_modif_data(CpG.df,cell_metadata,promoter_info=NULL),paste0("6_analysis/",samplei,"/CpG_level/1_all_cell_type.CpG_site.modif_level.csv"))
write.csv(merge_modif_data(CpG.df,cell_metadata=NULL,promoter_info),paste0("6_analysis/",samplei,"/gene_level/0_all_pixel.promoter.modif_level.csv"))
write.csv(merge_modif_data(CpG.df,cell_metadata,promoter_info),paste0("6_analysis/",samplei,"/gene_level/1_all_cell_type.promoter.modif_level.csv"))
write.csv(merge_modif_data(CpG.df,cell_metadata,windows_info2),paste0("6_analysis/",samplei,"/gene_level/2_all_cell_type.2kb_bin.modif_level.csv"))
write.csv(merge_modif_data(CpG.df,cell_metadata,windows_info2),paste0("6_analysis/",samplei,"/gene_level/3_germ_layer.2kb_bin.modif_level.csv"))
write.csv(merge_modif_data(CpG.df,cell_metadata=NULL,windows_info2),paste0("6_analysis/",samplei,"/gene_level/0_all_pixel.2kb_bin.modif_level.csv"))
write.csv(merge_modif_data(CpG.df,cell_metadata=NULL,repeat_info),paste0("6_analysis/",samplei,"/gene_level/0_all_pixel.repeat_family.modif_level.csv"))
write.csv(merge_modif_data(CpG.df,cell_metadata=NULL,gap_info),paste0("6_analysis/",samplei,"/gene_level/0_all_pixel.gap_family.modif_level.csv"))
write.csv(merge_modif_data(CpG.df,cell_metadata=NULL,body_info),paste0("6_analysis/",samplei,"/gene_level/0_all_pixel.genebody.modif_level.csv"))
write.csv(merge_modif_data(CpG.df,cell_metadata,body_info),paste0("6_analysis/",samplei,"/gene_level/1_all_cell_type.genebody.modif_level.csv"))

# 3 bed file
options(digits = 2)
ct_name = c("0"="0_hepatoblast","1"="1_spinalcord1","2"="2_spinalcord3","3"="3_hindlimb1",
            "4"="4_cardiacmyocyte","5"="5_hindbrain","6"="6_forebrain1","7"="7_mesenchyme",
            "8"="8_hindlimb2","9"="9_forebrain2","10"="10_forelimb1","11"="11_spinalcord2",
            "12"="12_forelimb2","13"="13_craniofacial","14"="14_pons","15"="15_midbrain1",
            "16"="16_midbrain2","17"="17_sinusvenosus")
ct_modif.df = read.csv(paste0("6_analysis/",samplei,"/CpG_level/1_all_cell_type.CpG_site.modif_level.csv"),row.names = 1)
ct_modif.df$seurat_clusters = gsub("^C(\\d+)_.*","\\1",ct_modif.df$id)
ct_modif.df[,c("V4","V5","V6")] = ct_modif.df[,c("V4","V5","V6")]/apply(ct_modif.df[,c("V4","V5","V6")],1,sum)
ct_modif.df[,c("V4","V5","V6")] = format(ct_modif.df[,c("V4","V5","V6")], digits = 1)
for (ct_i in names(ct_name)) {
  write.table(ct_modif.df[ct_modif.df$seurat_clusters==ct_i,c("V1","V2","V3","V5")],file = paste0("6_analysis/",samplei,"/CpG_level/track_new/",ct_name[ct_i],"_5mC.bdg"),
              row.names = F,col.names = F,quote = F,sep = "\t")
  write.table(ct_modif.df[ct_modif.df$seurat_clusters==ct_i,c("V1","V2","V3","V6")],file = paste0("6_analysis/",samplei,"/CpG_level/track_new/",ct_name[ct_i],"_5hmC.bdg"),
              row.names = F,col.names = F,quote = F,sep = "\t")
}


