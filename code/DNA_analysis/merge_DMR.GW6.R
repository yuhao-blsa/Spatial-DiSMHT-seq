library(ggplot2)
library(parallel)

samplei = "HEM23"
bin_size = "2kb"
sample.ls = c("GW6_1","GW6_2")
# 1 filter input
pixel_modif.df = read.csv(paste0("6_analysis/",samplei,"/gene_level/0_all_pixel.",bin_size,"_bin.modif_level.csv"),row.names = 1)
pixel_modif.df$coord = gsub("^(.*\\d+x\\d+)_.*","\\1",pixel_modif.df$id)
pixel_modif.df$windows = gsub("^.*\\d+x\\d+_(.*)","\\1",pixel_modif.df$id)
pixel_modif.df$V2 = gsub("^.*_(\\d+)_\\d+$","\\1",pixel_modif.df$windows)
pixel_modif.df$V3 = gsub("^.*_\\d+_(\\d+)$","\\1",pixel_modif.df$windows)
# cell_metadata = read.csv("../6_merge_GW6_RNA/Plot/1_SCT.cell_metadata.csv",stringsAsFactors = F,row.names = 1)
# cell_metadata$seurat_clusters = cell_metadata$C18
cell_metadata = read.csv("../6_merge_GW6_RNA/Plot_soft_merge/5_meta_data.ai_filter.add_grey.csv",stringsAsFactors = F)
cell_metadata$seurat_clusters = cell_metadata$C18_add_pred
DMR.df = read.csv(paste0("6_analysis/",samplei,"/DMR/2_Total.2kb_bin.DMR.csv"))
# Solve chr1_MU069434v1_random_4000_6000
split_region <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  len <- length(parts)
  chr <- paste(parts[1:(len - 2)], collapse = "_")
  start <- parts[len - 1]
  end <- parts[len]
  return(c(chr, start, end))
}
tmp = as.data.frame(do.call(rbind, lapply(DMR.df$windows,split_region)))
DMR.df = cbind(DMR.df,tmp)

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

# 4 DMR pixel average methylation level
DMR_mlevel.total = data.frame()
for (cti in sort(unique(DMR.df$seurat_clusters))) {
  DMR_hypoi = DMR.df[DMR.df$seurat_clusters==cti & DMR.df$diff_mlevel<0,c("windows","V1","V2","V3")]
  DMR_hyeri = DMR.df[DMR.df$seurat_clusters==cti & DMR.df$diff_mlevel>0,c("windows","V1","V2","V3")]
  DMR_hypoi$V8 = paste0("C",cti,"_hypo")
  DMR_hyeri$V8 = paste0("C",cti,"_hyer")
  inputi = pixel_modif.df[pixel_modif.df$windows %in% DMR_hypoi$windows, c("V1","V2","V3","V4","V5","V6","coord")]
  colnames(inputi) = c("V1","V2","V3","V4","V5","V6","V7")
  DMR_mlevel.hypoi = merge_modif_data(inputi,cell_metadata=NULL,DMR_hypoi)
  inputi = pixel_modif.df[pixel_modif.df$windows %in% DMR_hyeri$windows, c("V1","V2","V3","V4","V5","V6","coord")]
  colnames(inputi) = c("V1","V2","V3","V4","V5","V6","V7")
  DMR_mlevel.hyeri = merge_modif_data(inputi,cell_metadata=NULL,DMR_hyeri)
  
  DMR_mlevel.total = rbind(DMR_mlevel.total,DMR_mlevel.hypoi)
  DMR_mlevel.total = rbind(DMR_mlevel.total,DMR_mlevel.hyeri)
  }

DMR_mlevel.total$coord = gsub("_C\\d+_(hypo|hyer)$","",DMR_mlevel.total$id)
DMR_mlevel.total = DMR_mlevel.total[DMR_mlevel.total$coord %in% cell_metadata$coord[cell_metadata$isfinal==1],] # filter pixel
DMR_mlevel.total$V1 = cell_metadata$orig.ident[match(DMR_mlevel.total$coord,cell_metadata$coord)]
DMR_mlevel.total$V2 = cell_metadata$nA[match(DMR_mlevel.total$coord,cell_metadata$coord)]
DMR_mlevel.total$V3 = cell_metadata$nB[match(DMR_mlevel.total$coord,cell_metadata$coord)]
DMR_mlevel.total$m_level = DMR_mlevel.total$V5/apply(DMR_mlevel.total[,c("V4","V5","V6")],1,sum)
DMR_mlevel.total$h_level = DMR_mlevel.total$V6/apply(DMR_mlevel.total[,c("V4","V5","V6")],1,sum)
DMR_mlevel.total$group = gsub("^.*_\\d+x\\d+_(.*)$","\\1",DMR_mlevel.total$id)
colnames(DMR_mlevel.total) = c("id","orig.ident","nA","nB","un","mC","hmC","coord","m_level","h_level","group")
write.csv(DMR_mlevel.total,paste0("6_analysis/",samplei,"/DMR/5_Total.pixel_level.DMR_hypo_hyper.csv"),row.names = F)
# DMR_mlevel.total = read.csv(paste0("6_analysis/",samplei,"/DMR/5_Total.pixel_level.DMR_hypo_hyper.csv"))

system(paste0("mkdir -p 6_analysis/",samplei,"/DMR_plot"))
library(paletteer)
for (groupi in sort(unique(DMR_mlevel.total$group))) {
  DMR_mlevel.dfi = DMR_mlevel.total[DMR_mlevel.total$group==groupi,]
  DMR_mlevel.dfi$m_level = pmin(pmax(DMR_mlevel.dfi$m_level,quantile(DMR_mlevel.dfi$m_level,0.1)),quantile(DMR_mlevel.dfi$m_level,0.9))
  DMR_mlevel.dfi$h_level = pmin(pmax(DMR_mlevel.dfi$h_level,quantile(DMR_mlevel.dfi$h_level,0.1)),quantile(DMR_mlevel.dfi$h_level,0.9))
  pdf(paste0("6_analysis/",samplei,"/DMR_plot/1.",groupi,".5mC_level_0.1_0.9.pdf"))
  for (si in sample.ls) {
    p1 = ggplot(DMR_mlevel.dfi[DMR_mlevel.dfi$orig.ident==si,],aes(x=nB,y=nA,fill=m_level)) +
      geom_tile(show.legend=T,width=0.9,height=0.9) + theme_void() +
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",si)) +
      expand_limits(y=c(0,96),x=c(0,96)) + coord_fixed() +
      scale_fill_gradientn(colours=c(paletteer_c("grDevices::Viridis", 30,direction = -1)))
    print(p1)
  }
  dev.off()
  pdf(paste0("6_analysis/",samplei,"/DMR_plot/1.",groupi,".5hmC_level_0.1_0.9.pdf"))
  for (si in sample.ls) {
    p1 = ggplot(DMR_mlevel.dfi[DMR_mlevel.dfi$orig.ident==si,],aes(x=nB,y=nA,fill=h_level)) +
      geom_tile(show.legend=T,width=0.9,height=0.9) + theme_void() +
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",si)) +
      expand_limits(y=c(0,96),x=c(0,96)) + coord_fixed() +
      scale_fill_gradientn(colours=c(paletteer_c("grDevices::Viridis", 30,direction = -1)))
    print(p1)
  }
  dev.off()
  }


