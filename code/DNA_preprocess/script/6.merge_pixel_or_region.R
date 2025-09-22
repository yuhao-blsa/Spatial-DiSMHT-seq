library(ggplot2)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]
parse_config <- function(path) {
  lines <- readLines(path)
  lines <- lines[!grepl("^\\s*#", lines)] 
  lines <- lines[nzchar(lines)]  
  
  config <- list()
  for(line in lines) {
    if(grepl("=", line)) {
      pair <- strsplit(line, "=")[[1]]
      key <- trimws(pair[1])
      value <- trimws(pair[2])
      config[[key]] <- gsub("^['\"]|['\"]$", "", value)  
    }
  }
  return(config)
}
config<-parse_config(config_file)

sample_id = config$sample_id
outpath=config$outpath

bedtools=config$bedtools_path
samplei = config$sample_id
# 1 filter input
CpG.df = read.table(paste0(outpath,"5_Plot/CpG_site.single.",samplei,".bed"),stringsAsFactors = F)
CpG.df$V7 = gsub("_(\\d+)_","_\\1x",CpG.df$V7) 

system(paste0("mkdir -p ",outpath,"6_analysis/",samplei,"/CpG_level/track"))
system(paste0("mkdir -p ",outpath,"6_analysis/",samplei,"/gene_level"))

# get CpG-promoter link
site.info = CpG.df[,1:3]
site.info$id = paste0(site.info$V1,"_",site.info$V2,"_",site.info$V3)
site.info = site.info[!duplicated(site.info),]
write.table(site.info[,1:3],file = paste0(outpath,"6_analysis/",samplei,"/CpG_level/CpG_site.bed"),row.names = F,col.names = F,quote = F,sep = "\t")
system(paste0("cd ",outpath,"6_analysis/",samplei,"/CpG_level;",bedtools," intersect -a CpG_site.bed -b ",config$promoter," -wa -wb > CpG_site_and_promoter.bed"))
system(paste0("cd ",outpath,"6_analysis/",samplei,"/CpG_level;",bedtools," intersect -a CpG_site.bed -b ",config$genome_2kb," -wa -wb > CpG_site_and_2kb_bin.bed"))
system(paste0("cd ",outpath,"6_analysis/",samplei,"/CpG_level;",bedtools," intersect -a CpG_site.bed -b ",config$famiy," -wa -wb > CpG_site_and_repeat_family.bed"))
system(paste0("cd ",outpath,"6_analysis/",samplei,"/CpG_level;",bedtools," intersect -a CpG_site.bed -b ",config$CGI," -wa -wb > CpG_site_and_other_family.bed"))
promoter_info = read.table(paste0(outpath,"6_analysis/",samplei,"/CpG_level/CpG_site_and_promoter.bed"))
windows_info2 = read.table(paste0(outpath,"6_analysis/",samplei,"/CpG_level/CpG_site_and_2kb_bin.bed"))
windows_info2$V8 = paste0(windows_info2$V4,"_",windows_info2$V5,"_",windows_info2$V6)
repeat_info = read.table(paste0(outpath,"6_analysis/",samplei,"/CpG_level/CpG_site_and_repeat_family.bed"))
gap_info = read.table(paste0(outpath,"6_analysis/",samplei,"/CpG_level/CpG_site_and_other_family.bed"))

library(dplyr)
merge_modif_data <- function(CpG.df,cell_metadata=NULL,promoter_info=NULL){
  CpG.df$site = paste0(CpG.df$V1,"_",CpG.df$V2,"_",CpG.df$V3)
  if(!is.null(cell_metadata) & is.null(promoter_info)){
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
    CpG.df$ct = cell_metadata$seurat_clusters[match(CpG.df$V7,cell_metadata$coord)]
    promoter_info$site = paste0(promoter_info$V1,"_",promoter_info$V2,"_",promoter_info$V3)
    CpG.df$gene = promoter_info$V8[match(CpG.df$site,promoter_info$site)]
    CpG.df = CpG.df[!is.na(CpG.df$gene),]
    CpG.df$id = paste0("C",CpG.df$ct,"_",CpG.df$gene)
  }
  output = CpG.df[,c("V1","V2","V3","V4","V5","V6","id")]

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
  return(output)
}
write.csv(merge_modif_data(CpG.df,cell_metadata=NULL,promoter_info),paste0(outpath,"6_analysis/",samplei,"/gene_level/0_all_pixel.promoter.modif_level.csv"))
write.csv(merge_modif_data(CpG.df,cell_metadata=NULL,windows_info2),paste0(outpath,"6_analysis/",samplei,"/gene_level/0_all_pixel.2kb_bin.modif_level.csv"))
write.csv(merge_modif_data(CpG.df,cell_metadata=NULL,repeat_info),paste0(outpath,"6_analysis/",samplei,"/gene_level/0_all_pixel.repeat_family.modif_level.csv"))
write.csv(merge_modif_data(CpG.df,cell_metadata=NULL,gap_info),paste0(outpath,"6_analysis/",samplei,"/gene_level/0_all_pixel.gap_family.modif_level.csv"))


sample.ls = read.table(config$index)$V2
element_mlevel.total = rbind(read.csv(paste0(outpath,"6_analysis/",samplei,"/gene_level/0_all_pixel.repeat_family.modif_level.csv"),row.names = 1),
                             read.csv(paste0(outpath,"6_analysis/",samplei,"/gene_level/0_all_pixel.gap_family.modif_level.csv"),row.names = 1))
element_mlevel.total$coord = gsub("^(.*\\d+x\\d+)_.*","\\1",element_mlevel.total$id)
element_mlevel.total$group = gsub("^.*\\d+x\\d+_(.*)","\\1",element_mlevel.total$id)
element_mlevel.total$orig.ident = gsub("(^.*)_\\d+x\\d+$","\\1",element_mlevel.total$coord)
element_mlevel.total$nA = as.numeric(gsub("^.*_(\\d+)x\\d+$","\\1",element_mlevel.total$coord))
element_mlevel.total$nB = as.numeric(gsub("^.*_\\d+x(\\d+)$","\\1",element_mlevel.total$coord))
element_mlevel.total$m_level = element_mlevel.total$V5/apply(element_mlevel.total[,c("V4","V5","V6")],1,sum)
element_mlevel.total$h_level = element_mlevel.total$V6/apply(element_mlevel.total[,c("V4","V5","V6")],1,sum)

system(paste0("mkdir -p ",outpath,"6_analysis/",samplei,"/element_plot"))


library(paletteer)
for (groupi in sort(unique(element_mlevel.total$group))) {
  element_mlevel.dfi = element_mlevel.total[element_mlevel.total$group==groupi,]
  element_mlevel.dfi$m_level = pmin(pmax(element_mlevel.dfi$m_level,quantile(element_mlevel.dfi$m_level,0.05)),quantile(element_mlevel.dfi$m_level,0.95))
  element_mlevel.dfi$h_level = pmin(pmax(element_mlevel.dfi$h_level,quantile(element_mlevel.dfi$h_level,0.05)),quantile(element_mlevel.dfi$h_level,0.95))
  pdf(paste0(outpath,"6_analysis/",samplei,"/element_plot/1.",groupi,".5mC_level_0.05_0.95.pdf"))
  for (si in sample.ls) {
    p1 = ggplot(element_mlevel.dfi[element_mlevel.dfi$orig.ident==si,],aes(x=nB,y=nA,fill=m_level)) +
      geom_raster(show.legend=T) + theme_void() +
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",si)) +
      expand_limits(y=c(0,96),x=c(0,96)) + coord_fixed() +
      scale_fill_gradientn(colours=c(paletteer_c("grDevices::Viridis", 30,direction = -1)))+
      ggtitle(si) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p1)
  }
  dev.off()
  pdf(paste0(outpath,"6_analysis/",samplei,"/element_plot/1.",groupi,".5hmC_level_0.05_0.95.pdf"))
  for (si in sample.ls) {
    p1 = ggplot(element_mlevel.dfi[element_mlevel.dfi$orig.ident==si,],aes(x=nB,y=nA,fill=h_level)) +
      geom_raster(show.legend=T) + theme_void() +
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",si)) +
      expand_limits(y=c(0,96),x=c(0,96)) + coord_fixed() +
      scale_fill_gradientn(colours=c(paletteer_c("grDevices::Viridis", 30,direction = -1)))+
      ggtitle(si) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p1)
  }
  dev.off()
}

write.csv("Succeed",paste0(outpath,"log/6.merge_pixel_or_region.log"))



