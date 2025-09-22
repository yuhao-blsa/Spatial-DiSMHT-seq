# conda activate plot_env

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
dicti7 = read.table(config$index)$V2
outpath=config$outpath
min_CpG = as.numeric(config$min_CpG)
nchannel = as.numeric(config$channel)

##CpG soft data
id.ls = read.table(paste0(outpath,"4_modkit_perpixel_CpG/pixel_id.",sample_id,"_soft_threshold.txt"))$V1
id.ls = gsub("^.*/split_|\\.bed$", "", id.ls)
id.m =  do.call(rbind, strsplit(id.ls, "_"))
pixel.df = data.frame("id" = id.ls,
                      "sample" = id.m[,1],
                      "nA" = id.m[,2],
                      "nB" = id.m[,3],
                      "CpG" = 0,
                      "C_level" = 0,
                      "m_level" = 0,
                      "h_level" = 0)
CpG.df = matrix(nrow = 0, ncol = 7)

tmp <- mclapply(id.ls, function(i){
  input.df = read.table(paste0(outpath,"4_modkit_perpixel_CpG/",sample_id,"_soft_threshold/split_",i,".bed"),stringsAsFactors = F)
  CpG.i = input.df[input.df$V4=="m",c(1:3,13,12)]
  colnames(CpG.i) = c("chr","start","end","C","m")
  CpG.i$h = input.df[input.df$V4=="h",12]
  CpG.i$pixel = i
  CpG.i[4:6] = CpG.i[4:6]/apply(CpG.i[4:6],1,sum)
  return(list(CpG.i, c(nrow(CpG.i),apply(CpG.i[4:6],2,mean))))
},mc.cores = 12)

CpG.df = do.call(rbind, lapply(tmp, function(x) x[[1]]))
pixel.df[,5:8] = do.call(rbind, lapply(tmp, function(x) x[[2]]))
write.table(CpG.df,paste0(outpath,"5_Plot/CpG_site.single.",sample_id,".bed"),row.names = F,col.names = F,quote = F,sep = "\t")
write.table(pixel.df,paste0(outpath,"5_Plot/CpG_site.pixel_summary.",sample_id,".tsv"),row.names = F,col.names = F,quote = F,sep = "\t")

# plot
library(paletteer)
pixel.df$nA = as.numeric(pixel.df$nA)
pixel.df$nB = as.numeric(pixel.df$nB)
pdf(paste0(outpath,"5_Plot/2.CpG_coverage.",sample_id,".pdf"))
for (i7_i in dicti7) {
  #i7_i<-"i702"
  p1 = ggplot(pixel.df[pixel.df$sample==i7_i,],aes(x=nB,y=nA,fill=CpG)) +
    geom_raster(show.legend=T) + theme_void() +
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",i7_i)) +
    expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() +
    scale_fill_gradientn(colours=c(paletteer_c("grDevices::Blues 3", 30,direction = -1),rep("black",3)))+
    ggtitle(i7_i) +
    theme(plot.title = element_text(hjust = 0.5))
  print(p1)
}
dev.off()

pixel.df$m_level = pmin(pmax(pixel.df$m_level,quantile(pixel.df$m_level,0.05)),quantile(pixel.df$m_level,0.95))
pixel.df$h_level = pmin(pmax(pixel.df$h_level,quantile(pixel.df$h_level,0.05)),quantile(pixel.df$h_level,0.95))
pdf(paste0(outpath,"5_Plot/2.5mC_level_0.05_0.95.",sample_id,".pdf"))
for (i7_i in dicti7) {
  p1 = ggplot(pixel.df[pixel.df$sample==i7_i & pixel.df$CpG>min_CpG,],aes(x=nB,y=nA,fill=m_level)) +
    geom_raster(show.legend=T) + theme_void() +
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",i7_i)) +
    expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() +
    scale_fill_gradientn(colours=c(paletteer_c("grDevices::Viridis", 30,direction = -1)))+
    ggtitle(i7_i) +
    theme(plot.title = element_text(hjust = 0.5))
  print(p1)
}
dev.off()
pdf(paste0(outpath,"5_Plot/2.5hmC_level_0.05_0.95.",sample_id,".pdf"))
for (i7_i in dicti7) {
  p1 = ggplot(pixel.df[pixel.df$sample==i7_i & pixel.df$CpG>min_CpG,],aes(x=nB,y=nA,fill=h_level)) +
    geom_raster(show.legend=T) + theme_void() +
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",i7_i)) +
    expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() +
    scale_fill_gradientn(colours=c(paletteer_c("grDevices::Viridis", 30,direction = -1)))+
    ggtitle(i7_i) +
    theme(plot.title = element_text(hjust = 0.5))
  print(p1)
}
dev.off()

tmp.df = pixel.df[pixel.df$CpG>min_CpG,]
cat(sample_id,"median CpG cover:\n", file = paste0(outpath,"Log.txt"), append = TRUE)
cat(names(tapply(tmp.df$CpG, tmp.df$sample, median)),"\n", file =paste0(outpath,"Log.txt"), append = TRUE)
cat(tapply(tmp.df$CpG, tmp.df$sample, median),"\n", file =paste0(outpath,"Log.txt"), append = TRUE)
cat(sample_id,"mean m_level:\n", file =paste0(outpath,"Log.txt"), append = TRUE)
cat(names(tapply(tmp.df$m_level, tmp.df$sample, mean)),"\n", file = paste0(outpath,"Log.txt"), append = TRUE)
cat(tapply(tmp.df$m_level, tmp.df$sample, mean),"\n", file = paste0(outpath,"Log.txt"), append = TRUE)
cat(sample_id,"mean h_level:\n", file = paste0(outpath,"Log.txt"), append = TRUE)
cat(names(tapply(tmp.df$h_level, tmp.df$sample, mean)),"\n", file = paste0(outpath,"Log.txt"), append = TRUE)
cat(tapply(tmp.df$h_level, tmp.df$sample, mean),"\n", file = paste0(outpath,"Log.txt"), append = TRUE)


CpG.df$sample = pixel.df$sample[match(CpG.df$pixel,pixel.df$id)]
cat(sample_id,"bulk m_level:\n", file = paste0(outpath,"Log.txt"), append = TRUE)
cat(names(tapply(CpG.df$m, CpG.df$sample, mean)),"\n", file = paste0(outpath,"Log.txt"), append = TRUE)
cat(tapply(CpG.df$m, CpG.df$sample, mean),"\n", file = paste0(outpath,"Log.txt"), append = TRUE)
cat(sample_id,"bulk h_level\n", file =paste0(outpath,"Log.txt"), append = TRUE)
cat(names(tapply(CpG.df$h, CpG.df$sample, mean)),"\n", file = paste0(outpath,"Log.txt"), append = TRUE)
cat(tapply(CpG.df$h, CpG.df$sample, mean),"\n", file =paste0(outpath,"Log.txt"), append = TRUE)


###nonCpG soft data
for (type in c("CHG","CHH")) {
  id.ls = read.table(paste0(outpath,"4_modkit_perpixel_nonCpG/pixel_id.",sample_id,"_",type,"_soft_threshold.txt"))$V1
  id.ls = gsub("^.*/split_|\\.bed$", "", id.ls)
  id.m =  do.call(rbind, strsplit(id.ls, "_"))
  pixel.df = data.frame("id" = id.ls,
                        "sample" = id.m[,1],
                        "nA" = id.m[,2],
                        "nB" = id.m[,3],
                        "CpG" = 0,
                        "C_level" = 0,
                        "m_level" = 0,
                        "h_level" = 0)
  CpG.df = matrix(nrow = 0, ncol = 7)
  
  tmp <- mclapply(id.ls, function(i){
    # fast read
    input.df = data.table::fread(paste0(outpath,"4_modkit_perpixel_nonCpG/",sample_id,"_",type,"_soft_threshold/split_",i,".bed"),stringsAsFactors = F)
    input.df$V4 = gsub(",.*$","",input.df$V4)
    CpG.i = input.df[input.df$V4=="m",c(1:3,13,12)]
    colnames(CpG.i) = c("chr","start","end","C","m")
    CpG.i$h = input.df[input.df$V4=="h",12]
    CpG.i$pixel = i
    CpG.i[,4:6] = CpG.i[,4:6]/apply(CpG.i[,4:6],1,sum)
    return(list(1, c(nrow(CpG.i),apply(CpG.i[,4:6],2,mean))))
  },mc.cores = 20)
  
  pixel.df[,5:8] = do.call(rbind, lapply(tmp, function(x) x[[2]]))
  write.table(pixel.df,paste0(outpath,"5_Plot/CpG_site.pixel_summary.",sample_id,"_",type,".tsv"),row.names = F,col.names = F,quote = F,sep = "\t")

  # plot
  pixel.df$nA = as.numeric(pixel.df$nA)
  pixel.df$nB = as.numeric(pixel.df$nB)
  pdf(paste0(outpath,"5_Plot/3.CpG_coverage.",sample_id,"_",type,".pdf"))
  for (i7_i in dicti7) {
    p1 = ggplot(pixel.df[pixel.df$sample==i7_i,],aes(x=nB,y=nA,fill=CpG)) +
      geom_raster(show.legend=T) + theme_void() +
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",i7_i)) +
      expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() +
      scale_fill_gradientn(colours=c(paletteer_c("grDevices::Blues 3", 30,direction = -1),rep("black",3)))+
      ggtitle(i7_i) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p1)
  }
  dev.off()
  
  pixel.df$m_level = pmin(pmax(pixel.df$m_level,quantile(pixel.df$m_level,0.05)),quantile(pixel.df$m_level,0.95))
  pixel.df$h_level = pmin(pmax(pixel.df$h_level,quantile(pixel.df$h_level,0.05)),quantile(pixel.df$h_level,0.95))
  pdf(paste0(outpath,"5_Plot/3.5mC_level_0.05_0.95.",sample_id,"_",type,".pdf"))
  for (i7_i in dicti7) {
    p1 = ggplot(pixel.df[pixel.df$sample==i7_i,],aes(x=nB,y=nA,fill=m_level)) +
      geom_raster(show.legend=T) + theme_void() +
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",i7_i)) +
      expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() +
      scale_fill_gradientn(colours=c(paletteer_c("grDevices::Viridis", 30,direction = -1)))+
      ggtitle(i7_i) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p1)
  }
  dev.off()
  pdf(paste0(outpath,"5_Plot/3.5hmC_level_0.05_0.95.",sample_id,"_",type,".pdf"))
  for (i7_i in dicti7) {
    p1 = ggplot(pixel.df[pixel.df$sample==i7_i,],aes(x=nB,y=nA,fill=h_level)) +
      geom_raster(show.legend=T) + theme_void() +
      scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
      scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",i7_i)) +
      expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() +
      scale_fill_gradientn(colours=c(paletteer_c("grDevices::Viridis", 30,direction = -1)))+
      ggtitle(i7_i) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p1)
  }
  dev.off()
}

##map ref length
data.df = read.table(paste0(outpath,"5_Plot/1_mapped_ref_length.summary.",sample_id,".tsv"))
id.m =  do.call(rbind, strsplit(data.df$V1, "_"))
tmp.df = read.table(paste0(outpath,"5_Plot/CpG_site.pixel_summary.",sample_id,".tsv"),stringsAsFactors = F,col.names = c("id","sample","nA","nB","CpG","C_level","m_level","h_level"))
tmp.df = tmp.df[tmp.df$CpG>min_CpG,]

pixel.df = data.frame("id" = data.df$V1,
                      "sample" = id.m[,1],
                      "nA" = id.m[,2],
                      "nB" = id.m[,3],
                      "bp" = data.df$V2)
# plot
pixel.df$nA = as.numeric(pixel.df$nA)
pixel.df$nB = as.numeric(pixel.df$nB)
pdf(paste0(outpath,"5_Plot/3_mapped_ref_length.",sample_id,".pdf"))
for (i7_i in dicti7) {
  p1 = ggplot(pixel.df[pixel.df$sample==i7_i,],aes(x=nB,y=nA,fill=bp)) +
    geom_raster(show.legend=T) + theme_void() +
    scale_y_continuous(trans = "reverse",breaks = seq(0,95,5)) +
    scale_x_continuous(trans = "reverse",breaks = seq(0,95,5),name = paste("nB",i7_i)) +
    expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() +
    scale_fill_gradientn(colours=c(paletteer_c("grDevices::Blues 3", 30,direction = -1),rep("black",3)))+
    ggtitle(i7_i) +
    theme(plot.title = element_text(hjust = 0.5))
  print(p1)
}
dev.off()

cat(sample_id,"total mapped length:\n", file =paste0(outpath,"Log.txt"), append = TRUE)
cat(names(tapply(pixel.df$bp, pixel.df$sample, sum)),"\n", file = paste0(outpath,"Log.txt"), append = TRUE)
cat(tapply(pixel.df$bp, pixel.df$sample, sum)/1e6,"\n", file = paste0(outpath,"Log.txt"), append = TRUE)

pixel.df = pixel.df[pixel.df$id %in% tmp.df$id,]
cat(sample_id,"median mapped length:\n", file = paste0(outpath,"Log.txt"), append = TRUE)
cat(names(tapply(pixel.df$bp, pixel.df$sample, median)),"\n", file =paste0(outpath,"Log.txt"), append = TRUE)
cat(tapply(pixel.df$bp, pixel.df$sample, median)/1e3,"\n", file =paste0(outpath,"Log.txt"), append = TRUE)

write.csv("Succeed",paste0(outpath,"log/5.summarize.log"))
