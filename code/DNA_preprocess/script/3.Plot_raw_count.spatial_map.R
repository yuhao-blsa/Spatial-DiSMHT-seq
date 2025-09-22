library(ggplot2)

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
nchannel = as.numeric(config$channel)
pixel.df = read.table(paste0(outpath,"3_corrected/3_out_id2coord.",sample_id,".tsv"),stringsAsFactors = F,col.names = c("id","ni7","nb1","nb2","coord"))
pixel.df = pixel.df[,2:5]

# summarize data
pixel.summary = by(pixel.df,pixel.df$coord,nrow)
pixel.summary.df = pixel.df[!duplicated(pixel.df),]
pixel.summary.df = pixel.summary.df[order(pixel.summary.df$coord),]
pixel.summary.df$num = as.vector(pixel.summary[pixel.summary.df$coord])

# plot
library(paletteer)
ref.df1 = read.table(config$index)
dicti7 = ref.df1$V2
pdf(paste0(outpath,"3_corrected/4.raw_read_distribution.",sample_id,".pdf"))
for (i7_i in dicti7) {
  p1 = ggplot(pixel.summary.df[pixel.summary.df$ni7==i7_i,],aes(x=nb2,y=nb1,color=num)) +
    geom_point(shape=15,size=1.9,show.legend=T) + theme_void() +
    scale_y_continuous(trans = "reverse") +
    scale_x_continuous(trans = "reverse",name = paste("nB",i7_i)) +
    expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() +
    scale_color_gradientn(colours=c(paletteer_c("grDevices::Blues 3", 30,direction = -1),rep("black",3)))+
    ggtitle(i7_i) +
    theme(plot.title = element_text(hjust = 0.5))
  print(p1)
}
dev.off()

pdf(paste0(outpath,"3_corrected/4.raw_read_distribution.geom_tile.",sample_id,".pdf"))
for (i7_i in dicti7) {
  p1 = ggplot(pixel.summary.df[pixel.summary.df$ni7==i7_i,],aes(x=nb2,y=nb1,fill=num)) +
    geom_tile(color = "white",lwd = 0,linetype = 0,show.legend=T) + theme_void() + 
    scale_y_continuous(trans = "reverse") +
    scale_x_continuous(trans = "reverse",name = paste("nB",i7_i)) +
    expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() +
    scale_fill_gradientn(colours=c(paletteer_c("grDevices::Blues 3", 30,direction = -1),rep("black",3)))+
    ggtitle(i7_i) +
    theme(plot.title = element_text(hjust = 0.5))
  print(p1)
}
dev.off()

pdf(paste0(outpath,"3_corrected/4.raw_read_distribution.geom_raster.",sample_id,".pdf")) 
for (i7_i in dicti7) {
  p1 = ggplot(pixel.summary.df[pixel.summary.df$ni7==i7_i,],aes(x=nb2,y=nb1,fill=num)) +
    geom_raster(show.legend=T) + theme_void() + 
    scale_y_continuous(trans = "reverse") +
    scale_x_continuous(trans = "reverse",name = paste("nB",i7_i)) +
    expand_limits(y=c(0,nchannel),x=c(0,nchannel)) + coord_fixed() +
    scale_fill_gradientn(colours=c(paletteer_c("grDevices::Blues 3", 30,direction = -1),rep("black",3)))+
    ggtitle(i7_i) +
    theme(plot.title = element_text(hjust = 0.5))
  print(p1)
}
dev.off()

write.csv("Succeed",paste0(outpath,"log/3.correct_barcode.log"))
