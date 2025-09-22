library(ShortRead)
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

config <- parse_config(config_file)

sample_id = config$sample_id
outpath=config$outpath
# Load data
fastq_file <- paste0(outpath,"2_skewer_trim/3.I7barcode2barcode1.",sample_id,".fq.gz")
fastq_data <- readFastq(fastq_file)

# Correction Run
char_seq <- as.character(sread(fastq_data)) 
char_seq.df = data.frame("seq"=char_seq, "i7"=substr(char_seq, 1, 8), "b2"=substr(char_seq, 9, 16), "b1"=substr(char_seq, 17,24))
ref.ls1 = read.table(config$index)$V1
ref.ls2 = read.table(config$barcode1)$V1
ref.ls3 = read.table(config$barcode2)$V1

# Use Hamming distance
# compare the number of different characters in two strings
getRefMinDist <- function(ref.ls){
  nref = length(ref.ls)
  dist.ls = c()
  for (i in 1:(nref-1)) {
    for (j in (i+1):nref) {
      dist.ls = c(dist.ls, sum(strsplit(ref.ls[i], NULL)[[1]] != strsplit(ref.ls[j], NULL)[[1]]))
    }
  }
  return(min(dist.ls))
}
getSimilarRef <- function(chri,ref.ls,max_d=1,expand=1){
  dist.ls = c()
  for (i in 1:length(ref.ls)) {
    dist.ls = c(dist.ls, sum(strsplit(chri, NULL)[[1]] != strsplit(ref.ls[i], NULL)[[1]]))
  }
  dist.sort = sort(dist.ls)
  if(min(dist.ls)<=max_d){
    return(ref.ls[which(dist.ls<=max_d)[1]])
  }else if(min(dist.ls)<=(max_d+expand) & dist.sort[1]<dist.sort[2]){
    return(ref.ls[which(dist.ls==dist.sort[1])])
  }else{
    return(strrep("N", nchar(chri)))
  }}

char_seq.df$i7_correct = unlist(mclapply(char_seq.df$i7, getSimilarRef, ref.ls = ref.ls1, max_d=as.numeric(config$i7_correct_max_d), expand=as.numeric(config$i7_correct_expand),mc.cores = as.numeric(config$threads)))
char_seq.df$b2_correct = unlist(mclapply(char_seq.df$b2, getSimilarRef, ref.ls = ref.ls3, max_d=as.numeric(config$b2_correct_max_d), expand=as.numeric(config$b2_correct_expand),mc.cores = as.numeric(config$threads)))
char_seq.df$b1_correct = unlist(mclapply(char_seq.df$b1, getSimilarRef, ref.ls = ref.ls2, max_d=as.numeric(config$b1_correct_max_d), expand=as.numeric(config$b1_correct_expand),mc.cores = as.numeric(config$threads)))

char_seq.df$seq_correct = paste0(char_seq.df$i7_correct, char_seq.df$b2_correct, char_seq.df$b1_correct)
fastq_data@sread  <- DNAStringSet(char_seq.df$seq_correct)


# Save
correct_fastq_file=paste0(outpath,"3_corrected/1.I7barcode2barcode1.extract.corrected.",sample_id,".fq.gz")
writeFastq(fastq_data[char_seq.df$i7_correct !="NNNNNNNN"  & char_seq.df$b2_correct!="NNNNNNNN" & char_seq.df$b1_correct!="NNNNNNNN"], correct_fastq_file,mode = "w", compress = TRUE) 
cat(sprintf(paste0(sample_id," barcode total: %d\n"), dim(char_seq.df)[1]), file =paste0(outpath,"Log.txt"), append = TRUE)
cat(sprintf(paste0(sample_id," barcode right: %d\n"), sum(char_seq.df$i7_correct !="NNNNNNNN"  & char_seq.df$b2_correct!="NNNNNNNN" & char_seq.df$b1_correct!="NNNNNNNN")), file = paste0(outpath,"Log.txt"), append = TRUE)

# Save id2coord
id.ls <- as.character(id(fastq_data)) 
id.ls = gsub("_[ATCG]{2}$","",id.ls)
char_seq.df$id = id.ls
pixel.df = char_seq.df[char_seq.df$i7_correct!="NNNNNNNN" & char_seq.df$b2_correct!="NNNNNNNN" & char_seq.df$b1_correct!="NNNNNNNN",c(9,5:7)]

dicti7 = read.table(config$index)$V2
names(dicti7) = ref.ls1
dictb1 = 1:length(ref.ls2)
names(dictb1) = ref.ls2
dictb2 = 1:length(ref.ls3)
names(dictb2) = ref.ls3

pixel.df$ni7 = dicti7[pixel.df$i7_correct]
pixel.df$nb1 = dictb1[pixel.df$b1_correct]
pixel.df$nb2 = dictb2[pixel.df$b2_correct]
pixel.df$coord = paste(pixel.df$ni7,pixel.df$nb1,pixel.df$nb2,sep = "_")
write.table(pixel.df[,c(1,5:8)],paste0(outpath,"3_corrected/3_out_id2coord.",sample_id,".tsv"),row.names = F,col.names = F,quote = F,sep = "\t")
