library(rGREAT)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

samplei = "HEM23"
bin_size = "2kb"
# Ke LoH + Adult
color.ls5a = c("#7CAE00","#ffdf00","#AD0000","#7F7F7F","#C09B00","#ff4500","#00BAE0","#A3A500","#000080","#90EE90",
               "#EA8331","#00B0F6","#FF0000","#704241","#FFFF00","#C77CFF","#daa520","#adff2f","#8D2282","#E76BF3",
               "#39B600","#9590FF","#FF6A98","#FA62DB","#4169E1","#FF9999","#B22222","#FF62BC","#FFD032","#A2FF00",
               "#006400","#c0c0c0","#004A00","#FFB676","#87CEFA","#B34E00","#63E463","#6A5ACD")
names(color.ls5a) = paste0("C",0:37)
## summary
total_result.df = data.frame()
for (cti in c(3,8)) {
  gr=import(paste0("6_analysis/",samplei,"/DMR/4_C",cti,".",bin_size,"_bin.DMR.hypo.bed"),format = "bed")
  gr <- gr[ !grepl("(_fix$|_alt$)", seqnames(gr)) ]
  job=submitGreatJob(gr,species="hg38")
  result = getEnrichmentTables(job)
  result = result$`GO Biological Process`
  result = result[result$Binom_Fold_Enrichment>2 & result$Binom_Adjp_BH<0.05,c(1,2,8)]
  result$logp = -log10(result$Binom_Raw_PValue)
  result$class = paste0("C",cti)
  colnames(result)=c("GO_ID","GO_Term","Pvalue","logP","Sample")
  total_result.df = rbind(total_result.df,result)
}
write.csv(total_result.df,file = paste0("6_analysis/",samplei,"/DMR/7_DMR_GREAT_result.hypo.csv"))

top10.df = Reduce(rbind,lapply(c(13,6,9,15,16,14,5,1,11,2,7,10,12,4,17,0,3,8),function(i){
  tmp = total_result.df[total_result.df$Sample==paste0("C",i),]
  print(nrow(tmp))
  return(tmp[1:min(10,nrow(tmp)),])
}))
top10.df$name = paste(top10.df$Sample,top10.df$GO_Term,sep = "-")
top10.df$name=factor(top10.df$name,levels = rev(unique(top10.df$name)))

pdf(paste0("6_analysis/",samplei,"/DMR/7_DMR_GREAT_result.top10.hypo.pdf"),width = 10,height = 20)
p2 <- ggplot(top10.df, aes(x = name, y = logP, fill = Sample)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  coord_cartesian(ylim = c(0,45)) + 
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) + scale_fill_manual(values = color.ls5a)
print(p2)
dev.off()
pdf(paste0("6_analysis/",samplei,"/DMR/7_DMR_GREAT_result.top10.hypo.sub2.pdf"),width = 10,height = 10)
p2 <- ggplot(top10.df[1:90,], aes(x = name, y = logP, fill = Sample)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  coord_cartesian(ylim = c(0,45)) + 
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) + scale_fill_manual(values = color.ls5a)
print(p2)
p2 <- ggplot(top10.df[91:180,], aes(x = name, y = logP, fill = Sample)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  coord_cartesian(ylim = c(0,45)) + 
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) + scale_fill_manual(values = color.ls5a)
print(p2)
dev.off()


## Manual selection
select_result.df = read.csv(paste0("6_analysis/HEM23/DMR/7_DMR_GREAT_result.hypo.final_filtered_top3.csv"))

select_result.df$name = paste(select_result.df$Sample,select_result.df$GO_Term,sep = "-")
select_result.df$name=factor(select_result.df$name,levels = rev(unique(select_result.df$name)))

pdf(paste0("6_analysis/",samplei,"/DMR/8_DMR_GREAT_result.top3_selected.hypo.pdf"),width = 8,height = 8)
p2 <- ggplot(select_result.df, aes(x = name, y = logP, fill = Sample)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  coord_cartesian(ylim = c(0,45)) + 
  coord_flip() +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) + scale_fill_manual(values = color.ls5a)
print(p2)
dev.off()







#remove the 冗余词条；
k=unique(merdat_1[,c(1,2)])
#input the GO ID  to http://revigo.irb.hr/ ; download the table, select the "NULL" in "representative column"
sele_term=read.table(file = "Filtered_GO_term.txt",header = F,stringsAsFactors = F,sep=" ")
merdat_1=merdat_1[merdat_1$ID %in% sele_term$V1,]

####
# 使用ggplot2绘制热图
ggplot(result, aes(x = Sample, y = GO_Term, fill = logP)) +
  geom_tile(color = "white") +                              # 绘制热图单元格，添加白色边框
  scale_fill_gradient(low = "white", high = "red", name = "-log10(P-value)") + # 颜色梯度
  labs(title = "GO Enrichment Heatmap", x = "Sample", y = "GO Biological Process") + # 标题和轴标签
  theme_minimal() +                                         # 使用简洁主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),      # X轴标签倾斜
    axis.text.y = element_text(size = 10),                 # Y轴标签字体大小
    plot.title = element_text(size = 14, face = "bold")    # 标题加粗
  )
