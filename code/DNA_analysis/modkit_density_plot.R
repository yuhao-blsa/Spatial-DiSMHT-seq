library(dplyr)
library(tidyr)

result_data <- read.delim("read_mod_counts_filter.tsv", stringsAsFactors=FALSE)
name_data <- read.delim("3_out_id2coord.tsv", header=FALSE,stringsAsFactors=FALSE)

merged_data <- result_data %>%
  left_join(
    name_data %>% select(read_id = V1, name = V5,sample=V2),
    by = "read_id"
  )
merged_data$total_sites = merged_data$X5mC + merged_data$X5hmC + merged_data$Unmodified
summary_stats <- merged_data %>%
  group_by(sample) %>%
  summarise(
    total_methylated = sum(X5mC),
    total_sites = sum(total_sites),
    methylation_rate = total_methylated / total_sites * 100 
  )

cluster_info <- read.csv("5_C5_metadata.csv",header=TRUE)
library(stringr)

cluster_info1$name <- str_replace(cluster_info1$coord, "testis_1", "i701") %>%
  str_replace("x","_")
cluster_info2$name <- str_replace(cluster_info2$X, "E9_1", "i703") %>%
  str_replace("E9_2", "i704")%>%
  str_replace("x","_")


merged_data_testis <- merged_data %>%
  left_join(
    cluster_info1 %>% select(name,C10),  
    by = "name"
  )

merged_data_E9 <- merged_data %>%
  left_join(
    cluster_info2 %>% select(name,seurat_clusters),  
    by = "name"
  )

testis_cluster <- data.frame(
  C10 =  0:9,  # C0, C1, ..., C9
  cluster_name = c(
    "Spermatogonial population",
    "Spermatocyte",
    "Intermediate",
    "Acrosoma",
    "Round spermatids",
    "Elongating spermatids1",
    "Elongating spermatids2",
    "Sperm",
    "Sertoli",
    "Stroma"
  )
)

mouse_cluster <- read.csv("mouse_cluster.txt",header=FALSE)
mouse_cluster$seurat_clusters <- sapply(strsplit(mouse_cluster$V1, "_"), `[`, 1)
mouse_cluster$cluster_name <- sapply(strsplit(mouse_cluster$V1, "_"), `[`, 2)
merged_data3 <- merged_data_E9  %>%
  left_join(
    mouse_cluster %>% 
      select(cluster_name, seurat_clusters) %>% 
      mutate(seurat_clusters = as.integer(seurat_clusters)),  # 强制转换为整数
    by = "seurat_clusters"
  )

merged_data3$ratio_5mC <- merged_data3$X5mC/merged_data3$total_sites
merged_data3$ratio_5hmC <- merged_data3$X5hmC/merged_data3$total_sites
write.table(merged_data3, "E9_read_mod_counts_with_info.tsv", sep = "\t", row.names = FALSE,quote = FALSE)
merged_data3<-read.csv("read_mod_counts_with_info.tsv", sep = "\t")
filtered_data <- merged_data3 %>% 
  filter(!is.na(cluster_name))

hindbrain <-read.csv("5_C5_metadata.csv",sep=",")
hindbrain$name <- str_replace(hindbrain$X, "GW6_1", "i702") %>%
 str_replace("GW6_2", "i703")%>%
  str_replace("x","_")

merged_data_hindbrain <- merged_data3 %>%
  filter(cluster_name == "hindbrain")%>%
  left_join(
    hindbrain %>%  dplyr::select(name,seurat_clusters),  
    by = "name")

top_hindbrainpoints <- merged_data_hindbrain%>%
  arrange(desc(total_sites)) %>%
  slice_head(n = 5000) %>%
  ungroup()

top_hindbrainpoints_2500 <- merged_data_hindbrain %>%
  filter(seurat_clusters %in% c("1", "2")) %>%
  group_by(seurat_clusters) %>%
  arrange(desc(total_sites)) %>%
  slice_head(n = 2500) %>%
  ungroup()

filtered_data <- merged_data_hindbrain %>% 
  filter(!is.na(cluster_name))

cluster_counts <- filtered_data %>%
  group_by(cluster_name) %>%
  summarise(
    n_reads = n(),
    mean_5mC = mean(ratio_5mC),
    mean_5hmC = mean(ratio_5hmC)
  ) %>%
  arrange(desc(n_reads)) 


color.ls5a = c("#7CAE00","#ffdf00","#AD0000","#7F7F7F","#C09B00","#ff4500","#00BAE0","#A3A500","#000080","#90EE90",
               "#EA8331","#00B0F6","#FF0000","#704241","#FFFF00","#C77CFF","#daa520","#adff2f","#8D2282","#E76BF3",
               "#39B600","#9590FF","#FF6A98","#FA62DB","#4169E1","#FF9999","#B22222","#FF62BC","#FFD032","#A2FF00",
               "#006400","#c0c0c0","#004A00","#FFB676","#87CEFA","#B34E00","#63E463","#6A5ACD")
filtered_data$cluster_name <- factor(filtered_data$cluster_name)

top_points <- filtered_data%>%
  group_by(cluster_name) %>%
  arrange(desc(total_sites)) %>%
  slice_head(n = 5000) %>%
  ungroup()

top_points <- filtered_data %>%
  filter(cluster_name %in% c(
    "Spermatogonial population",
    "Spermatocyte",
    "Acrosoma",
    "Round spermatids",
    "Elongating spermatids1",
#    "Elongating spermatids2",
    "Sperm"
  )) %>%
  group_by(cluster_name) %>%
  arrange(desc(total_sites)) %>%
  slice_head(n = 5000) %>%
  ungroup()

correct_order <- c(
  "Spermatogonial population",
  "Spermatocyte",
  "Acrosoma",
  "Round spermatids",
  "Elongating spermatids1",
  "Sperm"
)
top_points <- top_points %>%
  mutate(
    cluster_name = factor(cluster_name, levels = correct_order)
  ) %>%
  arrange(cluster_name)

top_points_hindbrain <- top_points[top_points$cluster_name=="hindbrain",]

set.seed(123)  # 保证可重复性
filtered_data_50_sub <- filtered_data_50 %>%
  group_by(cell_category) %>%
  slice_sample(n = 1247) 
library(MASS)
library(ggplot2)
{
  top_hindbrainpoints_2500 <- merged_data_hindbrain %>%
    filter(seurat_clusters %in% c("1", "2")) %>%
    group_by(seurat_clusters) %>%
    arrange(desc(total_sites)) %>%
    slice_head(n = 2000) %>%
    ungroup()
  h_x <- MASS::bandwidth.nrd(top_points$ratio_5mC)
  h_y <- MASS::bandwidth.nrd(top_points$ratio_5hmC)
  ggplot(top_hindbrainpoints_2500, aes(ratio_5mC, ratio_5hmC)) +
    geom_point(  # 仅高5mC的点
      color = "#4169E1",
      #aes(color = factor(seurat_clusters)), 
      alpha = 0.2,
      size=1
    ) +
  geom_density_2d(
    color = "black",
    #h=c(h_x,h_y),
    size = 0.1,
    n=300,
    bins = 10,
    alpha = 1
  ) +
  facet_wrap(~ seurat_clusters) +
  labs(
    title = "Methylation Contour Lines of hindbrain(2500)",
    x = "5mC Ratio",
    y = "5hmC Ratio"
  ) +
  theme_classic() +
  theme(
    panel.background = element_rect(fill = "white"),
    strip.background = element_blank())+
  coord_cartesian(ylim = c(0, 0.2)) 
  
  ggplot() +
    geom_point( data = top_hindbrainpoints_2500,  
      aes(x = ratio_5mC, y = ratio_5hmC, color = factor(seurat_clusters)),
      alpha = 0.2,
      size=1
    ) +
    geom_density_2d(
      data=top_hindbrainpoints,
      aes(x = ratio_5mC, y = ratio_5hmC),
      color = "black",
      #h=c(h_x,h_y),
      size = 0.1,
      n=300,
      bins = 10,  
      alpha = 1
    ) +
    labs(
      title = "Methylation Contour Lines of hindbrain",
      x = "5mC Ratio",
      y = "5hmC Ratio"
    ) +
    scale_color_manual(
      values = c("1" = "#FF6B6B", "2" = "#4ECDC4", "NA" = "gray80"),
      na.value = "gray80",  
      name = "Seurat Cluster", 
      labels = c("C1", "C2", "Unassigned")) +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = "white"), 
      strip.background = element_blank())+
    coord_cartesian(ylim = c(0, 0.2)) 
  
  
}

