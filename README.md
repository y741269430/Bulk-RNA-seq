# Bluk-RNA-seq

## 目录    
先挖个坑    

## 加载以下包，函数    
```r
library(DESeq2)
library(factoextra)
library(ggplot2)
library(pheatmap)

myNormal <- function(x, grouplist){
  colData <- data.frame(row.names = colnames(x), grouplist)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = x,
                                        colData = colData,
                                        design = ~ grouplist)
  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds)  
  nor <- DESeq2::counts(dds, normalized = T)
  nor2 <- as.data.frame(nor)
}
myPCA <- function(x, grouplist, label_option = "ind") {
  factoextra::fviz_pca_ind(prcomp(t(x),scale. = T), 
                           repel = TRUE, pointsize = 5, 
                           palette = "jco", label = ifelse(label_option == "none", "", "ind"), mean.point = F, 
                           col.ind = grouplist) + 
    ggplot2::coord_fixed(1) + 
    ggplot2::labs(title = element_blank())+ 
    ggplot2::theme_bw()+
    ggplot2::theme(text = element_text(family = "serif"),
                   plot.title = element_text(hjust = 0.5),
                   axis.title = element_text(size = 15),
                   axis.text = element_text(size = 15),
                   legend.title = element_blank(), 
                   legend.text = element_text(face = "bold", size = 15),
                   legend.position = "bottom") + 
    xlim(-250, 250) + ylim(-250, 250)
}
myHeat <- function(x, grouplist, show_rownames = F, 
                   angle_col = 45, cluster_cols = F, cluster_rows = T, main = NA){
  # 输入一个标准化的矩阵，进行热图绘制，以grouplist为准，默认关闭行名，列名45度
  # 输出是一个热图
  # 保存热图一定要加上pheatmap的包名: pheatmap::
  
  group <- data.frame(grouplist)
  rownames(group) <- colnames(x)
  plot <- pheatmap::pheatmap(log10(x + 1), 
                             scale = "row", 
                             clustering_distance_rows = "correlation",
                             #color = rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(256)),
                             color = colorRampPalette(c("navy", "white", "firebrick3"))(256),
                             # color = colorRampPalette(c('#03045E', '#0077B6', '#00B4D8', '#90E0EF', '#CAF0F8',
                             #                            '#FAE0E4', '#F9BEC7', '#FF99AC', '#FF7096', '#FF477E'))(256),
                             fontsize = 8, 
                             fontfamily="serif",
                             display_numbers = F,
                             border_color = NA,
                             gaps_row = F,
                             cluster_cols = cluster_cols, 
                             cluster_rows = cluster_rows,
                             treeheight_row = 50,treeheight_col = 50,
                             cellwidth = NA, 
                             cellheight = NA,
                             legend = T,
                             show_rownames = show_rownames, 
                             show_colnames = T,
                             annotation_legend = T,
                             annotation_names_col = F,
                             annotation_col = group, 
                             angle_col = angle_col, 
                             main = main,
                             silent = T) 
  plot_grid(plot$gtable)
}
```

## 合并count文件变成一个矩阵    
```r
path <- '文件路径'

fileNames <- dir(path, pattern = ".count$")
filePath <- sapply(fileNames, function(x){paste(path, x, sep = '')})

data <- lapply(filePath, function(x){
  read.table(x, header = F, sep = "\t", colClasses = c("character", "numeric"))})

rawdata <- do.call(cbind, lapply(data, function(x){x <- x[,2]}))

data <- cbind(data.frame('ENSEMBL' = substr(data[[1]]$V1, 1, 18), rawdata))

rownames(data) <- data$ENSEMBL

# 重命名列（样本）
colnames(data)[-1] <- c('CTRL_1', 'CTRL_2', 'CTRL_3', 'Treat_1', 'Treat_2', 'Treat_3')

data <- data[,-1]

write.csv(data, file = paste0(path, 'rawcount.csv))
```

## 进行质控（PCA+热图+皮尔森相关性分析）

```r
grouplist <- c(rep('CTRL', 3), rep('Treat', 3))
nor <- myNormal(data, grouplist)

# 画PCA图
plot_pca <- myPCA(nor, grouplist)


# 画热图
nor_data_selcet <- nor[sample(nrow(nor), size = 2000, replace = FALSE), ]

plot_heat <- myHeat(nor_data_selcet, grouplist, cluster_cols = T, cluster_rows = T)


# 画相关性图
cor_data <- round(cor(nor), 3)

plot_cor <- as.ggplot(pheatmap(cor_data, c("#6D9EC1", "white", "#E46726"),
                               name = "Pearson correlation",
                               fontsize = 10, display_numbers = T))

p <- plot_pca+plot_heat+plot_cor

ggplot2::ggsave(paste0(path, 'PCA_Heatmap_Pearson_cor.pdf'), plot = p
                height = 6, width = 18, dpi = 300, limitsize = FALSE)
```



