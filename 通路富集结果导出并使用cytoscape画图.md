# 通路富集结果导出并使用cytoscape画图
## 目录 ####
- 1.数据加载与基因ID转换    
- 2.使用clusterProfiler进行基因ID转换
- 3.合并转换结果并去除重复的Ensembl ID
- 4.差异基因分组
- 5.GO富集分析
- 6.关键通路基因提取
- 7.Cytoscape 网络数据构建
- 8.Cytoscape 网络图绘制  

## 1.数据加载与基因ID转换    
```r
load('test.RData')  # 加载原始数据
```
## 2.使用clusterProfiler进行基因ID转换（Ensembl转Symbol/EntrezID）
```r
gene_info <- clusterProfiler::bitr(data$ENSEMBL, 
                                   fromType = 'ENSEMBL', 
                                   toType = c('SYMBOL', 'ENTREZID','GENENAME'), 
                                   OrgDb = org.Mm.eg.db)
```
## 3.合并转换结果并去除重复的Ensembl ID
```r
data <- merge(gene_info, data, by = "ENSEMBL") |> 
  subset(!duplicated(ENSEMBL))
```

## 4.差异基因分组    
```r
deg <- list(
  UP   = subset(data, log2FC >= 0),  # 上调基因组
  DOWN = subset(data, log2FC <= 0)    # 下调基因组
)

# 提取ENTREZID列表用于富集分析
entrez_list <- lapply(deg, function(x){ x <- x$ENTREZID })
```
---

## 5.GO富集分析    
### 5.1 使用compareCluster进行批量GO分析)

```r
go_res <- clusterProfiler::compareCluster(
  geneClusters = entrez_list,
  fun = "enrichGO",
  ont = "BP",
  OrgDb = org.Mm.eg.db,
  keyType = 'ENTREZID',
  readable = TRUE,
  pAdjustMethod = "BH"  # 添加FDR校正
)

# 分离上下调通路结果
go_up <- subset(go_res@compareClusterResult, Cluster == "UP")
go_down <- subset(go_res@compareClusterResult, Cluster == "DOWN")
```

### 5.2 使用enrichGO进行GO分析,没有区分上下调基因)
```r
gp_res_all <- clusterProfiler::enrichGO(
  data$ENTREZID, 
  ont = "BP",
  OrgDb = org.Mm.eg.db,
  keyType = 'ENTREZID',
  readable = TRUE,
  pAdjustMethod = "BH"  # 添加FDR校正
)
```


## 6.关键通路基因提取
### 如果你是结果1得到的GO通路，则选择1，如果是2则选择2
```r
go_df <- go_up                  # 结果1的GO通路
go_df <- gp_res_all@result      # 结果2的GO通路
```

```r
go_terms <- sample(go_df$Description, 3)  # 随机选取3个GO term
go_terms <- c('term_a', 'term_b')   # 或选取2个GO term
```
### 获取目标term的基因列表
```r
target_genes <- lapply(go_terms, \(term) {
  ids <- unlist(strsplit(subset(go_df, Description == term)$geneID, "/"))
  merge(data.frame(SYMBOL = ids), data, by = "SYMBOL")
}) |> setNames(gsub(" ", "_", go_terms))  # 名称格式化
```

---

## 7.Cytoscape 网络数据构建
### 构建边列表（包含log2FC和通路节点）
```r
edge_data <- lapply(seq_along(target_genes), \(i) {
  transform(target_genes[[i]], node = go_terms[i])
}) |> do.call(what = rbind)

write.xlsx(edge_data, 'cyto边文件.xlsx')
```

### 构建节点列表（基因节点+通路节点）
```r
node_data <- rbind(
  edge_data[, c("SYMBOL", "log2FC")],
  data.frame(SYMBOL = unique(edge_data$node), 
             log2FC = 0)  # 通路节点log2FC设为0
)

write.xlsx(node_data, 'cyto节点文件.xlsx')
```

---
## 8.Cytoscape 网络图绘制  
- [下载软件](https://cytoscape.org/download.html)
- [操作手册](http://manual.cytoscape.org/en/stable/)    
步骤：123导入边文件，45导入节点文件    
<img src="https://github.com/y741269430/Bluk-RNA-seq/blob/main/figure/cyto1.png" width="600" />
<img src="https://github.com/y741269430/Bluk-RNA-seq/blob/main/figure/cyto2.png" width="600" />
<img src="https://github.com/y741269430/Bluk-RNA-seq/blob/main/figure/cyto3.png" width="600" />
