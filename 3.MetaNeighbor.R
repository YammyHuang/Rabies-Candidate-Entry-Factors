adult.inte$maintype_sample <- paste(adult.inte$Maintype, adult.inte$sample, sep = "|")
library(dplyr)
table(adult.inte$sample)
adult.inte$sample <- recode(adult.inte$sample, 
                            "rAAV_MOCK" = "control_PFC", 
                            "RV_MOCK" = "Control_V1")
Idents(adult.inte) <- 'maintype_sample'

#build the tree
adult.inte<- BuildClusterTree(object = adult.inte, dims=1:10)
phy <- Tool(object = adult.inte, slot = 'BuildClusterTree')
p5<-plot(phy)
p5
unique(adult.inte$maintype_sample)

#MetaNeighbor
library(SingleCellExperiment)
library(MetaNeighbor)

# 创建 colData
new_colData <- data.frame(
  study_id = adult.inte$sample,
  cell_type = adult.inte$Maintype # 确保 row.names 对齐
)
new_colData$study_id <- as.character(new_colData$study_id)
new_colData$study_id[which(new_colData$study_id %in% c("Control_V1","RV_infected"))] <- "VISp"
new_colData$study_id[which(new_colData$study_id %in% c("control_PFC","rAAV_infected"))] <- "PFC"

dat <- SingleCellExperiment(adult.inte@assays$integrated@data,
  colData = new_colData
)
var_genes <- VariableFeatures(adult.inte)  # 替换为实际的高变基因获取方式
# 运行 MetaNeighbor
celltype_NV = MetaNeighborUS(
  var_genes = var_genes, 
  dat = dat, # 确保 dat 是矩阵
  study_id = dat$study_id,
  cell_type = dat$cell_type,
  fast_version = TRUE
)

library(pheatmap)
pheatmap(celltype_NV, cluster_rows=T, cluster_cols=F, display_numbers=T)

cor_inter<-celltype_NV[grep('VISp|PFC',rownames(celltype_NV),value = T),
                       grep('VISp|PFC',rownames(celltype_NV),value = T)]
p<-pheatmap(cor_inter, cluster_rows=T, cluster_cols=T, display_numbers=T,fontsize_number = 5,silent = T)
cluster_order <- rownames(cor_inter)[p$tree_row$order] 

pheatmap(cor_inter, cluster_rows=T, cluster_cols=T, display_numbers=T, 
         fontsize_number=5, fontsize_row=8, fontsize_col=8)
p5<-pheatmap(cor_inter[cluster_order, cluster_order], 
         cluster_rows=F, cluster_cols=F, display_numbers=T, fontsize_number=5)

cluster_order<-c("PFC|Pvalb" ,"VISp|Pvalb","PFC|Vip/Lamp5","VISp|Vip/Lamp5","PFC|Sst" ,"VISp|Sst",
                 "PFC|Microglia","VISp|Microglia","PFC|Oligo"  ,"VISp|Oligo","PFC|Astro","VISp|Astro",
                  "VISp|L2/3IT" ,"PFC|L2/3IT","VISp|L4/5IT", "PFC|L4/5IT","PFC|L6IT","VISp|L6IT" , "VISp|L5ET","PFC|L5ET","PFC|L6CT", "VISp|L6CT","PFC|L5NP","VISp|L5NP") 
ggsave("H:/Project1_RV Receptor Projection/FIG1.皮层单细胞RV rAAV感染数据分析/p5_MetaNeighbor.pdf", 
       plot = p5, width =8, height = 8, dpi = 300)




# 提取所有共有细胞类型（PFC和VISp共有的类型）
common_types <- sub("PFC\\|", "", grep("PFC\\|", rownames(cor_inter), value = TRUE))
common_types <- common_types[common_types %in% sub("VISp\\|", "", grep("VISp\\|", rownames(cor_inter), value = TRUE))]

# 提取PFC与VISp同类型细胞的相似性得分（单向即可，因矩阵对称）
scores <- sapply(common_types, function(type) {
  cor_inter[paste0("PFC|", type), paste0("VISp|", type)]
})

# 输出所有配对得分
print(data.frame(CellType = common_types, AUROC = scores))

# 计算平均相似性得分（假设理想情况下应接近1）
mean_score <- mean(scores)
cat(sprintf("\nMean AUROC between matched cell types: %.2f\n", mean_score))

# 显著性检验：判断得分是否显著高于随机预期（假设随机预期为0.5）
t_test_result <- t.test(scores, mu = 0.5, alternative = "greater")
cat(sprintf("t-test p-value vs. random expectation (0.5): %.1e\n", t_test_result$p.value))
