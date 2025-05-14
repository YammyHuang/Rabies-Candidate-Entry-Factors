
marker_genes<-c('Gad1',"Sox6","Pvalb","Sst","Prox1","Vip"
               ,"Ctss","Slc1a3","Mog","Slc17a7","Cux2","Deptor","Etv1","Foxp2","Syt6","Bcl6","Pou3f1"
               ,"Bcl6","Col23a1","Rorb",'Cdh13','Nnat','Rab3c',"Rorm",'Pcp4','Kcnip1','Tshz2','Penk','Calb1','Npr3')
# 计算每个细胞群和区域的平均表达值
avg_exp <- AverageExpression(rAAV_RV.inte, group.by = c("region", "Maintype"), return.seurat = TRUE)

# 提取 RNA 表达数据
data_matrix <- as.matrix(avg_exp@assays$RNA@data)

# 筛选 marker 基因的表达数据
marker_matrix <- data_matrix[rownames(data_matrix) %in% marker_genes, ]
# 计算基于 marker 基因的相关性矩阵
cor_matrix <- cor(marker_matrix, method = "pearson")
cor_matrix <- cor(marker_matrix, method = "spearman")
library(pheatmap)
pheatmap(cor_matrix, clustering_method = "ward.D2")

library(MetaNeighbor)
common_genes <- rownames(arabidopsis)[rownames(arabidopsis) %in% rownames(rice)]