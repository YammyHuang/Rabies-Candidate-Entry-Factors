library(Seurat)
DimPlot(rAAV_RV.inte,)
table(rAAV_RV.inte$sample)
FeaturePlot(rAAV_RV.inte,
            reduction = 'tsne',
            c('rAAV_infected','RV_infected','Syt6','Tshz2'),
            ncol = 2,
            cols=c("lightgrey", "red"),
            min.cutoff='q1')
Idents(rAAV_RV.inte) <- rAAV_RV.inte$sample  # 设定细胞身份

DimPlot(rAAV_RV.inte, reduction = "tsne", group.by = "sample",
        cols = c("#00bfc4", "lightgrey","red", "grey"), pt.size = 0.5) +
  ggtitle("t-SNE Visualization of Infected Cells")

Idents(rAAV_RV.inte) <- rAAV_RV.inte$Maintype  # 设定细胞身份
DimPlot(rAAV_RV.inte, reduction = "umap", pt.size = 0.5) +
  ggtitle("t-SNE Visualization of Infected Cells")
table(rAAV_RV.inte$Maintype,rAAV_RV.inte$sample)
