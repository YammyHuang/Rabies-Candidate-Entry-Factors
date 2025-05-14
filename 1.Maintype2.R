#20230909尝试整合mPFC rAAV与V1 RV数据
install.packages("ape")
RV_Infected<-RV_infected_mock_intergreted[,RV_infected_mock_intergreted$sample%in%"RV_infected"]
RV_MOCK<-RV_infected_mock_intergreted[,RV_infected_mock_intergreted$sample%in%"MOCK"]

rAAV_Infected<-all.Adult[,all.Adult$BC_label%in%"Barcoded"]
rAAV_MOCK<-all.Adult[,all.Adult$BC_label%in%"Unbarcoded"]

VIRUS_LIST<-list(RV_Infected,RV_MOCK,rAAV_Infected,rAAV_MOCK)
virus_sample_name<- c('RV_infected', 'RV_MOCK',"rAAV_infected","rAAV_MOCK")

for (i in 1:length(VIRUS_LIST)){
  VIRUS_LIST[[i]]@meta.data$sample <- virus_sample_name[[i]]}

p <- list()
for (i in 1:4){
  p[[i]] <- DimPlot(VIRUS_LIST[[i]], reduction = 'tsne') + 
    labs(title = virus_sample_name[i]) +
    guides(colour=guide_legend(ncol=2, override.aes = list(size=2)))
}
plot_grid(plotlist = p, ncol =2)

#整合
features <- SelectIntegrationFeatures(object.list = VIRUS_LIST)
adult.anchors <- FindIntegrationAnchors(object.list = VIRUS_LIST, 
                                        anchor.features = features)
adult.inte <- IntegrateData(anchorset = adult.anchors)
adult.inte <- ScaleData(adult.inte, verbose = FALSE)
adult.inte <- RunPCA(adult.inte, npcs = 10, verbose = FALSE)
adult.inte <- FindNeighbors(adult.inte, reduction = "pca", dims = 1:20)
adult.inte <- FindClusters(adult.inte, resolution = 1)
adult.inte <- RunUMAP(adult.inte, reduction = "pca", dims = 1:10)
adult.inte <- RunTSNE(adult.inte, reduction = "pca", dims = 1:10)

adult.inte$seurat_clusters

VlnPlot(adult.inte,
        features = c("Snap25",'Gad1','Gad2',"Sox6",
                     "Pvalb","Sst","Prox1","Vip",
                     "Aldoc", "Slc1a3", "Aqp4", "Olig2", "Olig1","Pdgfra", "C1ql1","Fcrls", "Trem2",
                     "Slc17a7",'Calb1','Cux2','Rorb','Bdnf',
                     'Ptn',"Col23a1",
                     'Tshz2','Cbln2','Grp',"Syt6",
                     'Pou3f1','Etv1','Adamts2','Dlk1','Npr3'),
        stack = TRUE, flip = TRUE,  assay = 'RNA') + NoLegend()
level_maintype<-c( "Pvalb","Sst","Vip/Lamp5","Astro","Oligo","Microglia","L2/3IT","L4/5IT","L6IT","L5NP","L6CT","L5ET")
Idents(adult.inte) <- factor(Idents(adult.inte),levels=level_maintype)
levels(Idents(adult.inte))

Idents(adult.inte) <- 'seurat_clusters'
Idents(adult.inte) <- 'Maintype'
adult.inte<- BuildClusterTree(object = adult.inte, dims=1:10)
phy <- Tool(object = adult.inte, slot = 'BuildClusterTree')
plot(phy)

adult.inte$Maintype <- as.character(adult.inte$seurat_clusters)
adult.inte$Maintype[which(adult.inte$Maintype %in% c(8))] <- "Pvalb"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(9))] <- "Sst"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(13,27))] <- "Vip/Lamp5"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(2,14))] <- "L2/3IT"#cux2
adult.inte$Maintype[which(adult.inte$Maintype %in% c(0,6))] <- "L4/5IT"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(3,4,23,25,19))] <- "L6IT"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(1))] <- "L6CT"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(11))] <- "L5ET"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(16))] <- "L5NP"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(10,15,21))] <- "Microglia"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(5,17,24,7))] <- "Astro"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(12,18,20,22,26))] <- "Oligo"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(7))] <- "7"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(19))] <- "19"
Idents(adult.inte) <- 'Maintype'
VlnPlot(adult.inte,
        features = c('Gad1',"Sox6","Pvalb","Sst","Prox1","Vip","Lamp5","Slc17a7","Col23a1",
                     "Rorb","Foxp2","Bcl6","Ctss","Slc1a3","Mog"),group.by = "Maintype",stack = TRUE, flip = TRUE, fill.by="ident", assay = 'RNA') + 
  NoLegend()
DimPlot(adult.inte,  reduction = 'umap', label = T, ncol = 2)
DimPlot(adult.inte, split.by = 'sample', reduction = 'tsne', label = T, ncol = 2) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.text = element_blank(), plot.title = element_text(size = 30)) +
  labs(x='', y='')
saveRDS(adult.inte,"H:/Project1_RV Receptor Projection/FIG1.皮层单细胞RV rAAV感染数据分析/rAAV_RV.inte2.RDS")


table(SLQ_IPC$predicted.id)
SLQ_IPC$`ACB-I`
SLQ_IPC$`BLA-I`