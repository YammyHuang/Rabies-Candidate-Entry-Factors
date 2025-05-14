#20230909尝试整合mPFC rAAV与V1 RV数据

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
adult.inte <- FindNeighbors(adult.inte, reduction = "pca", dims = 1:10)
adult.inte <- FindClusters(adult.inte, resolution = 1)
adult.inte <- RunUMAP(adult.inte, reduction = "pca", dims = 1:10)
adult.inte <- RunTSNE(adult.inte, reduction = "pca", dims = 1:10)
adult.inte <-rAAV_RV.inte
adult.inte$seurat_clusters
VlnPlot(adult.inte,
        features = c("Snap25",'Gad1',"Sox6","Col23a1","Pvalb","Sst","Prox1","Vip","Lamp5","Slc17a7",
                     "Rorb","Foxp2","Bcl6","Ctss","Slc1a3","Mog","Cux2","Rorb","Deptor","Etv1","Pou3f1",'Oprk1','Foxp2',"Syt6","Tshz2"),stack = TRUE, flip = TRUE,  assay = 'RNA') + 
  NoLegend()
VlnPlot(adult.inte,
        features = c("Snap25",'Gad1','Gad2',"Bcl11b","Sox6",
                     "Col23a1","Pvalb","Sst","Prox1","Vip","Lamp5","Slc17a7",
                     "Rorb","Foxp2","Bcl6","Ctss","Slc1a3","Mog","Cux2",
                     "Rorb","Deptor","Etv1","Pou3f1",'Oprk1','Foxp2',"Syt6","Tshz2"),
        stack = TRUE, flip = TRUE,  assay = 'RNA') + NoLegend()

VlnPlot(adult.inte,
        features = c("Snap25",'Gad1','Gad2',"Sox6",
                     "Pvalb","Sst","Prox1","Vip","Lamp5","Aldoc", "Slc1a3", "Aqp4","Slc17a6","Slc17a7",
                     'Calb1','Otof','Hap1','Cux2','Rorb','Tnnc1','Rspo1','Bdnf',
                     'Ptn','Sema3d','Cpne7','Fstl5','Slc24a2','Oprk1','Penk',"Col23a1",
                     'Syt6','Pcp4','Pou3f1','S100b','Etv1','Scn4b','Adamts2','Dlk1','Npr3',
                     'Tshz2','Cbln2','Grp',"Pdgfra", "C1ql1", "Olig2", "Olig1","Fcrls", "Trem2"),
        stack = TRUE, flip = TRUE,  assay = 'RNA') + NoLegend()
level_maintype<-c( "Pvalb","Sst","Vip/Lamp5","Astro","Oligo","Microglia","L2/3IT","L4/5IT","L6IT","L5NP","L6CT","L5ET")
Idents(adult.inte) <- factor(Idents(adult.inte),levels=level_maintype)
levels(Idents(adult.inte))
adult.inte$Maintype <- as.character(adult.inte$seurat_clusters)
adult.inte$Maintype[which(adult.inte$Maintype %in% c(9))] <- "Pvalb"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(11))] <- "Sst"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(13))] <- "Vip/Lamp5"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(0,1,6,9))] <- "L2/3IT"#cux2
adult.inte$Maintype[which(adult.inte$Maintype %in% c(5,21,22))] <- "L4/5IT"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(3,18))] <- "L6IT"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(2,22))] <- "L6CT"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(12,17))] <- "L5ET"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(8,15))] <- "L5NP"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(10,14,19))] <- "Microglia"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(4,16,17))] <- "Astro"
adult.inte$Maintype[which(adult.inte$Maintype %in% c(7,20))] <- "Oligo"
Idents(adult.inte) <- 'Maintype'
VlnPlot(adult.inte,
        features = c('Gad1',"Sox6","Pvalb","Sst","Prox1","Vip","Lamp5","Slc17a7","Col23a1",
                     "Rorb","Foxp2","Bcl6","Ctss","Slc1a3","Mog"),group.by = "Maintype",stack = TRUE, flip = TRUE, fill.by="ident", assay = 'RNA') + 
  NoLegend()

DimPlot(adult.inte, split.by = 'sample', reduction = 'tsne', label = T, ncol = 2) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), 
        axis.text = element_blank(), plot.title = element_text(size = 30)) +
  labs(x='', y='')
saveRDS(adult.inte,"H:/Project1_RV Receptor Projection/FIG1.皮层单细胞RV rAAV感染数据分析/rAAV_RV.inte2.RDS")
 