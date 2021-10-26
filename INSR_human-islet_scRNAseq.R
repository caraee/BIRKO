library("plyr")
library("tidyverse")
library("here")
library("Seurat")
library("RColorBrewer")
#devtools::install_github('satijalab/seurat-data')
library("SeuratData")

options(future.globals.maxSize = 4000 * 1024^2)
#InstallData("panc8")
data("panc8")
pancreas.list <- SplitObject(panc8, split.by = "tech")
pancreas.list <- pancreas.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")]
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- SCTransform(pancreas.list[[i]], verbose = FALSE)
}

pancreas.features <- SelectIntegrationFeatures(object.list = pancreas.list, nfeatures = 3000)
pancreas.list <- PrepSCTIntegration(object.list = pancreas.list, anchor.features = pancreas.features, 
                                    verbose = FALSE)

pancreas.anchors <- FindIntegrationAnchors(object.list = pancreas.list, normalization.method = "SCT", 
                                           anchor.features = pancreas.features, verbose = FALSE)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

pancreas.integrated <- RunPCA(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30)
plots <- DimPlot(pancreas.integrated, group.by = c("tech", "celltype"))
plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
                                                                     override.aes = list(size = 3)))

DefaultAssay(pancreas.integrated) <- "RNA"
# Normalize RNA data for visualization purposes
pbmc.integrated <- NormalizeData(pancreas.integrated, verbose = FALSE)

cluster_cols2<-c("#343d00",
                 "#834efd",
                 "#d7ff44",
                 "#011b9e",
                 "#9acb00",
                 "#e953ff",
                 "#00df71",
                 "#ac00ba",
                 "#438500",
                 "#dc0086",
                 "#008e5a",
                 "#ff0072",
                 "#02edf7",
                 "#ff3401",
                 "#70a0ff",
                 "#ffdb44",
                 "#001756",
                 "#deffb0",
                 "#3b003a",
                 "#91ffe3",
                 "#790034",
                 "#a5dcff",
                 "#892f00",
                 "#a78cff",
                 "#827a00",
                 "#ff99ce",
                 "#002d00",
                 "#ff5f57",
                 "#029cb4",
                 "#ff9a62",
                 "#000f14",
                 "#ffecc1",
                 "#3a0019",
                 "#f9d2ff",
                 "#00536d")

Idents(pancreas.integrated)<-"celltype"

unique(pancreas.integrated$celltype)

COI<-c("gamma","acinar","alpha","delta","beta","ductal","endothelial","activated_stellate",
       "schwann","mast","macrophage","epsilon","quiescent_stellate")

png(filename = "./Figures/INSR_islets_ridge_all.png",width=20,height=15,units="cm",res=300)
RidgePlot(object = pancreas.integrated, cols = cluster_cols2, sort = "increasing",
          features = c("INSR"), 
          idents = COI,
          log = TRUE)+NoLegend()
dev.off()

png(filename = "./Figures/INSR_islets_ridge_islet_endothelial_integrated.png",width=20,height=15,units="cm",res=300)
RidgePlot(object = pancreas.integrated, cols = cluster_cols2, sort = "increasing",
          features = c("INSR"), 
          idents = COI[c(1:7,12)])+NoLegend() & theme( plot.title = element_text( face = "italic") )
dev.off()

citation("panc8.SeuratData")
citation("Seurat")

features.plot <- "INSR"
DefaultAssay(pancreas.integrated) <- "integrated"
plot1 <- FeaturePlot(pancreas.integrated, features.plot,
                     cells = which(Idents(pancreas.integrated)%in%COI[c(1:7,12)]),
                     reduction = "umap",label=T,repel = T,order = T) &
  theme( plot.title = element_text( face = "italic") )

plot2 <- DimPlot(islets.integrated, cols=cluster_cols2, reduction = "umap",label=F,pt.size = 0.2)

png(filename = "./Figures/INSR_umap.png",width=30,height=15,units="cm",res=300)
plot1 
dev.off()

summary(pancreas.integrated@assays$RNA@counts["INSR",which(Idents(pancreas.integrated)%in%"beta")])
sum(pancreas.integrated@assays$RNA@counts["INSR",which(Idents(pancreas.integrated)%in%"beta")]>1)

table(Idents(pancreas.integrated))

save(pancreas.integrated,file="for_sis_panc.RData")
