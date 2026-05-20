library(Seurat)
library(SeuratData)
library(SingleCellExperiment)
library(scDesign3)
library(scales)
library(ggplot2)
library(cowplot)

brain <- LoadData("stxBrain", type = "anterior1")

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

brain

# transform into SingleCellExperiment object
sce <- as.SingleCellExperiment(brain, assay = "SCT")

spatial_coords <- GetTissueCoordinates(brain, image = "anterior1")

# 3. Adicionar as coordenadas diretamente no colData do SCE
# Garantindo o alinhamento correto pelos nomes dos spots/barcodes
colData(sce)$spatial1 <- spatial_coords[colnames(sce), "x"]
colData(sce)$spatial2 <- spatial_coords[colnames(sce), "y"]

# remove MT genes
mt_idx<- grep("mt-",rownames(sce))
if(length(mt_idx)!=0){
    sce <- sce[-mt_idx,]
}

sce
sce@colData

# we first select reference genes based on moran'I to
# reduce fitting time for scDesign3
num_genes <- 200
loc <- colData(sce)[, c("spatial1","spatial2")]
features <- FindSpatiallyVariableFeatures(counts(sce), 
                                          spatial.location = loc, 
                                          selection.method = "moransi", 
                                          nfeatures = num_genes)

top.features <- features[order(features$p.value),]

top.features <- rownames(top.features[1:num_genes,])

de_idx <- which(rownames(sce) %in% top.features)
sce <- sce[de_idx, ]



plot_exp <- function(sce, gene, pt_size=1){
    df_loc <- colData(sce)[, c("spatial1","spatial2")]
    df_exp <- as.data.frame(counts(sce)[gene, ])
    colnames(df_exp) <- c('exp')
    df_exp$exp <- rescale(log1p(df_exp$exp))
    
    df <- cbind(df_exp, df_loc)
    
    p <- ggplot(data = df, aes(x = spatial1, y = spatial2)) +
    geom_point(aes(color = exp), size = pt_size) +
    scale_colour_gradientn(colors = viridis_pal(option = "magma")(10), limits=c(0, 1)) +
    theme_cowplot() +
    theme(axis.text = element_blank(), axis.ticks = element_blank()) +
    ggtitle(gene)
    
    return(p)
}

options(repr.plot.height = 4, repr.plot.width = 5)
par(mfrow=(c(2,2)))
for(gene in rownames(sce)[1:4]){
    p <- plot_exp(sce, gene = gene, pt_size = 1.2)
    print(p)
}


# run scDesign3

set.seed(123)
