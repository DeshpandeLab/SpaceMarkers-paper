bcs_merge <- readRDS("~/FertigLab/Visium/VisiumCogaps/Results/01_s182967barcodeInfo.rds")
matrixPath <- "~/FertigLab/Visium/VisiumCogaps/Data/182967_1168993F_2_raw_feature_bc_matrix.h5"
matrix <- as.data.frame(t(as.matrix(Read10X_h5(matrixPath))))
matrix <- matrix[bcs_merge$barcode[bcs_merge$sample=="s182967"],]
sample_names <- unique(bcs_merge$sample)[2]

load("/Users/atul/Library/CloudStorage/OneDrive-JohnsHopkins/Incoming/images_tibble_flipped.rdata")
plots <- list()
geneList <- c("ERBB2","ESR1","INSRR","PIK3R3","FGFR1","IGFBP3","TFAP2E","FGFR4","TGFB2")
geneList <- c("CD4")
geneList<- c("VIM", "JUND", "MYH9", "CD74", "GNAS", "MSLN", "LGALS9", "ACTB")
myPalette <- colorRampPalette(rev(brewer.pal(n = 11,"Spectral")))

pdf("Expression_plots_Breast_Cancer_2.pdf",width = 7,height = 5)
for (ii in 1:length(geneList)) {
    if (max(matrix[,geneList[ii]])>2)
    {
        plots <- bcs_merge %>% 
        bind_cols(as.data.table(matrix)[, geneList[ii], with=FALSE]) %>% 
        ggplot(aes_string(x="imagecol",y="imagerow",fill=geneList[ii])) +
        geom_spatial(data=images_tibble[1,], aes(grob=grob), x=0.5, y=0.5)+
        geom_point(shape = 21, colour = "black", size = 1, alpha = 0.9, stroke = 0.1)+
        coord_cartesian(expand=FALSE)+
        scale_fill_gradientn(colours = myPalette(100))+
        xlim(0,max(bcs_merge %>% 
                select(width)))+
        ylim(0,max(bcs_merge %>% 
                select(height)))+
        xlab("") +
        ylab("") +
        ggtitle(geneList[ii])+
        theme_set(theme_bw(base_size = 10))+
        theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.text = element_blank(),
            axis.ticks = element_blank())
        print(plots)
    }
}
dev.off()
