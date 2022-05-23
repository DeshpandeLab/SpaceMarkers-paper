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

pdf("Expression_plots_PDAC.pdf",width = 8.5,height = 5)
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
        ylim(max(bcs_merge %>% 
                select(height)),0)+
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

region <- SpaceMarkers$hotspotRegions[,6]
region <- ifelse(!is.na(region) & !is.na(SpaceMarkers$hotspotRegions[,9]),"Interacting",ifelse(!is.na(region),region,SpaceMarkers$hotspotRegions[,9]))
region <- factor(region)

pdf("Expression_boxplots_PDAC.pdf",width = 5,height = 5)
plist <- list()
mplot2 <- mplot[!is.na(region),]
for (ii in 1:length(geneList)){
plist[[ii]]<- mplot2 %>% ggplot( aes_string(x="region", y=geneList[ii], fill="region")) +
    geom_boxplot() +
        scale_fill_viridis(discrete = TRUE, alpha=0.6) +
        geom_jitter(color="black", size=0.4, alpha=0.9) +
        theme_ipsum() +
        theme(
             legend.position="none",
             plot.title = element_text(size=11)
         ) +
     ggtitle(paste0(geneList[ii]," Expression (Log)")) +
    xlab("")
}
n <- length(plist)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(plist, ncol=4))
dev.off()
