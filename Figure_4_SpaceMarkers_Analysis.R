## author: Atul Deshpande
## email: adeshpande@jhu.edu
rm(list = ls())
setwd('.')
source('plottingScripts/patternSpotter.R')
## Specify data folder paths here: 
## Expected structure: parent_folder
##                              |____ VisiumDir (10x output format)
##                              |           |____ patient_id1
##                              |           |____ patient_id2
##                              |           .   .   .   .   .
##                              |____ CoGAPS_Analysis
##                                          |____ patient_id1
##                                          |           |____ cogapsFilePattern1
##                                          |           |____ cogapsFilePattern2
##                                          |             .   .   .   .   .
##                                          |____ patient_id2
##                                                      |____ cogapsFilePattern3
##                                                      |____ cogapsFilePattern4
##                                                        .   .   .   .   .
patient_id1 <- 'BreastCancer'
visiumDir <- "VisiumData/"
cogapsDir <- "CoGAPS_Analysis/"
cogapsFilePattern <- "182967_1168993F_2_CogapsResult_20.rds"
rngtools::RNGseed(123)

## Set these parameters
# SpaceMarkersMode: defaut mode is "residual". You can also set "DE" mode for Differential Expression mode.
SpaceMarkersMode = "residual"  
# SpaceMarkersRefPattern is the pattern whose "interaction" with every other pattern we want to study. If refPattern is not explicitly assigned, the code assumes Pattern_1 to be refPattern. Here the immune pattern Pattern_1 is the reference.
SpaceMarkersRefPattern = "Pattern_1" 

## Loading data
pVisiumPath <- paste0(visiumDir,patient_id1)
fullMat <- load10XExpr(pVisiumPath)
good_gene_threshold <- 3
goodGenes <- rownames(fullMat)[apply(fullMat,1,function(x) sum(x>0)>=good_gene_threshold)]
fullMat <- fullMat[goodGenes,]
spCoords <- load10XCoords(pVisiumPath)
rownames(spCoords) <- spCoords$barcode
pCoGAPSPath <- paste0(cogapsDir,patient_id1)

## Matching formats
cogapsFilePath <- dir(pCoGAPSPath,cogapsFilePattern,full.names = T)
CoGAPS_Result <- readRDS(cogapsFilePath)
features <- intersect(rownames(fullMat),rownames(CoGAPS_Result@featureLoadings))
barcodes <- intersect(colnames(fullMat),rownames(CoGAPS_Result@sampleFactors))
fullMat <- fullMat[features,barcodes]
cgMat <- CoGAPS_Result@featureLoadings[features,] %*% t(CoGAPS_Result@sampleFactors[barcodes,])
spCoords <- spCoords[barcodes,]
spPatterns <- cbind(spCoords,CoGAPS_Result@sampleFactors[barcodes,])

## Running scripts
optParams <- getSpatialParameters(spPatterns)
SpaceMarkers <- getInteractingGenes(data = fullMat, reconstruction = cgMat, spatialPatterns = spPatterns, refPattern = SpaceMarkersRefPattern, mode = SpaceMarkersMode)
SpaceMarkers$optParams <- optParams
for (i in seq(1,length(SpaceMarkers$interacting_genes)))
{
    filename <- paste0(gsub(x = cogapsFilePath,pattern = '.rds',replacement = paste0('_SpaceMarkers_',SpaceMarkersMode,'_')),names(SpaceMarkers$interacting_genes[[i]])[2],gsub(":",".",gsub(c(" "), "_", date())),".csv")
    write.csv(SpaceMarkers$interacting_genes[[i]],file = filename,row.names = F, quote = F)
}  
saveRDS(SpaceMarkers,gsub(x = cogapsFilePath,pattern = '.rds',replacement = paste0('_SpaceMarkers_',SpaceMarkersMode,'_',gsub(c(" "), "_", date()),'.rds')))

    