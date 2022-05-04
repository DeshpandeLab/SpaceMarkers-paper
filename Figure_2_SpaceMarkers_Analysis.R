rm(list = ls())
library(Matrix)
library(SpatialExperiment)
library(SpaceMarkers)
setwd("~/FertigLab/SpaceMarkers-paper")
source('./patternSpotter.R', echo=TRUE)
dataFolder <- "PDAC_J"
rngtools::RNGseed(123)
## Set these parameters
# SpInMarkersMode: defaut mode is "residual". You can also set "DE" mode for Differential Expression mode.
SpaceMarkersMode = "DE"  
# SpinMarkersRefPattern is the pattern whose "interaction" with every other pattern we want to study. If refPattern is not explicitly assigned, the code assumes Pattern_1 to be refPattern. Here the met-PDAC pattern "Pattern_6" is treated as the reference.
SpaceMarkersRefPattern = "Pattern_6" 

countPath <- list.files(dataFolder, pattern = "countsMatrix",full.names = T)
imagePath <- list.files(dataFolder, pattern = "image",full.names = T)
barcodePath <- list.files(dataFolder, pattern = "barcode",full.names = T)

cogapsFilePath <- list.files("PDAC_J", pattern = "_8Pattern.rds",full.names = T)

resultsFolderNames <- paste0(gsub(".rds","",cogapsFilePath),"_Results")
sapply(resultsFolderNames, function(rName) if (!dir.exists(rName)){dir.create(rName)})
  
## Read the Full Count Matrix
fullMat <- t(readRDS(countPath))

CoGAPS_Result <- readRDS(cogapsFilePath);
cgMat <- CoGAPS_Result@featureLoadings %*% t(CoGAPS_Result@sampleFactors)

features <- intersect(rownames(fullMat),rownames(CoGAPS_Result@featureLoadings))
barcodes <- intersect(colnames(fullMat),rownames(CoGAPS_Result@sampleFactors))

fullMat <- fullMat[features,barcodes]
cgMat <- cgMat[features,barcodes]

barcodeInfo <- readRDS(barcodePath)
spCoords <- barcodeInfo %>% select(c("barcode","imagecol","imagerow"))
colnames(spCoords)[colnames(spCoords)=="imagerow"]="y"
colnames(spCoords)[colnames(spCoords)=="imagecol"]="x"
spCoords <- spCoords[barcodes,]
spPatterns <- cbind(spCoords,CoGAPS_Result@sampleFactors[barcodes,])


## patternMarkers

PMall <- patternMarkers(CoGAPS_Result)
PMcut <- patternMarkers(CoGAPS_Result, threshold = "cut")
names(PMcut$PatternMarkers) <- names(PMall$PatternMarkers)
maxCounts <- apply(fullMat,1,max)
PMfiltAll <- lapply(PMall$PatternMarkers, function(PM) PM[maxCounts[PM]>2])
PMfiltCut <- lapply(PMcut$PatternMarkers, function(PM) PM[maxCounts[PM]>2])

dfAll <- data.frame(lapply(PMfiltAll, "length<-",max(lengths(PMfiltAll))))
dfCut <- data.frame(lapply(PMfiltCut, "length<-",max(lengths(PMfiltCut))))
write.csv(dfAll, file = gsub(pattern = ".rds", replacement = "_patternMarkers_All.csv",cogapsFilePath),row.names = FALSE, na = "")
write.csv(dfCut, file = gsub(pattern = ".rds", replacement = "_patternMarkers_Cut.csv",cogapsFilePath),row.names = FALSE, na = "")

fullMat <- log2(1+fullMat)
## Running scripts
optParams <- getSpatialParameters(spPatterns)
SpaceMarkers <- getInteractingGenes(data = fullMat, reconstruction = cgMat, optParams = optParams, spatialPatterns = spPatterns, refPattern = SpaceMarkersRefPattern, mode = SpaceMarkersMode)
SpaceMarkers$optParams <- optParams

for (i in seq(1,length(SpaceMarkers$interacting_genes)))
{
  filename <- paste0(gsub(x = cogapsFilePath,pattern = '.rds',replacement = paste0('_SpaceMarkers_',SpaceMarkersMode,'_')),names(SpaceMarkers$interacting_genes[[i]])[2],gsub(":",".",gsub(c(" "), "_", date())),".csv")
  write.csv(SpaceMarkers$interacting_genes[[i]],file = filename, row.names = F, quote = F)
}  
saveRDS(SpaceMarkers,gsub(x = cogapsFilePath,pattern = '.rds',replacement = paste0('_SpaceMarkers_',SpaceMarkersMode,'_',gsub(c(" "), "_", date()),'.rds')))

## patternSpotter Visualization
pos <- spCoords[,c("x","y")]; pos$y = -pos$y
patternSpotter(obj = CoGAPS_Result, locs = pos, patternList = c("Pattern_6","Pattern_9"), spotIndices = NULL, threshold = 2, plotTitle = "metastatic PDAC vs Immune cells",radius = 1.5)
  