patternSpotter <- function(obj, locs, geneList = NULL, spotIndices = NULL, patternList = NULL, threshold = 1, plotTitle = "Generic title", radius = 2, lwd = 0.1, groups = NA, colors_arr = NULL)
{
    require(STdeconvolve)
    if (is.null(patternList)) patternList <- colnames(obj@featureLoadings)
    if (is.null(geneList)) geneList <- rownames(obj@featureLoadings)
    if (is.null(spotIndices)) spotIndices <- rownames(obj@sampleFactors)

    geneFE <- obj@featureLoadings[geneList,] %*% diag(colSums(obj@sampleFactors[spotIndices,]))
    geneFE <- geneFE/apply(geneFE,1,sum)
    colnames(geneFE) <- colnames(obj@featureLoadings)
    geneFE <- geneFE[,patternList]
    if (is.null(colors_arr)) colors_arr <- dittoSeq::dittoColors()[1:length(patternList)]
    locs <- locs[spotIndices,]
    spotFE <- obj@sampleFactors[spotIndices,] %*% diag(colSums(obj@featureLoadings[unique(geneList),]))
    colnames(spotFE) <- colnames(obj@featureLoadings)
    if (length(patternList)<length(obj@featureLoadings))
    {
        spotFE <- cbind(spotFE[,patternList], Others = apply(spotFE[,setdiff(colnames(spotFE),patternList)],1,sum))
        colors <- c(colors,"#FFFFFF")
    }   
    present <- apply(spotFE,1,sum)>threshold
    vizAllTopics(spotFE[present,], locs[present,],topicCols = colors_arr, r=radius, lwd = lwd, plotTitle = plotTitle, groups = groups)
    #return(list(spotFE = spotFE,geneFE = geneFE))
}