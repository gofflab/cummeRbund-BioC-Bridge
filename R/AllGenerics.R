#New for cummeRbund 0.1.2
#setGeneric("loadGenes",function(fpkmFile,..) standardGeneric("loadGenes"))
#setGeneric("loadIsoforms",function(fpkmFile,..) standardGeneric("loadIsoforms"))
#setGeneric("loadTSS",function(fpkmFile,..) standardGeneric("loadTSS"))
#setGeneric("loadCDS",function(fpkmFile,..) standardGeneric("loadCDS"))

##################
#CuffSet
#################3

#CuffSet generics
setGeneric("samples",function(object) standardGeneric("samples"))
setGeneric("getGene",function(object,geneId,sampleIdList=NULL) standardGeneric("getGene"))
setGeneric("getGenes",function(object,geneIdList,sampleIdList=NULL) standardGeneric("getGenes"))
setGeneric("addFeatures",function(object, features, ...) standardGeneric("addFeatures"))
setGeneric("findSimilar",function(object,x,n) standardGeneric("findSimilar"))

###############
#CuffData
###############

#CuffData generics
setGeneric("features",function(object) standardGeneric("features"))
setGeneric("featureNames",function(object) standardGeneric("featureNames"))
setGeneric("fpkm",function(object, features=FALSE) standardGeneric("fpkm"))
setGeneric("fpkmMatrix",function(object) standardGeneric("fpkmMatrix"))
setGeneric("diffData",function(object, x, y, features=FALSE, lnFcCutoff=20, ...) standardGeneric("diffData"))
setGeneric("getLevels",function(object) standardGeneric("getLevels"))


#CuffDist generics
setGeneric("distValues",function(object,...) standardGeneric("distValues"))

#CuffData plotting
setGeneric("csDensity",function(object, logMode=TRUE, pseudocount=1.0, labels, features=FALSE, ...) standardGeneric("csDensity"))
setGeneric("csScatter",function(object, x, y, logMode=TRUE, pseudocount=1.0, labels, smooth=FALSE, ...) standardGeneric("csScatter"))
setGeneric("csVolcano",function(object, x, y, features=F, ...) standardGeneric("csVolcano"))
setGeneric("csBoxplot",function(object, logMode=T, ...) standardGeneric("csBoxplot"))

###################
#CuffGeneSet
###################



#################
#CuffFeatureSet
#################
setGeneric("csHeatmap",function(object,rescaling='none', clustering='none', labCol=T, labRow=T, logMode=T, pseudocount=1.0, border=FALSE, heatscale= c(low='darkred',mid='orange',high='white'), heatMidpoint=NULL, ...) standardGeneric("csHeatmap"))
setGeneric("csCluster",function(object, k, iter.max=100, ...) standardGeneric("csCluster"))

##################
#CuffGene
##################

##############
#CuffFeature
##############

#CuffFeature plotting
setGeneric("expressionBarplot",function(object, logMode=FALSE, pseudocount=1.0, showErrorbars=TRUE, ...) standardGeneric("expressionBarplot"))
setGeneric("expressionPlot",function(object, logMode=FALSE, pseudocount=1.0, drawSummary=FALSE, sumFun=mean_cl_boot, showErrorbars=TRUE, ...) standardGeneric("expressionPlot"))