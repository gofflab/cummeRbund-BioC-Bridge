#setGeneric("fpkm",		function(object) standardGeneric("fpkm"))
#setGeneric("fpkm<-",	function(object,value) standardGeneric("fpkm<-"))
#setGeneric("conf_lo",		function(object) standardGeneric("conf_lo"))
#setGeneric("conf_lo<-",	function(object,value) standardGeneric("conf_lo<-"))
#setGeneric("conf_hi",		function(object) standardGeneric("conf_hi"))
#setGeneric("conf_hi<-",	function(object,value) standardGeneric("conf_hi<-"))
#setGeneric("csApply", function(object, ...) standardGeneric("csApply"))
#setGeneric("csmelt", function(object) standardGeneric("csmelt"))
#setGeneric("write.fpkm", function(object,...) standardGeneric("write.fpkm"))
#
#setGeneric("expressionPlot",function(object, ...) standardGeneric("expressionPlot"))
##setGeneric("melt", function(data) standardGeneric("melt"))
#setGeneric("csHeatmap",function(object,...) standardGeneric("csHeatmap"))
#setGeneric("csDensity",function(object, ...) standardGeneric("csDensity"))
#setGeneric("csHist",function(object ,...) standardGeneric("csHist"))
#setGeneric("csBoxplot", function(object, ...) standardGeneric("csBoxplot"))
#setGeneric("csScatter",function(object,x,y,...) standardGeneric("csScatter"))
#
#setGeneric("diffData",       function(object) standardGeneric("diffData"))
#setGeneric("diffData<-",     function(object, value) standardGeneric("diffData<-"))

#New for cummeRbund 0.1.2
setGeneric("loadGenes",function(fpkmFile,..) standardGeneric("loadGenes"))
setGeneric("loadIsoforms",function(fpkmFile,..) standardGeneric("loadIsoforms"))
setGeneric("loadTSS",function(fpkmFile,..) standardGeneric("loadTSS"))
setGeneric("loadCDS",function(fpkmFile,..) standardGeneric("loadCDS"))

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