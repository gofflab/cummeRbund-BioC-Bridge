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
setGeneric("getGene",function(object,...) standardGeneric("getGene"))
setGeneric("getGenes",function(object,...) standardGeneric("getGenes"))
setGeneric("addFeatures",function(object,...) standardGeneric("addFeatures"))

###############
#CuffData
###############

#CuffData generics
setGeneric("features",function(object,...) standardGeneric("features"))
setGeneric("featureNames",function(object,...) standardGeneric("featureNames"))
setGeneric("fpkm",function(object,...) standardGeneric("fpkm"))
setGeneric("fpkmMatrix",function(object) standardGeneric("fpkmMatrix"))
setGeneric("diffData",function(object,x,y,...) standardGeneric("diffData"))
setGeneric("getLevels",function(object) standardGeneric("getLevels"))
setGeneric("getFeatures",function(object,...) standardGeneric("getFeatures"))

#CuffDist generics
setGeneric("values",function(object,...) standardGeneric("values"))

#CuffData plotting
setGeneric("csDensity",function(object, ...) standardGeneric("csDensity"))
setGeneric("csScatter",function(object,x,y,...) standardGeneric("csScatter"))
setGeneric("csVolcano",function(object,x,y,...) standardGeneric("csVolcano"))
setGeneric("csBoxplot",function(object,...) standardGeneric("csBoxplot"))

###################
#CuffGeneSet
####################



#################
#CuffFeatureSet
#################
setGeneric("csHeatmap",function(object,...) standardGeneric("csHeatmap"))
setGeneric("csCluster",function(object,...) standardGeneric("csCluster"))



##################
#CuffGene
###################

##############
#CuffFeature
##############

#CuffFeature plotting
setGeneric("expressionBarplot",function(object,...) standardGeneric("expressionBarplot"))
setGeneric("expressionPlot",function(object, ...) standardGeneric("expressionPlot"))