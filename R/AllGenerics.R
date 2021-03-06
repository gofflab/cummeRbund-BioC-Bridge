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
setGeneric("replicates",function(object) standardGeneric("replicates"))
setGeneric("conditions",function(object) standardGeneric("conditions"))
setGeneric("runInfo",function(object) standardGeneric("runInfo"))
setGeneric("varModel",function(object) standardGeneric("varModel"))
setGeneric("genes",function(object) standardGeneric("genes"))
setGeneric("isoforms",function(object) standardGeneric("isoforms"))
setGeneric("TSS",function(object) standardGeneric("TSS"))
setGeneric("CDS",function(object) standardGeneric("CDS"))
setGeneric("promoters",function(object) standardGeneric("promoters"))
setGeneric("splicing",function(object) standardGeneric("splicing"))
setGeneric("relCDS",function(object) standardGeneric("relCDS"))
setGeneric("getGene",function(object,geneId,sampleIdList=NULL) standardGeneric("getGene"))
setGeneric("getGenes",function(object,geneIdList,sampleIdList=NULL) standardGeneric("getGenes"))
setGeneric("getGeneId",function(object,idList) standardGeneric("getGeneId"))
setGeneric("findGene",function(object,query) standardGeneric("findGene"))
setGeneric("getFeatures",function(object,featureIdList,sampleIdList=NULL,...) standardGeneric("getFeatures"))
setGeneric("getSig", function(object,x,y,alpha=0.05,level="genes",method="BH",useCuffMTC=FALSE) standardGeneric("getSig"))
setGeneric("getSigTable", function(object,alpha=0.05,level='genes') standardGeneric("getSigTable"))
setGeneric("addFeatures",function(object, features, ...) standardGeneric("addFeatures"))
setGeneric("findSimilar",function(object,x,n,...) standardGeneric("findSimilar"))
setGeneric("makeGRanges",function(object,id,idField='transcript_id') standardGeneric("makeGRanges"))
setGeneric("sigMatrix",function(object,alpha=0.05,level='genes',orderByDist=F) standardGeneric("sigMatrix"))

###############
#CuffData
###############

#CuffData generics
setGeneric("features",function(object) standardGeneric("features"))
#setGeneric("featureNames",function(object) standardGeneric("featureNames"))
setGeneric("fpkm",function(object, features=FALSE,...) standardGeneric("fpkm"))
setGeneric("repFpkm",function(object,features=FALSE,...) standardGeneric("repFpkm"))
setGeneric("count",function(object,...) standardGeneric("count"))
setGeneric("filters",function(object) standardGeneric("filters"))
setGeneric("idField",function(object) standardGeneric("idField"))
setGeneric("tables",function(object) standardGeneric("tables"))
setGeneric("fpkmMatrix",function(object,...) standardGeneric("fpkmMatrix"))
setGeneric("repFpkmMatrix",function(object,...) standardGeneric("repFpkmMatrix"))
setGeneric("countMatrix",function(object,...) standardGeneric("countMatrix"))
setGeneric("repCountMatrix",function(object,...) standardGeneric("repCountMatrix"))
setGeneric("diffData",function(object, x, y, features=FALSE, lnFcCutoff=20, ...) standardGeneric("diffData"))
setGeneric("diffTable",function(object,...) standardGeneric("diffTable"))
setGeneric("getLevels",function(object) standardGeneric("getLevels"))
setGeneric("getRepLevels",function(object) standardGeneric("getRepLevels"))
setGeneric("MAplot",function(object,x,y,logMode=T,pseudocount=1,...) standardGeneric("MAplot"))
setGeneric("dispersionPlot",function(object) standardGeneric("dispersionPlot"))
setGeneric("makeRnk",function(object,x,y,filename,...) standardGeneric("makeRnk"))
setGeneric("csScatterMatrix",function(object,replicates=FALSE,logMode=TRUE,...) standardGeneric("csScatterMatrix"))
setGeneric("csVolcanoMatrix",function(object,replicates=FALSE,logMode=TRUE,...) standardGeneric("csVolcanoMatrix"))
setGeneric("MDSplot",function(object,replicates=FALSE,logMode=TRUE,...) standardGeneric("MDSplot"))
setGeneric("PCAplot",function(object,x="PC1", y="PC2",replicates=TRUE,pseudocount=1.0,scale=TRUE,...) standardGeneric("PCAplot"))
setGeneric("fpkmSCVPlot",function(object,FPKMLowerBound=1,...) standardGeneric("fpkmSCVPlot"))
setGeneric("csNMF",function(object,k,...) standardGeneric("csNMF"))

#CuffDist generics
setGeneric("distValues",function(object, x, y,...) standardGeneric("distValues"))
setGeneric("DB",function(object,...) standardGeneric("DB"))
#setGeneric("table",function(object,...) standardGeneric("table"))
setGeneric("type",function(object,...) standardGeneric("type"))
#setGeneric("testId",function(object,...) standardGeneric("testId"))

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
setGeneric("csFoldChangeHeatmap",function(object, control_condition, replicate_num=NULL, clustering='none', labCol=T, labRow=T, logMode=F, pseudocount=1.0, border=FALSE, heatscale=c(low='steelblue',mid='white',high='tomato'), heatMidpoint=0,fullnames=T,replicates=FALSE,method='none',heatRange=3, ...) standardGeneric("csFoldChangeHeatmap"))
setGeneric("csDistHeat",function(object, replicates=F, samples.not.genes=T, logMode=T, pseudocount=1.0, heatscale=c(low='lightyellow',mid='orange',high='darkred'), heatMidpoint=NULL, ...) standardGeneric("csDistHeat"))
setGeneric("csCluster",function(object, k, logMode=T, pseudocount=1,...) standardGeneric("csCluster"))
#setGeneric("csClusterPlot",function(clustering, pseudocount=1.0) standardGeneric("csClusterPlot"))
#setGeneric("diff",function(object) standardGeneric("diff"))
#setGeneric("annotation",function(object) standardGeneric("annotation"))
setGeneric("csSpecificity",function(object,logMode=T,pseudocount=1,relative=FALSE,...) standardGeneric("csSpecificity"))
setGeneric("csDendro",function(object,logMode=T,pseudocount=1,replicates=FALSE,...) standardGeneric("csDendro"))

##################
#CuffGene
##################
setGeneric("genePlot",function(object) standardGeneric("genePlot"))
setGeneric("makeGeneRegionTrack",function(object) standardGeneric("makeGeneRegionTrack"))
setGeneric("as.GRangesList",function(object,f="transcript") standardGeneric("as.GRangesList"))
setGeneric("csPie",function(object,...) standardGeneric("csPie"))


##############
#CuffFeature
##############
setGeneric("getGenome",function(object) standardGeneric	("getGenome"))

#CuffFeature plotting
setGeneric("expressionBarplot",function(object, logMode=FALSE, pseudocount=1.0, showErrorbars=TRUE, ...) standardGeneric("expressionBarplot"))
setGeneric("expressionPlot",function(object, logMode=FALSE, pseudocount=1.0, drawSummary=FALSE, sumFun=mean_cl_boot, showErrorbars=TRUE, ...) standardGeneric("expressionPlot"))
