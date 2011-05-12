# methods-CuffSet.R
# 
#Introduces the cuffSet Class for analysis,manipulation, and visualization of Cufflinks output data
#
# Author: Loyal A. Goff
#
###############################################################################

#Initialize
setMethod("initialize","cuffSet",
		function(.Object,
				assayData,
				phenoData,
				featureData,
				#diffData,
				fpkm = new("matrix"),
				conf_lo = new("matrix"),
				conf_hi = new("matrix"),
				...){
			if (missing(assayData)){
				if(missing(phenoData))
					phenoData<- annotatedDataFrameFrom(fpkm, byrow=FALSE)
				if(missing(featureData))
					featureData<- annotatedDataFrameFrom(fpkm,byrow=TRUE)
				.Object <-callNextMethod(.Object,
						phenoData = phenoData,
						featureData = featureData,
						fpkm = fpkm,
						conf_lo = conf_lo,
						conf_hi = conf_hi,
						...)
			} else if (missing(fpkm)) {
				if (missing(phenoData))
					phenoData<-annotatedDataFrameFrom(assayData,byrow=FALSE)
				if (missing(featureData)) 
					featureData<-annotatedDataFrameFrom(assayData,byrow=TRUE)
				.Object<-callNextMethod(.Object,
						assayData = assayData,
						phenoData = phenoData,
						featureData = featureData,
						...)
			} else stop("provide at most one of 'assayData' or 'fpkm' to initialize cuffSet", call.=FALSE)
			#Add diffData if available
							
			.harmonizeDimnames(.Object)
		}
)

#This was adopted from ExpressionSet definition.
.harmonizeDimnames <- function(object) {
	err <- function(conflicts)
		stop("assayData element dimnames conflict: ",
				paste(names(conflicts), collapse=", "))
	okNames <- list(featureNames(featureData(object)),
			sampleNames(phenoData(object)))
	#dimNames <- .assayDataDimnames(assayData(object))
	dimNames <- Biobase:::.assayDataDimnames(assayData(object))
	dimConflict <- function(dimNames, okNames, dim) {
		nm <- lapply(dimNames, "[[", dim)
		isConflict <- !sapply(nm, function(x, y) {
					is.null(x) || all.equal(x, y, check.attr=FALSE)
				}, okNames[[dim]])
		isNamed <- sapply(lapply(nm, names), length) > 0
		isNull <- sapply(nm, is.null)
		if (all(!isConflict & !isNamed & !isNull))
			return (FALSE)
		if (any(isConflict & !isNull))
			err(isConflict[!isNull])
		TRUE
	}
	if (dimConflict(dimNames, okNames, 1))
		featureNames(assayData(object)) <- okNames[[1]]
	if (dimConflict(dimNames, okNames, 2))
		sampleNames(assayData(object)) <- okNames[[2]]
	object
}

setValidity("cuffSet",function(object) {
			msg<-Biobase:::validMsg(NULL,Biobase:::isValidVersion(object,"cuffSet"))
			#msg<-validMsg(msg,assayDataValidMembers(assayData(object),c("fpkm")))
			msg<-validMsg(NULL,assayDataValidMembers(assayData(object),c("fpkm","conf_lo","conf_hi")))
			#msg<-Biobase:::validMsg(NULL,Biobase:::assayDataValidMembers(assayData(object),c("fpkm")))
			
			if (is.null(msg)) TRUE else msg
		}
)


setAs("cuffSet", "data.frame",
		function (from) data.frame(t(fpkm(from)), pData(from)))


as.data.frame.cuffSet <- function(x, row.names=NULL, optional=FALSE, ...)
	as(x, "data.frame")


setMethod("fpkm", signature(object="cuffSet"),
		function(object) assayDataElement(object,"fpkm")
)

setReplaceMethod("fpkm", signature(object="cuffSet",value="matrix"),
		function(object,value) assayDataElementReplace(object,"fpkm",value)
)

setMethod("conf_lo", signature(object="cuffSet"),
		function(object) assayDataElement(object,"conf_lo")
)

setReplaceMethod("conf_lo", signature(object="cuffSet",value="matrix"),
		function(object,value) assayDataElementReplace(object,"conf_lo",value)
)

setMethod("conf_hi", signature(object="cuffSet"),
		function(object) assayDataElement(object,"conf_hi")
)

setReplaceMethod("conf_hi", signature(object="cuffSet",value="matrix"),
		function(object,value) assayDataElementReplace(object,"conf_hi",value)
)

setMethod("geneNames",signature(object="cuffSet"), 
		function(object){
			.Defunct("featureNames","Biobase")
		}
)

setReplaceMethod("geneNames", signature(object="cuffSet",value="character"),
		function(object,value){
			.Defunct("featureNames","Biobase")
		}
)


.csApply <- function(object,MARGIN,FUN,...) {
	parent <-environment(FUN)
	if (is.null(parent))
		parent <- emptyenv()
	e1 <- new.env(parent=parent)
	multiassign(names(pData(object)),pData(object),env=e1)
	environment(FUN) <- e1
	apply(fpkm(object),MARGIN,FUN,...)
}


setMethod("csApply",signature = signature(object="cuffSet"), .csApply
)


setMethod("write.fpkm",
		signature(object="cuffSet"),
		function(object,file="tmp.fpkm",quote=FALSE,sep="\t", col.names=NA,...){
			write.table(fpkm(object), file=file,quote=quote,sep=sep,col.names=col.names,...)
		}
)

#########
#Data import and initialization
#########

readCufflinks <- function(fpkmFile,
		phenoDataFile,
		featureDataFile,
		experimentDataFile,
		notesFile,
		path,
		annotation,
		##arguments to read.* methods
		fpkmArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, na.string=na.string, ...),
		phenoDataArgs = list(sep=sep, header=header, row.names = row.names, quote = quote, na.string=na.string, stringsAsFactors = stringsAsFactors, ...),
		featureDataArgs = list(sep=sep, header=header, row.names = row.names, quote = quote, na.string=na.string, stringsAsFactors = stringsAsFactors, ...),
		experimentDataArgs = list(sep=sep, header=header, row.names=row.names, quote=quote, na.string=na.string, stringsAsFactors = stringsAsFactors, ...),
		sep = "\t",
		na.string = "-",
		header = TRUE,
		quote = "",
		stringsAsFactors = FALSE, 
		row.names=1L,
		...) {
	
	if (missing(fpkmFile))
		stop("fpkmFile can not be missing!")
	
	#Read Primary file (fpkmfFile) 
	fpkmArgs$file=fpkmFile
	full = as.data.frame(do.call(read.table,fpkmArgs))
	#print(head(full))
	fp <- as.matrix(full[,-c(1:9)])
	#print(head(fp))
	nCols <- dim(fp)[2]
	conf_lo <- fp[,seq(2,nCols,3)]
	conf_hi <- fp[,seq(3,nCols,3)]
	fp <- fp[,seq(1,nCols,3)]
	#
	#print(dim(fp))
	#print(dim(conf_lo))
	#print(dim(conf_hi))
	
	#fix colnames
	newColnames<-unlist(lapply(colnames(fp),function(x) substr(x,1,nchar(x)-5)))
	#print(newColnames)
	colnames(fp)<-newColnames
	colnames(conf_lo)<-newColnames
	colnames(conf_hi)<-newColnames
	
	featureInfo <- full[,c(1:9)]
	
	#Check phenoData
	if(!missing(phenoDataFile)) {
		phenoDataArgs$file = phenoDataFile
		pd = do.call(read.AnnotatedDataFrame, phenoDataArgs)
		if (!identical(sampleNames(pd), colnames(fp)))
			stop("Column names of fpkm matrix must be idential to\n",
					"the sample names of the phenodata table.\n",
					"You could use 'options(error=recover)' to compare the\n",
					"values of 'sampleNames(pd)' and 'colnames(fp)'.\n")
		
	} else {
		pd = annotatedDataFrameFrom(fp,byrow=FALSE)
	}
	
	#Check feature data
	if(!missing(featureDataFile)){
		featureDataArgs$file = featureDataFile
		fd = do.call(read.AnnotatedDataFrame,featureDataArgs)
	} else {
		fd = new("AnnotatedDataFrame", data=featureInfo)
	}
	
	#Create New object
	obj = new("cuffSet",fpkm=fp,conf_lo=conf_lo,conf_hi=conf_hi,phenoData=pd,featureData=fd)
	
	#set default columns for featureData and phenoData
	fData(obj)$symbol <- rownames(fData(obj))
	pData(obj)$sample <- rownames(pData(obj))
	
	#Additional experimental data
	##experimentData
	if(!missing(experimentDataFile))
		experimentDataArgs$file = experimentDataFile
	if(!is.null(experimentDataArgs$file))
		experimentData(obj)<-do.call(read.MIAME, experimentDataArgs)

	#annotation
	if(!missing(annotation))
		annotation(obj) <- annotation
	
	#notes
	if(!missing(notesFile))
		notes(obj) <- readLines(notesFile)
	
	#validObject(obj)
	
	obj
}

#Interpreted from Bioconductor
#readCufflinks <- function(fpkmFile,
#		phenoDataFile,
#		experimentDataFile,
#		notesFile,
#		path,
#		annotation,
#		##arguments to read.* methods
#		fpkmArgs = list(sep=sep, header=header, row.names = row.names, quote=quote, ...),
#		phenoDataArgs = list(sep=sep, header=header, row.names = row.names, quote = quote, stringsAsFactors = stringsAsFactors, ...),
#		experimentDataArgs = list(sep=sep, header=header, row.names=row.names, quote=quote, stringsAsFactors = stringsAsFactors, ...),
#		sep = "\t",
#		header = TRUE,
#		quote = "",
#		stringsAsFactors = FALSE, 
#		row.names=1L,
#		widget = getOption("BioC")$Base$use.widgets, 
#		...) {
#	
#	if (!missing(widget) && widget != FALSE)
#		stop("sorry, widgets are not yet available")
#	
#	##fpkm
#	if (missing(fpkmFile))
#		stop("fpkmFile can not be missing!")
#	fpkmArgs$file=fpkmFile
#	fp = as.matrix(do.call(read.table,fpkmArgs))
#	
#	##phenoData
#	if (!missing(phenoDataFile)) {
#		phenoDataArgs$file = phenoDataFile
#		pd = do.call(read.AnnotatedDataFrame, phenoDataArgs)
#		if (!identical(sampleNames(pd), colnames(fp)))
#			stop("Column names of fpkm matrix must be idential to\n",
#					"the sample names of the phenodata table.\n",
#					"You could use 'options(error=recover)' to compare the\n",
#					"values of 'sampleNames(pd)' and 'colnames(fp)'.\n"
#			)
#	} else {
#		pd = annotatedDataFramFrom(fp,byrow=FALSE)
#	}
#	
#	obj = new("cuffSet",fpkm=fp,phenoData=pd)
#	
#	##experimentData
#	if(!missing(experimentDataFile))
#		experimentDataArgs$file = experimentDataFile
#	if(!is.null(experimentDataArgs$file))
#		experimentData(obj)<-do.call(read.MIAME, experimentDataArgs)
#	
#	#annotation
#	if(!missing(annotation))
#		annotation(obj) <- annotation
#	
#	#notes
#	if(!missing(notesFile))
#		notes(obj) <- readLines(notesFile)
#	
#	validObject(obj)
#	obj
#}

############
#Data Manipulation
############

#.melt<-function(X){
#	datamelt<-melt(fpkm(X),varnames=c("symbol","sample"))
#	datamelt<-merge(datamelt,pData(X),by.x="sample",by.y="samp")
#	datamelt$trt<-factor(datamelt$trt)
#	datamelt
#}

.melt<-function(object){
	fpkmmelt<-melt(fpkm(object),varnames=c("symbol","sample"))
	colnames(fpkmmelt) = c("symbol","sample","fpkm")
	conf_lomelt<-melt(conf_lo(object),varnames=c("symbol","sample"))
	colnames(conf_lomelt)<-c("symbol","sample","conf_lo")
	conf_himelt<-melt(conf_hi(object),varnames=c("symbol","sample"))
	colnames(conf_himelt)<-c("symbol","sample","conf_hi")
	datamelt<-merge(pData(object),fpkmmelt,by.x="sample",by.y="sample")
	datamelt<-merge(datamelt,fData(object),by.x="symbol",by.y="symbol")
	datamelt<-merge(datamelt,conf_lomelt)
	datamelt<-merge(datamelt,conf_himelt)
	datamelt
}

setMethod("csmelt",signature=signature(object="cuffSet"),.melt)

#setMethod("melt",signature=signature(data="cuffSet"),.melt)


#################
#Plots
#################

.density<-function(object, logMode=TRUE, pseudocount=0.0001, labels, ...) {
	#Test class
	if (is(object, 'cuffSet')) {
		datamelt<-csmelt(object)
	} else {
		stop('Un-supported class of x.')
	}
	if(logMode) datamelt$fpkm<-datamelt$fpkm+pseudocount
	p <- ggplot(datamelt)
			if(logMode) {
				p<-p+geom_density(aes(x=log2(fpkm),group=sample,color=sample,fill=sample),alpha=I(1/3))
			} else {
				p<-p+geom_density(aes(x=fpkm,group=sample,color=sample,fill=sample),alpha=I(1/3))
			}
	if (!missing(labels)){
		labeldata<-subset(datamelt,gene_short_name %in% labels)
		if(logMode){
			labeldata$fpkm<-log2(labeldata$fpkm)
		}
		#print(labeldata)
		#print(str(p))
		p<-p+geom_point(data=labeldata,aes(x=fpkm,y=0,color=sample),size=2)
		p<-p+geom_text(data=labeldata,aes(x=fpkm,y=0,label=gene_short_name),hjust=0,vjust=0,size=2,angle=45)
	}
	p
	
}

setMethod('csDensity', signature(object='cuffSet'), .density)

setMethod('csHist', signature(object='cuffSet'), 
		function(object, ...) {
			csDensity(object, ...)
		}
)

.boxplot<-function(object,logMode=TRUE,main,...){
	tmp<-description(object)
	if (missing(main) && (is(tmp,"MIAME")))
		main <- tmp@title
	datamelt<-csmelt(object)
	p <- ggplot(datamelt)
	if(logMode) {
		p<-p+geom_boxplot(aes(x=sample,y=log2(fpkm),fill=sample),size=0.3,alpha=I(1/3))
	} else {
		p<-p+geom_boxplot(aes(x=sample,y=fpkm,fill=sample),alpha=I(1/3),size=0.3)
	}
	p
	
}

setMethod("csBoxplot",signature(object="cuffSet"),.boxplot)

.expressionPlot<-function(object,logMode=TRUE, pseudocount=0.0001, drawSummary=FALSE, sumFun=mean_cl_boot,...){
	datamelt<-csmelt(object)
	if(logMode){
		datamelt$fpkm<-log2(datamelt$fpkm+pseudocount)
	}
	p	<-	ggplot(datamelt) + 
			geom_line(aes(x=sample,y=fpkm,group=symbol),alpha=0.2) 
	
	#drawMean
	if(drawSummary){
		p <- p + stat_summary(aes(x=sample,y=fpkm,group=1),fun.data=sumFun,color="red",fill="red",alpha=0.2,size=1.1,geom="smooth")
	}
	
	p
}

setMethod("expressionPlot",signature(object="cuffSet"),.expressionPlot)

.barplot<-function(object,logMode=TRUE, pseudocount=0.0001, ...){
	datamelt<-csmelt(object)
	p	<-	ggplot(datamelt,aes(x=sample,y=fpkm,fill=sample)) +
			geom_bar() +
			geom_errorbar(aes(ymin=conf_lo,ymax=conf_hi),size=0.15) +
			facet_grid(symbol~.,scales="free") +
			opts(axis.text.x=theme_text(hjust=0,angle=-90))
	#This is a cool example plot for isoform expression faceted by gene_short_name, colored by sample
	#p+geom_bar(aes(fill=sample),position="dodge",stat="identity")+facet_wrap('gene_short_name',scales="free")+opts(legend.position="none",axis.text.x=theme_text(size=8,angle=-90,hjust=0))
	p
}

setMethod("expressionBarplot",signature(object="cuffSet"),.barplot)


.heatmap<-function(object,logMode=TRUE, ...){
	datamelt<-csmelt(object)
	p <-	ggplot(datamelt,aes(x=sample,y=symbol))
	if(logMode){ 
		p<-p+geom_tile(aes(fill=log2(fpkm+0.0001)))
	}else{
		p<-p+geom_tile(aes(fill=fpkm))
	}
	p<- p + scale_fill_gradient(low="white",high="red")
	p
}

#######################
#The following is borrowed from Malarkey and is not yet ready for prime time...
#I would like to replace the clustering here with JSdistance on rows and/or columns and make
#this package work with cuffSet objects by default. 
#There is no genericMethod yet, goal is to replace .heatmap with .ggheat for genericMethod 'csHeatmap'

.ggheat<-function(object, rescaling='none', clustering='none', labCol=T, labRow=T, logMode=T, pseudocount=0.001, 
		border=FALSE, heatscale= c(low='white',high='red'),...) {
	## the function can be be viewed as a two step process
	## 1. using the rehape package and other funcs the data is clustered, scaled, and reshaped
	## using simple options or by a user supplied function
	## 2. with the now resahped data the plot, the chosen labels and plot style are built
	
	#require(reshape)
	#require(ggplot2)
	
	m=fpkm(object)
	#remove genes with no expression in any condition
	m=m[!apply(m,1,sum)==0,]

	## you can either scale by row or column not both! 
	## if you wish to scale by both or use a different scale method then simply supply a scale
	## function instead NB scale is a base funct
	
	if(is.function(rescaling))
	{ 
		m=rescaling(m)
	} 
	else 
	{
		if(rescaling=='column') 
			m=scale(m, center=T)
		if(rescaling=='row') 
			m=t(scale(t(m),center=T))
	}
	
	## I have supplied the default cluster and euclidean distance (JSdist) - and chose to cluster after scaling
	## if you want a different distance/cluster method-- or to cluster and then scale
	## then you can supply a custom function 
	
	if(is.function(clustering)) 
	{
		m=clustering(m)
	}else
	{
		if(clustering=='row')
			m=m[hclust(JSdist(makeprobs(t(m))))$order, ]
		if(clustering=='column')  
			m=m[,hclust(JSdist(makeprobs(m)))$order]
		if(clustering=='both')
			m=m[hclust(JSdist(makeprobs(t(m))))$order ,hclust(JSdist(makeprobs(m)))$order]
	}
	## this is just reshaping into a ggplot format matrix and making a ggplot layer
	
	rows=dim(m)[1]
	cols=dim(m)[2]
	if(logMode) {
		melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows) ,melt(log2(m+pseudocount)))
	}else{
		melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows) ,melt(m))
	}
	g=ggplot(data=melt.m)
	
	## add the heat tiles with or without a white border for clarity
	
	if(border==TRUE)
		g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value),colour='white')
	if(border==FALSE)
		g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value))
	
	## add axis labels either supplied or from the colnames rownames of the matrix
	
	if(labCol==T) 
		g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=colnames(m))
	if(labCol==F) 
		g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=rep('',cols))
	
	if(labRow==T) 
		g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rownames(m))	
	if(labRow==F) 
		g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rep('',rows))	
	
	## get rid of grey panel background and gridlines
	
	g2=g2+opts(panel.grid.minor=theme_line(colour=NA), panel.grid.major=theme_line(colour=NA),
			panel.background=theme_rect(fill=NA, colour=NA))
	
	## finally add the fill colour ramp of your choice (default is blue to red)-- and return
	return(g2+scale_fill_continuous("", heatscale[1], heatscale[2]))
	
}

#####################

#setMethod("csHeatmap",signature(object="cuffSet"),.heatmap)
setMethod("csHeatmap",signature(object="cuffSet"),.ggheat)


.scatter<-function(object,x=1,y=2,logMode=TRUE,pseudocount=0.00001, labels, smooth=FALSE, ...){
	fp<-fpkm(object)
	fp<-as.data.frame(fp)
	if(logMode){
		fp<-fp+pseudocount
	}
	#Find column indices

	if(is.numeric(x)){
		xname<-colnames(fp)[x]
	}else{
		xname<-x
	}
	if(is.numeric(y)){
		yname<-colnames(fp)[y]
	}else{
		yname<-y
	}
	print(c(xname,yname))
	p <- ggplot(fp) 
	p <- p + aes_string(x=xname,y=yname)
	p <- p + geom_point(size=1.2,alpha=I(1/3)) + geom_abline(intercept=0,slope=1,linetype=2) + geom_rug(size=0.5,alpha=0.01)
	
	#add smoother
	if(smooth){
		p <- p + stat_smooth(method="lm",fill="blue",alpha=0.2)
	}
	
	#Add highlights from labels
	if(!missing(labels)){
		labelIdx<-fData(object)$gene_short_name %in% labels
		labelfp<-fp[labelIdx,]
		labelfp$gene_short_name<-fData(object)$gene_short_name[labelIdx]
		#print(head(labelfp))
		p <- p + geom_point(data=labelfp,size=1.2,color="red")
		p <- p + geom_text(data=labelfp,aes(label=gene_short_name),color="red",hjust=0,vjust=0,angle=45,size=2)
	}

	
	if(logMode){
		p <- p + scale_y_log2() + scale_x_log2()
	}
	p
}

setMethod("csScatter",signature(object="cuffSet"),.scatter)

##################
#Clustering
##################
#Find a quick way to generate JSdistance matrix.
#This can be used as input for pam (from 'cluster' library) for quick partioning around mediods and clustering
#Example:
#cuff.JS<-JSdist(makeprobs(t(fpkm(cuff)))) #Current row-based JS distance calculation in R (very slow for large matrices).
#cuff.pam<-pam(cuff.JS,12) #PAM can take a dissimilarity matrix as input and a value of 'k'.
#
#fData(cuff.known)$cluster<-as.factor(cuff.pam$clustering) #Add cluster information to cuffSet featureData.
#p<-expressionPlot(cuff.known)
#p+facet_wrap("cluster")+aes_string(color='cluster') #Expression plots by cluster.
#


#############
#Unit Tests
##############
#library(Biobase)
#library(reshape)
#library(ggplot2)
#library(cummeRbund)

#path<-"/Users/lgoff/Documents/workspace/cummeRbund/data/"
#fpkm_tracking_file <-"genes.fpkm_tracking"
#pData_file <- "pData.txt"
#fname<-paste(path,fpkm_tracking_file,sep="")
#pDataFile = paste(path,pData_file,sep="")
#a<-readCufflinks(fpkmFile=fname,phenoDataFile=pDataFile)

#Make fake lincRNA featureData
#fData(a)$lincRNA<-sample(c("lincRNA","pcGene"),size=999,replace=T)


#melt cuffSet
#datamelt<-csmelt(a)

#plots
#p<-ggplot(datamelt) + scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1")
#p + geom_density(aes(x=fpkm,group=sample,fill=sample),alpha=I(1/3)) + scale_x_log() +
#		opts(legend.position="right",
#				plot.title = theme_text(size = 17), 
#				axis.title.x = theme_text(size = 15), 
#				axis.title.y = theme_text(size = 15, angle = 90), 
#				strip.text.x = theme_text(size = 15), 
#				strip.text.y = theme_text(size = 15, angle = -90)) + 
		#facet_grid(reliability ~ .)
#		facet_grid(. ~ .)


