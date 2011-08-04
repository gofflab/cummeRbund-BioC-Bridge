########################
#methods-CuffFeatureSet.R
#
#Author: Loyal A. Goff
#
#Date created: 5-17-2011
#
#Description:
#########################

#################
#Initialize		#
#################
setMethod("initialize","CuffFeatureSet",
		function(.Object,
				annotation=data.frame(),
				fpkm=data.frame(),
				diff=data.frame(),
				... ){
			.Object<-callNextMethod(.Object,
					annotation=annotation,
					fpkm=fpkm,
					diff=diff,
					...)				
		}
)

#################
#Validate		#
#################
#TODO: Add validity constraints
setValidity("CuffFeatureSet",function(object){
			TRUE
		}
)		

#################
#Class Methods	#
#################
#setMethod("show","CufFFeatureSet",function(object){
#		cat(class(object),"instance")
#		}
#)

#################
#Subsetting		#
#################
#TODO: Add subset methods to return a CuffFeature object
#setMethod("[","CuffFeatureSet",function(object,featureID){
#			
#		}
#)

#################
#Accessors
##################
.samples<-function(object){
	res<-fpkm(object)$sample_name
	res
}

setMethod("samples","CuffFeatureSet",.samples)

.fpkm<-function(object,features=FALSE){
	if (features){
		return (merge(object@annotation,object@fpkm))
	}else{
		return(object@fpkm)
	}
}
setMethod("fpkm",signature(object="CuffFeatureSet"),.fpkm)

.featureNames<-function(object){
	data.frame(tracking_id=object@annotation[,1],gene_short_name=object@annotation$gene_short_name)
}

setMethod("featureNames",signature(object="CuffFeatureSet"),.featureNames)

.features<-function(object){
	object@annotation
}

setMethod("features",signature(object="CuffFeatureSet"),.features)

.fpkmMatrix<-function(object){
	res<-fpkm(object)
	colnames(res)[1]<-"tracking_id"
	res<-res[,c(1:3)]
	res<-melt(res)
	res<-cast(res,tracking_id~sample_name)
	res<-data.frame(res[,-1],row.names=res[,1])
}

setMethod("fpkmMatrix",signature(object="CuffFeatureSet"),.fpkmMatrix)

.diffData<-function(object){
	object@diff
}

setMethod("diffData",signature(object="CuffFeatureSet"),.diffData)


#################
#Plotting		#
#################
#Basic heatmap
.heatmap<-function(object,logMode=TRUE,pseudocount=0.0001){
	dat<-fpkm(object)
	if(logMode){
		dat$fpkm<- log10(dat$fpkm+pseudocount)
	}
	colnames(dat)[1] <- "tracking_id"
	p<-ggplot(dat)
	p <- p + geom_tile(aes(x=tracking_id,y=sample_name,fill=fpkm)) + scale_fill_gradient(low="white",high="red") + opts(axis.text.x=theme_text(angle=-90, hjust=0))
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
	
	m=fpkmMatrix(object)

	#remove genes with no expression in any condition
	m=m[!apply(m,1,sum)==0,]
	
	## you can either scale by row or column not both! 
	## if you wish to scale by both or use a different scale method then simply supply a scale
	## function instead NB scale is a base funct
	
	if(is.function(rescaling))
	{ 
		m=rescaling(m)
	} else {
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
	}else{
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
		melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), melt( log10(m+pseudocount)))
	}else{
		melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), melt(m))
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
	
	##adjust x-axis labels
	g2=g2+opts(axis.text.x=theme_text(angle=-90, hjust=0))
	
	## finally add the fill colour ramp of your choice (default is blue to red)-- and return
	return(g2+scale_fill_continuous("", heatscale[1], heatscale[2]))
	
}

setMethod("csHeatmap",signature("CuffFeatureSet"),.ggheat)

#Scatterplot
.scatter<-function(object,x,y,logMode=TRUE,pseudocount=0.0001,labels, smooth=FALSE,colorByStatus=FALSE,...){
	dat<-fpkmMatrix(object)
	samp<-samples(object)
	
	#check to make sure x and y are in samples
	if (!all(c(x,y) %in% samp)){
		stop("One or more values of 'x' or 'y' are not valid sample names!")
	}
	
	#add pseudocount if necessary
	if(logMode){
		for (i in samp){
			dat[[i]]<-dat[[i]]+pseudocount
		}
	}
	
	#make plot object
	p<-ggplot(dat)
	p<- p + aes_string(x=x,y=y)
	p<- p + geom_point(size=1.2,alpha=I(1/3)) + geom_abline(intercept=0,slope=1,linetype=2) + geom_rug(size=0.5,alpha=0.01)
	
	#add smoother
	if(smooth){
		p <- p + stat_smooth(method="lm",fill="blue",alpha=0.2)
	}
	
	#Add highlights from labels
#	if(!missing(labels)){
#		labelIdx<-fData(object)$gene_short_name %in% labels
#		labelfp<-fp[labelIdx,]
#		labelfp$gene_short_name<-fData(object)$gene_short_name[labelIdx]
#		#print(head(labelfp))
#		p <- p + geom_point(data=labelfp,size=1.2,color="red")
#		p <- p + geom_text(data=labelfp,aes(label=gene_short_name),color="red",hjust=0,vjust=0,angle=45,size=2)
#	}
#	
	#logMode
	if(logMode){
		p <- p + scale_y_log10() + scale_x_log10()
	}
	
	#Add title & Return value
	#p<- p + opts(title=object@tables$mainTable)
	p
}

setMethod("csScatter",signature(object="CuffFeatureSet"), .scatter)

#Volcano plot
.volcano<-function(object,x,y,xlimits=c(-20,20),...){
	dat<-diffData(object=object,x=x,y=y)
	s1<-unique(dat$sample_1)
	s2<-unique(dat$sample_2)
	
	p<-ggplot(dat)
	p<- p + geom_point(aes(x=ln_fold_change,y=-log10(p_value),color=significant),size=1,alpha=I(1/3))
	
	#Set axis limits
	p<- p + scale_x_continuous(limits=xlimits)
	p
}

setMethod("csVolcano",signature(object="CuffFeatureSet"), .volcano)

.barplot<-function(object,logMode=TRUE,pseudocount=0.0001,...){
	dat<-fpkm(object,features=T)
	#TODO: Test dat to ensure that there are >0 rows to plot.  If not, trap error and move on...
	
	colnames(dat)[1]<-"tracking_id"
	p<-ggplot(dat,aes(x=tracking_id,y=fpkm,fill=tracking_id))
	
	p<- p +
			geom_bar() +
			geom_errorbar(aes(ymin=conf_lo,ymax=conf_hi,group=1),size=0.15) +
			facet_wrap('sample_name') +
			opts(axis.text.x=theme_text(hjust=0,angle=-90))
	
	#This does not make immediate sense with the conf_hi and conf_lo values.  Need to figure out appropriate transformation for these
	#if(logMode)
	#p<-p+scale_y_ log10()
	p + opts(legend.position = "none")
	
}

setMethod("expressionBarplot",signature(object="CuffFeatureSet"),.barplot)

#################
#Clustering		#
#################
#Kmeans by expression profile using JSdist?
.cluster<-function(object,k,metric='euclidean',iter.max=100, ...){
	m=as.data.frame(fpkmMatrix(object))
	clusters<-kmeans(m,k,iter.max=iter.max)$cluster
	m$ids<-rownames(m)
	m$cluster<-factor(clusters)
	m.melt<-melt(m,id.vars=c("ids","cluster"))
	c<-ggplot(m.melt)
	c<-c+geom_line(aes(x=variable,y=value,color=cluster,group=ids)) + facet_wrap('cluster',scales='free')
	c
}

setMethod("csCluster",signature(object="CuffFeatureSet"),.cluster)
#################
#Misc			#
#################