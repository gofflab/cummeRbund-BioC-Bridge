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
				repFpkm=data.frame(),
				count=data.frame(),
				... ){
			.Object<-callNextMethod(.Object,
					annotation=annotation,
					fpkm=fpkm,
					diff=diff,
					repFpkm=repFpkm,
					count=count,
					...)
		}
)

#################
#Validate		#
#################
#TODO: Add validity constraints
setValidity("CuffFeatureSet",function(object){
			TRUE
			#Add test for genes with no expression
		}
)		

#################
#Class Methods	#
#################
setMethod("show","CuffFeatureSet",function(object){
		cat(class(object),"instance for ",length(object)," features\nSlots:
			\tannotation
			\tfpkm
			\trepFpkm
			\tdiff
			\tcount\n")
		}
)

setMethod("length","CuffFeatureSet",
		function(x){
			dim(x@annotation)[1]
		}
)

#Add dim method so you can find number of samples in CuffFeatureSet, not just number of genes


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

.replicates<-function(object){
	res<-repFpkm(object)$rep_name
	res
}

setMethod("replicates","CuffFeatureSet",.replicates)

.fpkm<-function(object,features=FALSE){
	myFPKM<-object@fpkm
	myFPKM$stdev<-(myFPKM$conf_hi-myFPKM$fpkm)/2
	if (features){
		return (merge(object@annotation,myFPKM,by=1))
	}else{
		return(myFPKM)
	}
}
setMethod("fpkm",signature(object="CuffFeatureSet"),.fpkm)

.repFpkm<-function(object,features=FALSE){
	myFPKM<-object@repFpkm
	#myFPKM$stdev<-(myFPKM$conf_hi-myFPKM$fpkm)/2
	if (features){
		return (merge(object@annotation,myFPKM,by=1))
	}else{
		return(myFPKM)
	}
}
setMethod("repFpkm",signature(object="CuffFeatureSet"),.repFpkm)

.featureNames<-function(object){
	data.frame(tracking_id=object@annotation[,1],gene_short_name=object@annotation$gene_short_name)
}

setMethod("featureNames",signature(object="CuffFeatureSet"),.featureNames)

.annotation<-function(object){
	object@annotation
}

setMethod(BiocGenerics::annotation,signature(object="CuffFeatureSet"),.annotation)

.fpkmMatrix<-function(object,fullnames=FALSE,sampleIdList){
	#Sample subsetting
	if(!missing(sampleIdList)){
		if (!all(sampleIdList %in% samples(object))){
			stop("Sample does not exist!")
		}else{
			mySamples<-sampleIdList
		}
	}else{
		mySamples<-samples(object)
	}
	if(fullnames){
		res<-fpkm(object,features=TRUE)
		res$tracking_id<-paste(res$gene_short_name,res[,1],sep="|")
	}else{
		res<-fpkm(object)
		colnames(res)[1]<-"tracking_id"	
	}
	selectedRows<-c('tracking_id','sample_name','fpkm')
	res<-res[,selectedRows]
	res<-melt(res)
	res<-dcast(res,tracking_id~sample_name)
	res<-data.frame(res[,-1],row.names=res[,1])
	if(!missing(sampleIdList)){
		res<-res[,mySamples]
	}
	res
}

setMethod("fpkmMatrix",signature(object="CuffFeatureSet"),.fpkmMatrix)

.repFpkmMatrix<-function(object,fullnames=FALSE,repIdList){
	#Sample subsetting
	if(!missing(repIdList)){
		if (!all(repIdList %in% replicates(object))){
			stop("Replicate does not exist!")
		}else{
			myReps<-repIdList
		}
	}else{
		myReps<-replicates(object)
	}
	if(fullnames){
		res<-repFpkm(object,features=TRUE)
		res$tracking_id<-paste(res$gene_short_name,res[,1],sep="|")
	}else{
		res<-repFpkm(object)
		colnames(res)[1]<-"tracking_id"	
	}
	selectedRows<-c('tracking_id','rep_name','fpkm')
	res<-res[,selectedRows]
	res<-melt(res)
	res<-dcast(res,tracking_id~rep_name)
	res<-data.frame(res[,-1],row.names=res[,1])
	if(!missing(repIdList)){
		res<-res[,myReps]
	}
	res
}

setMethod("repFpkmMatrix",signature(object="CuffFeatureSet"),.repFpkmMatrix)

.countMatrix<-function(object,fullnames=FALSE,sampleIdList){
	#Sample subsetting
	if(!missing(sampleIdList)){
		if (!all(sampleIdList %in% samples(object))){
			stop("Sample does not exist!")
		}else{
			mySamples<-sampleIdList
		}
	}else{
		mySamples<-samples(object)
	}
	if(fullnames){
		res<-count(object,features=TRUE)
		res$tracking_id<-paste(res$gene_short_name,res[,1],sep="|")
	}else{
		res<-count(object)
		colnames(res)[1]<-"tracking_id"	
	}
	selectedRows<-c('tracking_id','sample_name','count')
	res<-res[,selectedRows]
	res<-melt(res)
	res<-dcast(res,tracking_id~sample_name)
	res<-data.frame(res[,-1],row.names=res[,1])
	if(!missing(sampleIdList)){
		res<-res[,mySamples]
	}
	res
}

setMethod("countMatrix",signature(object="CuffFeatureSet"),.countMatrix)

.repCountMatrix<-function(object,fullnames=FALSE,repIdList){
	#Sample subsetting
	if(!missing(repIdList)){
		if (!all(repIdList %in% replicates(object))){
			stop("Replicate does not exist!")
		}else{
			myReps<-repIdList
		}
	}else{
		myReps<-replicates(object)
	}
	if(fullnames){
		res<-repFpkm(object,features=TRUE)
		res$tracking_id<-paste(res$gene_short_name,res[,1],sep="|")
	}else{
		res<-repFpkm(object)
		colnames(res)[1]<-"tracking_id"	
	}
	selectedRows<-c('tracking_id','rep_name','external_scaled_frags')
	res<-res[,selectedRows]
	res<-melt(res)
	res<-dcast(res,tracking_id~rep_name)
	res<-data.frame(res[,-1],row.names=res[,1])
	if(!missing(repIdList)){
		res<-res[,myReps]
	}
	res
}

setMethod("repCountMatrix",signature(object="CuffFeatureSet"),.repCountMatrix)

.diffData<-function(object){
	return(object@diff)
}

setMethod("diffData",signature(object="CuffFeatureSet"),.diffData)

.diffTable<-function(object,logCutoffValue=99999){
	measureVars<-c('status','value_1','value_2','log2_fold_change','test_stat','p_value','q_value','significant')
	all.diff<-diffData(object,features=TRUE)
	all.diff$log2_fold_change[all.diff$log2_fold_change>=logCutoffValue]<-Inf
	all.diff$log2_fold_change[all.diff$log2_fold_change<=-logCutoffValue]<--Inf
	all.diff.melt<-melt(all.diff,measure.vars=measureVars)
	#all.diff.melt<-all.diff.melt[!grepl("^value_",all.diff.melt$variable),]
	all.diff.cast<-dcast(all.diff.melt,formula=...~sample_2+sample_1+variable)
	all.diff.cast
}

setMethod("diffTable",signature(object="CuffFeatureSet"),.diffTable)


.count<-function(object,features=FALSE){
	if (features){
		return (merge(object@annotation,object@count,by=1))
	}else{
		return(object@count)
	}
}

setMethod("count",signature(object="CuffFeatureSet"),.count)

setMethod("annotation","CuffFeatureSet",function(object){
			return(object@annotation)
		})

#################
#Plotting		#
#################
#Basic heatmap
.heatmap<-function(object,logMode=TRUE,pseudocount=1.0){
	dat<-fpkm(object)
	if(logMode){
		dat$fpkm<- log10(dat$fpkm+pseudocount)
	}
	colnames(dat)[1] <- "tracking_id"
	p<-ggplot(dat)
	p <- p + geom_tile(aes(x=tracking_id,y=sample_name,fill=fpkm)) + scale_fill_gradient(low="white",high="red") + theme(axis.text.x=element_text(angle=-90, hjust=0))
	p
}

#######################
#The following is borrowed from Malarkey and is not yet ready for prime time...
#I would like to replace the clustering here with JSdistance on rows and/or columns and make
#this package work with cuffSet objects by default. 
#There is no genericMethod yet, goal is to replace .heatmap with .ggheat for genericMethod 'csHeatmap'

.ggheat<-function(object, rescaling='none', clustering='none', labCol=T, labRow=T, logMode=T, pseudocount=1.0, 
		border=FALSE, heatscale=c(low='lightyellow',mid='orange',high='darkred'), heatMidpoint=NULL,fullnames=T,replicates=FALSE,method='none',...) {
	## the function can be be viewed as a two step process
	## 1. using the rehape package and other funcs the data is clustered, scaled, and reshaped
	## using simple options or by a user supplied function
	## 2. with the now resahped data the plot, the chosen labels and plot style are built

	if(replicates){
		m=repFpkmMatrix(object,fullnames=fullnames)
	}else{
		m=fpkmMatrix(object,fullnames=fullnames)
	}
	#remove genes with no expression in any condition
	m=m[!apply(m,1,sum)==0,]
	
	## you can either scale by row or column not both! 
	## if you wish to scale by both or use a different scale method then simply supply a scale
	## function instead NB scale is a base funct
	
    if(logMode) 
    {
      m = log10(m+pseudocount)
    }
	
	
	if(is.function(rescaling))
	{ 
		m=rescaling(m)
	} else {
		if(rescaling=='column'){
			m=scale(m, center=T)
			m[is.nan(m)] = 0
		}
		if(rescaling=='row'){ 
			m=t(scale(t(m),center=T))
			m[is.nan(m)] = 0
		}
	}
	
	
	## I have supplied the default cluster and euclidean distance (JSdist) - and chose to cluster after scaling
	## if you want a different distance/cluster method-- or to cluster and then scale
	## then you can supply a custom function 
	
	if(!is.function(method)){
		method = function(mat){JSdist(makeprobs(t(mat)))}	
	}

	if(clustering=='row')
		m=m[hclust(method(m))$order, ]
	if(clustering=='column')  
		m=m[,hclust(method(t(m)))$order]
	if(clustering=='both')
		m=m[hclust(method(m))$order ,hclust(method(t(m)))$order]

	## this is just reshaping into a ggplot format matrix and making a ggplot layer
	

	rows=dim(m)[1]
	cols=dim(m)[2]
	
	
	
    # if(logMode) {
    #   melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), melt( log10(m+pseudocount)))
    # }else{
    #   melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), melt(m))
    # }
    


    melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), melt(m))

	g=ggplot(data=melt.m)
	
	## add the heat tiles with or without a white border for clarity
	
	if(border==TRUE)
		g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value),colour='grey')
	if(border==FALSE)
		g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value))
	
	## add axis labels either supplied or from the colnames rownames of the matrix
	
	if(labCol==T) 
	{
		g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=colnames(m))
	}
	if(labCol==F) 
	{
		g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=rep('',cols))
	}
	
	
	if(labRow==T) 
	{
		g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rownames(m))	
	}
	if(labRow==F)
	{ 
		g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rep('',rows))	
	}
	
	# Get rid of the ticks, they get way too dense with lots of rows
    g2 <- g2 + theme(axis.ticks = element_blank()) 

	## get rid of grey panel background and gridlines
	
	g2=g2+theme(panel.grid.minor=element_line(colour=NA), panel.grid.major=element_line(colour=NA),
			panel.background=element_rect(fill=NA, colour=NA))
	
	##adjust x-axis labels
	g2=g2+theme(axis.text.x=element_text(angle=-90, hjust=0))

    #write(paste(c("Length of heatscale is :", length(heatscale))), stderr())
	
	if (logMode)
	{
	   legendTitle <- bquote(paste(log[10]," FPKM + ",.(pseudocount),sep=""))
	   #legendTitle <- paste(expression(plain(log)[10])," FPKM + ",pseudocount,sep="")
	} else {
	   legendTitle <- "FPKM"
	}
	
	if (length(heatscale) == 2){
	    g2 <- g2 + scale_fill_gradient(low=heatscale[1], high=heatscale[2], name=legendTitle)
	} else if (length(heatscale) == 3) {
	    if (is.null(heatMidpoint))
	    {
	        heatMidpoint = (max(m) + min(m)) / 2.0
	        #write(heatMidpoint, stderr())
	    }

	    g2 <- g2 + scale_fill_gradient2(low=heatscale[1], mid=heatscale[2], high=heatscale[3], midpoint=heatMidpoint, name=legendTitle)
	}
	
	#g2<-g2+scale_x_discrete("",breaks=tracking_ids,labels=gene_short_names)
	
	
	## finally add the fill colour ramp of your choice (default is blue to red)-- and return
	return (g2)
	
}

setMethod("csHeatmap",signature("CuffFeatureSet"),.ggheat)


.fcheatmap<-function(object, control_condition, replicate_num=NULL, clustering='none', labCol=T, labRow=T, logMode=F, pseudocount=1.0, 
		border=FALSE, heatscale=c(low='steelblue',mid='white',high='tomato'), heatMidpoint=0,fullnames=T,replicates=FALSE,method='none',heatRange=3, ...) {
	## the function can be be viewed as a two step process
	## 1. using the rehape package and other funcs the data is clustered, scaled, and reshaped
	## using simple options or by a user supplied function
	## 2. with the now resahped data the plot, the chosen labels and plot style are built

	if(replicates){
	    if (is.null(replicate_num)){
	        print ("Error: if replicates == TRUE, you must specify both a control condition and a replicate number")
	        return()
	    }
		m=repFpkmMatrix(object,fullnames=fullnames)
		selected_rep <- paste(control_condition,replicate_num,sep="_")
		
		#remove genes with no expression in any condition
    	m=m[!apply(m,1,sum)==0,]
    	m=log2(m/m[,selected_rep])
		m=m[,names(m) != selected_rep]
	}else{
		m=fpkmMatrix(object,fullnames=fullnames)
		#remove genes with no expression in any condition
    	m=m[!apply(m,1,sum)==0,]
		m=log2(m/m[,control_condition])
		m=m[,names(m) != control_condition]
	}
	
	m_vals <- unlist(as.list(m))[is.finite(unlist(as.list(m)))]
	m_max <- max(m_vals)
	m_min <- max(m_vals)
	range_lim <- max(c(abs(m_max), abs(m_min)))
	range_lim <- min(c(heatRange, range_lim))
	m_max <- range_lim
	m_min <- -range_lim
	m[m <  m_min] <- m_min
	m[m >  m_max] <- m_max
	m[is.na(m)] <- 0
	
	## I have supplied the default cluster and euclidean distance (JSdist) - and chose to cluster after scaling
	## if you want a different distance/cluster method-- or to cluster and then scale
	## then you can supply a custom function 
	
	if(!is.function(method)){
		method = dist
	}

	if(clustering=='row')
		m=m[hclust(method(m))$order, ]
	if(clustering=='column')  
		m=m[,hclust(method(t(m)))$order]
	if(clustering=='both')
		m=m[hclust(method(m))$order ,hclust(method(t(m)))$order]
	
	rows=dim(m)[1]
	cols=dim(m)[2]

    melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows), melt(m))

	g=ggplot(data=melt.m)
	
	## add the heat tiles with or without a white border for clarity
	
	if(border==TRUE)
		g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value),colour='grey')
	if(border==FALSE)
		g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value))
	
	## add axis labels either supplied or from the colnames rownames of the matrix
	
	if(labCol==T) 
	{
		g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=colnames(m))
	}
	if(labCol==F) 
	{
		g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=rep('',cols))
	}
	
	
	if(labRow==T) 
	{
		g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rownames(m))	
	}
	if(labRow==F)
	{ 
		g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rep('',rows))	
	}
	
	# Get rid of the ticks, they get way too dense with lots of rows
    g2 <- g2 + theme(axis.ticks = element_blank()) 

	## get rid of grey panel background and gridlines
	
	g2=g2+theme(panel.grid.minor=element_line(colour=NA), panel.grid.major=element_line(colour=NA),
			panel.background=element_rect(fill=NA, colour=NA))
	
	##adjust x-axis labels
	g2=g2+theme(axis.text.x=element_text(angle=-90, hjust=0))

    #write(paste(c("Length of heatscale is :", length(heatscale))), stderr())
    if(replicates){
        legendTitle <- bquote(paste(log[2], frac("FPKM",.(selected_rep),sep="")))
    }else{
    	legendTitle <- bquote(paste(log[2], frac("FPKM",.(control_condition),sep="")))
	}
	#legendTitle <- paste(expression(plain(log)[10])," FPKM + ",pseudocount,sep="")

	if (length(heatscale) == 2){
	    g2 <- g2 + scale_fill_gradient(low=heatscale[1], high=heatscale[2], name=legendTitle)
	} else if (length(heatscale) == 3) {
	    if (is.null(heatMidpoint))
	    {
	        heatMidpoint = (max(m) + min(m)) / 2.0
	        #write(heatMidpoint, stderr())
	    }

	    g2 <- g2 + scale_fill_gradient2(low=heatscale[1], mid=heatscale[2], high=heatscale[3], midpoint=heatMidpoint, name=legendTitle)
	}
	
	#g2<-g2+scale_x_discrete("",breaks=tracking_ids,labels=gene_short_names)
	
	
	## finally add the fill colour ramp of your choice (default is blue to red)-- and return
	return (g2)
	
}

setMethod("csFoldChangeHeatmap",signature("CuffFeatureSet"),.fcheatmap)

# Distance Heatmaps
.distheat<-function(object, replicates=F, samples.not.genes=T, logMode=T, pseudocount=1.0, heatscale=c(low='lightyellow',mid='orange',high='darkred'), heatMidpoint=NULL, ...) {
  # get expression from a sample or gene perspective
	if(replicates){
		obj.fpkm<-repFpkmMatrix(object,fullnames=T)
	}else{
		obj.fpkm<-fpkmMatrix(object,fullnames=T)
	}

	if(samples.not.genes) {
    obj.fpkm.pos = obj.fpkm[rowSums(obj.fpkm)>0,]
  } else {
    obj.fpkm = t(obj.fpkm)
    obj.fpkm.pos = obj.fpkm[,colSums(obj.fpkm)>0]
  }
  
  if(logMode) {
    obj.fpkm.pos = log10(obj.fpkm.pos+pseudocount)
  }

  # compute distances
  obj.dists = JSdist(makeprobs(obj.fpkm.pos))
  
  # cluster to order
  obj.hc = hclust(obj.dists)

  # make data frame
  dist.df = melt(as.matrix(obj.dists),varnames=c("X1","X2"))

  # initialize
  g = ggplot(dist.df, aes(x=X1, y=X2, fill=value))

  # draw
  labels = labels(obj.dists)
  g = g + geom_tile() + scale_x_discrete("", limits=labels[obj.hc$order]) + scale_y_discrete("", limits=labels[obj.hc$order])

  # roll labels
  g = g + theme(axis.text.x=element_text(angle=-90, hjust=0), axis.text.y=element_text(angle=0, hjust=1))

  # drop grey panel background and gridlines
  g = g + theme(panel.grid.minor=element_line(colour=NA), panel.grid.major=element_line(colour=NA), panel.background=element_rect(fill=NA, colour=NA))

  # adjust heat scale
  if (length(heatscale) == 2) {
    g = g + scale_fill_gradient(low=heatscale[1], high=heatscale[3], name="JS Distance")
  }
  else if (length(heatscale) == 3) {
    if (is.null(heatMidpoint)) {
      heatMidpoint = max(obj.dists) / 2.0
    }
    g = g + scale_fill_gradient2(low=heatscale[1], mid=heatscale[2], high=heatscale[3], midpoint=heatMidpoint, name="JS Distance")
  }
  
  g <- g + geom_text(aes(label=format(value,digits=3)),size=6)
  
  # return
  g
}

setMethod("csDistHeat", signature("CuffFeatureSet"), .distheat)

#Scatterplot
.scatter<-function(object,x,y,logMode=TRUE,pseudocount=0.0,labels, smooth=FALSE,colorByStatus=FALSE,...){
	dat<-fpkmMatrix(object,fullnames=T)
	
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
	
	#Attach tracking_id and gene_short_name
	if(!missing(labels)){
		require(stringr)
		tracking<-str_split_fixed(rownames(dat),"\\|",2)
		dat$gene_short_name<-tracking[,1]
		dat$tracking_id<-tracking[,2]

		labeled.dat<-dat[dat$gene_short_name %in% labels,]
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
	if(!missing(labels)){
		p <- p + geom_point(data=labeled.dat,aes_string(x=x,y=y),size=1.3,color="red")
		p <- p + geom_text(data=labeled.dat,aes_string(x=x,y=y,label='gene_short_name'),color="red",hjust=0,vjust=0,angle=0,size=4)
	}
	
	#logMode
	if(logMode){
		p <- p + scale_y_log10() + scale_x_log10()
	}
	
	if (logMode)
    {
        p <- p + ylab(paste(y, "FPKM +",pseudocount))
        p <- p + xlab(paste(x, "FPKM +",pseudocount))
    } else {
        p <- p + ylab(paste(y, "FPKM"))
        p <- p + xlab(paste(x, "FPKM"))
    }
	
	#Add title & Return value
	#p<- p + labs(title=object@tables$mainTable)
	p
}

#.scatter(sigGenes,'P7_lincBrn1b_KO_Brain','P7_WT_Brain',labels=c('Arc'))

setMethod("csScatter",signature(object="CuffFeatureSet"), .scatter)

#Volcano plot
.volcano<-function(object,x,y,alpha=0.05,showSignificant=TRUE,xlimits=c(-20,20),...){
	samp<-samples(object)
	
	#check to make sure x and y are in samples
	if (!all(c(x,y) %in% samp)){
		stop("One or more values of 'x' or 'y' are not valid sample names!")
	}
	dat<-diffData(object=object)
	dat$significant<-'no'
	dat$significant[dat$q_value<=alpha]<-'yes'
	
	#subset dat for samples of interest
	mySamples<-c(x,y)
	dat<-dat[(dat$sample_1 %in% mySamples & dat$sample_2 %in% mySamples),]
	
	#Labels
	s1<-unique(dat$sample_1)
	s2<-unique(dat$sample_2)
	
	p<-ggplot(dat)
	if(showSignificant){
		p<- p + geom_point(aes(x=log2_fold_change,y=-log10(p_value),color=significant),alpha=I(1/3),size=1.2)
	}else{
		p<- p + geom_point(aes(x=log2_fold_change,y=-log10(p_value)),alpha=I(1/3),size=1.2)
	}
	
	p<- p + labs(title=paste(s2,"/",s1,sep=""))
	
	#Set axis limits
	p<- p + scale_x_continuous(limits=xlimits)
	
	#Default cummeRbund colorscheme
	p<-p + scale_colour_manual(values = c("black","red"))
	
	p
}

setMethod("csVolcano",signature(object="CuffFeatureSet"), .volcano)

.barplot<-function(object,logMode=TRUE,pseudocount=1.0,showErrorbars=TRUE,showStatus=TRUE,replicates=FALSE,...){
	quant_types<-c("OK","FAIL","LOWDATA","HIDATA","TOOSHORT")
	quant_types<-factor(quant_types,levels=quant_types)
	quant_colors<-c("black","red","blue","orange","green")
	names(quant_colors)<-quant_types
	
	dat<-fpkm(object,features=T)
	dat$warning<-""
	dat$warning[dat$quant_status!="OK"]<-"!"
	#TODO: Test dat to ensure that there are >0 rows to plot.  If not, trap error and move on...
	
	#Handle replicates
	if(replicates){
		repDat<-repFpkm(object)
		colnames(repDat)[1]<-"tracking_id"
	}
	
	colnames(dat)[1]<-"tracking_id"
	#tracking_ids<-dat$tracking_id
	obj_features <- annotation(object)
	tracking_ids <- obj_features[,1]
	
	gene_labels<-obj_features$gene_short_name
	#print(gene_labels)
	gene_labels[is.na(gene_labels)] = tracking_ids[is.na(gene_labels)]
	#print(gene_labels)
	
	dodge <- position_dodge(width=0.9) 
	
	if(logMode)
	{
	    dat$fpkm <- dat$fpkm + pseudocount
	    dat$conf_hi <- dat$conf_hi + pseudocount
	    dat$conf_lo <- dat$conf_lo + pseudocount
    }
	
	p<-ggplot(dat,aes(x=tracking_id,y=fpkm,fill=sample_name))
	p <- p + geom_bar(position=dodge,stat='identity')

	if (showErrorbars)
	{
	    p <- p +
		    geom_errorbar(aes(ymin=conf_lo,ymax=conf_hi),position=dodge,width=0.5)
	}
	
	if (logMode)
	{
	    p <- p + scale_y_log10()
		status_pos = 1
    }
	if(showStatus){
		if(logMode){
			p <- p + geom_text(aes(x=tracking_id,y=1,label=warning),position=dodge,stat='identity',color='red',vjust=1.5,size=3)
		}else{
			p <- p + geom_text(aes(x=tracking_id,y=0,label=warning),position=dodge,color='red',stat='identity',vjust=1.5,size=3)
		}
	}
	#gene_labels = dat$gene_short_name
	p <- p + scale_x_discrete("",breaks=tracking_ids,labels=gene_labels) + 
	    theme(axis.text.x=element_text(hjust=0,angle=-90))
    	
    # p<- p +
    #       geom_bar() +
    #       geom_errorbar(aes(ymin=conf_lo,ymax=conf_hi,group=1),size=0.15) +
    #       facet_wrap('sample_name') +
    #       theme(axis.text.x=element_text(hjust=0,angle=-90))
	
	#This does not make immediate sense with the conf_hi and conf_lo values.  Need to figure out appropriate transformation for these
	#if(logMode)
	#p<-p+scale_y_ log10()
	
    if (logMode)
    {
        p <- p + ylab(paste("FPKM +",pseudocount))
    } else {
        p <- p + ylab("FPKM")
    }
	
	#p <- p + theme(legend.position = "none")
	
	#Default cummeRbund colorscheme
	p<-p + scale_fill_hue(l=50,h.start=200) + scale_color_hue(l=50,h.start=200)
	
	#Recolor quant flags
	p<- p+ scale_colour_manual(name='quant_status',values=quant_colors)
	
	p
	
}

setMethod("expressionBarplot",signature(object="CuffFeatureSet"),.barplot)

.expressionPlot<-function(object,logMode=FALSE,pseudocount=1.0, drawSummary=FALSE, sumFun=mean_cl_boot, showErrorbars=TRUE,showStatus=TRUE,replicates=FALSE,...){
	quant_types<-c("OK","FAIL","LOWDATA","HIDATA","TOOSHORT")
	quant_types<-factor(quant_types,levels=quant_types)
	quant_colors<-c("black","red","blue","orange","green")
	names(quant_colors)<-quant_types
	
	dat<-fpkm(object)
	colnames(dat)[1]<-"tracking_id"
	
	if(replicates){
		repDat<-repFpkm(object)
		repDat$replicate<-as.factor(repDat$replicate)
		colnames(repDat)[1]<-"tracking_id"
	}
	
	if(logMode)
	{
		dat$fpkm <- dat$fpkm + pseudocount
		dat$conf_hi <- dat$conf_hi + pseudocount
		dat$conf_lo <- dat$conf_lo + pseudocount
		
		if(replicates){
			repDat$fpkm<-repDat$fpkm + pseudocount
		}
	}
	p <- ggplot(dat)
	#dat$fpkm<- log10(dat$fpkm+pseudocount)
	p <- p + geom_line(aes(x=sample_name,y=fpkm,group=tracking_id,color=tracking_id))
	
	if(replicates){
		p <- p + geom_point(aes(x=sample_name,y=fpkm,color=tracking_id),size=2.5,shape=18,data=repDat)
	}
	
	if (showErrorbars)
	{
		p <- p +
				geom_errorbar(aes(x=sample_name, ymin=conf_lo,ymax=conf_hi,color=tracking_id,group=tracking_id),width=0.25)
	}
	
	if(showStatus){
		p <- p+ geom_point(aes(x=sample_name,y=fpkm,shape=quant_status))
	}
	
	if (logMode)
	{
		p <- p + scale_y_log10()
	}
	
	#drawMean
	if(drawSummary){
		p <- p + stat_summary(aes(x=sample_name,y=fpkm,group=1),fun.data=sumFun,color="red",fill="red",alpha=0.2,size=1.1,geom="smooth")
	}
	
	if (logMode)
	{
		p <- p + ylab(paste("FPKM +",pseudocount))
	} else {
		p <- p + ylab("FPKM")
	}
	
	#Recolor quant flags
	#p<- p + scale_colour_manual(name='quant_status',values=quant_colors)
	
	p
}

setMethod("expressionPlot",signature(object="CuffFeatureSet"),.expressionPlot)

#TODO: Add csDensity plot for CuffFeatureSet objects

#TODO: Add ecdf plot for CuffFeatureSet and CuffData objects

#################
#Clustering		#
#################
#Kmeans by expression profile using JSdist
#TODO:Make this function return an object of type CuffClusterSet
#.cluster<-function(object, k, pseudocount=1, ...){
#	library(cluster)
#	m<-as.data.frame(fpkmMatrix(object))
#	m<-m[rowSums(m)>0,]
#	n<-JSdist(makeprobs(t(m)))
#	clusters<-pam(n,k)
#	clusters$fpkm<-m
#	m<-m+pseudocount
#	m$ids<-rownames(m)
#	m$cluster<-factor(clusters$clustering)
#	m.melt<-melt(m,id.vars=c("ids","cluster"))
#	c<-ggplot(m.melt)
#	c<-c+geom_line(aes(x=variable,y=value,color=cluster,group=ids)) + facet_wrap('cluster',scales='free')+scale_y_log10()
#	
#	#Default cummeRbund colorscheme
#	c<-c + scale_color_hue(l=50,h.start=200)
#	c
#}

.cluster<-function(object, k, logMode=T, method='none', pseudocount=1,...){
	require(cluster)
	m<-as.data.frame(fpkmMatrix(object))
	m<-m[rowSums(m)>0,]
	if(logMode){
		m<-log10(m+pseudocount)
	}
	
	if(!is.function(method)){
		method = function(mat){JSdist(makeprobs(t(m)))}	
	}		
	n<-method(m)
	clusters<-pam(n,k, ...)
	#clsuters<-pamk(n,krange=2:20)
	class(clusters)<-"list"
	clusters$fpkm<-m
	clusters
}

setMethod("csCluster",signature(object="CuffFeatureSet"),.cluster)

.ratioCluster<-function(object,k,ratioTo=NULL,pseudocount=0.0001,...){
	require(cluster)
	m<-as.data.frame(fpkmMatrix(object)+pseudocount)
	#TODO: ensure that ratioTo is in colnames(fpkmMatrix(object))
	m.ratio<-m/m[[ratioTo]]
	m.log.ratio<-log2(m.ratio)
	n<-dist(m.log.ratio)
	clusters<-pam(n,k,...)
	class(clusters)<-"list"
	clusters$fpkm<-m
	clusters
}

csClusterPlot<-function(clustering,pseudocount=1.0,logMode=FALSE,drawSummary=TRUE, sumFun=mean_cl_boot){
	m<-clustering$fpkm+pseudocount
	m$ids<-rownames(clustering$fpkm)
	m$cluster<-factor(clustering$clustering)
	m.melt<-melt(m,id.vars=c("ids","cluster"))
	c<-ggplot(m.melt)
	c<-c+geom_line(aes(x=variable,y=value,color=cluster,group=ids)) + theme_bw() + facet_wrap('cluster',scales='free_y')
	if(drawSummary){
		c <- c + stat_summary(aes(x=variable,y=value,group=1),fun.data=sumFun,color="black",fill="black",alpha=0.2,size=1.1,geom="smooth")
	}
	if(logMode){
		c<-c + scale_y_log10()
	}
	c<-c + scale_color_hue(l=50,h.start=200) + theme(axis.text.x=element_text(angle=-90,hjust=0))
	c
}

.findK<-function(object, k.range=c(2:20), logMode=T, pseudocount=1,...){
	require(cluster)
	m<-as.data.frame(fpkmMatrix(object))
	m<-m[rowSums(m)>0,]
	if(logMode){
		m<-log10(m+pseudocount)
	}
	n<-JSdist(makeprobs(t(m)))
	myWidths<-c()
	for (k in k.range){
		print(k)
		myWidths<-c(myWidths,pam(n,k,...)$silinfo$avg.width)
	}
	plot(k.range,myWidths)
}
	
	
	

######################
# 
######################

.dendro<-function(object,logMode=T,pseudocount=1,replicates=FALSE){
	if(replicates){
		fpkmMat<-repFpkmMatrix(object)
	}else{
		fpkmMat<-fpkmMatrix(object)
	}
	if(logMode){
		fpkmMat<-log10(fpkmMat+pseudocount)
	}
	res<-JSdist(makeprobs(fpkmMat))
	#colnames(res)<-colnames(fpkmMat)
	
	#res<-as.dist(res)
	res<-as.dendrogram(hclust(res))
	plot(res)
	res
}

setMethod("csDendro",signature(object="CuffFeatureSet"),.dendro)

#.csNMF<-function(object,rank,method,...){
#	#Non-negative Matrix Factorization
#	require(NMF)
#	#Retrieve FPKM values
#	fpkms<-fpkmMatrix(object)
#	
#	#Remove rows with sum()==0
#	fpkms<-fpkms[rowSums(fpkms)>0,]
#	
#}
#
#setMethod("csNMF",signature(object="CuffFeatureSet"),.csNMF)

##Takes as first argument the object returned from csCluster (a modified 'cluster' list)
#.clusterPlot<-function(clusters, pseudocount=1, ...){
#	m<-clusters$fpkm
#	m<-m+pseudocount
#	m$ids<-rownames(m)
#	m$cluster<-factor(clusters$clustering)
#	m.melt<-melt(m,id.vars=c("ids","cluster"))
#	c<-ggplot(m.melt)
#	c<-c+geom_line(aes(x=variable,y=value,color=cluster,group=ids)) + facet_wrap('cluster',scales='free')+scale_y_log10()
#	
#	#Default cummeRbund colorscheme
#	c<-c + scale_color_hue(l=50,h.start=200)
#	c
#}

######################
# Exploratory Analysis
######################
.nmf<-function(object,k,logMode=T,pseudocount=1,maxiter=1000,replicates=FALSE,fullnames=FALSE){
	require(NMFN)
	if(missing(k)) stop("Please provide a rank value for factorization (arg=k)")
	
	if(replicates){
		m=repFpkmMatrix(object,fullnames=fullnames)
	}else{
		m=fpkmMatrix(object,fullnames=fullnames)
	}
	
	if(logMode) 
	{
		m = log10(m+pseudocount)
	}
	
	myNMF<-nnmf(m,k=k,maxiter=maxiter)
	return (myNMF)
}

setMethod("csNMF",signature(object="CuffFeatureSet"),.nmf)

.density<-function(object, logMode = TRUE, pseudocount=0.0, labels, features=FALSE, replicates=FALSE,...){
	if(replicates){
		dat<-repFpkm(object,features=features)
		colnames(dat)[colnames(dat)=="rep_name"]<-"condition"
	}else{
		dat<-fpkm(object,features=features)
		colnames(dat)[colnames(dat)=="sample_name"]<-"condition"
	}
	if(logMode) dat$fpkm<-dat$fpkm+pseudocount
	p<-ggplot(dat)
	if(logMode) {
		p<-p+geom_density(aes(x= log10(fpkm),group=condition,color=condition,fill=condition),alpha=I(1/3))
	}else{
		p<-p+geom_density(aes(x=fpkm,group=condition,color=condition,fill=condition),alpha=I(1/3))
	}
	
	#p<-p + labs(title=object@tables$mainTable)
	
	#Default cummeRbund colorscheme
	p<-p + scale_fill_hue(l=50,h.start=200) + scale_color_hue(l=50,h.start=200)
	
	#TODO: Add label callout
	p
}

setMethod("csDensity",signature(object="CuffFeatureSet"),.density)

#################
# Dimensionality Reduction (use sparingly for Gene Sets)
#################

.MDSplot<-function(object,replicates=FALSE,logMode=TRUE,pseudocount=1.0){
	if(replicates){
		dat<-repFpkmMatrix(object)
		#repData<-sapply(replicates(object),function(x){strsplit(x,"_")[[1]][1]}) This is to color by condition and not replicate...
	}else{
		dat<-fpkmMatrix(object)
	}
	
	if(logMode){
		dat<-log10(dat+pseudocount)
	}
	
	d<-JSdist(makeprobs(dat))
	fit <- cmdscale(d,eig=TRUE, k=2)
	res<-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
	p <- ggplot(res)
	p <- p + geom_point(aes(x=M1,y=M2,color=names)) + geom_text(aes(x=M1,y=M2,label=names,color=names)) + theme_bw()
	p
}

setMethod("MDSplot",signature(object="CuffFeatureSet"),.MDSplot)

.PCAplot<-function(object,x="PC1", y="PC2",replicates=FALSE,pseudocount=1.0,scale=TRUE,showPoints=TRUE,...){
	if(replicates){
		fpkms<-repFpkmMatrix(object)
	}else{
		fpkms<-fpkmMatrix(object)
	}
	fpkms<-log10(fpkms+pseudocount)
	PC<-prcomp(fpkms,scale=scale,...)
	dat <- data.frame(obsnames=row.names(PC$x), PC$x)
	#dat$shoutout<-""
	#dat$shoutout[matchpt(PC$rotation,PC$x)$index]<-rownames(pca$x[matchpt(pca$rotation,pca$x)$index,])
	plot <- ggplot(dat, aes_string(x=x, y=y)) 
	if(showPoints){
		plot<- plot + geom_point(alpha=.4, size=0.8, aes(label=obsnames))
	}
	plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2) #+ geom_text(aes(label=shoutout),size=2,color="red")
	datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
	mult <- min(
			(max(dat[,y]) - min(dat[,y])/(max(datapc[,y])-min(datapc[,y]))),
			(max(dat[,x]) - min(dat[,x])/(max(datapc[,x])-min(datapc[,x])))
	)
	datapc <- transform(datapc,
			v1 = .7 * mult * (get(x)),
			v2 = .7 * mult * (get(y))
	)
	plot <- plot + 
			#coord_equal() + 
			geom_text(data=datapc, aes(x=v1, y=v2, label=varnames,color=varnames), vjust=1)
	plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2,color=varnames), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75) + theme_bw()
	plot
}

setMethod('PCAplot',signature(object="CuffFeatureSet"),.PCAplot)

#################
#Misc			#
#################
.makeIdentityMatrix<-function(sampleNames){
	d<-diag(length(sampleNames))
}

.specificity<-function(object,logMode=T,pseudocount=1,relative=FALSE,...){
	fpkms<-fpkmMatrix(object,...)
	if(logMode){
		fpkms<-log10(fpkms+pseudocount)
	}
	fpkms<-t(makeprobs(t(fpkms)))
	d<-diag(ncol(fpkms))
	res<-apply(d,MARGIN=1,function(q){
				JSdistFromP(fpkms,q)
			})
	colnames(res)<-paste(colnames(fpkms),"_spec",sep="")
	
	if(relative){
		res<-res/max(res)
	}
	1-res
}

setMethod("csSpecificity",signature(object="CuffFeatureSet"),.specificity)
