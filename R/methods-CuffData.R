#methods-CuffData.R
#
#Author: Loyal A. Goff
#
#
####################

##################
#Initialize
##################
setMethod("initialize","CuffData",
			function(.Object,
					DB,
					tables=list(mainTable = "",dataTable = "",expDiffTable = "",featureTable = "",countTable = "", replicateTable = ""),
					filters=list(),
					type = c("genes","isoforms","TSS","CDS"),
					idField,
					... ){
				.Object<-callNextMethod(.Object,
						DB = DB,
						tables = tables,
						filters = filters,
						type = type,
						idField = idField,
						...)				
		}
)

setValidity("CuffData",function(object){
		TRUE
		}
)			

################
#Class Methods
################
setMethod("show","CuffData",
		function(object){
			size<-dim(object)
			cat(class(object), "instance with:\n\t",size[1],"features and",size[2],"samples\n")
		}
)

setMethod("dim","CuffData",
		function(x){
			countQuery<-paste("SELECT COUNT(",x@idField,") as n FROM ",x@tables$mainTable)
			nIds<-dbGetQuery(x@DB,countQuery)
			sampleQuery<-("SELECT COUNT(sample_name) as n FROM samples")
			nSamples<-dbGetQuery(x@DB,sampleQuery)
			c(nIds$n,nSamples$n)
		}
)

#####################
#Feature Table
#####################

.addFeatures<-function(object,features,...){
	if(!is.data.frame(features)){
		stop("features must be a data.frame")
	}
	colnames(features)[1]<-object@idField
	colnames(features)<-make.db.names(object@DB,colnames(features),unique=T)
	dbWriteTable(object@DB,object@tables$featureTable,features,row.names=F,overwrite=T)
	indexQuery<-paste("CREATE INDEX ",object@idField," ON ", object@tables$featureTable," (",object@idField,")",sep="")
	res<-dbGetQuery(object@DB,indexQuery)
}

setMethod("addFeatures",signature="CuffData",.addFeatures)

###################
#Accessors
###################
.annotation<-function(object){
	featureQuery<-paste("SELECT * FROM ",object@tables$mainTable," x LEFT JOIN ",object@tables$featureTable," xf ON x.",object@idField,"=xf.",object@idField,sep="")
	dbGetQuery(object@DB, featureQuery)
}

setMethod("annotation","CuffData",.annotation)

.featureNames<-function(object){
	featureQuery<-paste("SELECT ",object@idField," FROM ",object@tables$mainTable, sep="")
	res<-dbGetQuery(object@DB,featureQuery)
	res[,1]
}

setMethod("featureNames","CuffData",.featureNames)

.samples<-function(object){
	res<-dbReadTable(object@DB,'samples')
	res<-res$sample_name
	res
}

setMethod("samples","CuffData",.samples)

.replicates<-function(object){
	res<-dbReadTable(object@DB,'replicates')
	res<-res$rep_name
	res
}

setMethod("replicates","CuffData",.replicates)

.fpkm<-function(object,features=FALSE,sampleIdList){
	#Sample subsetting
	if(!missing(sampleIdList)){
		if(.checkSamples(object@DB,sampleIdList)){
			myLevels<-sampleIdList
		}else{
			stop("Sample does not exist!")
		}
	}else{
		myLevels<-getLevels(object)
	}
	#Sample Search String (SQL)
	sampleString<-'('
	for (i in myLevels){
		sampleString<-paste(sampleString,"'",i,"',",sep="")
	}
	sampleString<-substr(sampleString,1,nchar(sampleString)-1)
	sampleString<-paste(sampleString,")",sep="")
	
	if(!features){
		FPKMQuery<-paste("SELECT * FROM ",object@tables$dataTable," WHERE sample_name IN ",sampleString,sep="")
	}else{
		FPKMQuery<-paste("SELECT xf.*,xm.*,x.sample_name,x.fpkm,x.conf_hi,x.conf_lo FROM ",object@tables$dataTable," x LEFT JOIN ",object@tables$featureTable," xf ON x.",object@idField,"=xf.",object@idField," LEFT JOIN ",object@tables$mainTable," xm ON x.",object@idField,"=xm.",object@idField," WHERE x.sample_name IN ",sampleString,sep="")
		#print(FPKMQuery)
	}
	res<-dbGetQuery(object@DB,FPKMQuery)
	res$sample_name<-factor(res$sample_name,levels=getLevels(object))
	res$stdev<-(res$conf_hi-res$fpkm)/2
	res
}

setMethod("fpkm","CuffData",.fpkm)

.repFpkm<-function(object,features=FALSE,repIdList){
	#Sample subsetting
	if(!missing(repIdList)){
		if(.checkReps(object@DB,repIdList)){
			myLevels<-repIdList
		}else{
			stop("Replicate does not exist!")
		}
	}else{
		myLevels<-getRepLevels(object)
	}
	#Sample Search String (SQL)
	sampleString<-'('
	for (i in myLevels){
		sampleString<-paste(sampleString,"'",i,"',",sep="")
	}
	sampleString<-substr(sampleString,1,nchar(sampleString)-1)
	sampleString<-paste(sampleString,")",sep="")
	
	if(!features){
		FPKMQuery<-paste("SELECT * FROM ",object@tables$replicateTable," WHERE rep_name IN ",sampleString,sep="")
	}else{
		FPKMQuery<-paste("SELECT xf.*,xm.*,x.rep_name,x.raw_frags,x.internal_scaled_frags,x.external_scaled_frags,x.fpkm,x.effective_length,x.status FROM ",object@tables$replicateTable," x LEFT JOIN ",object@tables$featureTable," xf on x.",object@idField,"=xf.",object@idField," LEFT JOIN ",object@tables$mainTable," xm ON x.",object@idField,"=xm.",object@idField," WHERE x.rep_name IN ",sampleString,sep="")
	}
	#print(FPKMQuery)
	res<-dbGetQuery(object@DB,FPKMQuery)
	res$rep_name<-factor(res$rep_name,levels=getRepLevels(object))
	#res$stdev<-(res$conf_hi-res$fpkm)/2 #Not really available yet since conf_hi and conf_lo are not 
	res
}

setMethod("repFpkm","CuffData",.repFpkm)

.count<-function(object,sampleIdList){
	#Sample subsetting
	if(!missing(sampleIdList)){
		if(.checkSamples(object@DB,sampleIdList)){
			myLevels<-sampleIdList
		}else{
			stop("Sample does not exist!")
		}
	}else{
		myLevels<-getLevels(object)
	}
	#Sample Search String (SQL)
	sampleString<-'('
	for (i in myLevels){
		sampleString<-paste(sampleString,"'",i,"',",sep="")
	}
	sampleString<-substr(sampleString,1,nchar(sampleString)-1)
	sampleString<-paste(sampleString,")",sep="")
	
	CountQuery<-FPKMQuery<-paste("SELECT * FROM ",object@tables$countTable," WHERE sample_name IN ",sampleString,sep="")
	
	res<-dbGetQuery(object@DB,CountQuery)
	res$sample_name<-factor(res$sample_name,levels=getLevels(object))
	res
	
}

setMethod("count","CuffData",.count)

.fpkmMatrix<-function(object,fullnames=FALSE,sampleIdList){
	#TODO: fix fullnames for CuffData::fpkmMatrix
	#Sample subsetting
	if(!missing(sampleIdList)){
		if(.checkSamples(object@DB,sampleIdList)){
			myLevels<-sampleIdList
		}else{
			stop("Sample does not exist!")
		}
	}else{
		myLevels<-getLevels(object)
	}
	
	samp<-samples(object)
	FPKMMatQuery<-paste("select x.",object@idField,", ",sep="")
	for (i in samp){
		FPKMMatQuery<-paste(FPKMMatQuery,"sum(case when xd.sample_name ='",i,"' then fpkm end) as ",i,",",sep="")
	}
	FPKMMatQuery<-substr(FPKMMatQuery, 1, nchar(FPKMMatQuery)-1)
	FPKMMatQuery<-paste(FPKMMatQuery," from ",object@tables$mainTable," x LEFT JOIN ",object@tables$dataTable," xd on x.",object@idField," = xd.",object@idField," group by x.",object@idField,sep="")
	res<-dbGetQuery(object@DB,FPKMMatQuery)
	res<-data.frame(res[,-1],row.names=res[,1])
	if(!missing(sampleIdList)){
		res<-data.frame(res[,sampleIdList],row.names=rownames(res))
		colnames(res)<-sampleIdList
	}
	res
}

setMethod("fpkmMatrix","CuffData",.fpkmMatrix)

.repFpkmMatrix<-function(object,fullnames=FALSE,repIdList){
	#Sample subsetting
	if(!missing(repIdList)){
		if(.checkReps(object@DB,repIdList)){
			myLevels<-repIdList
		}else{
			stop("Replicate does not exist!")
		}
	}else{
		myLevels<-getRepLevels(object)
	}
	
	samp<-replicates(object)
	FPKMMatQuery<-paste("select x.",object@idField,", ",sep="")
	for (i in samp){
		FPKMMatQuery<-paste(FPKMMatQuery,"sum(case when xd.rep_name ='",i,"' then fpkm end) as ",i,",",sep="")
	}
	FPKMMatQuery<-substr(FPKMMatQuery, 1, nchar(FPKMMatQuery)-1)
	FPKMMatQuery<-paste(FPKMMatQuery," from ",object@tables$mainTable," x LEFT JOIN ",object@tables$replicateTable," xd on x.",object@idField," = xd.",object@idField," group by x.",object@idField,sep="")
	res<-dbGetQuery(object@DB,FPKMMatQuery)
	res<-data.frame(res[,-1],row.names=res[,1])
	if(!missing(repIdList)){
		res<-data.frame(res[,repIdList],row.names=rownames(res))
		colnames(res)<-repIdList
	}
	res
}

setMethod("repFpkmMatrix","CuffData",.repFpkmMatrix)

.countMatrix<-function(object,fullnames=FALSE,sampleIdList){
	#Sample subsetting
	if(!missing(sampleIdList)){
		if(.checkSamples(object@DB,sampleIdList)){
			myLevels<-sampleIdList
		}else{
			stop("Sample does not exist!")
		}
	}else{
		myLevels<-getLevels(object)
	}
	
	samp<-samples(object)
	CountMatQuery<-paste("select x.",object@idField,", ",sep="")
	for (i in samp){
		CountMatQuery<-paste(CountMatQuery,"sum(case when xd.sample_name ='",i,"' then count end) as ",i,",",sep="")
	}
	CountMatQuery<-substr(CountMatQuery, 1, nchar(CountMatQuery)-1)
	CountMatQuery<-paste(CountMatQuery," from ",object@tables$mainTable," x LEFT JOIN ",object@tables$countTable," xd on x.",object@idField," = xd.",object@idField," group by x.",object@idField,sep="")
	res<-dbGetQuery(object@DB,CountMatQuery)
	res<-data.frame(res[,-1],row.names=res[,1])
	if(!missing(sampleIdList)){
		res<-data.frame(res[,sampleIdList],row.names=rownames(res))
		colnames(res)<-sampleIdList
	}
	res
}

setMethod("countMatrix","CuffData",.countMatrix)

.repCountMatrix<-function(object,fullnames=FALSE,repIdList){
	#Sample subsetting
	if(!missing(repIdList)){
		if(.checkReps(object@DB,repIdList)){
			myLevels<-repIdList
		}else{
			stop("Replicate does not exist!")
		}
	}else{
		myLevels<-getRepLevels(object)
	}
	reps<-replicates(object)
	repCountMatQuery<-paste("select x.",object@idField,", ",sep="")
	for (i in reps){
		repCountMatQuery<-paste(repCountMatQuery,"sum(case when xr.rep_name ='",i,"' then external_scaled_frags end) as ",i,",",sep="")
	}
	repCountMatQuery<-substr(repCountMatQuery, 1, nchar(repCountMatQuery)-1)
	repCountMatQuery<-paste(repCountMatQuery," from ",object@tables$mainTable," x LEFT JOIN ",object@tables$replicateTable," xr on x.",object@idField," = xr.",object@idField," group by x.",object@idField,sep="")
	res<-dbGetQuery(object@DB,repCountMatQuery)
	res<-data.frame(res[,-1],row.names=res[,1])
	if(!missing(repIdList)){
		res<-data.frame(res[,repIdList],row.names=rownames(res))
		colnames(res)<-repIdList
	}
	res
}

setMethod("repCountMatrix","CuffData",.repCountMatrix)

.statsMatrix<-function(object){
	statsQuery<-paste("SELECT xd.*, xc.count, xc.variance as count_variance , xc.uncertainty as count_uncertainty, xc.dispersion as count_dispersion, (xd.conf_hi-xd.fpkm)/2 as fpkm_stdev,((xd.conf_hi-xd.fpkm)/2)/xd.fpkm AS 'CV' FROM ",object@tables$dataTable," xd LEFT JOIN ",object@tables$countTable," xc ON xd.",object@idField,"=xc.",object@idField," AND xd.sample_name=xc.sample_name",sep="")
	res<-dbGetQuery(object@DB,statsQuery)
	res$sample_name<-factor(res$sample_name,levels=samples(object))
	res
}

#This needs a lot of work...
#TODO: Change this to remove lnFcCutoff but make sure that functions that rely on diffData have their own FC cutoff so that plotting doesn't suffer
.diffData<-function(object,x,y,features=FALSE){
	if(missing(x) && missing(y)){
		if(!features){
			diffQuery<-paste("SELECT * FROM ",object@tables$expDiffTable,sep="")
		}else{
			diffQuery<-paste("SELECT xm.*, xed.*, xf.* FROM ",object@tables$mainTable," xm LEFT JOIN ",object@tables$expDiffTable," xed ON xm.",object@idField,"=xed.",object@idField," LEFT JOIN ",object@tables$featureTable," xf ON xm.",object@idField,"=xf.",object@idField,sep="")
		}
	}else if (missing(x) || missing(y)){
		stop("You must supply both x and y or neither.")
	}else{
		if(!features){
			diffQuery<-paste("SELECT x.",object@idField,", xed.* FROM ",object@tables$mainTable," x LEFT JOIN ",object@tables$expDiffTable," xed on x.",object@idField," = xed.",object@idField," WHERE ((sample_1 = '",x,"' AND sample_2 = '",y,"') OR (sample_1 = '",y,"' AND sample_2 = '",x,"'))",sep="")
		}else{
			diffQuery<-paste("SELECT xm.*, xed.*, xf.* FROM ",object@tables$mainTable," xm LEFT JOIN ",object@tables$expDiffTable," xed on xm.",object@idField," = xed.",object@idField," LEFT JOIN ",object@tables$featureTable," xf ON xm.",object@idField,"=xf.",object@idField," WHERE ((sample_1 = '",x,"' AND sample_2 = '",y,"') OR (sample_1 = '",y,"' AND sample_2 = '",x,"'))",sep="")
		}
	}
	dat<-dbGetQuery(object@DB,diffQuery)
	#diffQuery
	dat
}

setMethod("diffData",signature(object="CuffData"),.diffData)

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

setMethod("diffTable",signature(object="CuffData"),.diffTable)

.QCTable<-function(object){
	qcSQL<-paste("SELECT d.*, (d.conf_hi-d.fpkm)/2 as stdev, c.count, c.variance, c.uncertainty, c.dispersion, (c.variance/d.fpkm) AS 'IOD', ((d.conf_hi-d.fpkm)/2)/d.fpkm AS 'CV' FROM ",object@tables$dataTable," d LEFT JOIN ",object@tables$countTable," c ON d.gene_id=c.gene_id AND d.sample_name=c.sample_name",sep="")
	res<-dbGetQuery(object@DB,qcSQL)
	res
}

.getMA<-function(object,x,y,logMode=T,pseudocount=1){
	if (missing(x) || missing(y)){
		stop("You must supply both x and y.")
	}else{
		sql<-paste("SELECT x.",object@idField,", sum(case when x.sample_name = '",x,"' then x.fpkm end) AS 'x', sum(case when x.sample_name = '",y,"' then x.fpkm end) AS 'y' FROM ",object@tables$dataTable," x GROUP BY x.",object@idField,";",sep="")
		#print(sql)
		dat<-dbGetQuery(object@DB,sql)
		
		if(logMode){
			dat$x<-log10(dat$x+pseudocount)
			dat$y<-log10(dat$y+pseudocount)
		}
		dat$A<-(dat$x+dat$y)/2
		dat$M<-dat$x/dat$y
		res<-dat[,c(1,4:5)]
		res
	}
}

.getCountMA<-function(object,x,y,logMode=T,pseudocount=1){
	if (missing(x) || missing(y)){
		stop("You must supply both x and y.")
	}else{
		sql<-paste("SELECT x.",object@idField,", sum(case when x.sample_name = '",x,"' then x.count end) AS 'x', sum(case when x.sample_name = '",y,"' then x.count end) AS 'y' FROM ",object@tables$countTable," x GROUP BY x.",object@idField,";",sep="")
		dat<-dbGetQuery(object@DB,sql)
		
		if(logMode){
			dat$x<-log10(dat$x+pseudocount)
			dat$y<-log10(dat$y+pseudocount)
		}
		dat$A<-(dat$x+dat$y)/2
		dat$M<-dat$x/dat$y
		res<-dat[,c(1,4:5)]
		res
	}
}

#.getRankOrder<-function(object,x,y,logMode=TRUE,pseudocount=1,ratio=TRUE){
#	if (missing(x) || missing(y)){
#		stop("You must supply both x and y.")
#	}else{
#		
#	}
#}

setMethod("DB","CuffData",function(object){
		return(object@DB)
		})

setMethod("tables","CuffData",function(object){
		return(object@tables)
		})

setMethod("filters","CuffData",function(object){
		return(object@filters)
		})

setMethod("type","CuffData",function(object){
		return(object@type)
		})

setMethod("idField","CuffData",function(object){
		return(object@idField)
		})

##################
#Setters
##################


##################
#Subsetting
##################
#Example query
#"SELECT * FROM genes WHERE gene_id in ('XLOC_000005','XLOC_000015','XLOC_000055','XLOC_000595','XLOC_005998','ucscCodingXLOC_018816')"
.getLevels<-function(object){
	levelsQuery<-'SELECT s.sample_name FROM samples s ORDER BY s.sample_index ASC'
	levels<-dbGetQuery(object@DB,levelsQuery)$sample_name
	levels
}

setMethod("getLevels",signature(object="CuffData"),.getLevels)

.getRepLevels<-function(object){
	levelsQuery<-'SELECT r.rep_name FROM replicates r JOIN samples s ON r.sample_name=s.sample_name ORDER BY s.sample_index ASC'
	levels<-dbGetQuery(object@DB,levelsQuery)$rep_name
	levels
}

setMethod("getRepLevels",signature(object="CuffData"),.getRepLevels)

#Useful SQL commands

#SELECT g.gene_id, g.class_code, g.nearest_ref_id, g.gene_short_name, g.locus, g.length, g.coverage, g.status, gd.sample_name, gd.fpkm, gd.conf_hi, gd.conf_lo FROM genes g LEFT JOIN geneData gd ON g.gene_id = gd.gene_id WHERE (g.gene_id = 'XLOC_000001');

#SELECT g.gene_id, ged.* FROM genes g LEFT JOIN geneExpDiffData ged on g.gene_id = ged.gene_id WHERE ((sample_1 = 'H1_hESC' AND sample_2 = 'Fibroblasts') OR (sample_1 = 'Fibroblasts' AND sample_2 = 'H1_hESC')) AND ged.log2_fold_change>-20 AND ged.log2_fold_change<20 ;

#Pivot table SQL for scatterplots
#select g.*, sum(case when gd.sample_name = 'Fibroblasts' then fpkm end) as Fibroblasts, sum(case when gd.sample_name = 'H1_hESC' then fpkm end) as H1_hESC from genes g LEFT JOIN geneData gd on g.gene_id = gd.gene_id group by g.gene_id;


##################
#Plotting
##################

.density<-function(object, logMode = TRUE, pseudocount=0.0, labels, features=FALSE, replicates=FALSE,...){
	if(is(object,'CuffData')) {
		if(replicates){
			dat<-repFpkm(object,features=features)
			colnames(dat)[colnames(dat)=="rep_name"]<-"condition"
		}else{
			dat<-fpkm(object,features=features)
			colnames(dat)[colnames(dat)=="sample_name"]<-"condition"
		}
	} else {
		stop('Un-supported class of object.')
	}
	if(logMode) dat$fpkm<-dat$fpkm+pseudocount
	p<-ggplot(dat)
		if(logMode) {
			p<-p+geom_density(aes(x= log10(fpkm),group=condition,color=condition,fill=condition),alpha=I(1/3))
		}else{
			p<-p+geom_density(aes(x=fpkm,group=condition,color=condition,fill=condition),alpha=I(1/3))
		}
	
	p<-p + labs(title=object@tables$mainTable)
	
	#Default cummeRbund colorscheme
	p<-p + scale_fill_hue(l=50,h.start=200) + scale_color_hue(l=50,h.start=200)
	
	#TODO: Add label callout
	p
}

setMethod("csDensity",signature(object="CuffData"),.density)

.scatter<-function(object,x,y,logMode=TRUE,pseudocount=1.0,labels,smooth=FALSE,colorByStatus=FALSE, drawRug=TRUE, ...){
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
	
	#Right now, this does nothing, because 'significant' is not returned from fpkmMatrix object so I don't have this as a feature to draw
	if(colorByStatus){
		p<- p + geom_point(size=1.2,alpha=I(1/3))
	}else{
		p<- p + geom_point(size=1.2,alpha=I(1/3))
	}
	#Add symmetry line
	p<- p + geom_abline(intercept=0,slope=1,linetype=2) 
	
	#Add rug
	if(drawRug){
		p<- p + geom_rug(size=0.8,alpha=0.01)
	}
	
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
	p<- p + labs(title=object@tables$mainTable)
	p
}

setMethod("csScatter",signature(object="CuffData"), .scatter)

.scatterMat<-function(object,replicates=FALSE,logMode=TRUE,pseudocount=1.0,hexbin=FALSE,useCounts=FALSE,...){
	if(replicates){
		if(useCounts){
			dat<-repCountMatrix(object)
		}else{
			dat<-repFpkmMatrix(object)
		}
	}else{
		if(useCounts){
			dat<-countMatrix(object)
		}else{
			dat<-fpkmMatrix(object)
		}
	}
	
	if(hexbin){
		dat<-dat+pseudocount
	}
	
	if(useCounts){
		myLab = "Normalized Counts"
	}else{
		myLab = "FPKM"
	}
	if(logMode){
		myLab = paste("log10 ",myLab,sep="")
		p <- .plotmatrix(log10(dat),hexbin=hexbin,...)
	}else{
		p <- .plotmatrix(dat,hexbin=hexbin,...)
	}
	
	p<- p + geom_abline(intercept=0,slope=1,linetype=2)
	
	p <- p + theme_bw() + ylab(myLab) + xlab(myLab)
	
	#p<- p + aes(alpha=0.01)
	
	p
	
}

setMethod("csScatterMatrix",signature(object="CuffData"),.scatterMat)

.volcano<-function(object,x,y,alpha=0.05,showSignificant=TRUE,features=FALSE,xlimits=c(-20,20),...){
	samp<-samples(object)
	
	#check to make sure x and y are in samples
	if (!all(c(x,y) %in% samp)){
		stop("One or more values of 'x' or 'y' are not valid sample names!")
	}
	
	dat<-diffData(object=object,features=features)
	
	#subset dat for samples of interest
	mySamples<-c(x,y)
	dat<-dat[(dat$sample_1 %in% mySamples & dat$sample_2 %in% mySamples),]
	dat$significant <- 'no'
	dat$significant[dat$q_value<=alpha]<-'yes'
	s1<-unique(dat$sample_1)
	s2<-unique(dat$sample_2)
	
	p<-ggplot(dat)
	if(showSignificant){
		p<- p + geom_point(aes(x=log2_fold_change,y=-log10(p_value),color=significant),size=0.8)
	}else{
		p<- p + geom_point(aes(x=log2_fold_change,y=-log10(p_value)),size=1.2)
	}
	#Add title and return
	p<- p + labs(title=paste(object@tables$mainTable,": ",s2,"/",s1,sep=""))
	p<- p + scale_colour_manual(values = c("black","red"))
	
	#Set axis limits
	p<- p + scale_x_continuous(limits=xlimits)
	
	p <- p + xlab(bquote(paste(log[2],"(fold change)",sep=""))) + 
	    ylab(bquote(paste(-log[10],"(p value)",sep="")))
	p
}

setMethod("csVolcano",signature(object="CuffData"), .volcano)

.volcanoMatrix<-function(object,alpha=0.05,xlimits=c(-20,20),mapping=aes(),...){
	dat<-diffData(object)
	part1<-dat[,c('sample_1','sample_2','value_1','value_2','test_stat','p_value','q_value')]
	part2<-data.frame(sample_1=part1$sample_2,sample_2=part1$sample_1,value_1=part1$value_2,value_2=part1$value_1,test_stat=-part1$test_stat,p_value=part1$p_value,q_value=part1$q_value)
	dat<-rbind(part1,part2)
	
	myLevels<-union(dat$sample_1,dat$sample_2)
	dat$sample_1<-factor(dat$sample_1,levels=myLevels)
	dat$sample_2<-factor(dat$sample_2,levels=myLevels)
	dat$log2_fold_change<-log2(dat$value_2/dat$value_1)
	
	#Set significance value
	dat$significant<-"no"
	dat$significant[dat$q_value<=alpha]<-"yes"
	
	#May need to expand filler for time-series data (when there aren't always all pairwise comparisons on which to facet
	filler<-data.frame(sample_1=factor(myLevels,levels=myLevels),sample_2=factor(myLevels,levels=myLevels),label="")
	filler$label<-as.character(filler$label)
	mapping <- defaults(mapping, aes_string(x = "log2_fold_change", y = "-log10(p_value)", color="significant"))
	class(mapping) <- "uneval"
	
	p <-ggplot(dat) + geom_point(mapping,na.rm=TRUE,size=0.8) + scale_colour_manual(values = c("black","red")) + geom_text(aes(x=0,y=15,label=label),data=filler) + facet_grid(sample_1~sample_2)
	
	p<- p + geom_vline(aes(x=0),linetype=2)
	
	p <- p + theme_bw() + xlab(bquote(paste(log[2],"(fold change)",sep=""))) + ylab(bquote(paste(-log[10],"(p value)",sep="")))

	p
	
}

setMethod("csVolcanoMatrix",signature(object="CuffData"),.volcanoMatrix)

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
	g = g + geom_tile(color="black") + scale_x_discrete("", limits=labels[obj.hc$order]) + scale_y_discrete("", limits=labels[obj.hc$order])
	
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
	if(samples.not.genes){
		g <- g + geom_text(aes(label=format(value,digits=3)))
	}
	# return
	g
}

setMethod("csDistHeat", signature("CuffData"), .distheat)


.boxplot<-function(object,logMode=TRUE,pseudocount=0.0001,replicates=FALSE,...){
	if(replicates){
		dat<-repFpkm(object)
		colnames(dat)[colnames(dat)=="rep_name"]<-"condition"
	}else{
		dat<-fpkm(object)
		colnames(dat)[colnames(dat)=="sample_name"]<-"condition"
	}
	if(logMode) {
		dat$fpkm<-dat$fpkm+pseudocount
		p<-ggplot(dat)
		p<-p+geom_boxplot(aes(x=condition,y=log10(fpkm),fill=condition),size=0.3,alpha=I(1/3))
	} else {
		p <- ggplot(dat)
		p<-p+geom_boxplot(aes(x=condition,y=fpkm,fill=condition),alpha=I(1/3),size=0.3)
	}
	p<- p + theme(axis.text.x=element_text(angle=-90, hjust=0))
	
	#Default cummeRbund colorscheme
	p<-p + scale_fill_hue(l=50,h.start=200)
	
	p
}

setMethod("csBoxplot",signature(object="CuffData"),.boxplot)

.dendro<-function(object,logMode=T,pseudocount=1,replicates=FALSE,...){
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
	plot(res,main=paste("All",deparse(substitute(object)),sep=" "),...)
	res
}

setMethod("csDendro",signature(object="CuffData"),.dendro)

.MAplot<-function(object,x,y,logMode=T,pseudocount=1,smooth=FALSE,useCount=FALSE){
	if(useCount){
		dat<-.getCountMA(object,x,y,logMode=logMode,pseudocount=pseudocount)
	}else{
		dat<-.getMA(object,x,y,logMode=logMode,pseudocount=pseudocount)
	}
	p<-ggplot(dat)
	p<-p+geom_point(aes(x=A,y=log2(M)),size=0.8)
	
	#add smoother
	if(smooth){
		p <- p + stat_smooth(aes(x=A,y=log2(M)),fill="blue")
	}
	#add baseline
	p <- p + geom_hline(yintercept=0)
	
	p
}

setMethod("MAplot",signature(object="CuffData"),.MAplot)

.dispersionPlot<-function(object){
	dat<-count(object)
	p<-ggplot(dat)
	p<-p+geom_point(aes(x=count,y=dispersion,color=sample_name)) + facet_wrap('sample_name') + scale_x_log10() + scale_y_log10()
	p
}

setMethod("dispersionPlot",signature(object="CuffData"),.dispersionPlot)

.MDSplot<-function(object,replicates=FALSE,logMode=TRUE,pseudocount=1.0){
	if(replicates){
		dat<-repFpkmMatrix(object)
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

setMethod("MDSplot",signature(object="CuffData"),.MDSplot)

#Not sure if I want to include this or not..
.PCAplot<-function(object,x="PC1", y="PC2",replicates=FALSE,pseudocount=1.0,scale=TRUE,...){
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
	plot <- ggplot(dat, aes_string(x=x, y=y)) + geom_point(alpha=.4, size=0.8, aes(label=obsnames))
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
			geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 3, vjust=1, color="red")
	plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red") + theme_bw()
	plot
}

setMethod('PCAplot',signature(object="CuffData"),.PCAplot)

#This is most definitely a work in progress
.confidencePlot<-function(object,percentCutoff=20){
	res<-.statsMatrix(object)
	#res$CV[is.na(res$CV)]<-0.0
	res$CV<-as.numeric(res$CV)
	p<-ggplot(res)
	p<-p + 	geom_point(aes(x=log10(fpkm),y=log10(CV)),size=1.5,alpha=0.3)
	p<-p + 	geom_smooth(aes(x=log10(fpkm),y=log10(CV),color=sample_name)) +
	#p<-p + 	geom_point(aes(x=log10(fpkm),y=log10(count),color=-log10(CV))) +
	#p<-p + 	geom_point(aes(x=log10(fpkm),y=log2(fpkm/count),color=-log10(CV))) + geom_hline(aes(0),linetype=2) +
			#facet_wrap('sample_name') + 
			geom_abline(intercept=0,slope=1,linetype=2,size=0.3) + 
			scale_color_gradient(name="%CV",low="darkblue",high="white",limits=c(0,percentCutoff), na.value = "white") +
			labs(title=object@type)
	p
}

.CVdensity<-function(object){
	dat<-.statsMatrix(object)
	p<-ggplot(dat)
	p <- p + geom_density(aes(x=CV,fill=sample_name),alpha=0.3,position="dodge") + scale_x_log10()
	p
}

.fpkmSCVPlot<-function(object,FPKMLowerBound=1){
	dat<-repFpkm(object)
	colnames(dat)[1]<-"tracking_id"
	dat<-dat[,c('tracking_id','sample_name','fpkm')]
	dat<-dat[dat$fpkm>0,]
	
	#Option 3 (tapply on log10(replicateFPKM) values)
	dat.means<-tapply(dat$fpkm,dat[,c('tracking_id','sample_name')],function(x){mean(x,na.rm=T)})
	dat.sd<-tapply(dat$fpkm,dat[,c('tracking_id','sample_name')],function(x){sd(x,na.rm=T)})
	#write("Calculating replicate fpkm mean...",stderr())
	dat.means<-melt(dat.means)
	#write("Calculating replicate fpkm stdev...",stderr())
	dat.sd<-melt(dat.sd)
	colnames(dat.means)[colnames(dat.means)=="value"]<-'fpkm'
	colnames(dat.sd)[colnames(dat.sd)=="value"]<-'stdev'
	dat<-cbind(dat.means,dat.sd$stdev)
	colnames(dat)[colnames(dat)=="dat.sd$stdev"]<-'stdev'
	dat<-dat[!is.na(dat$stdev) & !is.na(dat$fpkm),]
	dat<-dat[dat$fpkm>0 & dat$stdev>0,]
	dat$sample_name<-factor(dat$sample_name,levels=samples(object))
	
	p <-ggplot(dat,aes(x=fpkm,y=(stdev/fpkm)^2),na.rm=T)
	#p <-ggplot(dat,aes(x=log10(fpkm+1),y=log10(stdev)),na.rm=T)
	p <- p + #geom_point(aes(color=sample_name),size=1,na.rm=T) +
		stat_smooth(aes(color=sample_name,fill=sample_name),na.rm=T,method='auto',fullrange=T) + 
		scale_x_log10() +
		scale_y_continuous(name=bquote(CV^2)) +
		xlab(bquote(paste(log[10],"FPKM",sep=" "))) +
		theme_bw() + xlim(c(log10(FPKMLowerBound),max(log10(dat$fpkm)))) + labs(title=object@type)
	p
	
}

setMethod("fpkmSCVPlot",signature(object="CuffData"),.fpkmSCVPlot)

.sensitivityPlot<-function(object){
	return
}

#TODO:Log2FC vs Test-statistic

#TODO:log2FPKM vs log2(stdev) (color by sample)

#TODO: Index of dispersion alpha (FPKM vs ?)
#SELECT gd.*, (gd.conf_hi-gd.fpkm)/2 as fpkm_stdev, gc.count, gc.variance as count_variance , gc.uncertainty, gc.dispersion, ((gd.conf_hi-gd.fpkm)/2)/gd.fpkm AS 'CV' FROM geneData gd LEFT JOIN geneCount gc ON gd.gene_id=gc.gene_id AND gd.sample_name=gc.sample_name;



#############
# Other Methods
#############
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

setMethod("csSpecificity",signature(object="CuffData"),.specificity)

#############
# GSEA helper methods
#############

.makeRnk<-function(object,x,y,filename,...){
	#Creates a log2-fold change .rnk file for all genes given an x/y comparison.
	#While this method will work for any cuffData level, it really only makes sense for 'genes' as this is what is required for GSEA...
	#Must provide 'x' and 'y' so as to provide a log2 fold change.
	samp<-samples(object)
	#check to make sure x and y are in samples
	if (!all(c(x,y) %in% samp)){
		stop("One or more values of 'x' or 'y' are not valid sample names!")
	}
	query<-paste("SELECT gd.gene_id,g.gene_short_name,(sum(CASE WHEN gd.sample_name='",x,"' THEN gd.fpkm+1 END))/(sum(CASE WHEN gd.sample_name='",y,"' THEN gd.fpkm+1 END)) as 'ratio' FROM genes g LEFT JOIN geneData gd ON g.gene_id=gd.gene_id GROUP BY g.gene_id ORDER BY ratio DESC",sep="")
	res<-dbGetQuery(object@DB,query)
	res$ratio<-log2(res$ratio)
	#Remove gene_id field
	res<-res[,-1]
	#Remove rows with "NA" for gene_short_name
	res<-res[!is.na(res$gene_short_name),]
	#Write to file
	write.table(res,file=filename,sep="\t",quote=F,...,row.names=F,col.names=F)
	#res
}

setMethod("makeRnk",signature(object="CuffData"),.makeRnk)


#############
#Utility functions
#############
.checkSamples<-function(dbConn,sampleIdList){
	dbSamples<-dbReadTable(dbConn,"samples")
	if (all(sampleIdList %in% dbSamples$sample_name)){
		return(TRUE)
	}else{
		return(FALSE)
	}
}

.checkReps<-function(dbConn,repIdList){
	dbReps<-dbReadTable(dbConn,"replicates")
	if (all(repIdList %in% dbReps$rep_name)){
		return(TRUE)
	}else{
		return(FALSE)
	}
}

###################
# Coersion methods
###################
#As ExpressionSet

###################
#Utility functions
###################

#.volcanoMatrix <- function(data){
#	densities <- do.call('rbind',lapply(1:))
#	p <-ggplot(data) + facet_grid(sample1~sample2,scales="free") + geom_point(aes(x=log2_fold_change,y=-log10(p_value))) + stat_density(aes(x = log2_fold_change, 
#					y = ..scaled.. * diff(range(log2_fold_change)) + min(log2_fold_change)), data = densities, 
#			position = "identity", colour = "grey20", geom = "line")
#	
#}

