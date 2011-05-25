#methods-CuffData.R
#
#
#
#
####################

##################
#Initialize
##################
setMethod("initialize","CuffData",
			function(.Object,
					DB,
					tables=list(mainTable = "",dataTable = "",expDiffTable = "",featureTable = "", otherTable = ""),
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
.features<-function(object){
	featureQuery<-paste("SELECT * FROM ",object@tables$mainTable," x LEFT JOIN ",object@tables$featureTable," xf ON x.",object@idField,"=xf.",object@idField,sep="")
	dbGetQuery(object@DB, featureQuery)
}

setMethod("features","CuffData",.features)

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

.fpkm<-function(object){
	FPKMQuery<-paste("SELECT * FROM",object@tables$dataTable)
	dbGetQuery(object@DB,FPKMQuery)
}

setMethod("fpkm","CuffData",.fpkm)

.fpkmMatrix<-function(object){
	samp<-samples(object)
	FPKMMatQuery<-paste("select x.",object@idField,", ",sep="")
	for (i in samp){
		FPKMMatQuery<-paste(FPKMMatQuery,"sum(case when xd.sample_name ='",i,"' then fpkm end) as ",i,",",sep="")
	}
	FPKMMatQuery<-substr(FPKMMatQuery, 1, nchar(FPKMMatQuery)-1)
	FPKMMatQuery<-paste(FPKMMatQuery," from ",object@tables$mainTable," x LEFT JOIN ",object@tables$dataTable," xd on x.",object@idField," = xd.",object@idField," group by x.",object@idField,sep="")
	res<-dbGetQuery(object@DB,FPKMMatQuery)
	data.frame(res[,-1],row.names=res[,1])
}

setMethod("fpkmMatrix","CuffData",.fpkmMatrix)

.diffData<-function(object,x,y,lnFcCutoff=20){
	diffQuery<-paste("SELECT x.",object@idField,", xed.* FROM ",object@tables$mainTable," x LEFT JOIN ",object@tables$expDiffTable," xed on x.",object@idField," = xed.",object@idField," WHERE ((sample_1 = '",x,"' AND sample_2 = '",y,"') OR (sample_1 = '",y,"' AND sample_2 = '",x,"')) AND xed.ln_fold_change>",-lnFcCutoff," AND xed.ln_fold_change<",lnFcCutoff,sep="")
	dat<-dbGetQuery(object@DB,diffQuery)
	#diffQuery
	dat
}

setMethod("diffData",signature(object="CuffData"),.diffData)

##################
#Setters
##################


##################
#Subsetting
##################
#Example query
#"SELECT * FROM genes WHERE gene_id in ('XLOC_000005','XLOC_000015','XLOC_000055','XLOC_000595','XLOC_005998','ucscCodingXLOC_018816')"




#Useful SQL commands

#SELECT g.gene_id, g.class_code, g.nearest_ref_id, g.gene_short_name, g.locus, g.length, g.coverage, g.status, gd.sample_name, gd.fpkm, gd.conf_hi, gd.conf_lo FROM genes g LEFT JOIN geneData gd ON g.gene_id = gd.gene_id WHERE (g.gene_id = 'XLOC_000001');

#SELECT g.gene_id, ged.* FROM genes g LEFT JOIN geneExpDiffData ged on g.gene_id = ged.gene_id WHERE ((sample_1 = 'H1_hESC' AND sample_2 = 'Fibroblasts') OR (sample_1 = 'Fibroblasts' AND sample_2 = 'H1_hESC')) AND ged.ln_fold_change>-20 AND ged.ln_fold_change<20 ;

#Pivot table SQL for scatterplots
#select g.*, sum(case when gd.sample_name = 'Fibroblasts' then fpkm end) as Fibroblasts, sum(case when gd.sample_name = 'H1_hESC' then fpkm end) as H1_hESC from genes g LEFT JOIN geneData gd on g.gene_id = gd.gene_id group by g.gene_id;


##################
#Plotting
##################

.density<-function(object, logMode = TRUE, pseudocount=0.0001, labels, ...){
	if(is(object,'CuffData')) {
		dat<-fpkm(object)
	} else {
		stop('Un-supported class of object.')
	}
	if(logMode) dat$fpkm<-dat$fpkm+pseudocount
	p<-ggplot(dat)
		if(logMode) {
			p<-p+geom_density(aes(x=log2(fpkm),group=sample_name,color=sample_name,fill=sample_name),alpha=I(1/3))
		}else{
			p<-p+geom_density(aes(x=fpkm,group=sample_name,color=sample_name,fill=sample_name),alpha=I(1/3))
		}
	
	p<-p + opts(title=object@tables$mainTable)	
	#TODO: Add label callout
	p
}

setMethod("csDensity",signature(object="CuffData"),.density)

.boxplot<-function(){
	
}

.expressionPlot<-function(){
	
}

.barplot<-function(){
	
}

.heatmap<-function(){
	
}

.ggheat<-function(){
	
}

.scatter<-function(object,x,y,logMode=TRUE,pseudocount=0.0001,labels, smooth=FALSE,...){
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
		p <- p + scale_y_log2() + scale_x_log2()
	}
	
	#Add title & Return value
	p<- p + opts(title=object@tables$mainTable)
	p
}

setMethod("csScatter",signature(object="CuffData"), .scatter)

.volcano<-function(object,x,y){
	dat<-diffData(object=object,x=x,y=y)
	s1<-unique(dat$sample_1)
	s2<-unique(dat$sample_2)
	
	p<-ggplot(dat)
	p<- p + geom_point(aes(x=ln_fold_change,y=-log10(p_value),color=significant),size=1,alpha=I(1/3))
	
	#Add title and return
	p<- p + opts(title=paste(object@tables$mainTable,": ",s2,"/",s1,sep=""))
	p
}

setMethod("csVolcano",signature(object="CuffData"), .volcano)


