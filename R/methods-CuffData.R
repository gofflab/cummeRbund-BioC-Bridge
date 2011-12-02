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
					tables=list(mainTable = "",dataTable = "",expDiffTable = "",featureTable = ""),
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

.fpkm<-function(object,features=FALSE){
	if(!features){
		FPKMQuery<-paste("SELECT * FROM",object@tables$dataTable)
	}else{
		FPKMQuery<-paste("SELECT xf.*,x.sample_name,x.fpkm,x.conf_hi, x.conf_lo FROM ",object@tables$dataTable," x LEFT JOIN ",object@tables$featureTable," xf ON x.",object@idField,"=xf.",object@idField,sep="")
	}
	res<-dbGetQuery(object@DB,FPKMQuery)
	res$sample_name<-factor(res$sample_name,levels=getLevels(object))
	res
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


#This needs a lot of work...
#TODO: Change this to remove lnFcCutoff but make sure that functions that rely on diffData have their own FC cutoff so that plotting doesn't suffer
.diffData<-function(object,x,y,features=FALSE){
	if(missing(x) && missing(y)){
		if(!features){
			diffQuery<-paste("SELECT * FROM ",object@tables$expDiffTable,sep="")
		}else{
			diffQuery<-paste("SELECT * FROM ",object@tables$expDiffTable," x LEFT JOIN ",object@tables$featureTable," xf ON x.",object@idField,"=xf.",object@idField,sep="")
		}
	}else if (missing(x) || missing(y)){
		stop("You must supply both x and y or neither.")
	}else{
		if(!features){
			diffQuery<-paste("SELECT x.",object@idField,", xed.* FROM ",object@tables$mainTable," x LEFT JOIN ",object@tables$expDiffTable," xed on x.",object@idField," = xed.",object@idField," WHERE ((sample_1 = '",x,"' AND sample_2 = '",y,"') OR (sample_1 = '",y,"' AND sample_2 = '",x,"'))",sep="")
		}else{
			diffQuery<-paste("SELECT x.",object@idField,", xed.*, xf.* FROM ",object@tables$mainTable," x LEFT JOIN ",object@tables$expDiffTable," xed on x.",object@idField," = xed.",object@idField," LEFT JOIN ",object@tables$featureTable," xf ON x.",object@idField,"=xf.",object@idField," WHERE ((sample_1 = '",x,"' AND sample_2 = '",y,"') OR (sample_1 = '",y,"' AND sample_2 = '",x,"'))",sep="")
		}
	}
	dat<-dbGetQuery(object@DB,diffQuery)
	#diffQuery
	dat
}

setMethod("diffData",signature(object="CuffData"),.diffData)

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

#Useful SQL commands

#SELECT g.gene_id, g.class_code, g.nearest_ref_id, g.gene_short_name, g.locus, g.length, g.coverage, g.status, gd.sample_name, gd.fpkm, gd.conf_hi, gd.conf_lo FROM genes g LEFT JOIN geneData gd ON g.gene_id = gd.gene_id WHERE (g.gene_id = 'XLOC_000001');

#SELECT g.gene_id, ged.* FROM genes g LEFT JOIN geneExpDiffData ged on g.gene_id = ged.gene_id WHERE ((sample_1 = 'H1_hESC' AND sample_2 = 'Fibroblasts') OR (sample_1 = 'Fibroblasts' AND sample_2 = 'H1_hESC')) AND ged.ln_fold_change>-20 AND ged.ln_fold_change<20 ;

#Pivot table SQL for scatterplots
#select g.*, sum(case when gd.sample_name = 'Fibroblasts' then fpkm end) as Fibroblasts, sum(case when gd.sample_name = 'H1_hESC' then fpkm end) as H1_hESC from genes g LEFT JOIN geneData gd on g.gene_id = gd.gene_id group by g.gene_id;


##################
#Plotting
##################

.density<-function(object, logMode = TRUE, pseudocount=1.0, labels, features=FALSE, ...){
	if(is(object,'CuffData')) {
		dat<-fpkm(object,features=features)
	} else {
		stop('Un-supported class of object.')
	}
	if(logMode) dat$fpkm<-dat$fpkm+pseudocount
	p<-ggplot(dat)
		if(logMode) {
			p<-p+geom_density(aes(x= log10(fpkm),group=sample_name,color=sample_name,fill=sample_name),alpha=I(1/3))
		}else{
			p<-p+geom_density(aes(x=fpkm,group=sample_name,color=sample_name,fill=sample_name),alpha=I(1/3))
		}
	
	p<-p + opts(title=object@tables$mainTable)
	
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
		p<- p + geom_rug(size=0.5,alpha=0.01)
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
	p<- p + opts(title=object@tables$mainTable)
	p
}

setMethod("csScatter",signature(object="CuffData"), .scatter)

.volcano<-function(object,x,y,features=FALSE,xlimits=c(-20,20),...){
	samp<-samples(object)
	
	#check to make sure x and y are in samples
	if (!all(c(x,y) %in% samp)){
		stop("One or more values of 'x' or 'y' are not valid sample names!")
	}
	
	dat<-diffData(object=object,features=features)
	
	#subset dat for samples of interest
	mySamples<-c(x,y)
	dat<-dat[(dat$sample_1 %in% mySamples & dat$sample_2 %in% mySamples),]
	s1<-unique(dat$sample_1)
	s2<-unique(dat$sample_2)
	
	p<-ggplot(dat)
	p<- p + geom_point(aes(x=ln_fold_change,y=-log10(p_value),color=significant),size=0.8)
	
	#Add title and return
	p<- p + opts(title=paste(object@tables$mainTable,": ",s2,"/",s1,sep=""))
	p<- p + scale_colour_manual(values=c("red", "steelblue"))
	
	#Set axis limits
	p<- p + scale_x_continuous(limits=xlimits)
	
	p <- p + xlab(bquote(paste(log[2],"(fold change)",sep=""))) + 
	    ylab(bquote(paste(-log[10],"(p value)",sep="")))
	p
}

setMethod("csVolcano",signature(object="CuffData"), .volcano)

.boxplot<-function(object,logMode=TRUE,pseudocount=0.0001,...){
	dat<-fpkm(object)
	if(logMode) {
		dat$fpkm<-dat$fpkm+pseudocount
		p<-ggplot(dat)
		p<-p+geom_boxplot(aes(x=sample_name,y=log10(fpkm),fill=sample_name),size=0.3,alpha=I(1/3))
	} else {
		p <- ggplot(dat)
		p<-p+geom_boxplot(aes(x=sample_name,y=fpkm,fill=sample_name),alpha=I(1/3),size=0.3)
	}
	p<- p + opts(axis.text.x=theme_text(angle=-90, hjust=0))
	
	#Default cummeRbund colorscheme
	p<-p + scale_fill_hue(l=50,h.start=200)
	
	p
}

setMethod("csBoxplot",signature(object="CuffData"),.boxplot)

#############
# Other Methods
#############

