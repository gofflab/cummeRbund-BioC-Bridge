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
					tables,
					filters,
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

##################
#Database connectivity
##################




##################
#Data Retrieval
##################
.names<-function(object){
	
}

setMethod("names","CuffData",.names)

.features<-function(object){
	dbReadTable(object@DB, object@tables$mainTable)
}

setMethod("features","CuffData",.features)

.fpkm<-function(object){
	FPKMQuery<-paste("SELECT * FROM",object@tables$dataTable)
	dbGetQuery(object@DB,FPKMQuery)
}

setMethod("fpkm","CuffData",.fpkm)

.diff<-function(){
	
}

setMethod("diff","CuffData",.diff)

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

.scatter<-function(){
	
}

.volcano<-function(){
	
}