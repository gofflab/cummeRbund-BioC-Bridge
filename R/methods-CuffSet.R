##################
#methods-CuffSet.R
#
#Introduces the CuffSet Class for analysis, manipulation, and plotting of Cufflinks data
#
#Author: Loyal A. Goff
#
##################

#Initialize
setMethod("initialize","CuffSet",
		function(.Object,
				DB,
				conditions=data.frame(),
				genes,
				isoforms,
				TSS,
				CDS,
				...){
			.Object<-callNextMethod(.Object,
					DB = DB,
					conditions = conditions,
					genes = genes,
					isoforms = isoforms,
					TSS = TSS,
					CDS = CDS,
					...)				
		}
)

##################
#Class Methods
##################
setMethod("show","CuffSet",
		function(object){
			cat(class(object), "instance with:\n\t",
					dim(object@genes)[2],"samples\n\t",
					dim(object@genes)[1],"genes\n\t",
					dim(object@isoforms)[1],"isoforms\n\t",
					dim(object@TSS)[1],"TSS\n\t",
					dim(object@CDS)[1],"CDS\n"
					)
		}
)

#This does not subset appropriately yet
#TODO: Fix for multiple values of i
#
#Solution is to test i to determine if it is of type 'numeric' (index), list (multi-index), or 'character' (gene_ids)
#
#TODO: 	- Add 'j' to select on sampleNames as well
#		- Add ability to search on gene_short_name(s) or featureIDs
		
setMethod("[",signature(x="CuffSet"),function(x, i, ...){
			featureIDs<-featureNames(x@genes)[i]
			if(length(featureIDs)==1){
				res<-getGene(x,featureID)
			}else{
				res<-getGenes(x,featureIDs)
			}
			res
		}
)


setValidity("CuffSet",
		function(object){
		TRUE	
		}
)

############
#Accessors
############
.samples<-function(object){
	sampleQuery<-"SELECT * FROM samples s LEFT JOIN phenoData p on s.sample_name = p.sample_name"
	dbGetQuery(object@DB,sampleQuery)
}

setMethod("samples",signature(object="CuffSet"),.samples)

#make CuffGene objects from a gene_ids
.getGene<-function(object,geneId){
	
	#dbQueries
	geneAnnotationQuery<-paste("SELECT * from genes WHERE gene_id ='",geneId,"'",sep="")
	geneFPKMQuery<-paste("SELECT y.* from genes x JOIN geneData y ON x.gene_id=y.gene_id WHERE x.gene_id ='",geneId,"'",sep="")
	geneDiffQuery<-paste("SELECT y.* from genes x JOIN geneExpDiffData y ON x.gene_id=y.gene_id WHERE x.gene_id ='",geneId,"'",sep="")
	
	isoformAnnotationQuery<-paste("SELECT * from isoforms WHERE gene_id ='",geneId,"'",sep="")
	isoformFPKMQuery<-paste("SELECT y.* from isoforms x JOIN isoformData y ON x.isoform_id = y.isoform_id WHERE x.gene_id ='",geneId,"'",sep="")
	isoformDiffQuery<-paste("SELECT y.* from isoforms x JOIN isoformExpDiffData y ON x.isoform_id = y.isoform_id WHERE x.gene_id ='",geneId,"'",sep="")
	
	TSSAnnotationQuery<-paste("SELECT * from TSS WHERE gene_id ='",geneId,"'",sep="")
	TSSFPKMQuery<-paste("SELECT y.* from TSS x JOIN TSSData y ON x.TSS_group_id=y.TSS_group_id WHERE x.gene_id ='",geneId,"'",sep="")
	TSSDiffQuery<-paste("SELECT y.* from TSS x JOIN TSSExpDiffData y ON x.TSS_group_id=y.TSS_group_id WHERE x.gene_id ='",geneId,"'",sep="")
	
	CDSAnnotationQuery<-paste("SELECT * from CDS WHERE gene_id ='",geneId,"'",sep="")
	CDSFPKMQuery<-paste("SELECT y.* from CDS x JOIN CDSData y ON x.CDS_id = y.CDS_id WHERE x.gene_id ='",geneId,"'",sep="")
	CDSDiffQuery<-paste("SELECT y.* from CDS x JOIN CDSExpDiffData y ON x.CDS_id = y.CDS_id WHERE x.gene_id ='",geneId,"'",sep="")
	
	res<-new("CuffGene",
			id=geneId,
			annotation=dbGetQuery(object@DB,geneAnnotationQuery),
			fpkm=dbGetQuery(object@DB,geneFPKMQuery),
			diff=dbGetQuery(object@DB,geneDiffQuery),
			isoforms=new("CuffFeature",
					annotation=dbGetQuery(object@DB,isoformAnnotationQuery),
					fpkm=dbGetQuery(object@DB,isoformFPKMQuery),
					diff=dbGetQuery(object@DB,isoformDiffQuery)
					),
			TSS=new("CuffFeature",
					annotation=dbGetQuery(object@DB,TSSAnnotationQuery),
					fpkm=dbGetQuery(object@DB,TSSFPKMQuery),
					diff=dbGetQuery(object@DB,TSSDiffQuery)
			),
			CDS=new("CuffFeature",
					annotation=dbGetQuery(object@DB,CDSAnnotationQuery),
					fpkm=dbGetQuery(object@DB,CDSFPKMQuery),
					diff=dbGetQuery(object@DB,CDSDiffQuery)
			)

			
		)
		
		res
}

setMethod("getGene",signature(object="CuffSet"),.getGene)
	
.getGenes<-function(object,geneIdList){
	#Make WHERE Clause search string
	whereString<-'WHERE x.gene_id %in% ('
	for (i in geneIdList){
		whereString<-paste(whereString,"'",i,"',",sep="")
	}
	whereString<-substr(whereString,1,nchar(whereString)-1)
	whereString<-paste(whereString,")",sep="")
	
	#dbQueries
	geneAnnotationQuery<-paste("SELECT * from genes x ", whereString,sep="")
	geneFPKMQuery<-paste("SELECT y.* from genes x JOIN geneData y ON x.gene_id=y.gene_id ", whereString,sep="")
	geneDiffQuery<-paste("SELECT y.* from genes x JOIN geneExpDiffData y ON x.gene_id=y.gene_id ", whereString,sep="")
	
	isoformAnnotationQuery<-paste("SELECT * from isoforms ", whereString,sep="")
	isoformFPKMQuery<-paste("SELECT y.* from isoforms x JOIN isoformData y ON x.isoform_id = y.isoform_id ", whereString,sep="")
	isoformDiffQuery<-paste("SELECT y.* from isoforms x JOIN isoformExpDiffData y ON x.isoform_id = y.isoform_id ", whereString,sep="")
	
	TSSAnnotationQuery<-paste("SELECT * from TSS ", whereString,sep="")
	TSSFPKMQuery<-paste("SELECT y.* from TSS x JOIN TSSData y ON x.TSS_group_id=y.TSS_group_id ", whereString,sep="")
	TSSDiffQuery<-paste("SELECT y.* from TSS x JOIN TSSExpDiffData y ON x.TSS_group_id=y.TSS_group_id ", whereString,sep="")
	
	CDSAnnotationQuery<-paste("SELECT * from CDS ", whereString,sep="")
	CDSFPKMQuery<-paste("SELECT y.* from CDS x JOIN CDSData y ON x.CDS_id = y.CDS_id ", whereString,sep="")
	CDSDiffQuery<-paste("SELECT y.* from CDS x JOIN CDSExpDiffData y ON x.CDS_id = y.CDS_id ", whereString,sep="")

	res<-new("CuffGeneSet",
			ids=geneIdList,
			annotation=dbGetQuery(object@DB,geneAnnotationQuery),
			fpkm=dbGetQuery(object@DB,geneFPKMQuery),
			diff=dbGetQuery(object@DB,geneDiffQuery),
			isoforms=new("CuffFeatureSet",
					annotation=dbGetQuery(object@DB,isoformAnnotationQuery),
					fpkm=dbGetQuery(object@DB,isoformFPKMQuery),
					diff=dbGetQuery(object@DB,isoformDiffQuery)
					),
			TSS=new("CuffFeatureSet",
					annotation=dbGetQuery(object@DB,TSSAnnotationQuery),
					fpkm=dbGetQuery(object@DB,TSSFPKMQuery),
					diff=dbGetQuery(object@DB,TSSDiffQuery)
					),
			CDS=new("CuffFeatureSet",
					annotation=dbGetQuery(object@DB,CDSAnnotationQuery),
					fpkm=dbGetQuery(object@DB,CDSFPKMQuery),
					diff=dbGetQuery(object@DB,CDSDiffQuery)
					)
			)
	res
}

setMethod("getGenes",signature(object="CuffSet"),.getGenes)

############
#SQL access
############


#####################
#Add FeatureData    #
#####################
.addFeatures<-function(object,features,level="genes"){
	if(!is.data.frame(features)){
		stop("features must be a data.frame")
	}
	colnames(features)[1]<-slot(object,level)@idField
	colnames(features)<-make.db.names(object@DB,colnames(features),unique=T)
	dbWriteTable(object@DB,slot(object,level)@tables$featureTable,features,row.names=F,overwrite=T)
	indexQuery<-paste("CREATE INDEX ",slot(object,level)@idField," ON ", slot(object,level)@tables$featureTable," (",slot(object,level)@idField,")",sep="")
	res<-dbGetQuery(object@DB,indexQuery)
}

setMethod("addFeatures",signature(object="CuffSet"),.addFeatures)
