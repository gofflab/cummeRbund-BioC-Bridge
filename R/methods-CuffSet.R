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
				promoters,
				splicing,
				relCDS,
				...){
			.Object<-callNextMethod(.Object,
					DB = DB,
					conditions = conditions,
					genes = genes,
					isoforms = isoforms,
					TSS = TSS,
					CDS = CDS,
					promoters = promoters,
					splicing = splicing,
					relCDS = relCDS,
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
					dim(object@CDS)[1],"CDS\n\t",
					dim(object@promoters)[1],"promoters\n\t",
					dim(object@splicing)[1],"splicing\n\t",
					dim(object@relCDS)[1],"relCDS\n"
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
	
	#get sample_name levels
	myLevels<-getLevels(object)
	
	whereString = paste("WHERE x.gene_id ='",geneId,"' OR x.gene_short_name = '",geneId,"'",sep="")
	
	#dbQueries
	geneAnnotationQuery<-paste("SELECT * from genes x ",whereString,sep="")
	geneFPKMQuery<-paste("SELECT y.* from genes x JOIN geneData y ON x.gene_id=y.gene_id ",whereString,sep="")
	geneDiffQuery<-paste("SELECT y.* from genes x JOIN geneExpDiffData y ON x.gene_id=y.gene_id ",whereString,sep="")
	
	isoformAnnotationQuery<-paste("SELECT * from isoforms i JOIN genes x ON i.gene_id = x.gene_id ",whereString,sep="")
	isoformFPKMQuery<-paste("SELECT y.* from isoforms i JOIN isoformData y ON i.isoform_id = y.isoform_id JOIN genes x ON i.gene_id = x.gene_id ",whereString,sep="")
	isoformDiffQuery<-paste("SELECT y.* from isoforms i JOIN isoformExpDiffData y ON i.isoform_id = y.isoform_id JOIN genes x ON i.gene_id = x.gene_id ",whereString,sep="")
	
	TSSAnnotationQuery<-paste("SELECT * from TSS t JOIN genes x ON t.gene_id = x.gene_id ",whereString,sep="")
	TSSFPKMQuery<-paste("SELECT y.* from TSS t JOIN TSSData y ON t.TSS_group_id=y.TSS_group_id JOIN genes x ON t.gene_id = x.gene_id ",whereString,sep="")
	TSSDiffQuery<-paste("SELECT y.* from TSS t JOIN TSSExpDiffData y ON t.TSS_group_id=y.TSS_group_id JOIN genes x ON t.gene_id = x.gene_id ",whereString,sep="")
	
	CDSAnnotationQuery<-paste("SELECT * from CDS c JOIN genes x ON c.gene_id = x.gene_id ",whereString,sep="")
	CDSFPKMQuery<-paste("SELECT y.* from CDS c JOIN CDSData y ON c.CDS_id = y.CDS_id JOIN genes x ON c.gene_id = x.gene_id ",whereString,sep="")
	CDSDiffQuery<-paste("SELECT y.* from CDS c JOIN CDSExpDiffData y ON c.CDS_id = y.CDS_id JOIN genes x ON c.gene_id = x.gene_id ",whereString,sep="")
	
	begin<-dbSendQuery(object@DB,"BEGIN;")
	
	#fetch records
	#genes
	genes.fpkm<-dbGetQuery(object@DB,geneFPKMQuery)
	genes.fpkm$sample_name<-factor(genes.fpkm$sample_name,levels=myLevels)
	genes.diff<-dbGetQuery(object@DB,geneDiffQuery)
	genes.diff$sample_1<-factor(genes.diff$sample_1,levels=myLevels)
	genes.diff$sample_2<-factor(genes.diff$sample_2,levels=myLevels)
	
	#isoforms
	isoform.fpkm<-dbGetQuery(object@DB,isoformFPKMQuery)
	isoform.fpkm$sample_name<-factor(isoform.fpkm$sample_name,levels=myLevels)
	isoform.diff<-dbGetQuery(object@DB,isoformDiffQuery)
	isoform.diff$sample_1<-factor(isoform.diff$sample_1,levels=myLevels)
	isoform.diff$sample_2<-factor(isoform.diff$sample_2,levels=myLevels)
	
	#CDS
	CDS.fpkm<-dbGetQuery(object@DB,CDSFPKMQuery)
	CDS.fpkm$sample_name<-factor(CDS.fpkm$sample_name,levels=myLevels)
	CDS.diff<-dbGetQuery(object@DB,CDSDiffQuery)
	CDS.diff$sample_1<-factor(CDS.diff$sample_1,levels=myLevels)
	CDS.diff$sample_2<-factor(CDS.diff$sample_2,levels=myLevels)
	
	#TSS
	TSS.fpkm<-dbGetQuery(object@DB,TSSFPKMQuery)
	TSS.fpkm$sample_name<-factor(TSS.fpkm$sample_name,levels=myLevels)
	TSS.diff<-dbGetQuery(object@DB,TSSDiffQuery)
	TSS.diff$sample_1<-factor(TSS.diff$sample_1,levels=myLevels)
	TSS.diff$sample_2<-factor(TSS.diff$sample_2,levels=myLevels)
	
	res<-new("CuffGene",
			id=geneId,
			annotation=dbGetQuery(object@DB,geneAnnotationQuery),
			fpkm=genes.fpkm,
			diff=genes.diff,
			isoforms=new("CuffFeature",
					annotation=dbGetQuery(object@DB,isoformAnnotationQuery),
					fpkm=isoform.fpkm,
					diff=isoform.diff
					),
			TSS=new("CuffFeature",
					annotation=dbGetQuery(object@DB,TSSAnnotationQuery),
					fpkm=TSS.fpkm,
					diff=TSS.diff
			),
			CDS=new("CuffFeature",
					annotation=dbGetQuery(object@DB,CDSAnnotationQuery),
					fpkm=CDS.fpkm,
					diff=CDS.diff
			)

			
		)
	end<-dbSendQuery(object@DB,"END;")
		
		res
}

setMethod("getGene",signature(object="CuffSet"),.getGene)
	
.getGenes<-function(object,geneIdList){
	#Make WHERE Clause search string
#	whereString<-'WHERE x.gene_id IN ('
#	for (i in geneIdList){
#		whereString<-paste(whereString,"'",i,"',",sep="")
#	}
#	whereString<-substr(whereString,1,nchar(whereString)-1)
#	whereString<-paste(whereString,")",sep="")
	
	idString<-'('
	for (i in geneIdList){
		idString<-paste(idString,"'",i,"',",sep="")
	}
	idString<-substr(idString,1,nchar(idString)-1)
	idString<-paste(idString,")",sep="")
	
	whereStringGene<-paste('WHERE x.gene_id IN ',idString,' OR x.gene_short_name IN ',idString,sep="")
	whereString<-paste('WHERE x.gene_id IN ',idString,' OR g.gene_short_name IN ',idString,sep="")
	
	#get sample_name levels
	myLevels<-getLevels(object)
	
	#dbQueries
	idQuery<-paste("SELECT DISTINCT gene_id from genes x ",whereStringGene,sep="")
	
	geneAnnotationQuery<-paste("SELECT * from genes x ", whereStringGene,sep="")
	geneFPKMQuery<-paste("SELECT y.* from genes x JOIN geneData y ON x.gene_id=y.gene_id ", whereStringGene,sep="")
	geneDiffQuery<-paste("SELECT y.* from genes x JOIN geneExpDiffData y ON x.gene_id=y.gene_id ", whereStringGene,sep="")
	
	isoformAnnotationQuery<-paste("SELECT * from isoforms x JOIN genes g on x.gene_id=g.gene_id ", whereString,sep="")
	isoformFPKMQuery<-paste("SELECT y.* from isoforms x JOIN isoformData y ON x.isoform_id = y.isoform_id JOIN genes g on x.gene_id=g.gene_id ", whereString,sep="")
	isoformDiffQuery<-paste("SELECT y.* from isoforms x JOIN isoformExpDiffData y ON x.isoform_id = y.isoform_id JOIN genes g on x.gene_id=g.gene_id ", whereString,sep="")
	
	TSSAnnotationQuery<-paste("SELECT * from TSS x JOIN genes g on x.gene_id=g.gene_id ", whereString,sep="")
	TSSFPKMQuery<-paste("SELECT y.* from TSS x JOIN TSSData y ON x.TSS_group_id=y.TSS_group_id JOIN genes g on x.gene_id=g.gene_id ", whereString,sep="")
	TSSDiffQuery<-paste("SELECT y.* from TSS x JOIN TSSExpDiffData y ON x.TSS_group_id=y.TSS_group_id JOIN genes g on x.gene_id=g.gene_id ", whereString,sep="")
	
	CDSAnnotationQuery<-paste("SELECT * from CDS x JOIN genes g on x.gene_id=g.gene_id ", whereString,sep="")
	CDSFPKMQuery<-paste("SELECT y.* from CDS x JOIN CDSData y ON x.CDS_id = y.CDS_id JOIN genes g on x.gene_id=g.gene_id ", whereString,sep="")
	CDSDiffQuery<-paste("SELECT y.* from CDS x JOIN CDSExpDiffData y ON x.CDS_id = y.CDS_id JOIN genes g on x.gene_id=g.gene_id ", whereString,sep="")
	
	begin<-dbSendQuery(object@DB,"BEGIN;")
	
	#fetch records
	#genes
	genes.fpkm<-dbGetQuery(object@DB,geneFPKMQuery)
	genes.fpkm$sample_name<-factor(genes.fpkm$sample_name,levels=myLevels)
	genes.diff<-dbGetQuery(object@DB,geneDiffQuery)
	genes.diff$sample_1<-factor(genes.diff$sample_1,levels=myLevels)
	genes.diff$sample_2<-factor(genes.diff$sample_2,levels=myLevels)
	
	#isoforms
	isoform.fpkm<-dbGetQuery(object@DB,isoformFPKMQuery)
	isoform.fpkm$sample_name<-factor(isoform.fpkm$sample_name,levels=myLevels)
	isoform.diff<-dbGetQuery(object@DB,isoformDiffQuery)
	isoform.diff$sample_1<-factor(isoform.diff$sample_1,levels=myLevels)
	isoform.diff$sample_2<-factor(isoform.diff$sample_2,levels=myLevels)
	
	#CDS
	CDS.fpkm<-dbGetQuery(object@DB,CDSFPKMQuery)
	CDS.fpkm$sample_name<-factor(CDS.fpkm$sample_name,levels=myLevels)
	CDS.diff<-dbGetQuery(object@DB,CDSDiffQuery)
	CDS.diff$sample_1<-factor(CDS.diff$sample_1,levels=myLevels)
	CDS.diff$sample_2<-factor(CDS.diff$sample_2,levels=myLevels)
	
	#TSS
	TSS.fpkm<-dbGetQuery(object@DB,TSSFPKMQuery)
	TSS.fpkm$sample_name<-factor(TSS.fpkm$sample_name,levels=myLevels)
	TSS.diff<-dbGetQuery(object@DB,TSSDiffQuery)
	TSS.diff$sample_1<-factor(TSS.diff$sample_1,levels=myLevels)
	TSS.diff$sample_2<-factor(TSS.diff$sample_2,levels=myLevels)
	
	res<-new("CuffGeneSet",
			#TODO: Fix ids so that it only displays those genes in CuffGeneSet
			ids=as.character(dbGetQuery(object@DB,idQuery)),
			annotation=dbGetQuery(object@DB,geneAnnotationQuery),
			fpkm=genes.fpkm,
			diff=genes.diff,
			isoforms=new("CuffFeatureSet",
					annotation=dbGetQuery(object@DB,isoformAnnotationQuery),
					fpkm=isoform.fpkm,
					diff=isoform.diff
					),
			TSS=new("CuffFeatureSet",
					annotation=dbGetQuery(object@DB,TSSAnnotationQuery),
					fpkm=TSS.fpkm,
					diff=TSS.diff
					),
			CDS=new("CuffFeatureSet",
					annotation=dbGetQuery(object@DB,CDSAnnotationQuery),
					fpkm=CDS.fpkm,
					diff=CDS.diff
					)
			)
	end<-dbSendQuery(object@DB,"END;")		
	res
}

setMethod("getGenes",signature(object="CuffSet"),.getGenes)

#Test
#myGenes<-getGenes(a,sample(geneFeatures$gene_id,10))

############
#SQL access
############


################
#Misc Utilities
################
.getLevels<-function(object){
	levelsQuery<-'SELECT s.sample_name FROM samples s ORDER BY s.sample_index ASC'
	levels<-dbGetQuery(object@DB,levelsQuery)$sample_name
	levels
}

setMethod("getLevels",signature(object="CuffSet"),.getLevels)


#####################
#Add FeatureData    #
#####################
.addFeatures<-function(object,features,level="genes",...){
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

#TODO: Add method to purge existing feature data table to allow 'refresh' of feature level data
