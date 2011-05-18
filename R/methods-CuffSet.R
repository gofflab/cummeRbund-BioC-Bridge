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

#make CuffGene objects from a list of gene_ids
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
	


############
#SQL access
############