########################
#methods-CuffGene.R
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


#################
#Validate		#
#################
#TODO: Add validity constraints
setValidity("CuffGene",function(object){
			objLen = length(object)
			if(objLen==0){
				write("No gene set returned (Gene might not be in database?)",stderr())
				return(FALSE)
			}
			if(objLen>1){
				write("Warning: Possibly more than one gene returned",stderr())
				return(TRUE)
			}
		}
)		

#################
#Class Methods	#
#################
setMethod("show","CuffGene",function(object){
		cat(class(object),"instance for gene",object@id,"\nShort name:\t",unique(object@annotation$gene_short_name),
						"\nSlots:\n\t annotation\n\t features\n\t fpkm\n\t repFpkm\n\t diff\n\t count\n\t",
						"isoforms\t",class(object@isoforms),"instance of size",length(object@isoforms),"\n\t",
						"TSS\t\t",class(object@TSS),"instance of size",length(object@TSS),"\n\t",
						"CDS\t\t",class(object@CDS),"instance of size",length(object@CDS),"\n"
						)			
		}
)

setMethod("length","CuffGene",
		function(x){
			dim(x@annotation)[1]
		}
)

#################
#Subsetting		#
#################

#################
#Accessors
#################
#isoforms
setMethod("isoforms","CuffGene",function(object){
		return(object@isoforms)	
		})
#TSS
setMethod("TSS","CuffGene",function(object){
		return(object@TSS)
		})

#CDS
setMethod("CDS","CuffGene",function(object){
		return(object@CDS)
		})

#promoters
setMethod("promoters","CuffGene",function(object){
			return(object@promoters)
		})
#splicing
setMethod("splicing","CuffGene",function(object){
			return(object@splicing)
		})
#relCDS
setMethod("relCDS","CuffGene",function(object){
			return(object@relCDS)
		})

#################
#Plotting		#
#################
.makeGeneRegionTrack<-function(object){
	featCols<-c('seqnames','start','end','source','gene_id','exon_number','transcript_id','transcript_id','exon_number','strand')
	feats<-object@features[featCols]
	newColnames<-c('seqnames','start','end','feature','gene','exon','transcript','symbol','rank','strand')
	mychr<-unique(feats$seqnames)
	colnames(feats)<-newColnames
	feats<-feats[,-1]
	#print(feats)
	feats$symbol[is.na(feats$symbol)]<-"NA"
	#THIS NEEDS TO BE MADE GENERIC
	genetrack<-GeneRegionTrack(feats,genome="hg19",chromosome=mychr,name='CuffDiff',showId=T,stacking="pack")
	genetrack
}

.plot<-function(object){
	trackList<-list()
	myStart<-min(object@features$start)
	myEnd<-max(object@features$end)
	#Make the following conditional on network connectivity
	ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = unique(object@features$seqnames))
	trackList<-c(trackList,ideoTrack)
	
	axtrack<-GenomeAxisTrack()
	trackList<-c(trackList,axtrack)
	
	genetrack<-.makeGeneRegionTrack(object)
	
	trackList<-c(trackList,genetrack)
	
	biomTrack<-BiomartGeneRegionTrack(genome="hg19",chromosome=as.character(unique(object@features$seqnames)),start=myStart,end=myEnd,name="ENSEMBL",showId=T)
	trackList<-c(trackList,biomTrack)
	
	
	plotTracks(trackList,from=myStart-2000,to=myEnd+2000)
}

setMethod("genePlot",signature(object="CuffGene"),.plot)

#################
#Coersion methods
#################
#As GRanges

#As GRangesLit