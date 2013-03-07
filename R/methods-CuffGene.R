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

#features
setMethod("features","CuffGene",function(object){
			return(object@features)
		})

#################
#Plotting		#
#################
.makeGeneRegionTrack<-function(object){
	featCols<-c('seqnames','start','end','source','gene_id','exon_number','isoform_id','isoform_id','exon_number','strand')
	feats<-features(object)[,featCols]
	newColnames<-c('seqnames','start','end','feature','gene','exon','transcript','symbol','rank','strand')
	mychr<-unique(feats$seqnames)
	colnames(feats)<-newColnames
	feats<-feats[,-1]
	#print(feats)
	feats$symbol[is.na(feats$symbol)]<-"NA"
	#THIS NEEDS TO BE MADE GENERIC
	genetrack<-GeneRegionTrack(feats,genome=object@genome,chromosome=mychr,name='CuffDiff',showId=T,stacking="pack")
	genetrack
}

setMethod("makeGeneRegionTrack",signature(object="CuffGene"),.makeGeneRegionTrack)

.plot<-function(object){
	trackList<-list()
	myStart<-min(object@features$start)
	myEnd<-max(object@features$end)
	#Make the following conditional on network connectivity
	ideoTrack <- IdeogramTrack(genome = object@genome, chromosome = unique(object@features$seqnames))
	trackList<-c(trackList,ideoTrack)
	
	axtrack<-GenomeAxisTrack()
	trackList<-c(trackList,axtrack)
	
	genetrack<-.makeGeneRegionTrack(object)
	
	trackList<-c(trackList,genetrack)
	
	biomTrack<-BiomartGeneRegionTrack(genome=object@genome,chromosome=as.character(unique(object@features$seqnames)),start=myStart,end=myEnd,name="ENSEMBL",showId=T)
	trackList<-c(trackList,biomTrack)
	
	
	plotTracks(trackList,from=myStart-2000,to=myEnd+2000)
}

setMethod("genePlot",signature(object="CuffGene"),.plot)



#################
#Feature Plotting
#################
.ideogram<-function(object){
	myStart<-min(object@features$start)
	myEnd<-max(object@features$end)
	mychr<-unique(object@features$seqnames)
	p<-plotIdeogram(genome=object@genome,subchr=mychr,zoom.region=c(myStart,myEnd))
	p
	
}

.plot2<-function(object,...){
	#Ideogram
	ideoTrack<-.ideogram(object)
	
	#Expression levels
	expressionTrack<-expressionPlot(isoforms(object),facet=T,...)+theme_bw() + theme(legend.position='none')
	hasAxis(expressionTrack)<-TRUE
	
	#Transcript Models
	#modelTrack<-autoplot(.asGRangesList(object),aes(fill=transcript,group=transcript),gap.geom="arrow") + theme_bw() + scale_fill_hue(l=50,h.start=200) + scale_color_hue(l=50,h.start=200)
	modelTrack<-ggplot(as.GRangesList(object),)
	
	hasAxis(modelTrack)<-TRUE
	
	#Plot it all...
	tracks(ideoTrack,modelTrack,expressionTrack,heights=c(1,3,3),fixed=c(TRUE,TRUE,FALSE),main=unique(object@annotation$gene_short_name))
	
}

.pie<-function(object,level="isoforms",pseudocount=0.0001,...){
	dat<-fpkm(slot(object,level))
	colnames(dat)[1]<-'tracking_id'
	#dat$fpkm<-dat$fpkm+pseudocount
	#print(dat)
	p<-ggplot(dat,aes(x="",y=fpkm,fill=tracking_id))
	
	p<- p + geom_bar(stat="identity",position="fill",line="black")
	
	p<- p + coord_polar(theta='y')
	
	p<-p + scale_fill_hue(l=50,h.start=200) + scale_color_hue(l=50,h.start=200) + theme_bw()
	
	p<-p + facet_wrap('sample_name') + theme(axis.text.x = element_blank(),aspect.ratio=1)
	
	p
}

setMethod("csPie",signature(object="CuffGene"),.pie)

#################
#Coersion methods
#################
#As GRanges
.as.GRanges<-function(from){
	#featCols<-c('seqnames','start','end','source','gene_id','exon_number','isoform_id','isoform_id','exon_number','strand')
	feats<-from@features
	#newColnames<-c('seqnames','start','end','feature','gene','exon','transcript','symbol','rank','strand')
	#ExpressionValues (transcript)
	fpkm<-cbind(isoform_id=rownames(fpkmMatrix(isoforms(from))),fpkmMatrix(isoforms(from)))
	feats<-merge(feats,fpkm)
	colnames(feats)[colnames(feats)=='isoform_id']<-'transcript'
	colnames(feats)[colnames(feats)=='gene_id']<-'gene'
	colnames(feats)[colnames(feats)=='exon_number']<-'exon'
	colnames(feats)[colnames(feats)=='source']<-'feature'
	feats$symbol<-feats$transcript
	corCols<-c('seqnames','start','end','strand','width')
	myGR<-GRanges(Rle(feats$seqnames),ranges=IRanges(feats$start,end=feats$end),strand=Rle(feats$strand),elementMetadata=feats[,!colnames(feats) %in% corCols])
	colnames(elementMetadata(myGR))<-colnames(feats[,!colnames(feats) %in% corCols])
	myGR
}

setAs("CuffGene","GRanges",.as.GRanges)

#As GRangesList
.as.GRangesList<-function(object,f="transcript"){
	gr<-as(object,"GRanges")
	grl<-split(gr,f)
	grl
}
setMethod("as.GRangesList",signature(object="CuffGene"),.as.GRangesList)


#######################
# ggbio integration
#######################
#setMethod("ggplot", "CuffGene", function(data, ...){
#			df <- mold(as.GRangesList(data))
#			g <- ggplot(df, ...)
#			g$.data <- as.GRangesList(data)
#			g <- ggbio(g)
#			g
#		})