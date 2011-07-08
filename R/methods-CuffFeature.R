########################
#methods-CuffFeature.R
#
#Author: Loyal A. Goff
#
#Date created: 5-17-2011
#
#Description: A 'data' object for a collection of cufflinks features, irrespective of type ('genes','isoforms','TSS','CDS')
#########################

#################
#Initialize		#
#################
setMethod("initialize","CuffFeature",
		function(.Object,
				annotation=data.frame(),
				fpkm=data.frame(),
				diff=data.frame(),
				... ){
			.Object<-callNextMethod(.Object,
					annotation=annotation,
					fpkm=fpkm,
					diff=diff,
					...)				
		}
)

#################
#Validate		#
#################
#TODO: Add validity constraints
setValidity("CuffFeature",function(object){
			TRUE
		}
)		

#################
#Class Methods	#
#################
setMethod("show","CuffFeature",
		function(object){
			cat(class(object), "instance with ",length(object),"elements\n")
		}
)

setMethod("length","CuffFeature",
		function(x){
			dim(x@annotation)[1]
		}
)
#################
#Subsetting		#
#################


#################
#Accessors		#
#################
.fpkm<-function(object){
	object@fpkm
}
setMethod("fpkm",signature="CuffFeature",.fpkm)

#################
#Setters		#
#################


#################
#Plotting		#
#################
.barplot<-function(object,logMode=TRUE,pseudocount=0.0001,...){
	dat<-fpkm(object)
	#TODO: Test dat to ensure that there are >0 rows to plot.  If not, trap error and move on...
	
	colnames(dat)[1]<-"tracking_id"
	p<-ggplot(dat,aes(x=sample_name,y=fpkm,fill=sample_name))
	
	p<- p +
		geom_bar() +
		geom_errorbar(aes(ymin=conf_lo,ymax=conf_hi,group=1),size=0.15) +
		facet_wrap('tracking_id') +
		opts(title=object@annotation$gene_short_name,axis.text.x=theme_text(hjust=0,angle=-90))
	
	#This does not make immediate sense with the conf_hi and conf_lo values.  Need to figure out appropriate transformation for these
	#if(logMode)
		#p<-p+scale_y_log2()
	p + opts(legend.position = "none")
	
}

setMethod("expressionBarplot",signature(object="CuffFeature"),.barplot)


.expressionPlot<-function(object,logMode=TRUE,pseudocount=0.0001, drawSummary=FALSE, sumFun=mean_cl_boot,...){
	dat<-fpkm(object)
	colnames(dat)[1]<-"tracking_id"
	if(logMode){
		dat$fpkm<-log2(dat$fpkm+pseudocount)
	}
	p	<-	ggplot(dat) + 
			geom_line(aes(x=sample_name,y=fpkm,color=tracking_id,group=tracking_id)) 
	
	#drawMean
	if(drawSummary){
		p <- p + stat_summary(aes(x=sample_name,y=fpkm,group=1),fun.data=sumFun,color="red",fill="red",alpha=0.2,size=1.1,geom="smooth")
	}
	
	p
}

setMethod("expressionPlot",signature(object="CuffFeature"),.expressionPlot)
#################
#Misc			#
#################