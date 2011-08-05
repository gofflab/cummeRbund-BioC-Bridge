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
			TRUE #length(object)==1
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
.barplot<-function(object,logMode=FALSE,pseudocount=0.0001,showErrorbars=TRUE,...){
	dat<-fpkm(object)
	#TODO: Test dat to ensure that there are >0 rows to plot.  If not, trap error and move on...
	
	colnames(dat)[1]<-"tracking_id"
	
	if(logMode)
	{
	    dat$fpkm <- dat$fpkm + pseudocount
	    dat$conf_hi <- dat$conf_hi + pseudocount
	    dat$conf_lo <- dat$conf_lo + pseudocount
    }

    p<-ggplot(dat,aes(x=sample_name,y=fpkm,fill=sample_name))
    
	#dat$fpkm<- log10(dat$fpkm+pseudocount)
	p <- p + 
	    geom_bar()
	if (showErrorbars)
	{
	    p <- p +
		    geom_errorbar(aes(ymin=conf_lo,ymax=conf_hi,group=1))
	}
	
	if (logMode)
	{
	    p <- p + scale_y_log10()
    }
	
	p <- p + facet_wrap('tracking_id') +
    	opts(title=object@annotation$gene_short_name,axis.text.x=theme_text(hjust=0,angle=-90))
		
    # p <- p + ylim(min(dat$conf_lo), max(dat$conf_hi))
	
	p <- p + opts(legend.position = "none")
	p <- p + scale_fill_brewer(palette="Set1")
	p
}

setMethod("expressionBarplot",signature(object="CuffFeature"),.barplot)


.expressionPlot<-function(object,logMode=FALSE,pseudocount=0.0001, drawSummary=FALSE, sumFun=mean_cl_boot, showErrorbars=T,...){
	dat<-fpkm(object)
	colnames(dat)[1]<-"tracking_id"
	p <- ggplot(dat)
	if(logMode)
	{
	    dat$fpkm <- dat$fpkm + pseudocount
	    dat$conf_hi <- dat$conf_hi + pseudocount
	    dat$conf_lo <- dat$conf_lo + pseudocount
    }

	#dat$fpkm<- log10(dat$fpkm+pseudocount)
	p <- p + 
	    geom_line(aes(x=sample_name,y=fpkm,color=tracking_id,group=tracking_id))
	if (showErrorbars)
	{
	    p <- p +
		    geom_errorbar(aes(x=sample_name, ymin=conf_lo,ymax=conf_hi, color=tracking_id, group=tracking_id))
	}
	
	if (logMode)
	{
	    p <- p + scale_y_log10()
    }

	
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