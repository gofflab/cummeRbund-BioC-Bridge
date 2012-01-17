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

.fpkmMatrix<-function(object){
	res<-fpkm(object)
	colnames(res)[1]<-"tracking_id"
	res<-res[,c(1:3)]
	res<-melt(res)
	res<-cast(res,tracking_id~sample_name)
	res<-data.frame(res[,-1],row.names=res[,1])
}

setMethod("fpkmMatrix",signature(object="CuffFeature"),.fpkmMatrix)

#setMethod("diff","CuffFeature",function(object){
#		return(object@diff)
#		})

.diffData<-function(object){
	object@diff
}

setMethod("diffData",signature(object="CuffFeature"),.diffData)

setMethod("annotation","CuffFeature",function(object){
		return(object@annotation)
		})

#################
#Setters		#
#################


#################
#Plotting		#
#################
.barplot<-function(object,logMode=FALSE,pseudocount=1.0,showErrorbars=TRUE,...){
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
		    geom_errorbar(aes(ymin=conf_lo,ymax=conf_hi,group=1),width=0.5)
	}
	
	if (logMode)
	{
	    p <- p + scale_y_log10()
    }
	
    p <- p + facet_wrap('tracking_id') +
          opts(title=object@annotation$gene_short_name,axis.text.x=theme_text(hjust=0,angle=-90))
		
    # p <- p + facet_wrap('tracking_id')
    #     # gene_labels <- object@annotation$gene_short_name
    #     # gene_labels[is.na(object@annotation$gene_short_name)] = dat[is.na(object@annotation$gene_short_name),1]
    #     # print("gene_labels:")
    #     # print(str(gene_labels))
    #     # p <- p + opts(title=gene_labels,axis.text.x=theme_text(hjust=0,angle=-90))   
	
    # p <- p + ylim(min(dat$conf_lo), max(dat$conf_hi))
	
    if (logMode)
    {
        p <- p + ylab(paste("FPKM +",pseudocount))
    } else {
        p <- p + ylab("FPKM")
    }
	
	
	p <- p + opts(legend.position="none")
	
	#Default cummeRbund colorscheme
	p<-p + scale_fill_hue(l=50,h.start=200)
	
	p
}

setMethod("expressionBarplot",signature(object="CuffFeature"),.barplot)


.expressionPlot<-function(object,logMode=FALSE,pseudocount=1.0, drawSummary=FALSE, sumFun=mean_cl_boot, showErrorbars=TRUE,...){
	dat<-fpkm(object)
	colnames(dat)[1]<-"tracking_id"
	if(logMode)
	{
	    dat$fpkm <- dat$fpkm + pseudocount
	    dat$conf_hi <- dat$conf_hi + pseudocount
	    dat$conf_lo <- dat$conf_lo + pseudocount
    }
	p <- ggplot(dat)
	#dat$fpkm<- log10(dat$fpkm+pseudocount)
	p <- p + 
	    geom_line(aes(x=sample_name,y=fpkm,color=tracking_id,group=tracking_id))
	if (showErrorbars)
	{
	    p <- p +
		    geom_errorbar(aes(x=sample_name, ymin=conf_lo,ymax=conf_hi, color=tracking_id, group=tracking_id),width=0.25)
	}
	
	if (logMode)
	{
	    p <- p + scale_y_log10()
    }

	
	#drawMean
	if(drawSummary){
		p <- p + stat_summary(aes(x=sample_name,y=fpkm,group=1),fun.data=sumFun,color="red",fill="red",alpha=0.2,size=1.1,geom="smooth")
	}
	
	if (logMode)
    {
        p <- p + ylab(paste("Log10 FPKM + ",pseudocount))
    } else {
        p <- p + ylab("FPKM")
    }
	
	#Default cummeRbund colorscheme
	p<-p + scale_fill_hue(l=50,h.start=200) + scale_color_hue(l=50,h.start=200)
	
	#Add Title
	p<-p + opts(title=object@annotation$gene_short_name,axis.text.x=theme_text(hjust=0,angle=-90))
	
	p
}

setMethod("expressionPlot",signature(object="CuffFeature"),.expressionPlot)
#################
#Misc			#
#################