########################
#methods-CuffFeatureSet.R
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


#################
#Class Methods	#
#################
#setMethod("show","CufFFeatureSet",function(object){
#		cat(class(object),"instance")
#		}
#)

#################
#Subsetting		#
#################

#setMethod("[","CuffFeatureSet",function(object,featureID){
#			
#		}
#)

#################
#Accessors
##################
.fpkm<-function(object){
	object@fpkm
}
setMethod("fpkm",signature="CuffFeatureSet",.fpkm)


#################
#Plotting		#
#################
.heatmap<-function(object,logMode=TRUE,pseudocount=0.0001){
	dat<-fpkm(object)
	if(logMode){
		dat$fpkm<-log2(dat$fpkm+pseudocount)
	}
	colnames(dat)[1] <- "tracking_id"
	p<-ggplot(dat)
	p <- p + geom_tile(aes(x=tracking_id,y=sample_name,fill=fpkm)) + scale_fill_gradient(low="white",high="red") + opts(axis.text.x=theme_text(angle=-90, hjust=0))
	p
}

setMethod("csHeatmap",signature("CuffFeatureSet"),.heatmap)

#################
#Misc			#
#################