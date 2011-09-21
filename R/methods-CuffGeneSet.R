########################
#methods-CuffGeneSet.R
#
#Author: Loyal A. Goff
#
#Date created: 5-17-2011
#
#Description: Defines a class of cufflinks data for multiple genes
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
setMethod("show","CuffGeneSet",function(object){
			cat(class(object),"instance for genes",object@ids,"\nShort name:\t",unique(object@annotation$gene_short_name),
					"\nSlots:\n\t annotation\n\t fpkm\n\t diff\n\t",
					"isoforms\t",class(object@isoforms),"instance of size",length(object@isoforms),"\n\t",
					"TSS\t\t",class(object@TSS),"instance of size",length(object@TSS),"\n\t",
					"CDS\t\t",class(object@CDS),"instance of size",length(object@CDS),"\n"
			)			
		}
)

#################
#Accessors
#################
#isoforms
setMethod("isoforms","CuffGeneSet",function(object){
			return(object@isoforms)	
		})
#TSS
setMethod("TSS","CuffGeneSet",function(object){
			return(object@TSS)
		})

#CDS
setMethod("CDS","CuffGeneSet",function(object){
			return(object@CDS)
		})


#################
#Subsetting		#
#################


#################
#Plotting		#
#################


#################
#Misc			#
#################