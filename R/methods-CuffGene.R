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


#################
#Class Methods	#
#################
setMethod("show","CuffGene",function(object){
		cat(class(object),"instance for gene",object@id,"\nShort name:\t",unique(object@annotation$gene_short_name),
						"\nSlots:\n\t annotation\n\t fpkm\n\t diff\n\t",
						"isoforms\t",class(object@isoforms),"instance of size",length(object@isoforms),"\n\t",
						"TSS\t\t",class(object@TSS),"instance of size",length(object@TSS),"\n\t",
						"CDS\t\t",class(object@CDS),"instance of size",length(object@CDS),"\n"
						)			
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


#################
#Plotting		#
#################


#################
#Misc			#
#################