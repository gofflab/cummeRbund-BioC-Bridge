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
		function(.Object
				){
					
		}
)

##################
#Class Methods
##################
setMethod("show","CuffSet",
		function(object){
			cat(class(object), "instance with:\n")
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
sampleNames<-function(object){
	
}

############
#SQL access
############