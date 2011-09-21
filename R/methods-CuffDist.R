#methods-CuffDist.R
#
#Author: Loyal A. Goff
#
#
####################

##################
#Initialize
##################
setMethod("initialize","CuffDist",
			function(.Object,
					DB,
					table="",
					type = c("promoter","splicing","relCDS"),
					testId = c("gene_id","tss_group_id"),
					... ){
				.Object<-callNextMethod(.Object,
						DB = DB,
						table = table,
						type = type,
						testId = testId,
						...)				
		}
)

setValidity("CuffDist",function(object){
		TRUE
		}
)			

################
#Class Methods
################
setMethod("show","CuffDist",
		function(object){
			size<-dim(object)
			cat(class(object), "instance with:\n\t",size[1]," ",object@type," records\n")
		}
)

setMethod("dim","CuffDist",
		function(x){
			countQuery<-paste("SELECT COUNT(",x@testId,") as n FROM ",x@table)
			nIds<-dbGetQuery(x@DB,countQuery)
			c(nIds$n)
		}
)

###################
#Accessors
###################
.values<-function(object){
	valueQuery<-paste("SELECT * FROM ",object@table,sep="")
	dbGetQuery(object@DB, valueQuery)
}

setMethod("distValues","CuffDist",.values)

setMethod("DB","CuffDist",function(object){
		return(object@DB)
		})

#setMethod("table","CuffDist",function(object){
#		return(object@table)
#		})

setMethod("type","CuffDist",function(object){
		return(object@type)
		})

setMethod("testId","CuffDist",function(object){
		return(object@testId)
		})

##################
#Setters
##################


##################
#Subsetting
##################


##################
#Plotting
##################

