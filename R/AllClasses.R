# AllClasses.R
# 
# Author: lgoff
###############################################################################


#CuffData Class
setClass("CuffData",
		representation(	DBfile = "character",
						
						mainTable = "character",
						diffTable = "character",
						expDiffTable = "character",
						filterList = "list"
						)
		)

#CuffSet Class
setClass("CuffSet",
		representation(	DB = "dbConnect",
						conditions = "data.frame",
						genes = "CuffData",
						isoforms = "CuffData",
						TSS = "CuffData",
						CDS = "CuffData"
						),
		prototype = prototype(
				)

		)