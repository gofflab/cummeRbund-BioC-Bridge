# AllClasses.R
# 
# Author: lgoff
###############################################################################

#TODO: I get the distinct feeling that these two should be nested environments, but I don't really know what that means.

#CuffData Class
setClass("CuffData",
		representation(DB = "SQLiteConnection",
						tables = "list",
						filters = "list",
						type = "character",
						idField = "character"
						)
		)

#CuffSet Class
setClass("CuffSet",
		representation(DB = "SQLiteConnection",
						conditions = "data.frame",
						genes = "CuffData",
						isoforms = "CuffData",
						TSS = "CuffData",
						CDS = "CuffData"
						)#,
#		prototype = prototype(
#				conditions = data.frame(),
#				genes = new("CuffData", DB=dbConnect(dbDriver("SQLite"),"cuffData.db"), tables = list(mainTable = "genes",dataTable = "geneData",expDiffTable = "geneExpDiff",otherTable = "promoterDiffData"), filters = list(),type = "genes",idField = "gene_id"),
#				isoforms = new("CuffData", DB=dbConnect(dbDriver("SQLite"),"cuffData.db"), tables = list(mainTable = "isoforms",dataTable = "isoformData",expDiffTable = "isoformExpDiff",otherTable = ""), filters = list(),type="isoforms",idField = "isoform_id"),
#				TSS = new("CuffData", DB = DB=dbConnect(dbDriver("SQLite"),"cuffData.db"), tables = list(mainTable = "TSS",dataTable = "TSSData",expDiffTable = "TSSExpDiff", otherTable = "splicingDiffData"), filters = list(),type = "TSS",idField = "TSS_id"),
#				CDS = new("CuffData", DB = DB=dbConnect(dbDriver("SQLite"),"cuffData.db"), tables = list(mainTable = "CDS",dataTable = "CDSData",expDiffTable = "CDSExpDiff", otherTable = "CDSDiffData"), filters = list(),type = "CDS",idField = "CDS_id")
#				)
)