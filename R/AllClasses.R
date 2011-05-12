# AllClasses.R
# 
# Author: lgoff
###############################################################################

##diffData Class
#setClassUnion("DiffData",c("list","environment"))
#
##CuffSet Class (extends eSet) - Handles "fpkm_tracking" output files from cufflinks
#setClass("cuffSet",
#		representation(diffData = "DiffData"),
#		contains = "eSet",
#		prototype = prototype(
#				diffData=list(),
#				new("VersionedBiobase",
#						versions=c(classVersion("eSet"),cuffSet="0.1")
#				)
#		)
#)

#CuffData class
setClassUnion("CuffData", c("list", "environment"))

#Main class "CuffSet"
setClass("CuffSet",
		representation(	geneData = "CuffData",
						transcriptData = "CuffData",
						TSSData = "CuffData",
						CDSData = "CuffData",
						geneDiffData = "CuffData",
						transcriptDiffData = "CuffData",
						TSSDiffData = "CuffData",
						CDSDiffData = "CuffData",
						phenoData = "AnnotatedDataFrame",
						featureData = "AnnotatedDataFrame",
						experimentData = "MIAME",
						annotation = "character",
						protocolData = "AnnotatedDataFrame",
						biasData = "AnnotatedDataFrame"
						),
		contains="VersionedBiobase",
		prototype = prototype(
						new("VersionedBiobase", versions=c(CuffSet="0.2.0")),
						geneData = list(),
						transcriptData = list(),
						TSSData = list(),
						CDSData = list(),
						geneDiffData = list(),
						transcriptDiffData = list(),
						TSSDiffData = list(),
						CDSDiffData = list(),
						biasData = list(),
						phenoData = new("AnnotatedDataFrame",
										dimLabels=c("sampleNames","sampleColumns")),
						featureData = new("AnnotatedDataFrame",
								dimLabels=c("sampleNames","sampleColumns")),
						experimentData = new("MIAME"),
						annotation = character(),
						protocolData = new("AnnotatedDataFrame",
								dimLabels=c("sampleNames","sampleColumns")),
						)
				
		)

