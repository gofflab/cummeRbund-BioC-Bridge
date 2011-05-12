# TODO: Add comment
# 
# Author: lgoff
###############################################################################

#diffData Class
setClassUnion("DiffData",c("list","environment"))

#CuffSet Class (extends eSet) - Handles "fpkm_tracking" output files from cufflinks
setClass("cuffSet",
		representation(diffData = "DiffData"),
		contains = "eSet",
		prototype = prototype(
				diffData=list(),
				new("VersionedBiobase",
						versions=c(classVersion("eSet"),cuffSet="0.1")
				)
		)
)
