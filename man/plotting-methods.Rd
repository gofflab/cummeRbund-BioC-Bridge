\name{plotting-methods}

%CuffData aliases
\alias{csBoxplot,CuffData-method}
\alias{csDensity,CuffData-method}
\alias{csScatter,CuffData-method}
\alias{csVolcano,CuffData-method

%CuffFeatureSet aliases
\alias{csScatter,CuffFeatureSet-method}
\alias{csVolcano,CuffFeatureSet-method}
\alias{csHeatmap,CuffFeatureSet-method}


\title{Plotting methods for cummeRbund classes}

\description{Various plotting methods for cuffdiff data}

\usage{
\S4method{csBoxplot}{CuffData}(object, logMode=TRUE, ...)

\S4method{csDensity}{CuffData}(object, logMode=TRUE, pseudocount=1.0, labels, features=FALSE, ...)

\S4method{csScatter}{CuffData}(object, x, y, logMode=TRUE, pseudocount=1.0, labels, smooth=FALSE, colorByStatus = FALSE, ...)

\S4method{csVolcano}{CuffData}(object, x, y, features=FALSE, ...)

\S4method{csHeatmap}{CuffFeatureSet}(object, rescaling='none', clustering='none', labCol=T, labRow=T, logMode=T, pseudocount=1.0, 
		border=FALSE, heatscale= c(low='darkred',mid='orange',high='white'), heatMidpoint=NULL,...)

}

\arguments{
	\item{object}{}

}