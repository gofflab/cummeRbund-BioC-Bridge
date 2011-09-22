\name{plotting-methods}

%CuffData aliases
\alias{csBoxplot,CuffData-method}
\alias{csDensity,CuffData-method}
\alias{csScatter,CuffData-method}
\alias{csVolcano,CuffData-method}

%CuffFeatureSet aliases
\alias{csScatter,CuffFeatureSet-method}
\alias{csVolcano,CuffFeatureSet-method}
\alias{csHeatmap,CuffFeatureSet-method}
\alias{csCluster,CuffFeatureSet-method}
\alias{expressionBarplot,CuffFeatureSet-method}
\alias{expressionPlot,CuffFeatureSet-method}


%CuffFeature aliases
\alias{expressionBarplot,CuffFeature-method}
\alias{expressionPlot,CuffFeature-method}


\title{Plotting methods for cummeRbund classes}

\description{Various plotting methods for cuffdiff data}

\usage{
\S4method{csBoxplot}{CuffData}(object, logMode=TRUE, ...)

\S4method{csDensity}{CuffData}(object, logMode=TRUE, pseudocount=1.0, labels, features=FALSE, ...)

\S4method{csScatter}{CuffData}(object, x, y, logMode=TRUE, pseudocount=1.0, labels, smooth=FALSE, colorByStatus = FALSE, ...)

\S4method{csVolcano}{CuffData}(object, x, y, features=FALSE, ...)

\S4method{csHeatmap}{CuffFeatureSet}(object, rescaling='none', clustering='none', labCol=T, labRow=T, logMode=T, pseudocount=1.0, border=FALSE, heatscale= c(low='darkred',mid='orange',high='white'), heatMidpoint=NULL,...)

\S4method{csCluster}{CuffFeatureSet}(object,k,iter.max=100,...)

\S4method{expressionPlot}{CuffFeature}(object, logMode=FALSE, pseudocount=1.0, drawSummary=FALSE, sumFun=mean_cl_boot, showErrorbars=TRUE, ...)

\S4method{expressionBarplot}{CuffFeatureSet}(object, logMode=TRUE, pseudocount=1.0, showErrorbars=TRUE, ...)

}

\arguments{
	\item{object}{An object of class CuffData, CuffFeature, or CuffFeatureSet}
	
	\item{x}{A character argument for the first condition (sample) in a pairwise comparison}

	\item{y}{A character argument for the second condition (sample) in a pairwise comparison}
	
	\item{logMode}{Logical argument whether or not to log10 transform FPKM values}
	
	\item{pseudocount}{Numerical value added to all FPKMs avoid errors in log10 transfromation of true zero values}
	
	\item{labels}{Not implemented yet}
	
	\item{features}{Logical argument whether or not to include attributes from features table in output}
	
	\item{smooth}{Logical argument whether or not to draw a smoother summary line on a returned scatterplot}
	
	\item{colorByStatus}{Logical argument whether or not to shade significant points a different color}
	
	\item{\dots}{Additional arguments to plotting function(s)}
	
	\item{rescaling}{One of c('none','row','column','both') to rescale data in heatmap}
	
	\item{clustering}{One of c('none','row','column','both') or a function that returns a dist object.  This argument indicates which dimensions of a heatmap are clustered.  If you pass
	a function that returns a 'dist' object, this function will be used to reorder rows and columns.  By default clustering uses the Jensen-Shannon distance as a similarity metric.}
	
	\item{labRow}{Logical argument whether or not to draw labels for Row}
	
	\item{labCol}{See labRow}
	
	\item{drawSummary}{Logical argument whether or not to draw a summary line (using sumFun) for all features in FeatureSet or GeneSet}
	
	\item{sumFun}{Function used to plot summary line}
	
	\item{showErrorbars}{Logical argument whether or not to draw error bars on expression plot [default: TRUE]}
	
	\item{border}{Logical argument whether or not to draw a border around heatmap}
	
	\item{heatscale}{Vector of length 2 or 3 to indicate colors for heatmap scale}
	
	\item{heatMidpoint}{Numerical value as midpoint of heatscale}

	\item{k}{Numerical argument for number of clusters in K-means clustering}
	
	\item{iter.max}{Maximum number of iterations for K-means clustering}

}

\details{
None
}
\value{
A ggplot2 plot object.
}
\references{
None
}
\author{
Loyal A. Goff
}
\note{
None
}
\examples{
	a<-readCufflinks(system.file("extdata", package="cummeRbund")) #Read cufflinks data and create CuffSet object
	genes<-a@genes #CuffData object for all genes
	
	b<-csBoxplot(genes) #Create csBoxplot
	b #Render plot
	
	d<-csDensity(genes) #Create csDensity plot
	d #Render plot
	
	s<-csScatter(genes,'hESC','Fibroblasts',smooth=TRUE) #Create plot object
	s #render plot object
	
	data(sampleData)
	h<-csHeatmap(sampleGeneSet) #Create csHeatmap plot
	h #Render plot
	
	v<-csVolcano(genes,"hESC","Fibroblasts") # Volcano plot of all genes for conditions x='hESC' and y='Fibroblast'
	v #print plot
	
	csCluster(sampleGeneSet,4)
	
}