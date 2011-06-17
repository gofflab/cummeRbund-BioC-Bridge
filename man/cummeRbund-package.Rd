\name{cummeRbund-package}
\alias{cummeRbund-package}
\alias{cummeRbund}
\docType{package}
\title{
cummeRbund: The finishing touch on your Tuxedo workflow. Analysis, manipulation, and visualization of Cufflinks HTS data.
}
\description{
Allows for persistent storage, access, and manipulation of Cufflinks high-throughput sequencing data.  In addition, provides numerous plotting functions for commonly used visualizations.
}
\details{
\tabular{ll}{
Package: \tab cummeRbund\cr
Version: \tab 0.1.2\cr
Suggests: \tab \cr
Depends: \tab R (>= 2.7.0), RSQLite, reshape, ggplot2, methods\cr
License: \tab MIT License\cr
Collate: \tab AllGenerics.R
AllClasses.R
database-setup.R
methods-CuffSet.R
methods-CuffData.R
methods-CuffGeneSet.R
methods-CuffFeatureSet.R
methods-CuffGene.R
methods-CuffFeature.R
tools.R\cr
LazyLoad: \tab yes\cr
biocViews: \tab Infrastructure, Bioinformatics, Visualization\cr
Packaged: \tab 2011-06-06 17:28:38 UTC; lgoff\cr
Built: \tab R 2.12.2; ; 2011-06-06 17:53:09 UTC; unix\cr
}

Index:
\preformatted{
JSdist                  Jensen-Shannon distance on columns
csDensity               Density plot of CuffData
csScatter               Scatter Plot
csVolcano               Volcano Plot
diffData                Differential comparison data
expressionBarplot       Barplot
featureNames            Feature names
features                Features
fpkm                    Retrieve FPKM values
fpkmMatrix              Retrieve FPKM values as matrix
makeprobs               Transform a matrix into probabilities by
                        columns
samples                 Get sample list from CuffData object
shannon.entropy         Shannon entropy
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{cummeRbund-manual} \tab An R package for visualization and analysis of Cufflinks high-throughput sequencing data (source, pdf)\cr
}
}
\author{
L. Goff, C. Trapnell

Maintainer: Loyal A. Goff <lgoff@csail.mit.edu>
}
\keyword{ package }
