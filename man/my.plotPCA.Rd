\name{my.plotPCA}
\alias{my.plotPCA}
\title{Modeling after plotPCA of DESeq}
\description{Modeling after plotPCA of DESeq}
\usage{
my.plotPCA(x, intgroup, ablne = 0,
           colours = c("red","hotpink","darkblue", "lightblue","red3","hotpink3",
           "royalblue3","lightskyblue3"),
           LINE.V = TRUE)
}
\arguments{
  \item{x}{An ExpressionSet object, with matrix data (x) in ‘assay(x)’, produced for example by ExpressionSet(assayData=Test.matrix, phenoData=phenoData)}
  \item{intgroup}{Subtype condition: a character vector of names in ‘colData(x)’ to use for grouping.}
  \item{ablne}{An x-axis coordinate for drawing a vertical line. Default is 0.}
  \item{colours}{Colors for subtypes present in the condition. }
  \item{LINE.V}{Determines whether or not to draw line. Default is \code{TRUE}.}
}
\value{
  Returns an image containing:
  \item{pcafig}{The plot.}
}
\examples{
  library("Biobase")  
  
  data_path <- system.file("extdata", "Sample_IHC_PAM_Mat.Rdat", package = "PCAPAM50")
  load(data_path) # Loads Test.ihc and Test.matrix
  
  pData = data.frame(condition=Test.ihc$IHC)
  rownames(pData) = Test.ihc$PatientID
  phenoData = new("AnnotatedDataFrame", data=pData)#, varMetadata=Metadata
  XSet      = ExpressionSet(assayData=Test.matrix, phenoData=phenoData)
  my.plotPCA(XSet, intgroup=pData$condition, ablne=2.4,
  		colours = c("hotpink","darkblue","lightblue","lightblue3","red"),
  		LINE.V = TRUE)

}
\seealso{
  \code{\link[stats]{prcomp}}, \code{\link[base]{merge}}, \code{\link[base]{set.seed}}
}
