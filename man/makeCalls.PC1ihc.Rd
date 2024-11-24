\name{makeCalls.PC1ihc}
\alias{makeCalls.PC1ihc}
\title{Make intermediate intrinsic subtype calls}
\description{
  This function processes clinical IHC subtyping data and preprocessed PAM50 gene expression data
  to form a gene expression-guided ER-balanced set.This set is created by combining IHC classification
  information and using principal component 1 (PC1) to guide the separation.The function computes
  the median for each gene in this ER-balanced set, updates a calibration file, and runs subtype
  prediction algorithms to generate intermediate intrinsic subtype calls based on the PAM50 method.
  Various diagnostics and subtyping results are returned.
}
\usage{
makeCalls.PC1ihc(df.cln, seed = 118, mat, outDir=NULL)
}
\arguments{
  \item{df.cln}{Data frame of clinical data; It should include the columns 'PatientID' and 'IHC'.}
  \item{seed}{Seed for random number generation to ensure reproducibility. Default is 118.}
  \item{mat}{Matrix of preprocessed PAM50 expression data.}
  \item{outDir}{Directory for output files.If \code{NULL}, a subdirectory named "Calls.PC1ihc" within the temporary directory will be used.}
}
\value{
  Returns a list containing:
  \item{Int.sbs}{Data frame with integrated subtype and clinical data.}
  \item{score.fl}{Data frame with scores from subtype predictions.}
  \item{mdns.fl}{Data frame with median values for each gene in the ER-balanced set.}
  \item{SBS.colr}{Colors associated with each subtype from the prediction results.}
  \item{outList}{Detailed results from subtype prediction functions.}
  \item{PC1cutoff}{Cutoff values for PC1 used in subsetting.}
  \item{DF.PC1}{Data frame of initial PCA results merged with clinical data.}
}
\examples{
\donttest{
  data_path <- system.file("extdata", "Sample_IHC_PAM_Mat.Rdat", package = "PCAPAM50")
  load(data_path) # Loads Test.ihc and Test.matrix


  # Prepare the data
  Test.ihc$ER_status <- rep("NA", length(Test.ihc$PatientID))
  Test.ihc$ER_status[grep("^L",Test.ihc$IHC)] = "pos"
  Test.ihc$ER_status[-grep("^L",Test.ihc$IHC)] = "neg"
  Test.ihc <- Test.ihc[order(Test.ihc$ER_status, decreasing = TRUE),]
  Test.matrix <- Test.matrix[, Test.ihc$PatientID]


  df.cln <- data.frame(PatientID = Test.ihc$PatientID, IHC = Test.ihc$IHC, stringsAsFactors = FALSE)


  

  # Call the function
  result <- makeCalls.PC1ihc(df.cln=df.cln, seed = 118, mat = Test.matrix, outDir=NULL)
}
}

\seealso{
  \code{\link[stats]{prcomp}}, \code{\link[base]{merge}}, \code{\link[base]{set.seed}}
}
