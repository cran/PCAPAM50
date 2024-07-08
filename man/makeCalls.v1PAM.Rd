\name{makeCalls.v1PAM}
\alias{makeCalls.v1PAM}
\title{Make PCAPAM50 calls}
\description{
  This function uses the intermediate intrinsic subtype calls and preprocessed PAM50 gene expression
  data to create an ER-balanced set and produces PCAPAM50 Calls.
}
\usage{
makeCalls.v1PAM(df.pam, seed = 118, mat, inputDir=NULL)
}
\arguments{
  \item{df.pam}{Data frame of PAM data; It should include the columns 'PatientID' and 'PAM50'.}
  \item{seed}{Seed for random number generation to ensure reproducibility.}
  \item{mat}{Matrix of preprocessed PAM50 expression data.}
  \item{inputDir}{Directory for output files.If \code{NULL}, a subdirectory named "Calls.PCAPAM50" within the temporary directory will be used.}
}
\value{
  Returns a list containing:
  \item{Int.sbs}{Data frame with integrated subtype and clinical data.}
  \item{score.fl}{Data frame with scores from subtype predictions.}
  \item{mdns.fl}{Data frame with median values for each gene in the ER-balanced set.}
  \item{SBS.colr}{Colors associated with each subtype from the prediction results.}
  \item{outList}{Detailed results from subtype prediction functions.}
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
  inputDir <- "Call.PC1"

  # Make a secondary ER-balanced subset and derive intermediate intrinsic subtype calls
  result <- makeCalls.PC1ihc(df.cln=df.cln, seed = 118, mat = Test.matrix, inputDir=inputDir)

  df.pc1pam = data.frame(PatientID=result$Int.sbs$PatientID,
  			PAM50=result$Int.sbs$Int.SBS.Mdns.PC1ihc,
			IHC=result$Int.sbs$IHC,
			stringsAsFactors=FALSE) ### IHC column is optional
  
  

  # Make a tertiary ER-balanced set and PCAPAM50 calls
  res <- makeCalls.v1PAM(df.pam = df.pc1pam, seed = 118, mat = Test.matrix, inputDir=NULL)
}
}
\seealso{
  \code{\link[stats]{prcomp}}, \code{\link[base]{merge}}, \code{\link[base]{set.seed}}
}
