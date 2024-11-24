library(testthat)
library(PCAPAM50)
library(ComplexHeatmap)

test_that("makeCalls.PC1ihc works correctly with test data", {
  # Load the test data
  data_path <- system.file("extdata", "Sample_IHC_PAM_Mat.Rdat", package = "PCAPAM50")
  load(data_path) # Loads Test.ihc and Test.matrix
  expect_true(file.exists(data_path), info = "Data file not found!")
  
    
  # Prepare the data
  Test.ihc$ER_status <- rep("NA", length(Test.ihc$PatientID))
  Test.ihc$ER_status[grep("^L",Test.ihc$IHC)] = "pos"
  Test.ihc$ER_status[-grep("^L",Test.ihc$IHC)] = "neg"
  Test.ihc <- Test.ihc[order(Test.ihc$ER_status, decreasing = TRUE),]  
   
  df.cln <- data.frame(PatientID = Test.ihc$PatientID, IHC = Test.ihc$IHC, stringsAsFactors = FALSE)
  
  
  #outDir <- "Call.PC1"
  outDir <- file.path(tempdir(), "Call.PC1")
  dir.create(outDir, showWarnings = FALSE)  # Ensure the directory exists
  
   
  
  
   
  # Call the function
  result <- makeCalls.PC1ihc(df.cln=df.cln, seed = 118, mat = Test.matrix, outDir=outDir)
  
  # Clean up the temporary directory for testcase
  unlink(outDir, recursive = TRUE)
    
  # Test the result
  expect_type(result, "list")
  expect_true("Int.sbs" %in% names(result))
  expect_true("score.fl" %in% names(result))
  expect_true("mdns.fl" %in% names(result))
  expect_true("SBS.colr" %in% names(result))
  expect_true("outList" %in% names(result))
  expect_true("PC1cutoff" %in% names(result))
  expect_true("DF.PC1" %in% names(result))
})
