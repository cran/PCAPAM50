#' Sample IHC and Matrix Data for Testing
#'
#' These files are used for testing purposes in the PCA-PAM50 approach.
#' - `Sample_IHC_PAM-Mat.Rdat`: Contains `Test.ihc` and `Test.matrix`.
#' - `TEST_pam_mat.txt`: Contains test data for PAM50 predictions.
#'
#' @name TestData
#' @docType data
#' @keywords data
#' @examples
#' # Example for loading sample test data
#' data_path <- system.file("extdata", "Sample_IHC_PAM-Mat.Rdat", package = "pcapam50")
#' load(data_path)
#' 
#' # Example for loading TCGA test data
#' tcga_data_path <- system.file("extdata", "TCGA.712BC_IHC_PAM-Mat.Rdat", package = "pcapam50")
#' load(tcga_data_path)
NULL

