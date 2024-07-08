#--last updated 5/29/2024 3:57 PM

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("
========================================
PCAPAM50 version 1.0.0
Documentation: https://www.wriwindber.org/tools-portal/pca-pam50/

If you use it in published research, please cite:
PCAPAM50: Raj-Kumar PK, et al.Scientific reports. 2019 May 28;9(1):7956.
PAM50: Parker JS, et al. Journal of clinical oncology. 2009 Mar 3;27(8):1160.

This message can be suppressed by:
  suppressPackageStartupMessages(library(PCAPAM50))

========================================
")
}


###### 
##
my.plotPCA = function(x, intgroup, ablne = 0,colours = c("red","hotpink","darkblue","lightblue","red3","hotpink3","royalblue3","lightskyblue3"),LINE.V=TRUE)
{


    # Assuming 'exprs' is correctly defined or imported; using Biobase::exprs explicitly if from Biobase package

    fac = factor(intgroup)
    
    if (nlevels(fac) == 0) {
        stop("intgroup has no valid levels. Please check your 'intgroup' data.")
    }

  
  # Check if number of colors matches number of factor levels in intgroup
    if (length(colours) < nlevels(fac)) {
        stop("Please provide enough colors to match the number of levels in the intgroup factor.")
    }
    
    # Perform PCA
    pca = prcomp(t(Biobase::exprs(x)))
    
    # Define key text based on factor levels
    key_text <- levels(fac)
  
    
   # Generate the PCA plot using lattice::xyplot
    pcafig = lattice::xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x), xlab = "PC1",
                             panel = function(x, y, ...) {
                                 lattice::panel.xyplot(x, y, ...)
                                 if (LINE.V) {
                                     lattice::panel.abline(v = ablne, lty = 3, lwd = 1.5)
                                 }
                             },
                             key = list(space = "right",
                                        text = list(labels = levels(fac)),  
                                        points = list(pch = 16, cex = 2, col = colours),
                                        columns = min(nlevels(fac), 2)),
                                        pch = 16, 
                             cex = 2, 
                             aspect = "iso", 
                             col = colours)



  return(pcafig)
}


#### function to form a ER-balance subset and derive its median---write it to PAM50 dir
makeCalls.ihc = function(df.cln, seed=118, mat, inputDir=NULL){
  #message("###clinical subtype data.frame should have a column --PatientID-- with which mat cols are also named")
  #message("##IHC subtype column should be named ---IHC---")
  if (is.null(inputDir)) {
    inputDir <- file.path(tempdir(), "Calls.PAM50")
  }

  message("inputDir: ", inputDir)
  # Check and create directory if it doesn't exist
  if (!dir.exists(inputDir)) {
    dir.create(inputDir, recursive = TRUE)
    message("Directory created: ", inputDir)
  } else {
    message("Directory already exists: ", inputDir)
  }

  # Initial checks for 'df.cln' and 'mat'
  if (is.null(df.cln) || !'PatientID' %in% colnames(df.cln) || !'IHC' %in% colnames(df.cln)) {
    stop("Clinical data 'df.cln' is missing or does not contain required columns 'PatientID' and 'IHC'. Refer to the vignette to create 'test.clinical' data.")
  }

  if (is.null(mat) || !is.matrix(mat)) {
    stop("Gene expression matrix 'mat' is missing or not correctly formatted. Refer to the vignette to create 'test.matrix' data.")
  }

   # Continue with function if checks pass
  message("All necessary input data are correctly formatted and present. Grab a cup of coffee while we do the heavy lifting for you!")
  
    # Save the current working directory and ensure it is reset upon function exit
  oldwd <- getwd()
  on.exit(setwd(oldwd))
  
  # Save the current par settings
  oldpar <- par(no.readonly = TRUE)
  
  # Ensure the par settings are reset when the function exits
  on.exit(par(oldpar))
  


  #--Prepare the required input parameters for PAM50 method
  predFiles              = file.path(inputDir, "TEST_pam_mat.txt")#provided inputFile for predictions
  short                  = "Ihc.Mdns"# short name that will be used for output files
  calibrationParameters  = "Mdns.PAM50"# suffix for naming intrinsic subtype column in o/p
  Out.mdns.fl            = "Mdns.PAM50.txt"## where new medians will be written
  calibrationFile        = file.path(inputDir, Out.mdns.fl)
  hasClinical            = F #--PCAPAM50 donot use clinical info module of PAM50 classifier

  # Write the Test.matrix to the file
  write.table(mat, predFiles, sep = "\t", col.names = NA)

  # Verify file creation
  if (!file.exists(predFiles)) {
    stop("File not created: ", predFiles)
  }

  # Define parameters
  paramDir <- system.file( "PAM50/bioclassifier_R", package = "PCAPAM50")#"extdata",
    if (!file.exists(paramDir)) {
    stop("Parameter directory not found: ", paramDir)
  }


  # Print paths for debugging
  message("paramater Dir path: ", paramDir, "\n")
  message("Input File path: ", predFiles, "\n")
  message("Output Median File path: ", calibrationFile, "\n")


  # Convert IHC column to uppercase to handle case insensitivity
  df.cln$IHC = toupper(df.cln$IHC)

  # Get ER- samples data.frame by excluding those starting with "L"
  ERN.ihc = df.cln[!grepl("^L", df.cln$IHC), ]
  #ERN.ihc = df.cln[which(df.cln$IHC %in% c("TN","Her2+")),] ### get ER- samples data.frame
  dim(ERN.ihc)	#[1] 153   9
  
  # Get ER+ samples data.frame (subtypes starting with "L")
  ERP.ihc <- df.cln[grepl("^L", df.cln$IHC), ]
  #ERP.ihc = df.cln[which(df.cln$IHC %in% c("LA","LB1","LB2")),]
  dim(ERP.ihc) #[1] 559   9
  
  # Determine the smaller size between ER+ and ER-
  sample_size <- min(dim(ERN.ihc)[1], dim(ERP.ihc)[1])

  # Sample equal number of ER+ and ER- samples
  set.seed(seed)  # Set seed for reproducibility
  if (dim(ERN.ihc)[1] <= dim(ERP.ihc)[1]) {
    i <- sample(dim(ERP.ihc)[1], sample_size)
    ERP.ihc_sampled <- ERP.ihc[i, ]
    ERN.ihc_sampled <- ERN.ihc
  } else {
    i <- sample(dim(ERN.ihc)[1], sample_size)
    ERN.ihc_sampled <- ERN.ihc[i, ]
    ERP.ihc_sampled <- ERP.ihc
  }
  
  # Subset PAM50 matrix with the IDs corresponding to balanced ER+ and ER-
  mbal.ihc <- mat[, c(ERP.ihc_sampled$PatientID, ERN.ihc_sampled$PatientID)]
  

  #set.seed(seed);i = sample(dim(ERP.ihc)[1],dim(ERN.ihc)[1]) # take equal number of ER+ and ER- samples
  #length(ERP.ihc$PatientID[i]) # ER positive samples
  #length(ERN.ihc$PatientID)    # ER negative samples

  ######=== subset pam50 matrix with the IDs corresponding to balanced ER+ and ER-
  #mbal.ihc = mat[,c(ERP.ihc$PatientID[i],ERN.ihc$PatientID)]

  #dim(mbal.ihc) #[1]  50 306

  # Find median
  mdns      = apply(mbal.ihc,1,median,na.rm=T) # compute median of each row i.e gene
  mdns.df   = as.data.frame(mdns)
  df.mdns   = data.frame(X=rownames(mdns.df),mdns.ihc=mdns.df$mdns) # ER-blanced set based on IHC status alone---
  colnames(df.mdns) = c("X",calibrationParameters)


  # merge this median to PAM50-medians file "mediansPerDataset_v2.txt"
  fl.nm     = paste(paramDir,"mediansPerDataset_v2.txt",sep="/")
  message("Reading file: ", fl.nm)
  if (!file.exists(fl.nm)) {
    stop("File not found: ", fl.nm)
  }



  fl.mdn    = read.table(fl.nm,sep="\t",header=T,stringsAsFactors=F)

  df.al     = merge(fl.mdn,df.mdns,by="X")

  #save(df.al,file="./Rdat/mediansPerDataset_15Nov17.Rdat")
  write.table(df.al, file=calibrationFile, sep="\t", col.names=T,row.names=F,quote=F)

  ####
  # run the assignment algorithm
  # all the global variable required for assignment algorithm are set outside of function where it is been called
  ####

  message("Sourcing: ", paste(paramDir, "subtypePrediction_distributed_PNmodified.R", sep = "/"))
  message("Sourcing: ", paste(paramDir, "subtypePrediction_functions_PNmodified.R", sep = "/"))
  source(paste(paramDir,"subtypePrediction_functions_PNmodified.R",sep="/"))
  source(paste(paramDir,"subtypePrediction_distributed_PNmodified.R",sep="/"))

  # Call the modified function
  res = subtypePrediction_distributed(paramDir, inputDir, short, predFiles, hasClinical, calibrationFile, calibrationParameters)

  # Verify the expected output file paths
  ot.fl = paste(inputDir, "/", short, "_pam50scores.txt", sep = "")
  message("Expected output file: ", ot.fl)

  kl        = read.table(ot.fl,sep="\t",header=T,stringsAsFactors=F)#note md -- for median centered

  df.kl     = data.frame(PatientID=kl$X,Int.SBS.ihcMd=kl$Call,stringsAsFactors=FALSE)
  colnames(df.kl) = c("PatientID",paste("Int.SBS",calibrationParameters,sep="."))
  df.kl     = merge(df.kl,df.cln,by="PatientID")

  return(list(Int.sbs=df.kl, score.fl=kl, mdns.fl=df.al, SBS.colr = res$subtypeColors, outList=res$out))
}


#### function to form a ER-balance subet and derive its median---write it to PAM50 dir
makeCalls.PC1ihc = function(df.cln, seed=118, mat, inputDir=NULL){

  #message("###clinical subtype data.frame should have a column --PatientID-- with which mat cols are also named")
  #message("##IHC subtype column should be named ---IHC---")
  # Set the default inputDir to a temporary directory if not provided
  if (is.null(inputDir)) {
    inputDir <- file.path(tempdir(), "Calls.PC1ihc")
  }

  
  message("inputDir: ", inputDir)

  # Check and create directory if it doesn't exist
  if (!dir.exists(inputDir)) {
    dir.create(inputDir, recursive = TRUE)
    message("Directory created: ", inputDir)
  } else {
    message("Directory already exists: ", inputDir)
  }

  # Initial checks for 'df.cln' and 'mat'
  if (is.null(df.cln) || !'PatientID' %in% colnames(df.cln) || !'IHC' %in% colnames(df.cln)) {
    stop("Clinical data 'df.cln' is missing or does not contain required columns 'PatientID' and 'IHC'. Refer to the vignette to create 'test.clinical' data.")
  }

  if (is.null(mat) || !is.matrix(mat)) {
    stop("Gene expression matrix 'mat' is missing or not correctly formatted. Refer to the vignette to create 'test.matrix' data.")
  }

   # Continue with function if checks pass
  message("All necessary input data are correctly formatted and present. Grab a cup of coffee while we do the heavy lifting for you!")
  
    # Save the current working directory and ensure it is reset upon function exit
  oldwd <- getwd()
  on.exit(setwd(oldwd))
  
  # Save the current par settings
  oldpar <- par(no.readonly = TRUE)
  
  # Ensure the par settings are reset when the function exits
  on.exit(par(oldpar))
  


  #--Prepare the required input parameters for PAM50 method
  predFiles              = file.path(inputDir, "TEST_pam_mat.txt")#provided inputFile for predictions
  short                  = "PC1ihc.Mdns"# short name that will be used for output files
  calibrationParameters  = "Mdns.PC1ihcPAM50"# suffix for naming intrinsic subtype column in o/p
  Out.mdns.fl            = "Mdns.PC1ihcPAM50.txt"## where new medians will be written
  calibrationFile        = file.path(inputDir, Out.mdns.fl)
  hasClinical            = F #--PCAPAM50 donot use clinical info module of PAM50 classifier

  # Write the Test.matrix to the file
  write.table(mat, predFiles, sep = "\t", col.names = NA)

  # Verify file creation
  if (!file.exists(predFiles)) {
    stop("File not created: ", predFiles)
  }

  # Define parameters
  paramDir <- system.file( "PAM50/bioclassifier_R", package = "PCAPAM50")#"extdata",
    if (!file.exists(paramDir)) {
    stop("Parameter directory not found: ", paramDir)
  }


  # Print paths for debugging
  message("paramater Dir path: ", paramDir, "\n")
  message("Input File path: ", predFiles, "\n")
  message("Output Median File path: ", calibrationFile, "\n")


  # Pull the PCA components


  #rv      = rowVars(mat)
  #select  = order(rv, decreasing = TRUE)[seq_len(dim(mat)[1])] # the input is PAM50 matrix --50 genes -- get from dimension
  pca     = prcomp(t(mat))#[select,]
  pc12    = pca$x[,1:2] #get two principal
  df.pc1  = data.frame(PatientID=rownames(pc12),PC1 = pc12[,1],stringsAsFactors=F)
  df.pca1 = merge(df.cln,df.pc1,by="PatientID")
  
  #--our function works best if majority of ER- cases fall in the positive PC1 axis--check
  # Identify ER-negative cases
  er_negative <- !grepl("^L", df.pca1$IHC)
  
  # Determine if the majority of ER-negative cases fall in the negative axis of PC1
  if (sum(df.pca1$PC1[er_negative] < 0) > sum(df.pca1$PC1[er_negative] >= 0)) {
    #print("yes")
    df.pca1$PC1 <- -df.pca1$PC1
  }
  

  # Ensure that IHC is not a factor or has all necessary levels defined
    if (is.factor(df.pca1$IHC)) {
        df.pca1$IHC = as.character(df.pca1$IHC)
    }
    
    # Convert IHC column to uppercase to handle case insensitivity
    df.pca1$IHC <- toupper(df.pca1$IHC)
    
    # Function to count the number of misclassified cases on a given PC1 point ---find the cutoff
    getno = function(x) {
    	p.rgt = length(which(grepl("^L", df.pca1$IHC) & df.pca1$PC1 > x)) / length(which(grepl("^L", df.pca1$IHC)))
    	n.lft = length(which(!grepl("^L", df.pca1$IHC) & df.pca1$PC1 < x)) / length(which(!grepl("^L", df.pca1$IHC)))
    	tot = (p.rgt + n.lft) * 100
    	return(list(PC1 = x, Mis = tot))
    }

  ### function to count the number of mis-classified cases on a given PC1 point ---find the cutoff
  #getno = function(x){
  #  p.rgt = length(which(df.pca1$IHC %in% c("LB1","LB2","LA") & df.pca1$PC1 >x))/length(which(df.pca1$IHC %in% c("LB1","LB2","LA") ))
  #  n.lft = length(which(df.pca1$IHC %in% c("Her2+","TN") & df.pca1$PC1 <x))/ length(which(df.pca1$IHC %in% c("Her2+","TN")))
  #  tot   = (p.rgt + n.lft) * 100
  #  return(list(PC1=x,Mis=tot))
  #}



  df.mis  = do.call(rbind.data.frame,lapply(seq(-20,20,by=0.1),getno))
  num.min = df.mis$PC1[which(df.mis$Mis == min(df.mis$Mis))]
  
  plt.fl     = paste(inputDir,"PC1_misclassified_cases.png",sep="/")
  png(filename = plt.fl, width = 1400, height = 1500, units = "px", pointsize = 10, bg = "white", type = c("cairo"),res=200)
  par(mar=c(5,5,1,1))
  plot(x=df.mis$PC1, y=df.mis$Mis,xlab="PC1",ylab="% of misclassified cases",type="l",lwd=2,col="red",lty=1,cex.axis=2,cex.lab=2) 
  abline(v=num.min, lty=3,lwd=3)
  dev.off()


  #ERP.pc1ihc = df.pca1[which(df.pca1$IHC %in% c("LB1","LB2","LA") & df.pca1$PC1 <= mean(num.min)),] # used mean to overcome situation where there are two minimum
  #ERN.pc1ihc = df.pca1[which(df.pca1$IHC %in% c("Her2+","TN") & df.pca1$PC1 > mean(num.min)),]
  
  ERP.pc1ihc = df.pca1[which(grepl("^L", df.pca1$IHC) & df.pca1$PC1 <= mean(num.min)), ] # used mean to overcome situation where there are two minimum
  ERN.pc1ihc = df.pca1[which(!grepl("^L", df.pca1$IHC) & df.pca1$PC1 > mean(num.min)), ]

  dim(ERP.pc1ihc)#  543   3
  dim(ERN.pc1ihc)#  118   3
  
  # Determine the smaller size between ER+ and ER-
  sample_size <- min(dim(ERP.pc1ihc)[1], dim(ERN.pc1ihc)[1])

  # Sample equal number of ER+ and ER- samples
  set.seed(seed)  # Set seed for reproducibility
  if (dim(ERN.pc1ihc)[1] <= dim(ERP.pc1ihc)[1]) {
    i <- sample(dim(ERP.pc1ihc)[1], sample_size)
    ERP.pc1ihc_sampled <- ERP.pc1ihc[i, ]
    ERN.pc1ihc_sampled <- ERN.pc1ihc
  } else {
    i <- sample(dim(ERN.pc1ihc)[1], sample_size)
    ERN.pc1ihc_sampled <- ERN.pc1ihc[i, ]
    ERP.pc1ihc_sampled <- ERP.pc1ihc
  }
  
  # Subset PAM50 matrix with the IDs corresponding to balanced ER+ and ER-
  mbal.pc1ihc <- mat[, c(ERP.pc1ihc_sampled$PatientID, ERN.pc1ihc_sampled$PatientID)]

  #set.seed(seed);i = sample(dim(ERP.pc1ihc)[1],dim(ERN.pc1ihc)[1]) # take equal number of ER+ and ER- samples
  #length(ERP.pc1ihc$PatientID[i]) # ER positive samples
  #length(ERN.pc1ihc$PatientID)    # ER negative samples

  ######=== subset pam50 matrix with the IDs corresponding to balanced ER+ and ER-
  #mbal.pc1ihc = mat[,c(ERP.pc1ihc$PatientID[i],ERN.pc1ihc$PatientID)]

  dim(mbal.pc1ihc) #[1]  50 236

  # Find median
  mdns      = apply(mbal.pc1ihc,1,median,na.rm=T) # compute median of each row i.e gene
  mdns.df   = as.data.frame(mdns)
  df.mdns   = data.frame(X=rownames(mdns.df),mdns.pc1ihc=mdns.df$mdns) # ER-blanced set based on IHC status alone---
  colnames(df.mdns) = c("X",calibrationParameters)


  # merge this median to PAM50-medians file "mediansPerDataset_v2.txt"
  fl.nm     = paste(paramDir,"mediansPerDataset_v2.txt",sep="/")
  message("Reading file: ", fl.nm)
  if (!file.exists(fl.nm)) {
    stop("File not found: ", fl.nm)
  }


  fl.mdn    = read.table(fl.nm,sep="\t",header=T,stringsAsFactors=F)

  df.al     = merge(fl.mdn,df.mdns,by="X")

  #save(df.al,file="./Rdat/mediansPerDataset_15Nov17.Rdat")
  write.table(df.al, file=calibrationFile, sep="\t", col.names=T,row.names=F,quote=F)

  ####
  # run the assignment algorithm
  ####

  message("Sourcing: ", paste(paramDir, "subtypePrediction_distributed_PNmodified.R", sep = "/"))
  message("Sourcing: ", paste(paramDir, "subtypePrediction_functions_PNmodified.R", sep = "/"))
  source(paste(paramDir,"subtypePrediction_functions_PNmodified.R",sep="/"))
  source(paste(paramDir,"subtypePrediction_distributed_PNmodified.R",sep="/"))

  # Call the modified function
  res = subtypePrediction_distributed(paramDir, inputDir, short, predFiles, hasClinical, calibrationFile, calibrationParameters)

  # Verify the expected output file paths
  ot.fl = paste(inputDir, "/", short, "_pam50scores.txt", sep = "")
  message("Expected output file: ", ot.fl)
  kl        = read.table(ot.fl,sep="\t",header=T,stringsAsFactors=F)#note md -- for median centered

  df.kl     = data.frame(PatientID=kl$X,Int.SBS.ihcMd=kl$Call,stringsAsFactors=FALSE)
  colnames(df.kl) = c("PatientID",paste("Int.SBS",calibrationParameters,sep="."))
  df.kl     = merge(df.kl,df.cln,by="PatientID")

  return(list(Int.sbs=df.kl, score.fl=kl, mdns.fl=df.al, SBS.colr = res$subtypeColors, outList=res$out, PC1cutoff=df.mis, DF.PC1=df.pca1))
}


makeCalls.v1PAM = function(df.pam, seed=118, mat, inputDir=NULL){

  #message("###df.pam data.frame should have a column --PatientID-- with which mat cols are also named")
  #message("##v1PAM subtype column should be named ---PAM50---")
  if (is.null(inputDir)) {
    inputDir <- file.path(tempdir(), "Calls.PCAPAM50")
  }

  message("inputDir: ", inputDir)

  # Check and create directory if it doesn't exist
  if (!dir.exists(inputDir)) {
    dir.create(inputDir, recursive = TRUE)
    message("Directory created: ", inputDir)
  } else {
    message("Directory already exists: ", inputDir)
  }

  # Initial checks for 'df.cln' and 'mat'
  if (is.null(df.pam) || !'PatientID' %in% colnames(df.pam) || !'PAM50' %in% colnames(df.pam)) {
    stop("Clinical data 'df.pam' is missing or does not contain required columns 'PatientID' and 'PAM50'. Refer to the vignette to create 'test.intermediate subtype' data.")
  }

  if (is.null(mat) || !is.matrix(mat)) {
    stop("Gene expression matrix 'mat' is missing or not correctly formatted. Refer to the vignette to create 'test.matrix' data.")
  }

   # Continue with function if checks pass
  message("All necessary input data are correctly formatted and present. Grab a cup of coffee while we do the heavy lifting for you!")
  
    # Save the current working directory and ensure it is reset upon function exit
  oldwd <- getwd()
  on.exit(setwd(oldwd))
  
  # Save the current par settings
  oldpar <- par(no.readonly = TRUE)
  
  # Ensure the par settings are reset when the function exits
  on.exit(par(oldpar))
  


  #--Prepare the required input parameters for PAM50 method
  predFiles              = file.path(inputDir, "TEST_pam_mat.txt")#provided inputFile for predictions
  short                  = "PCAPAM50.Mdns"# short name that will be used for output files
  calibrationParameters  = "Mdns.PCAPAM50"# suffix for naming intrinsic subtype column in o/p
  Out.mdns.fl            = "Mdns.PCAPAM50.txt"## where new medians will be written
  calibrationFile        = file.path(inputDir, Out.mdns.fl)
  hasClinical            = F #--PCAPAM50 donot use clinical info module of PAM50 classifier

  # Write the Test.matrix to the file
  write.table(mat, predFiles, sep = "\t", col.names = NA)

  # Verify file creation
  if (!file.exists(predFiles)) {
    stop("File not created: ", predFiles)
  }

  # Define parameters
  paramDir <- system.file( "PAM50/bioclassifier_R", package = "PCAPAM50")#"extdata",
    if (!file.exists(paramDir)) {
    stop("Parameter directory not found: ", paramDir)
  }


  # Print paths for debugging
  message("paramater Dir path: ", paramDir, "\n")
  message("Input File path: ", predFiles, "\n")
  message("Output Median File path: ", calibrationFile, "\n")


  ERN.pam = df.pam[which(df.pam$PAM50 %in% c("Basal")),] ### get ER- samples data.frame
  dim(ERN.pam)	#[1] 119   2

  ERP.pam = df.pam[which(df.pam$PAM50 %in% c("LumA")),]
  dim(ERP.pam) #[1] 352   2
  
  # Determine the smaller size between ER+ and ER-
  sample_size <- min(dim(ERP.pam)[1], dim(ERN.pam)[1])

  # Sample equal number of ER+ and ER- samples
  set.seed(seed)  # Set seed for reproducibility
  if (dim(ERN.pam)[1] <= dim(ERP.pam)[1]) {
    i <- sample(dim(ERP.pam)[1], sample_size)
    ERP.pam_sampled <- ERP.pam[i, ]
    ERN.pam_sampled <- ERN.pam
  } else {
    i <- sample(dim(ERN.pam)[1], sample_size)
    ERN.pam_sampled <- ERN.pam[i, ]
    ERP.pam_sampled <- ERP.pam
  }
  
  # Subset PAM50 matrix with the IDs corresponding to balanced ER+ and ER-
  mbal.pam <- mat[, c(ERP.pam_sampled$PatientID, ERN.pam_sampled$PatientID)]

  #set.seed(seed);i = sample(dim(ERP.pam)[1],dim(ERN.pam)[1]) # take equal number of ER+ and ER- samples
  #length(ERP.pam$PatientID[i]) # ER positive samples
  #length(ERN.pam$PatientID)    # ER negative samples

  ######=== subset pam50 matrix with the IDs corresponding to balanced ER+ and ER-
  #mbal.pam = mat[,c(ERP.pam$PatientID[i],ERN.pam$PatientID)]

  dim(mbal.pam) #[1]  50 306

  # Find median
  mdns      = apply(mbal.pam,1,median,na.rm=T) # compute median of each row i.e gene
  mdns.df   = as.data.frame(mdns)
  df.mdns   = data.frame(X=rownames(mdns.df),mdns.pam=mdns.df$mdns) # ER-blanced set based on IHC status alone---
  colnames(df.mdns) = c("X",calibrationParameters)


  # merge this median to PAM50-medians file "mediansPerDataset_v2.txt"
  fl.nm     = paste(paramDir,"mediansPerDataset_v2.txt",sep="/")
  message("Reading file: ", fl.nm)
  if (!file.exists(fl.nm)) {
    stop("File not found: ", fl.nm)
  }



  fl.mdn    = read.table(fl.nm,sep="\t",header=T,stringsAsFactors=F)

  df.al     = merge(fl.mdn,df.mdns,by="X")

  #save(df.al,file="./Rdat/mediansPerDataset_15Nov17.Rdat")
  write.table(df.al, file=calibrationFile, sep="\t", col.names=T,row.names=F,quote=F)#paste(paramDir,mdns.outFile,sep="/")

  ####
  # run the assignment algorithm
  # all the global variable required for assignment algorithm are set outside of function where it is been called
  ####

  message("Sourcing: ", paste(paramDir, "subtypePrediction_distributed_PNmodified.R", sep = "/"))
  message("Sourcing: ", paste(paramDir, "subtypePrediction_functions_PNmodified.R", sep = "/"))
  source(paste(paramDir,"subtypePrediction_functions_PNmodified.R",sep="/"))
  source(paste(paramDir,"subtypePrediction_distributed_PNmodified.R",sep="/"))

  # Call the modified function
  res = subtypePrediction_distributed(paramDir, inputDir, short, predFiles, hasClinical, calibrationFile, calibrationParameters)

  # Verify the expected output file paths
  ot.fl = paste(inputDir, "/", short, "_pam50scores.txt", sep = "")
  message("Expected output file: ", ot.fl)

  kl        = read.table(ot.fl,sep="\t",header=T,stringsAsFactors=F)#note md -- for median centered

  df.kl     = data.frame(PatientID=kl$X,Int.SBS.pamMd=kl$Call,stringsAsFactors=FALSE)
  colnames(df.kl) = c("PatientID",paste("Int.SBS",calibrationParameters,sep="."))
  df.kl     = merge(df.kl,df.pam,by="PatientID") 

  return(list(Int.sbs=df.kl, score.fl=kl, mdns.fl=df.al, SBS.colr = res$subtypeColors, outList=res$out))
}



