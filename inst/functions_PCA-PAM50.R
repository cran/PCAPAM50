#--last updated 5/29/2024 3:57 PM

###### modeling after plotPCA of DESeq
##
my.plotPCA = function(x, intgroup, ablne = 0,colours = c("red","hotpink","darkblue","lightblue","red3","hotpink3","royalblue3","lightskyblue3"),LINE.V=T)
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
makeCalls.ihc = function(df.cln, seed=118, mat, inputDir="Calls.pam50"){
  #message("###clinical subtype data.frame should have a column --PatientID-- with which mat cols are also named")
  #message("##IHC subtype column should be named ---IHC---")
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

  #--Prepare the required input parameters for PAM50 method
  predFiles              = file.path(inputDir, "TEST_pam_mat.txt")#provided inputFile for predictions
  short                  = "Ihc.Mdns"# short name that will be used for output files
  calibrationParameters  = "Mdns.Ihc"# suffix for naming intrinsic subtype column in o/p
  Out.mdns.fl            = "Mdns.Ihc.txt"## where new medians will be written
  calibrationFile        = file.path(inputDir, Out.mdns.fl)
  hasClinical            = F #--pcapam50 donot use clinical info module of PAM50 classifier

  # Write the Test.matrix to the file
  write.table(mat, predFiles, sep = "\t", col.names = NA)

  # Verify file creation
  if (!file.exists(predFiles)) {
    stop("File not created: ", predFiles)
  }

  # Define parameters
  paramDir <- system.file( "PAM50/bioclassifier_R", package = "pcapam50")#"extdata",
  if (!file.exists(paramDir)) {
    stop("Parameter directory not found: ", paramDir)
  }


  # Print paths for debugging
  cat("paramater Dir path: ", paramDir, "\n")
  cat("Input File path: ", predFiles, "\n")
  cat("Output Median File path: ", calibrationFile, "\n")




  ERN.ihc = df.cln[which(df.cln$IHC %in% c("TN","Her2+")),] ### get ER- samples data.frame
  dim(ERN.ihc)	#[1] 153   9

  ERP.ihc = df.cln[which(df.cln$IHC %in% c("LA","LB1","LB2")),]
  dim(ERP.ihc) #[1] 559   9

  set.seed(seed);i = sample(dim(ERP.ihc)[1],dim(ERN.ihc)[1]) # take equal number of ER+ and ER- samples
  length(ERP.ihc$PatientID[i]) # ER positive samples
  length(ERN.ihc$PatientID)    # ER negative samples

  ######=== subset pam50 matrix with the IDs corresponding to balanced ER+ and ER-
  mbal.ihc = mat[,c(ERP.ihc$PatientID[i],ERN.ihc$PatientID)]

  dim(mbal.ihc) #[1]  50 306

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
makeCalls.PC1ihc = function(df.cln, seed=118, mat, inputDir="Calls.PC1ihc"){

  #message("###clinical subtype data.frame should have a column --PatientID-- with which mat cols are also named")
  #message("##IHC subtype column should be named ---IHC---")
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

  #--Prepare the required input parameters for PAM50 method
  predFiles              = file.path(inputDir, "TEST_pam_mat.txt")#provided inputFile for predictions
  short                  = "PC1ihc.Mdns"# short name that will be used for output files
  calibrationParameters  = "Mdns.PC1ihc"# suffix for naming intrinsic subtype column in o/p
  Out.mdns.fl            = "Mdns.PC1ihc.txt"## where new medians will be written
  calibrationFile        = file.path(inputDir, Out.mdns.fl)
  hasClinical            = F #--pcapam50 donot use clinical info module of PAM50 classifier

  # Write the Test.matrix to the file
  write.table(mat, predFiles, sep = "\t", col.names = NA)

  # Verify file creation
  if (!file.exists(predFiles)) {
    stop("File not created: ", predFiles)
  }

  # Define parameters
  paramDir <- system.file( "PAM50/bioclassifier_R", package = "pcapam50")#"extdata",
  if (!file.exists(paramDir)) {
    stop("Parameter directory not found: ", paramDir)
  }


  # Print paths for debugging
  cat("paramater Dir path: ", paramDir, "\n")
  cat("Input File path: ", predFiles, "\n")
  cat("Output Median File path: ", calibrationFile, "\n")


  # Pull the PCA components


  #rv      = rowVars(mat)
  #select  = order(rv, decreasing = TRUE)[seq_len(dim(mat)[1])] # the input is PAM50 matrix --50 genes -- get from dimension
  pca     = prcomp(t(mat))#[select,]
  pc12    = pca$x[,1:2] #get two principal
  df.pc1  = data.frame(PatientID=rownames(pc12),PC1 = pc12[,1],stringsAsFactors=F)
  df.pca1 = merge(df.cln,df.pc1,by="PatientID")

  # Ensure that IHC is not a factor or has all necessary levels defined
  if (is.factor(df.pca1$IHC)) {
    df.pca1$IHC = as.character(df.pca1$IHC)
  }

  ### function to count the number of mis-classified cases on a given PC1 point ---find the cutoff
  getno = function(x){
    p.rgt = length(which(df.pca1$IHC %in% c("LB1","LB2","LA") & df.pca1$PC1 >x))/length(which(df.pca1$IHC %in% c("LB1","LB2","LA") ))
    n.lft = length(which(df.pca1$IHC %in% c("Her2+","TN") & df.pca1$PC1 <x))/ length(which(df.pca1$IHC %in% c("Her2+","TN")))
    tot   = (p.rgt + n.lft) * 100
    return(list(PC1=x,Mis=tot))
  }



  df.mis  = do.call(rbind.data.frame,lapply(seq(-20,20,by=0.1),getno))
  num.min = df.mis$PC1[which(df.mis$Mis == min(df.mis$Mis))]

  ERP.pc1ihc = df.pca1[which(df.pca1$IHC %in% c("LB1","LB2","LA") & df.pca1$PC1 <= mean(num.min)),] # used mean to overcome situation where there are two minimum
  ERN.pc1ihc = df.pca1[which(df.pca1$IHC %in% c("Her2+","TN") & df.pca1$PC1 > mean(num.min)),]

  dim(ERP.pc1ihc)#  543   3
  dim(ERN.pc1ihc)#  118   3

  set.seed(seed);i = sample(dim(ERP.pc1ihc)[1],dim(ERN.pc1ihc)[1]) # take equal number of ER+ and ER- samples
  length(ERP.pc1ihc$PatientID[i]) # ER positive samples
  length(ERN.pc1ihc$PatientID)    # ER negative samples

  ######=== subset pam50 matrix with the IDs corresponding to balanced ER+ and ER-
  mbal.pc1ihc = mat[,c(ERP.pc1ihc$PatientID[i],ERN.pc1ihc$PatientID)]

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


makeCalls.v1PAM = function(df.pam, seed=118, mat, inputDir="Calls.pcapam50"){

  #message("###df.pam data.frame should have a column --PatientID-- with which mat cols are also named")
  #message("##v1PAM subtype column should be named ---PAM50---")
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

  #--Prepare the required input parameters for PAM50 method
  predFiles              = file.path(inputDir, "TEST_pam_mat.txt")#provided inputFile for predictions
  short                  = "pcapam50.Mdns"# short name that will be used for output files
  calibrationParameters  = "Mdns.pcapam50"# suffix for naming intrinsic subtype column in o/p
  Out.mdns.fl            = "Mdns.pcapam50.txt"## where new medians will be written
  calibrationFile        = file.path(inputDir, Out.mdns.fl)
  hasClinical            = F #--pcapam50 donot use clinical info module of PAM50 classifier

  # Write the Test.matrix to the file
  write.table(mat, predFiles, sep = "\t", col.names = NA)

  # Verify file creation
  if (!file.exists(predFiles)) {
    stop("File not created: ", predFiles)
  }

  # Define parameters
  paramDir <- system.file( "PAM50/bioclassifier_R", package = "pcapam50")#"extdata",
  if (!file.exists(paramDir)) {
    stop("Parameter directory not found: ", paramDir)
  }


  # Print paths for debugging
  cat("paramater Dir path: ", paramDir, "\n")
  cat("Input File path: ", predFiles, "\n")
  cat("Output Median File path: ", calibrationFile, "\n")


  ERN.pam = df.pam[which(df.pam$PAM50 %in% c("Basal")),] ### get ER- samples data.frame
  dim(ERN.pam)	#[1] 119   2

  ERP.pam = df.pam[which(df.pam$PAM50 %in% c("LumA")),]
  dim(ERP.pam) #[1] 352   2

  set.seed(seed);i = sample(dim(ERP.pam)[1],dim(ERN.pam)[1]) # take equal number of ER+ and ER- samples
  length(ERP.pam$PatientID[i]) # ER positive samples
  length(ERN.pam$PatientID)    # ER negative samples

  ######=== subset pam50 matrix with the IDs corresponding to balanced ER+ and ER-
  mbal.pam = mat[,c(ERP.pam$PatientID[i],ERN.pam$PatientID)]

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
  df.kl     = merge(df.kl,df.pam[,-which(names(df.pam) %in% c("PAM50"))],by="PatientID") # avoiding PAM50 column

  return(list(Int.sbs=df.kl, score.fl=kl, mdns.fl=df.al, SBS.colr = res$subtypeColors, outList=res$out))
}
