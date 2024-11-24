subtypePrediction_distributed <- function(paramDir, inputDir, short, predFiles, hasClinical, calibrationFile, calibrationParameters) {
  # Download normalized expression data
  
  # Load the published centroids for classification
  trainCentroids <- paste(paramDir, "pam50_centroids.txt", sep = "/") # _PNmdfd
  trainFile <- paste(paramDir, "220arrays_nonUBCcommon+12normal_50g.txt", sep = "/")
  proliferationGenes <- c("CCNB1", "UBE2C", "BIRC5", "KNTC2", "CDC20", "PTTG1", "RRM2", "MKI67", "TYMS", "CEP55", "CDCA1")
  stdArray <- TRUE # just for visualization, and only set to F if many missing genes
  # predFiles <- paste(inputDir, inputFile, sep = "/")
  
  # Some constants
  glthreshold <- -0.15
  ghthreshold <- 0.1
  gplthreshold <- -0.25
  gphthreshold <- 0.1
  clthreshold <- -0.1
  chthreshold <- 0.2
  cplthreshold <- -0.2
  cphthreshold <- 0.2
  
  # Begin analyses
  
  # Only need train data for visualizations
  x <- readarray(trainFile, hr = 2)
  x$xd <- standardize(medianCtr(x$xd))
  
  # Load the published centroids for classification
  pamout.centroids <- read.table(trainCentroids, sep = "\t", header = TRUE, row.names = 1)
  
  pdfname1 <- paste(inputDir, paste("predictionScores_pam50RankCorrelation_1_", short, ".pdf", sep = ""), sep = "/")
  pdfname2 <- paste(inputDir, paste("predictionScores_pam50RankCorrelation_2_", short, ".pdf", sep = ""), sep = "/")
  clustername <- paste(inputDir, paste(short, "_PAM50_normalized_heatmap", sep = ""), sep = "/")
  outFile <- paste(inputDir, paste(short, "_pam50scores.txt", sep = ""), sep = "/")
  
  # Read in the data file
  if (hasClinical) {
    xhr <- 2
  } else {
    xhr <- 1
  }
  y <- readarray(predFiles, hr = xhr, method = collapseMethod,impute=FALSE)
  
  # Normalization
  if (is.na(calibrationParameters)) {
    y$xd <- medianCtr(y$xd)
  } else {
    if (calibrationParameters != -1) {
      medians <- readarray(calibrationFile, hr = 1)
      print(paste("calibration to:", calibrationParameters)) # dimnames(medians$xd)[[2]][calibrationParameters]
      tm <- overlapSets(medians$xd, y$xd)
      y$xd <- (tm$y - tm$x[, calibrationParameters])
    }
  }
  
  num.missing <- NA
  
  if (stdArray) {
    y$xd <- standardize(y$xd)
  }
  
  erScore <- as.vector(t(y$xd["ESR1", ]))
  her2Score <- as.vector(t(y$xd["ERBB2", ]))
  
  # Assign the subtype scores and calculate the proliferation score
  this.proliferationGenes <- dimnames(y$xd)[[1]] %in% proliferationGenes
  
  prolifScore <- apply(y$xd[this.proliferationGenes, ], 2, mean, na.rm = TRUE)
  
  out <- sspPredict(pamout.centroids, classes = "", y$xd, std = FALSE, distm = "spearman", centroids = TRUE)
  out$distances <- -1 * out$distances
  
  call.conf <- c()
  for (j in 1:length(out$predictions)) {
    call.conf[j] <- 1 - cor.test(out$testData[, j], out$centroids[, which(colnames(pamout.centroids) == out$predictions[j])], method = "spearman")$p.value
  }
  call.conf <- round(call.conf, 2)
  
  # Calculate the risk scores
  genomic <- 0.04210193 * out$distances[, 1] + 0.12466938 * out$distances[, 2] + -0.35235561 * out$distances[, 3] + 0.14213283 * out$distances[, 4]
  genomicWprolif <- -0.0009299747 * out$distances[, 1] + 0.0692289192 * out$distances[, 2] + -0.0951505484 * out$distances[, 3] + 0.0493487685 * out$distances[, 4] + 0.3385116381 * prolifScore
  if (hasClinical) {
    xT <- as.numeric(as.vector(y$classes$T))
    combined <- 0.0442770 * out$distances[, 1] + 0.1170297 * out$distances[, 2] + -0.2608388 * out$distances[, 3] + 0.1055908 * out$distances[, 4] + 0.1813751 * xT
    combinedWprolif <- -0.009383416 * out$distances[, 1] + 0.073725503 * out$distances[, 2] + -0.090436516 * out$distances[, 3] + 0.053013865 * out$distances[, 4] + 0.131605960 * xT + 0.327259375 * prolifScore
  }
  
  # Threshold the risk score
  griskgroups <- genomic
  griskgroups[genomic > ghthreshold] <- "high"
  griskgroups[genomic > glthreshold & genomic < ghthreshold] <- "med"
  griskgroups[genomic < glthreshold] <- "low"
  gpriskgroups <- genomicWprolif
  gpriskgroups[genomicWprolif > gphthreshold] <- "high"
  gpriskgroups[genomicWprolif > gplthreshold & genomicWprolif < gphthreshold] <- "med"
  gpriskgroups[genomicWprolif < gplthreshold] <- "low"
  
  genomic <- 100 * (genomic + 0.35) / 0.85
  genomicWprolif <- 100 * (genomicWprolif + 0.35) / 0.85
  
  # Write output files
  if (hasClinical) {
    criskgroups <- combined
    criskgroups[combined > chthreshold] <- "high"
    criskgroups[combined > clthreshold & combined < chthreshold] <- "med"
    criskgroups[combined < clthreshold] <- "low"
    cpriskgroups <- combinedWprolif
    cpriskgroups[combinedWprolif > cphthreshold] <- "high"
    cpriskgroups[combinedWprolif > cplthreshold & combinedWprolif < cphthreshold] <- "med"
    cpriskgroups[combinedWprolif < cplthreshold] <- "low"
    
    combined <- 100 * (combined + 0.35) / 0.85
    combinedWprolif <- 100 * (combinedWprolif + 0.35) / 0.85
    
    outtable <- cbind(out$distances, out$predictions, call.conf, genomic, griskgroups, prolifScore, genomicWprolif, gpriskgroups, combined, criskgroups, combinedWprolif, cpriskgroups, erScore, her2Score)
    dimnames(outtable)[[2]] <- c("Basal", "Her2", "LumA", "LumB", "Normal", "Call", "Confidence",
                                 "ROR-S (Subtype Only)", "ROR-S Group (Subtype Only)", "Proliferation Score",
                                 "ROR-P (Subtype + Proliferation)", "ROR-P Group (Subtype + Proliferation)",
                                 "ROR-C (Subtype + Clinical)", "ROR-C Group (Subtype + Clinical)",
                                 "ROR-PC (Subtype + Clinical + Proliferation)", "ROR-PC Group (Subtype + Clinical + Proliferation)",
                                 "ER", "Her2")
  } else {
    outtable <- cbind(out$distances, out$predictions, call.conf, genomic, griskgroups, prolifScore, genomicWprolif, gpriskgroups, erScore, her2Score)
    dimnames(outtable)[[2]] <- c("Basal", "Her2", "LumA", "LumB", "Normal", "Call", "Confidence",
                                 "ROR-S (Subtype Only)", "ROR-S Group (Subtype Only)", "Proliferation Score",
                                 "ROR-P (Subtype + Proliferation)", "ROR-P Group (Subtype + Proliferation)",
                                 "ER", "Her2")
  }
  write.table(outtable, outFile, sep = "\t", col.names = NA)
  
  # Make some plots for evaluation
  print(paste("ER range:", quantile(erScore, .9, na.rm = TRUE) - quantile(erScore, .1, na.rm = TRUE)))
  
  subtypeColors <- as.character(out$predictions)
  subtypeColors[subtypeColors == "Basal"] <- "red"
  subtypeColors[subtypeColors == "Her2"] <- "hotpink"
  subtypeColors[subtypeColors == "LumA"] <- "darkblue"
  subtypeColors[subtypeColors == "LumB"] <- "skyblue"
  subtypeColors[subtypeColors == "Normal"] <- "green"
  conf.colors <- as.character(call.conf)
  conf.colors[call.conf >= 0.95] <- "black"
  conf.colors[call.conf < 0.95] <- "red"
  
  df.sub = as.data.frame(out$predictions)
  sub.nm = paste("Subtype",gsub('.*\\.','',calibrationParameters),sep=".")
  colnames(df.sub) = sub.nm
  df.sub$Confidence = call.conf
  
  pdf(paste(clustername, ".pdf", sep = ""), 
    width = 8.5,#max(8.5, ncol(out$testData) * 0.13),  # Adjust width based on the number of columns
    height = 11)#max(11, nrow(out$testData) * 0.7)) # Adjust height based on the number of rows
  ht = myComplexHeatmap(x= out$testData, t.colors = df.sub, contrast = 2) 
  draw(ht,heatmap_legend_side="left", annotation_legend_side="top")
  dev.off()
  
  pdf(pdfname1, height = 10, width = 12)
  pars <- par(no.readonly = TRUE)
  myplot(out, short, prolifScore)
  dev.off()
  
  tm<-overlapSets(x$xd,y$xd)
  tm$x<-tm$x[,!is.na(x$classes$subtype)]
  tm<-cbind(tm$x,	impute::impute.knn(as.matrix(tm$y))$data)
  classes<-matrix(nrow=4,ncol=dim(tm)[2])
  nTrainSamples<-length(x$classes$subtype[!is.na(x$classes$subtype)])
  
  classes[1, ] <- c(rep("train", nTrainSamples), rep(short, dim(y$xd)[2]))
  classes[2, ] <- c(x$classes$subtype[!is.na(x$classes$subtype)], rep(NA, dim(tm)[2] - dim(x$xd[, !is.na(x$classes$subtype)])[2]))
  classes[3, ] <- c(rep(NA, dim(tm)[2] - length(out$predictions)), out$predictions)
  
  pdf(pdfname2, height = 6, width = 12)
  tm <- scale(tm, center = FALSE)
  par(mfrow = c(1, 3))
  pcaEA(tm, classes[1, ], mainStr = "Training and Test sets", showNames = FALSE, showClasses = FALSE)
  pcaEA(tm[, !is.na(as.vector(t(classes[2, ])))], classes[2, !is.na(as.vector(t(classes[2, ])))], mainStr = "Training cases", showNames = FALSE, showClasses = FALSE, groupColors = c("red", "hotpink", "darkblue", "skyblue", "green"))
  pcaEA(tm[, !is.na(as.vector(t(classes[3, ])))], classes[3, !is.na(as.vector(t(classes[3, ])))], mainStr = "Test cases", showNames = FALSE, showClasses = FALSE, groupColors = c("red", "hotpink", "darkblue", "skyblue", "green"))
  par(pars)
  dev.off()
  
  return(list(subtypeColors = subtypeColors, out = out))
}
