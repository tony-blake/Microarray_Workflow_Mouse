makebeta1_pathwaygenetable <- function(indata){
  
  # Get list of genes in bph1
  bp1test <- make_beta1matrix2d(indata)
  
  testgenesinGOBP1 <- bp1test[[2]]
  bp1matrix <- bp1test[[1]]
  
  # create data frame from list
  testGG = lapply(testgenesinGOBP1, ldply)
  testDDBP1 = ldply(testGG)
  
  library(annotate)
  
  #convert all entrez id's in "genesinGO" list to gene symbols
  
  tsymbolsforGOBP1 <- list()
  
  for(i in 1:length(testgenesinGOBP1)){
    
    tsymbolsforGOBP1[[i]] <- mget(testgenesinGOBP1[[i]], org.Mm.egSYMBOL)
    
    tsymbolsforGOBP1
  }
  
  # create dataframe from list of lists
  tSG1BP1 = lapply(tsymbolsforGOBP1, ldply)
  tddbp1 = ldply(tSG1BP1)
  
  #Add symbols column to genes in pathways dataframe
  
  colnames(tddbp1)[2] <- "symbols"
  
  testDDBP1$SYMBOLS <- tddbp1$symbols
  
  colnames(testDDBP1)[1] <- "pathways"
  colnames(testDDBP1)[2] <- "ENTREZID"
  
  tDDbp1 <- testDDBP1[!duplicated(testDDBP1),]
  tDDbp1e <- tDDbp1[,c(1,3)]
  tddbp1e <- aggregate(SYMBOLS ~ pathways, FUN = "c", data=tDDbp1e)
  
  
  
  
  
  elements <- list()
  
  elements[[1]] <- bp1matrix
  elements[[2]] <- tddbp1e
  elements[[3]] <- testDDBP1
  
  elements
  
}

#################################################### beta2 ######################################################

makebeta2_pathwaygenetable <- function(indata){
  
  # Get list of genes in bph1
  bp1test <- make_beta2matrix2d(indata)
  
  testgenesinGOBP1 <- bp1test[[2]]
  bp1matrix <- bp1test[[1]]
  
  # create data frame from list
  testGG = lapply(testgenesinGOBP1, ldply)
  testDDBP1 = ldply(testGG)
  
  library(annotate)
  
  #convert all entrez id's in "genesinGO" list to gene symbols
  
  tsymbolsforGOBP1 <- list()
  
  for(i in 1:length(testgenesinGOBP1)){
    
    tsymbolsforGOBP1[[i]] <- mget(testgenesinGOBP1[[i]], org.Mm.egSYMBOL)
    
    tsymbolsforGOBP1
  }
  
  # create dataframe from list of lists
  tSG1BP1 = lapply(tsymbolsforGOBP1, ldply)
  tddbp1 = ldply(tSG1BP1)
  
  #Add symbols column to genes in pathways dataframe
  
  colnames(tddbp1)[2] <- "symbols"
  
  testDDBP1$SYMBOLS <- tddbp1$symbols
  
  colnames(testDDBP1)[1] <- "pathways"
  colnames(testDDBP1)[2] <- "ENTREZID"
  
  tDDbp1 <- testDDBP1[!duplicated(testDDBP1),]
  tDDbp1e <- tDDbp1[,c(1,3)]
  tddbp1e <- aggregate(SYMBOLS ~ pathways, FUN = "c", data=tDDbp1e)
  
  
  
  
  
  elements <- list()
  
  elements[[1]] <- bp1matrix
  elements[[2]] <- tddbp1e
  elements[[3]] <- testDDBP1
  
  elements
  
}

################################################# beta1 down ##############################################

makebeta1down_pathwaygenetable <- function(indata){
  
  # Get list of genes in bph1
  bp1test <- makedown_beta1matrix2d(indata)
  
  testgenesinGOBP1 <- bp1test[[2]]
  bp1matrix <- bp1test[[1]]
  
  # create data frame from list
  testGG = lapply(testgenesinGOBP1, ldply)
  testDDBP1 = ldply(testGG)
  
  library(annotate)
  
  #convert all entrez id's in "genesinGO" list to gene symbols
  
  tsymbolsforGOBP1 <- list()
  
  for(i in 1:length(testgenesinGOBP1)){
    
    tsymbolsforGOBP1[[i]] <- mget(testgenesinGOBP1[[i]], org.Mm.egSYMBOL)
    
    tsymbolsforGOBP1
  }
  
  # create dataframe from list of lists
  tSG1BP1 = lapply(tsymbolsforGOBP1, ldply)
  tddbp1 = ldply(tSG1BP1)
  
  #Add symbols column to genes in pathways dataframe
  
  colnames(tddbp1)[2] <- "symbols"
  
  testDDBP1$SYMBOLS <- tddbp1$symbols
  
  colnames(testDDBP1)[1] <- "pathways"
  colnames(testDDBP1)[2] <- "ENTREZID"
  
  tDDbp1 <- testDDBP1[!duplicated(testDDBP1),]
  tDDbp1e <- tDDbp1[,c(1,3)]
  tddbp1e <- aggregate(SYMBOLS ~ pathways, FUN = "c", data=tDDbp1e)
  
  
  
  
  
  elements <- list()
  
  elements[[1]] <- bp1matrix
  elements[[2]] <- tddbp1e
  elements[[3]] <- testDDBP1
  
  elements
  
}

#################################################### beta2 ######################################################

makebeta2down_pathwaygenetable <- function(indata){
  
  # Get list of genes in bph1
  bp1test <- makedown_beta2matrix2d(indata)
  
  testgenesinGOBP1 <- bp1test[[2]]
  bp1matrix <- bp1test[[1]]
  
  # create data frame from list
  testGG = lapply(testgenesinGOBP1, ldply)
  testDDBP1 = ldply(testGG)
  
  library(annotate)
  
  #convert all entrez id's in "genesinGO" list to gene symbols
  
  tsymbolsforGOBP1 <- list()
  
  for(i in 1:length(testgenesinGOBP1)){
    
    tsymbolsforGOBP1[[i]] <- mget(testgenesinGOBP1[[i]], org.Mm.egSYMBOL)
    
    tsymbolsforGOBP1
  }
  
  # create dataframe from list of lists
  tSG1BP1 = lapply(tsymbolsforGOBP1, ldply)
  tddbp1 = ldply(tSG1BP1)
  
  #Add symbols column to genes in pathways dataframe
  
  colnames(tddbp1)[2] <- "symbols"
  
  testDDBP1$SYMBOLS <- tddbp1$symbols
  
  colnames(testDDBP1)[1] <- "pathways"
  colnames(testDDBP1)[2] <- "ENTREZID"
  
  tDDbp1 <- testDDBP1[!duplicated(testDDBP1),]
  tDDbp1e <- tDDbp1[,c(1,3)]
  tddbp1e <- aggregate(SYMBOLS ~ pathways, FUN = "c", data=tDDbp1e)
  
  
  
  
  
  elements <- list()
  
  elements[[1]] <- bp1matrix
  elements[[2]] <- tddbp1e
  elements[[3]] <- testDDBP1
  
  elements
  
}

