make_beta1matrix2d <- function(indata){
  
  dim(indata <- indata[indata[,3] <= 0.05,])
  library(org.Mm.eg.db); library(gplots); library(bioDist)
  # get gene in GO cat
  gGOcat <- list()
  for(i in 1:dim(indata)[1]){
    gGOcat[[indata[,1][i]]] <- as.vector(base::unlist(mget(indata[,1][i], org.Mm.egGO2ALLEGS)))
  }
  
  # filter each GO cat by the genes found to be DE
  f <- function(vecGO, vecDE){
    return(vecGO[which(vecGO %in% vecDE)])
  }
  
  gGOcat2 <- list()
  
  gGOcat <- lapply(gGOcat, FUN=function(x){f(x,GOstats_upgenes[[1]])})
  ar <- array()
  ar <- array(ar, c(length(gGOcat),length(gGOcat)), dimnames=list(names(gGOcat),names(gGOcat)))
  for(i in 1:length(gGOcat)) {
    for(j in 1:length(gGOcat)) {
      
      # intersect / union
      ar[i,j] <- (length(Reduce(intersect, list(gGOcat[[i]], gGOcat[[j]]))) / length(Reduce(union, list(gGOcat[[i]], gGOcat[[j]]))))
      
    }
  }
  dimnames(ar)[[1]] <- names(gGOcat)
  dimnames(ar)[[2]] <- names(gGOcat)
  
  
  goNames <- as.vector(base::unlist(mget(names(gGOcat), GOTERM)))
  goNames <- as.vector(base::unlist(lapply(goNames, FUN=function(x){x@Term})))
  
  goNames <- paste(goNames, names(gGOcat), sep=" - ")
  dimnames(ar)[[1]] <- goNames
  dimnames(ar)[[2]] <- goNames
  
  names(gGOcat) <- goNames
  
  
  argGOcat <- list()
  
  argGOcat[[1]] <- ar
  
  argGOcat[[2]] <- gGOcat
  
  argGOcat
}



############################################# beta2 ################################################################


make_beta2matrix2d <- function(indata){
  
  dim(indata <- indata[indata[,3] <= 0.05,])
  library(org.Mm.eg.db); library(gplots); library(bioDist)
  # get gene in GO cat
  gGOcat <- list()
  for(i in 1:dim(indata)[1]){
    gGOcat[[indata[,1][i]]] <- as.vector(base::unlist(mget(indata[,1][i], org.Mm.egGO2ALLEGS)))
  }
  
  # filter each GO cat by the genes found to be DE
  f <- function(vecGO, vecDE){
    return(vecGO[which(vecGO %in% vecDE)])
  }
  
  gGOcat2 <- list()
  
  gGOcat <- lapply(gGOcat, FUN=function(x){f(x,GOstats_upgenes[[2]])})
  ar <- array()
  ar <- array(ar, c(length(gGOcat),length(gGOcat)), dimnames=list(names(gGOcat),names(gGOcat)))
  for(i in 1:length(gGOcat)) {
    for(j in 1:length(gGOcat)) {
      
      # intersect / union
      ar[i,j] <- (length(Reduce(intersect, list(gGOcat[[i]], gGOcat[[j]]))) / length(Reduce(union, list(gGOcat[[i]], gGOcat[[j]]))))
      
    }
  }
  dimnames(ar)[[1]] <- names(gGOcat)
  dimnames(ar)[[2]] <- names(gGOcat)
  
  
  goNames <- as.vector(base::unlist(mget(names(gGOcat), GOTERM)))
  goNames <- as.vector(base::unlist(lapply(goNames, FUN=function(x){x@Term})))
  
  goNames <- paste(goNames, names(gGOcat), sep=" - ")
  dimnames(ar)[[1]] <- goNames
  dimnames(ar)[[2]] <- goNames
  
  names(gGOcat) <- goNames
  
  
  argGOcat <- list()
  
  argGOcat[[1]] <- ar
  
  argGOcat[[2]] <- gGOcat
  
  argGOcat
}

#################################### beta1 down ##########################################

makedown_beta1matrix2d <- function(indata){
  
  dim(indata <- indata[indata[,3] <= 0.05,])
  library(org.Mm.eg.db); library(gplots); library(bioDist)
  # get gene in GO cat
  gGOcat <- list()
  for(i in 1:dim(indata)[1]){
    gGOcat[[indata[,1][i]]] <- as.vector(base::unlist(mget(indata[,1][i], org.Mm.egGO2ALLEGS)))
  }
  
  # filter each GO cat by the genes found to be DE
  f <- function(vecGO, vecDE){
    return(vecGO[which(vecGO %in% vecDE)])
  }
  
  gGOcat2 <- list()
  
  gGOcat <- lapply(gGOcat, FUN=function(x){f(x,GOstats_downgenes[[1]])})
  ar <- array()
  ar <- array(ar, c(length(gGOcat),length(gGOcat)), dimnames=list(names(gGOcat),names(gGOcat)))
  for(i in 1:length(gGOcat)) {
    for(j in 1:length(gGOcat)) {
      
      # intersect / union
      ar[i,j] <- (length(Reduce(intersect, list(gGOcat[[i]], gGOcat[[j]]))) / length(Reduce(union, list(gGOcat[[i]], gGOcat[[j]]))))
      
    }
  }
  dimnames(ar)[[1]] <- names(gGOcat)
  dimnames(ar)[[2]] <- names(gGOcat)
  
  
  goNames <- as.vector(base::unlist(mget(names(gGOcat), GOTERM)))
  goNames <- as.vector(base::unlist(lapply(goNames, FUN=function(x){x@Term})))
  
  goNames <- paste(goNames, names(gGOcat), sep=" - ")
  dimnames(ar)[[1]] <- goNames
  dimnames(ar)[[2]] <- goNames
  
  names(gGOcat) <- goNames
  
  
  argGOcat <- list()
  
  argGOcat[[1]] <- ar
  
  argGOcat[[2]] <- gGOcat
  
  argGOcat
}



############################################# beta2 down ################################################################


makedown_beta2matrix2d <- function(indata){
  
  dim(indata <- indata[indata[,3] <= 0.05,])
  library(org.Mm.eg.db); library(gplots); library(bioDist)
  # get gene in GO cat
  gGOcat <- list()
  for(i in 1:dim(indata)[1]){
    gGOcat[[indata[,1][i]]] <- as.vector(base::unlist(mget(indata[,1][i], org.Mm.egGO2ALLEGS)))
  }
  
  # filter each GO cat by the genes found to be DE
  f <- function(vecGO, vecDE){
    return(vecGO[which(vecGO %in% vecDE)])
  }
  
  gGOcat2 <- list()
  
  gGOcat <- lapply(gGOcat, FUN=function(x){f(x,GOstats_downgenes[[2]])})
  ar <- array()
  ar <- array(ar, c(length(gGOcat),length(gGOcat)), dimnames=list(names(gGOcat),names(gGOcat)))
  for(i in 1:length(gGOcat)) {
    for(j in 1:length(gGOcat)) {
      
      # intersect / union
      ar[i,j] <- (length(Reduce(intersect, list(gGOcat[[i]], gGOcat[[j]]))) / length(Reduce(union, list(gGOcat[[i]], gGOcat[[j]]))))
      
    }
  }
  dimnames(ar)[[1]] <- names(gGOcat)
  dimnames(ar)[[2]] <- names(gGOcat)
  
  
  goNames <- as.vector(base::unlist(mget(names(gGOcat), GOTERM)))
  goNames <- as.vector(base::unlist(lapply(goNames, FUN=function(x){x@Term})))
  
  goNames <- paste(goNames, names(gGOcat), sep=" - ")
  dimnames(ar)[[1]] <- goNames
  dimnames(ar)[[2]] <- goNames
  
  names(gGOcat) <- goNames
  
  
  argGOcat <- list()
  
  argGOcat[[1]] <- ar
  
  argGOcat[[2]] <- gGOcat
  
  argGOcat
}