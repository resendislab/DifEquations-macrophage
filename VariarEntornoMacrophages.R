rm(list = ls())

library(deSolve) # If odesolve is not already loaded
library(rootSolve)
library(BoolNetPerturb)

source("ODE_MacrophagePolarizationProbFinalComplete.R")
Pheno <- function(x){
  lbl <- ""
  #    if (!x["TCR"]){
  #        lbl <- "naive"}
  if (!(x["NFKB"]||x['STAT1']||x['STAT3']||x['STAT6']||x['HIF1A']||x['ERK']||x['Fra1']||x['AP1'])){
    lbl <- "M0"}
  if ((x['NFKB'] || x['STAT1'] || (x['TNFA'] && x['AP1']) || (x['TNFAe'] && x['AP1']) )){
    lbl <- "M1"}                                                    
  if ((x['STAT6'])){
    lbl <- "M2a"}        
  if (((x['IL1B'] && x['AP1']) || x['ERK'])){
    lbl <- "M2b"}
  if ((x['STAT3'])){
    lbl <- "M2c"}
  if (((x['TLR4'] && x['A2a']) || (x['Fra1'] && x['AP1']) || x['HIF1A'] || x['Fra1'])){
    lbl <- "M2d"}
  if (lbl == ""){
    lbl <- "NoLabel"}
  return(lbl)
}

dir.create("Results-M2d", showWarnings = FALSE)
#dir.create("Results/csv", showWarnings = FALSE)
#dir.create("Results/svg", showWarnings = FALSE)


Attractors <- read.csv("AttractoresPrueba.csv", header=TRUE)
rownames(Attractors) <- as.character(Attractors$Label)
Attractors$Label <- NULL
genes <- colnames(Attractors)
Attractors

labels.rules <- read.csv("Labels.csv")
#labels.rules$rules <- 
#    unlist(lapply(labels.rules$rules, function(x) 
#    gsub("([:|:]|[:&:]|[:!:]|[:(:]|[:):])"," ",x)))

#    unlist(lapply(labels.rules$rules, function(x) 
#        gsub("([:|:]|[:&:]|[:!:]).*","\\1",x)))

AttractorsMat <- Attractors
transitionMatrix <- data.frame(matrix(NA, nrow = sum(Attractors), ncol = 7))
t <- 1

Kpars <- rep(1,length(genes))
names(Kpars) <- paste("alpha", genes, sep="")
Parms <- c(h = 25, Kpars)
Parms<-c(b=0.5,Parms)


# Genes to modify
EnvGen <- c("HIF1A","STAT3")

for(Attri in 1:nrow(AttractorsMat)) {
  InitialState <- as.numeric(AttractorsMat[Attri,])
  InitialLabel <- rownames(AttractorsMat[Attri,])
  names(InitialState) <- genes
  
  Ks <- seq(0,1, length=11)
  Ks <- expand.grid( rep(list(Ks), length(EnvGen)) )
  colnames(Ks)<- EnvGen
  
  AttrSums <- numeric( dim(Ks)[1] )
  AttrLabels <- character( dim(Ks)[1] )
  AttrMatrix <- matrix(0, dim(Ks)[1], length(InitialState)+2)
  #rownames(AttrMatrix) <- as.character(Ks)
  colnames(AttrMatrix) <- c(names(InitialState), "Sum", "Label" )
  
  for (i in 1:dim(Ks)[1]) {
    State <- InitialState
    State[EnvGen] <- as.double(Ks[i,])
    StateF <- runsteady(y = State, fun = MacrophagePolarizationProbFinalComplete, parms = Parms, times = c(0, 1e5))$y
    # Final - Initial - Modified Extrinsic
    Sum <- sum(abs(StateF - InitialState)) - sum(Ks[i,])
    StateF.logi <- (StateF >= 0.75)
    label <- Pheno(StateF.logi)
    AttrMatrix[i,] <- c(StateF,Sum,label)
  }
  
  write.csv(AttrMatrix, paste("Results-M2d/", InitialLabel, "_",
                              paste(EnvGen, collapse = "_"), ".csv", sep="") )
  
}