library(BoolNet)
library(BoolNetPerturb)


macrophage<-loadNetwork("macrophage2.txt")
Macrophage.Phenotypes<-data.frame(labels=c('M0', 'M1', 'M2a','M2b','M2c', 'M2d'), rules=c('!(NFKB | STAT1 |STAT3 |STAT6|HIF1A|ERK|Fra1|AP1)', 'NFKB|STAT1|TNFA & AP1| TNFAe & AP1', 'STAT6', '  (IL1B & AP1)| ERK','STAT3', '(TLR4 & A2a) | (Fra1 & AP1) | HIF1A| Fra1'))
#We Obtain the attractors
attr.s <- getAttractors(macrophage, method = "sat.exhaustive", type="synchronous", returnTable=TRUE)
attr.table.s <- attractorToDataframe(attr.s, Boolean=TRUE)
#Label Attractors 
attr.labels.s<-labelAttractors(attr.s, Macrophage.Phenotypes, macrophage$genes)
attr.table.s<-merge(x=attr.table.s, y=attr.labels.s, by.x='attractor', by.y=0)
colnames(attr.table.s)[colnames(attr.table.s)=="y"]<-"label"
attr.table.s$size<-stringr::str_count(attr.table.s$label,'/')+1
attr.table.s<-attr.table.s[c(c("label", "attractor","state","size"),macrophage$genes)]
write.csv(attr.table.s, "FinalAttractors.csv", row.names=FALSE)

## We apply the transformation to the following boolean model 
####
macrophage<-loadNetwork("macrophage2.txt")

cont_simulation <- function(macrophage, LOGIC="Probabilistic", EQ="Villarreal"){
  ODEnetwork <- booleanToODE(macrophage, logic=LOGIC, eq = EQ)
}
macrophageODE<-cont_simulation(macrophage)
##We obtain the set of differential equations
capture.output(macrophageODE$func, file= "Output_func.txt")
###########################################################################################################################
###########################################################################################################################
library(deSolve)
library(rootSolve)
library(plyr)
# load ode network
source("ODE_MacrophagePolarizationProbFinal.R")
attr.bool <- read.csv("MacrophagePolarization_Final.csv", stringsAsFactors=FALSE)
node.names <- names(attr.bool)[-length(attr.bool)]
parms <- c(rep(1,length(node.names)-15), rep(0,15))
names(parms) <- paste("alpha", node.names, sep="")
parms <- c(h=25, parms)  
parms<-c(b=0.5,parms)

###########################################################################################################################
###########################################################################################################################
## Verify attractors 
df <- apply(attr.bool, 1, function(s0) {
    label <- s0[length(s0)]
    s0 <- s0[-length(s0)]
    s0 <- gsub("\\*", "0", s0)
    s0 <- sapply(s0, as.numeric)
    sf <- runsteady(y = s0, fun = MacrophagePolarizationProbFinal, parms = parms, times = c(0, 1e5))$y
    dif <- sum(abs(sf-s0))
    c(label,sf,difference=dif)
})
write.csv(t(df), 'Diff_bool_contProb.csv')
###########################################################################################################################
###########################################################################################################################
###M0 ONE-NODE
M0 = rep(0, length(node.names))
names(M0) <- node.names
values <- seq(0,1,by=.025)  


dir.create("Results-M0-1node", showWarnings = FALSE)
dir.create("Results-M0-1node/csv", showWarnings = FALSE)

extrinsic <- node.names[16:length(node.names)]

for (env in extrinsic) {
    env.index <- which(node.names==env)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M0
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProbFinal, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M0-1node/csv/M0_",env,'.csv') )
    write.csv(res, paste0("Results-M0-1node/csv/M0Prob_",env,'.csv'))
}
###########################################################################################################################
###########################################################################################################################
dir.create("Results-M0-env", showWarnings = FALSE)
dir.create("Results-M0-env/csv", showWarnings = FALSE)

environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M0
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProbFinal, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M0-env/csv/M0_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M0-env/csv/M0Prob_",paste0(env,collapse='+'),'.csv') )
}
###########################################################################################################################
###########################################################################################################################
##Perturbating the interleukins but in the opposite direction 
###M0 ONE-NODE
M0 = rep(0, length(node.names))
values <- seq(1,0,by=-.025)


dir.create("Results-M0-1node", showWarnings = FALSE)
dir.create("Results-M0-1node/csvopposite", showWarnings = FALSE)

extrinsic <- node.names[15:length(node.names)]

for (env in extrinsic) {
    env.index <- which(node.names==env)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M0
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProb, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M0-1node/csv/M0_",env,'.csv') )
    write.csv(res, paste0("Results-M0-1node/csvopposite/M0_",env,'.csv'))
}
###########################################################################################################################
###########################################################################################################################
##Perturbing The M1 phenotype and how resilent it is to change
###But first we must adjust the system of differential equations 
# load ode network
source("ODE_MacrophagePolarizationProbFinalComplete.R")
attr.bool <- read.csv("MacrophagePolarization_Final.csv", stringsAsFactors=FALSE)
node.names <- names(attr.bool)[-length(attr.bool)]
parms <- c(rep(1,length(node.names)-29), rep(1,29))
names(parms) <- paste("alpha", node.names, sep="")
parms <- c(h=25, parms)  
parms<-c(b=0.5,parms)
dir.create("Results-M1-env", showWarnings = FALSE)
dir.create("Results-M1-env/csv", showWarnings = FALSE)
M1 = c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M1) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values-1), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    res[1,7] <- 1
    res[1,8] <- 1
    res[1,11] <- 1
    for (i in seq(length(values))) {
        s0 <- M1
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProbFinalComplete, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M1-env/csv/M1_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M1-env/csv/M1Prob_",paste0(env,collapse='+'),'.csv') )
}
###########################################################################################################################
###########################################################################################################################
##Perturbing The M2aM2cM2d phenotype and how resilent it is to change
dir.create("Results-M2aM2cM2d-env", showWarnings = FALSE)
dir.create("Results-M2aM2cM2d-env/csv", showWarnings = FALSE)
M2aM2cM2d = c(1,1,1,0,1,0,0,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M2aM2cM2d) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M2aM2cM2d
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProbFinalComplete, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M2aM2cM2d-env/csv/M2aM2cM2d_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M2aM2cM2d-env/csv/M2aM2cM2dProb_",paste0(env,collapse='+'),'.csv') )
}
###########################################################################################################################
###########################################################################################################################
##Perturbing The M1M2bM2cM2d phenotype and how resilent it is to change
dir.create("Results-M1M2bM2cM2d-env", showWarnings = FALSE)
dir.create("Results-M1M2bM2cM2d-env/csv", showWarnings = FALSE)
M1M2bM2cM2d = c(1,1,0,0,1,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M1M2bM2cM2d) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M1M2bM2cM2d
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProbFinalComplete, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M1M2bM2cM2d-env/csv/M1M2bM2cM2d_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M1M2bM2cM2d-env/csv/M1M2bM2cM2dProb_",paste0(env,collapse='+'),'.csv') )
}

###########################################################################################################################
###########################################################################################################################
##Perturbing The M2cM2d phenotype and how resilent it is to change
dir.create("Results-M2cM2d-env", showWarnings = FALSE)
dir.create("Results-M2cM2d-env/csv", showWarnings = FALSE)
M2cM2d = c(1,1,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M2cM2d) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M2cM2d
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProbFinalComplete, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M2cM2d-env/csv/M2cM2d_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M2cM2d-env/csv/M2cM2dProb_",paste0(env,collapse='+'),'.csv') )
}

###########################################################################################################################
###########################################################################################################################
##Perturbing The M1M2bM2d phenotype and how resilent it is to change
dir.create("Results-M1M2bM2d-env", showWarnings = FALSE)
dir.create("Results-M1M2bM2d-env/csv", showWarnings = FALSE)
M1M2bM2d = c(0,0,0,0,1,1,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M1M2bM2d) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M1M2bM2d
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProbFinalComplete, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M1M2bM2d-env/csv/M1M2bM2d_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M1M2bM2d-env/csv/M1M2bM2dProb_",paste0(env,collapse='+'),'.csv') )
}
###########################################################################################################################
###########################################################################################################################
##Perturbing The M2aM2d phenotype and how resilent it is to change
dir.create("Results-M2aM2d-env", showWarnings = FALSE)
dir.create("Results-M2aM2d-env/csv", showWarnings = FALSE)
M2aM2d = c(0,1,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M2aM2d) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M2aM2d
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProbFinalComplete, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M2aM2d-env/csv/M2aM2d_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M2aM2d-env/csv/M2aM2dProb_",paste0(env,collapse='+'),'.csv') )
}


###########################################################################################################################
###########################################################################################################################
##Perturbing The M2aM2bM2d phenotype and how resilent it is to change
dir.create("Results-M2aM2bM2d-env", showWarnings = FALSE)
dir.create("Results-M2aM2bM2d-env/csv", showWarnings = FALSE)
M2aM2bM2d = c(0,1,1,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M2aM2bM2d) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M2aM2bM2d
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProbFinalComplete, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M2aM2bM2d-env/csv/M2aM2bM2d_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M2aM2bM2d-env/csv/M2aM2bM2dProb_",paste0(env,collapse='+'),'.csv') )
}

###########################################################################################################################
###########################################################################################################################
##Perturbing The M2bM2d phenotype and how resilent it is to change
dir.create("Results-M2bM2d-env", showWarnings = FALSE)
dir.create("Results-M2bM2d-env/csv", showWarnings = FALSE)
M2bM2d = c(0,1,0,0,0,0,0,0,0,0,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M2bM2d) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M2bM2d
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProbFinalComplete, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M2bM2d-env/csv/M2bM2d_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M2bM2d-env/csv/M2bM2dProb_",paste0(env,collapse='+'),'.csv') )
}

###########################################################################################################################
###########################################################################################################################
##Perturbing The M1M2b phenotype and how resilent it is to change
dir.create("Results-M1M2b-env", showWarnings = FALSE)
dir.create("Results-M1M2b-env/csv", showWarnings = FALSE)
M1M2b = c(0,0,0,0,0,1,1,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M1M2b) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M1M2b
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProbFinalComplete, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M1M2b-env/csv/M1M2b_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M1M2b-env/csv/M1M2bProb_",paste0(env,collapse='+'),'.csv') )
}

###########################################################################################################################
###########################################################################################################################
##Perturbing The M1M2d phenotype and how resilent it is to change
dir.create("Results-M1M2d-env", showWarnings = FALSE)
dir.create("Results-M1M2d-env/csv", showWarnings = FALSE)
M1M2d = c(0,0,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M1M2d) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M1M2d
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProbFinalComplete, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M1M2d-env/csv/M1M2d_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M1M2d-env/csv/M1M2dProb_",paste0(env,collapse='+'),'.csv') )
}

###########################################################################################################################
###########################################################################################################################
##Perturbing The M2d phenotype and how resilent it is to change
dir.create("Results-M2d-env", showWarnings = FALSE)
dir.create("Results-M2d-env/csv", showWarnings = FALSE)
M2d = c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M2d) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M2d
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProbFinalComplete, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M2d-env/csv/M2d_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M2d-env/csv/M2dProb_",paste0(env,collapse='+'),'.csv') )
}

###########################################################################################################################
###########################################################################################################################
##Perturbing The M2b phenotype and how resilent it is to change
dir.create("Results-M2b-env", showWarnings = FALSE)
dir.create("Results-M2b-env/csv", showWarnings = FALSE)
M2b = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M2b) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M2b
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophagePolarizationProbFinalComplete, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M2b-env/csv/M2b_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M2b-env/csv/M2bProb_",paste0(env,collapse='+'),'.csv') )
}



#######################################################################################################################
#              Transforming our TGEM in a set of differential equations and evaluating
#                the behaviour in 4 types of tumor microenvironments
#######################################################################################################################
NfkbNoHif1a<-fixGenes(macrophage, c('NFKB', 'HIF1A'), c(1,0))
cont_simulation <- function(NfkbNoHif1a, LOGIC="Probabilistic", EQ="Villarreal"){
  ODEnetwork <- booleanToODE(NfkbNoHif1a, logic=LOGIC, eq = EQ)
}
macrophageTGEM<-cont_simulation(NfkbNoHif1a)
##We obtain the set of differential equations
capture.output(macrophageTGEM$func, file= "Output_funcTGEM.txt")
###################################################################
#Obtain the attractors of our TGEM. 
 #We Obtain the attractors
attr.TGEM <- getAttractors(NfkbNoHif1a, method = "sat.exhaustive", type="synchronous", returnTable=TRUE)
attr.table.TGEM <- attractorToDataframe(attr.TGEM, Boolean=TRUE)
#Label Attractors 
attr.labels.TGEM<-labelAttractors(attr.TGEM, Macrophage.Phenotypes, NfkbNoHif1a$genes)
attr.table.TGEM<-merge(x=attr.table.TGEM, y=attr.labels.TGEM, by.x='attractor', by.y=0)
colnames(attr.table.TGEM)[colnames(attr.table.TGEM)=="y"]<-"label"
attr.table.TGEM$size<-stringr::str_count(attr.table.TGEM$label,'/')+1
attr.table.TGEM<-attr.table.TGEM[c(c("label", "attractor","state","size"),NfkbNoHif1a$genes)]
write.csv(attr.table.TGEM, "FinalAttractorsTGEM.csv", row.names=FALSE)
################################################################################################
#Lets simulate our TGEM to see how it behaves by modifying each extracellular node             #
#                                                                                              #
###############################################################################################
dir.create("Results-TGEMM0-1node", showWarnings = FALSE)
dir.create("Results-TGEMM0-1node/csv", showWarnings = FALSE)

extrinsic <- node.names[16:length(node.names)]

for (env in extrinsic) {
  env.index <- which(node.names==env)
  res <- matrix(0, length(values), length(node.names))
  colnames(res) <- node.names
  rownames(res) <- values
  for (i in seq(length(values))) {
    s0 <- M0
    s0[env.index] <- values[i]
    res[i,] <- runsteady(y = s0, fun = MacrophageTGEM, parms = parms, times = c(0, 1e3))$y
    #dif <- sum(abs(sf-s0))
    #print(dif)
  }
  print( paste0("Results-TGEMM0-1node/csv/M0_",env,'.csv') )
  write.csv(res, paste0("Results-TGEMM0-1node/csv/M0Prob_",env,'.csv'))
}
#######################################################################################################
#             Could we obtain desired phenotypes with our TGEM based on the                           #
#                               specific Macrophage phenotypes Microenvironments?                     #
######################################################################################################
dir.create("Results-TGEMMicro-env", showWarnings = FALSE)
dir.create("Results-TGEMMicro-env/csv", showWarnings = FALSE)

environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
  env.index <- match(env,node.names)
  res <- matrix(0, length(values), length(node.names))
  colnames(res) <- node.names
  rownames(res) <- values
  for (i in seq(length(values))) {
    s0 <- M0
    s0[env.index] <- values[i]
    res[i,] <- runsteady(y = s0, fun = MacrophageTGEM, parms = parms, times = c(0, 1e3))$y
    #dif <- sum(abs(sf-s0))
    #print(dif)
  }
  print( paste0("Results-TGEMMicro-env/csv/TGEM_",paste0(env,collapse='+'),'.csv') )
  write.csv(res, paste0("Results-TGEMMicro-env/csv/TGEMProb_",paste0(env,collapse='+'),'.csv') )
  
}

#######################################################################################################
#             Could we obtain desired phenotypes with our TGEM based on the                           #
#                               specific Breast Cancer Microenvironments?                             #
######################################################################################################
library(deSolve)
library(rootSolve)
library(plyr)
# load ode network
source("ODE_MacrophageTGEM.R")
attr.bool <- read.csv("MacrophagePolarizationTGEM_Final.csv", stringsAsFactors=FALSE)
node.names <- names(attr.bool)[-length(attr.bool)]
parms <- c(rep(1,length(node.names)-29), rep(1,29))
names(parms) <- paste("alpha", node.names, sep="")
parms <- c(h=25, parms)  
parms<-c(b=0.5,parms)


dir.create("Results-TGEM-env", showWarnings = FALSE)
dir.create("Results-TGEM-env/csv", showWarnings = FALSE)

environments2 <- read.csv("TumorMicroenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments2 <- apply(environments2, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments2) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M0
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophageTGEM, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-TGEM-env/csv/TGEM_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-TGEM-env/csv/TGEMProb_",paste0(env,collapse='+'),'.csv') )
    
    }
    

############################################################################################################################
#From our TGEM we will perturb specifically the phenotypes obtained based on how well they will performed on the 4 tumor microenvironments
##Perturbing The M1M2bM2d phenotype and how resilent it is to change
dir.create("Results-M1M2bM2dTGEM-env", showWarnings = FALSE)
dir.create("Results-M1M2bM2dTGEM-env/csv", showWarnings = FALSE)
M1M2bM2d = c(0,1,0,0,1,1,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M1M2bM2d) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M1M2bM2d
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophageTGEM, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M1M2bM2dTGEM-env/csv/M1M2bM2dTGEM_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M1M2bM2dTGEM-env/csv/M1M2bM2dTGEM_",paste0(env,collapse='+'),'.csv') )
}

############################################################################################################################
#From our TGEM we will perturb specifically the phenotypes obtained based on how well they will performed on the 4 tumor microenvironments
##Perturbing The M1M2b phenotype and how resilent it is to change
dir.create("Results-M1M2bTGEM-env", showWarnings = FALSE)
dir.create("Results-M1M2bTGEM-env/csv", showWarnings = FALSE)
M1M2b = c(0,1,0,0,0,1,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M1M2b) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M1M2b
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophageTGEM, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M1M2bTGEM-env/csv/M1M2bTGEM_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M1M2bTGEM-env/csv/M1M2bTGEM_",paste0(env,collapse='+'),'.csv') )
}


############################################################################################################################
#From our TGEM we will perturb specifically the phenotypes obtained based on how well they will performed on the 4 tumor microenvironments
##Perturbing The M1 phenotype and how resilent it is to change
dir.create("Results-M1TGEM-env", showWarnings = FALSE)
dir.create("Results-M1TGEM-env/csv", showWarnings = FALSE)
M1 = c(0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M1) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M1
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophageTGEM, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M1TGEM-env/csv/M1TGEM_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M1TGEM-env/csv/M1TGEM_",paste0(env,collapse='+'),'.csv') )
}

############################################################################################################################
#From our TGEM we will perturb specifically the phenotypes obtained based on how well they will performed on the 4 tumor microenvironments
##Perturbing The M1M2d phenotype and how resilent it is to change
dir.create("Results-M1M2dTGEM-env", showWarnings = FALSE)
dir.create("Results-M1M2dTGEM-env/csv", showWarnings = FALSE)
M1M2d = c(0,1,0,0,0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
names(M1M2d) <- node.names
values <- seq(0,1,by=.025) 
environments <- read.csv("Microenvironments.csv", header=F, row.names=1, stringsAsFactors=FALSE)
environments <- apply(environments, 1, function(s0) {s0 <- s0[s0!=""]})
values <- seq(0,1,by=.025) 

for (env in environments) {
    env.index <- match(env,node.names)
    res <- matrix(0, length(values), length(node.names))
    colnames(res) <- node.names
    rownames(res) <- values
    for (i in seq(length(values))) {
        s0 <- M1M2d
        s0[env.index] <- values[i]
        res[i,] <- runsteady(y = s0, fun = MacrophageTGEM, parms = parms, times = c(0, 1e3))$y
        #dif <- sum(abs(sf-s0))
        #print(dif)
    }
    print( paste0("Results-M1M2dTGEM-env/csv/M1M2dTGEM_",paste0(env,collapse='+'),'.csv') )
    write.csv(res, paste0("Results-M1M2dTGEM-env/csv/M1M2dTGEM_",paste0(env,collapse='+'),'.csv') )
}