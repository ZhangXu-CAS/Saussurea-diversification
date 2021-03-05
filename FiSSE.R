##Used by Jacob Landis to look at speciation rates of Polemoniaceae between annual and perennial species
##Used by Xu Zhang for Saussurea binary traits dependent diversification

#perform a FiSSE analysis with binary traits
library(phytools)
library(ape)
library(geiger)
library(nlme)
library(diversitree)

#setwd("~/FiSSE")
#need the R script traitDependent_functions.R to complete the run
source("traitDependent_functions.R")

#starting files of the tree and the data file
mydata <- read.table("Stem.txt",row.names=1,header = TRUE)
mytree <- read.tree("Saussurea.tre")

#compares names between the tree and the data to list any discrepancies
comparison <- name.check(phy=mytree,data=mydata)
comparison
# prune taxa that don't have data but are present in the tree
mytree <- drop.tip(mytree,comparison$tree_not_data)

#double check to make sure that taxa all match with tree and data
name.check(phy=mytree,data=mydata)

#create vector for traits
# if the vector names are changed, then commands following this will need to be altered to maintain the procedure
states <- mydata[,2]
names(states) <- row.names(mydata)

FISSE.res <- FISSE.binary(mytree, states)
FISSE.res
