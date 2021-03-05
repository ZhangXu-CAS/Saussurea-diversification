#script used by Jacob Landis to perform a HiSSE analysis with habit coded as Stem vs Stemless
# This script tests 25 different models in the HiSSE, BiSSE and null framework; then 
# produces boxplots to observe differences. Last ran January 2019
# Modified by Xu Zhang in July 2020

#R libraries needed
library(phytools)
library(geiger)
library(nlme)
library(hisse)
#library(devtools)
#devtools::install_github("thej022214/hisse")
#install.packages("hisse")

#starting files of the tree and the data file
mydata <- read.table("Stem.txt",row.names=1,header = TRUE)
mytree <- read.tree("Saussurea.tre")

#compares names between the tree and the data to list any discrepancies
comparison <- name.check(phy=mytree,data=mydata)
comparison

# prune taxa that don't have data but are present in the tree
mytree <- drop.tip(mytree,comparison$tree_not_data)
#write.tree(mytree,"Sau_drop.tre")
#double check to make sure that taxa all match with tree and data
name.check(phy=mytree,data=mydata)
#comparison <- name.check(phy=mytree,data=mydata)


#output types can be any of the following: "turnover", "net.div", or "raw"


###################################################################################
#####run the full model with different turnover.anc, eps.anc, and transition rates run 1
#full hisse model
cat("Running Full HiSSE Model 1 \n")
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual

pp1 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual,output.type="raw")

pp1 ### added the run numbers for model testing
###################################################################################
###Bissee model all free, run 2
cat("Bissee model all free, run 2 \n")
turnover.anc = c(1,2,0,0)
eps.anc = c(1,2,0,0)

#if one wanted to run a Bisse model in Hisse
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
trans.rates.bisse

pp2 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse)

pp2
###################################################################################
###Bissee model, extinction 0=1, run 3
cat("Bissee model, extinction 0=1, run 3 \n")

turnover.anc = c(1,2,0,0)
eps.anc = c(1,1,0,0)

#if one wanted to run a Bisse model in Hisse
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
trans.rates.bisse

pp3 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse)

pp3
###################################################################################
###Bissee model with equal q's, run 4
cat("Bissee model with equal q's, run 4 \n")

turnover.anc = c(1,2,0,0)
eps.anc = c(1,2,0,0)

#if one wanted to run a Bisse model in Hisse
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
trans.rates.bisse
trans.rates.bisse.equal = ParEqual(trans.rates.bisse, c(1,2))
trans.rates.bisse.equal

pp4 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse.equal)

pp4
###################################################################################
###Bissee model with equal q's and e0=e1, run 5
cat("Bissee model with equal q's and e0=e1, run 5 \n")

turnover.anc = c(1,2,0,0)
eps.anc = c(1,1,0,0)

#if one wanted to run a Bisse model in Hisse
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
trans.rates.bisse
trans.rates.bisse.equal = ParEqual(trans.rates.bisse, c(1,2))
trans.rates.bisse.equal

pp5 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse.equal)

pp5
###################################################################################

#2-state character independent CID-2 model, equal q's but different transition rates, run 6
#null two model, traits A and B have a diversification rate
cat("2-state character independent CID-2 model, equal q's but different transition rates, run 6 \n")

turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,2,2)

#full 8 transition model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))

#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp6 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp6
###################################################################################
#2-state character independent CID-2 model, equal q's and e's but different transition rates, run 7
#null two model, traits A and B have a diversification rate
cat("2-state character independent CID-2 model, equal q's and e's but different transition rates, run 7 \n")

turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,1,1)

#full 8 transition model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))

#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp7 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp7
###################################################################################
######CD-4 model, q's equal, run 8

#full 8 transition model
cat("full 8 transition model \n")

trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))

#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp8.hisse.null4 <- hisse(mytree, mydata, f=c(0.5,0.5), turnover.anc=rep(c(1,2,3,4),2),eps.anc=rep(c(1,2,3,4),2), trans.rate=trans.rates.nodual.allequal)

pp8.hisse.null4
###################################################################################
######CD-4 model, q's and e's equal, run 9
cat("CD-4 model, q's and e's equal, run9 \n")

eps.anc = c(1,1,1,1)
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))

#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp9.hisse.null4 = hisse(mytree, mydata, f=c(0.5,0.5), turnover.anc=rep(c(1,2,3,4),2),eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal)

pp9.hisse.null4
###################################################################################
#####run the full model with transition rates equal, run 10
#full hisse model
cat("run the full model with transition rates equal, run 10 \n")
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp10 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp10
###################################################################################
#####run the full model with transition rates and e's equal, run 11
cat("run the full model with transition rates and e's equal, run 11 \n")

turnover.anc = c(1,2,3,4)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp11 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp11
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=1a=0b, and extinction 0a=1a=0b run 12
cat("run the full model with transition rates equal, and turnover 0a=1a=0b, and extinction 0a=1a=0b run 12 \n")

turnover.anc = c(1,1,1,2)
eps.anc = c(1,1,1,2)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp12 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp12
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=1a=0b, and extinction equal run 13
cat("run the full model with transition rates equal, and turnover 0a=1a=0b, and extinction equal run 13 \n")

turnover.anc = c(1,1,1,2)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp13 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp13
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=0b, and e0a=e0b, run 14
cat("run the full model with transition rates equal, and turnover 0a=0b, and e0a=e0b, run 14 \n")

turnover.anc = c(1,2,1,3)
eps.anc = c(1,2,1,3)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp14 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp14
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=0b, and extinction equal, run 15
cat("run the full model with transition rates equal, and turnover 0a=0b, and extinction equal, run 15 \n")

turnover.anc = c(1,2,1,3)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp15 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp15
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=1a, and e0a=e1a, run 16
cat("run the full model with transition rates equal, and turnover 0a=1a, and e0a=e1a, run 16 \n")

turnover.anc = c(1,1,2,3)
eps.anc = c(1,1,2,3)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp16 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp16
###################################################################################
#####run the full model with transition rates equal, and turnover 0a=1a, and extinction equal, run 17
cat("run the full model with transition rates equal, and turnover 0a=1a, and extinction equal, run 17 \n")

turnover.anc = c(1,1,2,3)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

pp17 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp17
###################################################################################
#####run the full model with different turnover.anc, eps.anc, and but rates q0b1b=0,q1b0b=0, all other equals, run 18
cat("run the full model with different turnover.anc, eps.anc, and but rates q0b1b=0,q1b0b=0, all other equals, run 18 \n")
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp18 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp18
###################################################################################
#####run the full model with different e's equal and but rates q0b1b=0,q1b0b=0, all other equals, run 19
cat("run the full model with different e's equal and but rates q0b1b=0,q1b0b=0, all other equals, run 19 \n")

turnover.anc = c(1,2,3,4)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp19 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp19
###################################################################################
#####run the full model with r0a=r1a=r0b, e0a=e0b=e0b and rates q0b1b=0,q1b0b=0, all other equals, run 20
cat("run the full model with r0a=r1a=r0b, e0a=e0b=e0b and rates q0b1b=0,q1b0b=0, all other equals, run 20 \n")

turnover.anc = c(1,1,1,2)
eps.anc = c(1,1,1,2)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp20 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp20
###################################################################################
#####run the full model with r0a=r1a=r0b, extinction equals and rates q0b1b=0,q1b0b=0, all other equals, run 21
cat("run the full model with r0a=r1a=r0b, extinction equals and rates q0b1b=0,q1b0b=0, all other equals, run 21 \n")

turnover.anc = c(1,1,1,2)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp21 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp21
###################################################################################
#####run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 22
cat("run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 22 \n")

turnover.anc = c(1,2,1,3)
eps.anc = c(1,2,1,3)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp22 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp22
###################################################################################
#####run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 23
cat("run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 23 \n")

turnover.anc = c(1,2,1,3)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp23 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp23
###################################################################################
#####run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 24
cat("run the full model withr0a=r0b, e0a=e0b, and rates q0b1b=0,q1b0b=0, all other equals, run 24 \n")

turnover.anc = c(1,1,2,3)
eps.anc = c(1,1,2,3)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp24 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp24
###################################################################################
#####run the full model withr0a=r0b, extinction equal, and rates q0b1b=0,q1b0b=0, all other equals, run 25
cat("run the full model withr0a=r0b, extinction equal, and rates q0b1b=0,q1b0b=0, all other equals, run 25 \n")

turnover.anc = c(1,1,2,3)
eps.anc = c(1,1,1,1)

#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,9,10,12))
trans.rates.nodual
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6))
trans.rates.nodual.allequal

pp25 = hisse(mytree, mydata, f=c(0.5,0.5), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")

pp25
###################################################################################
##Write a dataframe for comparison of modles

Model_test_leave_col<- data.frame(loglik= c(pp1[["loglik"]],pp2[["loglik"]],pp3[["loglik"]],pp4[["loglik"]],pp5[["loglik"]],pp6[["loglik"]],
								  pp7[["loglik"]],pp8.hisse.null4[["loglik"]],pp9.hisse.null4[["loglik"]],pp10[["loglik"]],pp11[["loglik"]],pp12[["loglik"]],
								  pp13[["loglik"]],pp14[["loglik"]],pp15[["loglik"]],pp16[["loglik"]],pp17[["loglik"]],pp18[["loglik"]],
								  pp19[["loglik"]],pp20[["loglik"]],pp21[["loglik"]],pp22[["loglik"]],pp23[["loglik"]],pp24[["loglik"]],pp25[["loglik"]]),
			        AIC= c(pp1[["AIC"]],pp2[["AIC"]],pp3[["AIC"]],pp4[["AIC"]],pp5[["AIC"]],pp6[["AIC"]],
								  pp7[["AIC"]],pp8.hisse.null4[["AIC"]],pp9.hisse.null4[["AIC"]],pp10[["AIC"]],pp11[["AIC"]],pp12[["AIC"]],
								  pp13[["AIC"]],pp14[["AIC"]],pp15[["AIC"]],pp16[["AIC"]],pp17[["AIC"]],pp18[["AIC"]],
								  pp19[["AIC"]],pp20[["AIC"]],pp21[["AIC"]],pp22[["AIC"]],pp23[["AIC"]],pp24[["AIC"]],pp25[["AIC"]]))

write.csv(Model_test_leave_col,"Model_test_leave_col.csv")

#Run MarginRecon after rerunning the best fit model
pp.recon_leave_col <- MarginRecon(mytree, mydata, f=c(0.5,0.5), pars=pp1$solution, hidden.states=TRUE)

save(pp.recon_leave_col, file="full_HiSSE_recon_leave_col.RData")


rates_leave_col<- GetModelAveRates(pp.recon_leave_col, type = "tips")

save(rates_leave_col, file="Tip_rates_fullHiSSE_model_leave_col.RData")

write.table(rates_leave_col, file="Tip_rates_fullHiSSE_model_leave_col.txt",quote=FALSE,sep="\t",row.names=FALSE)

#after combinging state calls into one new column called states
rates <- read.table("Tip_rates_fullHiSSE_model_leave_col.txt", header=TRUE)

pdf(file="leave_col_recon.pdf",20,20)
plot.hisse.states(pp.recon_leave_col, rate.param="net.div", show.tip.label=TRUE, legend="none",
                  fsize=0.5,do.observed.only=TRUE, rate.colors=c("steelblue","tomato"))
legend("topleft", legend=c("Same","Different"), fill=c("steelblue","tomato"),cex=2,title="States")
dev.off()

pdf(file="Boxpot_netdiv_leave_col.pdf")
boxplot(rates$net.div~rates$state, boxwex=0.5, 
        notch = TRUE,main="Net Diversification",
        col = c("steelblue","tomato"), 
        xlab="Same vs. Different", ylab="Diversification rate")
dev.off()

pdf(file="Boxplot_speciation_leave_col.pdf")
boxplot(rates$speciation~rates$state, boxwex=0.5,
        notch = TRUE,main="Speciation",
        col = c("steelblue","tomato"), 
        xlab="Same vs. Different", ylab="Speciation rate")
dev.off()

pdf(file="Boxplot_extinction_leave_col.pdf")
boxplot(rates$extinction~rates$state,boxwex=0.5, 
        notch = TRUE,main="Extinction", 
        col = c("steelblue","tomato"),
        xlab="Same vs. Different", ylab="Extinction rate")
dev.off()


