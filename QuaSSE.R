####QuaSSE analysis --xu zhang 2021/1/9
#Remove all data

rm(list=ls())
#perform a QuaSSE analysis for scored traits in Polemoniaceae
library(diversitree)
library(phytools)
library(ape)
library(geiger)
library(nlme)

#starting files of the tree and the data file
mydata <- read.table("niche_breadth.txt",row.names=1,header = T)
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


#standard deviation of traits, for corolla length its 1.186, width is 0.522
states.sd <- sd(states)

######make models and analysis
####
p <- starting.point.quasse(mytree,states)
p

###Create a piecewise “linear” function
xr <- range(states) + c(-1,1) * 20 * p["diffusion"]
linear.x <- make.linear.x(xr[1], xr[2])

###make models, only fitted speciation functions to simplify the analysis 
make.models <- function(lambda, mu)
make.quasse(mytree, states, states.sd, lambda, mu)

###make.corolla
nodrift <- function(f)
constrain(f, drift ~0)
###Create the likelihood functions where speciation is a constant, linear, sigmoidal, or hump-shaped
f.c <- make.models(constant.x, constant.x)
f.l <- make.models(linear.x, constant.x)
f.s <- make.models(sigmoid.x, constant.x)
f.h <- make.models(noroptimal.x, constant.x)

###fit the constant model
control <- list(parscale=.1, reltol=0.001)
mle.c <- find.mle(nodrift(f.c), p, lower=0, control=control, verbose=0)

###Starting points for the constrained analyses based on this constrained fit.
p.c <- mle.c$par
p.l <- c(p.c[1], l.m=0, p.c[2:3])
p.s <- p.h <- c(p.c[1], p.c[1], mean(xr), 1, p.c[2:3])
names(p.s) <- argnames(nodrift(f.s))
names(p.h) <- argnames(nodrift(f.h))
mle.l <- find.mle(nodrift(f.l), p.l, control=control, verbose=0)
mle.s <- find.mle(nodrift(f.s), p.s, control=control, verbose=0)
#mle.h <- find.mle(nodrift(f.h), p.h, control=control, verbose=0)	

###Run the fits with the drift parameter added, starting from the constrained model’s ML parameters:

mle.d.l <- find.mle(f.l, coef(mle.l, TRUE), control=control, verbose=0)
mle.d.s <- find.mle(f.s, coef(mle.s, TRUE), control=control, verbose=0)
#mle.d.h <- find.mle(f.h, coef(mle.h, TRUE), control=control, verbose=0)

model.test<- anova(mle.c, linear=mle.l, sigmoidal=mle.s,drift.linear=mle.d.l, drift.sigmoidal=mle.d.s)
write.csv(model.test,"QuaSSE.niche_breadth.modeltest.csv")

coef(mle.l)
coef(mle.d.s)
