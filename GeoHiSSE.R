#script used by Jacob Landis to perform a GeoHiSSE analysis with with species in the CA-FP
# This script tests 4 different models including the standard GeoSSE and GeoHiSSE framework 
# then produces boxplots to observe differences. 
# modified bu Xu Zhang

library( devtools )
#install_github(repo = "thej022214/hisse", ref = "master")
library(hisse)
library(diversitree)
library(ape)
library(geiger)

#setwd("~/GeoHiSSE")

#starting files of the tree and the data file
mydata <- read.table("habitat2.txt",row.names=1,header=TRUE)
#mytree <- read.nexus("MCC.fixed.nex")
mytree <- read.tree("Saussurea.tre")

#compares names between the tree and the data to list any discrepancies
comparison <- name.check(phy=mytree,data=mydata)
comparison

# prune taxa that don't have data but are present in the tree
mytree <- drop.tip(mytree,comparison$tree_not_data)

#double check to make sure that taxa all match with tree and data
name.check(phy=mytree,data=mydata)


f=c(0.5,0.5,0.5)

## Model 1 - Dispersal parameters vary only, no range-dependent diversification.
turnover <- c(1,1,0)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas = 0)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod1 <-  GeoHiSSE(phy = mytree, data = mydata, f=f,
                  speciation=turnover, extirpation=eps,
                  hidden.areas =FALSE, trans.rate=trans.rate.mod,
                   speciation.upper=100, trans.upper=10)

## Model 2. Canonical GeoSSE model, range effect on diversification 
turnover <- c(1,2,3)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod2 <- GeoHiSSE(phy = mytree, data = mydata, f=f,
                 speciation=turnover, extirpation=eps,
                 hidden.areas=FALSE, trans.rate=trans.rate.mod,
                 speciation.upper=100, trans.upper=10)
## Model 3. GeoHiSSE model with 1 hidden area, no range-dependent diversification.
## Note below how parameters vary among hidden classes but are the same within each 
##      hidden class.
turnover <- c(1,1,0,2,2,0)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1, make.null=TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod3 <- GeoHiSSE(phy = mytree, data = mydata, f=f,
                 speciation=turnover, extirpation=eps,
                 hidden.areas=TRUE, trans.rate=trans.rate.mod,
                 speciation.upper=100, trans.upper=10)

## Model 4. GeoHiSSE model with 1 hidden area, no range-dependent diversification.

turnover <- c(1,2,3,4,5,6)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=1)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod4 <- GeoHiSSE(phy = mytree, data =mydata, f=f,
                 speciation=turnover, extirpation=eps,
                 hidden.areas=TRUE, trans.rate=trans.rate.mod,
                 speciation.upper=100, trans.upper=10)

###Model 5. Heterogeneous diversification, not tied to range evolution. Assumes 5 distinct diversification rates.
## Model 5. MuSSE-like model with no hidden trait, no cladogenetic effects.

#turnover <- c(1,2,0)
#eps <- c(1,1)
#trans.rate <- TransMatMakerGeoHiSSE(hidden.areas=0, make.null=FALSE,
#                                    separate.extirpation = TRUE)
#trans.rate.mod <- ParEqual(trans.rate, c(1,2))
#trans.rate.mod <- ParEqual(trans.rate.mod, c(2,3))
#mod5 <- GeoHiSSE(phy = mytree, data = sim.dat, f=f,
#                 speciation=turnover, extirpation=eps,
#                 hidden.areas=FALSE, trans.rate=trans.rate.mod,
#                 speciation.upper=100, trans.upper=10, sann=FALSE,
#                 assume.cladogenetic = FALSE)
#
GetAICWeights(list(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4), criterion="AIC")
    
    
#GetModelWeight(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4)
geohiss_model_test<- data.frame(loglik=c(mod1[["loglik"]],mod2[["loglik"]],mod3[["loglik"]],mod4[["loglik"]]),
                                AIC= c(mod1[["AIC"]],mod2[["AIC"]],mod3[["AIC"]],mod4[["AIC"]]))
geohiss_model_test
write.csv(geohiss_model_test,"geohiss_model_test.csv")
## As the number of models in the set grows, naming each model in the set can become hard.
## So one can use a list (created by some automated code) as an imput also:


recon.mod <- MarginReconGeoSSE(phy = mod4$phy, data = mod4$data, f = mod4$f,
                                pars = mod4$solution,
                                root.type = mod4$root.type, root.p = mod4$root.p,
                                aic = mod4$AIC)


rates <- GetModelAveRates(recon.mod, type = "tips")

save(rates, file="Tip_rates_best_model.Rdata")

write.csv(rates, file="Tip_rates.csv")

#after combinging state calls into one new column called states
rates <- read.csv("Tip_rates.csv", header=TRUE)

pdf(file="Boxpot_netdiv.pdf")

boxplot(rates$net.div~rates$state, boxwex=0.5, 
        notch = F,main="Net Diversification",
        col = c("steelblue","tomato","gold"), 
        xlab=" Widespread vs. Alpine vs. Lowland ", ylab="Diversification rate")
dev.off()

pdf(file="Boxplot_speciation.pdf")
boxplot(rates$speciation~rates$state, boxwex=0.5, 
        notch = F,main="Speciation",
        col = c("steelblue","tomato","gold"), 
        xlab=" Widespread vs. Alpine vs. Lowland ", ylab="Speciation rate")
dev.off()

pdf(file="Boxplot_extinction.pdf")
boxplot(rates$extinction~rates$state, boxwex=0.5, 
        notch = F,main="Extinction",
        col = c("steelblue","tomato","gold"), 
        xlab=" Widespread vs. Alpine vs. Lowland ", ylab="Extinction rate")
dev.off()

#plot help is from plot.hisse.states function
#plot tree with geographic states only

dev.off()