rm(list=ls())
library(ape)
library (diversitree)
library(phytools)
library(geiger)

phy<-read.tree("Saussurea.tre")

X <- read.table("Phyrows.txt", row.names = 1,header = T)

trait <- X$Phyrows
names(trait)<- row.names(X)
comparison <- name.check(phy=phy,data=trait)
comparison
phy <- drop.tip(phy,comparison$tree_not_data)
name.check(phy=phy,data=trait)
samplingf<-c(0.5, 0.5, 0.5,0.5)
p <- starting.point.musse(phy, k=4)
p

lik.musse<-make.musse(phy,trait,4,sampling.f=samplingf)
lik.musse
lik.null<-constrain(lik.musse,lambda2 ~ lambda1, lambda3 ~ lambda1,lambda4 ~ lambda1,
                    mu2 ~ mu1,  mu3 ~ mu1,mu4 ~ mu1)
lik.null
fit.null<-find.mle(lik.null,x.init=p[argnames(lik.null)])
fit.null
fit.full<-find.mle(lik.musse,x.init=p[argnames(lik.musse)])
fit.full

lik.lambda <- constrain(lik.musse, mu2 ~ mu1, mu3 ~ mu1,mu4 ~ mu1)
fit.lambda <- find.mle(lik.lambda, p[argnames(lik.lambda)])
fit.lambda

lik.mu <- constrain(lik.musse, lambda2 ~ lambda1, lambda3 ~ lambda1,lambda4 ~ lambda1)
fit.mu <- find.mle(lik.mu, p[argnames(lik.mu)])
fit.mu

AnovaResults <- anova(fit.null, 
                      all.different=fit.full, 
                      free.lambda=fit.lambda, 
                      free.mu=fit.mu)
AnovaResults
write.csv(AnovaResults,"Musse_modeltest_Phyrows.csv")
aicw(setNames(AnovaResults$AIC,rownames(AnovaResults)))
#############
prior <- make.prior.exponential(1/2)
prior
prelim <- mcmc(lik.lambda, fit.lambda$par, nsteps=100, prior=prior, 
               w=1, print.every=0) #### where the error occurs##
head(prelim)
## the following is how Fitzjohn recommends setting w
w <- diff(sapply(prelim[2:(ncol(prelim)-1)], quantile, c(0.05, 0.95)))
w
mcmc.fit.full<- mcmc(lik.lambda, p[colnames(w)], nsteps=5000, 
                       prior=prior, w=w, print.every=100)
head(mcmc.fit.full)
plot(mcmc.fit.full$i, mcmc.fit.full$p, type="l", xlab="generation",
     ylab="log(L)")
mcmc.fit.full<- mcmc.fit.full[501:4500,]
colMeans(mcmc.fit.full)[2:ncol(mcmc.fit.full)]
lambda<- mcmc.fit.full[,grep("lambda",colnames(mcmc.fit.full))]
mu<-mcmc.fit.full[,grep("mu",colnames(mcmc.fit.full))]
#write.csv(mcmc.fit.full[,grep("mu",colnames(mcmc.fit.full))],"mu.csv")
#mu<- read.csv("mu.csv",header = T,row.names = 1)

pdf("Musse_Phyrows.pdf")
colors<-setNames(c("#E62A8A","#7570B3", "#667532","#DBA62B"),1:4)
profiles.plot(lambda, col.line=colors, las=1, legend.pos="topright")
#profiles.plot(mu,col.line=colors, las=1, legend.pos="topright")
net.div<-lambda-mu
colnames(net.div)<-paste("lambda(",1:4,")",sep="")
profiles.plot(net.div, 
              xlab="Net diversification rate", ylab="Probability density",
              legend.pos="topleft",col.line=setNames(colors,colnames(net.div)),
              lty=1)
dev.off()

