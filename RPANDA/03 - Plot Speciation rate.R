###Retrieved from Condamine et al. 2018 -Systematic Biology)
### Modified by Xu Zhang for Saussurea paleoenviroment dependent diversification anlysis
##########################
## Plot Speciation rate ##
##########################

library(picante)
library(pspline)

pdf("Saussurea diversification rate as a function of paleoenvironment.pdf")

par(mfrow=c(2,2), mar=c(4,4,2,2))

##
## a) Paleotemprature
##

InfTemp<-read.table("12Ma_Temperature.txt", header=T)

plot(-InfTemp[,"Age"],InfTemp[,"Temperature"],  xlab="", ylab="Temperatures", type="p", lwd=0.5, xlim=c(-12,0), col="grey", las=1, cex.axis=0.8, cex.main=1, bty = "n")
temp.spl<-smooth.spline(-InfTemp[,"Age"], InfTemp[,"Temperature"],df=100) #df=207.891
abline(v=c(-11.608,-5.33,-2.58),col="grey",lty="dotted",lwd="1")
lines(temp.spl, col="cornflowerblue", lwd=2)


legend(-12,16, bty="n",c("Late Mio."),text.col="black",cex=0.8)
legend(-7,16, bty="n",c("Pli."),text.col="black",cex=0.8)
legend(-4,16, bty="n",c("Ple."),text.col="black",cex=0.8)



##
## b) Speciation rate and Temperature
##

res<-sm.spline(InfTemp[,"Age"],InfTemp[,"Temperature"],df=20) #df=207.891
Temp_fun<-function(x){predict(res,x)}

parnass_TEMPDEPlamb_par1<-0.7585
parnass_TEMPDEPlamb_par2<--0.0933

f.lamb.mean<-function(x){parnass_TEMPDEPlamb_par1*exp(parnass_TEMPDEPlamb_par2*Temp_fun(x))}
plot(-InfTemp[,"Age"], f.lamb.mean(InfTemp[,"Age"]), ty="l",col="cornflowerblue",xlim=c(-12,0), ylim=c(0,1), lwd=2, yaxt="n", xlab="",ylab="Speciation rate (event/lineage/Myr)", cex.axis=0.8, cex.main=1, bty = "n")
axis(2, at = seq(0, 1, by = 0.2), las=1, cex.axis=0.8)
abline(v=c(-11.608,-5.33,-2.58),col="grey",lty="dotted",lwd="1")

f.lamb.low<-function(x){(parnass_TEMPDEPlamb_par1-0.05258)*exp((parnass_TEMPDEPlamb_par2-0.0708)*Temp_fun(x))}
lines(-InfTemp[,"Age"], f.lamb.low(InfTemp[,"Age"]),ty="l",col="cornflowerblue",xlim=c(-12,0),lwd=1,lty="dotted",yaxt="n")

f.lamb.high<-function(x){(parnass_TEMPDEPlamb_par1+0.0525)*exp((parnass_TEMPDEPlamb_par2+0.0708)*Temp_fun(x))}
lines(-InfTemp[,"Age"], f.lamb.high(InfTemp[,"Age"]),ty="l",col="cornflowerblue",xlim=c(-12,0),lwd=1,lty="dotted",yaxt="n")


dev.off()
