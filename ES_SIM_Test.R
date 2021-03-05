#--------------------------------------------------
# Trait dependent tests using ES-SIM
#  Harvey and Rabosky 2017, Methods Ecol. Evol. 
#  https://doi.org/10.1111/2041-210X.12949
# 
# zx 2021/1/17
rm(list=ls())

library("ape")
library("mvtnorm")
library("phytools")
library("geiger")
library("tidyverse")

# import function
source("essim_DR.R")
###read tree
tree<- read.tree("Saussurea.tre")
###read ecology factors and rates
Eco_rates <- read.table("Eco_rates.txt",header = T)

row.names(Eco_rates) <- Eco_rates$Taxa

Eco_rates.vector <-  treedata(tree, Eco_rates)
###DR statistic rates
DR <- as.numeric(Eco_rates.vector$data[,6])
names(DR) <- rownames(Eco_rates.vector$data)
##climate pc1
pc1 <- as.numeric(Eco_rates.vector$data[,2])
names(pc1) <- rownames(Eco_rates.vector$data)

###climate pc2
pc2 <- as.numeric(Eco_rates.vector$data[,3])
names(pc2) <- rownames(Eco_rates.vector$data)

###niche breadth
NB <- as.numeric(Eco_rates.vector$data[,4])
names(NB) <- rownames(Eco_rates.vector$data)

rangsize <- as.numeric(Eco_rates.vector$data[,5])
names(rangsize) <- rownames(Eco_rates.vector$data)

essim_DR(Eco_rates.vector$phy, pc1, nsim= 3000, is = DR)
#rho   P Value 
#0.1695933 0.3593333  
essim_DR(Eco_rates.vector$phy, pc2, nsim= 3000, is = DR)
##      rho   P Value 
##0.09784959 0.64866667
essim_DR(Eco_rates.vector$phy, NB, nsim= 3000, is = DR)
##      rho   P Value 
### 0.36296127 0.02733333 
essim_DR(Eco_rates.vector$phy, rangsize, nsim= 3000, is = DR)
##   rho  P Value 
##0.39874943 0.0180000_clade

###default inverse equal splits statistic  
source("essim.R")
essim(Eco_rates.vector$phy, pc1, nsim= 3000)
##      rho   P Value 
## 0.1877708 0.3346667
essim(Eco_rates.vector$phy, pc2, nsim= 3000)
##      rho   P Value 
##0.09546193 0.63466667 
essim(Eco_rates.vector$phy, NB, nsim= 3000)
##      rho   P Value 
### 0.38663457 0.01933333 
essim(Eco_rates.vector$phy, rangsize, nsim= 3000)
##   rho  P Value 
##0.41113213 0.01133333

###graphy 
plot.new()
pdf("eco_rates_corr.pdf")
ggplot(Eco_rates, aes(Niche_breadth, DR_rates))+
  geom_point(aes(color = clade),size=2.0)+ 
  theme_minimal()+ 
  geom_smooth( method = lm,color = "darkgrey")
ggplot(Eco_rates, aes(Range_size, DR_rates))+
  geom_point(aes(color = clade),size=2.0)+ 
  theme_minimal()+ 
  geom_smooth(method = lm,color = "darkgrey")
dev.off()
