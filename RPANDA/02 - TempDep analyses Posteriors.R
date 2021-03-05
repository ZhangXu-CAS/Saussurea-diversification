
## Macroevolutionary perspectives to environmental change (Condamine et al. 2013 - Ecol. Lett.)
## Retrieved from Condamine et al. 2018 -Systematic Biology)
### Modified by Xu Zhang for Saussurea paleoenviroment dependent diversification anlysis

#Required R packages
library(RPANDA)
library(picante)
library(pspline)
source("tables.summary.R")

#R codes and functions
source("fit_bd.R")
source("fit_env_bd.R")
source("likelihood_bd.R")
source("Phi.R")
source("Psi.R")
source("integrate.R")


####################
## Fit the models ##
####################
env_data<-read.table("12Ma_Temperature.txt",header=T)

Sassurea<-read.nexus("Sassurea.trees")
posteriors<-sample(Sassurea, 500)
missing.lineages<-0.5
names<-c("Sassurea")

finalSaussurea_Temp<-list()
Saussurea_res_tem<-c()

for (i in 1:length(posteriors))
{
	print(i)
	phyloi<-posteriors[[i]]
	tot_time<-max(node.age(phyloi)$ages)
	f<-Ntip(phyloi)/(Ntip(phyloi)+missing.lineages)
	cond="crown"
	
#####################################################
###### Temp Independence (exponential variation) ####
#####################################################
# BCST (Pure birth)
print("BCST")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){0}
lamb_par<-c(0.1)
mu_par<-c()
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

	treei_BCST<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BCST)
	
# BCST DCST (constant Birth-death)
print("BCST DCST")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]}
lamb_par<-c(treei_BCST$lamb_par[1])
mu_par<-c(0.01)
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

	treei_BCSTDCST<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BCSTDCST)

	
	
#####################################################
###### Temp Dependence (exponential variation) ######
#####################################################

	print(i)
	
# BTempVar EXPO
print("BTempVar EXPO")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.1,0.01)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

	treei_BTempVar_EXPO<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTempVar_EXPO)


# BTempVar DCST EXPO
print("BTempVar DCST EXPO")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]}
lamb_par<-c(abs(treei_BTempVar_EXPO$lamb_par[1]),treei_BTempVar_EXPO$lamb_par[2])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

	treei_BTempVarDCST_EXPO<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTempVarDCST_EXPO)


# BCST DTempVar EXPO
print("BCST DTempVar EXPO")
f.lamb<-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(treei_BCSTDCST$lamb_par[1])
mu_par<-c(0.01,0.001)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

	treei_BCSTDTempVar_EXPO<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BCSTDTempVar_EXPO)


# BTempVar DTempVar EXPO
print("BTempVar DTempVar EXPO")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(abs(treei_BTempVarDCST_EXPO$lamb_par[1]),treei_BTempVarDCST_EXPO$lamb_par[2])
mu_par<-c(0.05,0.01)
cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

	treei_BTempVarDTempVar_EXPO<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTempVarDTempVar_EXPO)

################################################
###### Temp Dependence (linear variation) ######
################################################

	print(i)
	
# BTempVar LIN
print("BTempVar LIN")
f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.01,0)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

	treei_BTempVar_LIN<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTempVar_LIN)


# BTempVar DCST LIN
print("BTempVar DCST LIN")
f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){y[1]}
lamb_par<-c(abs(treei_BTempVar_LIN$lamb_par[1]),treei_BTempVar_LIN$lamb_par[2])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

	treei_BTempVarDCST_LIN<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTempVarDCST_LIN)


# BCST DTempVar LIN
print("BCST DTempVar LIN")
f.lamb<-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]+y[2]*x}
lamb_par<-c(treei_BCSTDCST$lamb_par[1])
#mu_par<-c(0.01,0.001)
mu_par<-c(0.01,0)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

	treei_BCSTDTempVar_LIN<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BCSTDTempVar_LIN)


# BTempVar DTempVar LIN
print("BTempVar DTempVar LIN")
f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){y[1]+y[2]*x}
lamb_par<-c(abs(treei_BTempVarDCST_LIN$lamb_par[1]),0)
mu_par<-c(0.05,0)
cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

	treei_BTempVarDTempVar_LIN<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTempVarDTempVar_LIN)

############# Temp_results ###########################################

	results<-matrix(NA,7,8)
	colnames(results)<-c("Models","Parameters","logL","AICc","Lambda","AlphaTemp","Mu","BetaTemp")

#Models
	results[,1]<-c("BCSTDCST","BTempVarDCST_EXPO","BCSTDTempVar_EXPO","BTempVarDTempVar_EXPO",
	"BTempVarDCST_LIN","BCSTDTempVar_LIN","BTempVarDTempVar_LIN")

#Parameters
	results[1,2]<-2
	results[2,2]<-3
	results[3,2]<-3
	results[4,2]<-4
	results[5,2]<-3
	results[6,2]<-3
	results[7,2]<-4


#logL
	results[1,3]<-round(treei_BCSTDCST$LH,3)
	results[2,3]<-round(treei_BTempVarDCST_EXPO$LH,3)
	results[3,3]<-round(treei_BCSTDTempVar_EXPO$LH,3)
	results[4,3]<-round(treei_BTempVarDTempVar_EXPO$LH,3)
	results[5,3]<-round(treei_BTempVarDCST_LIN$LH,3)
	results[6,3]<-round(treei_BCSTDTempVar_LIN$LH,3)
	results[7,3]<-round(treei_BTempVarDTempVar_LIN$LH,3)
				  
				  
#AICc             
	results[1,4]<-round(treei_BCSTDCST$aicc,3)
	results[2,4]<-round(treei_BTempVarDCST_EXPO$aicc,3)
	results[3,4]<-round(treei_BCSTDTempVar_EXPO$aicc,3)
	results[4,4]<-round(treei_BTempVarDTempVar_EXPO$aicc,3)
	results[5,4]<-round(treei_BTempVarDCST_LIN$aicc,3)
	results[6,4]<-round(treei_BCSTDTempVar_LIN$aicc,3)
	results[7,4]<-round(treei_BTempVarDTempVar_LIN$aicc,3)
				  
				  
#Lambda0          
	results[1,5]<-round((treei_BCSTDCST$lamb_par[1]),3)
	results[2,5]<-round((treei_BTempVarDCST_EXPO$lamb_par[1]),3)
	results[3,5]<-round((treei_BCSTDTempVar_EXPO$lamb_par[1]),3)
	results[4,5]<-round((treei_BTempVarDTempVar_EXPO$lamb_par[1]),3)
	results[5,5]<-round((treei_BTempVarDCST_LIN$lamb_par[1]),3)
	results[6,5]<-round((treei_BCSTDTempVar_LIN$lamb_par[1]),3)
	results[7,5]<-round((treei_BTempVarDTempVar_LIN$lamb_par[1]),3)
				 
				 
#Alpha Temp      
				 
	results[2,6]<-round(treei_BTempVarDCST_EXPO$lamb_par[2],4)
	results[4,6]<-round(treei_BTempVarDTempVar_EXPO$lamb_par[2],4)
	results[5,6]<-round(treei_BTempVarDCST_LIN$lamb_par[2],4)
	results[7,6]<-round(treei_BTempVarDTempVar_LIN$lamb_par[2],4)
				  
				  
#Mu0              
	results[1,7]<-round((treei_BCSTDCST$mu_par[1]),3)
	results[2,7]<-round((treei_BTempVarDCST_EXPO$mu_par[1]),3)
	results[3,7]<-round((treei_BCSTDTempVar_EXPO$mu_par[1]),3)
	results[4,7]<-round((treei_BTempVarDTempVar_EXPO$mu_par[1]),3)
	results[5,7]<-round((treei_BTempVarDCST_LIN$mu_par[1]),3)
	results[6,7]<-round((treei_BCSTDTempVar_LIN$mu_par[1]),3)
	results[7,7]<-round((treei_BTempVarDTempVar_LIN$mu_par[1]),3)
				  
				  
#Beta Temp        
	results[3,8]<-round(treei_BCSTDTempVar_EXPO$mu_par[2],4)
	results[4,8]<-round(treei_BTempVarDTempVar_EXPO$mu_par[2],4)
	results[6,8]<-round(treei_BCSTDTempVar_LIN$mu_par[2],4)
	results[7,8]<-round(treei_BTempVarDTempVar_LIN$mu_par[2],4)

	
	finalSaussurea_Temp[[i]]<-results
	
resi_tem<-list("Clade"=names,"Clade_age"=tot_time,"Taxon_sampling"=Ntip(phyloi),"Clade_size"=Ntip(phyloi)+missing.lineages,"Sampling_fraction"=f,
"BCSTDCST"=treei_BCSTDCST,"BTempVarDCST_EXPO"=treei_BTempVarDCST_EXPO,
"BCSTDTempVar_EXPO"=treei_BCSTDTempVar_EXPO,
"BTempVarDTempVar_EXPO"=treei_BTempVarDTempVar_EXPO,
"BTempVarDCST_LIN"=treei_BTempVarDCST_LIN,
"BCSTDTempVar_LIN"=treei_BCSTDTempVar_LIN,"BTempVarDTempVar_LIN"=treei_BTempVarDTempVar_LIN)

Saussurea_res_tem<-c(Saussurea_res_tem,list(resi_tem))
}

final_table_Temp<-tables.summary(finalSaussurea_Temp)

write.csv(final_table_Temp,file="Temp_results_Saussurea_TempDep_posteriors.csv")
save(final_table_Temp,file="final_table_TempSaussurea_TempDep_posteriors.Rdata")
save(finalSaussurea_Temp,file="finalSaussurea_Temp_TempDep_posteriors.Rdata")
save(Saussurea_res_tem,file="Saussurea_TempDep_posteriors.Rdata")
