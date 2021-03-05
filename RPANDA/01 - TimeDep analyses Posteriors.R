##
## Macroevolutionary perspectives to environmental change (Condamine et al. 2013 - Ecol. Lett.)
## Retrieved from Condamine et al. 2018 -Systematic Biology)
###Modified by Xu Zhang 2021/1/4 for Saussurea paleoenviroment dependent diversification anlysis
#Required R packages
library(RPANDA)
source("tables.summary.R")

#R codes and functions
source("fit_bd.R")
source("likelihood_bd.R")
source("Phi.R")
source("Psi.R")
source("integrate.R")


####################
## Fit the models ##
####################



#Sassurea<-read.nexus("Sassurea.trees")
posteriors<-sample(Sassurea, 200)
missing.lineages<-0.5
names<-c("Sassurea")

finalSaussurea_time<-list()
Saussurea_res<-c()

for (i in 1:length(posteriors))
{
	print(i)
	phyloi<-posteriors[[i]]
	tot_time<-max(node.age(phyloi)$ages)
	f<-Ntip(phyloi)/(Ntip(phyloi)+missing.lineages)
	cond="crown"
	
#####################################################
###### Time Independence (exponential variation) ####
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
###### Time Dependence (exponential variation) ######
#####################################################

	print(i)
	
# BTimeVar EXPO
print("BTimeVar EXPO")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){0}
lamb_par<-c(0.1,0.01)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=T; expo.mu=F; fix.mu=T

	treei_BTimeVar_EXPO<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVar_EXPO)
	

# BTimeVar DCST EXPO
print("BTimeVar DCST EXPO")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){y[1]}
lamb_par<-c(treei_BTimeVar_EXPO$lamb_par[1],treei_BTimeVar_EXPO$lamb_par[2])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=T; expo.mu=F; fix.mu=F

	treei_BTimeVarDCST_EXPO<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVarDCST_EXPO)


# BCST DTimeVar EXPO
print("BCST DTimeVar EXPO")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(treei_BCSTDCST$lamb_par[1])
mu_par<-c(0.01,0.01)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=T; fix.mu=F

	treei_BCSTDTimeVar_EXPO<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BCSTDTimeVar_EXPO)


# BTimeVar DTimeVar EXPO
print("BTimeVar DTimeVar EXPO")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(treei_BTimeVarDCST_EXPO$lamb_par[1],0.01)
mu_par<-c(0.05,0.01)
cst.lamb=F; cst.mu=F; expo.lamb=T; expo.mu=T; fix.mu=F

	treei_BTimeVarDTimeVar_EXPO<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVarDTimeVar_EXPO)


	################################################
	###### Time Dependence (linear variation) ######
	################################################
	
	print(i)
	
	# BTimeVar LIN
	print("BTimeVar LIN")
	f.lamb<-function(x,y){y[1]+(y[2]*x)}
	f.mu<-function(x,y){0}
	lamb_par<-c(0.01,0)
	mu_par<-c()
	cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
	
	treei_BTimeVar_LIN<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVar_LIN)
	
	
	# BTimeVar DCST LIN
	print("BTimeVar DCST LIN")
	f.lamb<-function(x,y){y[1]+(y[2]*x)}
	f.mu<-function(x,y){y[1]}
	lamb_par<-c(abs(treei_BTimeVar_LIN$lamb_par[1]),treei_BTimeVar_LIN$lamb_par[2])
	mu_par<-c(0.01)
	cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F
	
	treei_BTimeVarDCST_LIN<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVarDCST_LIN)
	
	
	# BCST DTimeVar LIN
	print("BCST DTimeVar LIN")
	f.lamb<-function(x,y){y[1]}
	f.mu<-function(x,y){y[1]+(y[2]*x)}
	lamb_par<-c(treei_BCSTDCST$lamb_par[1])
	mu_par<-c(0.01,0)
	cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
	
	treei_BCSTDTimeVar_LIN<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BCSTDTimeVar_LIN)
	
	
	# BTimeVar DTimeVar LIN
	print("BTimeVar DTimeVar LIN")
	f.lamb<-function(x,y){y[1]+(y[2]*x)}
	f.mu<-function(x,y){y[1]+(y[2]*x)}
	lamb_par<-c(abs(treei_BTimeVarDCST_LIN$lamb_par[1]),0)
	mu_par<-c(0.05,0)
	cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
	
	treei_BTimeVarDTimeVar_LIN<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVarDTimeVar_LIN)
	
###################################
###########Time_results############
###################################

	Time_results<-matrix(NA,7,8)
	colnames(Time_results)<-c("Models","df","logL","AICc","Lambda","AlphaTime","Mu","BetaTime")
	
	#Models
	Time_results[,1]<-c("BCSTDCST","BTimeVarDCST_EXPO","BCSTDTimeVar_EXPO","BTimeVarDTimeVar_EXPO","BTimeVarDCST_LIN","BCSTDTimeVar_LIN","BTimeVarDTimeVar_LIN")
	
	#Parameters
	Time_results[1,2]<-2
	Time_results[2,2]<-3
	Time_results[3,2]<-3
	Time_results[4,2]<-4
	Time_results[5,2]<-3
	Time_results[6,2]<-3
	Time_results[7,2]<-4

	
	#logL

	Time_results[1,3]<-round(treei_BCSTDCST$LH,3)
	Time_results[2,3]<-round(treei_BTimeVarDCST_EXPO$LH,3)
	Time_results[3,3]<-round(treei_BCSTDTimeVar_EXPO$LH,3)
	Time_results[4,3]<-round(treei_BTimeVarDTimeVar_EXPO$LH,3)
	Time_results[5,3]<-round(treei_BTimeVarDCST_LIN$LH,3)
	Time_results[6,3]<-round(treei_BCSTDTimeVar_LIN$LH,3)
	Time_results[7,3]<-round(treei_BTimeVarDTimeVar_LIN$LH,3)
	
	#AICc
	Time_results[1,4]<-round(treei_BCSTDCST$aicc,3)
	Time_results[2,4]<-round(treei_BTimeVarDCST_EXPO$aicc,3)
	Time_results[3,4]<-round(treei_BCSTDTimeVar_EXPO$aicc,3)
	Time_results[4,4]<-round(treei_BTimeVarDTimeVar_EXPO$aicc,3)
	Time_results[5,4]<-round(treei_BTimeVarDCST_LIN$aicc,3)
	Time_results[6,4]<-round(treei_BCSTDTimeVar_LIN$aicc,3)
	Time_results[7,4]<-round(treei_BTimeVarDTimeVar_LIN$aicc,3)
	
	
	#Lambda0

	Time_results[1,5]<-round(abs(treei_BCSTDCST$lamb_par[1]),3)
	Time_results[2,5]<-round(abs(treei_BTimeVarDCST_EXPO$lamb_par[1]),3)
	Time_results[3,5]<-round(abs(treei_BCSTDTimeVar_EXPO$lamb_par[1]),3)
	Time_results[4,5]<-round(abs(treei_BTimeVarDTimeVar_EXPO$lamb_par[1]),3)
	Time_results[5,5]<-round(abs(treei_BTimeVarDCST_LIN$lamb_par[1]),3)
	Time_results[6,5]<-round(abs(treei_BCSTDTimeVar_LIN$lamb_par[1]),3)
	Time_results[7,5]<-round(abs(treei_BTimeVarDTimeVar_LIN$lamb_par[1]),3)
	
	#Alpha Time

	Time_results[2,6]<-round(treei_BTimeVarDCST_EXPO$lamb_par[2],4)
	Time_results[4,6]<-round(treei_BTimeVarDTimeVar_EXPO$lamb_par[2],4)
	Time_results[5,6]<-round(treei_BTimeVarDCST_LIN$lamb_par[2],4)
	Time_results[7,6]<-round(treei_BTimeVarDTimeVar_LIN$lamb_par[2],4)
	
	#Mu0
	Time_results[1,7]<-round(abs(treei_BCSTDCST$mu_par[1]),3)
	Time_results[2,7]<-round(abs(treei_BTimeVarDCST_EXPO$mu_par[1]),3)
	Time_results[3,7]<-round(abs(treei_BCSTDTimeVar_EXPO$mu_par[1]),3)
	Time_results[4,7]<-round(abs(treei_BTimeVarDTimeVar_EXPO$mu_par[1]),3)
	Time_results[5,7]<-round(abs(treei_BTimeVarDCST_LIN$mu_par[1]),3)
	Time_results[6,7]<-round(abs(treei_BCSTDTimeVar_LIN$mu_par[1]),3)
	Time_results[7,7]<-round(abs(treei_BTimeVarDTimeVar_LIN$mu_par[1]),3)
	
	#Beta Time
	Time_results[3,8]<-round(treei_BCSTDTimeVar_EXPO$mu_par[2],4)
	Time_results[4,8]<-round(treei_BTimeVarDTimeVar_EXPO$mu_par[2],4)
	Time_results[6,8]<-round(treei_BCSTDTimeVar_LIN$mu_par[2],4)
	Time_results[7,8]<-round(treei_BTimeVarDTimeVar_LIN$mu_par[2],4)
	
	finalSaussurea_time[[i]]<-Time_results

	resi<-list("Clade"=names,"Clade_age"=tot_time,"Taxon_sampling"=Ntip(phyloi),"Clade_size"=Ntip(phyloi)+missing.lineages,"Sampling_fraction"=f,
	           "BCSTDCST"=treei_BCSTDCST, "BTimeVarDCST_EXPO"=treei_BTimeVarDCST_EXPO,
			   "BCSTDTimeVar_EXPO"=treei_BCSTDTimeVar_EXPO,"BTimeVarDTimeVar_EXPO"=treei_BTimeVarDTimeVar_EXPO,
	           "BTimeVarDCST_LIN"=treei_BTimeVarDCST_LIN,"BCSTDTimeVar_LIN"=treei_BCSTDTimeVar_LIN,
			   "BTimeVarDTimeVar_LIN"=treei_BTimeVarDTimeVar_LIN)
	
Saussurea_res<-c(Saussurea_res,list(resi))
}

final_table_Time<-tables.summary(finalSaussurea_time)


write.csv(final_table_Time,file="Time_results_Saussurea_TimeDep_posteriors.csv")
save(final_table_Time,file="final_table_Time_Saussurea_TimeDep_posteriors.Rdata")
save(finalSaussurea_time,file="finalSaussurea_time_TimeDep_posteriors.Rdata")
save(Saussurea_res,file="Saussurea_TimeDep_posteriors.Rdata")

