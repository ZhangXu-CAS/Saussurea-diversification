## Retrieved from Condamine et al. 2018 -Systematic Biology)
## Used for diversification rates dynamics of Saussurea -- Xu Zhang

library(TESS)
tree <- read.tree("Saussurea.tre")
treefile <- "Saussurea"
is.ultrametric(tree)
##False
library(phytools)
tree<- force.ultrametric(tree)
is.ultrametric(tree)
##True
analysis_name <- sprintf("CoMET_%s",treefile)
CDT <- "survival"
EXPECTED_NUM_EVENTS <- 2
MCMC_ITERATIONS <- 10000000
ALLOW_MASS_EXTINCTION <- TRUE
if ( ALLOW_MASS_EXTINCTION == TRUE ) {
analysis_name <- sprintf("CoMET_%s_ME",treefile)
} else {
analysis_name <- sprintf("CoMET_%s",treefile)
}
rho <- 0.5

#priorForms <- c("lognormal","normal","gamma")
priorForms <- c("lognormal")
tess.analysis(tree=tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS, 
				numExpectedMassExtinctions=EXPECTED_NUM_EVENTS, 
				initialSpeciationRate=2.0, initialExtinctionRate=exp(1.0), 
				empiricalHyperPriorInflation = 10.0,  
				samplingProbability=rho, 
				estimateMassExtinctionTimes = ALLOW_MASS_EXTINCTION, 
				estimateNumberMassExtinctions = ALLOW_MASS_EXTINCTION, 
				MAX_ITERATIONS = MCMC_ITERATIONS, THINNING = 100,  
				MAX_TIME = Inf, MIN_ESS = 1000, 
				CONDITION=CDT, dir = analysis_name)
out <- tess.process.output(analysis_name, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS)
NUM_FIGS <- 6
FIG_TYPES <- c("speciation rates","extinction rates","speciation shift times",
				"extinction shift times","mass extinction times","mass extinction Bayes factors")
if ( ALLOW_MASS_EXTINCTION == FALSE ) {
NUM_FIGS <- 6
FIG_TYPES <- c("speciation rates","extinction rates","speciation shift times",
		"extinction shift times","mass extinction times","mass extinction Bayes factors")
}
pdf(sprintf("%s.pdf",analysis_name))
layout.mat <- matrix(1:NUM_FIGS,nrow=2,ncol=NUM_FIGS / 2)
layout(layout.mat)
tess.plot.output(out,fig.types=FIG_TYPES,las=1)
dev.off()
treefile <- "Saussurea"
analysis_name <- sprintf("CoMET_%s",treefile)
CDT <- "survival"
EXPECTED_NUM_EVENTS <- 2
MCMC_ITERATIONS <- 5000000
ALLOW_MASS_EXTINCTION <- TRUE
if ( ALLOW_MASS_EXTINCTION == TRUE ) {
analysis_name <- sprintf("CoMET_%s_ME",treefile)
} else {
analysis_name <- sprintf("CoMET_%s",treefile)
}
tree <- read.nexus(file=sprintf("%s.tre",treefile) )
tree<- force.ultrametric(tree)
rho <- 0.5

#priorForms <- c("lognormal","normal","gamma")
priorForms <- c("lognormal")
tess.analysis(tree=tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS, 
				numExpectedMassExtinctions=EXPECTED_NUM_EVENTS, 
				initialSpeciationRate=2.0, initialExtinctionRate=1.0, 
				empiricalHyperPriorInflation = 10.0, 
				empiricalHyperPriorForm = priorForms, 
				samplingProbability=rho, 
				estimateMassExtinctionTimes = ALLOW_MASS_EXTINCTION, 
				estimateNumberMassExtinctions = ALLOW_MASS_EXTINCTION, 
				MAX_ITERATIONS = MCMC_ITERATIONS, THINNING = 100,  
				MAX_TIME = Inf, MIN_ESS = 1000, CONDITION=CDT, 
				dir = analysis_name,estimateNumberRateChanges = T,BURNIN = 5000)
out <- tess.process.output(analysis_name, tree, numExpectedRateChanges=EXPECTED_NUM_EVENTS)
NUM_FIGS <- 6
FIG_TYPES <- c("speciation rates","extinction rates",
			"speciation shift times","extinction shift times","mass extinction times","mass extinction Bayes factors")
if ( ALLOW_MASS_EXTINCTION == FALSE ) {
NUM_FIGS <- 6
FIG_TYPES <- c("speciation rates","extinction rates","speciation shift times",
				"extinction shift times","mass extinction times","mass extinction Bayes factors")
}
pdf(sprintf("%s.pdf",analysis_name))
layout.mat <- matrix(1:NUM_FIGS,nrow=2,ncol=NUM_FIGS / 2)
layout(layout.mat)
tess.plot.output(out,fig.types=FIG_TYPES,las=1)
dev.off()
