# Saussurea diversification
The project aims to gain insights into the drivers of radiating diversification in biodiversity hotspots using Saussurea (Asteraceae) as a case. 
## The manuscript is preprinted in bioRxiv:

Zhang X, Landis JB, Sun Y, Zhang H, Feng T, Lin N, Tiamiyu BB, Huang X, Deng T, Wang H, Sun H. 2021. Insights into the drivers of radiating diversification in biodiversity hotspots using Saussurea (Asteraceae) as a case. bioRxiv: 2021.2003.2015.435394. https://www.biorxiv.org/content/10.1101/2021.03.15.435394v1 

The R codes include:

a: RPANDA script, modified from the study of Condamine et al. (2018), fiting a series of time- and temperature-dependent likelihood diversification birth-death (BD) models;

b: BAMMtools, provided by Jacob B. Landis, used for diversification rate dynamics after BAMM analysis;

c: DR_statistic, a script to calculate DR statistics from Jetz et al. (2012);

d: ES_SIM_Test, trait dependent tests using ES-SIM (Harvey and Rabosky 2017), including essim.R and essim_DR.R;

e: FiSSE, used for binary traits dependent diversification;

f: GeoHiSSE, provided by Jacob B. Landis for  testing hypotheses about range-dependent diversification processes (Caetano et al., 2018); 

g: HiSSE, a script testing 25 different models in the HiSSE, BiSSE and null framework;

h: Musse, preparing to run MuSSE (Multi-State Speciation and Extinction) on a phylogenetic tree and multi-state character distribution; 

i: QuaSSE, preparing to run QuaSSE (Quantitative State Speciation and Extinction) on a phylogenetic tree and quantitative character distribution; 

j: Rgbif, calculating species area and niche breadth for each species, used for ecological niche driving diversification;

k: TESS_analysis, detecting the abrupt changes in speciation and extinction rates, Retrieved from Condamine et al. (2018);

l: traitDependent_functions

## References

Caetano DS, O'Meara BC, Beaulieu JM. 2018. Hidden state models improve state-dependent diversification approaches, including biogeographical models. Evolution 72(11): 2308-2324.

Condamine FL, Rolland J, Höhna S, Sperling FAH, Sanmartín I. 2018. Testing the Role of the Red Queen and Court Jester as Drivers of the Macroevolution of Apollo Butterflies. Systematic Biology 67(6): 940-964.

Harvey MG, Rabosky DL. 2018. Continuous traits and speciation rates: Alternatives to state-dependent diversification models. Methods in Ecology and Evolution 9(4): 984-993.

Jetz W, Thomas GH, Joy JB, Hartmann K, Mooers AO. 2012. The global diversity of birds in space and time. Nature 491(7424): 444-448.

Howard CC, Landis JB, Beaulieu JM, Cellinese N. 2020. Geophytism in monocots leads to higher rates of diversification. New Phytologist 225(2): 1023-1032.
