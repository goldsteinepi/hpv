#################
# HPV simulation study
# Citation: Goldstein ND, LeVasseur MT, Tran NK, Purtle J, Welles SL, Eppes SC. Modeling HPV vaccination scale-up among urban young men who have sex with men in the context of HIV. Manuscript in preparation.
# Note: Simulation datasets may be downloaded from: https://drive.google.com/file/d/1GwEIS989ah49DPK5zaWjFhNB4SoKmCPx/view?usp=sharing
# 1/24/17 -- Neal Goldstein
#################


### FUNCTIONS ###

prevention_paradigm = function(rand_vector, prevention_var, dataset, intervention_level=NA)
{
  if (prevention_var=="TAPVL") {
    #TAPVL (viral load under TAP, 1=suppressed, 2=not suppressed)
    starting_vals = ifelse(dataset$HIV_STATUS==3 & rand_vector<=0.42, 1, ifelse(dataset$HIV_STATUS==3 & rand_vector>0.42, 2, ifelse(dataset$HIV_STATUS==2 & dataset$HIV==1, 2, -1)))
  } else if (prevention_var=="PROBCOND") {
    #PROBCOND (condom probability)
    starting_vals = ifelse(dataset$RACE==1 & rand_vector<=0.23, 1, ifelse(dataset$RACE==1 & rand_vector<=0.51, 0, ifelse(dataset$RACE==1 & rand_vector<=0.71, 0.75, ifelse(dataset$RACE==1 & rand_vector<=0.86, 0.50, ifelse(dataset$RACE==1, 0.25, ifelse(dataset$RACE==2 & rand_vector<=0.27, 1, ifelse(dataset$RACE==2 & rand_vector<=0.47, 0, ifelse(dataset$RACE==2 & rand_vector<=0.69, 0.75, ifelse(dataset$RACE==2 & rand_vector<=0.86, 0.50, ifelse(dataset$RACE==2, 0.25, NA))))))))))
  } else if (prevention_var=="SEROTYPE") {
    #SEROTYPE (1=seropositioning/insertive, 2=seropositioning/receptive, 3=serosorting, 4=no seroadaptive behavior)
    starting_vals = ifelse(dataset$HIV_STATUS==3 & rand_vector<=0.133, 2, ifelse(dataset$HIV_STATUS==3 & rand_vector<=0.347, 3, ifelse(dataset$HIV_STATUS==3 & rand_vector>0.347, 4, ifelse(dataset$HIV_STATUS==2 & rand_vector<=0.053, 1, ifelse(dataset$HIV_STATUS==2 & rand_vector<=0.374, 3, ifelse(dataset$HIV_STATUS==2 & rand_vector>0.374, 4, ifelse(dataset$HIV_STATUS==1 & rand_vector<=0.095, 1, ifelse(dataset$HIV_STATUS==1 & rand_vector<=0.402, 3, ifelse(dataset$HIV_STATUS==1 & rand_vector>0.402, 4, NA)))))))))
  } else if (prevention_var=="ONPREP") {
    #ONPREP (0=not on prep, 1=on prep)
    starting_vals = ifelse(dataset$HIV_STATUS==1 & rand_vector<=(intervention_level/((1-0.19)-((1-0.19)*0.369))), 1, 0)
  } else if (prevention_var=="HPVDOSES") {
    #HPVDOSES (number of doses of HPV vaccine)
    starting_vals = ifelse(rand_vector<=(1-intervention_level), 0, ifelse(rand_vector<=(intervention_level*(4/13)+(1-intervention_level)), 1, ifelse(rand_vector<=(intervention_level*((4+2)/13)+(1-intervention_level)), 2, 3)))
  }
  
  return(starting_vals)
}

#network modeling
library("EpiModel") 


### STEP0: GLOBAL and NETWORK INITIALIZATION ###

#benchmarking
#time0 = Sys.time()

#ensure that all runs start equivalently
set.seed(777)

#number of simulations
nsims=50

#number of prevention paradigm scenarios evaluated per simulation
nparadigms = 6

#ensure that each simulation is unique
seed_val_population = sample(1:(nsims*10),nsims,replace=F)

#number of agents in model/initial YMSM population size in Philadelphia
nagents_start = 5329

#number of agents in model total including migration (value set later)
nagents_end = NA

#length of simulation (in days); change year multiplier to change length
length_sim = round(365 * 10, 0)

#DEBUG code
#nsims = nparadigms = 1
#length_sim = 365*1
#sims = day = current_data = 1

#initialize container for results of all simulations
abm_sims_50 = list(NA) #abm_sims is simulation data
long_sims_50 = list(NA) #long_sims is longitudinal summary data
network_sims_50 = list(NA) #network_sims is sexual network data

#initialize the network with number of agents
nw = network.initialize(n=nagents_start, directed=F)

#add assortative mixing characteristics
#race distribution: 1=White, 2=Non-white
nw = set.vertex.attribute(nw, attrname="race", value=ifelse(runif(nagents_start,0,1)<=0.34, 1, 2))

#define network parameters and estimate the network
est = suppressWarnings(netest(nw, formation=~edges+nodefactor("race")+nodematch("race"), target.stats=c(round(nagents_start*0.76/2), round(nagents_start*0.4/2), round(nagents_start*0.6/2/2)), coef.diss=dissolution_coefs(dissolution = ~offset(edges), duration=100)))
ns = netsim(est, param=param.net(inf.prob=0, act.rate=0, b.rate=(0.128/365), ds.rate=(0.128/365), di.rate=(0.128/365)), init=init.net(i.num=0), control=control.net(type="SI", nsteps=(length_sim-2), nsims=nsims, ncores=20, verbose=T))

#diagnostics
#dx = netdx(est, nsims=nsims, nsteps=(length_sim-2), ncores=20, keep.tedgelist=T, nwstats.formula=~edges+nodefactor("race")+nodematch("race"))
#summary(dx)
#plot(dx)

#cleanup
rm(nw,est)
gc()

#benchmarking
#time1 = Sys.time()


### SAVE INITIALIZATION: NETWORK ###

#save.image("simulation initialization network.RData")
#load("simulation initialization network.RData")

#write individual networkDynamic objects for optimization
for (sims in 1:nsims) {
  simnw = ns$network[[sims]]
  save(simnw, file=paste("simulation initialization network/sim", sims, ".RData", sep=""))
}

rm(simnw,sims,ns)
save.image("simulation initialization network/global.RData")


### STEP1-3: POPULATION INITIALIZATION ###

library("EpiModel") 
load("simulation initialization network/global.RData")

for (sims in 1:nsims)
{
  cat("\n\n************** ","Simulation: ",sims," **************\n",sep="")
  
  
  ### STEP1: EXTRACT SEXUAL NETWORK ###
  
  load(file=paste("simulation initialization network/sim", sims, ".RData", sep=""))
  
  #extract edge activity (i.e. sexual network)
  nw = get.edge.activity(simnw, as.spellList=T)
  nw$onset[nw$onset==-Inf] = 0
  nw$terminus[nw$terminus==Inf] = (length_sim-1)
  
  #extract vertex activity (i.e. migration)
  migration = get.vertex.activity(simnw, as.spellList=T)
  migration$onset[1:nagents_start] = 0 #all original agents start at day 0 to align to edgelist
  nagents_end = c(nagents_end, nrow(migration))
  
  #build sexual network dataset
  sex_network = data.frame(matrix(data=NA, nrow=nagents_end[sims+1], ncol=(6+length_sim)), stringsAsFactors=F)
  colnames(sex_network) = c("ID","ENTRY","EXIT","N_PARTNERS","RELATIONSHIP_DAYS","N_SEX_ACTS", paste("DAY",1:length_sim,sep=""))
  sex_network$ID = migration$vertex.id
  
  #add migration data
  sex_network$ENTRY = migration$onset + 1
  sex_network$EXIT = migration$terminus + 1
  
  for (i in 1:nrow(sex_network))
  {
    #check for partnerships
    partners = nw[nw$tail==i | nw$head==i, ]
    
    if (nrow(partners)>0) {
      #add each partner to daily network
      #note: concurrent relationships will be dissolved so that most recent partner overrides concurrency
      for (j in 1:nrow(partners))
      {
        ids = c(partners$tail[j], partners$head[j])
        sex_network[i, (7+partners$onset[j]):(7+partners$terminus[j])] = ids[!i==ids]
      }
    }
    
    #count partners, relationship days, sex acts (computed based on number of months in a relationship)
    sex_network$N_PARTNERS[i] = length(unique(na.omit(as.numeric(sex_network[i,7:ncol(sex_network)]))))
    sex_network$RELATIONSHIP_DAYS[i] = sum(!is.na(sex_network[i,7:ncol(sex_network)]))
    sex_network$N_SEX_ACTS[i] = round((rpois(1,80.6)/12 * sex_network$RELATIONSHIP_DAYS[i]/30))
  }
  rm(i,j,ids,partners,nw,migration)
  
  
  ### STEP2: INITIALIZING BASE POPULATION ###
  
  #counterfactual within each simulation: ensures that characteristics of the population are the same within each simulation
  set.seed(seed_val_population[sims])
  
  #generate an initial population of MSM, with same characteristics for all prevention scenarios
  baseline_pop = data.frame("ID"=sex_network$ID, "ACTIVE"=NA, "AGE"=NA, "RACE"=NA, "HIV"=NA, "HIV_VL"=NA, "HIV_STATUS"=NA, "HPV_ANAL"=NA, "HPV_ORAL"=NA, "HPV6_ANAL"=NA, "HPV6_ORAL"=NA, "HPV11_ANAL"=NA, "HPV11_ORAL"=NA, "HPV16_ANAL"=NA, "HPV16_ORAL"=NA, "HPV18_ANAL"=NA, "HPV18_ORAL"=NA, "HPV31_ANAL"=NA, "HPV31_ORAL"=NA, "HPV33_ANAL"=NA, "HPV33_ORAL"=NA, "HPV45_ANAL"=NA, "HPV45_ORAL"=NA, "HPV52_ANAL"=NA, "HPV52_ORAL"=NA, "HPV58_ANAL"=NA, "HPV58_ORAL"=NA, "HPV6_ANAL_CLEARED"=0, "HPV6_ORAL_CLEARED"=0, "HPV11_ANAL_CLEARED"=0, "HPV11_ORAL_CLEARED"=0, "HPV16_ANAL_CLEARED"=0, "HPV16_ORAL_CLEARED"=0, "HPV18_ANAL_CLEARED"=0, "HPV18_ORAL_CLEARED"=0, "HPV31_ANAL_CLEARED"=0, "HPV31_ORAL_CLEARED"=0, "HPV33_ANAL_CLEARED"=0, "HPV33_ORAL_CLEARED"=0, "HPV45_ANAL_CLEARED"=0, "HPV45_ORAL_CLEARED"=0, "HPV52_ANAL_CLEARED"=0, "HPV52_ORAL_CLEARED"=0, "HPV58_ANAL_CLEARED"=0, "HPV58_ORAL_CLEARED"=0, "CIRC"=NA, "SEX_ACT"=NA, "SEX_POSITION_ANAL_PREF"=NA, "SEX_POSITION_ORAL_PREF"=NA, "SEX_POSITION_ANAL"=NA, "SEX_POSITION_ORAL"=NA, "STI_TEST_DAY"=NA, "MAX_SEX"=NA, "SEX_COUNT"=0, "PARTNERS"=NA, "ORIGINAL_STATUS_HIV"=NA, "ORIGINAL_STATUS_HPV_ANAL"=NA, "ORIGINAL_STATUS_HPV_ORAL"=NA, "ORIGINAL_STATUS_HPV6_ANAL"=NA, "ORIGINAL_STATUS_HPV6_ORAL"=NA, "ORIGINAL_STATUS_HPV11_ANAL"=NA, "ORIGINAL_STATUS_HPV11_ORAL"=NA, "ORIGINAL_STATUS_HPV16_ANAL"=NA, "ORIGINAL_STATUS_HPV16_ORAL"=NA, "ORIGINAL_STATUS_HPV18_ANAL"=NA, "ORIGINAL_STATUS_HPV18_ORAL"=NA, "ORIGINAL_STATUS_HPV31_ANAL"=NA, "ORIGINAL_STATUS_HPV31_ORAL"=NA, "ORIGINAL_STATUS_HPV33_ANAL"=NA, "ORIGINAL_STATUS_HPV33_ORAL"=NA, "ORIGINAL_STATUS_HPV45_ANAL"=NA, "ORIGINAL_STATUS_HPV45_ORAL"=NA, "ORIGINAL_STATUS_HPV52_ANAL"=NA, "ORIGINAL_STATUS_HPV52_ORAL"=NA, "ORIGINAL_STATUS_HPV58_ANAL"=NA, "ORIGINAL_STATUS_HPV58_ORAL"=NA, "PREP_PREVENT"=0, "TAP_PREVENT"=0, "SERO_PREVENT"=0, "COND_PREVENT_ANAL"=0, "COND_PREVENT_ORAL"=0, "VAX_PREVENT_ANAL"=0, "VAX_PREVENT_ORAL"=0,"OVERALL_PREVENT_ANAL"=0,"OVERALL_PREVENT_ORAL"=0, "DISCORDANT"=0, "CAUSE_INFECT"=0, "DAYS_KNOWN_HIV_POSITIVE"=-1, "INCIDENCE_DAYS_HIV"=-1, "INCIDENCE_DAYS_HPV6_ANAL"=-1, "INCIDENCE_DAYS_HPV6_ORAL"=-1, "INCIDENCE_DAYS_HPV11_ANAL"=-1, "INCIDENCE_DAYS_HPV11_ORAL"=-1, "INCIDENCE_DAYS_HPV16_ANAL"=-1, "INCIDENCE_DAYS_HPV16_ORAL"=-1, "INCIDENCE_DAYS_HPV18_ANAL"=-1, "INCIDENCE_DAYS_HPV18_ORAL"=-1, "INCIDENCE_DAYS_HPV31_ANAL"=-1, "INCIDENCE_DAYS_HPV31_ORAL"=-1, "INCIDENCE_DAYS_HPV33_ANAL"=-1, "INCIDENCE_DAYS_HPV33_ORAL"=-1, "INCIDENCE_DAYS_HPV45_ANAL"=-1, "INCIDENCE_DAYS_HPV45_ORAL"=-1, "INCIDENCE_DAYS_HPV52_ANAL"=-1, "INCIDENCE_DAYS_HPV52_ORAL"=-1, "INCIDENCE_DAYS_HPV58_ANAL"=-1, "INCIDENCE_DAYS_HPV58_ORAL"=-1, stringsAsFactors=F)
  
  #age distribution: 1=<=18-19, 2=>=20-26
  baseline_pop$AGE = ifelse(runif(nagents_end[sims+1],0,1)<=0.20, 1, 2)
  
  #race distribution: 1=White, 2=Non-white
  #note: inherit from epimodel based on assortative mixing characteristic
  baseline_pop$RACE = get.vertex.attribute(simnw, "race")
  
  #infect initial population with HIV contingent upon race
  baseline_pop$HIV = ifelse(baseline_pop$AGE==1 & baseline_pop$RACE==1 & runif(nagents_end[sims+1],0,1)<=0.00497, 1, ifelse(baseline_pop$AGE==2 & baseline_pop$RACE==1 & runif(nagents_end[sims+1],0,1)<=0.12069, 1, ifelse(baseline_pop$AGE==1 & baseline_pop$RACE==2 & runif(nagents_end[sims+1],0,1)<=0.048387, 1, ifelse(baseline_pop$AGE==2 & baseline_pop$RACE==2 & runif(nagents_end[sims+1],0,1)<=0.224543, 1, 0))))
  
  #HIV viral load category: 1=chronic, 2=acute
  baseline_pop$HIV_VL = ifelse(baseline_pop$HIV==1 & runif(nagents_end[sims+1],0,1)<=0.025, 2, ifelse(baseline_pop$HIV==1, 1, -1))
  
  #knowledge of HIV status: 1=HIV negative and aware, 2=unknown, 3=HIV positive and aware
  baseline_pop$HIV_STATUS = ifelse(baseline_pop$HIV==0 & runif(nagents_end[sims+1],0,1)<=0.369, 2, ifelse(baseline_pop$HIV==0, 1, NA))
  baseline_pop$HIV_STATUS = ifelse(is.na(baseline_pop$HIV_STATUS) & baseline_pop$HIV==1 & runif(nagents_end[sims+1],0,1)<=0.440, 2, ifelse(baseline_pop$HIV==1, 3, baseline_pop$HIV_STATUS))
  
  #infect intial population with HPV contingent upon HIV
  baseline_pop$HPV6_ANAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.329, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.154, 1, 0))
  baseline_pop$HPV6_ORAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.034, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.008, 1, 0))
  baseline_pop$HPV11_ANAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.182, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.074, 1, 0))
  baseline_pop$HPV11_ORAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.011, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0, 1, 0))
  baseline_pop$HPV16_ANAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.364, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.125, 1, 0))
  baseline_pop$HPV16_ORAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.023, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.011, 1, 0))
  baseline_pop$HPV18_ANAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.329, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.073, 1, 0))
  baseline_pop$HPV18_ORAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.034, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.005, 1, 0))
  baseline_pop$HPV31_ANAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.141, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.036, 1, 0))
  baseline_pop$HPV31_ORAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.025, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.002, 1, 0))
  baseline_pop$HPV33_ANAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.069, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.064, 1, 0))
  baseline_pop$HPV33_ORAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.051, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.010, 1, 0))
  baseline_pop$HPV45_ANAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.152, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.121, 1, 0))
  baseline_pop$HPV45_ORAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.018, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.002, 1, 0))
  baseline_pop$HPV52_ANAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.242, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.090, 1, 0))
  baseline_pop$HPV52_ORAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.018, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.010, 1, 0))
  baseline_pop$HPV58_ANAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.211, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0.098, 1, 0))
  baseline_pop$HPV58_ORAL = ifelse(baseline_pop$HIV==1, ifelse(runif(nagents_end[sims+1],0,1)<=0.007, 1, 0), ifelse(runif(nagents_end[sims+1],0,1)<=0, 1, 0))
  baseline_pop$HPV_ANAL = ifelse(baseline_pop$HPV6_ANAL==1 | baseline_pop$HPV11_ANAL==1 | baseline_pop$HPV16_ANAL==1 | baseline_pop$HPV18_ANAL==1 | baseline_pop$HPV31_ANAL==1 | baseline_pop$HPV33_ANAL==1 | baseline_pop$HPV45_ANAL==1 | baseline_pop$HPV52_ANAL==1 | baseline_pop$HPV58_ANAL==1, 1, 0)
  baseline_pop$HPV_ORAL = ifelse(baseline_pop$HPV6_ORAL==1 | baseline_pop$HPV11_ORAL==1 | baseline_pop$HPV16_ORAL==1 | baseline_pop$HPV18_ORAL==1 | baseline_pop$HPV31_ORAL==1 | baseline_pop$HPV33_ORAL==1 | baseline_pop$HPV45_ORAL==1 | baseline_pop$HPV52_ORAL==1 | baseline_pop$HPV58_ORAL==1, 1, 0)
  
  #set proportion of prevalent cases as incident
  baseline_pop$INCIDENCE_DAYS_HIV = ifelse(baseline_pop$HIV_VL==2, round(runif(sum(baseline_pop$HIV_VL==2),0,63)), baseline_pop$INCIDENCE_DAYS_HIV)
  baseline_pop$INCIDENCE_DAYS_HPV6_ANAL = ifelse(baseline_pop$HPV6_ANAL==1 & runif(nagents_end[sims+1],0,1)<=0.19, round(runif(sum(baseline_pop$HPV6_ANAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV6_ANAL)
  baseline_pop$INCIDENCE_DAYS_HPV6_ORAL = ifelse(baseline_pop$HPV6_ORAL==1 & runif(nagents_end[sims+1],0,1)<=0.16, round(runif(sum(baseline_pop$HPV6_ORAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV6_ORAL)
  baseline_pop$INCIDENCE_DAYS_HPV11_ANAL = ifelse(baseline_pop$HPV11_ANAL==1 & runif(nagents_end[sims+1],0,1)<=0.19, round(runif(sum(baseline_pop$HPV11_ANAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV11_ANAL)
  baseline_pop$INCIDENCE_DAYS_HPV11_ORAL = ifelse(baseline_pop$HPV11_ORAL==1 & runif(nagents_end[sims+1],0,1)<=0.16, round(runif(sum(baseline_pop$HPV11_ORAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV11_ORAL)
  baseline_pop$INCIDENCE_DAYS_HPV16_ANAL = ifelse(baseline_pop$HPV16_ANAL==1 & runif(nagents_end[sims+1],0,1)<=0.19, round(runif(sum(baseline_pop$HPV16_ANAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV16_ANAL)
  baseline_pop$INCIDENCE_DAYS_HPV16_ORAL = ifelse(baseline_pop$HPV16_ORAL==1 & runif(nagents_end[sims+1],0,1)<=0.16, round(runif(sum(baseline_pop$HPV16_ORAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV16_ORAL)
  baseline_pop$INCIDENCE_DAYS_HPV18_ANAL = ifelse(baseline_pop$HPV18_ANAL==1 & runif(nagents_end[sims+1],0,1)<=0.19, round(runif(sum(baseline_pop$HPV18_ANAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV18_ANAL)
  baseline_pop$INCIDENCE_DAYS_HPV18_ORAL = ifelse(baseline_pop$HPV18_ORAL==1 & runif(nagents_end[sims+1],0,1)<=0.16, round(runif(sum(baseline_pop$HPV18_ORAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV18_ORAL)
  baseline_pop$INCIDENCE_DAYS_HPV31_ANAL = ifelse(baseline_pop$HPV31_ANAL==1 & runif(nagents_end[sims+1],0,1)<=0.19, round(runif(sum(baseline_pop$HPV31_ANAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV31_ANAL)
  baseline_pop$INCIDENCE_DAYS_HPV31_ORAL = ifelse(baseline_pop$HPV31_ORAL==1 & runif(nagents_end[sims+1],0,1)<=0.16, round(runif(sum(baseline_pop$HPV31_ORAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV31_ORAL)
  baseline_pop$INCIDENCE_DAYS_HPV33_ANAL = ifelse(baseline_pop$HPV33_ANAL==1 & runif(nagents_end[sims+1],0,1)<=0.19, round(runif(sum(baseline_pop$HPV33_ANAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV33_ANAL)
  baseline_pop$INCIDENCE_DAYS_HPV33_ORAL = ifelse(baseline_pop$HPV33_ORAL==1 & runif(nagents_end[sims+1],0,1)<=0.16, round(runif(sum(baseline_pop$HPV33_ORAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV33_ORAL)
  baseline_pop$INCIDENCE_DAYS_HPV45_ANAL = ifelse(baseline_pop$HPV45_ANAL==1 & runif(nagents_end[sims+1],0,1)<=0.19, round(runif(sum(baseline_pop$HPV45_ANAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV45_ANAL)
  baseline_pop$INCIDENCE_DAYS_HPV45_ORAL = ifelse(baseline_pop$HPV45_ORAL==1 & runif(nagents_end[sims+1],0,1)<=0.16, round(runif(sum(baseline_pop$HPV45_ORAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV45_ORAL)
  baseline_pop$INCIDENCE_DAYS_HPV52_ANAL = ifelse(baseline_pop$HPV52_ANAL==1 & runif(nagents_end[sims+1],0,1)<=0.19, round(runif(sum(baseline_pop$HPV52_ANAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV52_ANAL)
  baseline_pop$INCIDENCE_DAYS_HPV52_ORAL = ifelse(baseline_pop$HPV52_ORAL==1 & runif(nagents_end[sims+1],0,1)<=0.16, round(runif(sum(baseline_pop$HPV52_ORAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV52_ORAL)
  baseline_pop$INCIDENCE_DAYS_HPV58_ANAL = ifelse(baseline_pop$HPV58_ANAL==1 & runif(nagents_end[sims+1],0,1)<=0.19, round(runif(sum(baseline_pop$HPV58_ANAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV58_ANAL)
  baseline_pop$INCIDENCE_DAYS_HPV58_ORAL = ifelse(baseline_pop$HPV58_ORAL==1 & runif(nagents_end[sims+1],0,1)<=0.16, round(runif(sum(baseline_pop$HPV58_ORAL),0,228)), baseline_pop$INCIDENCE_DAYS_HPV58_ORAL)  
  
  #circumcision status
  #baseline_pop$CIRC = ifelse(runif(nagents_end[sims+1],0,1)<=0.762, 1, 0)
  baseline_pop$CIRC = 0 #set everyone as uncircumcised as inconclusive evidence to benefit in US population
  
  #sex act preference, 1=anal only, 2=oral only, 3=both
  randact = runif(nagents_end[sims+1],0,1)
  baseline_pop$SEX_ACT = ifelse(randact<=0.22, 1, ifelse(randact<=0.66, 2, 3))
  rm(randact)
  
  #sex_position preference, 1=insertive, 2=receptive, 3=both
  analpref = runif(nagents_end[sims+1],0,1)
  oralpref = runif(nagents_end[sims+1],0,1)
  baseline_pop$SEX_POSITION_ANAL_PREF = ifelse(analpref<=0.24, 1, ifelse(analpref<=0.56, 2, 3))
  baseline_pop$SEX_POSITION_ORAL_PREF = ifelse(oralpref<=0.16, 1, ifelse(oralpref<=0.30, 2, 3))
  rm(analpref,oralpref)
  
  #sex position (scale of dominance), if preference cannot be honored
  baseline_pop$SEX_POSITION_ANAL = runif(nagents_end[sims+1],0,1)
  baseline_pop$SEX_POSITION_ORAL = runif(nagents_end[sims+1],0,1)
  
  #day of year will be tested for HIV (and other STIs), based on a future test probability of 1=tested 1/3 of year, 2=tested 2/3 of year, 3=tested 3/3 of year, 4=not tested
  testprob = runif(nagents_end[sims+1],0,1)
  futuretest = ifelse(baseline_pop$HIV_STATUS==3, 4, ifelse(testprob<=0.19, 1, ifelse(testprob<=0.30, 2, ifelse(testprob<=0.70, 3, 4))))
  probdaytest = runif(nagents_end[sims+1],0,1)
  baseline_pop$STI_TEST_DAY = ifelse(futuretest==1, round(1+(110-1)*probdaytest), ifelse(futuretest==2, round(111+(220-111)*probdaytest), ifelse(futuretest==3, round(221+(330-221)*probdaytest), NA)))
  rm(testprob,probdaytest,futuretest)
  
  #original status of infections
  baseline_pop$ORIGINAL_STATUS_HIV = baseline_pop$HIV
  baseline_pop$ORIGINAL_STATUS_HPV6_ANAL = baseline_pop$HPV6_ANAL
  baseline_pop$ORIGINAL_STATUS_HPV6_ORAL = baseline_pop$HPV6_ORAL
  baseline_pop$ORIGINAL_STATUS_HPV11_ANAL = baseline_pop$HPV11_ANAL
  baseline_pop$ORIGINAL_STATUS_HPV11_ORAL = baseline_pop$HPV11_ORAL
  baseline_pop$ORIGINAL_STATUS_HPV16_ANAL = baseline_pop$HPV16_ANAL
  baseline_pop$ORIGINAL_STATUS_HPV16_ORAL = baseline_pop$HPV16_ORAL
  baseline_pop$ORIGINAL_STATUS_HPV18_ANAL = baseline_pop$HPV18_ANAL
  baseline_pop$ORIGINAL_STATUS_HPV18_ORAL = baseline_pop$HPV18_ORAL
  baseline_pop$ORIGINAL_STATUS_HPV31_ANAL = baseline_pop$HPV31_ANAL
  baseline_pop$ORIGINAL_STATUS_HPV31_ORAL = baseline_pop$HPV31_ORAL
  baseline_pop$ORIGINAL_STATUS_HPV33_ANAL = baseline_pop$HPV33_ANAL
  baseline_pop$ORIGINAL_STATUS_HPV33_ORAL = baseline_pop$HPV33_ORAL
  baseline_pop$ORIGINAL_STATUS_HPV45_ANAL = baseline_pop$HPV45_ANAL
  baseline_pop$ORIGINAL_STATUS_HPV45_ORAL = baseline_pop$HPV45_ORAL
  baseline_pop$ORIGINAL_STATUS_HPV52_ANAL = baseline_pop$HPV52_ANAL
  baseline_pop$ORIGINAL_STATUS_HPV52_ORAL = baseline_pop$HPV52_ORAL
  baseline_pop$ORIGINAL_STATUS_HPV58_ANAL = baseline_pop$HPV58_ANAL
  baseline_pop$ORIGINAL_STATUS_HPV58_ORAL = baseline_pop$HPV58_ORAL
  baseline_pop$ORIGINAL_STATUS_HPV_ANAL = baseline_pop$HPV_ANAL
  baseline_pop$ORIGINAL_STATUS_HPV_ORAL = baseline_pop$HPV_ORAL
  
  #sex and partner distribution
  baseline_pop$MAX_SEX = sex_network$N_SEX_ACTS
  baseline_pop$PARTNERS = sex_network$N_PARTNERS
  
  
  ### STEP3: ASSIGNMENT OF PREVENTION PARADIGMS ###
  
  #baseline population
  randprep = runif(nagents_end[sims+1],0,1)
  randtap = runif(nagents_end[sims+1],0,1)
  randcond = runif(nagents_end[sims+1],0,1)
  randsero = runif(nagents_end[sims+1],0,1)
  randvax =  runif(nagents_end[sims+1],0,1)
  
  #create prevention paradigm datasets: abcde where a=prep, b=treatment as prevention, c=condom use, d=seroadaption, e=vaccination
  paradigm11111 = baseline_pop
  
  #PrEP, TasP, condoms, seroadaption, vaccination
  paradigm11111$PREP = 1 #turn off for calibration to historic data
  paradigm11111$TAP = 1
  paradigm11111$COND = 1
  paradigm11111$SERO = 1
  paradigm11111$VAX = 1
  paradigm11111$ONPREP = prevention_paradigm(randprep,"ONPREP",paradigm11111,0.122)
  paradigm11111$TAPVL = prevention_paradigm(randtap,"TAPVL",paradigm11111)
  paradigm11111$PROBCOND = prevention_paradigm(randcond,"PROBCOND",paradigm11111)
  paradigm11111$SEROTYPE = prevention_paradigm(randsero,"SEROTYPE",paradigm11111)
  
  #assign into lists by vaccination coverage
  abm_sims = list("Paradigm11111Vax0"=cbind(paradigm11111,"HPVDOSES"=prevention_paradigm(randvax,"HPVDOSES",paradigm11111,0)),
                  "Paradigm11111Vax13"=cbind(paradigm11111,"HPVDOSES"=prevention_paradigm(randvax,"HPVDOSES",paradigm11111,0.13)),
                  "Paradigm11111Vax25"=cbind(paradigm11111,"HPVDOSES"=prevention_paradigm(randvax,"HPVDOSES",paradigm11111,0.25)),
                  "Paradigm11111Vax50"=cbind(paradigm11111,"HPVDOSES"=prevention_paradigm(randvax,"HPVDOSES",paradigm11111,0.50)),
                  "Paradigm11111Vax80"=cbind(paradigm11111,"HPVDOSES"=prevention_paradigm(randvax,"HPVDOSES",paradigm11111,0.80)),
                  "Paradigm11111Vax100"=cbind(paradigm11111,"HPVDOSES"=prevention_paradigm(randvax,"HPVDOSES",paradigm11111,1)))
  
  #determine vaccine immunogenicity
  immune = runif(nagents_end[sims+1],0,1)<=0.594
  abm_sims$Paradigm11111Vax0$IMMUNE = ifelse(abm_sims$Paradigm11111Vax0$HPVDOSES>=1 & immune, 1, 0)
  abm_sims$Paradigm11111Vax13$IMMUNE = ifelse(abm_sims$Paradigm11111Vax13$HPVDOSES>=1 & immune, 1, 0)
  abm_sims$Paradigm11111Vax25$IMMUNE = ifelse(abm_sims$Paradigm11111Vax25$HPVDOSES>=1 & immune, 1, 0)
  abm_sims$Paradigm11111Vax50$IMMUNE = ifelse(abm_sims$Paradigm11111Vax50$HPVDOSES>=1 & immune, 1, 0)
  abm_sims$Paradigm11111Vax80$IMMUNE = ifelse(abm_sims$Paradigm11111Vax80$HPVDOSES>=1 & immune, 1, 0)
  abm_sims$Paradigm11111Vax100$IMMUNE = ifelse(abm_sims$Paradigm11111Vax100$HPVDOSES>=1 & immune, 1, 0)
  rm(immune)
  
  #initialize a list for longitudinal tracking of daily outcomes; must match length of abm_sims
  long_sims = list("Paradigm11111Vax0"=data.frame("DAY"=1:length_sim, "HIV"=NA, "HPV6_ANAL"=NA, "HPV6_ORAL"=NA, "HPV11_ANAL"=NA, "HPV11_ORAL"=NA, "HPV16_ANAL"=NA, "HPV16_ORAL"=NA, "HPV18_ANAL"=NA, "HPV18_ORAL"=NA, "HPV31_ANAL"=NA, "HPV31_ORAL"=NA, "HPV33_ANAL"=NA, "HPV33_ORAL"=NA, "HPV45_ANAL"=NA, "HPV45_ORAL"=NA, "HPV52_ANAL"=NA, "HPV52_ORAL"=NA, "HPV58_ANAL"=NA, "HPV58_ORAL"=NA, stringsAsFactors=F),
                   "Paradigm11111Vax13"=data.frame("DAY"=1:length_sim, "HIV"=NA, "HPV6_ANAL"=NA, "HPV6_ORAL"=NA, "HPV11_ANAL"=NA, "HPV11_ORAL"=NA, "HPV16_ANAL"=NA, "HPV16_ORAL"=NA, "HPV18_ANAL"=NA, "HPV18_ORAL"=NA, "HPV31_ANAL"=NA, "HPV31_ORAL"=NA, "HPV33_ANAL"=NA, "HPV33_ORAL"=NA, "HPV45_ANAL"=NA, "HPV45_ORAL"=NA, "HPV52_ANAL"=NA, "HPV52_ORAL"=NA, "HPV58_ANAL"=NA, "HPV58_ORAL"=NA, stringsAsFactors=F),
                   "Paradigm11111Vax25"=data.frame("DAY"=1:length_sim, "HIV"=NA, "HPV6_ANAL"=NA, "HPV6_ORAL"=NA, "HPV11_ANAL"=NA, "HPV11_ORAL"=NA, "HPV16_ANAL"=NA, "HPV16_ORAL"=NA, "HPV18_ANAL"=NA, "HPV18_ORAL"=NA, "HPV31_ANAL"=NA, "HPV31_ORAL"=NA, "HPV33_ANAL"=NA, "HPV33_ORAL"=NA, "HPV45_ANAL"=NA, "HPV45_ORAL"=NA, "HPV52_ANAL"=NA, "HPV52_ORAL"=NA, "HPV58_ANAL"=NA, "HPV58_ORAL"=NA, stringsAsFactors=F),
                   "Paradigm11111Vax50"=data.frame("DAY"=1:length_sim, "HIV"=NA, "HPV6_ANAL"=NA, "HPV6_ORAL"=NA, "HPV11_ANAL"=NA, "HPV11_ORAL"=NA, "HPV16_ANAL"=NA, "HPV16_ORAL"=NA, "HPV18_ANAL"=NA, "HPV18_ORAL"=NA, "HPV31_ANAL"=NA, "HPV31_ORAL"=NA, "HPV33_ANAL"=NA, "HPV33_ORAL"=NA, "HPV45_ANAL"=NA, "HPV45_ORAL"=NA, "HPV52_ANAL"=NA, "HPV52_ORAL"=NA, "HPV58_ANAL"=NA, "HPV58_ORAL"=NA, stringsAsFactors=F),
                   "Paradigm11111Vax80"=data.frame("DAY"=1:length_sim, "HIV"=NA, "HPV6_ANAL"=NA, "HPV6_ORAL"=NA, "HPV11_ANAL"=NA, "HPV11_ORAL"=NA, "HPV16_ANAL"=NA, "HPV16_ORAL"=NA, "HPV18_ANAL"=NA, "HPV18_ORAL"=NA, "HPV31_ANAL"=NA, "HPV31_ORAL"=NA, "HPV33_ANAL"=NA, "HPV33_ORAL"=NA, "HPV45_ANAL"=NA, "HPV45_ORAL"=NA, "HPV52_ANAL"=NA, "HPV52_ORAL"=NA, "HPV58_ANAL"=NA, "HPV58_ORAL"=NA, stringsAsFactors=F),
                   "Paradigm11111Vax100"=data.frame("DAY"=1:length_sim, "HIV"=NA, "HPV6_ANAL"=NA, "HPV6_ORAL"=NA, "HPV11_ANAL"=NA, "HPV11_ORAL"=NA, "HPV16_ANAL"=NA, "HPV16_ORAL"=NA, "HPV18_ANAL"=NA, "HPV18_ORAL"=NA, "HPV31_ANAL"=NA, "HPV31_ORAL"=NA, "HPV33_ANAL"=NA, "HPV33_ORAL"=NA, "HPV45_ANAL"=NA, "HPV45_ORAL"=NA, "HPV52_ANAL"=NA, "HPV52_ORAL"=NA, "HPV58_ANAL"=NA, "HPV58_ORAL"=NA, stringsAsFactors=F))
  
  #cleanup
  rm(randtap,randcond,randsero,randprep,randvax,baseline_pop,paradigm11111,simnw)
  gc()
  
  #add this initialization to the lists
  names(abm_sims) = paste("simulation",sims,names(abm_sims),sep="_")
  abm_sims_50 = c(abm_sims_50,abm_sims)
  names(long_sims) = paste("simulation",sims,names(long_sims),sep="_")
  long_sims_50 = c(long_sims_50,long_sims)
  sex_network = list(sex_network)
  names(sex_network) = paste("simulation",sims,sep="_")
  network_sims_50 = c(network_sims_50,sex_network)
}
rm(abm_sims,long_sims,sex_network,seed_val_population,prevention_paradigm,sims)
abm_sims_50[[1]] = NULL
long_sims_50[[1]] = NULL
network_sims_50[[1]] = NULL
nagents_end = nagents_end[-1]
gc()

#benchmarking
#time2 = Sys.time()

#migration: per time step, % of population expected to enter/exit
#sum(network_sims_50[[1]]$ENTRY==1 & network_sims_50[[1]]$EXIT==length_sim)
#sum(network_sims_50[[1]]$ENTRY!=1 | network_sims_50[[1]]$EXIT!=length_sim)
#sum(network_sims_50[[1]]$ENTRY!=1 | network_sims_50[[1]]$EXIT!=length_sim)/length_sim
#sum(network_sims_50[[1]]$ENTRY!=1 | network_sims_50[[1]]$EXIT!=length_sim)/length_sim/nagents_start


### SAVE INITIALIZATION: POPULATION ###

save.image("simulation initialization population.RData")
load("simulation initialization population.RData")


### STEP4: SIMULATE ###

# #sensitivity analysis for improved vaccination efficacy
# for (sims in 0:(nsims-1))
# {
#   
#   immune = runif(nagents_end[sims+1],0,1)<=0.949
#   abm_sims_50[[sims*6+1]]$IMMUNE = ifelse(abm_sims_50[[sims*6+1]]$HPVDOSES>=1 & immune, 1, 0)
#   abm_sims_50[[sims*6+2]]$IMMUNE = ifelse(abm_sims_50[[sims*6+2]]$HPVDOSES>=1 & immune, 1, 0)
#   abm_sims_50[[sims*6+3]]$IMMUNE = ifelse(abm_sims_50[[sims*6+3]]$HPVDOSES>=1 & immune, 1, 0)
#   abm_sims_50[[sims*6+4]]$IMMUNE = ifelse(abm_sims_50[[sims*6+4]]$HPVDOSES>=1 & immune, 1, 0)
#   abm_sims_50[[sims*6+5]]$IMMUNE = ifelse(abm_sims_50[[sims*6+5]]$HPVDOSES>=1 & immune, 1, 0)
#   abm_sims_50[[sims*6+6]]$IMMUNE = ifelse(abm_sims_50[[sims*6+6]]$HPVDOSES>=1 & immune, 1, 0)
#   rm(immune)
#   
# }
# rm(sims)

#ensure that all runs start equivalently
set.seed(777)

for (sims in 1:nsims)
{
  #load previously initialized objects
  abm_sims = abm_sims_50[((sims-1)*nparadigms+1):(((sims-1)*nparadigms)+nparadigms)]
  long_sims = long_sims_50[((sims-1)*nparadigms+1):(((sims-1)*nparadigms)+nparadigms)]
  sex_network = network_sims_50[[sims]]
  
  #ensure that each day is unique
  seed_val_day = sample(1:(length_sim*10),length_sim,replace=F)
  
  for (day in 1:length_sim)
  {
    cat("\n\n************** ","Simulation: ",sims," *** Day: ",day," **************\n",sep="")
    
    #active agents (i.e. those in Philadelphia today)
    active = sex_network$ID[sex_network$ENTRY<=day & sex_network$EXIT>=day]
    
    #go through each prevention paradigm scenario
    for (current_data in 1:nparadigms)
    {
      #cat("\n\n************** ","Scenario: ",current_data," **************\n",sep="")
      
      ## SET DAILY VARIABLES ##
      
      #ensure each scenario is identical for a given day
      set.seed(seed_val_day[day])
      
      #set active agents
      abm_sims[[current_data]]$ACTIVE = 0
      abm_sims[[current_data]]$ACTIVE[active] = 1
      
      #check for annual STI testing today
      abm_sims[[current_data]]$HIV_STATUS = ifelse((!is.na(abm_sims[[current_data]]$STI_TEST_DAY) & ((day-abm_sims[[current_data]]$STI_TEST_DAY) %% 365)==0) & (abm_sims[[current_data]]$HIV==1), 3, abm_sims[[current_data]]$HIV_STATUS)
      
      #check if testing identify as positive, used for TAP, after 30 days you may be virally suppressed
      abm_sims[[current_data]]$DAYS_KNOWN_HIV_POSITIVE = ifelse((!is.na(abm_sims[[current_data]]$STI_TEST_DAY)& (day==abm_sims[[current_data]]$STI_TEST_DAY) & (abm_sims[[current_data]]$HIV==1)), 1, abm_sims[[current_data]]$DAYS_KNOWN_HIV_POSITIVE)
      
      #increment days positive, if positive
      abm_sims[[current_data]]$DAYS_KNOWN_HIV_POSITIVE = ifelse(abm_sims[[current_data]]$DAYS_KNOWN_HIV_POSITIVE>=0, abm_sims[[current_data]]$DAYS_KNOWN_HIV_POSITIVE+1, abm_sims[[current_data]]$DAYS_KNOWN_HIV_POSITIVE)
      
      #in a TAP scenario, if put on treatment, after 30 days viral load goes to undetectable
      abm_sims[[current_data]]$TAPVL = ifelse(abm_sims[[current_data]]$TAP==1 & abm_sims[[current_data]]$DAYS_KNOWN_HIV_POSITIVE==30 & abm_sims[[current_data]]$HIV==1 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.42, 1, abm_sims[[current_data]]$TAPVL)
      
      #increment HIV incident days
      abm_sims[[current_data]]$INCIDENCE_DAYS_HIV = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HIV>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HIV+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HIV)
      
      #resolve acute phase of HIV infection
      abm_sims[[current_data]]$HIV_VL = ifelse(abm_sims[[current_data]]$HIV_VL==2 & abm_sims[[current_data]]$INCIDENCE_DAYS_HIV==63, 1, abm_sims[[current_data]]$HIV_VL)
      
      #if HIV positive and aware and set as receptive
      abm_sims[[current_data]]$SEROTYPE = ifelse(abm_sims[[current_data]]$HIV_STATUS==3 & abm_sims[[current_data]]$SEROTYPE==1, 2, abm_sims[[current_data]]$SEROTYPE) 
      
      #increment HPV incident days
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ORAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ORAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ORAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ORAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ORAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ORAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ORAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ORAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ORAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ORAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ORAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ORAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ORAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ORAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ORAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ORAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ORAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ORAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ORAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ORAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ORAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ORAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ORAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ORAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ORAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ORAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ORAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ORAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ORAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ORAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ORAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ORAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL)
      abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ORAL = ifelse(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ORAL>=0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ORAL+1, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ORAL)
      
      #resolve spontaneous HPV clearance based on HIV status
      abm_sims[[current_data]]$HPV6_ANAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV6_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV6_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV6_ANAL))
      abm_sims[[current_data]]$HPV6_ORAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV6_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV6_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV6_ORAL))
      abm_sims[[current_data]]$HPV11_ANAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV11_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV11_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV11_ANAL))
      abm_sims[[current_data]]$HPV11_ORAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV11_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV11_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV11_ORAL))
      abm_sims[[current_data]]$HPV16_ANAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV16_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV16_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV16_ANAL))
      abm_sims[[current_data]]$HPV16_ORAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV16_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV16_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV16_ORAL))
      abm_sims[[current_data]]$HPV18_ANAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV18_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV18_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV18_ANAL))
      abm_sims[[current_data]]$HPV18_ORAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV18_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV18_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV18_ORAL))
      abm_sims[[current_data]]$HPV31_ANAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV31_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV31_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV31_ANAL))
      abm_sims[[current_data]]$HPV31_ORAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV31_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV31_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV31_ORAL))
      abm_sims[[current_data]]$HPV33_ANAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV33_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV33_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV33_ANAL))
      abm_sims[[current_data]]$HPV33_ORAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV33_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV33_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV33_ORAL))
      abm_sims[[current_data]]$HPV45_ANAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV45_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV45_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV45_ANAL))
      abm_sims[[current_data]]$HPV45_ORAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV45_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV45_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV45_ORAL))
      abm_sims[[current_data]]$HPV52_ANAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV52_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV52_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV52_ANAL))
      abm_sims[[current_data]]$HPV52_ORAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV52_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV52_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV52_ORAL))
      abm_sims[[current_data]]$HPV58_ANAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV58_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV58_ANAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV58_ANAL))
      abm_sims[[current_data]]$HPV58_ORAL = ifelse(abm_sims[[current_data]]$HIV==0 & abm_sims[[current_data]]$HPV58_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.66, -1, ifelse(abm_sims[[current_data]]$HIV==1 & abm_sims[[current_data]]$HPV58_ORAL==1 & abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ORAL==228 & runif(nrow(abm_sims[[current_data]]),0,1)<=0.59, -1, abm_sims[[current_data]]$HPV58_ORAL))
      abm_sims[[current_data]]$HPV6_ANAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV6_ANAL==-1, abm_sims[[current_data]]$HPV6_ANAL_CLEARED+1, abm_sims[[current_data]]$HPV6_ANAL_CLEARED)
      abm_sims[[current_data]]$HPV6_ANAL[abm_sims[[current_data]]$HPV6_ANAL==-1] = 0
      abm_sims[[current_data]]$HPV6_ORAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV6_ORAL==-1, abm_sims[[current_data]]$HPV6_ORAL_CLEARED+1, abm_sims[[current_data]]$HPV6_ORAL_CLEARED)
      abm_sims[[current_data]]$HPV6_ORAL[abm_sims[[current_data]]$HPV6_ORAL==-1] = 0
      abm_sims[[current_data]]$HPV11_ANAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV11_ANAL==-1, abm_sims[[current_data]]$HPV11_ANAL_CLEARED+1, abm_sims[[current_data]]$HPV11_ANAL_CLEARED)
      abm_sims[[current_data]]$HPV11_ANAL[abm_sims[[current_data]]$HPV11_ANAL==-1] = 0
      abm_sims[[current_data]]$HPV11_ORAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV11_ORAL==-1, abm_sims[[current_data]]$HPV11_ORAL_CLEARED+1, abm_sims[[current_data]]$HPV11_ORAL_CLEARED)
      abm_sims[[current_data]]$HPV11_ORAL[abm_sims[[current_data]]$HPV11_ORAL==-1] = 0
      abm_sims[[current_data]]$HPV16_ANAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV16_ANAL==-1, abm_sims[[current_data]]$HPV16_ANAL_CLEARED+1, abm_sims[[current_data]]$HPV16_ANAL_CLEARED)
      abm_sims[[current_data]]$HPV16_ANAL[abm_sims[[current_data]]$HPV16_ANAL==-1] = 0
      abm_sims[[current_data]]$HPV16_ORAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV16_ORAL==-1, abm_sims[[current_data]]$HPV16_ORAL_CLEARED+1, abm_sims[[current_data]]$HPV16_ORAL_CLEARED)
      abm_sims[[current_data]]$HPV16_ORAL[abm_sims[[current_data]]$HPV16_ORAL==-1] = 0
      abm_sims[[current_data]]$HPV18_ANAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV18_ANAL==-1, abm_sims[[current_data]]$HPV18_ANAL_CLEARED+1, abm_sims[[current_data]]$HPV18_ANAL_CLEARED)
      abm_sims[[current_data]]$HPV18_ANAL[abm_sims[[current_data]]$HPV18_ANAL==-1] = 0
      abm_sims[[current_data]]$HPV18_ORAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV18_ORAL==-1, abm_sims[[current_data]]$HPV18_ORAL_CLEARED+1, abm_sims[[current_data]]$HPV18_ORAL_CLEARED)
      abm_sims[[current_data]]$HPV18_ORAL[abm_sims[[current_data]]$HPV18_ORAL==-1] = 0
      abm_sims[[current_data]]$HPV31_ANAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV31_ANAL==-1, abm_sims[[current_data]]$HPV31_ANAL_CLEARED+1, abm_sims[[current_data]]$HPV31_ANAL_CLEARED)
      abm_sims[[current_data]]$HPV31_ANAL[abm_sims[[current_data]]$HPV31_ANAL==-1] = 0
      abm_sims[[current_data]]$HPV31_ORAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV31_ORAL==-1, abm_sims[[current_data]]$HPV31_ORAL_CLEARED+1, abm_sims[[current_data]]$HPV31_ORAL_CLEARED)
      abm_sims[[current_data]]$HPV31_ORAL[abm_sims[[current_data]]$HPV31_ORAL==-1] = 0
      abm_sims[[current_data]]$HPV33_ANAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV33_ANAL==-1, abm_sims[[current_data]]$HPV33_ANAL_CLEARED+1, abm_sims[[current_data]]$HPV33_ANAL_CLEARED)
      abm_sims[[current_data]]$HPV33_ANAL[abm_sims[[current_data]]$HPV33_ANAL==-1] = 0
      abm_sims[[current_data]]$HPV33_ORAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV33_ORAL==-1, abm_sims[[current_data]]$HPV33_ORAL_CLEARED+1, abm_sims[[current_data]]$HPV33_ORAL_CLEARED)
      abm_sims[[current_data]]$HPV33_ORAL[abm_sims[[current_data]]$HPV33_ORAL==-1] = 0
      abm_sims[[current_data]]$HPV45_ANAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV45_ANAL==-1, abm_sims[[current_data]]$HPV45_ANAL_CLEARED+1, abm_sims[[current_data]]$HPV45_ANAL_CLEARED)
      abm_sims[[current_data]]$HPV45_ANAL[abm_sims[[current_data]]$HPV45_ANAL==-1] = 0
      abm_sims[[current_data]]$HPV45_ORAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV45_ORAL==-1, abm_sims[[current_data]]$HPV45_ORAL_CLEARED+1, abm_sims[[current_data]]$HPV45_ORAL_CLEARED)
      abm_sims[[current_data]]$HPV45_ORAL[abm_sims[[current_data]]$HPV45_ORAL==-1] = 0
      abm_sims[[current_data]]$HPV52_ANAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV52_ANAL==-1, abm_sims[[current_data]]$HPV52_ANAL_CLEARED+1, abm_sims[[current_data]]$HPV52_ANAL_CLEARED)
      abm_sims[[current_data]]$HPV52_ANAL[abm_sims[[current_data]]$HPV52_ANAL==-1] = 0
      abm_sims[[current_data]]$HPV52_ORAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV52_ORAL==-1, abm_sims[[current_data]]$HPV52_ORAL_CLEARED+1, abm_sims[[current_data]]$HPV52_ORAL_CLEARED)
      abm_sims[[current_data]]$HPV52_ORAL[abm_sims[[current_data]]$HPV52_ORAL==-1] = 0
      abm_sims[[current_data]]$HPV58_ANAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV58_ANAL==-1, abm_sims[[current_data]]$HPV58_ANAL_CLEARED+1, abm_sims[[current_data]]$HPV58_ANAL_CLEARED)
      abm_sims[[current_data]]$HPV58_ANAL[abm_sims[[current_data]]$HPV58_ANAL==-1] = 0
      abm_sims[[current_data]]$HPV58_ORAL_CLEARED = ifelse(abm_sims[[current_data]]$HPV58_ORAL==-1, abm_sims[[current_data]]$HPV58_ORAL_CLEARED+1, abm_sims[[current_data]]$HPV58_ORAL_CLEARED)
      abm_sims[[current_data]]$HPV58_ORAL[abm_sims[[current_data]]$HPV58_ORAL==-1] = 0
      abm_sims[[current_data]]$HPV_ANAL = ifelse(abm_sims[[current_data]]$HPV6_ANAL==1 | abm_sims[[current_data]]$HPV11_ANAL==1 | abm_sims[[current_data]]$HPV16_ANAL==1 | abm_sims[[current_data]]$HPV18_ANAL==1 | abm_sims[[current_data]]$HPV31_ANAL==1 | abm_sims[[current_data]]$HPV33_ANAL==1 | abm_sims[[current_data]]$HPV45_ANAL==1 | abm_sims[[current_data]]$HPV52_ANAL==1 | abm_sims[[current_data]]$HPV58_ANAL==1, 1, 0)
      abm_sims[[current_data]]$HPV_ORAL = ifelse(abm_sims[[current_data]]$HPV6_ORAL==1 | abm_sims[[current_data]]$HPV11_ORAL==1 | abm_sims[[current_data]]$HPV16_ORAL==1 | abm_sims[[current_data]]$HPV18_ORAL==1 | abm_sims[[current_data]]$HPV31_ORAL==1 | abm_sims[[current_data]]$HPV33_ORAL==1 | abm_sims[[current_data]]$HPV45_ORAL==1 | abm_sims[[current_data]]$HPV52_ORAL==1 | abm_sims[[current_data]]$HPV58_ORAL==1, 1, 0)
      
      ## ENGAGE IN SEX ##
      
      #active partnerships today
      sex_selection = data.frame("Ego"=sex_network$ID[which(!is.na(sex_network[,6+day]))], "Partner"=as.numeric(na.omit(sex_network[,6+day])), stringsAsFactors=F)
      
      #sex today, multiplier is to make the algorithm more conservative from calibration (i.e. less sex)
      sex_selection$Ego_sex = (sex_network$N_SEX_ACTS[sex_selection$Ego] / sex_network$RELATIONSHIP_DAYS[sex_selection$Ego]) <= (runif(nrow(sex_selection),0,1)*0.45)
      sex_selection$Partner_sex = (sex_network$N_SEX_ACTS[sex_selection$Partner] / sex_network$RELATIONSHIP_DAYS[sex_selection$Partner]) <= (runif(nrow(sex_selection),0,1)*0.45)
      sex_selection = sex_selection[-1*which(sex_selection$Ego_sex==F | sex_selection$Partner_sex==F), ]
      
      #remove the duplicates (if ego has sex with partner, partner has sex with ego)
      dedup=T; i=1
      while(dedup) {
        dup = which(sex_selection$Ego[i]==sex_selection$Partner)
        if (length(dup)>0) { sex_selection = sex_selection[-dup, ] }
        i=i+1
        if (i>nrow(sex_selection)) { dedup=F }
      }
      rm(dedup,i,dup)
      
      #ensuring at least one partnership
      if (nrow(sex_selection)>0)
      {
        #create the exposure matrix for each pathogen
        sex_exposed = data.frame("Ego"=sex_selection$Ego, "Partner"=sex_selection$Partner, "Exposed"=0, "HIV"=0, "HPV6_ANAL"=0, "HPV6_ORAL"=0, "HPV11_ANAL"=0, "HPV11_ORAL"=0, "HPV16_ANAL"=0, "HPV16_ORAL"=0, "HPV18_ANAL"=0, "HPV18_ORAL"=0, "HPV31_ANAL"=0, "HPV31_ORAL"=0, "HPV33_ANAL"=0, "HPV33_ORAL"=0, "HPV45_ANAL"=0, "HPV45_ORAL"=0, "HPV52_ANAL"=0, "HPV52_ORAL"=0, "HPV58_ANAL"=0, "HPV58_ORAL"=0, stringsAsFactors=F)
        sex_exposed$HIV[which((abm_sims[[current_data]]$HIV[sex_selection$Ego] + abm_sims[[current_data]]$HIV[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV6_ANAL[which((abm_sims[[current_data]]$HPV6_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV6_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV11_ANAL[which((abm_sims[[current_data]]$HPV11_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV11_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV16_ANAL[which((abm_sims[[current_data]]$HPV16_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV16_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV18_ANAL[which((abm_sims[[current_data]]$HPV18_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV18_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV31_ANAL[which((abm_sims[[current_data]]$HPV31_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV31_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV33_ANAL[which((abm_sims[[current_data]]$HPV33_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV33_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV45_ANAL[which((abm_sims[[current_data]]$HPV45_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV45_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV52_ANAL[which((abm_sims[[current_data]]$HPV52_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV52_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV58_ANAL[which((abm_sims[[current_data]]$HPV58_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV58_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV6_ORAL[which((abm_sims[[current_data]]$HPV6_ORAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV6_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV6_ORAL[which((abm_sims[[current_data]]$HPV6_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV6_ORAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV11_ORAL[which((abm_sims[[current_data]]$HPV11_ORAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV11_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV11_ORAL[which((abm_sims[[current_data]]$HPV11_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV11_ORAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV16_ORAL[which((abm_sims[[current_data]]$HPV16_ORAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV16_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV16_ORAL[which((abm_sims[[current_data]]$HPV16_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV16_ORAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV18_ORAL[which((abm_sims[[current_data]]$HPV18_ORAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV18_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV18_ORAL[which((abm_sims[[current_data]]$HPV18_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV18_ORAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV31_ORAL[which((abm_sims[[current_data]]$HPV31_ORAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV31_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV31_ORAL[which((abm_sims[[current_data]]$HPV31_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV31_ORAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV33_ORAL[which((abm_sims[[current_data]]$HPV33_ORAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV33_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV33_ORAL[which((abm_sims[[current_data]]$HPV33_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV33_ORAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV45_ORAL[which((abm_sims[[current_data]]$HPV45_ORAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV45_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV45_ORAL[which((abm_sims[[current_data]]$HPV45_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV45_ORAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV52_ORAL[which((abm_sims[[current_data]]$HPV52_ORAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV52_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV52_ORAL[which((abm_sims[[current_data]]$HPV52_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV52_ORAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV58_ORAL[which((abm_sims[[current_data]]$HPV58_ORAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV58_ANAL[sex_selection$Partner])==1)] = 1
        sex_exposed$HPV58_ORAL[which((abm_sims[[current_data]]$HPV58_ANAL[sex_selection$Ego] + abm_sims[[current_data]]$HPV58_ORAL[sex_selection$Partner])==1)] = 1
        
        #collapse HPV types into single HPV variable
        #sex_exposed$HPV_ANAL = ifelse(sex_exposed$HPV6_ANAL==1 | sex_exposed$HPV11_ANAL==1 | sex_exposed$HPV16_ANAL==1 | sex_exposed$HPV18_ANAL==1 | sex_exposed$HPV31_ANAL==1 | sex_exposed$HPV33_ANAL==1 | sex_exposed$HPV45_ANAL==1 | sex_exposed$HPV52_ANAL==1 | sex_exposed$HPV58_ANAL==1, 1, 0)
        #sex_exposed$HPV_ORAL = ifelse(sex_exposed$HPV6_ORAL==1 | sex_exposed$HPV11_ORAL==1 | sex_exposed$HPV16_ORAL==1 | sex_exposed$HPV18_ORAL==1 | sex_exposed$HPV31_ORAL==1 | sex_exposed$HPV33_ORAL==1 | sex_exposed$HPV45_ORAL==1 | sex_exposed$HPV52_ORAL==1 | sex_exposed$HPV58_ORAL==1, 1, 0)
        
        #HIV/HPV anal coinfections (used for seroadaption scenario)
        #sex_exposed$HIV_HPV_ANAL = sex_exposed$HIV * sex_exposed$HPV_ANAL
        sex_exposed$HIV_HPV6_ANAL = sex_exposed$HIV * sex_exposed$HPV6_ANAL
        sex_exposed$HIV_HPV11_ANAL = sex_exposed$HIV * sex_exposed$HPV11_ANAL
        sex_exposed$HIV_HPV16_ANAL = sex_exposed$HIV * sex_exposed$HPV16_ANAL
        sex_exposed$HIV_HPV18_ANAL = sex_exposed$HIV * sex_exposed$HPV18_ANAL
        sex_exposed$HIV_HPV31_ANAL = sex_exposed$HIV * sex_exposed$HPV31_ANAL
        sex_exposed$HIV_HPV33_ANAL = sex_exposed$HIV * sex_exposed$HPV33_ANAL
        sex_exposed$HIV_HPV45_ANAL = sex_exposed$HIV * sex_exposed$HPV45_ANAL
        sex_exposed$HIV_HPV52_ANAL = sex_exposed$HIV * sex_exposed$HPV52_ANAL
        sex_exposed$HIV_HPV58_ANAL = sex_exposed$HIV * sex_exposed$HPV58_ANAL
        
        #check for unexposed
        sex_exposed$Exposed = rowSums(sex_exposed[,c("HIV","HPV6_ANAL","HPV6_ORAL","HPV11_ANAL","HPV11_ORAL","HPV16_ANAL","HPV16_ORAL","HPV18_ANAL","HPV18_ORAL","HPV31_ANAL","HPV31_ORAL","HPV33_ANAL","HPV33_ORAL","HPV45_ANAL","HPV45_ORAL","HPV52_ANAL","HPV52_ORAL","HPV58_ANAL","HPV58_ORAL")])
        sex_unexposed = c(sex_exposed$Ego[sex_exposed$Exposed==0], sex_exposed$Partner[sex_exposed$Exposed==0])
        sex_exposed = sex_exposed[sex_exposed$Exposed!=0, ]
        rm(sex_selection)
        
        #assign for anal intercourse (AI) sex position, 1=ego receptive/partner insertive, 2=ego insertive/partner receptive, 3=flip fuck 
        #based on preference first, then scale of dominance if preference cannot be honored
        sex_exposed$simAI = ifelse(abm_sims[[current_data]]$SEX_POSITION_ANAL_PREF[sex_exposed$Ego]==1 & abm_sims[[current_data]]$SEX_POSITION_ANAL_PREF[sex_exposed$Partner]==2, 2, ifelse(abm_sims[[current_data]]$SEX_POSITION_ANAL_PREF[sex_exposed$Ego]==2 & abm_sims[[current_data]]$SEX_POSITION_ANAL_PREF[sex_exposed$Partner]==1, 1, ifelse(abm_sims[[current_data]]$SEX_POSITION_ANAL_PREF[sex_exposed$Ego]==3 & abm_sims[[current_data]]$SEX_POSITION_ANAL_PREF[sex_exposed$Partner]==3, 3, NA)))
        sex_exposed$simAI = ifelse(is.na(sex_exposed$simAI), ifelse(abs(abm_sims[[current_data]]$SEX_POSITION_ANAL[sex_exposed$Ego]-abm_sims[[current_data]]$SEX_POSITION_ANAL[sex_exposed$Partner])<=0.15, 3, ifelse(abm_sims[[current_data]]$SEX_POSITION_ANAL[sex_exposed$Ego]<abm_sims[[current_data]]$SEX_POSITION_ANAL[sex_exposed$Partner], 1, 2)), sex_exposed$simAI)
        
        #assign for oral intercourse (OI) sex position, 1=ego receptive/partner insertive, 2=ego insertive/partner receptive, 3=flip fuck 
        #based on preference first, then scale of dominance if preference cannot be honored
        sex_exposed$simOI = ifelse(abm_sims[[current_data]]$SEX_POSITION_ORAL_PREF[sex_exposed$Ego]==1 & abm_sims[[current_data]]$SEX_POSITION_ORAL_PREF[sex_exposed$Partner]==2, 2, ifelse(abm_sims[[current_data]]$SEX_POSITION_ORAL_PREF[sex_exposed$Ego]==2 & abm_sims[[current_data]]$SEX_POSITION_ORAL_PREF[sex_exposed$Partner]==1, 1, ifelse(abm_sims[[current_data]]$SEX_POSITION_ORAL_PREF[sex_exposed$Ego]==3 & abm_sims[[current_data]]$SEX_POSITION_ORAL_PREF[sex_exposed$Partner]==3, 3, NA)))
        sex_exposed$simOI = ifelse(is.na(sex_exposed$simOI), ifelse(abs(abm_sims[[current_data]]$SEX_POSITION_ORAL[sex_exposed$Ego]-abm_sims[[current_data]]$SEX_POSITION_ORAL[sex_exposed$Partner])<=0.15, 3, ifelse(abm_sims[[current_data]]$SEX_POSITION_ORAL[sex_exposed$Ego]<abm_sims[[current_data]]$SEX_POSITION_ORAL[sex_exposed$Partner], 1, 2)), sex_exposed$simOI)
        
        #determine sex type, 1=anal (both consent to anal only), 2=oral (one does not consent to anal), 3=both (both consent to anal but one, or both, would also like oral)
        sex_exposed$sex = ifelse(abm_sims[[current_data]]$SEX_ACT[sex_exposed$Ego]==1 & abm_sims[[current_data]]$SEX_ACT[sex_exposed$Partner]==1, 1, ifelse(abm_sims[[current_data]]$SEX_ACT[sex_exposed$Ego]==2 & abm_sims[[current_data]]$SEX_ACT[sex_exposed$Partner]==2, 2, ifelse(abm_sims[[current_data]]$SEX_ACT[sex_exposed$Ego]==3 & abm_sims[[current_data]]$SEX_ACT[sex_exposed$Partner]==3, 3, ifelse(abm_sims[[current_data]]$SEX_ACT[sex_exposed$Ego]==3 & abm_sims[[current_data]]$SEX_ACT[sex_exposed$Partner]==2, 3, ifelse(abm_sims[[current_data]]$SEX_ACT[sex_exposed$Ego]==2 & abm_sims[[current_data]]$SEX_ACT[sex_exposed$Partner]==3, 3, 2)))))
        
        if (sum(sex_exposed$Exposed)>0)
        {
          #pathogen discordance in at least one partnership
          
          ## DETERMINE INFECTION STATUS ##
          
          #create exposure lists of ego/partner IDs
          exposed_hiv_ego = sex_exposed$Ego[sex_exposed$HIV==1]
          exposed_hiv_partner = sex_exposed$Partner[sex_exposed$HIV==1]
          #exposed_hpv_anal_ego = sex_exposed$Ego[sex_exposed$HPV_ANAL==1]
          #exposed_hpv_anal_partner = sex_exposed$Partner[sex_exposed$HPV_ANAL==1]
          exposed_hpv6_anal_ego = sex_exposed$Ego[sex_exposed$HPV6_ANAL==1]
          exposed_hpv6_anal_partner = sex_exposed$Partner[sex_exposed$HPV6_ANAL==1]
          exposed_hpv11_anal_ego = sex_exposed$Ego[sex_exposed$HPV11_ANAL==1]
          exposed_hpv11_anal_partner = sex_exposed$Partner[sex_exposed$HPV11_ANAL==1]
          exposed_hpv16_anal_ego = sex_exposed$Ego[sex_exposed$HPV16_ANAL==1]
          exposed_hpv16_anal_partner = sex_exposed$Partner[sex_exposed$HPV16_ANAL==1]
          exposed_hpv18_anal_ego = sex_exposed$Ego[sex_exposed$HPV18_ANAL==1]
          exposed_hpv18_anal_partner = sex_exposed$Partner[sex_exposed$HPV18_ANAL==1]
          exposed_hpv31_anal_ego = sex_exposed$Ego[sex_exposed$HPV31_ANAL==1]
          exposed_hpv31_anal_partner = sex_exposed$Partner[sex_exposed$HPV31_ANAL==1]
          exposed_hpv33_anal_ego = sex_exposed$Ego[sex_exposed$HPV33_ANAL==1]
          exposed_hpv33_anal_partner = sex_exposed$Partner[sex_exposed$HPV33_ANAL==1]
          exposed_hpv45_anal_ego = sex_exposed$Ego[sex_exposed$HPV45_ANAL==1]
          exposed_hpv45_anal_partner = sex_exposed$Partner[sex_exposed$HPV45_ANAL==1]
          exposed_hpv52_anal_ego = sex_exposed$Ego[sex_exposed$HPV52_ANAL==1]
          exposed_hpv52_anal_partner = sex_exposed$Partner[sex_exposed$HPV52_ANAL==1]
          exposed_hpv58_anal_ego = sex_exposed$Ego[sex_exposed$HPV58_ANAL==1]
          exposed_hpv58_anal_partner = sex_exposed$Partner[sex_exposed$HPV58_ANAL==1]
          exposed_hpv6_oral_ego = sex_exposed$Ego[sex_exposed$HPV6_ORAL==1]
          exposed_hpv6_oral_partner = sex_exposed$Partner[sex_exposed$HPV6_ORAL==1]
          exposed_hpv11_oral_ego = sex_exposed$Ego[sex_exposed$HPV11_ORAL==1]
          exposed_hpv11_oral_partner = sex_exposed$Partner[sex_exposed$HPV11_ORAL==1]
          exposed_hpv16_oral_ego = sex_exposed$Ego[sex_exposed$HPV16_ORAL==1]
          exposed_hpv16_oral_partner = sex_exposed$Partner[sex_exposed$HPV16_ORAL==1]
          exposed_hpv18_oral_ego = sex_exposed$Ego[sex_exposed$HPV18_ORAL==1]
          exposed_hpv18_oral_partner = sex_exposed$Partner[sex_exposed$HPV18_ORAL==1]
          exposed_hpv31_oral_ego = sex_exposed$Ego[sex_exposed$HPV31_ORAL==1]
          exposed_hpv31_oral_partner = sex_exposed$Partner[sex_exposed$HPV31_ORAL==1]
          exposed_hpv33_oral_ego = sex_exposed$Ego[sex_exposed$HPV33_ORAL==1]
          exposed_hpv33_oral_partner = sex_exposed$Partner[sex_exposed$HPV33_ORAL==1]
          exposed_hpv45_oral_ego = sex_exposed$Ego[sex_exposed$HPV45_ORAL==1]
          exposed_hpv45_oral_partner = sex_exposed$Partner[sex_exposed$HPV45_ORAL==1]
          exposed_hpv52_oral_ego = sex_exposed$Ego[sex_exposed$HPV52_ORAL==1]
          exposed_hpv52_oral_partner = sex_exposed$Partner[sex_exposed$HPV52_ORAL==1]
          exposed_hpv58_oral_ego = sex_exposed$Ego[sex_exposed$HPV58_ORAL==1]
          exposed_hpv58_oral_partner = sex_exposed$Partner[sex_exposed$HPV58_ORAL==1]
          
          #check for who is exposed, 1=ego, 2=partner, 3=both (oral)
          sex_exposed$exposure_hiv[sex_exposed$HIV==1] = ifelse(abm_sims[[current_data]]$HIV[exposed_hiv_ego]<abm_sims[[current_data]]$HIV[exposed_hiv_partner], 1, 2)
          #sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1] = ifelse(abm_sims[[current_data]]$HPV_ANAL[exposed_hpv_anal_ego]<abm_sims[[current_data]]$HPV_ANAL[exposed_hpv_anal_partner], 1, ifelse(abm_sims[[current_data]]$HPV_ANAL[exposed_hpv_anal_ego]>abm_sims[[current_data]]$HPV_ANAL[exposed_hpv_anal_partner], 2, 3))
          sex_exposed$exposure_hpv6_anal[sex_exposed$HPV6_ANAL==1] = ifelse(abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_anal_ego]<abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_anal_partner], 1, 2)
          sex_exposed$exposure_hpv11_anal[sex_exposed$HPV11_ANAL==1] = ifelse(abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_anal_ego]<abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_anal_partner], 1, 2)
          sex_exposed$exposure_hpv16_anal[sex_exposed$HPV16_ANAL==1] = ifelse(abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_anal_ego]<abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_anal_partner], 1, 2)
          sex_exposed$exposure_hpv18_anal[sex_exposed$HPV18_ANAL==1] = ifelse(abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_anal_ego]<abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_anal_partner], 1, 2)
          sex_exposed$exposure_hpv31_anal[sex_exposed$HPV31_ANAL==1] = ifelse(abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_anal_ego]<abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_anal_partner], 1, 2)
          sex_exposed$exposure_hpv33_anal[sex_exposed$HPV33_ANAL==1] = ifelse(abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_anal_ego]<abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_anal_partner], 1, 2)
          sex_exposed$exposure_hpv45_anal[sex_exposed$HPV45_ANAL==1] = ifelse(abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_anal_ego]<abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_anal_partner], 1, 2)
          sex_exposed$exposure_hpv52_anal[sex_exposed$HPV52_ANAL==1] = ifelse(abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_anal_ego]<abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_anal_partner], 1, 2)
          sex_exposed$exposure_hpv58_anal[sex_exposed$HPV58_ANAL==1] = ifelse(abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_anal_ego]<abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_anal_partner], 1, 2)
          sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1] = ifelse((abm_sims[[current_data]]$HPV6_ORAL[exposed_hpv6_oral_ego]<abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_oral_partner]) & (abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_oral_ego]>abm_sims[[current_data]]$HPV6_ORAL[exposed_hpv6_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV6_ORAL[exposed_hpv6_oral_ego]>abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_oral_partner]) & (abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_oral_ego]<abm_sims[[current_data]]$HPV6_ORAL[exposed_hpv6_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV6_ORAL[exposed_hpv6_oral_ego]<abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_oral_partner]) | (abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_oral_ego]<abm_sims[[current_data]]$HPV6_ORAL[exposed_hpv6_oral_partner]), 1, ifelse((abm_sims[[current_data]]$HPV6_ORAL[exposed_hpv6_oral_ego]>abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_oral_partner]) | (abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_oral_ego]>abm_sims[[current_data]]$HPV6_ORAL[exposed_hpv6_oral_partner]), 2, NA))))
          sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1] = ifelse((abm_sims[[current_data]]$HPV11_ORAL[exposed_hpv11_oral_ego]<abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_oral_partner]) & (abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_oral_ego]>abm_sims[[current_data]]$HPV11_ORAL[exposed_hpv11_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV11_ORAL[exposed_hpv11_oral_ego]>abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_oral_partner]) & (abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_oral_ego]<abm_sims[[current_data]]$HPV11_ORAL[exposed_hpv11_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV11_ORAL[exposed_hpv11_oral_ego]<abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_oral_partner]) | (abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_oral_ego]<abm_sims[[current_data]]$HPV11_ORAL[exposed_hpv11_oral_partner]), 1, ifelse((abm_sims[[current_data]]$HPV11_ORAL[exposed_hpv11_oral_ego]>abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_oral_partner]) | (abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_oral_ego]>abm_sims[[current_data]]$HPV11_ORAL[exposed_hpv11_oral_partner]), 2, NA))))
          sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1] = ifelse((abm_sims[[current_data]]$HPV16_ORAL[exposed_hpv16_oral_ego]<abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_oral_partner]) & (abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_oral_ego]>abm_sims[[current_data]]$HPV16_ORAL[exposed_hpv16_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV16_ORAL[exposed_hpv16_oral_ego]>abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_oral_partner]) & (abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_oral_ego]<abm_sims[[current_data]]$HPV16_ORAL[exposed_hpv16_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV16_ORAL[exposed_hpv16_oral_ego]<abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_oral_partner]) | (abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_oral_ego]<abm_sims[[current_data]]$HPV16_ORAL[exposed_hpv16_oral_partner]), 1, ifelse((abm_sims[[current_data]]$HPV16_ORAL[exposed_hpv16_oral_ego]>abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_oral_partner]) | (abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_oral_ego]>abm_sims[[current_data]]$HPV16_ORAL[exposed_hpv16_oral_partner]), 2, NA))))
          sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1] = ifelse((abm_sims[[current_data]]$HPV18_ORAL[exposed_hpv18_oral_ego]<abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_oral_partner]) & (abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_oral_ego]>abm_sims[[current_data]]$HPV18_ORAL[exposed_hpv18_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV18_ORAL[exposed_hpv18_oral_ego]>abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_oral_partner]) & (abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_oral_ego]<abm_sims[[current_data]]$HPV18_ORAL[exposed_hpv18_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV18_ORAL[exposed_hpv18_oral_ego]<abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_oral_partner]) | (abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_oral_ego]<abm_sims[[current_data]]$HPV18_ORAL[exposed_hpv18_oral_partner]), 1, ifelse((abm_sims[[current_data]]$HPV18_ORAL[exposed_hpv18_oral_ego]>abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_oral_partner]) | (abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_oral_ego]>abm_sims[[current_data]]$HPV18_ORAL[exposed_hpv18_oral_partner]), 2, NA))))
          sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1] = ifelse((abm_sims[[current_data]]$HPV31_ORAL[exposed_hpv31_oral_ego]<abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_oral_partner]) & (abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_oral_ego]>abm_sims[[current_data]]$HPV31_ORAL[exposed_hpv31_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV31_ORAL[exposed_hpv31_oral_ego]>abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_oral_partner]) & (abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_oral_ego]<abm_sims[[current_data]]$HPV31_ORAL[exposed_hpv31_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV31_ORAL[exposed_hpv31_oral_ego]<abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_oral_partner]) | (abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_oral_ego]<abm_sims[[current_data]]$HPV31_ORAL[exposed_hpv31_oral_partner]), 1, ifelse((abm_sims[[current_data]]$HPV31_ORAL[exposed_hpv31_oral_ego]>abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_oral_partner]) | (abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_oral_ego]>abm_sims[[current_data]]$HPV31_ORAL[exposed_hpv31_oral_partner]), 2, NA))))
          sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1] = ifelse((abm_sims[[current_data]]$HPV33_ORAL[exposed_hpv33_oral_ego]<abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_oral_partner]) & (abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_oral_ego]>abm_sims[[current_data]]$HPV33_ORAL[exposed_hpv33_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV33_ORAL[exposed_hpv33_oral_ego]>abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_oral_partner]) & (abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_oral_ego]<abm_sims[[current_data]]$HPV33_ORAL[exposed_hpv33_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV33_ORAL[exposed_hpv33_oral_ego]<abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_oral_partner]) | (abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_oral_ego]<abm_sims[[current_data]]$HPV33_ORAL[exposed_hpv33_oral_partner]), 1, ifelse((abm_sims[[current_data]]$HPV33_ORAL[exposed_hpv33_oral_ego]>abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_oral_partner]) | (abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_oral_ego]>abm_sims[[current_data]]$HPV33_ORAL[exposed_hpv33_oral_partner]), 2, NA))))
          sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1] = ifelse((abm_sims[[current_data]]$HPV45_ORAL[exposed_hpv45_oral_ego]<abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_oral_partner]) & (abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_oral_ego]>abm_sims[[current_data]]$HPV45_ORAL[exposed_hpv45_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV45_ORAL[exposed_hpv45_oral_ego]>abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_oral_partner]) & (abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_oral_ego]<abm_sims[[current_data]]$HPV45_ORAL[exposed_hpv45_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV45_ORAL[exposed_hpv45_oral_ego]<abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_oral_partner]) | (abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_oral_ego]<abm_sims[[current_data]]$HPV45_ORAL[exposed_hpv45_oral_partner]), 1, ifelse((abm_sims[[current_data]]$HPV45_ORAL[exposed_hpv45_oral_ego]>abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_oral_partner]) | (abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_oral_ego]>abm_sims[[current_data]]$HPV45_ORAL[exposed_hpv45_oral_partner]), 2, NA))))
          sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1] = ifelse((abm_sims[[current_data]]$HPV52_ORAL[exposed_hpv52_oral_ego]<abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_oral_partner]) & (abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_oral_ego]>abm_sims[[current_data]]$HPV52_ORAL[exposed_hpv52_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV52_ORAL[exposed_hpv52_oral_ego]>abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_oral_partner]) & (abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_oral_ego]<abm_sims[[current_data]]$HPV52_ORAL[exposed_hpv52_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV52_ORAL[exposed_hpv52_oral_ego]<abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_oral_partner]) | (abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_oral_ego]<abm_sims[[current_data]]$HPV52_ORAL[exposed_hpv52_oral_partner]), 1, ifelse((abm_sims[[current_data]]$HPV52_ORAL[exposed_hpv52_oral_ego]>abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_oral_partner]) | (abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_oral_ego]>abm_sims[[current_data]]$HPV52_ORAL[exposed_hpv52_oral_partner]), 2, NA))))
          sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1] = ifelse((abm_sims[[current_data]]$HPV58_ORAL[exposed_hpv58_oral_ego]<abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_oral_partner]) & (abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_oral_ego]>abm_sims[[current_data]]$HPV58_ORAL[exposed_hpv58_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV58_ORAL[exposed_hpv58_oral_ego]>abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_oral_partner]) & (abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_oral_ego]<abm_sims[[current_data]]$HPV58_ORAL[exposed_hpv58_oral_partner]), 3, ifelse((abm_sims[[current_data]]$HPV58_ORAL[exposed_hpv58_oral_ego]<abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_oral_partner]) | (abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_oral_ego]<abm_sims[[current_data]]$HPV58_ORAL[exposed_hpv58_oral_partner]), 1, ifelse((abm_sims[[current_data]]$HPV58_ORAL[exposed_hpv58_oral_ego]>abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_oral_partner]) | (abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_oral_ego]>abm_sims[[current_data]]$HPV58_ORAL[exposed_hpv58_oral_partner]), 2, NA))))
          
          #indicators for various metrics
          sex_exposed$condct_hiv = 0            #HIV infection blocked from condom use, ego
          #sex_exposed$condct_hpv_anal = 0      #HPV anal infection blocked from condom use, ego
          sex_exposed$condct_hpv6_anal = 0      #HPV6 anal infection blocked from condom use, ego
          sex_exposed$condct_hpv11_anal = 0     #HPV11 anal infection blocked from condom use, ego
          sex_exposed$condct_hpv16_anal = 0     #HPV16 anal infection blocked from condom use, ego
          sex_exposed$condct_hpv18_anal = 0     #HPV18 anal infection blocked from condom use, ego
          sex_exposed$condct_hpv31_anal = 0     #HPV31 anal infection blocked from condom use, ego
          sex_exposed$condct_hpv33_anal = 0     #HPV33 anal infection blocked from condom use, ego
          sex_exposed$condct_hpv45_anal = 0     #HPV45 anal infection blocked from condom use, ego
          sex_exposed$condct_hpv52_anal = 0     #HPV52 anal infection blocked from condom use, ego
          sex_exposed$condct_hpv58_anal = 0     #HPV58 anal infection blocked from condom use, ego
          sex_exposed$condct_hpv6_oral = 0      #HPV6 oral infection blocked from condom use, ego
          sex_exposed$condct_hpv11_oral = 0     #HPV11 oral infection blocked from condom use, ego
          sex_exposed$condct_hpv16_oral = 0     #HPV16 oral infection blocked from condom use, ego
          sex_exposed$condct_hpv18_oral = 0     #HPV18 oral infection blocked from condom use, ego
          sex_exposed$condct_hpv31_oral = 0     #HPV31 oral infection blocked from condom use, ego
          sex_exposed$condct_hpv33_oral = 0     #HPV33 oral infection blocked from condom use, ego
          sex_exposed$condct_hpv45_oral = 0     #HPV45 oral infection blocked from condom use, ego
          sex_exposed$condct_hpv52_oral = 0     #HPV52 oral infection blocked from condom use, ego
          sex_exposed$condct_hpv58_oral = 0     #HPV58 oral infection blocked from condom use, ego
          sex_exposed$prepct = 0                #infection blocked from PrEP, ego
          sex_exposed$tapct = 0                 #infection blocked from TAP, ego
          sex_exposed$seroct_hiv = 0            #HIV infection blocked from seroadapation, ego
          sex_exposed$seroct_hpv = 0            #HPV anal infection blocked from seroadapation, ego
          sex_exposed$seroct_hpv6 = 0           #HPV6 anal infection blocked from seroadapation, ego
          sex_exposed$seroct_hpv11 = 0          #HPV11 anal infection blocked from seroadapation, ego
          sex_exposed$seroct_hpv16 = 0          #HPV16 anal infection blocked from seroadapation, ego
          sex_exposed$seroct_hpv18 = 0          #HPV18 anal infection blocked from seroadapation, ego
          sex_exposed$seroct_hpv31 = 0          #HPV31 anal infection blocked from seroadapation, ego
          sex_exposed$seroct_hpv33 = 0          #HPV33 anal infection blocked from seroadapation, ego
          sex_exposed$seroct_hpv45 = 0          #HPV45 anal infection blocked from seroadapation, ego
          sex_exposed$seroct_hpv52 = 0          #HPV52 anal infection blocked from seroadapation, ego
          sex_exposed$seroct_hpv58 = 0          #HPV58 anal infection blocked from seroadapation, ego
          #sex_exposed$vaxct_hpv_anal = 0       #HPV anal infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv6_anal = 0       #HPV6 anal infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv11_anal = 0      #HPV11 anal infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv16_anal = 0      #HPV16 anal infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv18_anal = 0      #HPV18 anal infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv31_anal = 0      #HPV31 anal infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv33_anal = 0      #HPV33 anal infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv45_anal = 0      #HPV45 anal infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv52_anal = 0      #HPV52 anal infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv58_anal = 0      #HPV58 anal infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv6_oral = 0       #HPV6 oral infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv11_oral = 0      #HPV11 oral infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv16_oral = 0      #HPV16 oral infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv18_oral = 0      #HPV18 oral infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv31_oral = 0      #HPV31 oral infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv33_oral = 0      #HPV33 oral infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv45_oral = 0      #HPV45 oral infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv52_oral = 0      #HPV52 oral infection blocked from vaccination, ego
          sex_exposed$vaxct_hpv58_oral = 0      #HPV58 oral infection blocked from vaccination, ego
          #sex_exposed$prevct_hpv_anal = 0      #HPV anal infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv6_anal = 0      #HPV6 anal infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv11_anal = 0     #HPV11 anal infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv16_anal = 0     #HPV16 anal infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv18_anal = 0     #HPV18 anal infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv31_anal = 0     #HPV31 anal infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv33_anal = 0     #HPV33 anal infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv45_anal = 0     #HPV45 anal infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv52_anal = 0     #HPV52 anal infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv58_anal = 0     #HPV58 anal infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv6_oral = 0      #HPV6 oral infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv11_oral = 0     #HPV11 oral infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv16_oral = 0     #HPV16 oral infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv18_oral = 0     #HPV18 oral infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv31_oral = 0     #HPV31 oral infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv33_oral = 0     #HPV33 oral infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv45_oral = 0     #HPV45 oral infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv52_oral = 0     #HPV52 oral infection blocked from any prevention strategy, ego
          sex_exposed$prevct_hpv58_oral = 0     #HPV58 oral infection blocked from any prevention strategy, ego
          sex_exposed$infectct = 0              #ego caused a new HPV anal infection to partner
          sex_exposed$newHIV = 0                #new HIV infection, ego
          sex_exposed$newVL = -1                #set HIV_VL to acute if new infection, ego
          #sex_exposed$newHPV_ANAL = 0          #new HPV anal infection, ego
          sex_exposed$newHPV6_ANAL = 0          #new HPV6 anal infection, ego
          sex_exposed$newHPV11_ANAL = 0         #new HPV11 anal infection, ego
          sex_exposed$newHPV16_ANAL = 0         #new HPV16 anal infection, ego
          sex_exposed$newHPV18_ANAL = 0         #new HPV18 anal infection, ego
          sex_exposed$newHPV31_ANAL = 0         #new HPV31 anal infection, ego
          sex_exposed$newHPV33_ANAL = 0         #new HPV33 anal infection, ego
          sex_exposed$newHPV45_ANAL = 0         #new HPV45 anal infection, ego
          sex_exposed$newHPV52_ANAL = 0         #new HPV52 anal infection, ego
          sex_exposed$newHPV58_ANAL = 0         #new HPV58 anal infection, ego
          sex_exposed$newHPV6_ORAL = 0          #new HPV6 oral infection, ego
          sex_exposed$newHPV11_ORAL = 0         #new HPV11 oral infection, ego
          sex_exposed$newHPV16_ORAL = 0         #new HPV16 oral infection, ego
          sex_exposed$newHPV18_ORAL = 0         #new HPV18 oral infection, ego
          sex_exposed$newHPV31_ORAL = 0         #new HPV31 oral infection, ego
          sex_exposed$newHPV33_ORAL = 0         #new HPV33 oral infection, ego
          sex_exposed$newHPV45_ORAL = 0         #new HPV45 oral infection, ego
          sex_exposed$newHPV52_ORAL = 0         #new HPV52 oral infection, ego
          sex_exposed$newHPV58_ORAL = 0         #new HPV58 oral infection, ego
          sex_exposed$partcondct_hiv = 0        #HIV infection blocked from condom use, partner
          #sex_exposed$partcondct_hpv_anal = 0  #HPV anal infection blocked from condom use, partner
          sex_exposed$partcondct_hpv6_anal = 0  #HPV6 anal infection blocked from condom use, partner
          sex_exposed$partcondct_hpv11_anal = 0 #HPV11 anal infection blocked from condom use, partner
          sex_exposed$partcondct_hpv16_anal = 0 #HPV16 anal infection blocked from condom use, partner
          sex_exposed$partcondct_hpv18_anal = 0 #HPV18 anal infection blocked from condom use, partner
          sex_exposed$partcondct_hpv31_anal = 0 #HPV31 anal infection blocked from condom use, partner
          sex_exposed$partcondct_hpv33_anal = 0 #HPV33 anal infection blocked from condom use, partner
          sex_exposed$partcondct_hpv45_anal = 0 #HPV45 anal infection blocked from condom use, partner
          sex_exposed$partcondct_hpv52_anal = 0 #HPV52 anal infection blocked from condom use, partner
          sex_exposed$partcondct_hpv58_anal = 0 #HPV58 anal infection blocked from condom use, partner
          sex_exposed$partcondct_hpv6_oral = 0  #HPV6 oral infection blocked from condom use, partner
          sex_exposed$partcondct_hpv11_oral = 0 #HPV11 oral infection blocked from condom use, partner
          sex_exposed$partcondct_hpv16_oral = 0 #HPV16 oral infection blocked from condom use, partner
          sex_exposed$partcondct_hpv18_oral = 0 #HPV18 oral infection blocked from condom use, partner
          sex_exposed$partcondct_hpv31_oral = 0 #HPV31 oral infection blocked from condom use, partner
          sex_exposed$partcondct_hpv33_oral = 0 #HPV33 oral infection blocked from condom use, partner
          sex_exposed$partcondct_hpv45_oral = 0 #HPV45 oral infection blocked from condom use, partner
          sex_exposed$partcondct_hpv52_oral = 0 #HPV52 oral infection blocked from condom use, partner
          sex_exposed$partcondct_hpv58_oral = 0 #HPV58 oral infection blocked from condom use, partner
          sex_exposed$partprepct = 0            #infection blocked from PrEP, partner
          sex_exposed$parttapct = 0             #infection blocked from TAP, partner
          sex_exposed$partseroct_hiv = 0        #HIV infection blocked from seroadapation, partner
          sex_exposed$partseroct_hpv = 0        #HPV anal infection blocked from seroadapation, partner
          sex_exposed$partseroct_hpv6 = 0       #HPV6 anal infection blocked from seroadapation, partner
          sex_exposed$partseroct_hpv11 = 0      #HPV11 anal infection blocked from seroadapation, partner
          sex_exposed$partseroct_hpv16 = 0      #HPV16 anal infection blocked from seroadapation, partner
          sex_exposed$partseroct_hpv18 = 0      #HPV18 anal infection blocked from seroadapation, partner
          sex_exposed$partseroct_hpv31 = 0      #HPV31 anal infection blocked from seroadapation, partner
          sex_exposed$partseroct_hpv33 = 0      #HPV33 anal infection blocked from seroadapation, partner
          sex_exposed$partseroct_hpv45 = 0      #HPV45 anal infection blocked from seroadapation, partner
          sex_exposed$partseroct_hpv52 = 0      #HPV52 anal infection blocked from seroadapation, partner
          sex_exposed$partseroct_hpv58 = 0      #HPV58 anal infection blocked from seroadapation, partner
          #sex_exposed$partvaxct_hpv_anal = 0   #HPV anal infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv6_anal = 0   #HPV6 anal infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv11_anal = 0  #HPV11 anal infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv16_anal = 0  #HPV16 anal infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv18_anal = 0  #HPV18 anal infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv31_anal = 0  #HPV31 anal infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv33_anal = 0  #HPV33 anal infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv45_anal = 0  #HPV45 anal infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv52_anal = 0  #HPV52 anal infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv58_anal = 0  #HPV58 anal infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv6_oral = 0   #HPV6 oral infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv11_oral = 0  #HPV11 oral infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv16_oral = 0  #HPV16 oral infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv18_oral = 0  #HPV18 oral infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv31_oral = 0  #HPV31 oral infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv33_oral = 0  #HPV33 oral infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv45_oral = 0  #HPV45 oral infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv52_oral = 0  #HPV52 oral infection blocked from vaccination, partner
          sex_exposed$partvaxct_hpv58_oral = 0  #HPV58 oral infection blocked from vaccination, partner
          #sex_exposed$partprevct_hpv_anal = 0  #HPV anal infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv6_anal = 0  #HPV6 anal infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv11_anal = 0 #HPV11 anal infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv16_anal = 0 #HPV16 anal infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv18_anal = 0 #HPV18 anal infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv31_anal = 0 #HPV31 anal infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv33_anal = 0 #HPV33 anal infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv45_anal = 0 #HPV45 anal infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv52_anal = 0 #HPV52 anal infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv58_anal = 0 #HPV58 anal infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv6_oral = 0  #HPV6 oral infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv11_oral = 0 #HPV11 oral infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv16_oral = 0 #HPV16 oral infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv18_oral = 0 #HPV18 oral infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv31_oral = 0 #HPV31 oral infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv33_oral = 0 #HPV33 oral infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv45_oral = 0 #HPV45 oral infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv52_oral = 0 #HPV52 oral infection blocked from any prevention strategy, partner
          sex_exposed$partprevct_hpv58_oral = 0 #HPV58 oral infection blocked from any prevention strategy, partner
          sex_exposed$partinfectct = 0          #partner caused a new HPV anal infection to ego
          sex_exposed$partnewHIV = 0            #new HIV infection, partner
          sex_exposed$partnewVL = -1            #set HIV_VL to acute if new infection, partner
          #sex_exposed$partnewHPV_ANAL = 0      #new HPV anal infection, partner
          sex_exposed$partnewHPV6_ANAL = 0      #new HPV6 anal infection, partner
          sex_exposed$partnewHPV11_ANAL = 0     #new HPV11 anal infection, partner
          sex_exposed$partnewHPV16_ANAL = 0     #new HPV16 anal infection, partner
          sex_exposed$partnewHPV18_ANAL = 0     #new HPV18 anal infection, partner
          sex_exposed$partnewHPV31_ANAL = 0     #new HPV31 anal infection, partner
          sex_exposed$partnewHPV33_ANAL = 0     #new HPV33 anal infection, partner
          sex_exposed$partnewHPV45_ANAL = 0     #new HPV45 anal infection, partner
          sex_exposed$partnewHPV52_ANAL = 0     #new HPV52 anal infection, partner
          sex_exposed$partnewHPV58_ANAL = 0     #new HPV58 anal infection, partner
          sex_exposed$partnewHPV6_ORAL = 0      #new HPV6 oral infection, partner
          sex_exposed$partnewHPV11_ORAL = 0     #new HPV11 oral infection, partner
          sex_exposed$partnewHPV16_ORAL = 0     #new HPV16 oral infection, partner
          sex_exposed$partnewHPV18_ORAL = 0     #new HPV18 oral infection, partner
          sex_exposed$partnewHPV31_ORAL = 0     #new HPV31 oral infection, partner
          sex_exposed$partnewHPV33_ORAL = 0     #new HPV33 oral infection, partner
          sex_exposed$partnewHPV45_ORAL = 0     #new HPV45 oral infection, partner
          sex_exposed$partnewHPV52_ORAL = 0     #new HPV52 oral infection, partner
          sex_exposed$partnewHPV58_ORAL = 0     #new HPV58 oral infection, partner
          
          #calculate probabilities of infection
          probinfect1_hiv = runif(nrow(sex_exposed),0,1)
          probinfect2_hiv = runif(nrow(sex_exposed),0,1)
          #probinfect1_hpv_anal = runif(nrow(sex_exposed),0,1)
          #probinfect2_hpv_anal = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv6_anal = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv6_anal = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv11_anal = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv11_anal = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv16_anal = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv16_anal = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv18_anal = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv18_anal = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv31_anal = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv31_anal = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv33_anal = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv33_anal = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv45_anal = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv45_anal = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv52_anal = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv52_anal = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv58_anal = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv58_anal = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv6_oral = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv6_oral = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv11_oral = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv11_oral = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv16_oral = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv16_oral = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv18_oral = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv18_oral = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv31_oral = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv31_oral = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv33_oral = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv33_oral = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv45_oral = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv45_oral = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv52_oral = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv52_oral = runif(nrow(sex_exposed),0,1)
          probinfect1_hpv58_oral = runif(nrow(sex_exposed),0,1)
          probinfect2_hpv58_oral = runif(nrow(sex_exposed),0,1)
          
          #determine infection under no prevention for each pathogen
          
          #hiv anal, ego then partner for chronic then acute
          sex_exposed$infection_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0134, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0010, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==0 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0059, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0134, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0010, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==0 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0059, 1, 0))))))
          sex_exposed$infection_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.1284, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0101, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==0 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0569, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.1284, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0101, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==1 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==0 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0569, 1, sex_exposed$infection_hiv[sex_exposed$HIV==1]))))))
          sex_exposed$partinfection_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0134, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0010, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==0 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0059, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0134, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0010, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==0 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0059, 1, 0))))))
          sex_exposed$partinfection_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.1284, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0101, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==0 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.0569, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] <= 0.1284, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0101, 1, ifelse(sex_exposed$sex[sex_exposed$HIV==1]!=2 & sex_exposed$exposure_hiv[sex_exposed$HIV==1]==2 & sex_exposed$simAI[sex_exposed$HIV==1]==3 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==0 & probinfect2_hiv[sex_exposed$HIV==1] <= 0.0569, 1, sex_exposed$partinfection_hiv[sex_exposed$HIV==1]))))))
          
          #hpv anal, ego then partner then both
          #sex_exposed$infection_hpv_anal[sex_exposed$HPV_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==1 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==1 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==0 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==3 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==1 & probinfect2_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==0 & probinfect2_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.071, 1, 0))))))
          #sex_exposed$partinfection_hpv_anal[sex_exposed$HPV_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==1 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==1 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==0 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==3 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==1 & probinfect2_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==0 & probinfect2_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.071, 1, 0))))))
          #sex_exposed$infection_hpv_anal[sex_exposed$HPV_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==3 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==1 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==3 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==1 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==3 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==0 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==3 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==3 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==3 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==1 & probinfect2_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==3 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==0 & probinfect2_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.071, 1, sex_exposed$infection_hpv_anal[sex_exposed$HPV_ANAL==1]))))))
          #sex_exposed$partinfection_hpv_anal[sex_exposed$HPV_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==3 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==1 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==3 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==1 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==3 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==0 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==3 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==3 & probinfect1_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==3 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==1 & probinfect2_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV_ANAL==1]!=2 & sex_exposed$exposure_hpv_anal[sex_exposed$HPV_ANAL==1]==3 & sex_exposed$simAI[sex_exposed$HPV_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv_anal_ego]==0 & probinfect2_hpv_anal[sex_exposed$HPV_ANAL==1] <= 0.071, 1, sex_exposed$partinfection_hpv_anal[sex_exposed$HPV_ANAL==1]))))))
          
          #hpv type specific anal, ego then partner
          sex_exposed$infection_hpv6_anal[sex_exposed$HPV6_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV6_ANAL==1]!=2 & sex_exposed$exposure_hpv6_anal[sex_exposed$HPV6_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV6_ANAL==1]==1 & probinfect1_hpv6_anal[sex_exposed$HPV6_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV6_ANAL==1]!=2 & sex_exposed$exposure_hpv6_anal[sex_exposed$HPV6_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV6_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv6_anal_ego]==1 & probinfect1_hpv6_anal[sex_exposed$HPV6_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV6_ANAL==1]!=2 & sex_exposed$exposure_hpv6_anal[sex_exposed$HPV6_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV6_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv6_anal_ego]==0 & probinfect1_hpv6_anal[sex_exposed$HPV6_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV6_ANAL==1]!=2 & sex_exposed$exposure_hpv6_anal[sex_exposed$HPV6_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV6_ANAL==1]==3 & probinfect1_hpv6_anal[sex_exposed$HPV6_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV6_ANAL==1]!=2 & sex_exposed$exposure_hpv6_anal[sex_exposed$HPV6_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV6_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv6_anal_ego]==1 & probinfect2_hpv6_anal[sex_exposed$HPV6_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV6_ANAL==1]!=2 & sex_exposed$exposure_hpv6_anal[sex_exposed$HPV6_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV6_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv6_anal_ego]==0 & probinfect2_hpv6_anal[sex_exposed$HPV6_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$partinfection_hpv6_anal[sex_exposed$HPV6_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV6_ANAL==1]!=2 & sex_exposed$exposure_hpv6_anal[sex_exposed$HPV6_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV6_ANAL==1]==1 & probinfect1_hpv6_anal[sex_exposed$HPV6_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV6_ANAL==1]!=2 & sex_exposed$exposure_hpv6_anal[sex_exposed$HPV6_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV6_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv6_anal_partner]==1 & probinfect1_hpv6_anal[sex_exposed$HPV6_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV6_ANAL==1]!=2 & sex_exposed$exposure_hpv6_anal[sex_exposed$HPV6_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV6_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv6_anal_partner]==0 & probinfect1_hpv6_anal[sex_exposed$HPV6_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV6_ANAL==1]!=2 & sex_exposed$exposure_hpv6_anal[sex_exposed$HPV6_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV6_ANAL==1]==3 & probinfect1_hpv6_anal[sex_exposed$HPV6_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV6_ANAL==1]!=2 & sex_exposed$exposure_hpv6_anal[sex_exposed$HPV6_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV6_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv6_anal_partner]==1 & probinfect2_hpv6_anal[sex_exposed$HPV6_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV6_ANAL==1]!=2 & sex_exposed$exposure_hpv6_anal[sex_exposed$HPV6_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV6_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv6_anal_partner]==0 & probinfect2_hpv6_anal[sex_exposed$HPV6_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$infection_hpv11_anal[sex_exposed$HPV11_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV11_ANAL==1]!=2 & sex_exposed$exposure_hpv11_anal[sex_exposed$HPV11_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV11_ANAL==1]==1 & probinfect1_hpv11_anal[sex_exposed$HPV11_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV11_ANAL==1]!=2 & sex_exposed$exposure_hpv11_anal[sex_exposed$HPV11_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV11_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv11_anal_ego]==1 & probinfect1_hpv11_anal[sex_exposed$HPV11_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV11_ANAL==1]!=2 & sex_exposed$exposure_hpv11_anal[sex_exposed$HPV11_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV11_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv11_anal_ego]==0 & probinfect1_hpv11_anal[sex_exposed$HPV11_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV11_ANAL==1]!=2 & sex_exposed$exposure_hpv11_anal[sex_exposed$HPV11_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV11_ANAL==1]==3 & probinfect1_hpv11_anal[sex_exposed$HPV11_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV11_ANAL==1]!=2 & sex_exposed$exposure_hpv11_anal[sex_exposed$HPV11_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV11_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv11_anal_ego]==1 & probinfect2_hpv11_anal[sex_exposed$HPV11_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV11_ANAL==1]!=2 & sex_exposed$exposure_hpv11_anal[sex_exposed$HPV11_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV11_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv11_anal_ego]==0 & probinfect2_hpv11_anal[sex_exposed$HPV11_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$partinfection_hpv11_anal[sex_exposed$HPV11_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV11_ANAL==1]!=2 & sex_exposed$exposure_hpv11_anal[sex_exposed$HPV11_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV11_ANAL==1]==1 & probinfect1_hpv11_anal[sex_exposed$HPV11_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV11_ANAL==1]!=2 & sex_exposed$exposure_hpv11_anal[sex_exposed$HPV11_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV11_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv11_anal_partner]==1 & probinfect1_hpv11_anal[sex_exposed$HPV11_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV11_ANAL==1]!=2 & sex_exposed$exposure_hpv11_anal[sex_exposed$HPV11_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV11_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv11_anal_partner]==0 & probinfect1_hpv11_anal[sex_exposed$HPV11_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV11_ANAL==1]!=2 & sex_exposed$exposure_hpv11_anal[sex_exposed$HPV11_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV11_ANAL==1]==3 & probinfect1_hpv11_anal[sex_exposed$HPV11_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV11_ANAL==1]!=2 & sex_exposed$exposure_hpv11_anal[sex_exposed$HPV11_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV11_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv11_anal_partner]==1 & probinfect2_hpv11_anal[sex_exposed$HPV11_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV11_ANAL==1]!=2 & sex_exposed$exposure_hpv11_anal[sex_exposed$HPV11_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV11_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv11_anal_partner]==0 & probinfect2_hpv11_anal[sex_exposed$HPV11_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$infection_hpv16_anal[sex_exposed$HPV16_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV16_ANAL==1]!=2 & sex_exposed$exposure_hpv16_anal[sex_exposed$HPV16_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV16_ANAL==1]==1 & probinfect1_hpv16_anal[sex_exposed$HPV16_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV16_ANAL==1]!=2 & sex_exposed$exposure_hpv16_anal[sex_exposed$HPV16_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV16_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv16_anal_ego]==1 & probinfect1_hpv16_anal[sex_exposed$HPV16_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV16_ANAL==1]!=2 & sex_exposed$exposure_hpv16_anal[sex_exposed$HPV16_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV16_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv16_anal_ego]==0 & probinfect1_hpv16_anal[sex_exposed$HPV16_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV16_ANAL==1]!=2 & sex_exposed$exposure_hpv16_anal[sex_exposed$HPV16_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV16_ANAL==1]==3 & probinfect1_hpv16_anal[sex_exposed$HPV16_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV16_ANAL==1]!=2 & sex_exposed$exposure_hpv16_anal[sex_exposed$HPV16_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV16_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv16_anal_ego]==1 & probinfect2_hpv16_anal[sex_exposed$HPV16_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV16_ANAL==1]!=2 & sex_exposed$exposure_hpv16_anal[sex_exposed$HPV16_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV16_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv16_anal_ego]==0 & probinfect2_hpv16_anal[sex_exposed$HPV16_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$partinfection_hpv16_anal[sex_exposed$HPV16_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV16_ANAL==1]!=2 & sex_exposed$exposure_hpv16_anal[sex_exposed$HPV16_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV16_ANAL==1]==1 & probinfect1_hpv16_anal[sex_exposed$HPV16_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV16_ANAL==1]!=2 & sex_exposed$exposure_hpv16_anal[sex_exposed$HPV16_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV16_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv16_anal_partner]==1 & probinfect1_hpv16_anal[sex_exposed$HPV16_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV16_ANAL==1]!=2 & sex_exposed$exposure_hpv16_anal[sex_exposed$HPV16_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV16_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv16_anal_partner]==0 & probinfect1_hpv16_anal[sex_exposed$HPV16_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV16_ANAL==1]!=2 & sex_exposed$exposure_hpv16_anal[sex_exposed$HPV16_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV16_ANAL==1]==3 & probinfect1_hpv16_anal[sex_exposed$HPV16_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV16_ANAL==1]!=2 & sex_exposed$exposure_hpv16_anal[sex_exposed$HPV16_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV16_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv16_anal_partner]==1 & probinfect2_hpv16_anal[sex_exposed$HPV16_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV16_ANAL==1]!=2 & sex_exposed$exposure_hpv16_anal[sex_exposed$HPV16_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV16_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv16_anal_partner]==0 & probinfect2_hpv16_anal[sex_exposed$HPV16_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$infection_hpv18_anal[sex_exposed$HPV18_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV18_ANAL==1]!=2 & sex_exposed$exposure_hpv18_anal[sex_exposed$HPV18_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV18_ANAL==1]==1 & probinfect1_hpv18_anal[sex_exposed$HPV18_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV18_ANAL==1]!=2 & sex_exposed$exposure_hpv18_anal[sex_exposed$HPV18_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV18_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv18_anal_ego]==1 & probinfect1_hpv18_anal[sex_exposed$HPV18_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV18_ANAL==1]!=2 & sex_exposed$exposure_hpv18_anal[sex_exposed$HPV18_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV18_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv18_anal_ego]==0 & probinfect1_hpv18_anal[sex_exposed$HPV18_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV18_ANAL==1]!=2 & sex_exposed$exposure_hpv18_anal[sex_exposed$HPV18_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV18_ANAL==1]==3 & probinfect1_hpv18_anal[sex_exposed$HPV18_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV18_ANAL==1]!=2 & sex_exposed$exposure_hpv18_anal[sex_exposed$HPV18_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV18_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv18_anal_ego]==1 & probinfect2_hpv18_anal[sex_exposed$HPV18_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV18_ANAL==1]!=2 & sex_exposed$exposure_hpv18_anal[sex_exposed$HPV18_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV18_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv18_anal_ego]==0 & probinfect2_hpv18_anal[sex_exposed$HPV18_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$partinfection_hpv18_anal[sex_exposed$HPV18_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV18_ANAL==1]!=2 & sex_exposed$exposure_hpv18_anal[sex_exposed$HPV18_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV18_ANAL==1]==1 & probinfect1_hpv18_anal[sex_exposed$HPV18_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV18_ANAL==1]!=2 & sex_exposed$exposure_hpv18_anal[sex_exposed$HPV18_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV18_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv18_anal_partner]==1 & probinfect1_hpv18_anal[sex_exposed$HPV18_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV18_ANAL==1]!=2 & sex_exposed$exposure_hpv18_anal[sex_exposed$HPV18_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV18_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv18_anal_partner]==0 & probinfect1_hpv18_anal[sex_exposed$HPV18_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV18_ANAL==1]!=2 & sex_exposed$exposure_hpv18_anal[sex_exposed$HPV18_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV18_ANAL==1]==3 & probinfect1_hpv18_anal[sex_exposed$HPV18_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV18_ANAL==1]!=2 & sex_exposed$exposure_hpv18_anal[sex_exposed$HPV18_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV18_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv18_anal_partner]==1 & probinfect2_hpv18_anal[sex_exposed$HPV18_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV18_ANAL==1]!=2 & sex_exposed$exposure_hpv18_anal[sex_exposed$HPV18_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV18_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv18_anal_partner]==0 & probinfect2_hpv18_anal[sex_exposed$HPV18_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$infection_hpv31_anal[sex_exposed$HPV31_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV31_ANAL==1]!=2 & sex_exposed$exposure_hpv31_anal[sex_exposed$HPV31_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV31_ANAL==1]==1 & probinfect1_hpv31_anal[sex_exposed$HPV31_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV31_ANAL==1]!=2 & sex_exposed$exposure_hpv31_anal[sex_exposed$HPV31_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV31_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv31_anal_ego]==1 & probinfect1_hpv31_anal[sex_exposed$HPV31_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV31_ANAL==1]!=2 & sex_exposed$exposure_hpv31_anal[sex_exposed$HPV31_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV31_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv31_anal_ego]==0 & probinfect1_hpv31_anal[sex_exposed$HPV31_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV31_ANAL==1]!=2 & sex_exposed$exposure_hpv31_anal[sex_exposed$HPV31_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV31_ANAL==1]==3 & probinfect1_hpv31_anal[sex_exposed$HPV31_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV31_ANAL==1]!=2 & sex_exposed$exposure_hpv31_anal[sex_exposed$HPV31_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV31_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv31_anal_ego]==1 & probinfect2_hpv31_anal[sex_exposed$HPV31_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV31_ANAL==1]!=2 & sex_exposed$exposure_hpv31_anal[sex_exposed$HPV31_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV31_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv31_anal_ego]==0 & probinfect2_hpv31_anal[sex_exposed$HPV31_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$partinfection_hpv31_anal[sex_exposed$HPV31_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV31_ANAL==1]!=2 & sex_exposed$exposure_hpv31_anal[sex_exposed$HPV31_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV31_ANAL==1]==1 & probinfect1_hpv31_anal[sex_exposed$HPV31_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV31_ANAL==1]!=2 & sex_exposed$exposure_hpv31_anal[sex_exposed$HPV31_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV31_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv31_anal_partner]==1 & probinfect1_hpv31_anal[sex_exposed$HPV31_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV31_ANAL==1]!=2 & sex_exposed$exposure_hpv31_anal[sex_exposed$HPV31_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV31_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv31_anal_partner]==0 & probinfect1_hpv31_anal[sex_exposed$HPV31_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV31_ANAL==1]!=2 & sex_exposed$exposure_hpv31_anal[sex_exposed$HPV31_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV31_ANAL==1]==3 & probinfect1_hpv31_anal[sex_exposed$HPV31_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV31_ANAL==1]!=2 & sex_exposed$exposure_hpv31_anal[sex_exposed$HPV31_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV31_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv31_anal_partner]==1 & probinfect2_hpv31_anal[sex_exposed$HPV31_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV31_ANAL==1]!=2 & sex_exposed$exposure_hpv31_anal[sex_exposed$HPV31_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV31_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv31_anal_partner]==0 & probinfect2_hpv31_anal[sex_exposed$HPV31_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$infection_hpv33_anal[sex_exposed$HPV33_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV33_ANAL==1]!=2 & sex_exposed$exposure_hpv33_anal[sex_exposed$HPV33_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV33_ANAL==1]==1 & probinfect1_hpv33_anal[sex_exposed$HPV33_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV33_ANAL==1]!=2 & sex_exposed$exposure_hpv33_anal[sex_exposed$HPV33_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV33_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv33_anal_ego]==1 & probinfect1_hpv33_anal[sex_exposed$HPV33_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV33_ANAL==1]!=2 & sex_exposed$exposure_hpv33_anal[sex_exposed$HPV33_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV33_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv33_anal_ego]==0 & probinfect1_hpv33_anal[sex_exposed$HPV33_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV33_ANAL==1]!=2 & sex_exposed$exposure_hpv33_anal[sex_exposed$HPV33_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV33_ANAL==1]==3 & probinfect1_hpv33_anal[sex_exposed$HPV33_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV33_ANAL==1]!=2 & sex_exposed$exposure_hpv33_anal[sex_exposed$HPV33_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV33_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv33_anal_ego]==1 & probinfect2_hpv33_anal[sex_exposed$HPV33_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV33_ANAL==1]!=2 & sex_exposed$exposure_hpv33_anal[sex_exposed$HPV33_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV33_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv33_anal_ego]==0 & probinfect2_hpv33_anal[sex_exposed$HPV33_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$partinfection_hpv33_anal[sex_exposed$HPV33_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV33_ANAL==1]!=2 & sex_exposed$exposure_hpv33_anal[sex_exposed$HPV33_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV33_ANAL==1]==1 & probinfect1_hpv33_anal[sex_exposed$HPV33_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV33_ANAL==1]!=2 & sex_exposed$exposure_hpv33_anal[sex_exposed$HPV33_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV33_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv33_anal_partner]==1 & probinfect1_hpv33_anal[sex_exposed$HPV33_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV33_ANAL==1]!=2 & sex_exposed$exposure_hpv33_anal[sex_exposed$HPV33_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV33_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv33_anal_partner]==0 & probinfect1_hpv33_anal[sex_exposed$HPV33_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV33_ANAL==1]!=2 & sex_exposed$exposure_hpv33_anal[sex_exposed$HPV33_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV33_ANAL==1]==3 & probinfect1_hpv33_anal[sex_exposed$HPV33_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV33_ANAL==1]!=2 & sex_exposed$exposure_hpv33_anal[sex_exposed$HPV33_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV33_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv33_anal_partner]==1 & probinfect2_hpv33_anal[sex_exposed$HPV33_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV33_ANAL==1]!=2 & sex_exposed$exposure_hpv33_anal[sex_exposed$HPV33_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV33_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv33_anal_partner]==0 & probinfect2_hpv33_anal[sex_exposed$HPV33_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$infection_hpv45_anal[sex_exposed$HPV45_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV45_ANAL==1]!=2 & sex_exposed$exposure_hpv45_anal[sex_exposed$HPV45_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV45_ANAL==1]==1 & probinfect1_hpv45_anal[sex_exposed$HPV45_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV45_ANAL==1]!=2 & sex_exposed$exposure_hpv45_anal[sex_exposed$HPV45_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV45_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv45_anal_ego]==1 & probinfect1_hpv45_anal[sex_exposed$HPV45_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV45_ANAL==1]!=2 & sex_exposed$exposure_hpv45_anal[sex_exposed$HPV45_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV45_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv45_anal_ego]==0 & probinfect1_hpv45_anal[sex_exposed$HPV45_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV45_ANAL==1]!=2 & sex_exposed$exposure_hpv45_anal[sex_exposed$HPV45_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV45_ANAL==1]==3 & probinfect1_hpv45_anal[sex_exposed$HPV45_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV45_ANAL==1]!=2 & sex_exposed$exposure_hpv45_anal[sex_exposed$HPV45_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV45_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv45_anal_ego]==1 & probinfect2_hpv45_anal[sex_exposed$HPV45_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV45_ANAL==1]!=2 & sex_exposed$exposure_hpv45_anal[sex_exposed$HPV45_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV45_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv45_anal_ego]==0 & probinfect2_hpv45_anal[sex_exposed$HPV45_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$partinfection_hpv45_anal[sex_exposed$HPV45_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV45_ANAL==1]!=2 & sex_exposed$exposure_hpv45_anal[sex_exposed$HPV45_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV45_ANAL==1]==1 & probinfect1_hpv45_anal[sex_exposed$HPV45_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV45_ANAL==1]!=2 & sex_exposed$exposure_hpv45_anal[sex_exposed$HPV45_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV45_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv45_anal_partner]==1 & probinfect1_hpv45_anal[sex_exposed$HPV45_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV45_ANAL==1]!=2 & sex_exposed$exposure_hpv45_anal[sex_exposed$HPV45_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV45_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv45_anal_partner]==0 & probinfect1_hpv45_anal[sex_exposed$HPV45_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV45_ANAL==1]!=2 & sex_exposed$exposure_hpv45_anal[sex_exposed$HPV45_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV45_ANAL==1]==3 & probinfect1_hpv45_anal[sex_exposed$HPV45_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV45_ANAL==1]!=2 & sex_exposed$exposure_hpv45_anal[sex_exposed$HPV45_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV45_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv45_anal_partner]==1 & probinfect2_hpv45_anal[sex_exposed$HPV45_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV45_ANAL==1]!=2 & sex_exposed$exposure_hpv45_anal[sex_exposed$HPV45_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV45_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv45_anal_partner]==0 & probinfect2_hpv45_anal[sex_exposed$HPV45_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$infection_hpv52_anal[sex_exposed$HPV52_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV52_ANAL==1]!=2 & sex_exposed$exposure_hpv52_anal[sex_exposed$HPV52_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV52_ANAL==1]==1 & probinfect1_hpv52_anal[sex_exposed$HPV52_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV52_ANAL==1]!=2 & sex_exposed$exposure_hpv52_anal[sex_exposed$HPV52_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV52_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv52_anal_ego]==1 & probinfect1_hpv52_anal[sex_exposed$HPV52_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV52_ANAL==1]!=2 & sex_exposed$exposure_hpv52_anal[sex_exposed$HPV52_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV52_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv52_anal_ego]==0 & probinfect1_hpv52_anal[sex_exposed$HPV52_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV52_ANAL==1]!=2 & sex_exposed$exposure_hpv52_anal[sex_exposed$HPV52_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV52_ANAL==1]==3 & probinfect1_hpv52_anal[sex_exposed$HPV52_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV52_ANAL==1]!=2 & sex_exposed$exposure_hpv52_anal[sex_exposed$HPV52_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV52_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv52_anal_ego]==1 & probinfect2_hpv52_anal[sex_exposed$HPV52_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV52_ANAL==1]!=2 & sex_exposed$exposure_hpv52_anal[sex_exposed$HPV52_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV52_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv52_anal_ego]==0 & probinfect2_hpv52_anal[sex_exposed$HPV52_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$partinfection_hpv52_anal[sex_exposed$HPV52_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV52_ANAL==1]!=2 & sex_exposed$exposure_hpv52_anal[sex_exposed$HPV52_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV52_ANAL==1]==1 & probinfect1_hpv52_anal[sex_exposed$HPV52_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV52_ANAL==1]!=2 & sex_exposed$exposure_hpv52_anal[sex_exposed$HPV52_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV52_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv52_anal_partner]==1 & probinfect1_hpv52_anal[sex_exposed$HPV52_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV52_ANAL==1]!=2 & sex_exposed$exposure_hpv52_anal[sex_exposed$HPV52_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV52_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv52_anal_partner]==0 & probinfect1_hpv52_anal[sex_exposed$HPV52_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV52_ANAL==1]!=2 & sex_exposed$exposure_hpv52_anal[sex_exposed$HPV52_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV52_ANAL==1]==3 & probinfect1_hpv52_anal[sex_exposed$HPV52_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV52_ANAL==1]!=2 & sex_exposed$exposure_hpv52_anal[sex_exposed$HPV52_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV52_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv52_anal_partner]==1 & probinfect2_hpv52_anal[sex_exposed$HPV52_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV52_ANAL==1]!=2 & sex_exposed$exposure_hpv52_anal[sex_exposed$HPV52_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV52_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv52_anal_partner]==0 & probinfect2_hpv52_anal[sex_exposed$HPV52_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$infection_hpv58_anal[sex_exposed$HPV58_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV58_ANAL==1]!=2 & sex_exposed$exposure_hpv58_anal[sex_exposed$HPV58_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV58_ANAL==1]==1 & probinfect1_hpv58_anal[sex_exposed$HPV58_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV58_ANAL==1]!=2 & sex_exposed$exposure_hpv58_anal[sex_exposed$HPV58_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV58_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv58_anal_ego]==1 & probinfect1_hpv58_anal[sex_exposed$HPV58_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV58_ANAL==1]!=2 & sex_exposed$exposure_hpv58_anal[sex_exposed$HPV58_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV58_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv58_anal_ego]==0 & probinfect1_hpv58_anal[sex_exposed$HPV58_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV58_ANAL==1]!=2 & sex_exposed$exposure_hpv58_anal[sex_exposed$HPV58_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV58_ANAL==1]==3 & probinfect1_hpv58_anal[sex_exposed$HPV58_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV58_ANAL==1]!=2 & sex_exposed$exposure_hpv58_anal[sex_exposed$HPV58_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV58_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv58_anal_ego]==1 & probinfect2_hpv58_anal[sex_exposed$HPV58_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV58_ANAL==1]!=2 & sex_exposed$exposure_hpv58_anal[sex_exposed$HPV58_ANAL==1]==1 & sex_exposed$simAI[sex_exposed$HPV58_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv58_anal_ego]==0 & probinfect2_hpv58_anal[sex_exposed$HPV58_ANAL==1] <= 0.071, 1, 0))))))
          sex_exposed$partinfection_hpv58_anal[sex_exposed$HPV58_ANAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV58_ANAL==1]!=2 & sex_exposed$exposure_hpv58_anal[sex_exposed$HPV58_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV58_ANAL==1]==1 & probinfect1_hpv58_anal[sex_exposed$HPV58_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV58_ANAL==1]!=2 & sex_exposed$exposure_hpv58_anal[sex_exposed$HPV58_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV58_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv58_anal_partner]==1 & probinfect1_hpv58_anal[sex_exposed$HPV58_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV58_ANAL==1]!=2 & sex_exposed$exposure_hpv58_anal[sex_exposed$HPV58_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV58_ANAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv58_anal_partner]==0 & probinfect1_hpv58_anal[sex_exposed$HPV58_ANAL==1] <= 0.071, 1, ifelse(sex_exposed$sex[sex_exposed$HPV58_ANAL==1]!=2 & sex_exposed$exposure_hpv58_anal[sex_exposed$HPV58_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV58_ANAL==1]==3 & probinfect1_hpv58_anal[sex_exposed$HPV58_ANAL==1] <= 0.546, 1, ifelse(sex_exposed$sex[sex_exposed$HPV58_ANAL==1]!=2 & sex_exposed$exposure_hpv58_anal[sex_exposed$HPV58_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV58_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv58_anal_partner]==1 & probinfect2_hpv58_anal[sex_exposed$HPV58_ANAL==1] <= 0.038, 1, ifelse(sex_exposed$sex[sex_exposed$HPV58_ANAL==1]!=2 & sex_exposed$exposure_hpv58_anal[sex_exposed$HPV58_ANAL==1]==2 & sex_exposed$simAI[sex_exposed$HPV58_ANAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv58_anal_partner]==0 & probinfect2_hpv58_anal[sex_exposed$HPV58_ANAL==1] <= 0.071, 1, 0))))))
          
          #hpv type specific oral, ego then partner, 1=new oral HPV infection, 2=new genital HPV infection, 3=new oral & genital HPV infection
          sex_exposed$infection_hpv6_oral[sex_exposed$HPV6_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==1 & probinfect1_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv6_oral_ego]==1 & probinfect1_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv6_oral_ego]==0 & probinfect1_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==3 & probinfect1_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv6_oral_ego]==1 & probinfect2_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv6_oral_ego]==0 & probinfect2_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==3 & probinfect1_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv6_oral_ego]==1 & probinfect2_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv6_oral_ego]==0 & probinfect2_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$partinfection_hpv6_oral[sex_exposed$HPV6_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==2 & probinfect1_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv6_oral_partner]==1 & probinfect1_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv6_oral_partner]==0 & probinfect1_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==3 & probinfect1_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv6_oral_partner]==1 & probinfect2_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv6_oral_partner]==0 & probinfect2_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==3 & probinfect1_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv6_oral_partner]==1 & probinfect2_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$exposure_hpv6_oral[sex_exposed$HPV6_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV6_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv6_oral_partner]==0 & probinfect2_hpv6_oral[sex_exposed$HPV6_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$infection_hpv11_oral[sex_exposed$HPV11_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==1 & probinfect1_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv11_oral_ego]==1 & probinfect1_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv11_oral_ego]==0 & probinfect1_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==3 & probinfect1_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv11_oral_ego]==1 & probinfect2_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv11_oral_ego]==0 & probinfect2_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==3 & probinfect1_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv11_oral_ego]==1 & probinfect2_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv11_oral_ego]==0 & probinfect2_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$partinfection_hpv11_oral[sex_exposed$HPV11_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==2 & probinfect1_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv11_oral_partner]==1 & probinfect1_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv11_oral_partner]==0 & probinfect1_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==3 & probinfect1_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv11_oral_partner]==1 & probinfect2_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv11_oral_partner]==0 & probinfect2_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==3 & probinfect1_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv11_oral_partner]==1 & probinfect2_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$exposure_hpv11_oral[sex_exposed$HPV11_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV11_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv11_oral_partner]==0 & probinfect2_hpv11_oral[sex_exposed$HPV11_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$infection_hpv16_oral[sex_exposed$HPV16_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==1 & probinfect1_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv16_oral_ego]==1 & probinfect1_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv16_oral_ego]==0 & probinfect1_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==3 & probinfect1_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv16_oral_ego]==1 & probinfect2_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv16_oral_ego]==0 & probinfect2_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==3 & probinfect1_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv16_oral_ego]==1 & probinfect2_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv16_oral_ego]==0 & probinfect2_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$partinfection_hpv16_oral[sex_exposed$HPV16_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==2 & probinfect1_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv16_oral_partner]==1 & probinfect1_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv16_oral_partner]==0 & probinfect1_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==3 & probinfect1_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv16_oral_partner]==1 & probinfect2_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv16_oral_partner]==0 & probinfect2_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==3 & probinfect1_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv16_oral_partner]==1 & probinfect2_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$exposure_hpv16_oral[sex_exposed$HPV16_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV16_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv16_oral_partner]==0 & probinfect2_hpv16_oral[sex_exposed$HPV16_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$infection_hpv18_oral[sex_exposed$HPV18_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==1 & probinfect1_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv18_oral_ego]==1 & probinfect1_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv18_oral_ego]==0 & probinfect1_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==3 & probinfect1_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv18_oral_ego]==1 & probinfect2_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv18_oral_ego]==0 & probinfect2_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==3 & probinfect1_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv18_oral_ego]==1 & probinfect2_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv18_oral_ego]==0 & probinfect2_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$partinfection_hpv18_oral[sex_exposed$HPV18_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==2 & probinfect1_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv18_oral_partner]==1 & probinfect1_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv18_oral_partner]==0 & probinfect1_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==3 & probinfect1_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv18_oral_partner]==1 & probinfect2_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv18_oral_partner]==0 & probinfect2_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==3 & probinfect1_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv18_oral_partner]==1 & probinfect2_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$exposure_hpv18_oral[sex_exposed$HPV18_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV18_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv18_oral_partner]==0 & probinfect2_hpv18_oral[sex_exposed$HPV18_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$infection_hpv31_oral[sex_exposed$HPV31_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==1 & probinfect1_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv31_oral_ego]==1 & probinfect1_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv31_oral_ego]==0 & probinfect1_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==3 & probinfect1_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv31_oral_ego]==1 & probinfect2_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv31_oral_ego]==0 & probinfect2_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==3 & probinfect1_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv31_oral_ego]==1 & probinfect2_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv31_oral_ego]==0 & probinfect2_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$partinfection_hpv31_oral[sex_exposed$HPV31_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==2 & probinfect1_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv31_oral_partner]==1 & probinfect1_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv31_oral_partner]==0 & probinfect1_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==3 & probinfect1_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv31_oral_partner]==1 & probinfect2_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv31_oral_partner]==0 & probinfect2_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==3 & probinfect1_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv31_oral_partner]==1 & probinfect2_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$exposure_hpv31_oral[sex_exposed$HPV31_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV31_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv31_oral_partner]==0 & probinfect2_hpv31_oral[sex_exposed$HPV31_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$infection_hpv33_oral[sex_exposed$HPV33_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==1 & probinfect1_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv33_oral_ego]==1 & probinfect1_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv33_oral_ego]==0 & probinfect1_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==3 & probinfect1_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv33_oral_ego]==1 & probinfect2_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv33_oral_ego]==0 & probinfect2_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==3 & probinfect1_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv33_oral_ego]==1 & probinfect2_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv33_oral_ego]==0 & probinfect2_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$partinfection_hpv33_oral[sex_exposed$HPV33_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==2 & probinfect1_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv33_oral_partner]==1 & probinfect1_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv33_oral_partner]==0 & probinfect1_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==3 & probinfect1_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv33_oral_partner]==1 & probinfect2_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv33_oral_partner]==0 & probinfect2_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==3 & probinfect1_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv33_oral_partner]==1 & probinfect2_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$exposure_hpv33_oral[sex_exposed$HPV33_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV33_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv33_oral_partner]==0 & probinfect2_hpv33_oral[sex_exposed$HPV33_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$infection_hpv45_oral[sex_exposed$HPV45_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==1 & probinfect1_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv45_oral_ego]==1 & probinfect1_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv45_oral_ego]==0 & probinfect1_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==3 & probinfect1_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv45_oral_ego]==1 & probinfect2_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv45_oral_ego]==0 & probinfect2_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==3 & probinfect1_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv45_oral_ego]==1 & probinfect2_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv45_oral_ego]==0 & probinfect2_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$partinfection_hpv45_oral[sex_exposed$HPV45_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==2 & probinfect1_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv45_oral_partner]==1 & probinfect1_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv45_oral_partner]==0 & probinfect1_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==3 & probinfect1_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv45_oral_partner]==1 & probinfect2_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv45_oral_partner]==0 & probinfect2_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==3 & probinfect1_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv45_oral_partner]==1 & probinfect2_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$exposure_hpv45_oral[sex_exposed$HPV45_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV45_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv45_oral_partner]==0 & probinfect2_hpv45_oral[sex_exposed$HPV45_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$infection_hpv52_oral[sex_exposed$HPV52_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==1 & probinfect1_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv52_oral_ego]==1 & probinfect1_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv52_oral_ego]==0 & probinfect1_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==3 & probinfect1_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv52_oral_ego]==1 & probinfect2_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv52_oral_ego]==0 & probinfect2_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==3 & probinfect1_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv52_oral_ego]==1 & probinfect2_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv52_oral_ego]==0 & probinfect2_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$partinfection_hpv52_oral[sex_exposed$HPV52_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==2 & probinfect1_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv52_oral_partner]==1 & probinfect1_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv52_oral_partner]==0 & probinfect1_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==3 & probinfect1_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv52_oral_partner]==1 & probinfect2_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv52_oral_partner]==0 & probinfect2_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==3 & probinfect1_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv52_oral_partner]==1 & probinfect2_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$exposure_hpv52_oral[sex_exposed$HPV52_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV52_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv52_oral_partner]==0 & probinfect2_hpv52_oral[sex_exposed$HPV52_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$infection_hpv58_oral[sex_exposed$HPV58_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==1 & probinfect1_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv58_oral_ego]==1 & probinfect1_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==2 & abm_sims[[current_data]]$CIRC[exposed_hpv58_oral_ego]==0 & probinfect1_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==3 & probinfect1_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv58_oral_ego]==1 & probinfect2_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv58_oral_ego]==0 & probinfect2_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==3 & probinfect1_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv58_oral_ego]==1 & probinfect2_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=2 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv58_oral_ego]==0 & probinfect2_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.014, 2, 0)))))))
          sex_exposed$partinfection_hpv58_oral[sex_exposed$HPV58_ORAL==1] = ifelse(sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==2 & probinfect1_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv58_oral_partner]==1 & probinfect1_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==1 & abm_sims[[current_data]]$CIRC[exposed_hpv58_oral_partner]==0 & probinfect1_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.014, 2, ifelse((sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==3 & probinfect1_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.032) & ((sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv58_oral_partner]==1 & probinfect2_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.008) | (sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv58_oral_partner]==0 & probinfect2_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.014)), 3, ifelse(sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==3 & probinfect1_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.032, 1, ifelse(sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv58_oral_partner]==1 & probinfect2_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.008, 2, ifelse(sex_exposed$sex[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$exposure_hpv58_oral[sex_exposed$HPV58_ORAL==1]!=1 & sex_exposed$simOI[sex_exposed$HPV58_ORAL==1]==3 & abm_sims[[current_data]]$CIRC[exposed_hpv58_oral_partner]==0 & probinfect2_hpv58_oral[sex_exposed$HPV58_ORAL==1] <= 0.014, 2, 0)))))))
          
          #determine infection under seroadaption scenarios
          if (unique(abm_sims[[current_data]]$SERO)==1)
          {
            #check if avoiding sex: both insertive seroadapters, or serosorter (ego/partner) and partner is aware HIV+
            avoid = ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==3 & abm_sims[[current_data]]$HIV_STATUS[exposed_hiv_partner]==3, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==3 & abm_sims[[current_data]]$HIV_STATUS[exposed_hiv_partner]==3, 1, 0)))
            partavoid = ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==3 & abm_sims[[current_data]]$HIV_STATUS[exposed_hiv_ego]==3, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==3 & abm_sims[[current_data]]$HIV_STATUS[exposed_hiv_ego]==3, 1, 0)))
            
            #seroadaption strategies for HIV
            seroblock_hiv = ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0010, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0059, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0101, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0569, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0010, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0059, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0101, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0569, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0134, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.1284, 1, 0))))))))))
            partseroblock_hiv = ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0010, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0059, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0101, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0569, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0010, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0059, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0101, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0569, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==1 & probinfect1_hiv[sex_exposed$HIV==1] > 0.0134, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego]==2 & probinfect1_hiv[sex_exposed$HIV==1] > 0.1284, 1, 0))))))))))
            
            #check if infection blocked for HIV
            sex_exposed$seroct_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==1, 1, ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & avoid==0 & seroblock_hiv==1, 1, 0))
            sex_exposed$partseroct_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==1, 1, ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & partavoid==0 & partseroblock_hiv==1, 1, 0))
            
            #seroadaption strategies for HIV, HPV implications; NAs (and possibly 0s) get returned for HIV infections without correspoding HPV infection
            #seroblock_hpv = ifelse(sex_exposed$infection_hpv_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & probinfect1_hpv_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            #partseroblock_hpv = ifelse(sex_exposed$partinfection_hpv_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & probinfect1_hpv_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            seroblock_hpv6 = ifelse(sex_exposed$infection_hpv6_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv6_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv6_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv6_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv6_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv6_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv6_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv6_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv6_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & probinfect1_hpv6_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            partseroblock_hpv6 = ifelse(sex_exposed$partinfection_hpv6_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv6_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv6_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv6_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv6_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv6_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv6_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv6_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv6_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & probinfect1_hpv6_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            seroblock_hpv11 = ifelse(sex_exposed$infection_hpv11_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv11_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv11_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv11_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv11_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv11_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv11_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv11_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv11_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & probinfect1_hpv11_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            partseroblock_hpv11 = ifelse(sex_exposed$partinfection_hpv11_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv11_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv11_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv11_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv11_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv11_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv11_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv11_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv11_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & probinfect1_hpv11_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            seroblock_hpv16 = ifelse(sex_exposed$infection_hpv16_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv16_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv16_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv16_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv16_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv16_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv16_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv16_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv16_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & probinfect1_hpv16_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            partseroblock_hpv16 = ifelse(sex_exposed$partinfection_hpv16_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv16_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv16_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv16_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv16_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv16_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv16_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv16_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv16_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & probinfect1_hpv16_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            seroblock_hpv18 = ifelse(sex_exposed$infection_hpv18_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv18_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv18_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv18_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv18_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv18_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv18_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv18_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv18_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & probinfect1_hpv18_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            partseroblock_hpv18 = ifelse(sex_exposed$partinfection_hpv18_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv18_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv18_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv18_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv18_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv18_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv18_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv18_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv18_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & probinfect1_hpv18_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            seroblock_hpv31 = ifelse(sex_exposed$infection_hpv31_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv31_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv31_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv31_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv31_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv31_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv31_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv31_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv31_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & probinfect1_hpv31_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            partseroblock_hpv31 = ifelse(sex_exposed$partinfection_hpv31_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv31_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv31_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv31_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv31_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv31_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv31_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv31_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv31_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & probinfect1_hpv31_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            seroblock_hpv33 = ifelse(sex_exposed$infection_hpv33_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv33_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv33_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv33_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv33_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv33_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv33_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv33_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv33_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & probinfect1_hpv33_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            partseroblock_hpv33 = ifelse(sex_exposed$partinfection_hpv33_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv33_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv33_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv33_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv33_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv33_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv33_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv33_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv33_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & probinfect1_hpv33_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            seroblock_hpv45 = ifelse(sex_exposed$infection_hpv45_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv45_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv45_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv45_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv45_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv45_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv45_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv45_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv45_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & probinfect1_hpv45_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            partseroblock_hpv45 = ifelse(sex_exposed$partinfection_hpv45_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv45_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv45_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv45_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv45_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv45_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv45_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv45_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv45_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & probinfect1_hpv45_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            seroblock_hpv52 = ifelse(sex_exposed$infection_hpv52_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv52_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv52_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv52_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv52_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv52_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv52_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv52_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv52_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & probinfect1_hpv52_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            partseroblock_hpv52 = ifelse(sex_exposed$partinfection_hpv52_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv52_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv52_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv52_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv52_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv52_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv52_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv52_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv52_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & probinfect1_hpv52_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            seroblock_hpv58 = ifelse(sex_exposed$infection_hpv58_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv58_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv58_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv58_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv58_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==1 & probinfect1_hpv58_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$infection_hpv58_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_ego]==2 & probinfect1_hpv58_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$infection_hpv58_anal[sex_exposed$HIV==1]==1 & avoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & probinfect1_hpv58_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            partseroblock_hpv58 = ifelse(sex_exposed$partinfection_hpv58_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv58_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv58_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_partner]==1 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv58_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv58_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==1 & probinfect1_hpv58_anal[sex_exposed$HIV==1] > 0.038, 1, ifelse(sex_exposed$partinfection_hpv58_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==2 & abm_sims[[current_data]]$CIRC[exposed_hiv_partner]==2 & probinfect1_hpv58_anal[sex_exposed$HIV==1] > 0.071, 1, ifelse(sex_exposed$partinfection_hpv58_anal[sex_exposed$HIV==1]==1 & partavoid==0 & abm_sims[[current_data]]$SEROTYPE[exposed_hiv_ego]==1 & probinfect1_hpv58_anal[sex_exposed$HIV==1] > 0.546, 1, 0)))))
            
            #check if infection blocked for HPV anal, removing NAs/0s from earlier that are not coinfections
            #if (sum(sex_exposed$HIV_HPV_ANAL)>0) { sex_exposed$seroct_hpv[sex_exposed$HIV_HPV_ANAL==1] = (ifelse(sex_exposed$infection_hpv_anal[sex_exposed$HIV==1]==1 & avoid==1, 1, ifelse(sex_exposed$infection_hpv_anal[sex_exposed$HIV==1]==1 & avoid==0 & seroblock_hpv==1, 1, 0)))[which(!is.na(sex_exposed$infection_hpv_anal[sex_exposed$HIV==1]))] }
            #if (sum(sex_exposed$HIV_HPV_ANAL)>0) { sex_exposed$partseroct_hpv[sex_exposed$HIV_HPV_ANAL==1] = (ifelse(sex_exposed$partinfection_hpv_anal[sex_exposed$HIV==1]==1 & partavoid==1, 1, ifelse(sex_exposed$partinfection_hpv_anal[sex_exposed$HIV==1]==1 & partavoid==0 & partseroblock_hpv==1, 1, 0)))[which(!is.na(sex_exposed$partinfection_hpv_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV6_ANAL)>0) { sex_exposed$seroct_hpv6[sex_exposed$HIV_HPV6_ANAL==1] = (ifelse(sex_exposed$infection_hpv6_anal[sex_exposed$HIV==1]==1 & avoid==1, 1, ifelse(sex_exposed$infection_hpv6_anal[sex_exposed$HIV==1]==1 & avoid==0 & seroblock_hpv6==1, 1, 0)))[which(!is.na(sex_exposed$infection_hpv6_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV6_ANAL)>0) { sex_exposed$partseroct_hpv6[sex_exposed$HIV_HPV6_ANAL==1] = (ifelse(sex_exposed$partinfection_hpv6_anal[sex_exposed$HIV==1]==1 & partavoid==1, 1, ifelse(sex_exposed$partinfection_hpv6_anal[sex_exposed$HIV==1]==1 & partavoid==0 & partseroblock_hpv6==1, 1, 0)))[which(!is.na(sex_exposed$partinfection_hpv6_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV11_ANAL)>0) { sex_exposed$seroct_hpv11[sex_exposed$HIV_HPV11_ANAL==1] = (ifelse(sex_exposed$infection_hpv11_anal[sex_exposed$HIV==1]==1 & avoid==1, 1, ifelse(sex_exposed$infection_hpv11_anal[sex_exposed$HIV==1]==1 & avoid==0 & seroblock_hpv11==1, 1, 0)))[which(!is.na(sex_exposed$infection_hpv11_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV11_ANAL)>0) { sex_exposed$partseroct_hpv11[sex_exposed$HIV_HPV11_ANAL==1] = (ifelse(sex_exposed$partinfection_hpv11_anal[sex_exposed$HIV==1]==1 & partavoid==1, 1, ifelse(sex_exposed$partinfection_hpv11_anal[sex_exposed$HIV==1]==1 & partavoid==0 & partseroblock_hpv11==1, 1, 0)))[which(!is.na(sex_exposed$partinfection_hpv11_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV16_ANAL)>0) { sex_exposed$seroct_hpv16[sex_exposed$HIV_HPV16_ANAL==1] = (ifelse(sex_exposed$infection_hpv16_anal[sex_exposed$HIV==1]==1 & avoid==1, 1, ifelse(sex_exposed$infection_hpv16_anal[sex_exposed$HIV==1]==1 & avoid==0 & seroblock_hpv16==1, 1, 0)))[which(!is.na(sex_exposed$infection_hpv16_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV16_ANAL)>0) { sex_exposed$partseroct_hpv16[sex_exposed$HIV_HPV16_ANAL==1] = (ifelse(sex_exposed$partinfection_hpv16_anal[sex_exposed$HIV==1]==1 & partavoid==1, 1, ifelse(sex_exposed$partinfection_hpv16_anal[sex_exposed$HIV==1]==1 & partavoid==0 & partseroblock_hpv16==1, 1, 0)))[which(!is.na(sex_exposed$partinfection_hpv16_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV18_ANAL)>0) { sex_exposed$seroct_hpv18[sex_exposed$HIV_HPV18_ANAL==1] = (ifelse(sex_exposed$infection_hpv18_anal[sex_exposed$HIV==1]==1 & avoid==1, 1, ifelse(sex_exposed$infection_hpv18_anal[sex_exposed$HIV==1]==1 & avoid==0 & seroblock_hpv18==1, 1, 0)))[which(!is.na(sex_exposed$infection_hpv18_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV18_ANAL)>0) { sex_exposed$partseroct_hpv18[sex_exposed$HIV_HPV18_ANAL==1] = (ifelse(sex_exposed$partinfection_hpv18_anal[sex_exposed$HIV==1]==1 & partavoid==1, 1, ifelse(sex_exposed$partinfection_hpv18_anal[sex_exposed$HIV==1]==1 & partavoid==0 & partseroblock_hpv18==1, 1, 0)))[which(!is.na(sex_exposed$partinfection_hpv18_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV31_ANAL)>0) { sex_exposed$seroct_hpv31[sex_exposed$HIV_HPV31_ANAL==1] = (ifelse(sex_exposed$infection_hpv31_anal[sex_exposed$HIV==1]==1 & avoid==1, 1, ifelse(sex_exposed$infection_hpv31_anal[sex_exposed$HIV==1]==1 & avoid==0 & seroblock_hpv31==1, 1, 0)))[which(!is.na(sex_exposed$infection_hpv31_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV31_ANAL)>0) { sex_exposed$partseroct_hpv31[sex_exposed$HIV_HPV31_ANAL==1] = (ifelse(sex_exposed$partinfection_hpv31_anal[sex_exposed$HIV==1]==1 & partavoid==1, 1, ifelse(sex_exposed$partinfection_hpv31_anal[sex_exposed$HIV==1]==1 & partavoid==0 & partseroblock_hpv31==1, 1, 0)))[which(!is.na(sex_exposed$partinfection_hpv31_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV33_ANAL)>0) { sex_exposed$seroct_hpv33[sex_exposed$HIV_HPV33_ANAL==1] = (ifelse(sex_exposed$infection_hpv33_anal[sex_exposed$HIV==1]==1 & avoid==1, 1, ifelse(sex_exposed$infection_hpv33_anal[sex_exposed$HIV==1]==1 & avoid==0 & seroblock_hpv33==1, 1, 0)))[which(!is.na(sex_exposed$infection_hpv33_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV33_ANAL)>0) { sex_exposed$partseroct_hpv33[sex_exposed$HIV_HPV33_ANAL==1] = (ifelse(sex_exposed$partinfection_hpv33_anal[sex_exposed$HIV==1]==1 & partavoid==1, 1, ifelse(sex_exposed$partinfection_hpv33_anal[sex_exposed$HIV==1]==1 & partavoid==0 & partseroblock_hpv33==1, 1, 0)))[which(!is.na(sex_exposed$partinfection_hpv33_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV45_ANAL)>0) { sex_exposed$seroct_hpv45[sex_exposed$HIV_HPV45_ANAL==1] = (ifelse(sex_exposed$infection_hpv45_anal[sex_exposed$HIV==1]==1 & avoid==1, 1, ifelse(sex_exposed$infection_hpv45_anal[sex_exposed$HIV==1]==1 & avoid==0 & seroblock_hpv45==1, 1, 0)))[which(!is.na(sex_exposed$infection_hpv45_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV45_ANAL)>0) { sex_exposed$partseroct_hpv45[sex_exposed$HIV_HPV45_ANAL==1] = (ifelse(sex_exposed$partinfection_hpv45_anal[sex_exposed$HIV==1]==1 & partavoid==1, 1, ifelse(sex_exposed$partinfection_hpv45_anal[sex_exposed$HIV==1]==1 & partavoid==0 & partseroblock_hpv45==1, 1, 0)))[which(!is.na(sex_exposed$partinfection_hpv45_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV52_ANAL)>0) { sex_exposed$seroct_hpv52[sex_exposed$HIV_HPV52_ANAL==1] = (ifelse(sex_exposed$infection_hpv52_anal[sex_exposed$HIV==1]==1 & avoid==1, 1, ifelse(sex_exposed$infection_hpv52_anal[sex_exposed$HIV==1]==1 & avoid==0 & seroblock_hpv52==1, 1, 0)))[which(!is.na(sex_exposed$infection_hpv52_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV52_ANAL)>0) { sex_exposed$partseroct_hpv52[sex_exposed$HIV_HPV52_ANAL==1] = (ifelse(sex_exposed$partinfection_hpv52_anal[sex_exposed$HIV==1]==1 & partavoid==1, 1, ifelse(sex_exposed$partinfection_hpv52_anal[sex_exposed$HIV==1]==1 & partavoid==0 & partseroblock_hpv52==1, 1, 0)))[which(!is.na(sex_exposed$partinfection_hpv52_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV58_ANAL)>0) { sex_exposed$seroct_hpv58[sex_exposed$HIV_HPV58_ANAL==1] = (ifelse(sex_exposed$infection_hpv58_anal[sex_exposed$HIV==1]==1 & avoid==1, 1, ifelse(sex_exposed$infection_hpv58_anal[sex_exposed$HIV==1]==1 & avoid==0 & seroblock_hpv58==1, 1, 0)))[which(!is.na(sex_exposed$infection_hpv58_anal[sex_exposed$HIV==1]))] }
            if (sum(sex_exposed$HIV_HPV58_ANAL)>0) { sex_exposed$partseroct_hpv58[sex_exposed$HIV_HPV58_ANAL==1] = (ifelse(sex_exposed$partinfection_hpv58_anal[sex_exposed$HIV==1]==1 & partavoid==1, 1, ifelse(sex_exposed$partinfection_hpv58_anal[sex_exposed$HIV==1]==1 & partavoid==0 & partseroblock_hpv58==1, 1, 0)))[which(!is.na(sex_exposed$partinfection_hpv58_anal[sex_exposed$HIV==1]))] }
            
            #decrement sex counter for those who have avoided sex (because it later gets incremented irrespective of seroadaption outcome)
            abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Ego[sex_exposed$HIV==1]] = ifelse(avoid==1, abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Ego[sex_exposed$HIV==1]] - 1, abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Ego[sex_exposed$HIV==1]])
            abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Partner[sex_exposed$HIV==1]] = ifelse(partavoid==1, abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Partner[sex_exposed$HIV==1]] - 1, abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Partner[sex_exposed$HIV==1]])
            
            rm(avoid,partavoid,seroblock_hiv,partseroblock_hiv,seroblock_hpv6,partseroblock_hpv6,seroblock_hpv11,partseroblock_hpv11,seroblock_hpv16,partseroblock_hpv16,seroblock_hpv18,partseroblock_hpv18,seroblock_hpv31,partseroblock_hpv31,seroblock_hpv33,partseroblock_hpv33,seroblock_hpv45,partseroblock_hpv45,seroblock_hpv52,partseroblock_hpv52,seroblock_hpv58,partseroblock_hpv58)
          }
          
          #determine infection under condom scenarios
          if (unique(abm_sims[[current_data]]$COND)==1)
          {
            #check if a condom was worn by either partner
            randcond = runif(nrow(sex_exposed),0,1)
            condom_worn_anal = (randcond<=abm_sims[[current_data]]$PROBCOND[sex_exposed$Ego]) | (randcond<=abm_sims[[current_data]]$PROBCOND[sex_exposed$Partner])
            condom_worn_oral = (randcond<=0.04)
            
            #check if the condom failed
            randfail = runif(nrow(sex_exposed),0,1)
            condom_success_anal = condom_worn_anal & (randfail<=0.705)
            condom_success_oral = condom_worn_oral & (randfail<=0.705)
            
            #check if infection blocked for HIV
            sex_exposed$condct_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & condom_success_anal[sex_exposed$HIV==1], 1, 0)
            sex_exposed$partcondct_hiv[sex_exposed$HIV==1] = ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & condom_success_anal[sex_exposed$HIV==1], 1, 0)
            
            #sensitivity analysis for lower condom effectiveness for HPV
            #condom_success_anal = condom_worn_anal & (randfail<=0.5)
            #condom_success_oral = condom_worn_oral & (randfail<=0.5)
            
            #check if infection blocked for HPV anal
            #sex_exposed$condct_hpv_anal[sex_exposed$HPV_ANAL==1] = ifelse(sex_exposed$infection_hpv_anal[sex_exposed$HPV_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV_ANAL==1], 1, 0)
            #sex_exposed$partcondct_hpv_anal[sex_exposed$HPV_ANAL==1] = ifelse(sex_exposed$partinfection_hpv_anal[sex_exposed$HPV_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV_ANAL==1], 1, 0)
            sex_exposed$condct_hpv6_anal[sex_exposed$HPV6_ANAL==1] = ifelse(sex_exposed$infection_hpv6_anal[sex_exposed$HPV6_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV6_ANAL==1], 1, 0)
            sex_exposed$partcondct_hpv6_anal[sex_exposed$HPV6_ANAL==1] = ifelse(sex_exposed$partinfection_hpv6_anal[sex_exposed$HPV6_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV6_ANAL==1], 1, 0)
            sex_exposed$condct_hpv11_anal[sex_exposed$HPV11_ANAL==1] = ifelse(sex_exposed$infection_hpv11_anal[sex_exposed$HPV11_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV11_ANAL==1], 1, 0)
            sex_exposed$partcondct_hpv11_anal[sex_exposed$HPV11_ANAL==1] = ifelse(sex_exposed$partinfection_hpv11_anal[sex_exposed$HPV11_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV11_ANAL==1], 1, 0)
            sex_exposed$condct_hpv16_anal[sex_exposed$HPV16_ANAL==1] = ifelse(sex_exposed$infection_hpv16_anal[sex_exposed$HPV16_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV16_ANAL==1], 1, 0)
            sex_exposed$partcondct_hpv16_anal[sex_exposed$HPV16_ANAL==1] = ifelse(sex_exposed$partinfection_hpv16_anal[sex_exposed$HPV16_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV16_ANAL==1], 1, 0)
            sex_exposed$condct_hpv18_anal[sex_exposed$HPV18_ANAL==1] = ifelse(sex_exposed$infection_hpv18_anal[sex_exposed$HPV18_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV18_ANAL==1], 1, 0)
            sex_exposed$partcondct_hpv18_anal[sex_exposed$HPV18_ANAL==1] = ifelse(sex_exposed$partinfection_hpv18_anal[sex_exposed$HPV18_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV18_ANAL==1], 1, 0)
            sex_exposed$condct_hpv31_anal[sex_exposed$HPV31_ANAL==1] = ifelse(sex_exposed$infection_hpv31_anal[sex_exposed$HPV31_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV31_ANAL==1], 1, 0)
            sex_exposed$partcondct_hpv31_anal[sex_exposed$HPV31_ANAL==1] = ifelse(sex_exposed$partinfection_hpv31_anal[sex_exposed$HPV31_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV31_ANAL==1], 1, 0)
            sex_exposed$condct_hpv33_anal[sex_exposed$HPV33_ANAL==1] = ifelse(sex_exposed$infection_hpv33_anal[sex_exposed$HPV33_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV33_ANAL==1], 1, 0)
            sex_exposed$partcondct_hpv33_anal[sex_exposed$HPV33_ANAL==1] = ifelse(sex_exposed$partinfection_hpv33_anal[sex_exposed$HPV33_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV33_ANAL==1], 1, 0)
            sex_exposed$condct_hpv45_anal[sex_exposed$HPV45_ANAL==1] = ifelse(sex_exposed$infection_hpv45_anal[sex_exposed$HPV45_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV45_ANAL==1], 1, 0)
            sex_exposed$partcondct_hpv45_anal[sex_exposed$HPV45_ANAL==1] = ifelse(sex_exposed$partinfection_hpv45_anal[sex_exposed$HPV45_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV45_ANAL==1], 1, 0)
            sex_exposed$condct_hpv52_anal[sex_exposed$HPV52_ANAL==1] = ifelse(sex_exposed$infection_hpv52_anal[sex_exposed$HPV52_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV52_ANAL==1], 1, 0)
            sex_exposed$partcondct_hpv52_anal[sex_exposed$HPV52_ANAL==1] = ifelse(sex_exposed$partinfection_hpv52_anal[sex_exposed$HPV52_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV52_ANAL==1], 1, 0)
            sex_exposed$condct_hpv58_anal[sex_exposed$HPV58_ANAL==1] = ifelse(sex_exposed$infection_hpv58_anal[sex_exposed$HPV58_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV58_ANAL==1], 1, 0)
            sex_exposed$partcondct_hpv58_anal[sex_exposed$HPV58_ANAL==1] = ifelse(sex_exposed$partinfection_hpv58_anal[sex_exposed$HPV58_ANAL==1]==1 & condom_success_anal[sex_exposed$HPV58_ANAL==1], 1, 0)
            
            #check if infection blocked for HPV oral
            sex_exposed$condct_hpv6_oral[sex_exposed$HPV6_ORAL==1] = ifelse(sex_exposed$infection_hpv6_oral[sex_exposed$HPV6_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV6_ORAL==1], 1, 0)
            sex_exposed$partcondct_hpv6_oral[sex_exposed$HPV6_ORAL==1] = ifelse(sex_exposed$partinfection_hpv6_oral[sex_exposed$HPV6_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV6_ORAL==1], 1, 0)
            sex_exposed$condct_hpv11_oral[sex_exposed$HPV11_ORAL==1] = ifelse(sex_exposed$infection_hpv11_oral[sex_exposed$HPV11_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV11_ORAL==1], 1, 0)
            sex_exposed$partcondct_hpv11_oral[sex_exposed$HPV11_ORAL==1] = ifelse(sex_exposed$partinfection_hpv11_oral[sex_exposed$HPV11_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV11_ORAL==1], 1, 0)
            sex_exposed$condct_hpv16_oral[sex_exposed$HPV16_ORAL==1] = ifelse(sex_exposed$infection_hpv16_oral[sex_exposed$HPV16_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV16_ORAL==1], 1, 0)
            sex_exposed$partcondct_hpv16_oral[sex_exposed$HPV16_ORAL==1] = ifelse(sex_exposed$partinfection_hpv16_oral[sex_exposed$HPV16_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV16_ORAL==1], 1, 0)
            sex_exposed$condct_hpv18_oral[sex_exposed$HPV18_ORAL==1] = ifelse(sex_exposed$infection_hpv18_oral[sex_exposed$HPV18_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV18_ORAL==1], 1, 0)
            sex_exposed$partcondct_hpv18_oral[sex_exposed$HPV18_ORAL==1] = ifelse(sex_exposed$partinfection_hpv18_oral[sex_exposed$HPV18_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV18_ORAL==1], 1, 0)
            sex_exposed$condct_hpv31_oral[sex_exposed$HPV31_ORAL==1] = ifelse(sex_exposed$infection_hpv31_oral[sex_exposed$HPV31_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV31_ORAL==1], 1, 0)
            sex_exposed$partcondct_hpv31_oral[sex_exposed$HPV31_ORAL==1] = ifelse(sex_exposed$partinfection_hpv31_oral[sex_exposed$HPV31_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV31_ORAL==1], 1, 0)
            sex_exposed$condct_hpv33_oral[sex_exposed$HPV33_ORAL==1] = ifelse(sex_exposed$infection_hpv33_oral[sex_exposed$HPV33_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV33_ORAL==1], 1, 0)
            sex_exposed$partcondct_hpv33_oral[sex_exposed$HPV33_ORAL==1] = ifelse(sex_exposed$partinfection_hpv33_oral[sex_exposed$HPV33_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV33_ORAL==1], 1, 0)
            sex_exposed$condct_hpv45_oral[sex_exposed$HPV45_ORAL==1] = ifelse(sex_exposed$infection_hpv45_oral[sex_exposed$HPV45_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV45_ORAL==1], 1, 0)
            sex_exposed$partcondct_hpv45_oral[sex_exposed$HPV45_ORAL==1] = ifelse(sex_exposed$partinfection_hpv45_oral[sex_exposed$HPV45_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV45_ORAL==1], 1, 0)
            sex_exposed$condct_hpv52_oral[sex_exposed$HPV52_ORAL==1] = ifelse(sex_exposed$infection_hpv52_oral[sex_exposed$HPV52_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV52_ORAL==1], 1, 0)
            sex_exposed$partcondct_hpv52_oral[sex_exposed$HPV52_ORAL==1] = ifelse(sex_exposed$partinfection_hpv52_oral[sex_exposed$HPV52_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV52_ORAL==1], 1, 0)
            sex_exposed$condct_hpv58_oral[sex_exposed$HPV58_ORAL==1] = ifelse(sex_exposed$infection_hpv58_oral[sex_exposed$HPV58_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV58_ORAL==1], 1, 0)
            sex_exposed$partcondct_hpv58_oral[sex_exposed$HPV58_ORAL==1] = ifelse(sex_exposed$partinfection_hpv58_oral[sex_exposed$HPV58_ORAL==1]>=1 & condom_success_oral[sex_exposed$HPV58_ORAL==1], 1, 0)
            
            rm(randcond,randfail,condom_worn_anal,condom_success_anal,condom_worn_oral,condom_success_oral)
          }
          
          #determine infection under TAP scenarios
          if (unique(abm_sims[[current_data]]$TAP)==1)
          {
            #check if infection blocked
            sex_exposed$tapct[sex_exposed$HIV==1] = ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$TAPVL[exposed_hiv_partner]==1, 1, 0)
            sex_exposed$parttapct[sex_exposed$HIV==1] = ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$TAPVL[exposed_hiv_ego]==1, 1, 0)
          }
          
          #determine infection under PrEP scenarios
          if (unique(abm_sims[[current_data]]$PREP)==1)
          {
            #check if infection blocked
            sex_exposed$prepct[sex_exposed$HIV==1] = ifelse(sex_exposed$infection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$ONPREP[exposed_hiv_ego]==1, 1, 0)
            sex_exposed$partprepct[sex_exposed$HIV==1] = ifelse(sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1 & abm_sims[[current_data]]$ONPREP[exposed_hiv_partner]==1, 1, 0)
          }
          
          #determine infection under vaccination scenarios
          if (unique(abm_sims[[current_data]]$VAX)==1)
          {
            #check if infection blocked for HPV anal
            #sex_exposed$vaxct_hpv_anal[sex_exposed$HPV_ANAL==1] = ifelse(sex_exposed$infection_hpv_anal[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv_anal_ego]==1, 1, 0)
            #sex_exposed$partvaxct_hpv_anal[sex_exposed$HPV_ANAL==1] = ifelse(sex_exposed$partinfection_hpv_anal[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv_anal_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv6_anal[sex_exposed$HPV6_ANAL==1] = ifelse(sex_exposed$infection_hpv6_anal[sex_exposed$HPV6_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv6_anal_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv6_anal[sex_exposed$HPV6_ANAL==1] = ifelse(sex_exposed$partinfection_hpv6_anal[sex_exposed$HPV6_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv6_anal_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv11_anal[sex_exposed$HPV11_ANAL==1] = ifelse(sex_exposed$infection_hpv11_anal[sex_exposed$HPV11_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv11_anal_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv11_anal[sex_exposed$HPV11_ANAL==1] = ifelse(sex_exposed$partinfection_hpv11_anal[sex_exposed$HPV11_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv11_anal_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv16_anal[sex_exposed$HPV16_ANAL==1] = ifelse(sex_exposed$infection_hpv16_anal[sex_exposed$HPV16_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv16_anal_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv16_anal[sex_exposed$HPV16_ANAL==1] = ifelse(sex_exposed$partinfection_hpv16_anal[sex_exposed$HPV16_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv16_anal_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv18_anal[sex_exposed$HPV18_ANAL==1] = ifelse(sex_exposed$infection_hpv18_anal[sex_exposed$HPV18_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv18_anal_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv18_anal[sex_exposed$HPV18_ANAL==1] = ifelse(sex_exposed$partinfection_hpv18_anal[sex_exposed$HPV18_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv18_anal_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv31_anal[sex_exposed$HPV31_ANAL==1] = ifelse(sex_exposed$infection_hpv31_anal[sex_exposed$HPV31_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv31_anal_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv31_anal[sex_exposed$HPV31_ANAL==1] = ifelse(sex_exposed$partinfection_hpv31_anal[sex_exposed$HPV31_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv31_anal_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv33_anal[sex_exposed$HPV33_ANAL==1] = ifelse(sex_exposed$infection_hpv33_anal[sex_exposed$HPV33_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv33_anal_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv33_anal[sex_exposed$HPV33_ANAL==1] = ifelse(sex_exposed$partinfection_hpv33_anal[sex_exposed$HPV33_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv33_anal_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv45_anal[sex_exposed$HPV45_ANAL==1] = ifelse(sex_exposed$infection_hpv45_anal[sex_exposed$HPV45_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv45_anal_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv45_anal[sex_exposed$HPV45_ANAL==1] = ifelse(sex_exposed$partinfection_hpv45_anal[sex_exposed$HPV45_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv45_anal_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv52_anal[sex_exposed$HPV52_ANAL==1] = ifelse(sex_exposed$infection_hpv52_anal[sex_exposed$HPV52_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv52_anal_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv52_anal[sex_exposed$HPV52_ANAL==1] = ifelse(sex_exposed$partinfection_hpv52_anal[sex_exposed$HPV52_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv52_anal_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv58_anal[sex_exposed$HPV58_ANAL==1] = ifelse(sex_exposed$infection_hpv58_anal[sex_exposed$HPV58_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv58_anal_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv58_anal[sex_exposed$HPV58_ANAL==1] = ifelse(sex_exposed$partinfection_hpv58_anal[sex_exposed$HPV58_ANAL==1]==1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv58_anal_partner]==1, 1, 0)
            
            #check if infection blocked for HPV oral
            sex_exposed$vaxct_hpv6_oral[sex_exposed$HPV6_ORAL==1] = ifelse(sex_exposed$infection_hpv6_oral[sex_exposed$HPV6_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv6_oral_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv6_oral[sex_exposed$HPV6_ORAL==1] = ifelse(sex_exposed$partinfection_hpv6_oral[sex_exposed$HPV6_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv6_oral_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv11_oral[sex_exposed$HPV11_ORAL==1] = ifelse(sex_exposed$infection_hpv11_oral[sex_exposed$HPV11_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv11_oral_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv11_oral[sex_exposed$HPV11_ORAL==1] = ifelse(sex_exposed$partinfection_hpv11_oral[sex_exposed$HPV11_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv11_oral_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv16_oral[sex_exposed$HPV16_ORAL==1] = ifelse(sex_exposed$infection_hpv16_oral[sex_exposed$HPV16_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv16_oral_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv16_oral[sex_exposed$HPV16_ORAL==1] = ifelse(sex_exposed$partinfection_hpv16_oral[sex_exposed$HPV16_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv16_oral_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv18_oral[sex_exposed$HPV18_ORAL==1] = ifelse(sex_exposed$infection_hpv18_oral[sex_exposed$HPV18_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv18_oral_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv18_oral[sex_exposed$HPV18_ORAL==1] = ifelse(sex_exposed$partinfection_hpv18_oral[sex_exposed$HPV18_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv18_oral_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv31_oral[sex_exposed$HPV31_ORAL==1] = ifelse(sex_exposed$infection_hpv31_oral[sex_exposed$HPV31_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv31_oral_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv31_oral[sex_exposed$HPV31_ORAL==1] = ifelse(sex_exposed$partinfection_hpv31_oral[sex_exposed$HPV31_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv31_oral_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv33_oral[sex_exposed$HPV33_ORAL==1] = ifelse(sex_exposed$infection_hpv33_oral[sex_exposed$HPV33_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv33_oral_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv33_oral[sex_exposed$HPV33_ORAL==1] = ifelse(sex_exposed$partinfection_hpv33_oral[sex_exposed$HPV33_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv33_oral_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv45_oral[sex_exposed$HPV45_ORAL==1] = ifelse(sex_exposed$infection_hpv45_oral[sex_exposed$HPV45_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv45_oral_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv45_oral[sex_exposed$HPV45_ORAL==1] = ifelse(sex_exposed$partinfection_hpv45_oral[sex_exposed$HPV45_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv45_oral_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv52_oral[sex_exposed$HPV52_ORAL==1] = ifelse(sex_exposed$infection_hpv52_oral[sex_exposed$HPV52_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv52_oral_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv52_oral[sex_exposed$HPV52_ORAL==1] = ifelse(sex_exposed$partinfection_hpv52_oral[sex_exposed$HPV52_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv52_oral_partner]==1, 1, 0)
            sex_exposed$vaxct_hpv58_oral[sex_exposed$HPV58_ORAL==1] = ifelse(sex_exposed$infection_hpv58_oral[sex_exposed$HPV58_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv58_oral_ego]==1, 1, 0)
            sex_exposed$partvaxct_hpv58_oral[sex_exposed$HPV58_ORAL==1] = ifelse(sex_exposed$partinfection_hpv58_oral[sex_exposed$HPV58_ORAL==1]>=1 & abm_sims[[current_data]]$IMMUNE[exposed_hpv58_oral_partner]==1, 1, 0)
          }
          
          #resolve infection statistics, HIV
          sex_exposed$newHIV[sex_exposed$HIV==1] = ifelse((sex_exposed$condct_hiv[sex_exposed$HIV==1] + sex_exposed$prepct[sex_exposed$HIV==1] + sex_exposed$tapct[sex_exposed$HIV==1] + sex_exposed$seroct_hiv[sex_exposed$HIV==1])==0 & sex_exposed$infection_hiv[sex_exposed$HIV==1]==1, 1, sex_exposed$newHIV[sex_exposed$HIV==1])
          sex_exposed$newVL[sex_exposed$HIV==1] = ifelse((sex_exposed$condct_hiv[sex_exposed$HIV==1] + sex_exposed$prepct[sex_exposed$HIV==1] + sex_exposed$tapct[sex_exposed$HIV==1] + sex_exposed$seroct_hiv[sex_exposed$HIV==1])==0 & sex_exposed$infection_hiv[sex_exposed$HIV==1]==1, 2, sex_exposed$newVL[sex_exposed$HIV==1])
          sex_exposed$partnewHIV[sex_exposed$HIV==1] = ifelse((sex_exposed$partcondct_hiv[sex_exposed$HIV==1] + sex_exposed$partprepct[sex_exposed$HIV==1] + sex_exposed$parttapct[sex_exposed$HIV==1] + sex_exposed$partseroct_hiv[sex_exposed$HIV==1])==0 & sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1, 1, sex_exposed$partnewHIV[sex_exposed$HIV==1])
          sex_exposed$partnewVL[sex_exposed$HIV==1] = ifelse((sex_exposed$partcondct_hiv[sex_exposed$HIV==1] + sex_exposed$partprepct[sex_exposed$HIV==1] + sex_exposed$parttapct[sex_exposed$HIV==1] + sex_exposed$partseroct_hiv[sex_exposed$HIV==1])==0 & sex_exposed$partinfection_hiv[sex_exposed$HIV==1]==1, 2, sex_exposed$partnewVL[sex_exposed$HIV==1])
          
          #resolve infection statistics, HPV anal
          sex_exposed$prevct_hpv6_anal[sex_exposed$HPV6_ANAL==1] = ifelse((sex_exposed$condct_hpv6_anal[sex_exposed$HPV6_ANAL==1] + sex_exposed$seroct_hpv6[sex_exposed$HPV6_ANAL==1] + sex_exposed$vaxct_hpv6_anal[sex_exposed$HPV6_ANAL==1])>0 & sex_exposed$infection_hpv6_anal[sex_exposed$HPV6_ANAL==1]==1, 1, sex_exposed$prevct_hpv6_anal[sex_exposed$HPV6_ANAL==1])
          sex_exposed$partprevct_hpv6_anal[sex_exposed$HPV6_ANAL==1] = ifelse((sex_exposed$partcondct_hpv6_anal[sex_exposed$HPV6_ANAL==1] + sex_exposed$partseroct_hpv6[sex_exposed$HPV6_ANAL==1] + sex_exposed$partvaxct_hpv6_anal[sex_exposed$HPV6_ANAL==1])>0 & sex_exposed$partinfection_hpv6_anal[sex_exposed$HPV6_ANAL==1]==1, 1, sex_exposed$partprevct_hpv6_anal[sex_exposed$HPV6_ANAL==1])
          sex_exposed$prevct_hpv11_anal[sex_exposed$HPV11_ANAL==1] = ifelse((sex_exposed$condct_hpv11_anal[sex_exposed$HPV11_ANAL==1] + sex_exposed$seroct_hpv11[sex_exposed$HPV11_ANAL==1] + sex_exposed$vaxct_hpv11_anal[sex_exposed$HPV11_ANAL==1])>0 & sex_exposed$infection_hpv11_anal[sex_exposed$HPV11_ANAL==1]==1, 1, sex_exposed$prevct_hpv11_anal[sex_exposed$HPV11_ANAL==1])
          sex_exposed$partprevct_hpv11_anal[sex_exposed$HPV11_ANAL==1] = ifelse((sex_exposed$partcondct_hpv11_anal[sex_exposed$HPV11_ANAL==1] + sex_exposed$partseroct_hpv11[sex_exposed$HPV11_ANAL==1] + sex_exposed$partvaxct_hpv11_anal[sex_exposed$HPV11_ANAL==1])>0 & sex_exposed$partinfection_hpv11_anal[sex_exposed$HPV11_ANAL==1]==1, 1, sex_exposed$partprevct_hpv11_anal[sex_exposed$HPV11_ANAL==1])
          sex_exposed$prevct_hpv16_anal[sex_exposed$HPV16_ANAL==1] = ifelse((sex_exposed$condct_hpv16_anal[sex_exposed$HPV16_ANAL==1] + sex_exposed$seroct_hpv16[sex_exposed$HPV16_ANAL==1] + sex_exposed$vaxct_hpv16_anal[sex_exposed$HPV16_ANAL==1])>0 & sex_exposed$infection_hpv16_anal[sex_exposed$HPV16_ANAL==1]==1, 1, sex_exposed$prevct_hpv16_anal[sex_exposed$HPV16_ANAL==1])
          sex_exposed$partprevct_hpv16_anal[sex_exposed$HPV16_ANAL==1] = ifelse((sex_exposed$partcondct_hpv16_anal[sex_exposed$HPV16_ANAL==1] + sex_exposed$partseroct_hpv16[sex_exposed$HPV16_ANAL==1] + sex_exposed$partvaxct_hpv16_anal[sex_exposed$HPV16_ANAL==1])>0 & sex_exposed$partinfection_hpv16_anal[sex_exposed$HPV16_ANAL==1]==1, 1, sex_exposed$partprevct_hpv16_anal[sex_exposed$HPV16_ANAL==1])
          sex_exposed$prevct_hpv18_anal[sex_exposed$HPV18_ANAL==1] = ifelse((sex_exposed$condct_hpv18_anal[sex_exposed$HPV18_ANAL==1] + sex_exposed$seroct_hpv18[sex_exposed$HPV18_ANAL==1] + sex_exposed$vaxct_hpv18_anal[sex_exposed$HPV18_ANAL==1])>0 & sex_exposed$infection_hpv18_anal[sex_exposed$HPV18_ANAL==1]==1, 1, sex_exposed$prevct_hpv18_anal[sex_exposed$HPV18_ANAL==1])
          sex_exposed$partprevct_hpv18_anal[sex_exposed$HPV18_ANAL==1] = ifelse((sex_exposed$partcondct_hpv18_anal[sex_exposed$HPV18_ANAL==1] + sex_exposed$partseroct_hpv18[sex_exposed$HPV18_ANAL==1] + sex_exposed$partvaxct_hpv18_anal[sex_exposed$HPV18_ANAL==1])>0 & sex_exposed$partinfection_hpv18_anal[sex_exposed$HPV18_ANAL==1]==1, 1, sex_exposed$partprevct_hpv18_anal[sex_exposed$HPV18_ANAL==1])
          sex_exposed$prevct_hpv31_anal[sex_exposed$HPV31_ANAL==1] = ifelse((sex_exposed$condct_hpv31_anal[sex_exposed$HPV31_ANAL==1] + sex_exposed$seroct_hpv31[sex_exposed$HPV31_ANAL==1] + sex_exposed$vaxct_hpv31_anal[sex_exposed$HPV31_ANAL==1])>0 & sex_exposed$infection_hpv31_anal[sex_exposed$HPV31_ANAL==1]==1, 1, sex_exposed$prevct_hpv31_anal[sex_exposed$HPV31_ANAL==1])
          sex_exposed$partprevct_hpv31_anal[sex_exposed$HPV31_ANAL==1] = ifelse((sex_exposed$partcondct_hpv31_anal[sex_exposed$HPV31_ANAL==1] + sex_exposed$partseroct_hpv31[sex_exposed$HPV31_ANAL==1] + sex_exposed$partvaxct_hpv31_anal[sex_exposed$HPV31_ANAL==1])>0 & sex_exposed$partinfection_hpv31_anal[sex_exposed$HPV31_ANAL==1]==1, 1, sex_exposed$partprevct_hpv31_anal[sex_exposed$HPV31_ANAL==1])
          sex_exposed$prevct_hpv33_anal[sex_exposed$HPV33_ANAL==1] = ifelse((sex_exposed$condct_hpv33_anal[sex_exposed$HPV33_ANAL==1] + sex_exposed$seroct_hpv33[sex_exposed$HPV33_ANAL==1] + sex_exposed$vaxct_hpv33_anal[sex_exposed$HPV33_ANAL==1])>0 & sex_exposed$infection_hpv33_anal[sex_exposed$HPV33_ANAL==1]==1, 1, sex_exposed$prevct_hpv33_anal[sex_exposed$HPV33_ANAL==1])
          sex_exposed$partprevct_hpv33_anal[sex_exposed$HPV33_ANAL==1] = ifelse((sex_exposed$partcondct_hpv33_anal[sex_exposed$HPV33_ANAL==1] + sex_exposed$partseroct_hpv33[sex_exposed$HPV33_ANAL==1] + sex_exposed$partvaxct_hpv33_anal[sex_exposed$HPV33_ANAL==1])>0 & sex_exposed$partinfection_hpv33_anal[sex_exposed$HPV33_ANAL==1]==1, 1, sex_exposed$partprevct_hpv33_anal[sex_exposed$HPV33_ANAL==1])
          sex_exposed$prevct_hpv45_anal[sex_exposed$HPV45_ANAL==1] = ifelse((sex_exposed$condct_hpv45_anal[sex_exposed$HPV45_ANAL==1] + sex_exposed$seroct_hpv45[sex_exposed$HPV45_ANAL==1] + sex_exposed$vaxct_hpv45_anal[sex_exposed$HPV45_ANAL==1])>0 & sex_exposed$infection_hpv45_anal[sex_exposed$HPV45_ANAL==1]==1, 1, sex_exposed$prevct_hpv45_anal[sex_exposed$HPV45_ANAL==1])
          sex_exposed$partprevct_hpv45_anal[sex_exposed$HPV45_ANAL==1] = ifelse((sex_exposed$partcondct_hpv45_anal[sex_exposed$HPV45_ANAL==1] + sex_exposed$partseroct_hpv45[sex_exposed$HPV45_ANAL==1] + sex_exposed$partvaxct_hpv45_anal[sex_exposed$HPV45_ANAL==1])>0 & sex_exposed$partinfection_hpv45_anal[sex_exposed$HPV45_ANAL==1]==1, 1, sex_exposed$partprevct_hpv45_anal[sex_exposed$HPV45_ANAL==1])
          sex_exposed$prevct_hpv52_anal[sex_exposed$HPV52_ANAL==1] = ifelse((sex_exposed$condct_hpv52_anal[sex_exposed$HPV52_ANAL==1] + sex_exposed$seroct_hpv52[sex_exposed$HPV52_ANAL==1] + sex_exposed$vaxct_hpv52_anal[sex_exposed$HPV52_ANAL==1])>0 & sex_exposed$infection_hpv52_anal[sex_exposed$HPV52_ANAL==1]==1, 1, sex_exposed$prevct_hpv52_anal[sex_exposed$HPV52_ANAL==1])
          sex_exposed$partprevct_hpv52_anal[sex_exposed$HPV52_ANAL==1] = ifelse((sex_exposed$partcondct_hpv52_anal[sex_exposed$HPV52_ANAL==1] + sex_exposed$partseroct_hpv52[sex_exposed$HPV52_ANAL==1] + sex_exposed$partvaxct_hpv52_anal[sex_exposed$HPV52_ANAL==1])>0 & sex_exposed$partinfection_hpv52_anal[sex_exposed$HPV52_ANAL==1]==1, 1, sex_exposed$partprevct_hpv52_anal[sex_exposed$HPV52_ANAL==1])
          sex_exposed$prevct_hpv58_anal[sex_exposed$HPV58_ANAL==1] = ifelse((sex_exposed$condct_hpv58_anal[sex_exposed$HPV58_ANAL==1] + sex_exposed$seroct_hpv58[sex_exposed$HPV58_ANAL==1] + sex_exposed$vaxct_hpv58_anal[sex_exposed$HPV58_ANAL==1])>0 & sex_exposed$infection_hpv58_anal[sex_exposed$HPV58_ANAL==1]==1, 1, sex_exposed$prevct_hpv58_anal[sex_exposed$HPV58_ANAL==1])
          sex_exposed$partprevct_hpv58_anal[sex_exposed$HPV58_ANAL==1] = ifelse((sex_exposed$partcondct_hpv58_anal[sex_exposed$HPV58_ANAL==1] + sex_exposed$partseroct_hpv58[sex_exposed$HPV58_ANAL==1] + sex_exposed$partvaxct_hpv58_anal[sex_exposed$HPV58_ANAL==1])>0 & sex_exposed$partinfection_hpv58_anal[sex_exposed$HPV58_ANAL==1]==1, 1, sex_exposed$partprevct_hpv58_anal[sex_exposed$HPV58_ANAL==1])
          sex_exposed$newHPV6_ANAL[sex_exposed$HPV6_ANAL==1] = ifelse((sex_exposed$condct_hpv6_anal[sex_exposed$HPV6_ANAL==1] + sex_exposed$seroct_hpv6[sex_exposed$HPV6_ANAL==1] + sex_exposed$vaxct_hpv6_anal[sex_exposed$HPV6_ANAL==1])==0 & sex_exposed$infection_hpv6_anal[sex_exposed$HPV6_ANAL==1]==1, 1, sex_exposed$newHPV6_ANAL[sex_exposed$HPV6_ANAL==1])
          sex_exposed$partnewHPV6_ANAL[sex_exposed$HPV6_ANAL==1] = ifelse((sex_exposed$partcondct_hpv6_anal[sex_exposed$HPV6_ANAL==1] + sex_exposed$partseroct_hpv6[sex_exposed$HPV6_ANAL==1] + sex_exposed$partvaxct_hpv6_anal[sex_exposed$HPV6_ANAL==1])==0 & sex_exposed$partinfection_hpv6_anal[sex_exposed$HPV6_ANAL==1]==1, 1, sex_exposed$partnewHPV6_ANAL[sex_exposed$HPV6_ANAL==1])
          sex_exposed$newHPV11_ANAL[sex_exposed$HPV11_ANAL==1] = ifelse((sex_exposed$condct_hpv11_anal[sex_exposed$HPV11_ANAL==1] + sex_exposed$seroct_hpv11[sex_exposed$HPV11_ANAL==1] + sex_exposed$vaxct_hpv11_anal[sex_exposed$HPV11_ANAL==1])==0 & sex_exposed$infection_hpv11_anal[sex_exposed$HPV11_ANAL==1]==1, 1, sex_exposed$newHPV11_ANAL[sex_exposed$HPV11_ANAL==1])
          sex_exposed$partnewHPV11_ANAL[sex_exposed$HPV11_ANAL==1] = ifelse((sex_exposed$partcondct_hpv11_anal[sex_exposed$HPV11_ANAL==1] + sex_exposed$partseroct_hpv11[sex_exposed$HPV11_ANAL==1] + sex_exposed$partvaxct_hpv11_anal[sex_exposed$HPV11_ANAL==1])==0 & sex_exposed$partinfection_hpv11_anal[sex_exposed$HPV11_ANAL==1]==1, 1, sex_exposed$partnewHPV11_ANAL[sex_exposed$HPV11_ANAL==1])
          sex_exposed$newHPV16_ANAL[sex_exposed$HPV16_ANAL==1] = ifelse((sex_exposed$condct_hpv16_anal[sex_exposed$HPV16_ANAL==1] + sex_exposed$seroct_hpv16[sex_exposed$HPV16_ANAL==1] + sex_exposed$vaxct_hpv16_anal[sex_exposed$HPV16_ANAL==1])==0 & sex_exposed$infection_hpv16_anal[sex_exposed$HPV16_ANAL==1]==1, 1, sex_exposed$newHPV16_ANAL[sex_exposed$HPV16_ANAL==1])
          sex_exposed$partnewHPV16_ANAL[sex_exposed$HPV16_ANAL==1] = ifelse((sex_exposed$partcondct_hpv16_anal[sex_exposed$HPV16_ANAL==1] + sex_exposed$partseroct_hpv16[sex_exposed$HPV16_ANAL==1] + sex_exposed$partvaxct_hpv16_anal[sex_exposed$HPV16_ANAL==1])==0 & sex_exposed$partinfection_hpv16_anal[sex_exposed$HPV16_ANAL==1]==1, 1, sex_exposed$partnewHPV16_ANAL[sex_exposed$HPV16_ANAL==1])
          sex_exposed$newHPV18_ANAL[sex_exposed$HPV18_ANAL==1] = ifelse((sex_exposed$condct_hpv18_anal[sex_exposed$HPV18_ANAL==1] + sex_exposed$seroct_hpv18[sex_exposed$HPV18_ANAL==1] + sex_exposed$vaxct_hpv18_anal[sex_exposed$HPV18_ANAL==1])==0 & sex_exposed$infection_hpv18_anal[sex_exposed$HPV18_ANAL==1]==1, 1, sex_exposed$newHPV18_ANAL[sex_exposed$HPV18_ANAL==1])
          sex_exposed$partnewHPV18_ANAL[sex_exposed$HPV18_ANAL==1] = ifelse((sex_exposed$partcondct_hpv18_anal[sex_exposed$HPV18_ANAL==1] + sex_exposed$partseroct_hpv18[sex_exposed$HPV18_ANAL==1] + sex_exposed$partvaxct_hpv18_anal[sex_exposed$HPV18_ANAL==1])==0 & sex_exposed$partinfection_hpv18_anal[sex_exposed$HPV18_ANAL==1]==1, 1, sex_exposed$partnewHPV18_ANAL[sex_exposed$HPV18_ANAL==1])
          sex_exposed$newHPV31_ANAL[sex_exposed$HPV31_ANAL==1] = ifelse((sex_exposed$condct_hpv31_anal[sex_exposed$HPV31_ANAL==1] + sex_exposed$seroct_hpv31[sex_exposed$HPV31_ANAL==1] + sex_exposed$vaxct_hpv31_anal[sex_exposed$HPV31_ANAL==1])==0 & sex_exposed$infection_hpv31_anal[sex_exposed$HPV31_ANAL==1]==1, 1, sex_exposed$newHPV31_ANAL[sex_exposed$HPV31_ANAL==1])
          sex_exposed$partnewHPV31_ANAL[sex_exposed$HPV31_ANAL==1] = ifelse((sex_exposed$partcondct_hpv31_anal[sex_exposed$HPV31_ANAL==1] + sex_exposed$partseroct_hpv31[sex_exposed$HPV31_ANAL==1] + sex_exposed$partvaxct_hpv31_anal[sex_exposed$HPV31_ANAL==1])==0 & sex_exposed$partinfection_hpv31_anal[sex_exposed$HPV31_ANAL==1]==1, 1, sex_exposed$partnewHPV31_ANAL[sex_exposed$HPV31_ANAL==1])
          sex_exposed$newHPV33_ANAL[sex_exposed$HPV33_ANAL==1] = ifelse((sex_exposed$condct_hpv33_anal[sex_exposed$HPV33_ANAL==1] + sex_exposed$seroct_hpv33[sex_exposed$HPV33_ANAL==1] + sex_exposed$vaxct_hpv33_anal[sex_exposed$HPV33_ANAL==1])==0 & sex_exposed$infection_hpv33_anal[sex_exposed$HPV33_ANAL==1]==1, 1, sex_exposed$newHPV33_ANAL[sex_exposed$HPV33_ANAL==1])
          sex_exposed$partnewHPV33_ANAL[sex_exposed$HPV33_ANAL==1] = ifelse((sex_exposed$partcondct_hpv33_anal[sex_exposed$HPV33_ANAL==1] + sex_exposed$partseroct_hpv33[sex_exposed$HPV33_ANAL==1] + sex_exposed$partvaxct_hpv33_anal[sex_exposed$HPV33_ANAL==1])==0 & sex_exposed$partinfection_hpv33_anal[sex_exposed$HPV33_ANAL==1]==1, 1, sex_exposed$partnewHPV33_ANAL[sex_exposed$HPV33_ANAL==1])
          sex_exposed$newHPV45_ANAL[sex_exposed$HPV45_ANAL==1] = ifelse((sex_exposed$condct_hpv45_anal[sex_exposed$HPV45_ANAL==1] + sex_exposed$seroct_hpv45[sex_exposed$HPV45_ANAL==1] + sex_exposed$vaxct_hpv45_anal[sex_exposed$HPV45_ANAL==1])==0 & sex_exposed$infection_hpv45_anal[sex_exposed$HPV45_ANAL==1]==1, 1, sex_exposed$newHPV45_ANAL[sex_exposed$HPV45_ANAL==1])
          sex_exposed$partnewHPV45_ANAL[sex_exposed$HPV45_ANAL==1] = ifelse((sex_exposed$partcondct_hpv45_anal[sex_exposed$HPV45_ANAL==1] + sex_exposed$partseroct_hpv45[sex_exposed$HPV45_ANAL==1] + sex_exposed$partvaxct_hpv45_anal[sex_exposed$HPV45_ANAL==1])==0 & sex_exposed$partinfection_hpv45_anal[sex_exposed$HPV45_ANAL==1]==1, 1, sex_exposed$partnewHPV45_ANAL[sex_exposed$HPV45_ANAL==1])
          sex_exposed$newHPV52_ANAL[sex_exposed$HPV52_ANAL==1] = ifelse((sex_exposed$condct_hpv52_anal[sex_exposed$HPV52_ANAL==1] + sex_exposed$seroct_hpv52[sex_exposed$HPV52_ANAL==1] + sex_exposed$vaxct_hpv52_anal[sex_exposed$HPV52_ANAL==1])==0 & sex_exposed$infection_hpv52_anal[sex_exposed$HPV52_ANAL==1]==1, 1, sex_exposed$newHPV52_ANAL[sex_exposed$HPV52_ANAL==1])
          sex_exposed$partnewHPV52_ANAL[sex_exposed$HPV52_ANAL==1] = ifelse((sex_exposed$partcondct_hpv52_anal[sex_exposed$HPV52_ANAL==1] + sex_exposed$partseroct_hpv52[sex_exposed$HPV52_ANAL==1] + sex_exposed$partvaxct_hpv52_anal[sex_exposed$HPV52_ANAL==1])==0 & sex_exposed$partinfection_hpv52_anal[sex_exposed$HPV52_ANAL==1]==1, 1, sex_exposed$partnewHPV52_ANAL[sex_exposed$HPV52_ANAL==1])
          sex_exposed$newHPV58_ANAL[sex_exposed$HPV58_ANAL==1] = ifelse((sex_exposed$condct_hpv58_anal[sex_exposed$HPV58_ANAL==1] + sex_exposed$seroct_hpv58[sex_exposed$HPV58_ANAL==1] + sex_exposed$vaxct_hpv58_anal[sex_exposed$HPV58_ANAL==1])==0 & sex_exposed$infection_hpv58_anal[sex_exposed$HPV58_ANAL==1]==1, 1, sex_exposed$newHPV58_ANAL[sex_exposed$HPV58_ANAL==1])
          sex_exposed$partnewHPV58_ANAL[sex_exposed$HPV58_ANAL==1] = ifelse((sex_exposed$partcondct_hpv58_anal[sex_exposed$HPV58_ANAL==1] + sex_exposed$partseroct_hpv58[sex_exposed$HPV58_ANAL==1] + sex_exposed$partvaxct_hpv58_anal[sex_exposed$HPV58_ANAL==1])==0 & sex_exposed$partinfection_hpv58_anal[sex_exposed$HPV58_ANAL==1]==1, 1, sex_exposed$partnewHPV58_ANAL[sex_exposed$HPV58_ANAL==1])
          #sex_exposed$prevct[sex_exposed$HPV_ANAL==1] = ifelse((sex_exposed$condct_hpv_anal[sex_exposed$HPV_ANAL==1] + sex_exposed$seroct_hpv[sex_exposed$HPV_ANAL==1] + sex_exposed$vaxct_hpv_anal[sex_exposed$HPV_ANAL==1])>0 & sex_exposed$infection_hpv_anal[sex_exposed$HPV_ANAL==1]==1, 1, sex_exposed$prevct[sex_exposed$HPV_ANAL==1])
          #sex_exposed$partprevct[sex_exposed$HPV_ANAL==1] = ifelse((sex_exposed$partcondct_hpv_anal[sex_exposed$HPV_ANAL==1] + sex_exposed$partseroct_hpv[sex_exposed$HPV_ANAL==1] + sex_exposed$partvaxct_hpv_anal[sex_exposed$HPV_ANAL==1])>0 & sex_exposed$partinfection_hpv_anal[sex_exposed$HPV_ANAL==1]==1, 1, sex_exposed$partprevct[sex_exposed$HPV_ANAL==1])
          #sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1] = ifelse((sex_exposed$condct_hpv_anal[sex_exposed$HPV_ANAL==1] + sex_exposed$seroct_hpv[sex_exposed$HPV_ANAL==1] + sex_exposed$vaxct_hpv_anal[sex_exposed$HPV_ANAL==1])==0 & sex_exposed$infection_hpv_anal[sex_exposed$HPV_ANAL==1]==1, 1, sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1])
          #sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1] = ifelse((sex_exposed$partcondct_hpv_anal[sex_exposed$HPV_ANAL==1] + sex_exposed$partseroct_hpv[sex_exposed$HPV_ANAL==1] + sex_exposed$partvaxct_hpv_anal[sex_exposed$HPV_ANAL==1])==0 & sex_exposed$partinfection_hpv_anal[sex_exposed$HPV_ANAL==1]==1, 1, sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1])
          #sex_exposed$partinfectct[sex_exposed$HPV_ANAL==1] = ifelse((sex_exposed$condct_hpv_anal[sex_exposed$HPV_ANAL==1] + sex_exposed$seroct_hpv[sex_exposed$HPV_ANAL==1] + sex_exposed$vaxct_hpv_anal[sex_exposed$HPV_ANAL==1])==0 & sex_exposed$infection_hpv_anal[sex_exposed$HPV_ANAL==1]==1, 1, sex_exposed$partinfectct[sex_exposed$HPV_ANAL==1])
          #sex_exposed$infectct[sex_exposed$HPV_ANAL==1] = ifelse((sex_exposed$partcondct_hpv_anal[sex_exposed$HPV_ANAL==1] + sex_exposed$partseroct_hpv[sex_exposed$HPV_ANAL==1] + sex_exposed$partvaxct_hpv_anal[sex_exposed$HPV_ANAL==1])==0 & sex_exposed$partinfection_hpv_anal[sex_exposed$HPV_ANAL==1]==1, 1, sex_exposed$infectct[sex_exposed$HPV_ANAL==1])
          
          #resolve infection statistics, HPV oral
          sex_exposed$prevct_hpv6_oral[sex_exposed$HPV6_ORAL==1] = ifelse((sex_exposed$condct_hpv6_oral[sex_exposed$HPV6_ORAL==1] + sex_exposed$vaxct_hpv6_oral[sex_exposed$HPV6_ORAL==1])>0 & sex_exposed$infection_hpv6_oral[sex_exposed$HPV6_ORAL==1]>=1, 1, sex_exposed$prevct_hpv6_oral[sex_exposed$HPV6_ORAL==1])
          sex_exposed$partprevct_hpv6_oral[sex_exposed$HPV6_ORAL==1] = ifelse((sex_exposed$partcondct_hpv6_oral[sex_exposed$HPV6_ORAL==1] + sex_exposed$partvaxct_hpv6_oral[sex_exposed$HPV6_ORAL==1])>0 & sex_exposed$partinfection_hpv6_oral[sex_exposed$HPV6_ORAL==1]>=1, 1, sex_exposed$partprevct_hpv6_oral[sex_exposed$HPV6_ORAL==1])
          sex_exposed$prevct_hpv11_oral[sex_exposed$HPV11_ORAL==1] = ifelse((sex_exposed$condct_hpv11_oral[sex_exposed$HPV11_ORAL==1] + sex_exposed$vaxct_hpv11_oral[sex_exposed$HPV11_ORAL==1])>0 & sex_exposed$infection_hpv11_oral[sex_exposed$HPV11_ORAL==1]>=1, 1, sex_exposed$prevct_hpv11_oral[sex_exposed$HPV11_ORAL==1])
          sex_exposed$partprevct_hpv11_oral[sex_exposed$HPV11_ORAL==1] = ifelse((sex_exposed$partcondct_hpv11_oral[sex_exposed$HPV11_ORAL==1] + sex_exposed$partvaxct_hpv11_oral[sex_exposed$HPV11_ORAL==1])>0 & sex_exposed$partinfection_hpv11_oral[sex_exposed$HPV11_ORAL==1]>=1, 1, sex_exposed$partprevct_hpv11_oral[sex_exposed$HPV11_ORAL==1])
          sex_exposed$prevct_hpv16_oral[sex_exposed$HPV16_ORAL==1] = ifelse((sex_exposed$condct_hpv16_oral[sex_exposed$HPV16_ORAL==1] + sex_exposed$vaxct_hpv16_oral[sex_exposed$HPV16_ORAL==1])>0 & sex_exposed$infection_hpv16_oral[sex_exposed$HPV16_ORAL==1]>=1, 1, sex_exposed$prevct_hpv16_oral[sex_exposed$HPV16_ORAL==1])
          sex_exposed$partprevct_hpv16_oral[sex_exposed$HPV16_ORAL==1] = ifelse((sex_exposed$partcondct_hpv16_oral[sex_exposed$HPV16_ORAL==1] + sex_exposed$partvaxct_hpv16_oral[sex_exposed$HPV16_ORAL==1])>0 & sex_exposed$partinfection_hpv16_oral[sex_exposed$HPV16_ORAL==1]>=1, 1, sex_exposed$partprevct_hpv16_oral[sex_exposed$HPV16_ORAL==1])
          sex_exposed$prevct_hpv18_oral[sex_exposed$HPV18_ORAL==1] = ifelse((sex_exposed$condct_hpv18_oral[sex_exposed$HPV18_ORAL==1] + sex_exposed$vaxct_hpv18_oral[sex_exposed$HPV18_ORAL==1])>0 & sex_exposed$infection_hpv18_oral[sex_exposed$HPV18_ORAL==1]>=1, 1, sex_exposed$prevct_hpv18_oral[sex_exposed$HPV18_ORAL==1])
          sex_exposed$partprevct_hpv18_oral[sex_exposed$HPV18_ORAL==1] = ifelse((sex_exposed$partcondct_hpv18_oral[sex_exposed$HPV18_ORAL==1] + sex_exposed$partvaxct_hpv18_oral[sex_exposed$HPV18_ORAL==1])>0 & sex_exposed$partinfection_hpv18_oral[sex_exposed$HPV18_ORAL==1]>=1, 1, sex_exposed$partprevct_hpv18_oral[sex_exposed$HPV18_ORAL==1])
          sex_exposed$prevct_hpv31_oral[sex_exposed$HPV31_ORAL==1] = ifelse((sex_exposed$condct_hpv31_oral[sex_exposed$HPV31_ORAL==1] + sex_exposed$vaxct_hpv31_oral[sex_exposed$HPV31_ORAL==1])>0 & sex_exposed$infection_hpv31_oral[sex_exposed$HPV31_ORAL==1]>=1, 1, sex_exposed$prevct_hpv31_oral[sex_exposed$HPV31_ORAL==1])
          sex_exposed$partprevct_hpv31_oral[sex_exposed$HPV31_ORAL==1] = ifelse((sex_exposed$partcondct_hpv31_oral[sex_exposed$HPV31_ORAL==1] + sex_exposed$partvaxct_hpv31_oral[sex_exposed$HPV31_ORAL==1])>0 & sex_exposed$partinfection_hpv31_oral[sex_exposed$HPV31_ORAL==1]>=1, 1, sex_exposed$partprevct_hpv31_oral[sex_exposed$HPV31_ORAL==1])
          sex_exposed$prevct_hpv33_oral[sex_exposed$HPV33_ORAL==1] = ifelse((sex_exposed$condct_hpv33_oral[sex_exposed$HPV33_ORAL==1] + sex_exposed$vaxct_hpv33_oral[sex_exposed$HPV33_ORAL==1])>0 & sex_exposed$infection_hpv33_oral[sex_exposed$HPV33_ORAL==1]>=1, 1, sex_exposed$prevct_hpv33_oral[sex_exposed$HPV33_ORAL==1])
          sex_exposed$partprevct_hpv33_oral[sex_exposed$HPV33_ORAL==1] = ifelse((sex_exposed$partcondct_hpv33_oral[sex_exposed$HPV33_ORAL==1] + sex_exposed$partvaxct_hpv33_oral[sex_exposed$HPV33_ORAL==1])>0 & sex_exposed$partinfection_hpv33_oral[sex_exposed$HPV33_ORAL==1]>=1, 1, sex_exposed$partprevct_hpv33_oral[sex_exposed$HPV33_ORAL==1])
          sex_exposed$prevct_hpv45_oral[sex_exposed$HPV45_ORAL==1] = ifelse((sex_exposed$condct_hpv45_oral[sex_exposed$HPV45_ORAL==1] + sex_exposed$vaxct_hpv45_oral[sex_exposed$HPV45_ORAL==1])>0 & sex_exposed$infection_hpv45_oral[sex_exposed$HPV45_ORAL==1]>=1, 1, sex_exposed$prevct_hpv45_oral[sex_exposed$HPV45_ORAL==1])
          sex_exposed$partprevct_hpv45_oral[sex_exposed$HPV45_ORAL==1] = ifelse((sex_exposed$partcondct_hpv45_oral[sex_exposed$HPV45_ORAL==1] + sex_exposed$partvaxct_hpv45_oral[sex_exposed$HPV45_ORAL==1])>0 & sex_exposed$partinfection_hpv45_oral[sex_exposed$HPV45_ORAL==1]>=1, 1, sex_exposed$partprevct_hpv45_oral[sex_exposed$HPV45_ORAL==1])
          sex_exposed$prevct_hpv52_oral[sex_exposed$HPV52_ORAL==1] = ifelse((sex_exposed$condct_hpv52_oral[sex_exposed$HPV52_ORAL==1] + sex_exposed$vaxct_hpv52_oral[sex_exposed$HPV52_ORAL==1])>0 & sex_exposed$infection_hpv52_oral[sex_exposed$HPV52_ORAL==1]>=1, 1, sex_exposed$prevct_hpv52_oral[sex_exposed$HPV52_ORAL==1])
          sex_exposed$partprevct_hpv52_oral[sex_exposed$HPV52_ORAL==1] = ifelse((sex_exposed$partcondct_hpv52_oral[sex_exposed$HPV52_ORAL==1] + sex_exposed$partvaxct_hpv52_oral[sex_exposed$HPV52_ORAL==1])>0 & sex_exposed$partinfection_hpv52_oral[sex_exposed$HPV52_ORAL==1]>=1, 1, sex_exposed$partprevct_hpv52_oral[sex_exposed$HPV52_ORAL==1])
          sex_exposed$prevct_hpv58_oral[sex_exposed$HPV58_ORAL==1] = ifelse((sex_exposed$condct_hpv58_oral[sex_exposed$HPV58_ORAL==1] + sex_exposed$vaxct_hpv58_oral[sex_exposed$HPV58_ORAL==1])>0 & sex_exposed$infection_hpv58_oral[sex_exposed$HPV58_ORAL==1]>=1, 1, sex_exposed$prevct_hpv58_oral[sex_exposed$HPV58_ORAL==1])
          sex_exposed$partprevct_hpv58_oral[sex_exposed$HPV58_ORAL==1] = ifelse((sex_exposed$partcondct_hpv58_oral[sex_exposed$HPV58_ORAL==1] + sex_exposed$partvaxct_hpv58_oral[sex_exposed$HPV58_ORAL==1])>0 & sex_exposed$partinfection_hpv58_oral[sex_exposed$HPV58_ORAL==1]>=1, 1, sex_exposed$partprevct_hpv58_oral[sex_exposed$HPV58_ORAL==1])
          sex_exposed$newHPV6_ORAL[sex_exposed$HPV6_ORAL==1] = ifelse((sex_exposed$condct_hpv6_oral[sex_exposed$HPV6_ORAL==1] + sex_exposed$vaxct_hpv6_oral[sex_exposed$HPV6_ORAL==1])==0 & sex_exposed$infection_hpv6_oral[sex_exposed$HPV6_ORAL==1]>=1, sex_exposed$infection_hpv6_oral[sex_exposed$HPV6_ORAL==1], sex_exposed$newHPV6_ORAL[sex_exposed$HPV6_ORAL==1])
          sex_exposed$partnewHPV6_ORAL[sex_exposed$HPV6_ORAL==1] = ifelse((sex_exposed$partcondct_hpv6_oral[sex_exposed$HPV6_ORAL==1] + sex_exposed$partvaxct_hpv6_oral[sex_exposed$HPV6_ORAL==1])==0 & sex_exposed$partinfection_hpv6_oral[sex_exposed$HPV6_ORAL==1]>=1, sex_exposed$partinfection_hpv6_oral[sex_exposed$HPV6_ORAL==1], sex_exposed$partnewHPV6_ORAL[sex_exposed$HPV6_ORAL==1])
          sex_exposed$newHPV11_ORAL[sex_exposed$HPV11_ORAL==1] = ifelse((sex_exposed$condct_hpv11_oral[sex_exposed$HPV11_ORAL==1] + sex_exposed$vaxct_hpv11_oral[sex_exposed$HPV11_ORAL==1])==0 & sex_exposed$infection_hpv11_oral[sex_exposed$HPV11_ORAL==1]>=1, sex_exposed$infection_hpv11_oral[sex_exposed$HPV11_ORAL==1], sex_exposed$newHPV11_ORAL[sex_exposed$HPV11_ORAL==1])
          sex_exposed$partnewHPV11_ORAL[sex_exposed$HPV11_ORAL==1] = ifelse((sex_exposed$partcondct_hpv11_oral[sex_exposed$HPV11_ORAL==1] + sex_exposed$partvaxct_hpv11_oral[sex_exposed$HPV11_ORAL==1])==0 & sex_exposed$partinfection_hpv11_oral[sex_exposed$HPV11_ORAL==1]>=1, sex_exposed$partinfection_hpv11_oral[sex_exposed$HPV11_ORAL==1], sex_exposed$partnewHPV11_ORAL[sex_exposed$HPV11_ORAL==1])
          sex_exposed$newHPV16_ORAL[sex_exposed$HPV16_ORAL==1] = ifelse((sex_exposed$condct_hpv16_oral[sex_exposed$HPV16_ORAL==1] + sex_exposed$vaxct_hpv16_oral[sex_exposed$HPV16_ORAL==1])==0 & sex_exposed$infection_hpv16_oral[sex_exposed$HPV16_ORAL==1]>=1, sex_exposed$infection_hpv16_oral[sex_exposed$HPV16_ORAL==1], sex_exposed$newHPV16_ORAL[sex_exposed$HPV16_ORAL==1])
          sex_exposed$partnewHPV16_ORAL[sex_exposed$HPV16_ORAL==1] = ifelse((sex_exposed$partcondct_hpv16_oral[sex_exposed$HPV16_ORAL==1] + sex_exposed$partvaxct_hpv16_oral[sex_exposed$HPV16_ORAL==1])==0 & sex_exposed$partinfection_hpv16_oral[sex_exposed$HPV16_ORAL==1]>=1, sex_exposed$partinfection_hpv16_oral[sex_exposed$HPV16_ORAL==1], sex_exposed$partnewHPV16_ORAL[sex_exposed$HPV16_ORAL==1])
          sex_exposed$newHPV18_ORAL[sex_exposed$HPV18_ORAL==1] = ifelse((sex_exposed$condct_hpv18_oral[sex_exposed$HPV18_ORAL==1] + sex_exposed$vaxct_hpv18_oral[sex_exposed$HPV18_ORAL==1])==0 & sex_exposed$infection_hpv18_oral[sex_exposed$HPV18_ORAL==1]>=1, sex_exposed$infection_hpv18_oral[sex_exposed$HPV18_ORAL==1], sex_exposed$newHPV18_ORAL[sex_exposed$HPV18_ORAL==1])
          sex_exposed$partnewHPV18_ORAL[sex_exposed$HPV18_ORAL==1] = ifelse((sex_exposed$partcondct_hpv18_oral[sex_exposed$HPV18_ORAL==1] + sex_exposed$partvaxct_hpv18_oral[sex_exposed$HPV18_ORAL==1])==0 & sex_exposed$partinfection_hpv18_oral[sex_exposed$HPV18_ORAL==1]>=1, sex_exposed$partinfection_hpv18_oral[sex_exposed$HPV18_ORAL==1], sex_exposed$partnewHPV18_ORAL[sex_exposed$HPV18_ORAL==1])
          sex_exposed$newHPV31_ORAL[sex_exposed$HPV31_ORAL==1] = ifelse((sex_exposed$condct_hpv31_oral[sex_exposed$HPV31_ORAL==1] + sex_exposed$vaxct_hpv31_oral[sex_exposed$HPV31_ORAL==1])==0 & sex_exposed$infection_hpv31_oral[sex_exposed$HPV31_ORAL==1]>=1, sex_exposed$infection_hpv31_oral[sex_exposed$HPV31_ORAL==1], sex_exposed$newHPV31_ORAL[sex_exposed$HPV31_ORAL==1])
          sex_exposed$partnewHPV31_ORAL[sex_exposed$HPV31_ORAL==1] = ifelse((sex_exposed$partcondct_hpv31_oral[sex_exposed$HPV31_ORAL==1] + sex_exposed$partvaxct_hpv31_oral[sex_exposed$HPV31_ORAL==1])==0 & sex_exposed$partinfection_hpv31_oral[sex_exposed$HPV31_ORAL==1]>=1, sex_exposed$partinfection_hpv31_oral[sex_exposed$HPV31_ORAL==1], sex_exposed$partnewHPV31_ORAL[sex_exposed$HPV31_ORAL==1])
          sex_exposed$newHPV33_ORAL[sex_exposed$HPV33_ORAL==1] = ifelse((sex_exposed$condct_hpv33_oral[sex_exposed$HPV33_ORAL==1] + sex_exposed$vaxct_hpv33_oral[sex_exposed$HPV33_ORAL==1])==0 & sex_exposed$infection_hpv33_oral[sex_exposed$HPV33_ORAL==1]>=1, sex_exposed$infection_hpv33_oral[sex_exposed$HPV33_ORAL==1], sex_exposed$newHPV33_ORAL[sex_exposed$HPV33_ORAL==1])
          sex_exposed$partnewHPV33_ORAL[sex_exposed$HPV33_ORAL==1] = ifelse((sex_exposed$partcondct_hpv33_oral[sex_exposed$HPV33_ORAL==1] + sex_exposed$partvaxct_hpv33_oral[sex_exposed$HPV33_ORAL==1])==0 & sex_exposed$partinfection_hpv33_oral[sex_exposed$HPV33_ORAL==1]>=1, sex_exposed$partinfection_hpv33_oral[sex_exposed$HPV33_ORAL==1], sex_exposed$partnewHPV33_ORAL[sex_exposed$HPV33_ORAL==1])
          sex_exposed$newHPV45_ORAL[sex_exposed$HPV45_ORAL==1] = ifelse((sex_exposed$condct_hpv45_oral[sex_exposed$HPV45_ORAL==1] + sex_exposed$vaxct_hpv45_oral[sex_exposed$HPV45_ORAL==1])==0 & sex_exposed$infection_hpv45_oral[sex_exposed$HPV45_ORAL==1]>=1, sex_exposed$infection_hpv45_oral[sex_exposed$HPV45_ORAL==1], sex_exposed$newHPV45_ORAL[sex_exposed$HPV45_ORAL==1])
          sex_exposed$partnewHPV45_ORAL[sex_exposed$HPV45_ORAL==1] = ifelse((sex_exposed$partcondct_hpv45_oral[sex_exposed$HPV45_ORAL==1] + sex_exposed$partvaxct_hpv45_oral[sex_exposed$HPV45_ORAL==1])==0 & sex_exposed$partinfection_hpv45_oral[sex_exposed$HPV45_ORAL==1]>=1, sex_exposed$partinfection_hpv45_oral[sex_exposed$HPV45_ORAL==1], sex_exposed$partnewHPV45_ORAL[sex_exposed$HPV45_ORAL==1])
          sex_exposed$newHPV52_ORAL[sex_exposed$HPV52_ORAL==1] = ifelse((sex_exposed$condct_hpv52_oral[sex_exposed$HPV52_ORAL==1] + sex_exposed$vaxct_hpv52_oral[sex_exposed$HPV52_ORAL==1])==0 & sex_exposed$infection_hpv52_oral[sex_exposed$HPV52_ORAL==1]>=1, sex_exposed$infection_hpv52_oral[sex_exposed$HPV52_ORAL==1], sex_exposed$newHPV52_ORAL[sex_exposed$HPV52_ORAL==1])
          sex_exposed$partnewHPV52_ORAL[sex_exposed$HPV52_ORAL==1] = ifelse((sex_exposed$partcondct_hpv52_oral[sex_exposed$HPV52_ORAL==1] + sex_exposed$partvaxct_hpv52_oral[sex_exposed$HPV52_ORAL==1])==0 & sex_exposed$partinfection_hpv52_oral[sex_exposed$HPV52_ORAL==1]>=1, sex_exposed$partinfection_hpv52_oral[sex_exposed$HPV52_ORAL==1], sex_exposed$partnewHPV52_ORAL[sex_exposed$HPV52_ORAL==1])
          sex_exposed$newHPV58_ORAL[sex_exposed$HPV58_ORAL==1] = ifelse((sex_exposed$condct_hpv58_oral[sex_exposed$HPV58_ORAL==1] + sex_exposed$vaxct_hpv58_oral[sex_exposed$HPV58_ORAL==1])==0 & sex_exposed$infection_hpv58_oral[sex_exposed$HPV58_ORAL==1]>=1, sex_exposed$infection_hpv58_oral[sex_exposed$HPV58_ORAL==1], sex_exposed$newHPV58_ORAL[sex_exposed$HPV58_ORAL==1])
          sex_exposed$partnewHPV58_ORAL[sex_exposed$HPV58_ORAL==1] = ifelse((sex_exposed$partcondct_hpv58_oral[sex_exposed$HPV58_ORAL==1] + sex_exposed$partvaxct_hpv58_oral[sex_exposed$HPV58_ORAL==1])==0 & sex_exposed$partinfection_hpv58_oral[sex_exposed$HPV58_ORAL==1]>=1, sex_exposed$partinfection_hpv58_oral[sex_exposed$HPV58_ORAL==1], sex_exposed$partnewHPV58_ORAL[sex_exposed$HPV58_ORAL==1])
          
          ## UPDATE STATS IN MAIN DATASET ##
          
          #HIV
          #abm_sims[[current_data]]$PREP_PREVENT[exposed_hiv_ego] = abm_sims[[current_data]]$PREP_PREVENT[exposed_hiv_ego] + sex_exposed$prepct[sex_exposed$HIV==1]
          #abm_sims[[current_data]]$PREP_PREVENT[exposed_hiv_partner] = abm_sims[[current_data]]$PREP_PREVENT[exposed_hiv_partner] + sex_exposed$partprepct[sex_exposed$HIV==1]
          #abm_sims[[current_data]]$TAP_PREVENT[exposed_hpv_anal_ego] = abm_sims[[current_data]]$TAP_PREVENT[exposed_hpv_anal_ego] + sex_exposed$tapct[sex_exposed$HIV==1]
          #abm_sims[[current_data]]$TAP_PREVENT[exposed_hpv_anal_partner] = abm_sims[[current_data]]$TAP_PREVENT[exposed_hpv_anal_partner] + sex_exposed$parttapct[sex_exposed$HIV==1]
          #abm_sims[[current_data]]$CAUSE_INFECT[exposed_hiv_ego] = ifelse(sex_exposed$infectct[sex_exposed$HIV==1]==1, abm_sims[[current_data]]$CAUSE_INFECT[exposed_hiv_ego] + sex_exposed$infectct[sex_exposed$HIV==1], abm_sims[[current_data]]$CAUSE_INFECT[exposed_hiv_ego])
          #abm_sims[[current_data]]$CAUSE_INFECT[exposed_hiv_partner] = ifelse(sex_exposed$partinfectct[sex_exposed$HIV==1]==1, abm_sims[[current_data]]$CAUSE_INFECT[exposed_hiv_partner] + sex_exposed$partinfectct[sex_exposed$HIV==1], abm_sims[[current_data]]$CAUSE_INFECT[exposed_hiv_partner])
          abm_sims[[current_data]]$HIV[exposed_hiv_ego] = ifelse(sex_exposed$newHIV[sex_exposed$HIV==1]==1, 1, abm_sims[[current_data]]$HIV[exposed_hiv_ego])
          abm_sims[[current_data]]$HIV[exposed_hiv_partner] = ifelse(sex_exposed$partnewHIV[sex_exposed$HIV==1]==1, 1, abm_sims[[current_data]]$HIV[exposed_hiv_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HIV[exposed_hiv_ego] = ifelse(sex_exposed$newHIV[sex_exposed$HIV==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HIV[exposed_hiv_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HIV[exposed_hiv_partner] = ifelse(sex_exposed$partnewHIV[sex_exposed$HIV==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HIV[exposed_hiv_partner])
          abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego] = ifelse(sex_exposed$newVL[sex_exposed$HIV==1]==2, 2, abm_sims[[current_data]]$HIV_VL[exposed_hiv_ego])
          abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner] = ifelse(sex_exposed$partnewVL[sex_exposed$HIV==1]==2, 2, abm_sims[[current_data]]$HIV_VL[exposed_hiv_partner])
          #abm_sims[[current_data]]$DISCORDANT[exposed_hiv_ego] = abm_sims[[current_data]]$DISCORDANT[exposed_hiv_ego] + 1
          #abm_sims[[current_data]]$DISCORDANT[exposed_hiv_partner] = abm_sims[[current_data]]$DISCORDANT[exposed_hiv_partner] + 1
          
          #HPV anal
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv6_anal_ego] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv6_anal_ego] + sex_exposed$condct_hpv6_anal[sex_exposed$HPV6_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv6_anal_partner] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv6_anal_partner] + sex_exposed$partcondct_hpv6_anal[sex_exposed$HPV6_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv11_anal_ego] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv11_anal_ego] + sex_exposed$condct_hpv11_anal[sex_exposed$HPV11_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv11_anal_partner] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv11_anal_partner] + sex_exposed$partcondct_hpv11_anal[sex_exposed$HPV11_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv16_anal_ego] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv16_anal_ego] + sex_exposed$condct_hpv16_anal[sex_exposed$HPV16_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv16_anal_partner] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv16_anal_partner] + sex_exposed$partcondct_hpv16_anal[sex_exposed$HPV16_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv18_anal_ego] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv18_anal_ego] + sex_exposed$condct_hpv18_anal[sex_exposed$HPV18_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv18_anal_partner] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv18_anal_partner] + sex_exposed$partcondct_hpv18_anal[sex_exposed$HPV18_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv31_anal_ego] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv31_anal_ego] + sex_exposed$condct_hpv31_anal[sex_exposed$HPV31_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv31_anal_partner] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv31_anal_partner] + sex_exposed$partcondct_hpv31_anal[sex_exposed$HPV31_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv33_anal_ego] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv33_anal_ego] + sex_exposed$condct_hpv33_anal[sex_exposed$HPV33_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv33_anal_partner] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv33_anal_partner] + sex_exposed$partcondct_hpv33_anal[sex_exposed$HPV33_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv45_anal_ego] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv45_anal_ego] + sex_exposed$condct_hpv45_anal[sex_exposed$HPV45_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv45_anal_partner] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv45_anal_partner] + sex_exposed$partcondct_hpv45_anal[sex_exposed$HPV45_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv52_anal_ego] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv52_anal_ego] + sex_exposed$condct_hpv52_anal[sex_exposed$HPV52_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv52_anal_partner] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv52_anal_partner] + sex_exposed$partcondct_hpv52_anal[sex_exposed$HPV52_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv58_anal_ego] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv58_anal_ego] + sex_exposed$condct_hpv58_anal[sex_exposed$HPV58_ANAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv58_anal_partner] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv58_anal_partner] + sex_exposed$partcondct_hpv58_anal[sex_exposed$HPV58_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv6_anal_ego] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv6_anal_ego] + sex_exposed$seroct_hpv6[sex_exposed$HPV6_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv6_anal_partner] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv6_anal_partner] + sex_exposed$partseroct_hpv6[sex_exposed$HPV6_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv11_anal_ego] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv11_anal_ego] + sex_exposed$seroct_hpv11[sex_exposed$HPV11_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv11_anal_partner] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv11_anal_partner] + sex_exposed$partseroct_hpv11[sex_exposed$HPV11_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv16_anal_ego] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv16_anal_ego] + sex_exposed$seroct_hpv16[sex_exposed$HPV16_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv16_anal_partner] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv16_anal_partner] + sex_exposed$partseroct_hpv16[sex_exposed$HPV16_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv18_anal_ego] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv18_anal_ego] + sex_exposed$seroct_hpv18[sex_exposed$HPV18_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv18_anal_partner] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv18_anal_partner] + sex_exposed$partseroct_hpv18[sex_exposed$HPV18_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv31_anal_ego] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv31_anal_ego] + sex_exposed$seroct_hpv31[sex_exposed$HPV31_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv31_anal_partner] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv31_anal_partner] + sex_exposed$partseroct_hpv31[sex_exposed$HPV31_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv33_anal_ego] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv33_anal_ego] + sex_exposed$seroct_hpv33[sex_exposed$HPV33_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv33_anal_partner] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv33_anal_partner] + sex_exposed$partseroct_hpv33[sex_exposed$HPV33_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv45_anal_ego] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv45_anal_ego] + sex_exposed$seroct_hpv45[sex_exposed$HPV45_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv45_anal_partner] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv45_anal_partner] + sex_exposed$partseroct_hpv45[sex_exposed$HPV45_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv52_anal_ego] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv52_anal_ego] + sex_exposed$seroct_hpv52[sex_exposed$HPV52_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv52_anal_partner] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv52_anal_partner] + sex_exposed$partseroct_hpv52[sex_exposed$HPV52_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv58_anal_ego] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv58_anal_ego] + sex_exposed$seroct_hpv58[sex_exposed$HPV58_ANAL==1]
          abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv58_anal_partner] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv58_anal_partner] + sex_exposed$partseroct_hpv58[sex_exposed$HPV58_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv6_anal_ego] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv6_anal_ego] + sex_exposed$vaxct_hpv6_anal[sex_exposed$HPV6_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv6_anal_partner] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv6_anal_partner] + sex_exposed$partvaxct_hpv6_anal[sex_exposed$HPV6_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv11_anal_ego] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv11_anal_ego] + sex_exposed$vaxct_hpv11_anal[sex_exposed$HPV11_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv11_anal_partner] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv11_anal_partner] + sex_exposed$partvaxct_hpv11_anal[sex_exposed$HPV11_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv16_anal_ego] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv16_anal_ego] + sex_exposed$vaxct_hpv16_anal[sex_exposed$HPV16_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv16_anal_partner] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv16_anal_partner] + sex_exposed$partvaxct_hpv16_anal[sex_exposed$HPV16_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv18_anal_ego] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv18_anal_ego] + sex_exposed$vaxct_hpv18_anal[sex_exposed$HPV18_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv18_anal_partner] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv18_anal_partner] + sex_exposed$partvaxct_hpv18_anal[sex_exposed$HPV18_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv31_anal_ego] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv31_anal_ego] + sex_exposed$vaxct_hpv31_anal[sex_exposed$HPV31_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv31_anal_partner] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv31_anal_partner] + sex_exposed$partvaxct_hpv31_anal[sex_exposed$HPV31_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv33_anal_ego] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv33_anal_ego] + sex_exposed$vaxct_hpv33_anal[sex_exposed$HPV33_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv33_anal_partner] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv33_anal_partner] + sex_exposed$partvaxct_hpv33_anal[sex_exposed$HPV33_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv45_anal_ego] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv45_anal_ego] + sex_exposed$vaxct_hpv45_anal[sex_exposed$HPV45_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv45_anal_partner] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv45_anal_partner] + sex_exposed$partvaxct_hpv45_anal[sex_exposed$HPV45_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv52_anal_ego] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv52_anal_ego] + sex_exposed$vaxct_hpv52_anal[sex_exposed$HPV52_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv52_anal_partner] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv52_anal_partner] + sex_exposed$partvaxct_hpv52_anal[sex_exposed$HPV52_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv58_anal_ego] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv58_anal_ego] + sex_exposed$vaxct_hpv58_anal[sex_exposed$HPV58_ANAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv58_anal_partner] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv58_anal_partner] + sex_exposed$partvaxct_hpv58_anal[sex_exposed$HPV58_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv6_anal_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv6_anal_ego] + sex_exposed$prevct_hpv6_anal[sex_exposed$HPV6_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv6_anal_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv6_anal_partner] + sex_exposed$partprevct_hpv6_anal[sex_exposed$HPV6_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv11_anal_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv11_anal_ego] + sex_exposed$prevct_hpv11_anal[sex_exposed$HPV11_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv11_anal_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv11_anal_partner] + sex_exposed$partprevct_hpv11_anal[sex_exposed$HPV11_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv16_anal_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv16_anal_ego] + sex_exposed$prevct_hpv16_anal[sex_exposed$HPV16_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv16_anal_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv16_anal_partner] + sex_exposed$partprevct_hpv16_anal[sex_exposed$HPV16_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv18_anal_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv18_anal_ego] + sex_exposed$prevct_hpv18_anal[sex_exposed$HPV18_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv18_anal_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv18_anal_partner] + sex_exposed$partprevct_hpv18_anal[sex_exposed$HPV18_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv31_anal_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv31_anal_ego] + sex_exposed$prevct_hpv31_anal[sex_exposed$HPV31_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv31_anal_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv31_anal_partner] + sex_exposed$partprevct_hpv31_anal[sex_exposed$HPV31_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv33_anal_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv33_anal_ego] + sex_exposed$prevct_hpv33_anal[sex_exposed$HPV33_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv33_anal_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv33_anal_partner] + sex_exposed$partprevct_hpv33_anal[sex_exposed$HPV33_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv45_anal_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv45_anal_ego] + sex_exposed$prevct_hpv45_anal[sex_exposed$HPV45_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv45_anal_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv45_anal_partner] + sex_exposed$partprevct_hpv45_anal[sex_exposed$HPV45_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv52_anal_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv52_anal_ego] + sex_exposed$prevct_hpv52_anal[sex_exposed$HPV52_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv52_anal_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv52_anal_partner] + sex_exposed$partprevct_hpv52_anal[sex_exposed$HPV52_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv58_anal_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv58_anal_ego] + sex_exposed$prevct_hpv58_anal[sex_exposed$HPV58_ANAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv58_anal_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv58_anal_partner] + sex_exposed$partprevct_hpv58_anal[sex_exposed$HPV58_ANAL==1]
          abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_anal_ego] = ifelse(sex_exposed$newHPV6_ANAL[sex_exposed$HPV6_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_anal_ego])
          abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_anal_partner] = ifelse(sex_exposed$partnewHPV6_ANAL[sex_exposed$HPV6_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_anal_partner])
          abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_anal_ego] = ifelse(sex_exposed$newHPV11_ANAL[sex_exposed$HPV11_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_anal_ego])
          abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_anal_partner] = ifelse(sex_exposed$partnewHPV11_ANAL[sex_exposed$HPV11_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_anal_partner])
          abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_anal_ego] = ifelse(sex_exposed$newHPV16_ANAL[sex_exposed$HPV16_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_anal_ego])
          abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_anal_partner] = ifelse(sex_exposed$partnewHPV16_ANAL[sex_exposed$HPV16_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_anal_partner])
          abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_anal_ego] = ifelse(sex_exposed$newHPV18_ANAL[sex_exposed$HPV18_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_anal_ego])
          abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_anal_partner] = ifelse(sex_exposed$partnewHPV18_ANAL[sex_exposed$HPV18_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_anal_partner])
          abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_anal_ego] = ifelse(sex_exposed$newHPV31_ANAL[sex_exposed$HPV31_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_anal_ego])
          abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_anal_partner] = ifelse(sex_exposed$partnewHPV31_ANAL[sex_exposed$HPV31_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_anal_partner])
          abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_anal_ego] = ifelse(sex_exposed$newHPV33_ANAL[sex_exposed$HPV33_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_anal_ego])
          abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_anal_partner] = ifelse(sex_exposed$partnewHPV33_ANAL[sex_exposed$HPV33_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_anal_partner])
          abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_anal_ego] = ifelse(sex_exposed$newHPV45_ANAL[sex_exposed$HPV45_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_anal_ego])
          abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_anal_partner] = ifelse(sex_exposed$partnewHPV45_ANAL[sex_exposed$HPV45_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_anal_partner])
          abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_anal_ego] = ifelse(sex_exposed$newHPV52_ANAL[sex_exposed$HPV52_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_anal_ego])
          abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_anal_partner] = ifelse(sex_exposed$partnewHPV52_ANAL[sex_exposed$HPV52_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_anal_partner])
          abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_anal_ego] = ifelse(sex_exposed$newHPV58_ANAL[sex_exposed$HPV58_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_anal_ego])
          abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_anal_partner] = ifelse(sex_exposed$partnewHPV58_ANAL[sex_exposed$HPV58_ANAL==1]==1, 1, abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_anal_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL[exposed_hpv6_anal_ego] = ifelse(sex_exposed$newHPV6_ANAL[sex_exposed$HPV6_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL[exposed_hpv6_anal_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL[exposed_hpv6_anal_partner] = ifelse(sex_exposed$partnewHPV6_ANAL[sex_exposed$HPV6_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL[exposed_hpv6_anal_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL[exposed_hpv11_anal_ego] = ifelse(sex_exposed$newHPV11_ANAL[sex_exposed$HPV11_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL[exposed_hpv11_anal_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL[exposed_hpv11_anal_partner] = ifelse(sex_exposed$partnewHPV11_ANAL[sex_exposed$HPV11_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL[exposed_hpv11_anal_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL[exposed_hpv16_anal_ego] = ifelse(sex_exposed$newHPV16_ANAL[sex_exposed$HPV16_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL[exposed_hpv16_anal_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL[exposed_hpv16_anal_partner] = ifelse(sex_exposed$partnewHPV16_ANAL[sex_exposed$HPV16_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL[exposed_hpv16_anal_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL[exposed_hpv18_anal_ego] = ifelse(sex_exposed$newHPV18_ANAL[sex_exposed$HPV18_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL[exposed_hpv18_anal_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL[exposed_hpv18_anal_partner] = ifelse(sex_exposed$partnewHPV18_ANAL[sex_exposed$HPV18_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL[exposed_hpv18_anal_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL[exposed_hpv31_anal_ego] = ifelse(sex_exposed$newHPV31_ANAL[sex_exposed$HPV31_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL[exposed_hpv31_anal_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL[exposed_hpv31_anal_partner] = ifelse(sex_exposed$partnewHPV31_ANAL[sex_exposed$HPV31_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL[exposed_hpv31_anal_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL[exposed_hpv33_anal_ego] = ifelse(sex_exposed$newHPV33_ANAL[sex_exposed$HPV33_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL[exposed_hpv33_anal_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL[exposed_hpv33_anal_partner] = ifelse(sex_exposed$partnewHPV33_ANAL[sex_exposed$HPV33_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL[exposed_hpv33_anal_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL[exposed_hpv45_anal_ego] = ifelse(sex_exposed$newHPV45_ANAL[sex_exposed$HPV45_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL[exposed_hpv45_anal_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL[exposed_hpv45_anal_partner] = ifelse(sex_exposed$partnewHPV45_ANAL[sex_exposed$HPV45_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL[exposed_hpv45_anal_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL[exposed_hpv52_anal_ego] = ifelse(sex_exposed$newHPV52_ANAL[sex_exposed$HPV52_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL[exposed_hpv52_anal_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL[exposed_hpv52_anal_partner] = ifelse(sex_exposed$partnewHPV52_ANAL[sex_exposed$HPV52_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL[exposed_hpv52_anal_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL[exposed_hpv58_anal_ego] = ifelse(sex_exposed$newHPV58_ANAL[sex_exposed$HPV58_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL[exposed_hpv58_anal_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL[exposed_hpv58_anal_partner] = ifelse(sex_exposed$partnewHPV58_ANAL[sex_exposed$HPV58_ANAL==1]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL[exposed_hpv58_anal_partner])
          
          # abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv_anal_ego] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv_anal_ego] + sex_exposed$seroct_hpv[sex_exposed$HPV_ANAL==1]
          # abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv_anal_partner] = abm_sims[[current_data]]$SERO_PREVENT[exposed_hpv_anal_partner] + sex_exposed$partseroct_hpv[sex_exposed$HPV_ANAL==1]
          # abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv_anal_ego] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv_anal_ego] + sex_exposed$condct_hpv_anal[sex_exposed$HPV_ANAL==1]
          # abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv_anal_partner] = abm_sims[[current_data]]$COND_PREVENT_ANAL[exposed_hpv_anal_partner] + sex_exposed$partcondct_hpv_anal[sex_exposed$HPV_ANAL==1]
          # abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv_anal_ego] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv_anal_ego] + sex_exposed$vaxct_hpv_anal[sex_exposed$HPV_ANAL==1]
          # abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv_anal_partner] = abm_sims[[current_data]]$VAX_PREVENT_ANAL[exposed_hpv_anal_partner] + sex_exposed$partvaxct_hpv_anal[sex_exposed$HPV_ANAL==1]
          # abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv_anal_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv_anal_ego] + sex_exposed$prevct_hpv_anal[sex_exposed$HPV_ANAL==1]
          # abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv_anal_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ANAL[exposed_hpv_anal_partner] + sex_exposed$partprevct_hpv_anal[sex_exposed$HPV_ANAL==1]
          # abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv_anal_partner]==1, 1, abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv_anal_ego]==1, 1, abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv_anal_partner]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv_anal_ego]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv_anal_partner]==1, 1, abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv_anal_ego]==1, 1, abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv_anal_partner]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv_anal_ego]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv_anal_partner]==1, 1, abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv_anal_ego]==1, 1, abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv_anal_partner]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv_anal_ego]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv_anal_partner]==1, 1, abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv_anal_ego]==1, 1, abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv_anal_partner]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv_anal_ego]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv_anal_partner]==1, 1, abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv_anal_ego]==1, 1, abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv_anal_partner]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv_anal_ego]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv_anal_partner]==1, 1, abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv_anal_ego]==1, 1, abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv_anal_partner]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv_anal_ego]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv_anal_partner]==1, 1, abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv_anal_ego]==1, 1, abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv_anal_partner]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv_anal_ego]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv_anal_partner]==1, 1, abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv_anal_ego]==1, 1, abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv_anal_partner]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv_anal_ego]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv_anal_partner]==1, 1, abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv_anal_ego]==1, 1, abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv_anal_partner])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL[exposed_hpv_anal_ego] = ifelse(sex_exposed$newHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv_anal_partner]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL[exposed_hpv_anal_ego])
          # abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL[exposed_hpv_anal_partner] = ifelse(sex_exposed$partnewHPV_ANAL[sex_exposed$HPV_ANAL==1]==1 & abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv_anal_ego]==1, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL[exposed_hpv_anal_partner])
          
          #HPV oral
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv6_oral_ego] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv6_oral_ego] + sex_exposed$condct_hpv6_oral[sex_exposed$HPV6_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv6_oral_partner] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv6_oral_partner] + sex_exposed$partcondct_hpv6_oral[sex_exposed$HPV6_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv11_oral_ego] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv11_oral_ego] + sex_exposed$condct_hpv11_oral[sex_exposed$HPV11_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv11_oral_partner] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv11_oral_partner] + sex_exposed$partcondct_hpv11_oral[sex_exposed$HPV11_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv16_oral_ego] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv16_oral_ego] + sex_exposed$condct_hpv16_oral[sex_exposed$HPV16_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv16_oral_partner] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv16_oral_partner] + sex_exposed$partcondct_hpv16_oral[sex_exposed$HPV16_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv18_oral_ego] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv18_oral_ego] + sex_exposed$condct_hpv18_oral[sex_exposed$HPV18_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv18_oral_partner] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv18_oral_partner] + sex_exposed$partcondct_hpv18_oral[sex_exposed$HPV18_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv31_oral_ego] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv31_oral_ego] + sex_exposed$condct_hpv31_oral[sex_exposed$HPV31_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv31_oral_partner] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv31_oral_partner] + sex_exposed$partcondct_hpv31_oral[sex_exposed$HPV31_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv33_oral_ego] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv33_oral_ego] + sex_exposed$condct_hpv33_oral[sex_exposed$HPV33_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv33_oral_partner] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv33_oral_partner] + sex_exposed$partcondct_hpv33_oral[sex_exposed$HPV33_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv45_oral_ego] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv45_oral_ego] + sex_exposed$condct_hpv45_oral[sex_exposed$HPV45_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv45_oral_partner] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv45_oral_partner] + sex_exposed$partcondct_hpv45_oral[sex_exposed$HPV45_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv52_oral_ego] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv52_oral_ego] + sex_exposed$condct_hpv52_oral[sex_exposed$HPV52_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv52_oral_partner] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv52_oral_partner] + sex_exposed$partcondct_hpv52_oral[sex_exposed$HPV52_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv58_oral_ego] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv58_oral_ego] + sex_exposed$condct_hpv58_oral[sex_exposed$HPV58_ORAL==1]
          abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv58_oral_partner] = abm_sims[[current_data]]$COND_PREVENT_ORAL[exposed_hpv58_oral_partner] + sex_exposed$partcondct_hpv58_oral[sex_exposed$HPV58_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv6_oral_ego] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv6_oral_ego] + sex_exposed$vaxct_hpv6_oral[sex_exposed$HPV6_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv6_oral_partner] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv6_oral_partner] + sex_exposed$partvaxct_hpv6_oral[sex_exposed$HPV6_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv11_oral_ego] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv11_oral_ego] + sex_exposed$vaxct_hpv11_oral[sex_exposed$HPV11_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv11_oral_partner] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv11_oral_partner] + sex_exposed$partvaxct_hpv11_oral[sex_exposed$HPV11_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv16_oral_ego] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv16_oral_ego] + sex_exposed$vaxct_hpv16_oral[sex_exposed$HPV16_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv16_oral_partner] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv16_oral_partner] + sex_exposed$partvaxct_hpv16_oral[sex_exposed$HPV16_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv18_oral_ego] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv18_oral_ego] + sex_exposed$vaxct_hpv18_oral[sex_exposed$HPV18_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv18_oral_partner] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv18_oral_partner] + sex_exposed$partvaxct_hpv18_oral[sex_exposed$HPV18_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv31_oral_ego] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv31_oral_ego] + sex_exposed$vaxct_hpv31_oral[sex_exposed$HPV31_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv31_oral_partner] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv31_oral_partner] + sex_exposed$partvaxct_hpv31_oral[sex_exposed$HPV31_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv33_oral_ego] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv33_oral_ego] + sex_exposed$vaxct_hpv33_oral[sex_exposed$HPV33_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv33_oral_partner] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv33_oral_partner] + sex_exposed$partvaxct_hpv33_oral[sex_exposed$HPV33_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv45_oral_ego] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv45_oral_ego] + sex_exposed$vaxct_hpv45_oral[sex_exposed$HPV45_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv45_oral_partner] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv45_oral_partner] + sex_exposed$partvaxct_hpv45_oral[sex_exposed$HPV45_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv52_oral_ego] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv52_oral_ego] + sex_exposed$vaxct_hpv52_oral[sex_exposed$HPV52_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv52_oral_partner] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv52_oral_partner] + sex_exposed$partvaxct_hpv52_oral[sex_exposed$HPV52_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv58_oral_ego] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv58_oral_ego] + sex_exposed$vaxct_hpv58_oral[sex_exposed$HPV58_ORAL==1]
          abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv58_oral_partner] = abm_sims[[current_data]]$VAX_PREVENT_ORAL[exposed_hpv58_oral_partner] + sex_exposed$partvaxct_hpv58_oral[sex_exposed$HPV58_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv6_oral_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv6_oral_ego] + sex_exposed$prevct_hpv6_oral[sex_exposed$HPV6_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv6_oral_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv6_oral_partner] + sex_exposed$partprevct_hpv6_oral[sex_exposed$HPV6_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv11_oral_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv11_oral_ego] + sex_exposed$prevct_hpv11_oral[sex_exposed$HPV11_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv11_oral_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv11_oral_partner] + sex_exposed$partprevct_hpv11_oral[sex_exposed$HPV11_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv16_oral_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv16_oral_ego] + sex_exposed$prevct_hpv16_oral[sex_exposed$HPV16_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv16_oral_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv16_oral_partner] + sex_exposed$partprevct_hpv16_oral[sex_exposed$HPV16_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv18_oral_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv18_oral_ego] + sex_exposed$prevct_hpv18_oral[sex_exposed$HPV18_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv18_oral_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv18_oral_partner] + sex_exposed$partprevct_hpv18_oral[sex_exposed$HPV18_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv31_oral_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv31_oral_ego] + sex_exposed$prevct_hpv31_oral[sex_exposed$HPV31_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv31_oral_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv31_oral_partner] + sex_exposed$partprevct_hpv31_oral[sex_exposed$HPV31_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv33_oral_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv33_oral_ego] + sex_exposed$prevct_hpv33_oral[sex_exposed$HPV33_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv33_oral_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv33_oral_partner] + sex_exposed$partprevct_hpv33_oral[sex_exposed$HPV33_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv45_oral_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv45_oral_ego] + sex_exposed$prevct_hpv45_oral[sex_exposed$HPV45_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv45_oral_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv45_oral_partner] + sex_exposed$partprevct_hpv45_oral[sex_exposed$HPV45_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv52_oral_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv52_oral_ego] + sex_exposed$prevct_hpv52_oral[sex_exposed$HPV52_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv52_oral_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv52_oral_partner] + sex_exposed$partprevct_hpv52_oral[sex_exposed$HPV52_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv58_oral_ego] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv58_oral_ego] + sex_exposed$prevct_hpv58_oral[sex_exposed$HPV58_ORAL==1]
          abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv58_oral_partner] = abm_sims[[current_data]]$OVERALL_PREVENT_ORAL[exposed_hpv58_oral_partner] + sex_exposed$partprevct_hpv58_oral[sex_exposed$HPV58_ORAL==1]
          abm_sims[[current_data]]$HPV6_ORAL[exposed_hpv6_oral_ego] = ifelse(sex_exposed$newHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==1 | sex_exposed$newHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV6_ORAL[exposed_hpv6_oral_ego])
          abm_sims[[current_data]]$HPV6_ORAL[exposed_hpv6_oral_partner] = ifelse(sex_exposed$partnewHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==1 | sex_exposed$partnewHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV6_ORAL[exposed_hpv6_oral_partner])
          abm_sims[[current_data]]$HPV11_ORAL[exposed_hpv11_oral_ego] = ifelse(sex_exposed$newHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==1 | sex_exposed$newHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV11_ORAL[exposed_hpv11_oral_ego])
          abm_sims[[current_data]]$HPV11_ORAL[exposed_hpv11_oral_partner] = ifelse(sex_exposed$partnewHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==1 | sex_exposed$partnewHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV11_ORAL[exposed_hpv11_oral_partner])
          abm_sims[[current_data]]$HPV16_ORAL[exposed_hpv16_oral_ego] = ifelse(sex_exposed$newHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==1 | sex_exposed$newHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV16_ORAL[exposed_hpv16_oral_ego])
          abm_sims[[current_data]]$HPV16_ORAL[exposed_hpv16_oral_partner] = ifelse(sex_exposed$partnewHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==1 | sex_exposed$partnewHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV16_ORAL[exposed_hpv16_oral_partner])
          abm_sims[[current_data]]$HPV18_ORAL[exposed_hpv18_oral_ego] = ifelse(sex_exposed$newHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==1 | sex_exposed$newHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV18_ORAL[exposed_hpv18_oral_ego])
          abm_sims[[current_data]]$HPV18_ORAL[exposed_hpv18_oral_partner] = ifelse(sex_exposed$partnewHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==1 | sex_exposed$partnewHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV18_ORAL[exposed_hpv18_oral_partner])
          abm_sims[[current_data]]$HPV31_ORAL[exposed_hpv31_oral_ego] = ifelse(sex_exposed$newHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==1 | sex_exposed$newHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV31_ORAL[exposed_hpv31_oral_ego])
          abm_sims[[current_data]]$HPV31_ORAL[exposed_hpv31_oral_partner] = ifelse(sex_exposed$partnewHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==1 | sex_exposed$partnewHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV31_ORAL[exposed_hpv31_oral_partner])
          abm_sims[[current_data]]$HPV33_ORAL[exposed_hpv33_oral_ego] = ifelse(sex_exposed$newHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==1 | sex_exposed$newHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV33_ORAL[exposed_hpv33_oral_ego])
          abm_sims[[current_data]]$HPV33_ORAL[exposed_hpv33_oral_partner] = ifelse(sex_exposed$partnewHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==1 | sex_exposed$partnewHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV33_ORAL[exposed_hpv33_oral_partner])
          abm_sims[[current_data]]$HPV45_ORAL[exposed_hpv45_oral_ego] = ifelse(sex_exposed$newHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==1 | sex_exposed$newHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV45_ORAL[exposed_hpv45_oral_ego])
          abm_sims[[current_data]]$HPV45_ORAL[exposed_hpv45_oral_partner] = ifelse(sex_exposed$partnewHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==1 | sex_exposed$partnewHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV45_ORAL[exposed_hpv45_oral_partner])
          abm_sims[[current_data]]$HPV52_ORAL[exposed_hpv52_oral_ego] = ifelse(sex_exposed$newHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==1 | sex_exposed$newHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV52_ORAL[exposed_hpv52_oral_ego])
          abm_sims[[current_data]]$HPV52_ORAL[exposed_hpv52_oral_partner] = ifelse(sex_exposed$partnewHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==1 | sex_exposed$partnewHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV52_ORAL[exposed_hpv52_oral_partner])
          abm_sims[[current_data]]$HPV58_ORAL[exposed_hpv58_oral_ego] = ifelse(sex_exposed$newHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==1 | sex_exposed$newHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV58_ORAL[exposed_hpv58_oral_ego])
          abm_sims[[current_data]]$HPV58_ORAL[exposed_hpv58_oral_partner] = ifelse(sex_exposed$partnewHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==1 | sex_exposed$partnewHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV58_ORAL[exposed_hpv58_oral_partner])
          
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ORAL[exposed_hpv6_oral_ego] = ifelse(sex_exposed$newHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==1 | sex_exposed$newHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ORAL[exposed_hpv6_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ORAL[exposed_hpv6_oral_partner] = ifelse(sex_exposed$partnewHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==1 | sex_exposed$partnewHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ORAL[exposed_hpv6_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ORAL[exposed_hpv11_oral_ego] = ifelse(sex_exposed$newHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==1 | sex_exposed$newHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ORAL[exposed_hpv11_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ORAL[exposed_hpv11_oral_partner] = ifelse(sex_exposed$partnewHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==1 | sex_exposed$partnewHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ORAL[exposed_hpv11_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ORAL[exposed_hpv16_oral_ego] = ifelse(sex_exposed$newHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==1 | sex_exposed$newHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ORAL[exposed_hpv16_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ORAL[exposed_hpv16_oral_partner] = ifelse(sex_exposed$partnewHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==1 | sex_exposed$partnewHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ORAL[exposed_hpv16_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ORAL[exposed_hpv18_oral_ego] = ifelse(sex_exposed$newHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==1 | sex_exposed$newHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ORAL[exposed_hpv18_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ORAL[exposed_hpv18_oral_partner] = ifelse(sex_exposed$partnewHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==1 | sex_exposed$partnewHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ORAL[exposed_hpv18_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ORAL[exposed_hpv31_oral_ego] = ifelse(sex_exposed$newHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==1 | sex_exposed$newHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ORAL[exposed_hpv31_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ORAL[exposed_hpv31_oral_partner] = ifelse(sex_exposed$partnewHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==1 | sex_exposed$partnewHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ORAL[exposed_hpv31_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ORAL[exposed_hpv33_oral_ego] = ifelse(sex_exposed$newHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==1 | sex_exposed$newHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ORAL[exposed_hpv33_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ORAL[exposed_hpv33_oral_partner] = ifelse(sex_exposed$partnewHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==1 | sex_exposed$partnewHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ORAL[exposed_hpv33_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ORAL[exposed_hpv45_oral_ego] = ifelse(sex_exposed$newHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==1 | sex_exposed$newHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ORAL[exposed_hpv45_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ORAL[exposed_hpv45_oral_partner] = ifelse(sex_exposed$partnewHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==1 | sex_exposed$partnewHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ORAL[exposed_hpv45_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ORAL[exposed_hpv52_oral_ego] = ifelse(sex_exposed$newHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==1 | sex_exposed$newHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ORAL[exposed_hpv52_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ORAL[exposed_hpv52_oral_partner] = ifelse(sex_exposed$partnewHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==1 | sex_exposed$partnewHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ORAL[exposed_hpv52_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ORAL[exposed_hpv58_oral_ego] = ifelse(sex_exposed$newHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==1 | sex_exposed$newHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ORAL[exposed_hpv58_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ORAL[exposed_hpv58_oral_partner] = ifelse(sex_exposed$partnewHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==1 | sex_exposed$partnewHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ORAL[exposed_hpv58_oral_partner])
          
          abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_oral_ego] = ifelse(sex_exposed$newHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==2 | sex_exposed$newHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_oral_ego])
          abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_oral_partner] = ifelse(sex_exposed$partnewHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==2 | sex_exposed$partnewHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV6_ANAL[exposed_hpv6_oral_partner])
          abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_oral_ego] = ifelse(sex_exposed$newHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==2 | sex_exposed$newHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_oral_ego])
          abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_oral_partner] = ifelse(sex_exposed$partnewHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==2 | sex_exposed$partnewHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV11_ANAL[exposed_hpv11_oral_partner])
          abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_oral_ego] = ifelse(sex_exposed$newHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==2 | sex_exposed$newHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_oral_ego])
          abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_oral_partner] = ifelse(sex_exposed$partnewHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==2 | sex_exposed$partnewHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV16_ANAL[exposed_hpv16_oral_partner])
          abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_oral_ego] = ifelse(sex_exposed$newHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==2 | sex_exposed$newHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_oral_ego])
          abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_oral_partner] = ifelse(sex_exposed$partnewHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==2 | sex_exposed$partnewHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV18_ANAL[exposed_hpv18_oral_partner])
          abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_oral_ego] = ifelse(sex_exposed$newHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==2 | sex_exposed$newHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_oral_ego])
          abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_oral_partner] = ifelse(sex_exposed$partnewHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==2 | sex_exposed$partnewHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV31_ANAL[exposed_hpv31_oral_partner])
          abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_oral_ego] = ifelse(sex_exposed$newHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==2 | sex_exposed$newHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_oral_ego])
          abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_oral_partner] = ifelse(sex_exposed$partnewHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==2 | sex_exposed$partnewHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV33_ANAL[exposed_hpv33_oral_partner])
          abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_oral_ego] = ifelse(sex_exposed$newHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==2 | sex_exposed$newHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_oral_ego])
          abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_oral_partner] = ifelse(sex_exposed$partnewHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==2 | sex_exposed$partnewHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV45_ANAL[exposed_hpv45_oral_partner])
          abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_oral_ego] = ifelse(sex_exposed$newHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==2 | sex_exposed$newHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_oral_ego])
          abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_oral_partner] = ifelse(sex_exposed$partnewHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==2 | sex_exposed$partnewHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV52_ANAL[exposed_hpv52_oral_partner])
          abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_oral_ego] = ifelse(sex_exposed$newHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==2 | sex_exposed$newHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_oral_ego])
          abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_oral_partner] = ifelse(sex_exposed$partnewHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==2 | sex_exposed$partnewHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==3, 1, abm_sims[[current_data]]$HPV58_ANAL[exposed_hpv58_oral_partner])
          
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL[exposed_hpv6_oral_ego] = ifelse(sex_exposed$newHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==2 | sex_exposed$newHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL[exposed_hpv6_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL[exposed_hpv6_oral_partner] = ifelse(sex_exposed$partnewHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==2 | sex_exposed$partnewHPV6_ORAL[sex_exposed$HPV6_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL[exposed_hpv6_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL[exposed_hpv11_oral_ego] = ifelse(sex_exposed$newHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==2 | sex_exposed$newHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL[exposed_hpv11_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL[exposed_hpv11_oral_partner] = ifelse(sex_exposed$partnewHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==2 | sex_exposed$partnewHPV11_ORAL[sex_exposed$HPV11_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL[exposed_hpv11_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL[exposed_hpv16_oral_ego] = ifelse(sex_exposed$newHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==2 | sex_exposed$newHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL[exposed_hpv16_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL[exposed_hpv16_oral_partner] = ifelse(sex_exposed$partnewHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==2 | sex_exposed$partnewHPV16_ORAL[sex_exposed$HPV16_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL[exposed_hpv16_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL[exposed_hpv18_oral_ego] = ifelse(sex_exposed$newHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==2 | sex_exposed$newHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL[exposed_hpv18_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL[exposed_hpv18_oral_partner] = ifelse(sex_exposed$partnewHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==2 | sex_exposed$partnewHPV18_ORAL[sex_exposed$HPV18_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL[exposed_hpv18_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL[exposed_hpv31_oral_ego] = ifelse(sex_exposed$newHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==2 | sex_exposed$newHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL[exposed_hpv31_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL[exposed_hpv31_oral_partner] = ifelse(sex_exposed$partnewHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==2 | sex_exposed$partnewHPV31_ORAL[sex_exposed$HPV31_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL[exposed_hpv31_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL[exposed_hpv33_oral_ego] = ifelse(sex_exposed$newHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==2 | sex_exposed$newHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL[exposed_hpv33_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL[exposed_hpv33_oral_partner] = ifelse(sex_exposed$partnewHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==2 | sex_exposed$partnewHPV33_ORAL[sex_exposed$HPV33_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL[exposed_hpv33_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL[exposed_hpv45_oral_ego] = ifelse(sex_exposed$newHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==2 | sex_exposed$newHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL[exposed_hpv45_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL[exposed_hpv45_oral_partner] = ifelse(sex_exposed$partnewHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==2 | sex_exposed$partnewHPV45_ORAL[sex_exposed$HPV45_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL[exposed_hpv45_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL[exposed_hpv52_oral_ego] = ifelse(sex_exposed$newHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==2 | sex_exposed$newHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL[exposed_hpv52_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL[exposed_hpv52_oral_partner] = ifelse(sex_exposed$partnewHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==2 | sex_exposed$partnewHPV52_ORAL[sex_exposed$HPV52_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL[exposed_hpv52_oral_partner])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL[exposed_hpv58_oral_ego] = ifelse(sex_exposed$newHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==2 | sex_exposed$newHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL[exposed_hpv58_oral_ego])
          abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL[exposed_hpv58_oral_partner] = ifelse(sex_exposed$partnewHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==2 | sex_exposed$partnewHPV58_ORAL[sex_exposed$HPV58_ORAL==1]==3, 0, abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL[exposed_hpv58_oral_partner])
          
          #update number of sex acts
          abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Ego] = abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Ego] + 1
          abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Partner] = abm_sims[[current_data]]$SEX_COUNT[sex_exposed$Partner] + 1
          
          if (length(sex_unexposed)>0)
            abm_sims[[current_data]]$SEX_COUNT[sex_unexposed] = abm_sims[[current_data]]$SEX_COUNT[sex_unexposed] + 1
          
          ## CLEAN UP ##
          
          rm(sex_exposed,sex_unexposed,exposed_hiv_ego,exposed_hiv_partner,probinfect1_hiv,probinfect2_hiv,exposed_hpv11_anal_ego,exposed_hpv11_anal_partner,exposed_hpv16_anal_ego,exposed_hpv16_anal_partner,exposed_hpv18_anal_ego,exposed_hpv18_anal_partner,exposed_hpv31_anal_ego,exposed_hpv31_anal_partner,exposed_hpv33_anal_ego,exposed_hpv33_anal_partner,exposed_hpv45_anal_ego,exposed_hpv45_anal_partner,exposed_hpv52_anal_ego,exposed_hpv52_anal_partner,exposed_hpv58_anal_ego,exposed_hpv58_anal_partner,exposed_hpv6_anal_ego,exposed_hpv6_anal_partner,probinfect1_hpv11_anal,probinfect1_hpv16_anal,probinfect1_hpv18_anal,probinfect1_hpv31_anal,probinfect1_hpv33_anal,probinfect1_hpv45_anal,probinfect1_hpv52_anal,probinfect1_hpv58_anal,probinfect1_hpv6_anal,probinfect2_hpv11_anal,probinfect2_hpv16_anal,probinfect2_hpv18_anal,probinfect2_hpv31_anal,probinfect2_hpv33_anal,probinfect2_hpv45_anal,probinfect2_hpv52_anal,probinfect2_hpv58_anal,probinfect2_hpv6_anal,exposed_hpv33_oral_ego,exposed_hpv45_oral_ego,exposed_hpv52_oral_ego,exposed_hpv58_oral_ego,exposed_hpv6_oral_ego,exposed_hpv11_oral_ego,exposed_hpv16_oral_ego,exposed_hpv18_oral_ego,exposed_hpv31_oral_ego,exposed_hpv33_oral_partner,exposed_hpv45_oral_partner,exposed_hpv52_oral_partner,exposed_hpv58_oral_partner,exposed_hpv6_oral_partner,exposed_hpv11_oral_partner,exposed_hpv16_oral_partner,exposed_hpv18_oral_partner,exposed_hpv31_oral_partner,probinfect1_hpv11_oral,probinfect1_hpv18_oral,probinfect1_hpv33_oral,probinfect1_hpv52_oral,probinfect1_hpv6_oral,probinfect2_hpv16_oral,probinfect2_hpv31_oral,probinfect2_hpv45_oral,probinfect2_hpv58_oral,probinfect1_hpv16_oral,probinfect1_hpv31_oral,probinfect1_hpv45_oral,probinfect1_hpv58_oral,probinfect2_hpv11_oral,probinfect2_hpv18_oral,probinfect2_hpv33_oral,probinfect2_hpv52_oral,probinfect2_hpv6_oral)
          #gc()
          
        } else {
          
          #pathogen concordance in every partnership
          
          ## UPDATE STATS IN MAIN DATASET ##
          
          abm_sims[[current_data]]$SEX_COUNT[sex_unexposed] = abm_sims[[current_data]]$SEX_COUNT[sex_unexposed] + 1
          
          ## CLEAN UP ##
          
          rm(sex_exposed,sex_unexposed)
          #gc()
          
        }
        
      } else {
        
        #no partnerships this day
        rm(sex_selection)   
        
      }
      
      #track longitudinal data using incidence days as indicator of new infection
      long_sims[[current_data]]$HIV[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HIV==0)
      long_sims[[current_data]]$HPV6_ANAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ANAL==0)
      long_sims[[current_data]]$HPV6_ORAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV6_ORAL==0)
      long_sims[[current_data]]$HPV11_ANAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ANAL==0)
      long_sims[[current_data]]$HPV11_ORAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV11_ORAL==0)
      long_sims[[current_data]]$HPV16_ANAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ANAL==0)
      long_sims[[current_data]]$HPV16_ORAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV16_ORAL==0)
      long_sims[[current_data]]$HPV18_ANAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ANAL==0)
      long_sims[[current_data]]$HPV18_ORAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV18_ORAL==0)
      long_sims[[current_data]]$HPV31_ANAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ANAL==0)
      long_sims[[current_data]]$HPV31_ORAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV31_ORAL==0)
      long_sims[[current_data]]$HPV33_ANAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ANAL==0)
      long_sims[[current_data]]$HPV33_ORAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV33_ORAL==0)
      long_sims[[current_data]]$HPV45_ANAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ANAL==0)
      long_sims[[current_data]]$HPV45_ORAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV45_ORAL==0)
      long_sims[[current_data]]$HPV52_ANAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ANAL==0)
      long_sims[[current_data]]$HPV52_ORAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV52_ORAL==0)
      long_sims[[current_data]]$HPV58_ANAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ANAL==0)
      long_sims[[current_data]]$HPV58_ORAL[day] = sum(abm_sims[[current_data]]$INCIDENCE_DAYS_HPV58_ORAL==0)
      
    }
    
  }
  rm(day,current_data,active)
  
  #return this simulation to the list
  abm_sims_50[((sims-1)*nparadigms+1):(((sims-1)*nparadigms)+nparadigms)] = abm_sims
  long_sims_50[((sims-1)*nparadigms+1):(((sims-1)*nparadigms)+nparadigms)] = long_sims
  
}
rm(sims,abm_sims,sex_network,network_sims_50,long_sims,seed_val_day)
gc()

#benchmarking
#time3 = Sys.time()


### SAVE SIMULATION RESULTS ###

#save.image("simulation results 045 calibration.RData")
save.image("simulation results 045.RData")
#save.image("simulation results 045 sens vax95.RData")
#save.image("simulation results 045 sens condom50.RData")


### STEP5: CALIBRATION ###

load("simulation results 0.45 calibration.RData")

#create an annual HIV incidence dataframe
annual_hiv = data.frame(matrix(data=NA, nrow=0, ncol=(length_sim/365)), stringsAsFactors=F)
for (i in seq(1,nsims*nparadigms,by=nparadigms))
{
  long_sims_50[[i]]$YEAR = (floor((long_sims_50[[i]]$DAY-1)/365) + 1)
  annual_hiv = rbind(annual_hiv, matrix(data=as.numeric(by(long_sims_50[[i]]$HIV, long_sims_50[[i]]$YEAR, FUN=sum)), nrow=1, ncol=length(unique(long_sims_50[[i]]$YEAR))))
  
}
rm(i)

colMeans(annual_hiv)

#plot predicted vs surveillance data

years = 2013:2017
surveillance = c(149,138,136,124,126)
#output to tif for publication 
#tiff("Figure1.tif",height=6,width=10,units='in',res=1200) 
plot(x=years, y=surveillance, ylim=c(0,250), pch=18, cex=2, xlab="Surveillance Year", ylab="Reported Cases")
if (ncol(annual_hiv)>5)
{
  annual_hiv = annual_hiv[,(ncol(annual_hiv)-4):ncol(annual_hiv)]
}
for (i in 1:length(surveillance)) 
{
  arrows(x0=years[i], y0=min(annual_hiv[,i]), x1=years[i], y1=max(annual_hiv[,i]), angle=90, length=0.05, code=3, lwd=2)
}
legend("topright", c("surveillance","predicted"), pch=c(18,NA), lty=c(NA,1))
#close file 
#dev.off() 

#check amount of sex
summary(abm_sims_50[[1]]$MAX_SEX)
summary(abm_sims_50[[1]]$SEX_COUNT)
hist(abm_sims_50[[1]]$MAX_SEX,breaks="fd")
hist(abm_sims_50[[1]]$SEX_COUNT,breaks="fd")
sum(abm_sims_50[[1]]$SEX_COUNT <= abm_sims_50[[1]]$MAX_SEX)
sum(abm_sims_50[[1]]$SEX_COUNT > abm_sims_50[[1]]$MAX_SEX)
boxplot(abm_sims_50[[1]]$SEX_COUNT/abm_sims_50[[1]]$MAX_SEX)
hist(abm_sims_50[[1]]$PARTNERS,breaks="fd")


### STEP6: TALLY RESULTS for PAPER ###

load("simulation results 045.RData")
#load("simulation results 045 sens vax95.RData")
#load("simulation results 045 sens condom50.RData")

#scenario labels
hpv_type = c(6, 11, 16, 18, 31, 33, 45, 52, 58)
vax_level = c(0, 13, 25, 50, 80, 100)

#create results data frame
results_infections_anal = data.frame("HPV_type"=NA,"Vax_level"=NA,"Prev_mean"=NA,"Prev_SE"=NA,"Incidence_mean"=NA,"Incidence_SE"=NA,"New_mean"=NA,"New_SE"=NA,"Total_mean"=NA,"Total_SE"=NA,"Prevent_mean"=NA,"Prevent_SE"=NA,"Vax_prevent_mean"=NA,"Vax_prevent_SE"=NA,"Sero_prevent_mean"=NA,"Sero_prevent_SE"=NA,"Cond_prevent_mean"=NA,"Cond_prevent_SE"=NA,stringsAsFactors=F)
results_infections_oral = data.frame("HPV_type"=NA,"Vax_level"=NA,"Prev_mean"=NA,"Prev_SE"=NA,"Incidence_mean"=NA,"Incidence_SE"=NA,"New_mean"=NA,"New_SE"=NA,"Total_mean"=NA,"Total_SE"=NA,"Prevent_mean"=NA,"Prevent_SE"=NA,"Vax_prevent_mean"=NA,"Vax_prevent_SE"=NA,"Cond_prevent_mean"=NA,"Cond_prevent_SE"=NA,stringsAsFactors=F)
results_longitudinal = data.frame("HPV_type"=NA,"Vax_level"=NA,"Yr1_anal_new_mean"=NA,"Yr1_anal_new_SE"=NA,"Yr5_anal_new_mean"=NA,"Yr5_anal_new_SE"=NA,"Yr10_anal_new_mean"=NA,"Yr10_anal_new_SE"=NA,"Yr1_oral_new_mean"=NA,"Yr1_oral_new_SE"=NA,"Yr5_oral_new_mean"=NA,"Yr5_oral_new_SE"=NA,"Yr10_oral_new_mean"=NA,"Yr10_oral_new_SE"=NA,stringsAsFactors=F)
for (i in 1:length(hpv_type))
{
  for (j in 1:length(vax_level))
  {
    stat_prev_anal = NA         #point prevalance at post simulation
    stat_incidence_anal = NA    #cumulative incidence (new + spontaneous clearance)
    stat_new_anal = NA          #new + spontaneous clearance
    stat_total_anal = NA        #infections at post simulation
    stat_prevent_anal = NA      #total preventions (aggregate across all serotypes)
    stat_prevent_vax_anal = NA  #total vax preventions (aggregate across all serotypes)
    stat_prevent_sero_anal = NA #total sero preventions (aggregate across all serotypes)
    stat_prevent_cond_anal = NA #total condom preventions (aggregate across all serotypes)
    stat_yr1_new_anal = NA      #1 yr new infections
    stat_yr5_new_anal = NA      #5 yr new infections
    stat_yr10_new_anal = NA      #10 yr new infections
    stat_prev_oral = NA         #point prevalance at post simulation
    stat_incidence_oral = NA    #cumulative incidence (new + spontaneous clearance)
    stat_new_oral = NA          #new + spontaneous clearance
    stat_total_oral = NA        #infections at post simulation
    stat_prevent_oral = NA      #total preventions (aggregate across all serotypes)
    stat_prevent_vax_oral = NA  #total vax preventions (aggregate across all serotypes)
    stat_prevent_cond_oral = NA #total condom preventions (aggregate across all serotypes)
    stat_yr1_new_oral = NA      #1 yr new infections
    stat_yr5_new_oral = NA      #5 yr new infections
    stat_yr10_new_oral = NA      #10 yr new infections
    
    for (k in 0:(nsims-1))
    {
      current = k*length(vax_level)
      
      if (hpv_type[i]==6) {
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV6_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV6_ANAL) + sum(abm_sims_50[[j+current]]$HPV6_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV6_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV6_ANAL) + sum(abm_sims_50[[j+current]]$HPV6_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV6_ANAL)))
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV6_ANAL)))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
        stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV6_ANAL[1:(365*1)])))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV6_ANAL[1:(365*5)])))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV6_ANAL[1:(365*10)])))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV6_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV6_ORAL) + sum(abm_sims_50[[j+current]]$HPV6_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV6_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV6_ORAL) + sum(abm_sims_50[[j+current]]$HPV6_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV6_ORAL)))
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV6_ORAL)))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
        stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV6_ORAL[1:(365*1)])))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV6_ORAL[1:(365*5)])))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV6_ORAL[1:(365*10)])))
      } else if (hpv_type[i]==11) {
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV11_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV11_ANAL) + sum(abm_sims_50[[j+current]]$HPV11_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV11_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV11_ANAL) + sum(abm_sims_50[[j+current]]$HPV11_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV11_ANAL)))
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV11_ANAL)))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
        stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV11_ANAL[1:(365*1)])))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV11_ANAL[1:(365*5)])))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV11_ANAL[1:(365*10)])))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV11_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV11_ORAL) + sum(abm_sims_50[[j+current]]$HPV11_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV11_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV11_ORAL) + sum(abm_sims_50[[j+current]]$HPV11_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV11_ORAL)))
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV11_ORAL)))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
        stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV11_ORAL[1:(365*1)])))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV11_ORAL[1:(365*5)])))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV11_ORAL[1:(365*10)])))
      } else if (hpv_type[i]==16) {
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV16_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV16_ANAL) + sum(abm_sims_50[[j+current]]$HPV16_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV16_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV16_ANAL) + sum(abm_sims_50[[j+current]]$HPV16_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV16_ANAL)))
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV16_ANAL)))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
        stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV16_ANAL[1:(365*1)])))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV16_ANAL[1:(365*5)])))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV16_ANAL[1:(365*10)])))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV16_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV16_ORAL) + sum(abm_sims_50[[j+current]]$HPV16_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV16_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV16_ORAL) + sum(abm_sims_50[[j+current]]$HPV16_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV16_ORAL)))
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV16_ORAL)))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
        stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV16_ORAL[1:(365*1)])))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV16_ORAL[1:(365*5)])))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV16_ORAL[1:(365*10)])))
      } else if (hpv_type[i]==18) {
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV18_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV18_ANAL) + sum(abm_sims_50[[j+current]]$HPV18_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV18_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV18_ANAL) + sum(abm_sims_50[[j+current]]$HPV18_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV18_ANAL)))
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV18_ANAL)))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
        stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV18_ANAL[1:(365*1)])))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV18_ANAL[1:(365*5)])))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV18_ANAL[1:(365*10)])))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV18_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV18_ORAL) + sum(abm_sims_50[[j+current]]$HPV18_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV18_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV18_ORAL) + sum(abm_sims_50[[j+current]]$HPV18_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV18_ORAL)))
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV18_ORAL)))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
        stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV18_ORAL[1:(365*1)])))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV18_ORAL[1:(365*5)])))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV18_ORAL[1:(365*10)])))
      } else if (hpv_type[i]==31) {
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV31_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV31_ANAL) + sum(abm_sims_50[[j+current]]$HPV31_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV31_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV31_ANAL) + sum(abm_sims_50[[j+current]]$HPV31_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV31_ANAL)))
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV31_ANAL)))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
        stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV31_ANAL[1:(365*1)])))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV31_ANAL[1:(365*5)])))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV31_ANAL[1:(365*10)])))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV31_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV31_ORAL) + sum(abm_sims_50[[j+current]]$HPV31_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV31_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV31_ORAL) + sum(abm_sims_50[[j+current]]$HPV31_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV31_ORAL)))
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV31_ORAL)))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
        stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV31_ORAL[1:(365*1)])))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV31_ORAL[1:(365*5)])))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV31_ORAL[1:(365*10)])))
      } else if (hpv_type[i]==33) {
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV33_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV33_ANAL) + sum(abm_sims_50[[j+current]]$HPV33_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV33_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV33_ANAL) + sum(abm_sims_50[[j+current]]$HPV33_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV33_ANAL)))
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV33_ANAL)))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
        stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV33_ANAL[1:(365*1)])))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV33_ANAL[1:(365*5)])))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV33_ANAL[1:(365*10)])))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV33_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV33_ORAL) + sum(abm_sims_50[[j+current]]$HPV33_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV33_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV33_ORAL) + sum(abm_sims_50[[j+current]]$HPV33_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV33_ORAL)))
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV33_ORAL)))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
        stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV33_ORAL[1:(365*1)])))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV33_ORAL[1:(365*5)])))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV33_ORAL[1:(365*10)])))
      } else if (hpv_type[i]==45) {
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV45_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV45_ANAL) + sum(abm_sims_50[[j+current]]$HPV45_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV45_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV45_ANAL) + sum(abm_sims_50[[j+current]]$HPV45_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV45_ANAL)))
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV45_ANAL)))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
        stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV45_ANAL[1:(365*1)])))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV45_ANAL[1:(365*5)])))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV45_ANAL[1:(365*10)])))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV45_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV45_ORAL) + sum(abm_sims_50[[j+current]]$HPV45_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV45_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV45_ORAL) + sum(abm_sims_50[[j+current]]$HPV45_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV45_ORAL)))
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV45_ORAL)))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
        stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV45_ORAL[1:(365*1)])))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV45_ORAL[1:(365*5)])))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV45_ORAL[1:(365*10)])))
      } else if (hpv_type[i]==52) {
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV52_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV52_ANAL) + sum(abm_sims_50[[j+current]]$HPV52_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV52_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV52_ANAL) + sum(abm_sims_50[[j+current]]$HPV52_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV52_ANAL)))
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV52_ANAL)))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
        stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV52_ANAL[1:(365*1)])))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV52_ANAL[1:(365*5)])))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV52_ANAL[1:(365*10)])))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV52_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV52_ORAL) + sum(abm_sims_50[[j+current]]$HPV52_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV52_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV52_ORAL) + sum(abm_sims_50[[j+current]]$HPV52_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV52_ORAL)))
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV52_ORAL)))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
        stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV52_ORAL[1:(365*1)])))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV52_ORAL[1:(365*5)])))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV52_ORAL[1:(365*10)])))
      } else if (hpv_type[i]==58) {
        stat_prev_anal = c(stat_prev_anal, (sum(abm_sims_50[[j+current]]$HPV58_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_anal = c(stat_incidence_anal, (sum(abm_sims_50[[j+current]]$HPV58_ANAL) + sum(abm_sims_50[[j+current]]$HPV58_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV58_ANAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_anal = c(stat_new_anal, (sum(abm_sims_50[[j+current]]$HPV58_ANAL) + sum(abm_sims_50[[j+current]]$HPV58_ANAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV58_ANAL)))
        stat_total_anal = c(stat_total_anal, (sum(abm_sims_50[[j+current]]$HPV58_ANAL)))
        stat_prevent_anal = c(stat_prevent_anal, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ANAL)))
        stat_prevent_vax_anal = c(stat_prevent_vax_anal, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ANAL)))
        stat_prevent_sero_anal = c(stat_prevent_sero_anal, (sum(abm_sims_50[[j+current]]$SERO_PREVENT)))
        stat_prevent_cond_anal = c(stat_prevent_cond_anal, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ANAL)))
        stat_yr1_new_anal = c(stat_yr1_new_anal, (sum(long_sims_50[[j+current]]$HPV58_ANAL[1:(365*1)])))
        stat_yr5_new_anal = c(stat_yr5_new_anal, (sum(long_sims_50[[j+current]]$HPV58_ANAL[1:(365*5)])))
        stat_yr10_new_anal = c(stat_yr10_new_anal, (sum(long_sims_50[[j+current]]$HPV58_ANAL[1:(365*10)])))
        stat_prev_oral = c(stat_prev_oral, (sum(abm_sims_50[[j+current]]$HPV58_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_incidence_oral = c(stat_incidence_oral, (sum(abm_sims_50[[j+current]]$HPV58_ORAL) + sum(abm_sims_50[[j+current]]$HPV58_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV58_ORAL))/nrow(abm_sims_50[[j+current]]))
        stat_new_oral = c(stat_new_oral, (sum(abm_sims_50[[j+current]]$HPV58_ORAL) + sum(abm_sims_50[[j+current]]$HPV58_ORAL_CLEARED) - sum(abm_sims_50[[j+current]]$ORIGINAL_STATUS_HPV58_ORAL)))
        stat_total_oral = c(stat_total_oral, (sum(abm_sims_50[[j+current]]$HPV58_ORAL)))
        stat_prevent_oral = c(stat_prevent_oral, (sum(abm_sims_50[[j+current]]$OVERALL_PREVENT_ORAL)))
        stat_prevent_vax_oral = c(stat_prevent_vax_oral, (sum(abm_sims_50[[j+current]]$VAX_PREVENT_ORAL)))
        stat_prevent_cond_oral = c(stat_prevent_cond_oral, (sum(abm_sims_50[[j+current]]$COND_PREVENT_ORAL)))
        stat_yr1_new_oral = c(stat_yr1_new_oral, (sum(long_sims_50[[j+current]]$HPV58_ORAL[1:(365*1)])))
        stat_yr5_new_oral = c(stat_yr5_new_oral, (sum(long_sims_50[[j+current]]$HPV58_ORAL[1:(365*5)])))
        stat_yr10_new_oral = c(stat_yr10_new_oral, (sum(long_sims_50[[j+current]]$HPV58_ORAL[1:(365*10)])))
      }
    }
    
    stat_prev_anal = stat_prev_anal[-1]
    stat_incidence_anal = stat_incidence_anal[-1]
    stat_new_anal = stat_new_anal[-1]
    stat_total_anal = stat_total_anal[-1]
    stat_prevent_anal = stat_prevent_anal[-1]
    stat_prevent_vax_anal = stat_prevent_vax_anal[-1]
    stat_prevent_sero_anal = stat_prevent_sero_anal[-1]
    stat_prevent_cond_anal = stat_prevent_cond_anal[-1]
    stat_yr1_new_anal = stat_yr1_new_anal[-1]
    stat_yr5_new_anal = stat_yr5_new_anal[-1]
    stat_yr10_new_anal = stat_yr10_new_anal[-1]
    stat_prev_oral = stat_prev_oral[-1]
    stat_incidence_oral = stat_incidence_oral[-1]
    stat_new_oral = stat_new_oral[-1]
    stat_total_oral = stat_total_oral[-1]
    stat_prevent_oral = stat_prevent_oral[-1]
    stat_prevent_vax_oral = stat_prevent_vax_oral[-1]
    stat_prevent_cond_oral = stat_prevent_cond_oral[-1]
    stat_yr1_new_oral = stat_yr1_new_oral[-1]
    stat_yr5_new_oral = stat_yr5_new_oral[-1]
    stat_yr10_new_oral = stat_yr10_new_oral[-1]
    
    results_infections_anal = rbind(results_infections_anal, data.frame("HPV_type"=hpv_type[i],"Vax_level"=vax_level[j],"Prev_mean"=mean(stat_prev_anal),"Prev_SE"=(sd(stat_prev_anal)/sqrt(length(stat_prev_anal))),"Incidence_mean"=mean(stat_incidence_anal),"Incidence_SE"=(sd(stat_incidence_anal)/sqrt(length(stat_incidence_anal))),"New_mean"=mean(stat_new_anal),"New_SE"=(sd(stat_new_anal)/sqrt(length(stat_new_anal))),"Total_mean"=mean(stat_total_anal),"Total_SE"=(sd(stat_total_anal)/sqrt(length(stat_total_anal))),"Prevent_mean"=mean(stat_prevent_anal),"Prevent_SE"=(sd(stat_prevent_anal)/sqrt(length(stat_prevent_anal))),"Vax_prevent_mean"=mean(stat_prevent_vax_anal),"Vax_prevent_SE"=(sd(stat_prevent_vax_anal)/sqrt(length(stat_prevent_vax_anal))),"Sero_prevent_mean"=mean(stat_prevent_sero_anal),"Sero_prevent_SE"=(sd(stat_prevent_sero_anal)/sqrt(length(stat_prevent_sero_anal))),"Cond_prevent_mean"=mean(stat_prevent_cond_anal),"Cond_prevent_SE"=(sd(stat_prevent_cond_anal)/sqrt(length(stat_prevent_cond_anal))),stringsAsFactors=F))
    results_infections_oral = rbind(results_infections_oral, data.frame("HPV_type"=hpv_type[i],"Vax_level"=vax_level[j],"Prev_mean"=mean(stat_prev_oral),"Prev_SE"=(sd(stat_prev_oral)/sqrt(length(stat_prev_oral))),"Incidence_mean"=mean(stat_incidence_oral),"Incidence_SE"=(sd(stat_incidence_oral)/sqrt(length(stat_incidence_oral))),"New_mean"=mean(stat_new_oral),"New_SE"=(sd(stat_new_oral)/sqrt(length(stat_new_oral))),"Total_mean"=mean(stat_total_oral),"Total_SE"=(sd(stat_total_oral)/sqrt(length(stat_total_oral))),"Prevent_mean"=mean(stat_prevent_oral),"Prevent_SE"=(sd(stat_prevent_oral)/sqrt(length(stat_prevent_oral))),"Vax_prevent_mean"=mean(stat_prevent_vax_oral),"Vax_prevent_SE"=(sd(stat_prevent_vax_oral)/sqrt(length(stat_prevent_vax_oral))),"Cond_prevent_mean"=mean(stat_prevent_cond_oral),"Cond_prevent_SE"=(sd(stat_prevent_cond_oral)/sqrt(length(stat_prevent_cond_oral))),stringsAsFactors=F))
    results_longitudinal = rbind(results_longitudinal, data.frame("HPV_type"=hpv_type[i],"Vax_level"=vax_level[j],"Yr1_anal_new_mean"=mean(stat_yr1_new_anal),"Yr1_anal_new_SE"=(sd(stat_yr1_new_anal)/sqrt(length(stat_yr1_new_anal))),"Yr5_anal_new_mean"=mean(stat_yr5_new_anal),"Yr5_anal_new_SE"=(sd(stat_yr5_new_anal)/sqrt(length(stat_yr5_new_anal))),"Yr10_anal_new_mean"=mean(stat_yr10_new_anal),"Yr10_anal_new_SE"=(sd(stat_yr10_new_anal)/sqrt(length(stat_yr10_new_anal))),"Yr1_oral_new_mean"=mean(stat_yr1_new_oral),"Yr1_oral_new_SE"=(sd(stat_yr1_new_oral)/sqrt(length(stat_yr1_new_oral))),"Yr5_oral_new_mean"=mean(stat_yr5_new_oral),"Yr5_oral_new_SE"=(sd(stat_yr5_new_oral)/sqrt(length(stat_yr5_new_oral))),"Yr10_oral_new_mean"=mean(stat_yr10_new_oral),"Yr10_oral_new_SE"=(sd(stat_yr10_new_oral)/sqrt(length(stat_yr10_new_oral))),stringsAsFactors=F))
  }
}
rm(i,j,k,current,stat_prev_anal,stat_incidence_anal,stat_new_anal,stat_total_anal,stat_prevent_anal,stat_prevent_vax_anal,stat_prevent_sero_anal,stat_prevent_cond_anal,stat_yr1_new_anal,stat_yr5_new_anal,stat_yr10_new_anal,stat_prev_oral,stat_incidence_oral,stat_new_oral,stat_total_oral,stat_prevent_oral,stat_prevent_vax_oral,stat_prevent_cond_oral,stat_yr1_new_oral,stat_yr5_new_oral,stat_yr10_new_oral)
results_infections_anal = results_infections_anal[-1, ]
results_infections_oral = results_infections_oral[-1, ]
results_longitudinal = results_longitudinal[-1, ]

#overall starting prevalence (all simulations start out the same)
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV_ANAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV_ORAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HIV[1:nagents_start])/nagents_start

sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV6_ANAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV11_ANAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV16_ANAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV18_ANAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV31_ANAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV33_ANAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV45_ANAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV52_ANAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV58_ANAL[1:nagents_start])/nagents_start

sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV6_ORAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV11_ORAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV16_ORAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV18_ORAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV31_ORAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV33_ORAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV45_ORAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV52_ORAL[1:nagents_start])/nagents_start
sum(abm_sims_50[[1]]$ORIGINAL_STATUS_HPV58_ORAL[1:nagents_start])/nagents_start

#population at 10yr simulation end
mean(nagents_end); sd(nagents_end)

#quantify sex
sex = NA
partners = NA
for (i in 1:nsims) {
  sex = c(sex, abm_sims_50[[i*nparadigms]]$SEX_COUNT)
  partners = c(partners, abm_sims_50[[i*nparadigms]]$PARTNERS)
}
summary(partners, na.rm=T)
summary(sex, na.rm=T)
rm(i,sex,partners)

#% reduction in 10yr prevalence overall at each vaccination level compared to present day vaccination levels (13%)
mean((results_infections_anal$Total_mean[results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$Vax_level==25]) / results_infections_anal$Total_mean[results_infections_anal$Vax_level==13]) * 100
mean((results_infections_anal$Total_mean[results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$Vax_level==50]) / results_infections_anal$Total_mean[results_infections_anal$Vax_level==13]) * 100
mean((results_infections_anal$Total_mean[results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$Vax_level==13]) * 100
mean((results_infections_anal$Total_mean[results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$Vax_level==100]) / results_infections_anal$Total_mean[results_infections_anal$Vax_level==13]) * 100

mean((results_infections_oral$Total_mean[results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$Vax_level==25]) / results_infections_oral$Total_mean[results_infections_oral$Vax_level==13]) * 100
mean((results_infections_oral$Total_mean[results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$Vax_level==50]) / results_infections_oral$Total_mean[results_infections_oral$Vax_level==13]) * 100
mean((results_infections_oral$Total_mean[results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$Vax_level==13]) * 100
mean((results_infections_oral$Total_mean[results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$Vax_level==100]) / results_infections_oral$Total_mean[results_infections_oral$Vax_level==13]) * 100

#% reduction in 10yr prevalence by serotype at 80% vaccination (healthy people goal) compared to present day vaccination levels (13%)
mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==6 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==6 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==6 & results_infections_anal$Vax_level==13]) * 100
mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==11 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==11 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==11 & results_infections_anal$Vax_level==13]) * 100
mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==16 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==16 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==16 & results_infections_anal$Vax_level==13]) * 100
mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==18 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==18 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==18 & results_infections_anal$Vax_level==13]) * 100
mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==31 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==31 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==31 & results_infections_anal$Vax_level==13]) * 100
mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==33 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==33 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==33 & results_infections_anal$Vax_level==13]) * 100
mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==45 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==45 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==45 & results_infections_anal$Vax_level==13]) * 100
mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==52 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==52 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==52 & results_infections_anal$Vax_level==13]) * 100
mean((results_infections_anal$Total_mean[results_infections_anal$HPV_type==58 & results_infections_anal$Vax_level==13] - results_infections_anal$Total_mean[results_infections_anal$HPV_type==58 & results_infections_anal$Vax_level==80]) / results_infections_anal$Total_mean[results_infections_anal$HPV_type==58 & results_infections_anal$Vax_level==13]) * 100

mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==6 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==6 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==6 & results_infections_oral$Vax_level==13]) * 100
mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==11 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==11 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==11 & results_infections_oral$Vax_level==13]) * 100
mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==16 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==16 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==16 & results_infections_oral$Vax_level==13]) * 100
mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==18 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==18 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==18 & results_infections_oral$Vax_level==13]) * 100
mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==31 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==31 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==31 & results_infections_oral$Vax_level==13]) * 100
mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==33 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==33 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==33 & results_infections_oral$Vax_level==13]) * 100
mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==45 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==45 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==45 & results_infections_oral$Vax_level==13]) * 100
mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==52 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==52 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==52 & results_infections_oral$Vax_level==13]) * 100
mean((results_infections_oral$Total_mean[results_infections_oral$HPV_type==58 & results_infections_oral$Vax_level==13] - results_infections_oral$Total_mean[results_infections_oral$HPV_type==58 & results_infections_oral$Vax_level==80]) / results_infections_oral$Total_mean[results_infections_oral$HPV_type==58 & results_infections_oral$Vax_level==13]) * 100

#percent reduction in new infections at 1yr, 5yr, 10yr time points at each vaccination level compared to present day vaccination levels (13%)
summary(((results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==25])) / (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==50])) / (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==80])) / (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==100])) / (results_longitudinal$Yr1_anal_new_mean[results_longitudinal$Vax_level==13])) * 100

summary(((results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==25])) / (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==50])) / (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==80])) / (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==100])) / (results_longitudinal$Yr5_anal_new_mean[results_longitudinal$Vax_level==13])) * 100

summary(((results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==25])) / (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==50])) / (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==80])) / (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==100])) / (results_longitudinal$Yr10_anal_new_mean[results_longitudinal$Vax_level==13])) * 100

summary(((results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==25])) / (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==50])) / (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==80])) / (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==100])) / (results_longitudinal$Yr1_oral_new_mean[results_longitudinal$Vax_level==13])) * 100

summary(((results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==25])) / (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==50])) / (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==80])) / (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==100])) / (results_longitudinal$Yr5_oral_new_mean[results_longitudinal$Vax_level==13])) * 100

summary(((results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==25])) / (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==50])) / (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==80])) / (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13])) * 100
summary(((results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13]) - (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==100])) / (results_longitudinal$Yr10_oral_new_mean[results_longitudinal$Vax_level==13])) * 100

#% prevention attributed to each strategy, anal
mean(results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==0]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==0] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==0] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==0]) * 100
mean(results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==13]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==13] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==13] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==13]) * 100
mean(results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==25]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==25] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==25] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==25]) * 100
mean(results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==50]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==50] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==50] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==50]) * 100
mean(results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==80]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==80] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==80] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==80]) * 100
mean(results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==100]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==100] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==100] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==100]) * 100

mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==0]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==0] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==0] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==0]) * 100
mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==13]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==13] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==13] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==13]) * 100
mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==25]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==25] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==25] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==25]) * 100
mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==50]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==50] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==50] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==50]) * 100
mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==80]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==80] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==80] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==80]) * 100
mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==100]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==100] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==100] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==100]) * 100

mean(results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==0]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==0] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==0] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==0]) * 100
mean(results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==13]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==13] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==13] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==13]) * 100
mean(results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==25]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==25] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==25] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==25]) * 100
mean(results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==50]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==50] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==50] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==50]) * 100
mean(results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==80]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==80] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==80] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==80]) * 100
mean(results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==100]) / mean(results_infections_anal$Cond_prevent_mean[results_infections_anal$Vax_level==100] + results_infections_anal$Sero_prevent_mean[results_infections_anal$Vax_level==100] + results_infections_anal$Vax_prevent_mean[results_infections_anal$Vax_level==100]) * 100

#% prevention attributed to each strategy, oral
mean(results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==0]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==0] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==0]) * 100
mean(results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==13]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==13] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==13]) * 100
mean(results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==25]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==25] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==25]) * 100
mean(results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==50]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==50] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==50]) * 100
mean(results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==80]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==80] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==80]) * 100
mean(results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==100]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==100] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==100]) * 100

mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==0]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==0] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==0]) * 100
mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==13]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==13] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==13]) * 100
mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==25]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==25] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==25]) * 100
mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==50]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==50] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==50]) * 100
mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==80]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==80] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==80]) * 100
mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==100]) / mean(results_infections_oral$Cond_prevent_mean[results_infections_oral$Vax_level==100] + results_infections_oral$Vax_prevent_mean[results_infections_oral$Vax_level==100]) * 100
