# script for developing the elimination proposal
# Extension of Watkins paper to:
# - apply to the UK - each node is a UK local authority
# - * include enhanced enterovirus surveillance (ENT)
# - * include ENV within specifc locations
# - * include variable importation rate as specified in previous analysis
# - include clustering of infection by having region and unit design prevalence

# I'm giong to assume monthly time-window as this should be almost equivalent 
# to the duration of shedding
# weekly doesn't make sense because if someone sheds in week 1 they may still be detected in week 2 and
# the simple methods doesn't really account for this.

#library(rjags)
library(prevalence)  # betaPERT distribution

rm(list = ls())

writewd <- setwd("~/")

# load in LA data - for pop size and Pr(importation)
datLA <- read.csv("summary_polio_risk_forR_16Oct18.csv")
head(datLA)
numLA <- dim(datLA)[1]  # 348 LA...

# STEP 1: go through each of detection sensitivites....

Iter <- 1000  # 1000 iterations of each parameter

nR <- datLA$PopSize
nR[is.na(nR)] <- 80000   # made up
PrP <- nR/sum(nR)           # proportion in each node

P_age <- rep(1,Iter)             # Proportion within specific age category (note: there are no age categories here)

Dprev_region <- rep(1/200,Iter)   #   should detect 1 region in 200 (regional design prevalence)                              # Design prevalence
#Dprev_unit <- rep(1/100000,Iter)    # should detect 1 person in 100,000    (AFP threshold...but per year)
Dprev_unit <- rep(1/(100000/12),Iter)   # 18.9% aged 0-15 and this is annual. for the age bit just assume no people above 15 get AFP?
#   3.5e-7 - 20 infections in UK pop
mean(datLA$NUMBER0)  
# we want to find a case in a LA if it's there...

# (assume everyone infected sheds poliovirus)

# (A). AFP surveillance

tmp <- betaExpert(best=1/190,lower=1/250,upper=1/150,p=0.95)
P_AFPclinical_wild <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta)   # CLINICAL case (<1%) ...  # hist(Naf)
signif(quantile(P_AFPclinical_wild,probs=c(0.025,0.5,0.975)),3)[c(2,1,3)] 

tmp <- betaExpert(best=1/1886,lower=1/2200,upper=1/1000,p=0.95)
P_AFPclinical_vdpv <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta)   # CLINICAL case (<1%) ...  # hist(Naf)
signif(quantile(P_AFPclinical_vdpv,probs=c(0.025,0.5,0.975)),3)[c(2,1,3)] 

tmp <- betaExpert(best=0.9,lower=0.6,upper=0.99,p=0.95)
P_not_AFP <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta)   # Notification rate (this is still "uncertain") ...  # hist(Nno)

tmp <- betaExpert(best=0.8,lower=0.50,upper=0.95,p=0.95)
P_stool <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta)   # Stool/CSF collection 

tmp <- betaExpert(best=0.97,lower=0.95,upper=0.9999,p=0.95)
P_AFPtest <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta)   # Test sensitivity # hist(Nte)

# (B). Enhanced Enterovirus surveillance

tmp <- betaExpert(best=(1/190)*0.01,lower=(1/250)*0.01,upper=(1/150)*0.01,p=0.95)  # wild
#tmp <- betaExpert(best=(1/190)*0.01,lower=(1/250)*0.01,upper=(1/150)*0.01,p=0.95)  # VDPV
P_MENclinical <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta)   # CLINICAL case (1-5%) ...  # hist(Naf)
signif(quantile(P_MENclinical,probs=c(0.5,0.025,0.975)),3)

tmp <- betaExpert(best=0.9,lower=0.6,upper=0.99,p=0.95)
P_not_MEN <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta)   # Notification rate (reporting to hospital) - same as AFP
signif(quantile(P_AFPclinical_vdpv,probs=c(0.025,0.5,0.975)),3)[c(2,1,3)] 

tmp <- betaExpert(best=0.8,lower=0.50,upper=0.95,p=0.95)
P_CSF <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta)   # Stool/CSF collection  # hist(Nst) - as good as AFP

tmp <- betaExpert(best=0.97,lower=0.95,upper=0.9999,p=0.95)
P_MENtest <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta)   # Test sensitivity # hist(Nte)

# (C). ENV surveillance

P_ENVcatchB <- datLA$BecktonCatch*0.8    # beckton *0.8 to account for not sampling for some reason # Pr(of being within catchment) 
P_ENVcatchRB <- (datLA$BirMManBradCatch+datLA$BecktonCatch)*0.8    # *0.8 to account for not sampling for some reason # Pr(of being within catchment) 
P_ENVcatchR <- (datLA$BirMManBradCatch)*0.8    # *0.8 to account for not sampling for some reason # Pr(of being within catchment) 
P_ENVcatchE <- datLA$AllHighRisk*0.8 

#tmp <- betaExpert(best=24/30,lower=0.5,upper=0.9,p=0.95) # monthly
tmp <- betaExpert(best=0.99,lower=0.9,upper=0.999,p=0.95) # fortnightly
P_ENVsamp <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta) # shedding for 16 of 30.4 days and sampling for 7 of 30.4...     # assume monthly sampling - if we went to weekly we would multiply 'unit' by 4
# could/should out some variation this
hist(P_ENVsamp)

tmp <- betaExpert(best=0.97,lower=0.95,upper=0.9999,p=0.95)
P_ENVtest <- rbeta(n=Iter,shape1=tmp$alpha,shape2=tmp$beta)   # ES test sens (this might be low because of CSF) # hist(Nst)

# STEP 2: based on above, calculate Unit Sensitivity (USe) and Surveillance Sensitivity (SSe) for 1 month.

# *** Accounting for DIFFERENTIAL RISK
tmp <- datLA$risk_all+1e-5      #+1e-10 #(plus a bit so min risk >0.00) [1:100]          #  a log scale and not min=1
RRi <- tmp/min(tmp)             # make lowest risk = 1.00
#RRi[RRi>10] <- 10 #RRi[RRi>67151] <- 67151
RRi <- RRi/min(RRi)
 
PrPi <- PrP #tmp/sum(tmp)  #PrP#[1:100]/sum(PrP[1:100])  sum(PrPi)

sum(datLA$PopSize[datLA$BecktonCatch==1])/sum(datLA$PopSize,na.rm = T)

# hist(RRi)                    # so min is 1.00 but RR for some is now very large...
# hist(PrPi)
# sum(PrPi)
source("elim_functions_11Oct18.r")   # used to solve equations below
data <- data.frame(RR=RRi,PrP=PrPi)
answers <-  ARfunct2(data)
AR <- answers  #  Adjusted risk...
# check sum(PrPi*AR) # = 1.00 hurray!! hist(AR)
# plot(RRi,AR)
hist(AR)
 
#  ** AFP Surveillance sensitivity for each region

# probabilty of detection might be easier to articulate first...
# AFP
signif(quantile(P_AFPclinical_wild*P_not_AFP*P_stool*P_AFPtest,probs = c(0.5,0.025,0.975)),3) # 
signif(quantile(P_AFPclinical_vdpv*P_not_AFP*P_stool*P_AFPtest,probs = c(0.5,0.025,0.975)),3)
# MEN
signif(quantile(P_MENclinical*P_not_MEN*P_CSF*P_MENtest,probs = c(0.5,0.025,0.975)),3)
# ENV
signif(quantile(0.8*P_ENVsamp*P_ENVtest,probs = c(0.025,0.5,0.975)),3)

# combined AFP and ENV - wild. Non-overlapping so just add together.
signif(quantile((P_AFPclinical_wild*P_not_AFP*P_stool*P_AFPtest)+(P_MENclinical*P_not_MEN*P_CSF*P_MENtest),probs = c(0.5,0.025,0.975)),3) # 
# VDPV
signif(quantile((P_AFPclinical_vdpv*P_not_AFP*P_stool*P_AFPtest)+(P_MENclinical*P_not_MEN*P_CSF*P_MENtest),probs = c(0.5,0.025,0.975)),3) # 

summary(P_AFPtest)
summary(Dprev_unit)

# then we combine surveillance with the effective probabilty of infection (EPI);
# EPI = AR*Ph (if only thinking about Ph - 1/100,000 in UK % 12? for per month?)
# EPI = AR*Ph*Pr  (detecting 1 infected region with at least 1/100000 per region)
# *** the loop below needs to be changed according to whether its WILD or VDPVs...***

SeR_AFPi <- Pneg_AFPi <- matrix(0,nrow=Iter,ncol=numLA)
for(i in 1:numLA){
  SeR_AFPi[,i] <- 1-(1-(Dprev_unit*P_AFPclinical_wild*P_not_AFP*P_stool*P_AFPtest))^nR[i]  # AFP
  Pneg_AFPi[,i] <- 1-(SeR_AFPi[,i]*Dprev_region*AR[i])
}
hist(SeR_AFPi[1,])  # 1 iteration
hist(Pneg_AFPi[1,]) 

SSCSe_AFP <- rep(NA,Iter)
for(i in 1:Iter){
  SSCSe_AFP[i] <- 1-prod(Pneg_AFPi[i,])
}
hist(SSCSe_AFP)
signif(quantile(SSCSe_AFP,probs = c(0.025,0.5,0.975)),3)

#  ** ENT Surveillance sensitivity for each region
SeR_ENTi <- Pneg_ENTi <- matrix(0,nrow=Iter,ncol=numLA)
for(i in 1:numLA){
  SeR_ENTi[,i] <- 1-(1-(Dprev_unit*P_MENclinical*P_not_MEN*P_CSF*P_MENtest))^nR[i]  # AFP
  Pneg_ENTi[,i] <- 1-(SeR_ENTi[,i]*Dprev_region*AR[i])
}
hist(SeR_ENTi[1,])  # 1 iteration
hist(Pneg_ENTi[1,]) 

SSCSe_ENT <- rep(NA,Iter)
for(i in 1:Iter){
  SSCSe_ENT[i] <- 1-prod(Pneg_ENTi[i,])
}
hist(SSCSe_ENT)
signif(quantile(SSCSe_ENT,probs = c(0.025,0.5,0.975)),3)

# ENV surveillance - assuming different options for sampling sites
# *** updated ***
# ENV_R - risky sites; Bradford, Manchester, Birmingham (3 sites)
# ENV_B - becton only (1 site)
# ENV_T - Thames (4 sites)
# ENV_E - exhaustive (all sites in table - 10 sites)
SeR_ENV_Ri <- Pneg_ENV_Ri <- matrix(0,nrow=Iter,ncol=numLA)  # risky
SeR_ENV_Bi <- Pneg_ENV_Bi <- matrix(0,nrow=Iter,ncol=numLA)
SeR_ENV_RBi <- Pneg_ENV_RBi <- matrix(0,nrow=Iter,ncol=numLA) # risky and beckton
SeR_ENV_Ei <- Pneg_ENV_Ei <- matrix(0,nrow=Iter,ncol=numLA)
for(i in 1:numLA){
  # risky
  SeR_ENV_Ri[,i] <- 1-(1-(Dprev_unit*P_ENVcatchR[i]*P_ENVsamp*P_ENVtest))^nR[i]  # ENV
  Pneg_ENV_Ri[,i] <- 1-(SeR_ENV_Ri[,i]*Dprev_region*AR[i])
  # beckton
  SeR_ENV_Bi[,i] <- 1-(1-(Dprev_unit*P_ENVcatchB[i]*P_ENVsamp*P_ENVtest))^nR[i]  # ENV
  Pneg_ENV_Bi[,i] <- 1-(SeR_ENV_Bi[,i]*Dprev_region*AR[i])
  # risk and beckton
  SeR_ENV_RBi[,i] <- 1-(1-(Dprev_unit*P_ENVcatchRB[i]*P_ENVsamp*P_ENVtest))^nR[i]  # AFP
  Pneg_ENV_RBi[,i] <- 1-(SeR_ENV_RBi[,i]*Dprev_region*AR[i])
  # all high risk
  SeR_ENV_Ei[,i] <- 1-(1-(Dprev_unit*P_ENVcatchE[i]*P_ENVsamp*P_ENVtest))^nR[i]  # AFP
  Pneg_ENV_Ei[,i] <- 1-(SeR_ENV_Ei[,i]*Dprev_region*AR[i])
}

SSCSe_ENV_R <- SSCSe_ENV_B <- SSCSe_ENV_RB <- SSCSe_ENV_E <- rep(NA,Iter)
for(i in 1:Iter){
  SSCSe_ENV_R[i] <- 1-prod(Pneg_ENV_Ri[i,])
  SSCSe_ENV_B[i] <- 1-prod(Pneg_ENV_Bi[i,])
  SSCSe_ENV_RB[i] <- 1-prod(Pneg_ENV_RBi[i,])
  SSCSe_ENV_E[i] <- 1-prod(Pneg_ENV_Ei[i,])
}
hist(SSCSe_ENV_R)

signif(quantile(SSCSe_ENV_R,probs = c(0.5,0.025,0.975)),3)
signif(quantile(SSCSe_ENV_B,probs = c(0.5,0.025,0.975)),3)
signif(quantile(SSCSe_ENV_RB,probs = c(0.5,0.025,0.975)),3)
signif(quantile(SSCSe_ENV_E,probs = c(0.5,0.025,0.975)),3)

# combining the surveillance together....
SSe_AFP <- 1 - ((1-SSCSe_AFP))
signif(quantile(SSe_AFP,probs=c(0.025,0.5,0.975)),3)[c(2,1,3)] # - terrible on its own!

SSe_ENT <- 1 - ((1-SSCSe_ENT))
signif(quantile(SSe_ENT,probs=c(0.5,0.025,0.975)),2) # 0.01220 0.00581 0.02090 (wild)

SSe_AE <- 1 - ((1-SSCSe_AFP)*(1-SSCSe_ENT))
signif(quantile(SSe_AE,probs=c(0.5,0.025,0.975)),2) # 0.01220 0.00581 0.02090 (wild)

SSe_R <- 1 - ((1-SSCSe_AFP)*(1-SSCSe_ENT)*(1-SSCSe_ENV_R))  # ENV_R = risky
hist(SSe_R)   # 
signif(quantile(SSe_R,probs=c(0.025,0.5,0.975)),3)[c(2,1,3)]  # 0.328 0.324 0.334 

# both? 
SSe_RB <- 1 - ((1-SSCSe_AFP)*(1-SSCSe_ENT)*(1-SSCSe_ENV_R)*(1-SSCSe_ENV_B))  # risky and Beckton
hist(SSe_RB)
signif(quantile(SSe_RB,probs=c(0.5,0.025,0.975)),3)  #  0.202 0.197 0.209

SSe_E <- 1 - ((1-SSCSe_AFP)*(1-SSCSe_ENT)*(1-SSCSe_ENV_E))  # ENV_E = exhaustive (all sites in table - 10 sites)
hist(SSe_E)   # 
signif(quantile(SSe_E,probs=c(0.025,0.5,0.975)),3)[c(2,1,3)]  # 0.328 0.324 0.334 



# STEP 3: include time, and calculate the Pr(Free) - AFP alone

# when was the last case of polio?
# 1993 - reported in Salisbury (1997)
# 1993 to now 1993-2018 = 25 years
# when was ENT introduced? 1997-1993 = 4 years in
# assuming ENV is introduced in final year...ie. 2019...year 25

max_time <- 25*12       # 25 years of monthly surveillance
time <- 0:max_time

# should be able to import the number of movements in and out of UK
# Have seperate movements for VDPV risk and wild risk.
pops <- read.csv("migration_polio_data_12Dec19.csv")
head(pops)

# AFP and ENT surveillance - default
Post_freeA <- matrix(0,nrow=length(time),ncol=Iter)   # A is with no introduction risk
Post_freeB <- matrix(0,nrow=length(time),ncol=Iter)   # B is some introduction risk
Prior_infeA <- matrix(0,nrow=length(time),ncol=Iter)
Prior_infeB <- matrix(0,nrow=length(time),ncol=Iter)

# + ENV surveillance
Post_freeA_ENV <- matrix(0,nrow=length(time),ncol=Iter)  
Post_freeB_ENV <- matrix(0,nrow=length(time),ncol=Iter)
Prior_infeA_ENV <- matrix(0,nrow=length(time),ncol=Iter)
Prior_infeB_ENV <- matrix(0,nrow=length(time),ncol=Iter)

# Pr(introduction) - note this should be updated to be region-specific
# datLA$risk_all is on a log scale and strictly speaking is the combined risk of importation and circulation
# default : best=1/(1000),lower=1/5000,upper=1/500
# should change it really...but it's small mostly
PP_intro <-   betaExpert(best=1/(1000),lower=1/5000,upper=1/500,p=0.95)  # default
#PP_intro <-   betaExpert(best=1/(10000),lower=1/50000,upper=1/5000,p=0.95)  # low
#PP_intro <-   betaExpert(best=1/(100),lower=1/1000,upper=1/10,p=0.95)  # high risk
Pr_intro <- rbeta(n=Iter,shape1=PP_intro$alpha,shape2=PP_intro$beta)

hist(Pr_intro) 

dim(Prior_infeA)  

# initially, we are unsure about freedom - prior is 0.5
# let's say last case was 6 mths ago
# if as a country we assume a 50% chance, we need to split it into the different LAs - otherwise risk is high!

survi <- "wild"  #"vdpv"   #wild"
if(survi=="wild"){
  Prior_infeA[1,] <- Prior_infeB[1,] <- Post_freeA[1,] <- Post_freeB[1,] <- rbeta(Iter,2,2)
  Prior_infeA_ENV[1,] <- Prior_infeB_ENV[1,] <- Post_freeA_ENV[1,] <- Post_freeB_ENV[1,] <- rbeta(Iter,2,2)
}
if(survi!="wild"){
  # vdpv never existed here!
  Prior_infeA[1,] <- Prior_infeB[1,] <- Post_freeA[1,] <- Post_freeB[1,] <- rbeta(Iter,2,2)
  Prior_infeA_ENV[1,] <- Prior_infeB_ENV[1,] <- Post_freeA_ENV[1,] <- Post_freeB_ENV[1,] <- rbeta(Iter,2,2)
}
for(t in 2:(max_time)){
  timei <- floor(t/12)+1993  # time in years
  row <- which(pops$Year==timei)
  intro_wild <- 1/12*( ((pops$Visit.PAK[row]/10)+pops$Visitor.PAK[row])*(pops$wild.PAK[row]/138500000) + ((pops$Visit.IND[row]/10)+pops$Visitor.IND[row])*(pops$wild.IND[row]/1053000000) + ((pops$Visit.AFR[row]/10)+pops$Visitor.AFR[row])*(pops$wild.AFR[row]/811000000))
  intro_vdpv <- 1/12*( ((pops$Visit.PAK[row]/10)+pops$Visitor.PAK[row])*(pops$VDPV.PAK[row]/138500000) + ((pops$Visit.IND[row]/10)+pops$Visitor.IND[row])*(pops$VDPV.IND[row]/1053000000) + ((pops$Visit.AFR[row]/10)+pops$Visitor.AFR[row])*(pops$VDPV.AFR[row]/811000000))
  
  if(survi=="wild"){
    intro <- intro_wild
  }
  if(survi!="wild"){
    intro <- intro_vdpv
  }
  # *** CLINICAL SURVEILLANCE
  # loop through each month in t
  #                    PRIOR infected        PRIOR infected           DATA
  if(t > 12*4){ # + ENT
    Post_freeA[t,] <- (1-(1-Post_freeA[t-1,]))/(1-(1-Post_freeA[t-1,])*SSe_AE[sample(1:Iter,Iter)]) # 
    # loop through each month - including introduction risk  (scenario B)
    #                    PRIOR infected        PRIOR infected           DATA
    Post_freeB[t,] <- (1-(Prior_infeB[t-1,]))/(1-(Prior_infeB[t-1,])*SSe_AE[sample(1:Iter,Iter)])
    Prior_infeB[t,] <- (1-Post_freeB[t,])+Pr_intro[sample(1:Iter,Iter)]*intro_wild-((1-Post_freeB[t,])*Pr_intro[sample(1:Iter,Iter)]*intro)
  }else{  # AFP only
    Post_freeA[t,] <- (1-(1-Post_freeA[t-1,]))/(1-(1-Post_freeA[t-1,])*SSe_AFP[sample(1:Iter,Iter)]) # 
    # loop through each month - including introduction risk  (scenario B)
    #                    PRIOR infected        PRIOR infected           DATA
    Post_freeB[t,] <- (1-(Prior_infeB[t-1,]))/(1-(Prior_infeB[t-1,])*SSe_AFP[sample(1:Iter,Iter)])
    Prior_infeB[t,] <- (1-Post_freeB[t,])+Pr_intro[sample(1:Iter,Iter)]*intro_wild-((1-Post_freeB[t,])*Pr_intro[sample(1:Iter,Iter)]*intro)
  }
  
  # *** + ENV SURVEILLANCE
  # loop through each month - ignoring introduction risk  (scenario A)
  #                    PRIOR infected        PRIOR infected           DATA
  if(t > 12*4 & t < 12*25){ # AFP & ENT
    Post_freeA_ENV[t,] <- (1-(1-Post_freeA_ENV[t-1,]))/(1-(1-Post_freeA_ENV[t-1,])*SSe_AE[sample(1:Iter,Iter)]) # 
    # loop through each month - including introduction risk  (scenario B)
    #                    PRIOR infected        PRIOR infected           DATA
    Post_freeB_ENV[t,] <- (1-(Prior_infeB_ENV[t-1,]))/(1-(Prior_infeB_ENV[t-1,])*SSe_AE[sample(1:Iter,Iter)])
    Prior_infeB_ENV[t,] <- (1-Post_freeB_ENV[t,])+Pr_intro[sample(1:Iter,Iter)]*intro_wild-((1-Post_freeB_ENV[t,])*Pr_intro[sample(1:Iter,Iter)]*intro)
  }
  if(t > 12*4 & t >= 12*25){ #+ ENV
    Post_freeA_ENV[t,] <- (1-(1-Post_freeA_ENV[t-1,]))/(1-(1-Post_freeA_ENV[t-1,])*SSe_R[sample(1:Iter,Iter)]) # 
    # loop through each month - including introduction risk  (scenario B)
    #                    PRIOR infected        PRIOR infected           DATA
    Post_freeB_ENV[t,] <- (1-(Prior_infeB_ENV[t-1,]))/(1-(Prior_infeB_ENV[t-1,])*SSe_R[sample(1:Iter,Iter)])
    Prior_infeB_ENV[t,] <- (1-Post_freeB_ENV[t,])+Pr_intro[sample(1:Iter,Iter)]*intro_wild-((1-Post_freeB_ENV[t,])*Pr_intro[sample(1:Iter,Iter)]*intro)
  }
  if(t <= 12*4 ){ # AFP only
    Post_freeA_ENV[t,] <- (1-(1-Post_freeA_ENV[t-1,]))/(1-(1-Post_freeA_ENV[t-1,])*SSe_AFP[sample(1:Iter,Iter)]) # 
    # loop through each month - including introduction risk  (scenario B)
    #                    PRIOR infected        PRIOR infected           DATA
    Post_freeB_ENV[t,] <- (1-(Prior_infeB_ENV[t-1,]))/(1-(Prior_infeB_ENV[t-1,])*SSe_AFP[sample(1:Iter,Iter)])
    Prior_infeB_ENV[t,] <- (1-Post_freeB_ENV[t,])+Pr_intro[sample(1:Iter,Iter)]*intro_wild-((1-Post_freeB_ENV[t,])*Pr_intro[sample(1:Iter,Iter)]*intro)
  }
}

quantile(Post_freeB_ENV[t,],probs=c(0.025,0.5,0.975))

# plots
fig_folder <- "~/"

#free_medsims <- apply(Post_freeB,1,quantile,probs=c(0.025,0.5,0.975))
#free_lowsims <- apply(Post_freeB,1,quantile,probs=c(0.025,0.5,0.975))
#free_highsims <- apply(Post_freeB,1,quantile,probs=c(0.025,0.5,0.975))

free_wildsims <- apply(Post_freeB,1,quantile,probs=c(0.025,0.5,0.975))[,1:300]
free_vdpvsims <- apply(Post_freeB,1,quantile,probs=c(0.025,0.5,0.975))[,1:300]

pdf(paste0(fig_folder,"Poliofree_AFP_ENT_type_3Feb20.pdf"),height=6,width=8)
#par(mfrow=c(2,1),mai=c(0.5,0.8,0.5,0.5))
plot(time[1:300]/12,free_wildsims[2,],pch=19,col="orange4",
     ylab="Pr(poliovirus free)",
     xaxt='n',xlab="",
     ylim=c(0,1),xlim=c(0,25))
points(time[1:300]/12,free_wildsims[1,],pch=19,col=rgb(139/255,90/255,0,0.1))
points(time[1:300]/12,free_wildsims[3,],pch=19,col=rgb(139/255,90/255,0,0.1))
# vdpv
points(time[1:300]/12,free_vdpvsims[2,],pch=19,col=rgb(102/255,205/255,0,1))
points(time[1:300]/12,free_vdpvsims[3,],pch=19,col=rgb(102/255,205/255,0,0.1))
points(time[1:300]/12,free_vdpvsims[1,],pch=19,col=rgb(102/255,205/255,0,0.1))
points(time[1:300]/12,free_wildsims[2,],pch=19,col="orange4")

legend('bottomright',legend=c('wild','vdpv'),pch=c(19,19),
       col=c(rgb(139/255,90/255,0,1),rgb(102/255,205/255,0,1)))

arrows(x0=4, y0=0.15, x1=4, y1=0.3, lwd=2 )
text(x=c(2,2),y=c(0.05,0.1)-0.01,labels = c("ENT surveillance","Introduction of"),pos=4)
axis(side=1,at=seq(0,(max_time+1)/12,1),labels=seq(1993,1993+25,1))
title(main="AFP & ENT surveillance",xlab="Year")
abline(h=0.95,lty=2)
dev.off()


pdf(paste0(fig_folder,"Poliofree_AFP_MEN_12Dec19.pdf"),height=6,width=8)
#par(mfrow=c(2,1),mai=c(0.5,0.8,0.5,0.5))
plot(time/12,free_medsims[2,],pch=19,col="orange4",
     ylab="Pr(polio free)",
     xaxt='n',xlab="",
     ylim=c(0,1),xlim=c(0,25))
points(time/12,free_medsims[1,],pch=19,col=rgb(139/255,90/255,0,0.1))
points(time/12,free_medsims[3,],pch=19,col=rgb(139/255,90/255,0,0.1))
# low
#points(time/12,free_lowsims[2,],pch=19,col=rgb(205/255,51/255,51/255,1)) # 205 51 51
#points(time/12,free_lowsims[1,],pch=19,col=rgb(205/255,51/255,51/255,0.1))
# high
#points(time/12,free_highsims[2,],pch=19,col=rgb(110/255,139/255,61/255,1))  # 110 139 61
#points(time/12,free_highsims[1,],pch=19,col=rgb(110/255,139/255,61/255,0.1))
arrows(x0=4, y0=0.15, x1=4, y1=0.3, lwd=2 )
text(x=c(2,2),y=c(0.05,0.1)-0.01,labels = c("MEN surveillance","Introduction of"),pos=4)
axis(side=1,at=seq(0,(max_time+1)/12,1),labels=seq(1993,1993+25,1))
title(main="AFP & MEN surveillance",xlab="Year")
abline(h=0.95,lty=2)
dev.off()

oo <- which(free_medsims[1,]>0.95)
(218/12)+1993
free_medsims[,75]

tmp2 <- data.frame(t(tmp))
tmp2$time <- (time/12)+1993
write.csv(tmp2,file="ddump.csv",row.names = F)

# plus ENV
# main thing here should be impact of importation I think..?

tmp <- apply(Post_freeB_ENV,1,quantile,probs=c(0.025,0.5,0.975))
pdf(paste0(fig_folder,"Elimination_AFP_ENT_ENV_import_12Dec18.pdf"),height=6*2,width=8)
par(mfrow=c(3,1),mai=c(0.5,0.8,0.5,0.5))
plot(time/12,tmp[2,],pch=4,xlab="Year of Surveillance",ylab="Pr(elimination)",
     xaxt='n',
     ylim=c(0,1),xlim=c(0,max(time)/12))
polygon(c(time/12,rev(time/12)),c(tmp[3,],rev(tmp[1,])),col="chartreuse4",border="chartreuse4")
points(time/12,tmp[2,],pch=4)
arrows(x0=4, y0=0.1, x1=4, y1=0.3, lwd=2 )
text(x=c(2,2),y=c(0.05,0.09)-0.01,labels = c("ENT surveillance","Introduction of"),pos=4)
arrows(x0=25, y0=0.7, x1=25, y1=0.8, lwd=2 )
text(x=c(21,21),y=c(0.65,0.69,0.72)-0.01,labels = c("ENV surveillance","introduction of","Planned"),pos=4)
axis(side=1,at=seq(0,(max_time+1)/12,1),labels=seq(1993,1993+27,1))
abline(h=0.95,lty=2)
abline(h=0.99,lty=2)
title(main="AFP, ENT & ENV Surveillance (High Importation Assumption)")
dev.off()

tmp2 <- data.frame(t(tmp))
tmp2$time <- (time/12)+1993
write.csv(tmp2,file="ddump.csv",row.names = F)


pdf(paste0(fig_folder,"Elimination_est_12Dec18.pdf"),height=5,width=5)
par(mfrow=c(2,1),mai=c(0.5,0.8,0.5,0.5))
tmp <- apply(Post_freeA,1,quantile,probs=c(0.025,0.5,0.975))
plot(time/12,tmp[2,],pch=4,ylab="Pr(Free from poliovirus)",
     ylim=c(0,1),xlim=c(0,25))
polygon(c(time/12,rev(time/12)),c(tmp[3,],rev(tmp[1,])),col="skyblue",border="skyblue")
points(time/12,tmp[2,],pch=4)
title(main="AFP & ENT surveillance")
abline(h=0.95,lty=2)
# plus ENV
par(mai=c(0.8,0.8,0.5,0.5))
tmp <- apply(Post_freeA_ENV,1,quantile,probs=c(0.025,0.5,0.975))
plot(time/12,tmp[2,],pch=4,xlab="Time (years)",ylab="Pr(Free from poliovirus)",
     ylim=c(0,1),xlim=c(0,25))
polygon(c(time/12,rev(time/12)),c(tmp[3,],rev(tmp[1,])),col="chartreuse4",border="chartreuse4")
points(time/12,tmp[2,],pch=4)
abline(h=0.95,lty=2)
title(main="AFP, ENT & ENV Surveillance")
dev.off()





# end
