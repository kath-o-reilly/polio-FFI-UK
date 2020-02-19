# Risk framework for UK ES sampling
# --- Using movement data for UK residents and overseas visitors
# --- Using census data for location of foreign born nationals (by country)
# --- Using COVER data for ENGLAND and WALES (not all of the PENTA data is available at LA level - some approximations required)
# --- possible issues - fitting geographical locations for different data sets

#install.packages(pkgs=c("maptools","rgdal","tmap","sf","gridExtra","RColorBrewer"))

library(devtools)
install_github("r-spatial/sf")  # see here for additional things to install (takes a while) https://github.com/r-spatial/sf/

library(sf)
library(maptools)
library(rgdal)
library(tmap)   # nice new package for plotting maps? see above to get it working...
library(gridExtra)
library(RColorBrewer)

rm(list = ls(all = TRUE))  # clear all data

fig_folder <- "~/"
# load in polygons
shp <- readShapePoly("shp_county/County_UnitaryAuthority.shp")
shp0 <- readShapePoly("shp_unitary/Local_UnitaryAuthority.shp")

# Is there just one shp file that corresponds to the data?
# migration data are of LOCAL AUTHORITIES (England = 326, Wales = 22)

names(shp)
dim(shp)  # n = 144. Includes England, Wales, Scotland (kinda counties and London as 1 unit)
dim(shp0) # n = 380. Includes England, Wales, Scotland (LA?)
table(shp$NAME)  # mostly counties
table(shp$CODE)

# quick check that we can plot
pdf(paste0(fig_folder,"test.pdf"),height=6,width=4)
tm_shape(shp) + tm_fill(col="AREA") 
dev.off()

# England data - detail
names(shp0)
shp0$abbrev <- substring(shp0$CODE,1,1)
table(shp0$abbrev)
shp0EW <- shp0[shp0$abbrev!="S",]   # remove Scotland
dim(shp0EW)   # n = 348, which corresponds with the census data
head(shp0EW@data)
 
# STEP1: load in data of migrants and vaccination

# read in csv file of migrants
resid <- read.csv("Data/nomis_2011_long_term_residents_clean.csv",header=T,sep=",")
head(resid)
short <- read.csv("Data/nomis_2011_short_term_residents_clean.csv",header=T,sep=",")
head(short)
pent <- read.csv("Data/cover_20162_clean_15Aug18.csv",header=T,sep=",")  # this requires mutiple years
head(pent)

#plot(pent$penta2016,pent$pentaAvg,pch=19)

# STEP2: add migrant and cover data to the shapefile

# long-term residents
oo <- match(shp0EW$CODE,resid$Geography)  # length(oo)==dim(shp0EW)[1] yes!
shp0EW$PakLTR <- resid$Pakistan[oo]
shp0EW$PakLTR[is.na(shp0EW$PakLTR)] <- 0
shp0EW$NieLTR <- resid$Nigeria[oo]
shp0EW$NieLTR[is.na(shp0EW$NieLTR)] <- 0
shp0EW$WestLTR <- resid$Western.Africa[oo]
shp0EW$MeaLTR <- resid$Middle.East[oo]*0.0616  # accountin for Syria being v small
shp0EW$MeaLTR[is.na(shp0EW$MeaLTR)] <- 0
shp0EW$SasLTR <- resid$Southern.Asia[oo]*0.0189  # accounting for Afg being small
shp0EW$SasLTR[is.na(shp0EW$SasLTR)] <- 0

# short-term residents
oo <- match(shp0EW$CODE,short$Geography)  # length(oo)==dim(shp0EW)[1] yes!
shp0EW$PakSTR <- short$Pakistan[oo]
shp0EW$PakSTR[is.na(shp0EW$PakSTR)] <- 0
shp0EW$NieSTR <- short$Nigeria[oo]
shp0EW$NieSTR[is.na(shp0EW$NieSTR)] <- 0
shp0EW$WestSTR <- short$Western.Africa[oo]
shp0EW$MeaSTR <- short$Middle.East[oo]*0.0616  # accountin for Syria being v small
shp0EW$MeaSTR[is.na(shp0EW$MeaSTR)] <- 0
shp0EW$SasSTR <- short$Southern.Asia[oo]*0.0189  # accounting for Afg being small
shp0EW$SasSTR[is.na(shp0EW$SasSTR)] <- 0

plot(shp0EW$PakSTR,shp0EW$PakLTR)
plot(shp0EW$NieSTR,shp0EW$NieLTR)
tmp <- shp0EW@data

# vaccination
oo <- match(shp0EW$CODE,pent$Geography) 
shp0EW$pent <- pent$pentaAvg[oo]
# fill in ones without data from country level stuff
oo <- match(shp0EW$CODE,short$Geography)
shp0EW$CCode <- short$CountyCode[oo]
shp0EW$County <- short$County[oo]
aa <- match(pent$Geography,shp0EW$CCode)
empt <- which(is.na(shp0EW$pent))
# loop through...
for(i in 1:length(empt)){
  oo <- shp0EW$CCode[empt[i]]
  o1 <- which(as.character(pent$Geography)==as.character(oo))
  if(length(o1)==1){
    shp0EW$pent[empt[i]] <- pent$pentaAvg[o1]
  }
}
tmp <- shp0EW@data[is.na(shp0EW$pent),]
# just 2 left, make them the avg of everything else
shp0EW$pent[is.na(shp0EW$pent)] <- mean(shp0EW$pent,na.rm = T)

# lovely plots made easy! https://cran.r-project.org/web/packages/tmap/vignettes/tmap-nutshell.html

# STEP3: estimate Risk from data
# note that this has been updated 

# Risk = Circulation^(visiting abroad (= MA) + visitors from abroad (= MB))
# (considering annual risk from each relevant country)
# Parameters
r0 <- 3
a <- 1/3 # could easily change this parameter (IPV vaccinated are 1/3 less transmissible)
visit_pak <- 552833  # UK residents visiting pakistan
mig_pak <- 65776    # residents of pakistan visitng UK
n_pak <- 193.2*1000000 # population of Pakistan
incid_pak <- 0.15/10000000 # incidence 
visit_nie <- 183807  # UK residents visiting nigeria
mig_nie <- 100904    # residents of nigeria visitng UK
n_nie <- 186*1000000 # population of Nigeria
incid_nie <- 0.01/10000000 # incidence 
visit_afg <- 15351  # UK residents visiting Afghanistan
mig_afg <- 5127    # residents of nigeria visitng UK
n_afg <- 34.6*1000000 # population of Nigeria
incid_afg <- 0.45/10000000 # incidence 
visit_syr <- 1000  # UK residents visiting Syria
mig_syr <- 3000    # residents of syria visitng UK
n_syr <- 18.4*1000000 # population of Syria
incid_syr <- 1.34/10000000 #(half it to see what the impact is) # incidence 

shp0EW$circ <- (1/(r0*((shp0EW$pent/100)*(a-1)+1)))

# *** Risk from Pakistan
shp0EW$MA_pak <- visit_pak*(shp0EW$PakLTR/sum(shp0EW$PakLTR,na.rm=T))*(1-(shp0EW$pent/100))*incid_pak
shp0EW$MB_pak <- mig_pak*(shp0EW$PakSTR/sum(shp0EW$PakSTR,na.rm=T))*incid_pak*(1-0.75)
shp0EW$risk_pak <- 1-((1-shp0EW$circ)^(shp0EW$MA_pak+shp0EW$MB_pak))
hist(1-shp0EW$risk_pak)
plot(shp0EW$pent,shp0EW$risk_pak)
#xx <- seq(0.1,0.99,length.out=10)
#yy <- 1-((1-(1/(xx)))^10)
#plot(xx,yy)
shp0EW$risk_pakLN <- log(shp0EW$risk_pak)
shp0EW$risk_pakLN[is.infinite(shp0EW$risk_pakLN)] <- -30

#tm_shape(shp0EW) + tm_fill(col="risk_pakLN",palette="-Reds") 
plot(shp0EW$pent,shp0EW$risk_pak)    # how do the parameters influence estimated risk?
plot(shp0EW$PakLTR,shp0EW$risk_pak)
plot(shp0EW$PakSTR,shp0EW$risk_pakLN)  # location of visitors is skewed to specific locations

# *** risk from Nigeria
shp0EW$MA_nie <- visit_nie*(shp0EW$NieLTR/sum(shp0EW$NieLTR,na.rm=T))*(1-(shp0EW$pent/100))*incid_nie
shp0EW$MB_nie <- mig_nie*(shp0EW$NieSTR/sum(shp0EW$NieSTR,na.rm=T))*incid_nie*(1-0.42)
shp0EW$risk_nie <- 1-((1-shp0EW$circ)^(shp0EW$MA_nie+shp0EW$MB_nie))
hist((shp0EW$risk_nie))   # risk seems lower than from Pakistan?
shp0EW$risk_nieLN <- log(shp0EW$risk_nie)
shp0EW$risk_nieLN[is.infinite(shp0EW$risk_nieLN)] <- -30

#tm_shape(shp0EW) + tm_fill(col="risk_nieLN",palette="-Reds") 
plot(shp0EW$pent,shp0EW$risk_nie)    # how do the parameters influence estimated risk?
plot(shp0EW$NieLTR,shp0EW$risk_nie)
plot(shp0EW$NieSTR,shp0EW$risk_nie)  # location of visitors is skewed to specific locations

# *** risk from Afghanistan
# note that we're using "southern Asia" data but proportions so should be ok?
shp0EW$MA_afg <- visit_afg*(shp0EW$SasLTR/sum(shp0EW$SasLTR,na.rm=T))*(1-(shp0EW$pent/100))*incid_afg
shp0EW$MB_afg <- mig_afg*(shp0EW$SasSTR/sum(shp0EW$SasSTR,na.rm=T))*incid_afg*(1-0.65)
shp0EW$risk_afg <- 1-((1-shp0EW$circ)^(shp0EW$MA_afg+shp0EW$MB_afg))
hist(shp0EW$risk_afg)  # about same as Nie
shp0EW$risk_afgLN <- log(shp0EW$risk_afg)
shp0EW$risk_afgLN[is.infinite(shp0EW$risk_afgLN)] <- -30

#tm_shape(shp0EW) + tm_fill(col="risk_afgLN",palette="-Reds") 
plot(shp0EW$pent,shp0EW$risk_afg)    # how do the parameters influence estimated risk?
plot(shp0EW$SasLTR,shp0EW$risk_afg)
plot(shp0EW$SasSTR,shp0EW$risk_afg)  # location of visitors is skewed to specific locations

# compare risks between countries
plot(shp0EW$PakLTR,shp0EW$risk_pakLN,col="green",ylim=c(-26,-5))
points(shp0EW$SasLTR,shp0EW$risk_afgLN,col="darkgreen")
points(shp0EW$NieLTR,shp0EW$risk_nieLN,col="red")

# *** risk from Syria
# note that we're using "middle east" data but proportions so should be ok?
shp0EW$MA_syr <- visit_syr*(shp0EW$MeaLTR/sum(shp0EW$MeaLTR,na.rm=T))*(1-(shp0EW$pent/100))*incid_syr
shp0EW$MB_syr <- mig_syr*(shp0EW$MeaSTR/sum(shp0EW$MeaSTR,na.rm=T))*incid_syr*(1-0.48)
shp0EW$risk_syr <- 1-((1-shp0EW$circ)^(shp0EW$MA_syr+shp0EW$MB_syr))
hist(shp0EW$risk_syr)   # about same as Nie
shp0EW$risk_syrLN <- log(shp0EW$risk_syr)
shp0EW$risk_syrLN[is.infinite(shp0EW$risk_syrLN)] <- -30

# compare risks between countries
plot(log(shp0EW$PakLTR+1),shp0EW$risk_pakLN,col="green",ylim=c(-26,-5))
points(log(shp0EW$SasLTR+1),shp0EW$risk_afgLN,col="darkgreen")
points(log(shp0EW$NieLTR+1),shp0EW$risk_nieLN,col="red")
points(log(shp0EW$MeaLTR+1),shp0EW$risk_syrLN,col="blue")

# add them all together
# Note: this could be updated to have WPV1 and VDPV2 risk seperately
shp0EW$risk_all <- shp0EW$risk_pak + shp0EW$risk_afg + shp0EW$risk_nie + shp0EW$risk_syr  # normal scale

summary(shp0EW$risk_all)
hist(shp0EW$risk_all)

brk <- seq(from=0,to=2e-04,length.out=10)
tm_shape(shp0EW) + tm_fill(col="risk_all",palette="YlOrRd",breaks=brk)

# we want to know which factors are driving the risk - vacc or specific movements...
# would CDA be appropriate here? use R2 from linear model *but outcome should be logged and anova
dat <- shp0EW@data

# baseline
m1 <- lm(risk_all ~ 1,data=dat)
summary(m1)  #
m2 <- lm(risk_all ~ pent,data=dat)
summary(m2)  # r2 = 0.2081
 
anova(m1,m2)

m3 <- lm(risk_all ~ MeaLTR,data=dat)
summary(m3)  # r2 = 0.1393

m5 <- lm(risk_all ~ NieLTR,data=dat)
summary(m5)  # r2 = 0.1946
m6 <- lm(risk_all ~ MeaSTR,data=dat)
summary(m6)  # r2 = 0.2425 
m7 <- lm(risk_all ~ MeaLTR,data=dat)
summary(m7)  # r2 = 0.3127

  

# for plotting we want -30 to -10 (-30,-20,-15,-10) or something like that?
# Afg
# Pak
# Nie
# Syr
# Pent

# names of the big ones - label LAs with pebta below 90% (just london boroughs...)
shp0EW$pentaNAME <- ""
oo <- shp0EW$pent<90 & !is.na(shp0EW$pent)
shp0EW$pentaNAME[oo] <- as.character(shp0EW$NAME[oo])
shp0EW$pentaNAME <- gsub(" London Borough","",shp0EW$pentaNAME)

shp0EW$LDNcode <- substr(shp0EW$CODE,1,3)
shpLDN <- shp0EW[shp0EW$LDNcode=="E09",]

brk <- c(0,80,85,90,95,100)
pdf(paste0(fig_folder,"penta_all_Oct18.pdf"),height=6,width=4)
tm_shape(shp0EW) + tm_fill(col="pent",palette="-YlOrRd",breaks=brk) + tm_borders("grey30") #+ tm_text("pentaNAME",size=0.5) 
dev.off()
pdf(paste0(fig_folder,"penta_lnd_Oct18.pdf"),height=6,width=4)
tm_shape(shpLDN) + tm_fill(col="pent",palette="-YlOrRd",breaks=brk) + tm_text("pentaNAME",size=1) + tm_borders("grey30") 
dev.off()

summary(shp0EW$risk_all)
brk <- round(seq(from=0,to=2e-04,length.out=10),5) #seq(-130,-50,10) #c(-30,-20,-17.5,-15,-12.5,-10)
# add text
#        2.5%          50%          95% 
#1.741136e-08 1.122031e-06 3.191407e-05 
#n =                      18
shp0EW$riskNAME <- ""
shp0EW$riskNAME[shp0EW$risk_all>3.191407e-05] <- as.character(shp0EW$NAME[shp0EW$risk_all>3.191407e-05])
# but not London
shp0EW$riskNAME[shp0EW$County=="Outer London"] <- ""
shpLDN$riskNAME <- ""
shpLDN$riskNAME[shpLDN$risk_all>3.191407e-05] <- as.character(shpLDN$NAME[shpLDN$risk_all>3.191407e-05])

# remove extra text
#shp0EW$riskNAME <- gsub(" London Borough","",shp0EW$riskNAME)
shp0EW$riskNAME <- gsub(" District /(B/)","",shp0EW$riskNAME)
shp0EW$riskNAME <- gsub("/(B/)","",shp0EW$riskNAME)
shpLDN$riskNAME <- gsub(" London Borough","",shpLDN$riskNAME)


# use multiplot?
source("~/Dropbox (VERG)/Dropbox/Zika/Rfiles/multiplot.r")
p1 <- tm_shape(shp0EW) + tm_fill(col="risk_all",palette="YlOrRd",breaks=brk, legend.title = "Risk") +
  tm_layout(frame=FALSE) + tm_borders("grey30") + tm_text("riskNAME",size=1,col="black",auto.placement=T)
p2 <- tm_shape(shpLDN) + tm_fill(col="risk_all",palette="YlOrRd",breaks=brk) + 
  tm_legend(show=FALSE) + tm_layout(frame=FALSE) + tm_borders("grey30") + tm_text("riskNAME",size=1,col="black",auto.placement=T)

pdf(paste0(fig_folder,"risk_all_Oct18.pdf"),height=6,width=8)
multiplot(p1,p2,layout=matrix(c(1,1,2), ncol=3, byrow=TRUE))
dev.off()

pdf(paste0(fig_folder,"risk_allEW_Oct18.pdf"),height=6,width=4)
p1
dev.off()
pdf(paste0(fig_folder,"risk_allL_Oct18.pdf"),height=4,width=4)
p2
dev.off()

# vaccination
brk <- c(0,80,85,90,95,98,100)
p1 <- tm_shape(shp0EW) + tm_fill(col="pent",palette="-GnBu",breaks=brk, legend.title = "Risk") +
  tm_layout(frame=FALSE)  + tm_borders("grey30") 
p2 <- tm_shape(shpLDN) + tm_fill(col="pent",palette="-GnBu",breaks=brk) + 
  tm_legend(show=FALSE) + tm_layout(frame=FALSE)  + tm_borders("grey30") 

pdf(paste0(fig_folder,"pent_all_Oct18.pdf"),height=6,width=8)
multiplot(p1,p2,layout=matrix(c(1,1,2), ncol=3, byrow=TRUE))
dev.off()


brk <- seq(from=0,to=2e-04,length.out=10) #c(-30,-20,-17.5,-15,-12.5,-10)
pdf(paste0(fig_folder,"riskpak_all_Oct18.pdf"),height=6,width=4)
tm_shape(shp0EW) + tm_fill(col="risk_pak",palette="YlOrRd",breaks=brk)  
dev.off()

pdf(paste0(fig_folder,"riskpak_lnd_Oct18.pdf"),height=4,width=4)
tm_shape(shpLDN) + tm_fill(col="risk_pakLN",palette="YlOrRd",breaks=brk) # + tm_text("pentaNAME",size=1) 
dev.off()

pdf(paste0(fig_folder,"risknie_all_Aug18.pdf"),height=6,width=4)
tm_shape(shp0EW) + tm_fill(col="risk_nieLN",palette="-YlOrRd",breaks=brk)  
dev.off()
pdf(paste0(fig_folder,"risknie_lnd_Aug18.pdf"),height=4,width=4)
tm_shape(shpLDN) + tm_fill(col="risk_nieLN",palette="-YlOrRd",breaks=brk) # + tm_text("pentaNAME",size=1) 
dev.off()

pdf(paste0(fig_folder,"riskafg_all_Aug18.pdf"),height=6,width=4)
tm_shape(shp0EW) + tm_fill(col="risk_afgLN",palette="-YlOrRd",breaks=brk)  
dev.off()
pdf(paste0(fig_folder,"riskafg_lnd_Aug18.pdf"),height=4,width=4)
tm_shape(shpLDN) + tm_fill(col="risk_afgLN",palette="-YlOrRd",breaks=brk) # + tm_text("pentaNAME",size=1) 
dev.off()

pdf(paste0(fig_folder,"risksyr_all_Aug18.pdf"),height=6,width=4)
tm_shape(shp0EW) + tm_fill(col="risk_syrLN",palette="-YlOrRd",breaks=brk)  
dev.off()
pdf(paste0(fig_folder,"risksyr_lnd_Aug18.pdf"),height=4,width=4)
tm_shape(shpLDN) + tm_fill(col="risk_syrLN",palette="-YlOrRd",breaks=brk) # + tm_text("pentaNAME",size=1) 
dev.off()

# proportion of risk
shp0EW$risk_prop <- shp0EW$risk_all/sum(shp0EW$risk_all)
hist(shp0EW$risk_prop)

# add population size
oo <- match(shp0EW$CODE,resid$Geography)  # length(oo)==dim(shp0EW)[1] yes!
shp0EW$PopSize <- resid$TotalPop2011Census[oo]
shp0EW$PopSizePRCT <- shp0EW$PopSize/sum(shp0EW$PopSize,na.rm=T)

tmp <- shp0EW@data

# so we want to illustrate the impoprtance of targetted sampling
tmp <- tmp[order(as.numeric(tmp$risk_all),decreasing = T),]
tmp$order <- 1:dim(tmp)[1]
tmp$crisk <- cumsum(tmp$risk_prop)
tmp$cpop <- cumsum(tmp$PopSizePRCT)
summary(tmp$cpop)
tmp[1:10,]
tmp$crisk_text <- ""
tmp$crisk_text[tmp$crisk<0.2] <- as.character(tmp$NAME[tmp$crisk<0.2])
# remove text
tmp$crisk_text <- gsub(" London Borough","",tmp$crisk_text)

# targetted surveillance
pdf(paste0(fig_folder,"Targeted_risk_08Oct18.pdf"),height=5,width=5)
plot(tmp$cpop*100,tmp$crisk*100,pch=19,xlab="Culmulative Population of England and Wales",
     ylab="Culmulative Percentage of Estimated Risk")
text(tmp$cpop*100,tmp$crisk*100,labels=as.character(tmp$crisk_text),pos=4,cex=0.9)
abline(a=0,b=1,lty=2,col="grey50")
dev.off()

# write to table
write.csv(tmp,"summary_polio_risk_16Oct18.csv",row.names=F)

tmp <- rbind(tmp1,tmp2,tmp3)
out <- data.frame(rbind(tmp1,tmp2,tmp3))
names(out) <- c("CODE","name","penta","out","pak","nie","risk","level")
head(out)
out$risk <- as.numeric(as.character(out$risk))
# no obs have a risk > 0, essentially because all have high enough coverage
write.table(tmp,"risk_all_regions.csv",sep=",",row.names=F,col.names=T)

 
 
# end