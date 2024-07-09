# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +++                                                    +++
# +++        Binomial GLMM for age-0 & age-1 sole        +++
# +++                                                    +++
# +++ Charlotte Van Moorleghem & Jochen Depestele (ILVO) +++
# +++                    February 2022                   +++
# +++                                                    +++
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(data.table)
library(dplyr)
library(lme4)
library(sf)
library(AICcmodavg)
library(beepr)
library(DHARMa)
library(tidyr)

tosave = FALSE

### set LOCAL folder directory ----
locdir <- rstudioapi::getActiveProject()
inputdir <- paste0(locdir,"/input")
outputdir <- paste0(locdir,"/output")

### Data set-up  ----
sf_SolFinal <- readRDS("C:/Users/cvanmoorleghem/OneDrive - ILVO/Development/RStudio/SEAwise T5.4/seawise_t5-4_juvsol7afg/output/sf_Solea_solea_final_withAbs&Predictors_YFS+BTS+BTS-VIII+DYFS+SNS_1981til2021__in3a204a4b4c7a7d7e7f7g7h7j28a8b.RDS")

### Data preparation ----

# Gear needs to be a factor if it is going to be a random effect in the model
sf_SolFinal$Gear <- as.factor(as.character(sf_SolFinal$Gear))
sf_SolFinal$Folk_7cl_t <- as.factor(as.character(sf_SolFinal$Folk_7cl_t))

# set very large distance for na values in Nearest_OWF variables
# 901240.1 m is max distance to nearest fully commissioned OWF
sf_SolFinal$fc_NearestOWF <- ifelse(is.na(sf_SolFinal$fc_NearestOWF),9999999,sf_SolFinal$fc_NearestOWF)
# 1052155.3 m is max distance to nearest OWF under construction
sf_SolFinal$uc_NearestOWF <- ifelse(is.na(sf_SolFinal$uc_NearestOWF),9999999,sf_SolFinal$uc_NearestOWF)
#' [might not be correct approach!]

# reduce dataset so it doesn't contain NAs
subSol_minSex <- as.data.table(sf_SolFinal)
subSol_minSex[,geometry:=NULL]
subSol_minSex <- subSol_minSex %>%
  group_by(haulID) %>%
  summarise(Year = first(Year), LengthCatMin = first(LengthCatMin), LengthCatMax = first(LengthCatMax), TotalNo = first(TotalNo), Age0 = sum(Age0), Age1 = sum(Age1), Age2 = sum(Age2), Age3 = sum(Age3), Age4 = sum(Age4), Age5 = sum(Age5), Age6 = sum(Age6), 
            bathymetry = first(bathymetry), SurCurDir = first(SurCurDir), SurCurSpeed = first(SurCurSpeed), WindDir = first(WindDir), WindSpeed = first(WindSpeed), SwellDir = first(SwellDir), SwellHeight = first(SwellHeight), SurTemp = first(SurTemp), BotTemp = first(BotTemp), SurSal = first(SurSal), BotSal = first(BotSal), SecchiDepth = first(SecchiDepth), TidePhase = first(TidePhase),
            Survey = first(Survey), Country = first(Country), Ship = first(Ship), Gear = first(Gear), GearEx = first(GearEx), GearDescription = first(GearDescription), singlegearwidth_m = first(singlegearwidth_m), gearwidth_m = first(gearwidth_m), Rigging = first(Rigging), Tickler = first(Tickler), CodendMesh = first(CodendMesh), Quarter = first(Quarter), Month = first(Month), Day = first(Day), HaulNo = first(HaulNo), 
            StNo = first(StNo), Area_27 = first(Area_27), ShootLat = first(ShootLat), ShootLong = first(ShootLong), HaulLat = first(HaulLat), HaulLong = first(HaulLong), HaulDur = first(HaulDur), HaulDist = first(HaulDist), GroundSpeed = first(GroundSpeed), SampleArea_m2 = first(SampleArea_m2), Cruise = first(Cruise), haulID = first(haulID),
            HumPop10km = first(HumPop10km), HumPop25km = first(HumPop25km), HumPop50km = first(HumPop50km), fc_OWFcapacity = first(fc_OWFcapacity), fc_NearestOWF = first(fc_NearestOWF), production_Polychaeta_Peracarida = first(production_Polychaeta_Peracarida),uc_NoTurbin = first(uc_NoTurbin), uc_NearestOWF = first(uc_NearestOWF), Folk_7cl_t = first(Folk_7cl_t))


subSol_minSex$Age0_PA <- ifelse(subSol_minSex$Age0 == 0, 0, 1)
subSol_minSex$Age1_PA <- ifelse(subSol_minSex$Age1 == 0, 0, 1)


subSol_minSex <- subSol_minSex %>% 
  filter(!is.na(HumPop25km)) %>%
  filter(!is.na(production_Polychaeta_Peracarida)) %>%
  filter(!is.na(Folk_7cl_t)) %>% 
  filter(!is.na(Month))

# rescale variables
subSol_minSex <- subSol_minSex %>%
  mutate(HumPop25km_scal = as.numeric(scale(HumPop25km)),
         year_scal = as.numeric(scale(Year)),
         quarter_scal = as.numeric(scale(Quarter)),
         month_scal = as.numeric(scale(Month)),
         bathymetry_scal = as.numeric(scale(bathymetry)),
         uc_NoTurbin_scal = as.numeric(scale(uc_NoTurbin)),
         uc_NearestOWF_scal = as.numeric(scale(uc_NearestOWF)),
         production_Polychaeta_Peracarida_scal = as.numeric(scale(production_Polychaeta_Peracarida)),
         SampleArea_m2_scal = as.numeric(scale(SampleArea_m2)),
         fc_OWFcapacity_scal = as.numeric(scale(fc_OWFcapacity)),
         fc_NearestOWF_scal = as.numeric(scale(fc_NearestOWF)))


# figure 1
mapview(subSol_minSex[which(subSol_minSex$Age0_PA == 1),])
mapview(subSol_minSex[which(subSol_minSex$Age1_PA == 1),])

### modelling
# assess multicollinearity
x <- sf_SolFinal[c("year_scal","month_scal","bathymetry_scal","Median.grain.size_scal","dist2coast_scal","HumPop10km_scal","uc_NoTurbin_scal","uc_NearestOWF_scal","fc_OWFcapacity_scal","fc_NearestOWF_scal","production_Polychaeta_Peracarida_scal","SampleArea_m2_scal","CertaintyIdx_Age0")]
st_geometry(x) <- NULL
cor(x, method = "pearson", use = "complete.obs")
# generally ok, highest values are:
# -0.56121921 for year & fc_NearestOWF; the later the year the nearer the nearest OWF (logical, as there are more constructed over the years)
# 0.47170938 for bathymetry and dist2coast; the nearer to the coast the shallower the water
plot(sf_SolFinal$dist2coast_scal,sf_SolFinal$bathymetry_scal)
plot(sf_SolFinal[which(sf_SolFinal$fc_NearestOWF != 9999999),]$fc_NearestOWF_scal,sf_SolFinal[which(sf_SolFinal$fc_NearestOWF != 9999999),]$year_scal)

##set up candidate models for the final step within the stepwise addition procedure for zero-aged sole
Cand.mod <- list()
##global model
Cand.mod[[1]]<-glmer(Age0_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+month_scal+fc_OWFcapacity_scal+HumPop25km_scal+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
Cand.mod[[2]]<-glmer(Age0_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+month_scal+fc_OWFcapacity_scal+HumPop25km_scal+year_scal+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
Cand.mod[[3]]<-glmer(Age0_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+month_scal+fc_OWFcapacity_scal+HumPop25km_scal+uc_NoTurbin_scal+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
Cand.mod[[4]]<-glmer(Age0_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+month_scal+fc_OWFcapacity_scal+HumPop25km_scal+uc_NearestOWF_scal+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
Cand.mod[[5]]<-glmer(Age0_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+month_scal+fc_OWFcapacity_scal+HumPop25km_scal+fc_NearestOWF_scal+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
Cand.mod[[6]]<-glmer(Age0_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+month_scal+fc_OWFcapacity_scal+HumPop25km_scal+production_Polychaeta_Peracarida_scal+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

confset(cand.set = Cand.mod, second.ord = FALSE, method = "raw")

# Confidence set for the best model
# 
# Method:	 raw sum of model probabilities
# 
# 95% confidence set:
#       K      AIC Delta_AIC AICWt
# Mod2 14 14832.18      0.00  0.37
# Mod3 14 14832.85      0.67  0.26
# Mod4 14 14833.14      0.95  0.23
# Mod1 13 14835.71      3.53  0.06
# Mod5 14 14836.64      4.46  0.04
# 
# Model probabilities sum to 0.96 
# 
# Warning message:
#   In confset(cand.set = Cand.mod, second.ord = FALSE, method = "raw") : 
#   Model names have been supplied automatically in the table


##set up candidate models
Cand.mod <- list()
##global model
Cand.mod[[1]]<-glmer(Age1_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
Cand.mod[[2]]<-glmer(Age1_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+year_scal+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
Cand.mod[[3]]<-glmer(Age1_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+month_scal+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
Cand.mod[[4]]<-glmer(Age1_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+HumPop25km_scal+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
Cand.mod[[5]]<-glmer(Age1_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+uc_NoTurbin_scal+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
Cand.mod[[6]]<-glmer(Age1_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+uc_NearestOWF_scal+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
Cand.mod[[7]]<-glmer(Age1_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+fc_OWFcapacity_scal+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
Cand.mod[[8]]<-glmer(Age1_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+fc_NearestOWF_scal+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
Cand.mod[[9]]<-glmer(Age1_PA~SampleArea_m2_scal+bathymetry_scal+Folk_7cl_t+production_Polychaeta_Peracarida_scal+(1|Gear),data=subSol_minSex,family=binomial(link="logit"),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

confset(cand.set = Cand.mod, second.ord = FALSE, method = "raw")

# Confidence set for the best model
# 
# Method:	 raw sum of model probabilities
# 
# 95% confidence set:
#       K      AIC Delta_AIC AICWt
# Mod4 11 16648.83      0.00  0.46
# Mod6 11 16649.91      1.08  0.27
# Mod5 11 16649.96      1.13  0.26
# 
# Model probabilities sum to 0.98 

plot(simulateResiduals(Cand.mod[[1]]),cex =10)
#' [all models had significantly deviating residuals according to KS (Kolmogorov-Smirnov) test]

simulationOutput <- simulateResiduals(fittedModel = Cand.mod[[1]])
testUniformity(simulationOutput = simulationOutput)
# test heterogeneity

lsmeans(Cand.mod[[2]],pairwise ~ Folk_7cl_t, adjust = "bonferroni", type = "response")

# build figures for SEAwise deliverable T5.4
par(mfrow=c(2,2))	
par(mar=c(15, 4, 4, 2) + 0.1)	# spelen met margin van graph

figSubstrate <- data.frame(gemiddelde=tapply(subSol_minSex$Age0_PA,subSol_minSex$Folk_7cl_t,mean),
                           SD=tapply(subSol_minSex$Age0_PA,subSol_minSex$Folk_7cl_t,sd),
                           SE=tapply(subSol_minSex$Age0_PA,subSol_minSex$Folk_7cl_t,sd)/sqrt(tapply(subSol_minSex$Age0_PA,subSol_minSex$Folk_7cl_t,length)))

figSubstrate
attach(figSubstrate)
grafiek<-NULL
grafiek<-barplot(gemiddelde,main="0-aged Sole",las=2,ylim=c(0,max(gemiddelde)+max(SE)), names.arg = c("Mud","Sandy Mud", "Muddy Sand", "Sand", "Coarse-grained\n Sediment","Mixed Sediment", "Rock & Boulders"), col = "steelblue", cex = 2, cex.main = 2, cex.axis = 1.5)
arrows(x0=grafiek,y0=gemiddelde-SE,x1=grafiek,y1=gemiddelde+SE,lwd=1.5,length=0.07,code=3,angle=90)
box()
detach(figSubstrate)

figSubstrate <- data.frame(gemiddelde=tapply(subSol_minSex$Age1_PA,subSol_minSex$Folk_7cl_t,mean),
                           SD=tapply(subSol_minSex$Age1_PA,subSol_minSex$Folk_7cl_t,sd),
                           SE=tapply(subSol_minSex$Age1_PA,subSol_minSex$Folk_7cl_t,sd)/sqrt(tapply(subSol_minSex$Age1_PA,subSol_minSex$Folk_7cl_t,length)))
figSubstrate
attach(figSubstrate)
grafiek<-NULL
grafiek<-barplot(gemiddelde,main="1-aged Sole",las=2,ylim=c(0,max(gemiddelde)+max(SE)), names.arg = c("Mud","Sandy Mud", "Muddy Sand", "Sand", "Coarse-grained\n Sediment","Mixed Sediment", "Rock & Boulders"), col = "steelblue", cex = 2, cex.main = 2, cex.axis = 1.5)
arrows(x0=grafiek,y0=gemiddelde-SE,x1=grafiek,y1=gemiddelde+SE,lwd=1.5,length=0.07,code=3,angle=90)
box()
detach(figSubstrate)

par(mar=c(9, 7, 4, 2) + 0.1) # restore margins

range(subSol$bathymetry)
subSol_minSex$depthCat <- cut(subSol_minSex$bathymetry,
                              breaks=seq(0, 250, 10))

figDepth <- data.frame(gemiddelde=tapply(subSol_minSex$Age0_PA,subSol_minSex$depthCat,mean),
                       SD=tapply(subSol_minSex$Age0_PA,subSol_minSex$depthCat,sd),
                       SE=tapply(subSol_minSex$Age0_PA,subSol_minSex$depthCat,sd)/sqrt(tapply(subSol_minSex$Age0_PA,subSol_minSex$depthCat,length)))
figDepth <- figDepth[1:11,]
attach(figDepth)
grafiek<-NULL
grafiek<-barplot(gemiddelde,las=2,ylim=c(0,max(gemiddelde)+max(SE)), names.arg = c("[0-10[","[10-20[","[20-30[","[30-40[","[40-50[","[50-60[","[60-70[","[70-80[","[80-90[","[90-100[", "[100-110["), col = "steelblue", cex = 2, cex.main = 2, cex.axis = 1.5)
title(xlab = "water depth (m)", mgp=c(8,1,0),cex.lab=2)
arrows(x0=grafiek,y0=gemiddelde-SE,x1=grafiek,y1=gemiddelde+SE,lwd=1.5,length=0.07,code=3,angle=90)
box()
detach(figDepth)

figDepth <- data.frame(gemiddelde=tapply(subSol_minSex$Age1_PA,subSol_minSex$depthCat,mean),
                       SD=tapply(subSol_minSex$Age1_PA,subSol_minSex$depthCat,sd),
                       SE=tapply(subSol_minSex$Age1_PA,subSol_minSex$depthCat,sd)/sqrt(tapply(subSol_minSex$Age1_PA,subSol_minSex$depthCat,length)))
figDepth <- figDepth[1:13,]
attach(figDepth)
grafiek<-NULL
grafiek<-barplot(gemiddelde,las=2,ylim=c(0,max(gemiddelde)+max(SE)), names.arg = c("[0-10[","[10-20[","[20-30[","[30-40[","[40-50[","[50-60[","[60-70[","[70-80[","[80-90[","[90-100[", "[100-110[", "[110-120[","[120-130["), col = "steelblue", cex = 2, cex.main = 2, cex.axis = 1.5)
arrows(x0=grafiek,y0=gemiddelde-SE,x1=grafiek,y1=gemiddelde+SE,lwd=1.5,length=0.07,code=3,angle=90)
box()
detach(figDepth)

# Use model estimates for OWF capacity and human population:
effect_plot(Cand.mod[[1]], pred = fc_OWFcapacity_scal, main.title = "OWF capacity", interval = TRUE, y.label = "Probability of 0-aged sole catch", plot.points = TRUE, point.alpha = 0.2, x.label = "OWF capacity (scaled)")
effect_plot(Cand.mod[[1]], pred = HumPop25km_scal, main.title = "Human Population", interval = TRUE, y.label = "", plot.points = TRUE, point.alpha = 0.2, x.label = "Coastal Population (scaled)")


