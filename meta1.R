getwd()
#setwd("R/R2d2/CC_Meta")
library(dplyr)
library(metafor)
library(glmulti)
m1<-(read.csv("CoverCropMeta3.csv"))
str(m1)

####FINAL CODE GOES HERE
##1. ORGANIZE factors
##2. Estimate Bulk Density
##3. Estimate Error
##4. Escalc
##5. Interactions
##6. Graph

############################################
############################################
######                       ###############
###### Organize Factors  (Cstock)###########
######                       ###############
############################################
############################################
dim(m1)

#summary(m1$CC.functional.group)
levels(m1$CC.functional.group)<-c("forb","graminoid","legume","legume","multiple","multiple.leg")
#summary(m1$CC.functional.group)


summary(m1$tillage)
levels(m1$tillage)<-c(NA,"Till","Till","NoTill","Till")
summary(m1$tillage)
summary(m1$Nitrogen)
levels(m1$Nitrogen)<-c(NA,"high","low","None","low")
summary(m1$Nitrogen)
summary(m1$Crop.type)
m1$Crop1<-m1$Crop.type
levels(m1$Crop.type)<-c("forb","graminoid","legume","multi","multi.leg",NA,"graminoid","forb")

summary(m1$Depth.for.MBC)
levels(m1$Depth.for.MBC)<-c(NA,5,10,12,15,20,30,5,7.5)


######DEPTH
m1$stdd3<-ifelse(is.na(m1$stdd3),15,m1$stdd3)
m1$depth.1<-ifelse((m1$stdd3<10.0),'shallow','medium')
m1$depth.1<-ifelse((m1$stdd3>15.1),'deep',m1$depth.1)

#########WEIGHTING VECTOR

#To weight studies based on number of replicates and duration of study. (Terrar, 2016, deGraaf, 2006; vanGroenigen 2006)
m1$Fc<-((m1$reps*m1$reps)/(m1$reps+m1$reps)+(m1$Years.in.experiment*m1$Years.in.experiment)/(m1$Years.in.experiment+m1$Years.in.experiment))
m1$Fc2<-m1$Fc/2

####################

############################################
############################################
######                       ###############
#####  BuLkDeNsItY_estimate  ###############
######                       ###############
############################################
############################################
############################################

BD.est<-(-0.011905*((m1$CC.soil.carbon.g.kg.stdD+m1$Cont.soil.C.g.kg.stdD)/2)-0.004203*(m1$stdd3)+1.562830)
BDest2<-(1.662094-0.146275*(log(m1$CC.carbon.Mg.ha.stdD+m1$Cont.C.Mg.ha.stdD)/2)+0.006104*m1$stdd3)

CC.Ccon.est<-(m1$CC.carbon.Mg.ha.stdD/(0.1*BDest2*m1$stdd3))
Cont.Ccon.est<-m1$Cont.C.Mg.ha.stdD/(0.1*BDest2*m1$stdd3)
CCest.Mg.ha<-m1$CC.soil.carbon.g.kg.stdD*BD.est*m1$stdd3*0.1 #g.kg*BD*d/10`
Contest.Mg.ha<-m1$Cont.soil.C.g.kg.stdD*BD.est*m1$stdd3*0.1

############KEEP THESE
m1$CC.g.kg.BDest<-ifelse(is.na(m1$CC.soil.carbon.g.kg.stdD),CC.Ccon.est,m1$CC.soil.carbon.g.kg.stdD)
m1$Cont.g.kg.BDest<-ifelse(is.na(m1$Cont.soil.C.g.kg.stdD),Cont.Ccon.est,m1$Cont.soil.C.g.kg.stdD)
m1$CC.Cstock.BDest<-ifelse(is.na(m1$CC.carbon.Mg.ha.stdD),CCest.Mg.ha,m1$CC.carbon.Mg.ha.stdD)
m1$Cont.Cstock.BDest<-ifelse(is.na(m1$Cont.C.Mg.ha.stdD),Contest.Mg.ha,m1$Cont.C.Mg.ha.stdD)

############################################
############################################
######                       ###############
#####     StdErr_estimate  ###############
######                       ###############
############################################
############################################
############################################
######ERROR_EST (DUMMY VAR)###############
cccse<-(m1$CC.C.Mg.ha.stD.SE/(0.1*m1$BDest2*m1$stdd3))
concse<-(m1$Cont.C.Mg.ha.stdD.SE/(0.1*m1$BDest2*m1$stdd3))
ccsse<-m1$CC.C.Mg.ha.stD.SE*m1$BDest*m1$stdd3*0.1 
consse<-m1$Cont.C.Mg.ha.stdD.SE*m1$BDest*m1$stdd3*0.1


###########KEEP THESE
m1$CC.Ccon.SE.est<-ifelse(is.na(m1$CC.g.kg.StdD.SE),cccse,m1$CC.g.kg.StdD.SE)
m1$Cont.Ccon.SE.est<-ifelse(is.na(m1$ContC.k.kg.StdD.SE),concse,m1$ContC.k.kg.StdD.SE)
m1$CC.Cstock.SE.est<-ifelse(is.na(m1$CC.C.Mg.ha.stD.SE),ccsse,m1$CC.C.Mg.ha.stD.SE)
m1$Cont.Cstock.SE.est<-ifelse(is.na(m1$Cont.C.Mg.ha.stdD.SE),consse,m1$Cont.C.Mg.ha.stdD.SE)


#NOW estimate error using These four (8) ^
###########
############ESTIMATE ERROR - C Stock
##To estimate the variance
#What is the ratio of SD to soil C
nozero<-function(x) {
  ifelse(x<=0,0.0001,x)}
###########SD calc
m1$CC.stdD.Mg.ha.SD<-(m1$CC.Cstock.SE.est)*(sqrt(m1$reps))
m1$Cont.stdD.Mg.ha.SD<-(m1$Cont.Cstock.SE.est)*(sqrt(m1$reps))
###########SD estimate

summary(C<-c(m1$CC.Cstock.BDest,m1$Cont.Cstock.BDest))
length(C)/2
t<-c(rep("Cover",384),rep("control",384))
summary(SD<-c(m1$CC.stdD.Mg.ha.SD,m1$Cont.stdD.Mg.ha.SD))
m3<-cbind.data.frame(C,SD,t)
colnames(m3)<-c("C","SD","TrT")
hist(m3$logSD<-(log(m3$SD+2)))
hist(m3$logC<-(log(m3$C+7)))

(boxplot.stats(m3$C.SD<-m3$logSD/m3$logC))
summary(m3$C.SD)
summary(m3$logSD.est<-m3$logC*(boxplot.stats(m3$C.SD)$stats[3]))
m3$logSD.hi<-m3$logC*(boxplot.stats(m3$C.SD)$stats[2])
m3$logSD.lo<-m3$logC*(boxplot.stats(m3$C.SD)$stats[4])

m3$SD.est<-nozero(exp(m3$logSD.est)-2)
m3$SDest.hi<-nozero(exp(m3$logSD.hi)-2)
m3$SDest.low<-nozero(exp(m3$logSD.lo)-2)
summary(m3$SD.est)
dim(m3)/2
m1$CC.SD.est<-m3$SD.est[1:384]
m1$Cont.SD.est<-m3$SD.est[385:768]
m1$CC.SDest.hi<-m3$SDest.hi[1:384]
m1$Cont.SDest.hi<-m3$SDest.hi[385:768]
m1$CC.SDest.lo<-m3$SDest.low[1:384]
m1$Cont.SDest.lo<-m3$SDest.low[385:768]

#########################################
#########################################
###########################################
##########EFFECT SIZE CALCULATION##########
#########################################
#########################################
###########################################
#######################################

m2<-escalc(m1i=CC.Cstock.BDest ,#Treatment response
           sd1i=CC.SD.est,#Treatment variance
           n1i=reps,#number of reps in study
           m2i=Cont.Cstock.BDest,#Control response
           sd2i=Cont.SD.est,#Control variance
           n2i=reps,#number of control reps in study
           data=m1,
           measure="ROM")#"ROM" for the log transformed ratio of means (Hedges et al., 1999).
title.m2="Effect of cover crop on soil carbon stock (% difference)"
#########################

#Function to Extract model components
#```{r extract model results}
rma.extract <- function(m){
  estimate <- summary(m)$b
  ci.lb <- summary(m)$ci.lb
  ci.ub <- summary(m)$ci.ub
  p <- summary(m)$pval
  modinf <- data.frame(estimate=estimate,
                       ci.lb = ci.lb,
                       ci.ub = ci.ub,
                       pval=p)
  return(modinf)
}
###########################
str(m2)
############################Meta Analysis
my.rma<-rma.mv(yi,vi,Fc2,random = ~1|study.ID,data=m2)##Including Fc shows preference for well replicated, long term studies
#USING FC makes sense BECAUSE
#rma.mv(yi,vi,mods=~reps,random=~1|study.ID,data=m2)

#######################################Effect of Factors
rma.ccFgroup.Fc<-rma.mv(yi,vi,Fc2,mods= ~CC.functional.group-1,random = ~1|study.ID,data=m2)
rma.till.Fc<-rma.mv(yi,vi,Fc2,mods= ~tillage-1,random = ~1|study.ID,data=m2)

rma.nitrogen.Fc<-rma.mv(yi,vi,Fc2,mods= ~Nitrogen-1,random = ~1|study.ID,data=m2)
rma.crop.Fc<-rma.mv(yi,vi,Fc2,mods=~Crop.type-1,random=~1|study.ID,data=m2)

rma.mv(yi,vi,Fc2,mods=~Continent-1,random=~1|study.ID,data=m2)

#### Most of the continuous variables are significant without the intercept (factor - 1) but not with
rma.mv(yi,vi,Fc2,mods=~elevation-1,random=~1|study.ID,data=m2)#sig without intercept but not with
rma.mv(yi,vi,Fc2,mods=~MAT.C-1,random=~1|study.ID,data=m2)#sig without intercept but not with
rma.mv(yi,vi,Fc2,mods=~pH-1,random=~1|study.ID,data=m2)#sig without intercept but not with
rma.mv(yi,vi,Fc2,mods=~MAP-1,random=~1|study.ID,data=m2)#sig without intercept but not with
rma.mv(yi,vi,Fc2,mods=~lat.decimal-1,random=~1|study.ID,data=m2)#sig without intercept but not with

###### Control C stock is a significant predictor with or without the intercept
rma.mv(yi,vi,Fc2,mods=~Cont.Cstock.BDest-1,random=~1|study.ID,data=m2)
summary(m2$Cont.Cstock.BDest)

plot(m2$yi~m2$CC.Cstock.BDest,ylab=("Effect Size"),xlab=("Control C stock"))

###### Experimental duration is significant with or without the intercept
rma.mv(yi,vi,Fc2,mods=~Years.in.experiment-1,random=~1|study.ID,data=m2)

## Could be driven by long studies (>25 years)
## Studies < 25 years not significant with intercept but significant without
m4<-m2[!(m2$Years.in.experiment>25),]
rma.mv(yi,vi,Fc2,mods=~Years.in.experiment-1,random=~1|study.ID,data=m4)
plot(m2$yi~m2$Years.in.experiment)
abline(lm(m2$yi~m2$Years.in.experiment))

######Crop yield is not significant with or without the intercept
rma.mv(yi,vi,Fc2,mods=~ESTIMATED_CROP_YIELD_MG.ha-1,random=~1|study.ID,data=m2)
m2$log.yield<-log(m2$ESTIMATED_CROP_YIELD_MG.ha)

#######Transformed crop yield And Transformed Residue added are significant without the intercept but not with
rma.mv(yi,vi,Fc2,mods=~log.yield-1,random=~1|study.ID,data=m2)
#hist(m2$TOTAL_AG_RESIDUE_ADDITIONS_MINUS_CONTROL)
m2$log.residue<-(log(m2$TOTAL_AG_RESIDUE_ADDITIONS_MINUS_CONTROL+3))
rma.mv(yi,vi,Fc2,mods=~log.residue-1,random=~1|study.ID,data=m2)
rma.extract(rma.mv(yi,vi,Fc2,mods=~TOTAL_AG_RESIDUE_ADDITIONS_MINUS_CONTROL,data=m2))
plot(m2$yi~m2$TOTAL_AG_RESIDUE_ADDITIONS_MINUS_CONTROL)
abline(lm(m2$yi~m2$TOTAL_AG_RESIDUE_ADDITIONS_MINUS_CONTROL,na.action=na.omit))
summary(m2$TOTAL_AG_RESIDUE_ADDITIONS_MINUS_CONTROL) #######ARE THESE UNITS ALL THE SAME??????

##############MMR

m.sums<-rbind(rma.extract(my.rma),
              rma.extract(rma.ccFgroup.Fc),
              rma.extract(rma.crop.Fc),
              rma.extract(rma.till.Fc),
              rma.extract(rma.nitrogen.Fc))

###########################
############################Meta Analysis No WEIGHTING FACTOR
my.rma<-rma.mv(yi,vi,random = ~1|study.ID,data=m2)##Including Fc shows preference for long term, well replicated studies

#######################################Effect of Factors
rma.ccFgroup<-rma.mv(yi,vi,mods= ~CC.functional.group-1,random = ~1|study.ID,data=m2)
rma.till<-rma.mv(yi,vi,mods= ~tillage-1,random = ~1|study.ID,data=m2)

rma.nitrogen<-rma.mv(yi,vi,mods= ~Nitrogen-1,random = ~1|study.ID,data=m2)
rma.crop<-rma.mv(yi,vi,mods=~Crop.type-1,random=~1|study.ID,data=m2)

rma.mv(yi,vi,mods=~sand.g.kg-1,random=~1|study.ID,data=m2)

#rma.depth<-rma.mv(yi,vi,mods= ~depth.1-1,random = ~1|study.ID,data=m2)
#rma.cont<-rma.mv(yi,vi,Fc,mods= ~Continent-1,random = ~1|study.ID,data=m2)

###############################

###############################
########CSTOCK Summary
############################
m.sums<-rbind(rma.extract(my.rma),
              rma.extract(rma.ccFgroup),
              rma.extract(rma.crop),
              rma.extract(rma.till),
              rma.extract(rma.nitrogen))
##########################################              

m.sums$percent.change<-(exp(m.sums$estimate)-1)*100
m.sums$ci.lb.perc<-(exp(m.sums$ci.lb)-1)*100
m.sums$ci.ub.perc<-(exp(m.sums$ci.ub)-1)*100

m.sums$moderator<-as.factor(c(rownames(m.sums)))
length(m.sums$moderator)
m.sums$label<-c("full model",
                "forb cover","graminoid cover","legume cover","multiple cover","multiple cover \n with legume",
                "forb cash","graminoid cash","legume cash","multiple cash","multiple cash  \n with legume",
                "Till","No Till",
                "high N","low N","No N")
m.sums$moderator <- factor(c("All trials", 
                             rep("Cover crop \n group",5),
                             rep("Cash crop \n group",5),
                             rep("Tillage", 2),
                             rep("Nitrogen \nFertilization",3)),
                           
                           levels=c("All trials",
                                    "Cover crop \n group",
                                    "Cash crop \n group",
                                    "Tillage",
                                    "Nitrogen \nFertilization"))

m.sums$n<-as.character(c((dim(m2)[1]),
                         (as.data.frame(table(m2$CC.functional.group))[,2]),
                         (as.data.frame(table(m2$Crop.type))[,2]),
                         as.data.frame(table(m2$tillage))[,2],
                         as.data.frame(table(m2$Nitrogen))[,2]))
str(m.sums)
summary(m.sums$ci.ub.perc)


m.sums$label<-factor(m.sums$label,levels=m.sums$label[order(m.sums$moderator,m.sums$estimate)])

m.sums$p.value<-round(m.sums$pval,3)
m.sums$pval<-ifelse(m.sums$pval<0.001,"p < 0.001", m.sums$p.value)
m.sums$change<-round(m.sums$percent.change,1)


sig<-rep(NA,length(m.sums[,1]))
sig<-as.character(ifelse(m.sums$p.value<0.05,"*",""))
sig<-ifelse(m.sums$p.value<0.01,"**",sig)
sig<-ifelse(m.sums$p.value<0.001,"***",sig)

m.sums$sig<-sig

pdat<-cbind(m.sums[,c(10,12:13,8)])

##################################################
######################################GRAPH
library(ggplot2)

p <- ggplot(m.sums, aes(y=label)) + facet_grid(moderator ~ ., scales="free_y", space="free_y")

p <- p + geom_point(aes(x=percent.change), size=3)
p <- p + geom_errorbarh(aes(x=percent.change,xmin=ci.lb.perc,xmax=ci.ub.perc), height=0, size=0.35)
p <- p + geom_vline(xintercept=0, color="grey25", linetype="dotted")

p <- p + geom_text(aes(label=change,x=-47),size=3) #x + 17
p <- p+ geom_text(aes(label=rep(" %",16),x=-41),size=3) #x + 23
p <- p + geom_text(aes(label=n, x=-32,fontface="italic"), size=3)#x + 7
p<- p+ geom_text(aes(label="(      )",x=-31.5,fontface="italic"),size=3) #x+7.5
p <- p + geom_text(aes(label=sig, x= -16), fontface="bold", size= 3)#if this is x

p <- p + labs(x=title.m2, y=NULL,size=5)
p <- p + theme(panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),
               strip.text.y=element_text(angle=0,size=10,color="white"),
               strip.background = element_rect(fill="turquoise4"))

p <- p + scale_x_continuous(limits=c(-50,80
                                     ))


p

ggsave("cstockFc2.png",width=5,height=5,units="in",dpi=300)

m.sums
title

