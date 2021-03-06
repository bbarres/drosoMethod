##############################################################################/
##############################################################################/
#R code for the Drosophila suzukii pesticide resistance: Experiment 2
##############################################################################/
##############################################################################/

source("droso_data_load.R")


##############################################################################/
#What is the effect of the age of the flies on their LD50?####
##############################################################################/

#we select the data of phosmet test with the St Foy population, with different
#age classes and we split it between males and females
agedata<-dataDroz[dataDroz$age_comp==1,]
#aggregate the data
agedata<-aggregate(cbind(dead,total)~dose+sex+age,data=agedata,"sum")
agedata_f<-agedata[agedata$sex=="female",]
agedata_m<-agedata[agedata$sex=="male",]

#let's model the mortality rate for the females of the different classes of
#age
age_mod_f<-drm(dead/total~dose,weights=total,
               data=agedata_f,curveid=age,
               fct=LN.3u(),
               type="binomial")
#testing the goodness-of-fit of the model
modelFit(age_mod_f)
#comparing the LD50
compParm(age_mod_f,"e")
agerez_f<-ED(age_mod_f,50,interval="delta",reference="control")

#because it seems there is a problem to display the 95CI with models using 
#the 'curveid' parameter, we plot the different modalities separately
age_mod_f24<-drm(dead/total~dose,weights=total,
                 data=agedata_f[agedata_f$age=="0-24h",],
                 fct=LN.3u(),
                 type="binomial")
age_mod_f48<-drm(dead/total~dose,weights=total,
                 data=agedata_f[agedata_f$age=="24-48h",],
                 fct=LN.3u(),
                 type="binomial")
age_mod_f96<-drm(dead/total~dose,weights=total,
                 data=agedata_f[agedata_f$age=="72-96h",],
                 fct=LN.3u(),
                 type="binomial")

#let's model the mortality rate for the males of the different classes of
#age
age_mod_m<-drm(dead/total~dose,weights=total,
               data=agedata_m,curveid=age,
               fct=LN.3u(),
               type="binomial")
#testing the goodness-of-fit of the model
modelFit(age_mod_m)
#comparing the LD50
compParm(age_mod_m,"e")
agerez_m<-ED(age_mod_m,50,interval="delta",reference="control")

#comparison of LD50 allowing different slope and "natural death"
age_mod_me<-drm(dead/total~dose,weights=total,
               data=agedata_m,curveid=age,
               fct=LN.3u(),
               type="binomial",
               pmodels=list(~factor(age)-1, ~factor(age)-1, ~1))
plot(age_mod_me)
anova(age_mod_m,age_mod_me)

#because it seems there is a problem to display the 95CI with models using 
#the 'curveid' parameter, we plot the different modalities separately
age_mod_m24<-drm(dead/total~dose,weights=total,
                 data=agedata_m[agedata_m$age=="0-24h",],
                 fct=LN.3u(),
                 type="binomial")
age_mod_m48<-drm(dead/total~dose,weights=total,
                 data=agedata_m[agedata_m$age=="24-48h",],
                 fct=LN.3u(),
                 type="binomial")
age_mod_m96<-drm(dead/total~dose,weights=total,
                 data=agedata_m[agedata_m$age=="72-96h",],
                 fct=LN.3u(),
                 type="binomial")


##############################################################################/
#Code for plotting Figure 2####
##############################################################################/

#code for the plot comparing the different age categories of the same
#population, for male and female
op<-par(mar=c(0,5,6,1),mfrow=c(2,1))
#female plot for different age category
plot(age_mod_f48,type="confidence",col=rgb(0.4,0.2,0.6,1),
     bty="n",axes=FALSE,ann=FALSE,lwd=3)
plot(age_mod_f48,type="obs",add=TRUE,pch=21,cex=2,
     col=rgb(0.4,0.2,0.6,0.3),bg=rgb(0.4,0.2,0.6,0.3))
box(lwd=3,lty=1)
axis(1,at=c(1,10,100,500),labels=FALSE,
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("0","20","40","60","80","100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
plot(age_mod_f24,type="confidence",add=TRUE,
     col=rgb(0.6,0.2,0.2,1),lwd=3)
plot(age_mod_f24,type="obs",add=TRUE,pch=21,cex=2,
     col=rgb(0.6,0.2,0.2,0.3),bg=rgb(0.6,0.2,0.2,0.3))
plot(age_mod_f96,type="confidence",add=TRUE,
     col=rgb(0.4,0.5,0.4,1),lwd=3)
plot(age_mod_f96,type="obs",add=TRUE,pch=21,cex=2,
     col=rgb(0.4,0.5,0.4,0.3),bg=rgb(0.4,0.5,0.4,0.3))
text(1.5,y=0.85,labels='\\VE',vfont=c("sans serif","bold"),cex=5)
title(ylab="Mortality rate",cex.lab=2,font.lab=2)
#male plot for different age category
par(mar=c(5,5,1,1))
plot(age_mod_m48,type="confidence",col=rgb(0.4,0.2,0.6,1),
     bty="n",axes=FALSE,ann=FALSE,lwd=3,lty=2)
plot(age_mod_m48,type="obs",add=TRUE,pch=24,cex=2,
     col=rgb(0.4,0.2,0.6,1))
box(lwd=3,lty=1)
axis(1,at=c(1,10,100,500),labels=c("0","10","100","500"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("0","20","40","60","80","100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
plot(age_mod_m24,type="confidence",add=TRUE,
     col=rgb(0.6,0.2,0.2,1),lwd=3,lty=2)
plot(age_mod_m24,type="obs",add=TRUE,pch=24,cex=2,
     col=rgb(0.6,0.2,0.2,1))
plot(age_mod_m96,type="confidence",add=TRUE,
     col=rgb(0.4,0.5,0.4,1),lwd=3,lty=2)
plot(age_mod_m96,type="obs",add=TRUE,pch=24,cex=2,
     col=rgb(0.4,0.5,0.4,1))
title(xlab="Dose (mg/l)",ylab="Mortality rate",cex.lab=2,font.lab=2)
text(1.5,y=0.85,labels='\\MA',vfont=c("sans serif","bold"),cex=5)
par(op)
#export to pdf 10*14 inches


##############################################################################/
#END
##############################################################################/