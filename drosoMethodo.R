###############################################################################
###############################################################################
#R code for the Drosophila suzukii pesticide resistance tests paper
###############################################################################
###############################################################################

#loading the libraries
library(drc)
library(plotrix)
library(gdata)

#set the working directory
setwd("~/work/Rfichiers/Githuber/droso_data")
setwd("K:/Projets de recherche/2015-PPV-résistance/ARTICLE-Dsuzukii")

#load the dataset
dataDroz<-read.table("droso_data.txt",header=T,sep="\t")
# #we remove the two concentrations that were used at the beginning of the 
# #test when we were still adjusting the range of doses for the test
# dataDroz<-dataDroz[dataDroz$dose!=603.70 & dataDroz$dose!=301.85,]
#creation of variable to distinguish between male and female and time 
#of exposure to pesticide
dataDroz<-cbind(dataDroz,"repet"=paste(dataDroz$date,dataDroz$sex, 
                                       dataDroz$exposition))


#let's sum the total number of individual tested and total number of 
#dead individual for each date
checkdat<-aggregate(cbind(dead,total)~date+sex+repet,data=dataDroz,"sum")
checkdat<-checkdat[order(checkdat$repet),]
plot(checkdat)


###############################################################################
#Is there a difference of DL50 between male and female?
###############################################################################

#we select the data of phosmet tests with the St Foy population, with
#flies of 24-48 hours age
sexdata<-dataDroz[dataDroz$number_comp==1,]
sexdata_f<-sexdata[sexdata$sex=="female",]
sexdata_m<-sexdata[sexdata$sex=="male",]

#let's do a model for every repetition
sex_mod<-drm(dead/total~dose,weights=total,
             data=sexdata,curveid=sex,
             fct=LN.3u(),
             type="binomial")
plot(sex_mod,type="confidence")
plot(sex_mod,type="obs",add=TRUE)
EDcomp(sex_mod,c(50,50))
sexrez<-ED(sex_mod,50,interval="delta",reference="control")

#because there is a bug to display the 95CI with models using curveid,
#we plot the different modality separately
sex_mod_f<-drm(dead/total~dose,weights=total,
               data=sexdata_f,
               fct=LN.3u(),
               type="binomial")
plot(sex_mod_f,type="confidence")
plot(sex_mod_f,type="obs",add=TRUE)
ED(sex_mod_f,50)

sex_mod_m<-drm(dead/total~dose,weights=total,
               data=sexdata_m,
               fct=LN.3u(),
               type="binomial")
plot(sex_mod_m,type="confidence")
plot(sex_mod_m,type="obs",add=TRUE)
ED(sex_mod_m,50)

op<-par(mar=c(5,5,4,1))
plot(sex_mod_f,type="confidence",col="black",bty="n",axes=FALSE,ann=FALSE,
     lwd=3)
plot(sex_mod_f,type="confidence",col="black",add=TRUE)
box(lwd=3,lty=1)
axis(1,at=c(1,10,100,500),labels=c("0","10","100","500"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("0","20","40","60","80","100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
plot(sex_mod_f,type="obs",add=TRUE,pch=21,cex=2,
     col=rgb(0,0,0,0.5),bg=rgb(0,0,0,0.5))
plot(sex_mod_m,type="confidence",add=TRUE,col="grey40",lty=2,lwd=3)
plot(sex_mod_m,type="obs",add=TRUE,pch=24,cex=2,
     col=rgb(0,0,0,0.3),bg=rgb(0,0,0,0.0))
title(xlab="Dose (mg/L)",ylab="Mortality rate",cex.lab=2,font.lab=2)
par(op)
#export .pdf 10*7 inches


###############################################################################
#What is the effect of the age of the flies on their LD50?
###############################################################################

#we select the data of phosmet test with the St Foy population, with different
#age classes
agedata<-dataDroz[dataDroz$age_comp==1,]
agedata_f<-agedata[agedata$sex=="female",]
agedata_m<-agedata[agedata$sex=="male",]

#let's model the mortality rate for the females of the different classes of
#age
age_mod_f<-drm(dead/total~dose,weights=total,
               data=agedata_f,curveid=age,
               fct=LN.3u(),
               type="binomial")
plot(age_mod_f,type="confidence")
plot(age_mod_f,type="obs",add=TRUE)
EDcomp(age_mod_f,c(50,50))
agerez_f<-ED(age_mod_f,50,interval="delta",reference="control")

#because there is a bug to display the 95CI with models using curveid,
#we plot the different modality separately
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
plot(age_mod_m,type="confidence")
plot(age_mod_m,type="obs",add=TRUE)
EDcomp(age_mod_m,c(50,50))
agerez_m<-ED(age_mod_m,50,interval="delta",reference="control")

#because there is a bug to display the 95CI with models using curveid,
#we plot the different modality separately
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
title(xlab="Dose (mg/L)",ylab="Mortality rate",cex.lab=2,font.lab=2)
text(1.5,y=0.85,labels='\\MA',vfont=c("sans serif","bold"),cex=5)
par(op)
#export .pdf 10*14 inches

#in order to take into account both the gender and the age of the population
#at the same time, we perform a logistic regression. 
#we set the age category "72-96h" as the reference
agedata$age<-relevel(agedata$age,ref="72-96h")
LogReg_gene<-glm(cbind(agedata$alive,agedata$dead)~dose+age+sex,
                 family=binomial(link=probit),data=agedata)
summary(LogReg_gene)


###############################################################################
#What is the effect of the genetic diversity of the tested population on LD50?
###############################################################################

#we select the data of phosmet test with the St Foy & SF IsoA populations
genDdata<-dataDroz[dataDroz$genediv_comp==1,]
genDdata_f<-genDdata[genDdata$sex=="female",]
genDdata_m<-genDdata[genDdata$sex=="male",]

#let's model the mortality rate for the females of both populations
genD_mod_f<-drm(dead/total~dose,weights=total,
                data=genDdata_f,curveid=population,
                fct=LN.2(),
                type="binomial")
plot(genD_mod_f,type="confidence")
plot(genD_mod_f,type="obs",add=TRUE)
EDcomp(genD_mod_f,c(50,50))
sexrez_f<-ED(genD_mod_f,50,interval="delta",reference="control")

#because there is a bug to display the 95CI with models using curveid,
#we plot the different modality separately
genD_mod_f_stf<-drm(dead/total~dose,weights=total,
                    data=genDdata_f[genDdata_f$population=="ste-foy",],
                    fct=LN.2(),
                    type="binomial")
genD_mod_f_isa<-drm(dead/total~dose,weights=total,
                    data=genDdata_f[genDdata_f$population=="sf-isoa",],
                    fct=LN.2(),
                    type="binomial")

#let's model the mortality rate for the males of both populations
genD_mod_m<-drm(dead/total~dose,weights=total,
                data=genDdata_m,curveid=population,
                fct=LN.2(),
                type="binomial")
plot(genD_mod_m,type="confidence")
plot(genD_mod_m,type="obs",add=TRUE)
EDcomp(genD_mod_m,c(50,50))
sexrez_m<-ED(genD_mod_m,50,interval="delta",reference="control")

#because there is a bug to display the 95CI with models using curveid,
#we plot the different modality separately
genD_mod_m_stf<-drm(dead/total~dose,weights=total,
                    data=genDdata_m[genDdata_m$population=="ste-foy",],
                    fct=LN.2(),
                    type="binomial")
genD_mod_m_isa<-drm(dead/total~dose,weights=total,
                    data=genDdata_m[genDdata_m$population=="sf-isoa",],
                    fct=LN.2(),
                    type="binomial")

#code for the plot comparing the different populations with different
#levels of genetic diversity, for male and female
op<-par(mar=c(0,5,6,1),mfrow=c(2,1))
#female plot for different genetic diversity populations
plot(genD_mod_f_stf,type="confidence",col=rgb(0.4,0.2,0.6,1),
     bty="n",axes=FALSE,ann=FALSE,lwd=3)
plot(genD_mod_f_stf,type="obs",add=TRUE,pch=21,cex=2,
     col=rgb(0.4,0.2,0.6,0.3),bg=rgb(0.4,0.2,0.6,0.3))
box(lwd=3,lty=1)
axis(1,at=c(1,10,50,150),labels=FALSE,
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("0","20","40","60","80","100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
segments(ED(genD_mod_f_stf,50,interval="delta",reference="control")[1],-0.2,
         ED(genD_mod_f_stf,50,interval="delta",reference="control")[1],0.5,
         lwd=3,col=rgb(0.4,0.2,0.6,1),lty=1)
plot(genD_mod_f_isa,type="confidence",add=TRUE,
     col=rgb(0.6,0.2,0.2,1),lwd=3)
plot(genD_mod_f_isa,type="obs",add=TRUE,pch=21,cex=2,
     col=rgb(0.6,0.2,0.2,0.3),bg=rgb(0.6,0.2,0.2,0.3))
segments(ED(genD_mod_f_isa,50,interval="delta",reference="control")[1],-0.2,
         ED(genD_mod_f_isa,50,interval="delta",reference="control")[1],0.5,
         lwd=3,col=rgb(0.6,0.2,0.2,1),lty=1)
abline(h=0.5,lwd=3,col=grey(0.5))
text(1.5,y=0.85,labels='\\VE',vfont=c("sans serif","bold"),cex=5)
title(ylab="Mortality rate",cex.lab=2,font.lab=2)
#male plot for different genetic diversity populations
par(mar=c(5,5,1,1))
plot(genD_mod_m_stf,type="confidence",col=rgb(0.4,0.2,0.6,1),
     bty="n",axes=FALSE,ann=FALSE,lwd=3,lty=2)
plot(genD_mod_m_stf,type="obs",add=TRUE,pch=24,cex=2,
     col=rgb(0.4,0.2,0.6,1))
box(lwd=3,lty=1)
axis(1,at=c(1,10,50,150),labels=c("0","10","50","150"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0,0.2,0.4,0.6,0.8,1),labels=c("0","20","40","60","80","100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
segments(ED(genD_mod_m_stf,50,interval="delta",reference="control")[1],-0.2,
         ED(genD_mod_m_stf,50,interval="delta",reference="control")[1],0.5,
         lwd=3,col=rgb(0.4,0.2,0.6,1),lty=2)
plot(genD_mod_m_isa,type="confidence",add=TRUE,
     col=rgb(0.6,0.2,0.2,1),lwd=3,lty=2)
plot(genD_mod_m_isa,type="obs",add=TRUE,pch=24,cex=2,
     col=rgb(0.6,0.2,0.2,1))
segments(ED(genD_mod_m_isa,50,interval="delta",reference="control")[1],-0.2,
         ED(genD_mod_m_isa,50,interval="delta",reference="control")[1],0.5,
         lwd=3,col=rgb(0.6,0.2,0.2,1),lty=2)
abline(h=0.5,lwd=3,col=grey(0.5))
title(xlab="Dose (mg/L)",ylab="Mortality rate",cex.lab=2,font.lab=2)
text(1.5,y=0.85,labels='\\MA',vfont=c("sans serif","bold"),cex=5)
par(op)
#export .pdf 8*14 inches

#in order to take into account both the gender and the population at the 
#same time, we performed a logistic regression. 
LogReg_gene<-glm(cbind(genDdata$alive,genDdata$dead)~dose+population+sex,
                 family=binomial(link=probit),data=genDdata)
summary(LogReg_gene)


###############################################################################
#What is the effect of number of flies used for a test on LD50 evaluation?
###############################################################################

#we select the data of phosmet test with the St Foy population
numberdata<-dataDroz[dataDroz$number_comp==1,]
#because on the 2016-06-23, the behaviour of the flies was strange, 
#this repetition was removed before analysis
numberdata<-numberdata[numberdata$date!="2016-06-23",]

#first we remove unnecessary levels in the data frame
numberdata<-drop.levels(numberdata)
#we then create a dataframe for the results
REZdroz<-data.frame("repet"=as.character(),"ED50"=as.numeric(),
                    "IC_low"=as.numeric(),"IC_up"=as.numeric(),
                    "SE"=as.numeric())
#and here comes the loop
for (i in 1: length(levels(numberdata$repet))) {
  temp.m1<-drm(dead/total~dose,weights=total,
               data=numberdata[numberdata$repet==levels(numberdata$repet)[i],],
               fct=LN.2(),
               type="binomial")
  temp<-ED(temp.m1,50,interval="delta",reference="control")
  tempx<-data.frame("date"=names(table(numberdata$repet))[i],
                    "ED50"=temp[1],"IC_low"=temp[3],"IC_up"=temp[4],
                    "SE"=temp[2])
  REZdroz<-rbind(REZdroz,tempx)
}
REZdroz<-REZdroz[order(as.character(REZdroz$date)),]
results<- merge(checkdat,REZdroz,by.x="repet",by.y="date",all=FALSE)

#the scatter plot of the LD50 analysis with different number of fly per dose
op<-par(mar=c(5,5,1,1),mfrow=c(1,2))
plot(results$ED50[results$sex=="female"]~results$total[results$sex=="female"],
     xlab ="Mean number of D. suzukii per dose",ylab="LD50 (mg/L)",
     ylim=c(0,100),xlim=c(50,270),bty="n",ann=FALSE,axes=FALSE,
     cex=2,pch=21,col=rgb(0,0,0,0.0),bg=rgb(0,0,0,0.0))
box(lwd=3,lty=1)
axis(1,at=c(70,105,140,175,210,245),labels=c("10","15","20","25","30","35"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0,20,40,60,80,100),labels=c("0","20","40","60","80","100"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
title(xlab="Mean number of D. suzukii per dose",
      ylab="LD50 (mg/L)",cex.lab=2,font.lab=2)
text(230,y=90,labels='\\VE',vfont=c("sans serif","bold"),cex=5)
plotCI(results$total[results$sex=="female"],
       results$ED50[results$sex=="female"],
       ui=results$ED50[results$sex=="female"]+
         results$SE[results$sex=="female"],
       li=results$ED50[results$sex=="female"]-
         results$SE[results$sex=="female"],
       #ui=results$IC_up[results$sex=="female"],
       #li=results$IC_low[results$sex=="female"],
       add=TRUE,cex=2,pch=21,col=rgb(0,0,0,1),pt.bg=rgb(0,0,0,0.3),
       gap=0.01)
totfem<-drm(dead/total~dose,weights=total,
            data=numberdata[numberdata$sex=="female",],
            fct=LN.2(),
            type="binomial")
totfemREZ<-ED(totfem,50,interval="delta",reference="control")
abline(totfemREZ[1],0,col="red",lwd=2)
abline(totfemREZ[3],0,col="red",lwd=2,lty=2)
abline(totfemREZ[4],0,col="red",lwd=2,lty=2)

par(mar=c(5,0,1,6))
plot(results$ED50[results$sex=="male"]~results$total[results$sex=="male"],
     xlab ="Mean number of D. suzukii per dose",ylab="LD50 (mg/L)",
     ylim=c(0,100),xlim=c(50,270),bty="n",ann=FALSE,axes=FALSE,
     cex=2,pch=24,col=rgb(0,0,0,0.0),bg=rgb(0,0,0,0.0))
box(lwd=3,lty=1)
axis(1,at=c(70,105,140,175,210,245),labels=c("10","15","20","25","30","35"),
     cex.axis=1.5,font.axis=2,lwd.ticks=2)
axis(2,at=c(0,20,40,60,80,100),labels=FALSE,
     cex.axis=1.5,font.axis=2,lwd.ticks=2,las=1)
title(xlab="Mean number of D. suzukii per dose",
      ylab="",cex.lab=2,font.lab=2)
text(230,y=90,labels='\\MA',vfont=c("sans serif","bold"),cex=5)
plotCI(results$total[results$sex=="male"],
       results$ED50[results$sex=="male"],
       ui=results$ED50[results$sex=="male"]+
         results$SE[results$sex=="male"],
       li=results$ED50[results$sex=="male"]-
         results$SE[results$sex=="male"],
       #ui=results$IC_up[results$sex=="male"],
       #li=results$IC_low[results$sex=="male"],
       add=TRUE,cex=2,pch=24,col=rgb(0,0,0,1),pt.bg=rgb(0,0,0,0.0),
       gap=0.01)
totmal<-drm(dead/total~dose,weights=total,
            data=numberdata[numberdata$sex=="male",],
            fct=LN.2(),
            type="binomial")
totmalREZ<-ED(totmal,50,interval="delta",reference="control")
abline(totmalREZ[1],0,col="red",lwd=2)
abline(totmalREZ[3],0,col="red",lwd=2,lty=2)
abline(totmalREZ[4],0,col="red",lwd=2,lty=2)
par(op)
#export .pdf 15*7 inches


###############################################################################
#What is the effect of exposure time on the evaluation of LD50?
###############################################################################

#we select the data of lambda-cyhalothrin test with the St Foy population
expodata<-dataDroz[dataDroz$expo_comp==1,]

#let's do a model for every repetition
expo_mod<-drm(dead/total~dose,weights=total,
              data=expodata,curveid=repet,
              fct=LN.3u(),
              type="binomial")
plot(expo_mod,type="confidence")
plot(expo_mod,type="obs",add=TRUE)
EDcomp(expo_mod,c(50,50))
exporez<-ED(expo_mod,50,interval="delta",reference="control")


