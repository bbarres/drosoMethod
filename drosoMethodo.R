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

#we do the same thing, but before that we combine the effective of the 
#different repetition
sexdata_conc<-aggregate(cbind(dead,total)~dose+sex,data=sexdata,"sum")
sex_mod2<-drm(dead/total~dose,weights=total,
              data=sexdata_conc,curveid=sex,
              fct=LN.3u(),
              type="binomial")
plot(sex_mod2,type="confidence")
plot(sex_mod2,type="obs",add=TRUE)
EDcomp(sex_mod2,c(50,50))
sexrez2<-ED(sex_mod2,50,interval="delta",reference="control")

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

EDcomp(list("sex_mod_f","sex_mod_m"),c(50,50))


###############################################################################
#What is the effect of the age of the flies on the LD50?
###############################################################################

#we select the data of phosmet test with the St Foy population, with different
#age classes
agedata<-dataDroz[dataDroz$age_comp==1,]
agedata_f<-agedata[agedata$sex=="female",]
agedata_m<-agedata[agedata$sex=="male",]


#let's do a model for every repetition mixing male and female
age_mod<-drm(dead/total~dose,weights=total,
             data=agedata,curveid=age,
             fct=LN.3u(),
             type="binomial")
plot(age_mod,type="confidence")
plot(age_mod,type="obs",add=TRUE)
EDcomp(age_mod,c(50,50))
agerez<-ED(age_mod,50,interval="delta",reference="control")

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


###############################################################################
#What is the effect of the genetic diversity of the tested population on LD50?
###############################################################################

#we select the data of phosmet test with the St Foy & SF IsoA populations
genDdata<-dataDroz[dataDroz$genediv_comp==1,]
genDdata_f<-genDdata[genDdata$sex=="female",]
genDdata_m<-genDdata[genDdata$sex=="male",]

#let's model the mortality rate for both populations
genD_mod<-drm(dead/total~dose,weights=total,
              data=genDdata,curveid=population,
              fct=LN.3u(),
              type="binomial")
plot(genD_mod,type="confidence")
plot(genD_mod,type="obs",add=TRUE)
EDcomp(genD_mod,c(50,50))
sexrez<-ED(genD_mod,50,interval="delta",reference="control")

#let's model the mortality rate for the females of both populations
genD_mod_f<-drm(dead/total~dose,weights=total,
                data=genDdata_f,curveid=population,
                fct=LN.3u(),
                type="binomial")
plot(genD_mod_f,type="confidence")
plot(genD_mod_f,type="obs",add=TRUE)
EDcomp(genD_mod_f,c(50,50))
sexrez_f<-ED(genD_mod_f,50,interval="delta",reference="control")

#let's model the mortality rate for the males of both populations
genD_mod_m<-drm(dead/total~dose,weights=total,
                data=genDdata_m,curveid=population,
                fct=LN.3u(),
                type="binomial")
plot(genD_mod_m,type="confidence",broken=TRUE)
plot(genD_mod_m,type="obs",add=TRUE)
EDcomp(genD_mod_m,c(50,50))
sexrez_m<-ED(genD_mod_m,50,interval="delta",reference="control")

#a combined graph of male and female regressions
op<-par(mfrow=c(2,1),mar=c(1,1,1,1))
plot(genD_f_mod,type="confidence")
plot(genD_f_mod,type="obs",add=TRUE)
abline(v=39.6,col="red")

plot(genD_m_mod,type="confidence")
plot(genD_m_mod,type="obs",add=TRUE)
abline(v=19.5,col="red")

par(op)

#another solution would be to do a logistic regression...


###############################################################################
#What is the effect of number of flies used for a test on LD50 evaluation?
###############################################################################

#we select the data of phosmet test with the St Foy population
numberdata<-dataDroz[dataDroz$number_comp==1,]

#let's do a model for every repetition
number_mod<-drm(dead/total~dose,weights=total,
                data=numberdata,curveid=repet,
                fct=LN.3u(),
                type="binomial")
plot(number_mod,type="confidence")
plot(number_mod,type="obs",add=TRUE)
numberrez<-ED(number_mod,50,interval="delta",reference="control")
write.table(numberrez,file="numberrez.txt",quote=FALSE,sep="\t",row.names=TRUE)


#another way to do the regression for each repetition is to use a loop######

#first we remove unnecessary levels in the data frame
dataDroz<-drop.levels(dataDroz)
#we then create a dataframe for the results
REZdroz<-data.frame("repet"=as.character(),"ED50"=as.numeric(),
                    "IC_low"=as.numeric(),"IC_up"=as.numeric())
#and here comes the loop
for (i in 1: length(levels(dataDroz$repet))) {
  temp.m1<-drm(dead/total~dose,weights=total,
               data=dataDroz[dataDroz$repet==levels(dataDroz$repet)[i],],
               fct=LN.3u(),
               type="binomial")
  temp<-ED(temp.m1,50,interval="delta",reference="control")
  tempx<-data.frame("date"=names(table(dataDroz$repet))[i],
                    "ED50"=temp[1],"IC_low"=temp[3],"IC_up"=temp[4])
  REZdroz<-rbind(REZdroz,tempx)
}

REZdroz<-REZdroz[order(as.character(REZdroz$date)),]
results<-cbind(checkdat,REZdroz)


#the scatter plot of the LD50 analysis with different number of fly per dose

plot(results$ED50[results$sexe=="femelle"]~results$total[results$sexe=="femelle"],
     xlab =" Number of tested D. suzukii adults",ylab="LD50 (mg/L)",
     main="LD50 values function of the number of tested females")
plotCI(results$total[results$sexe=="femelle"],
       results$ED50[results$sexe=="femelle"],
       ui=results$IC_up[results$sexe=="femelle"],
       li=results$IC_low[results$sexe=="femelle"],
       add=TRUE)
abline(39.5964,0,col="red",lwd=2)
abline(41.9867,0,col="red",lwd=2,lty=2)
abline(37.2061,0,col="red",lwd=2,lty=2)


plot(results$ED50[results$sexe=="male"]~results$total[results$sexe=="male"],
     xlab =" Number of tested D. suzukii adults",ylab="LD50 (mg/L)",
     main="LD50 values function of the number of tested males")
plotCI(results$total[results$sexe=="male"],
       results$ED50[results$sexe=="male"],
       ui=results$IC_up[results$sexe=="male"],
       li=results$IC_low[results$sexe=="male"],
       add=TRUE)
abline(19.1237,0,col="red",lwd=2)
abline(20.887,0,col="red",lwd=2,lty=2)
abline(17.36,0,col="red",lwd=2,lty=2)




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


