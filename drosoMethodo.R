###############################################################################
###############################################################################
#R code for the Drosophila pesticide resistance tests paper
###############################################################################
###############################################################################

#loading the libraries
library(drc)
library(plotrix)
library(gdata)

#set the working directory
setwd("~/work/Rfichiers/Githuber/droso_data")


###############################################################################
#What is the effect of the number of tested individuals on DL50 evaluation ?
###############################################################################

#load the dataset
dataDroz<-read.table("Droso_SF 24-48h phosmet 2016-2017.txt",header=T,sep="\t")
#creation of variable to distinguish between male and female
dataDroz<-cbind(dataDroz,"repet"=paste(dataDroz$date,dataDroz$sexe))

#let's sum the total number of individual tested and total number of 
#dead individual for each date
checkdat<-aggregate(cbind(dead,total)~date+sexe+repet,data=dataDroz,"sum")
checkdat<-checkdat[order(checkdat$repet),]
plot(checkdat)


#XXXXXhere put the general models for pooled male and female


#let's do a model for every repetition
droz_mod<-drm(dead/total~dose,weights=total,
              data=dataDroz,curveid=repet,
              fct=LN.3u(),
              type="binomial")
plot(droz_mod,type="confidence")
plot(droz_mod,type="obs",add=TRUE)
rez<-ED(droz_mod,50,interval="delta",reference="control")
write.table(rez,file="rez.txt",quote=FALSE,sep="\t",row.names=TRUE)


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












