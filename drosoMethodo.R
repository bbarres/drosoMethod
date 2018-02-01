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
plot(checkdat)

#let's do a model for every repetition
droz_mod<-drm(dead/total~dose,weights=total,
              data=dataDroz,curveid=repet,
              fct=LN.3u(),
              type="binomial")
plot(droz_mod,type="confidence")
plot(droz_mod,type="obs",add=TRUE)
rez<-ED(droz_mod,50,interval="delta",reference="control")
write.table(rez,file="rez.txt",quote=FALSE,sep="\t",row.names=TRUE)


####après là c'est pas encore terminé

#another way to do the regression for each repetition is to use a loop
#first we remove unnecessary levels in the data frame
dataDroz<-drop.levels(dataDroz)
#we then create a dataframe for the results
REZdroz<-data.frame("repet"=as.character(),"ED50"=as.numeric(),
                    "IC"=as.numeric())
#and here comes the loop
for (i in 1: length(levels(dataDroz$repet))) {
  temp.m1<-drm(dead/total~dose,weights=total,
               data=dataDroz[dataDroz$repet==levels(dataDroz$repet)[i],],
               fct=LN.3u(),
               type="binomial")
  temp<-ED(temp.m1,50,interval="delta",reference="control")
  tempx<-data.frame("date"=names(table(dataDroz$repet))[i],
                    "ED50"=temp[1],"IC"=temp[2])
  REZdroz<-rbind(REZdroz,tempx)
}

results<-cbind(checkdat,REZdroz)

plot(results$ED50[results$sexe=="femelle"]~results$total[results$sexe=="femelle"])
abline(36,0,col="red",lwd=2)
abline(38,0,col="red",lwd=2,lty=2)
abline(34,0,col="red",lwd=2,lty=2)

plot(results$ED50[results$sexe=="male"]~results$total[results$sexe=="male"])
abline(18,0,col="red",lwd=2)
abline(20.6,0,col="red",lwd=2,lty=2)
abline(15.4,0,col="red",lwd=2,lty=2)


REZbos$ED50[REZbos$ED50>30]<-30
plot(REZbos$ED50[order(REZbos$ED50)],main="Boscalid")
abline(0.39,0,col="green3",lwd=2)
abline(3.9,0,col="orange3",lwd=2)
write.table(REZbos,file="REZbos.txt",quote=FALSE,sep="\t",row.names=FALSE)

