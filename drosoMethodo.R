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

#let's sum the total number of individual tested and total number of 
#dead individual for each date
checkdat<-aggregate(cbind(dead,total)~date,data=dataDroz,"sum")
plot(checkdat)

#let's do a model for every repetition
droz_mod<-drm(dead/total~dose,weights=total,
              data=dataDroz,curveid=date,
              fct=LN.2(),
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
REZdroz<-data.frame("date"=bosc_rez,"ED50"=30)
#and here comes the loop
for (i in 1: dim(table(bosc.dat$sample_ID))[1]) {
  temp.m1<-drm(perc_croiss~dose,
               data=bosc.dat[bosc.dat$sample_ID==names(table(bosc.dat$sample_ID))[i],],
               fct=LL.4())
  temp<-ED(temp.m1,50,type="absolute")
  tempx<-data.frame("sample_ID"=names(table(bosc.dat$sample_ID))[i],
                    "ED50"=temp[1])
  REZbos<-rbind(REZbos,tempx)
}

REZbos$ED50[REZbos$ED50>30]<-30
plot(REZbos$ED50[order(REZbos$ED50)],main="Boscalid")
abline(0.39,0,col="green3",lwd=2)
abline(3.9,0,col="orange3",lwd=2)
write.table(REZbos,file="REZbos.txt",quote=FALSE,sep="\t",row.names=FALSE)

