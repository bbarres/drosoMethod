##############################################################################/
##############################################################################/
#Data and package loading for the analyses and figures plotting
##############################################################################/
##############################################################################/

#loading the packages necessary for the analysis
library(drc)
library(medrc)
library(plotrix)
library(gdata)
library(fmsb)
library(tidyr)
library(LDheatmap)
library(grid)
library(RColorBrewer)


##############################################################################/
#loading the bioassay dataset####
##############################################################################/

#load the dataset
dataDroz<-read.table("data/droso_data.txt",header=T,sep="\t",
                     stringsAsFactors=TRUE)
# #we remove the two concentrations that were used at the beginning of the 
# #test when we were still adjusting the range of doses for the test
# dataDroz<-dataDroz[dataDroz$dose!=603.70 & dataDroz$dose!=301.85,]
#creation of variable to distinguish between male and female and time 
#of exposure to pesticide
dataDroz<-cbind(dataDroz,"repet"=paste(dataDroz$date,dataDroz$sex, 
                                       dataDroz$exposition))


#we compute the total number of individual tested and total number of 	
#dead individual for each date	
checkdat<-aggregate(cbind(dead,total)~date+sex+repet,data=dataDroz,"sum")	
checkdat<-checkdat[order(checkdat$repet),]	
plot(checkdat)


##############################################################################/
#loading the genetic results data for the radarplot####
##############################################################################/

#load the datasets
dataNa<-read.table("data/droso_NumAll.txt",header=TRUE,sep="\t")
dataHS<-read.table("data/droso_GeneDiv.txt",header=TRUE,sep="\t")

#defining the two colors used for the two populations tested
colpoly_bord<-c(rgb(0.2,0.1,0.8,0.9),rgb(0.7,0.3,0.3,0.9))
colpoly_area<-c(rgb(0.2,0.1,0.8,0.4),rgb(0.7,0.3,0.3,0.4))


##############################################################################/
#Writing info session for reproducibility####
##############################################################################/

sink("session_info.txt")
print(sessioninfo::session_info())
sink()
#inspired by an R gist of François Briatte: 
#https://gist.github.com/briatte/14e47fb0cfb8801f25c889edea3fcd9b


##############################################################################/
#END
##############################################################################/