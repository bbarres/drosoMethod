###############################################################################
###############################################################################
#R code for the Drosophila suzukii pesticide resistance tests paper
###############################################################################
###############################################################################

#loading the libraries
library(fmsb)

#load the datasets
dataNa<-read.table("data/droso_NumAll.txt",header=TRUE,sep="\t")
dataHS<-read.table("data/droso_GeneDiv.txt",header=TRUE,sep="\t")

#defining the two colors used for the two populations tested
colpoly_bord<-c(rgb(0.2,0.1,0.8,0.9),rgb(0.7,0.3,0.3,0.9))
colpoly_area<-c(rgb(0.2,0.1,0.8,0.4),rgb(0.7,0.3,0.3,0.4))


###############################################################################
#radarplot for the comparison of the Number of alleles between two line
###############################################################################

radarchart(dataNa,
           #axis parameters
           axistype=1,caxislabels=seq(1,5,1),axislabcol="grey50",calcex=1.5,
           #grid parameters
           cglty=1,cglwd=2,cglcol="grey70",
           #labels parameters
           vlcex=1.5,
           #polygones parameters
           pcol=colpoly_bord,pfcol=colpoly_area,plwd=5,plty=1)
#export .pdf 7*7 inches


###############################################################################
#radarplot for the comparison of the Number of alleles between two line
###############################################################################

radarchart(dataHS,
           #axis parameters
           axistype=1,caxislabels=seq(0.00,1.00,0.25),axislabcol="grey50",
           calcex=1.5,
           #grid parameters
           cglty=1,cglwd=2,cglcol="grey70",
           #labels parameters
           vlcex=1.5,
           #polygones parameters
           pcol=colpoly_bord,pfcol=colpoly_area,plwd=5,plty=1)
#export .pdf 7*7 inches


###############################################################################
#Figure combining the two radarplot
###############################################################################

op<-par(mfrow=c(1,2),mar=c(2,1,2,0))
radarchart(dataNa,
           #axis parameters
           axistype=1,caxislabels=seq(1,5,1),axislabcol="grey50",calcex=1.5,
           #grid parameters
           cglty=1,cglwd=2,cglcol="grey70",
           #labels parameters
           vlcex=1.5,
           #polygones parameters
           pcol=colpoly_bord,pfcol=colpoly_area,plwd=5,plty=1)
radarchart(dataHS,
           #axis parameters
           axistype=1,caxislabels=seq(0.00,1.00,0.25),axislabcol="grey50",
           calcex=1.5,
           #grid parameters
           cglty=1,cglwd=2,cglcol="grey70",
           #labels parameters
           vlcex=1.5,
           #polygones parameters
           pcol=colpoly_bord,pfcol=colpoly_area,plwd=5,plty=1)
par(op)
#export .pdf 14*7 inches

###############################################################################
#END
###############################################################################