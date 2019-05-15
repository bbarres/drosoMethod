[![DOI](https://zenodo.org/badge/119716271.svg)](https://zenodo.org/badge/latestdoi/119716271)
# How to create a reliable and reproducible insecticide resistance bioassay? An example on the worldwide invasive pest, *Drosophila suzukii*.
*This repository contains the R code used for the data analyses and production of the figures of the related article*


![alt text](https://j2ejmg.db.files.1drv.com/y4mfs0HpAp-0lm3RXzqAl_6ox6ANJQa-eeY3mIva0J6-lCC_iOKhirczqHbvFa1CbVb0zPHC62CYNYdRDSlUcYTQsepfEoC7Rmwm5mL_yKFWTqgLlbRiQ8RWuDxwEzTYUQqne5s6Sj7aI_ky82MSBhwN4rsbfdgoEmAVv7WUUCsUatxVesPePWVoVl-Sv0hMsYnAh5W2h4q5jLprGqbSMofWQ?width=1584&height=588&cropmode=none)


## Context


## Datasets
There are three datasets used in this study. 

The first dataset contains the general data for all the experiment conducted: 
+ droso_data.txt
  + *date*: the date of the experimentation
  + *sex*: the sex of the flies
  + *dose*: the dose of the active substance tested
  + *alive*: number of flies alive at the end of the test
  + *dead*: number of dead flies at the end of the test
  + *total*: the total number of flies
  + *age*: the age category of the flies in hours at the start of the test
  + *population*: population ID
  + *exposition*: duration of the exposition
  + *active substance*: active substance used for the test
  + *numb_comp*: experimentation used for assessing the effect of the number of flies per dose on the LD50 evaluation (0/1 = no/yes)
  + *genediv_comp*: experimentation used for assessing the effect of the genetic diversity of the population tested on the LD50 evaluation (0/1 = no/yes)
  + *expo_comp*: experimentation used for assessing the effect of the duration of exposure on the LD50 evaluation (0/1 = no/yes)
  + *age_comp*: experimentation used for assessing the effect of the age of the flies on the LD50 evaluation (0/1 = no/yes)

The second file contains a table used to plot the radarplot for the Heterozygote indice for each locus for the two populations
+ droso_GeneDiv.txt: there is one column for each of the 13 microsatellite markers used

The third file contains a table used to plot the radarplot for the Number of alleles for each locus for the two populations
+ droso_NumAll.txt: there is one column for each of the 13 microsatellite markers used

## R scripts
+ **droso_data_load.R:** the script to load the different datasets in the environment
+ **droso_ageeffect.R:** the script to analyse the effect of the age of the *D. suzukii* on the LD50 evaluation
+ **droso_effeffect:** the script to analyse the effect of the number of individuals of *D. suzukii* per dose tested on the LD50 evaluation
+ **droso_expoeffect.R:** the script to analyse the effect of the duration of exposure  of *D. suzukii* on the LD50 evaluation
+ **droso_geneteffect.R:** the script to analyse the effect of the level of genetic homogeneity of *D. suzukii* populations on the LD50 evaluation
+ **droso_sexeffect.R:** the script to analyse the effect of the sex of *D. suzukii*  on the LD50 evaluation

## Citation
You will soon be able (hopefully) to cite the related study as follow : 


If you want to use (some of) the code found on this page or if you want to cite this repository : 

