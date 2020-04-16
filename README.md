[![DOI](https://zenodo.org/badge/119716271.svg)](https://zenodo.org/badge/latestdoi/119716271)
# Supporting data and code for: How to create a reliable and reproducible insecticide resistance bioassay? An example on the worldwide invasive pest, *Drosophila suzukii*.
*This repository contains the R code used for the data analyses and production of the figures of the related article*

![alt text](https://j2ejmg.db.files.1drv.com/y4mfs0HpAp-0lm3RXzqAl_6ox6ANJQa-eeY3mIva0J6-lCC_iOKhirczqHbvFa1CbVb0zPHC62CYNYdRDSlUcYTQsepfEoC7Rmwm5mL_yKFWTqgLlbRiQ8RWuDxwEzTYUQqne5s6Sj7aI_ky82MSBhwN4rsbfdgoEmAVv7WUUCsUatxVesPePWVoVl-Sv0hMsYnAh5W2h4q5jLprGqbSMofWQ?width=1584&height=588&cropmode=none)


## Context
Bioassays are the golden standard for pesticide resistance assessment of pest. Developing or adapting a method for a new species or a new active substance can require a significant amount of work. Several parameters that can impact the reliability and precision of the bioassay shouldn't be overlooked during the setting of a new protocol. In this study, we exemplify the effect of three biological parameters (sex, age, and genetic diversity) and two technical parameters (number of individuals tested and duration of exposure) on the median lethal dose (LD50) evaluation of an invasive pest species, *Drosophila suzukii*. 


## Datasets
There are three datasets used in this study. The files can be found in the "data" folder. 

+ **droso_data.txt:** the first dataset contains the general data for all the experiment conducted
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

+ **droso_GeneDiv.txt:** the second file contains a table used to plot the radarplot for the Heterozygote indice for each locus for the two populations. There is one column for each of the 13 microsatellite markers used. 

+ **droso_NumAll.txt:** the third file contains a table used to plot the radarplot for the Number of alleles for each locus for the two populations. There is one column for each of the 13 microsatellite markers used


## R scripts
+ **droso_data_load.R:** the script to load the different datasets in the environment
+ **droso_ageeffect.R:** the script to analyse the effect of the age of the *D. suzukii* on the LD50 evaluation
+ **droso_effeffect:** the script to analyse the effect of the number of individuals of *D. suzukii* per dose tested on the LD50 evaluation
+ **droso_expoeffect.R:** the script to analyse the effect of the duration of exposure  of *D. suzukii* on the LD50 evaluation
+ **droso_geneteffect.R:** the script to analyse the effect of the level of genetic homogeneity of *D. suzukii* populations on the LD50 evaluation
+ **droso_sexeffect.R:** the script to analyse the effect of the sex of *D. suzukii*  on the LD50 evaluation


## Citation
You will soon be able (hopefully) to cite the related study as follow: 
+ Blouquy L., Mottet C., Olivares J., Plantamp C., Siegwart M. and Barrès B.
[How to create a reliable and reproducible insecticide resistance bioassay? An example on the worldwide invasive pest, *Drosophila suzukii*. *Submitted*]()

If you want to use (some of) the code found on this page or if you want to cite this repository:
+ Benoit Barrès. bbarres/drosoMethod: [Supporting data and code for: How to create a reliable and reproducible insecticide resistance bioassay? An example on the worldwide invasive pest, *Drosophila suzukii*. Zenodo; 2019.](https://zenodo.org/badge/latestdoi/119716271)
