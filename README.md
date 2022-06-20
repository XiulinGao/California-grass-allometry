# California-grass-allometry
This repo is for hosting data, scripts, and results (git ignored however) for allometry study of 12 common California grasses.
See each section for details on study design, data structure, and how to use each script for data analysis.

## Study design
This is a greenhouse experiment on grass allometry. Study species see Table 1. We grew all plants from seeds, all seeds are requested from the [USDA germplasm program](https://npgsweb.ars-grin.gov/gringlobal/search). Seeds were germinated on flats in Feb. 2022 with no stratification. Once there were at least 20 seedlings for each species with mean height greater than 5cm, we transplanted all secpeis into either nursery cones or 42cm (depth) * 15cm (diameter)pots. Plants in small cones were harvested in the first 2 months (March, 2022 - April, 2022) after transplant. Plants in deep pots were harvested between late April, 2022 and March, 2023. Nursery cones were used for early-stage, small-sized plants to allow for more samples across wide size range given limited greenhouse space. 


Table 1. Study species informaiton. 

  **Species**	       |  **USDA code**		  |  **Photosynthesis**   | **Life form**    |
  :------------------|:-------------------:|:---------------------:|:--------------- |
  *Avena barbata*	|AVBA	   				  |C3                     | annual           |
  *Bromus hordeaceus*|BRHO2               |C3                     | annual           |
  *Brachypodium distachyon*|BRDI2         |C3                     | annual           |
  *Vulpia myuros*   |VUMY                 |C3                     | annual           |
  *Elymus glaucus*  |ELGL                 |C3                     | perennial        |
  *Festuca californica*|FECA              |C3                     | perennial        |
  *Nassella pulchra*|NAPU4                |C3                     | perennial        |
  *Aristida oligantha*|AROL               |C4                     | annual           |
  *Setaria pumila*  |SEPU8                |C4                     | annual           |
  *Sporobolus airoides*|SPAI              |C4                     | perennial        |
  *Aristida purpurea*|ARPU9               |C4                     | perennial        |
  *Muhlenbergia rigens*|MURI2             |C4                     | perennial        | 
  
  												


A block design was applied. We set up 5 blocks in a glassroom at the Oxford facility at University of California, Berkeley. For each block, we assigned 7 (if total samples were 35) or 4 (if total samples were 20) individuals of the same species. Individuals within the same block were randomly located on the workbench. We randomized locations of all plants within the same block every other month. On each sample day, one plant from each sampled species was taken from each block. Given growth rate varied among species, we sampled all annual species in the spring and early summer of 2022, however, all perennial species were kept growing till March in the greenhouse to reach a reasonable size. Sampling dates and species were chosen to maximize the sampled size range of each species. 

All plants were watered frequently to maintain moist soil surface during wet season in California. However, watering frequency was reduced to every other month during June, 2022 - September, 2022 for mimicking dry season precipitation. This watering regime was chosen to drought stress the plants but also keep all plants alive and growing during summer so that they would reach a reasonal size at the end of the experiment. We apply fertilizer (100 ppm N 20-20-20 CaNO3) once each month to all plants. A 16-hour photoperiod was applied using medal halide. 


## Data structure
Data are stored in repo Dir/data. All files are organized in the same way. For instance, biomass.csv stores biomass measurement of different plant organs for each sample, while biomass-metadata.csv stores all the metadata information for biomass.csv
Relevant scripts used for data analysis are under repo Dir/scripts, all are written in R. 
the results dir meant to store figures and tables produced by scripts, but due to large size of the files, they are all git ignored here. 

## Usage of scripts 
inside each script, there is a brief description for what each script is used for and in-line comments. If there's any question, please contact the author for more information at xiulingao@lbl.gov

