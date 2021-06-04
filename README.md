# Age-structured Jolly-Seber model expands inference and improves parameter estimation from capture-recapture data

### Hostetter NJ, NJ Lunn, ES Richardson, EV Regehr, and SJ Converse

### PLOS ONE

### DOI: *Accepted*

### Please contact the first author for questions about the code or data: Nathan Hostetter (njhostet@ncsu.edu)
__________________________________________________________________________________________________________________________________________

## Citation:
Hostetter NJ, NJ Lunn, ES Richardson, EV Regehr, and SJ Converse. Age-structured Jolly-Seber model expands inference and improves parameter estimation from capture-recapture data. PLOS One. *Accepted*

## Abstract:
Understanding the influence of individual attributes on demographic processes is a key objective of wildlife population studies. Capture-recapture and age data are commonly collected to investigate hypotheses about survival, reproduction, and viability. We present a novel age-structured Jolly-Seber model that incorporates age and capture-recapture data to provide comprehensive information on population dynamics, including abundance, age-dependent survival, recruitment, age structure, and population growth rates. We applied our model to a multi-year capture-recapture study of polar bears (Ursus maritimus) in western Hudson Bay, Canada (2012–2018), where management and conservation require a detailed understanding of how polar bears respond to climate change and other factors. In simulation studies, the age-structured Jolly-Seber model improved precision of survival, recruitment, and annual abundance estimates relative to standard Jolly-Seber models that omit age information. Furthermore, incorporating age information improved precision of population growth rates, increased power to detect trends in abundance, and allowed direct estimation of age-dependent survival and changes in annual age structure. Our case study provided detailed evidence for senescence in polar bear survival. Median survival estimates were lower (<0.95) for individuals aged <5 years, remained high (>0.95) for individuals aged 7–22 years, and subsequently declined to near zero for individuals >30 years. We also detected cascading effects of large recruitment classes on population age structure, which created major shifts in age structure when these classes entered the population and then again when they reached prime breeding ages (10–15 years old). Overall, age-structured Jolly-Seber models provide a flexible means to investigate ecological and evolutionary processes that shape populations (e.g., via senescence, life expectancy, and lifetime reproductive success) while improving our ability to investigate population dynamics and forecast population changes from capture-recapture data.


## Code 
1) [SupplementS2_SimulationAgeStructuredJS.R](./simulation/): This folder contains the R script to generate and analyze age-structured mark-recapture data described in the manuscript.
2) [SupplementS4_CaseStudy.R](./CaseStudy/): This folder contains the R script to load, format, and run the models described in the case study.


## Data
Datasets used in this project are all found in the [data folder](./data):

1) [SupplementS3_CaseStudyData.txt](./data/): Formatted data to run the age-structured Jolly-Seber model. See manuscript for detailed description of data.
