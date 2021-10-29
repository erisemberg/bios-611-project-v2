# bios-611-project-v2
This is the version of my BIOS 611 project that uses my own thesis data. The shiny HW was done on a different set of data (see github.com/erisemberg/bios-611-project). 

Genetic analyses to determine loci associated with allergic reaction to peanut
==============================================================================

This data science project performs genetic association analysis on a data set including genotypes and phenotypes from a backcross between two inbred mouse strains: CC27, which has peanut allergy, and C3H, which does not have peanut allergy. By studying their F2 projeny, we may gain insight into which genetic loci are associated with variation in peanut allergy phenotype. We do this with quantitative trait loci (QTL) analysis.

Phenotype data includes temperature post-peanut exposure and symptom score. 

To run the code associated with this project, build the docker container like this: 

```
docker build . -t pnut 
```

First produce the spreadsheet for analysis from the raw genotype and phenotype data. This code filters the genetic markers to those that are informative and then integrates the genetic data with the phenotype data and produces a file in the correct format for analysis with the `R/qtl` package. 
```
make clean 
make derived_data/Rqtl_CC27xC3H_BC.csv
```

Then, perform the QTL analysis and produce the figures. This code currently calculates a summary statistic for the trajectory of temperature for each mouse ("area above the curve"). Then we perform QTL analysis, using that summary statistic as the phenotype. 

To produce a genetic marker map:
```
make figures/genetic_map.png
```

To produce figures summarizing temperature data:
```
make figures/temp_histograms.png
make figures/temp_time-series.png
make figures/temp_aac_hist.png
make figures/temp_aac_transformed_hist.png
```

To produce a genome scan, which is a plot showing the results of the QTL analysis (where peaks represent significant associations between genotype and phenotype):
```
make figures/temp_aac_genome_scan.png
``` 

Logs for each step will be saved to the logs/ folder. 

To start an RStudio server at any point, run: 
```
docker run -p 8787:8787 -e PASSWORD=pw123 -v $(pwd):/home/rstudio/project -t pnut 
```