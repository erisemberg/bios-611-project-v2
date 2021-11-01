# bios-611-project-v2
This is the version of my BIOS 611 project that uses my own thesis data. The shiny HW was done on a different set of data (see github.com/erisemberg/bios-611-project). 

Genetic analyses to determine loci associated with allergic reaction to peanut
==============================================================================

This data science project performs genetic association analysis to identify loci associated with peanut allergy in mice. 

To run the code associated with this project, build the docker container like this: 

```
docker build . -t pnut 
```

To run the docker container and open a terminal session within the container:

```
docker run -e PASSWORD=pw123 --rm -v $(pwd):/home/rstudio/work -p 8787:8787 -it pnut /bin/bash
```

Navigate to the working directory:

```
cd home/rstudio/work
```

First produce the spreadsheet for analysis from the raw genotype and phenotype data. This code filters the genetic markers to those that are informative and then integrates the genetic data with the phenotype data and produces a file in the correct format for analysis with the `R/qtl` package. 

```
make clean 
make derived_data/Rqtl_CC27xC3H_BC.csv
```

Then, perform the QTL analysis and produce the figures. This code currently calculates a summary statistic for the trajectory of temperature for each mouse ("area above the curve"). Then we perform QTL analysis, using that summary statistic as the phenotype. 

Note that running `make` on any of the following figures will produce all of the figures, so you only need to run one of them. 

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

To create the final report:
```
make report.pdf
```

Logs for each step will be saved to the `logs/` folder. 

To start an RStudio server at any point, go to `http://localhost:8787/` and login with user `rstudio` and password `pw123`. 