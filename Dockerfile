FROM rocker/verse
RUN R -e "install.packages(\"tidyverse\")"
RUN R -e "install.packages(\"readxl\")"
RUN R -e "install.packages(\"qtl\")" 
RUN R -e "install.packages(\"MASS\")" 
RUN R -e "install.packages(\"MESS\")"
RUN R -e "install.packages(\"lme4\")"