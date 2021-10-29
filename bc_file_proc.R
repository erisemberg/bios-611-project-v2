# This script processes genotype and phenotype data for R/qtl analysis. Performs
# marker filtering and recoding (ATCG -> AA/AB/BB) and integrates phenotype data 
# into a csv file in the appropriate format for analysis with R/qtl. 
#
# Input:
#   Genotype file - first three columns must be marker name, chr, pos; followed 
#       any number of marker metadata columns; followed by parent/F1 genotypes; 
#       followed by F2 genotypes 
#       ### F2 genotypes must be the last columns 
#       ### Change Position column name to "Position_b38"
#       Example (SARS data): 1-marker name; 2-chr; 3-pos; 4-ref; 5-alt; 6-cc_reliable; 
#       7-cc_sdp; 8 to 23-parent genotypes; 24 to 27-F1 genotypes (female); 
#       28 to 29-F1 genotypes (male); 30 to 819: genotype data
#   Phenotype file - each row is an F2 mouse, each column is a phenotype. Must 
#       have a "Geno_ID" column in the format "Cr_RB05_<sex>_<padded UNC ID>"
#
# Output: csv file in R/qtl format

library(tidyverse)
library(readxl)
source("utils.R")

ensure_directory("logs")
ensure_directory("derived_data")
log <- make_logger("logs/file_processing_notes.md")

#---------------------------------Parameters-----------------------------------#
geno_file <- "source_data/Cr_WB02_miniMUGA-06242021.csv"
pheno_file <- "source_data/CC27xC3H_backcross_pheno_clean.xlsx"
out_file <- "derived_data/Rqtl_CC27xC3H_BC.csv"
Y_MT_out_file <- "derived_data/Geno_CC27xC3H_Y_MT.csv"

num_F2s = 365 # number of genotyped F2 mice
pheno.names = c("Cage", "Batch", "Geno_ID", "Temp_0min", "Temp_15min", "Temp_30min",
                "Temp_45min", "Temp_60min", "Delta_15min", "Delta_30min", "Delta_45min",
                "Delta_60min", "Min_Temp", "Min_Temp_Time", "Symptom_score",
                "Symptom_bin", "Diarrhea_bin", "Diarrhea_time") # Phenotypes to include
# infection = "SARS-CoV" # infection type
A.ref = "CC027.GeniUnc" # parent mouse representing A genotype 
B.ref = "C3H.HeJ" # parent mouse representing B genotype 

#------------------------------------------------------------------------------#
#-------------------------------Genotype data----------------------------------#
#------------------------------------------------------------------------------#
geno <- read.csv(geno_file)
log("Genotype data downloaded...")

F2_start = ncol(geno)-num_F2s+1 # column index of first F2 mouse column
F2_end = ncol(geno) # column index of last F2 mouse column 

geno <- geno[geno$Chromosome != 0, ] # Remove chr0 rows
geno$Position <- geno$Position/1e6 # convert bp position to cM

#----------------------------Metric calculation--------------------------------#
# Get ref, alt, het and N calls in the F2 mice at each marker
geno$num_ref <- rowSums(geno[,F2_start:F2_end] == geno$reference)
geno$num_alt <- rowSums(geno[,F2_start:F2_end] == geno$alternate)
geno$num_N <- rowSums(geno[,F2_start:F2_end] == "N")
geno$pct_N <- geno$num_N/num_F2s
geno$num_H <- rowSums(geno[,F2_start:F2_end] == "H")
# Calculate ref relative to alt - ref/(ref+alt) 
geno$ref_alt <- geno$num_ref/(geno$num_ref + geno$num_alt)
# calculate het relative to good calls - het/(ref+alt+het)
geno$het_all <- geno$num_H/(geno$num_ref + geno$num_alt + geno$num_H)

#--------------------------------Filtering-------------------------------------#
### Separate out MT and Y markers // what to do with PAR? 
geno.Y.MT <- geno[which((geno$Chromosome == "Y") | geno$Chromosome == "MT"), ]
### Filter Y chromosome markers - all females should have N calls for each Y marker
#geno.Y.MT <- geno[which((geno$Chromosome == "Y") & (geno$pct_N == 0))]
geno <- geno[which((geno$Chromosome != "Y") & (geno$Chromosome != "MT")),]

# Filter out markers with failure rate > 5% N  
geno <- geno[which(geno$pct_N <= 0.05),]

### Pre-filtering het_all and ref_alt - bad and uninformative markers
### will filter more on these metrics later 

# Remove markers where alt was only call (>99% alt) 
geno <- geno[which((geno$num_alt/num_F2s) <= 0.99),]
# Remove markers where ref was only call (>99% ref) 
geno <- geno[which((geno$num_ref/num_F2s) <= 0.99),]
# Remove markers where het was only call (>97% het?) 
geno <- geno[which((geno$num_H/num_F2s) <= 0.97),]

# Plot x = ref/(ref+alt) and y = het/(ref+alt+het)
png("figures/marker_qc_plot.png")
plot(x=geno$het_all, y=geno$ref_alt, main="Autosome/X markers", 
     ylab="Ref/(Ref + Alt)", xlab="Het/(Het + Ref + Alt)")
dev.off()

# Looking for het:hom ratio to be ~0.5 (backcross, so all should be AA/AB)
# Remove anything where het_all is <0.4 or >0.6   
geno <- geno[which((geno$het_all >= 0.4) & (geno$het_all <= 0.6)),]
#geno <- geno[which((geno$ref_alt <= 0.125) | (geno$ref_alt >= 0.875)),]
# Need hard cutoff to satisfy R/qtl 
geno <- geno[which((geno$ref_alt == 0) | (geno$ref_alt == 1)),]

log("Genotype data filtered...")

#---------------------------------Recoding-------------------------------------#
recoded.genos <- matrix(NA, nrow=dim(geno[,F2_start:F2_end])[1], 
                        ncol=dim(geno[,F2_start:F2_end])[2])
for (i in 1:ncol(recoded.genos)){ # For each F2 column 
  recoded.genos[(geno[,i+F2_start-1] == geno[[A.ref]]), i] <- "AA"
  recoded.genos[(geno[,i+F2_start-1] == geno[[B.ref]]), i] <- "BB"
  recoded.genos[(geno[,i+F2_start-1] == "H"),i] <- "AB"
}
geno[,F2_start:F2_end] <- recoded.genos 

geno[geno=="N"] <- "-" # Replace all "N" with "-"

log("Genotype data recoded...")

#------------------------------Prep for R/qtl----------------------------------#
rqtl <- geno[,-c(F2_end+1:ncol(geno))] # Remove marker stats 
rqtl <- rqtl[,-c(4:(F2_start-1))] # Remove marker metadata (other than name, chr, 
# and pos) and marker data for parent and F1 mice 

rqtl <- t(rqtl) # transpose

# Replace numbered colnames with marker names 
colnames(rqtl) <- rqtl[1,]
rqtl <- rqtl[-1,]

rqtl <- as.data.frame(rqtl)
rqtl <- rownames_to_column(rqtl, var="mouse_ID") %>% as_tibble

#------------------------------------------------------------------------------#
#-------------------------------Phenotype data---------------------------------#
#------------------------------------------------------------------------------#
pheno <- read_xlsx(pheno_file)
pheno$Geno_ID <- toupper(pheno$Geno_ID)
log("Phenotype data downloaded...")

# Before adding phenotype data, filter down to genotyped mice that we have 
# phenotype data for 
rqtl$mouse_ID <- c(rqtl$mouse_ID[1:2], toupper(rqtl$mouse_ID[3:nrow(rqtl)]))
rqtl <- rqtl %>% filter(mouse_ID == "Chromosome" | mouse_ID == "Position" | 
                          mouse_ID %in% pheno$Geno_ID)

# Empty phenotype column (with first two rows empty, as required for R/qtl)
new.col = c(rep("",2), rep(NA, nrow(rqtl)-2))

# Add phenotype data to genotype data
col.pos = 2
for (k in 1:length(pheno.names)){
  pheno.name = pheno.names[k]
  rqtl <- add_pheno(pheno.name, new.col, col.pos)
  col.pos = col.pos+1
}

log("Genotype and phenotype data integrated...")
#------------------------------------------------------------------------------#
#--------------------------------Output data-----------------------------------#
#------------------------------------------------------------------------------#
rqtl[1:2,"mouse_ID"] <- "" # Remove chr/pos row labels

write.csv(rqtl, out_file, row.names=F)
#write.csv(geno.Y.MT, Y_MT_out_file)

log("Output to csv")

