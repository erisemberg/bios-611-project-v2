### Functions for QTL analysis 

#-----------------------------Generic functions--------------------------------#

### Function to ensure directory exists before creating files in it 
ensure_directory <- function(directory){
  if(!dir.exists(directory)){
    dir.create(directory);
  }
}

### Function to create logger function 
make_logger <- function(filename, sep="\n"){
  if(file.exists(filename)){
    file.remove(filename);
  }
  function(...){
    text <- sprintf(...);
    cat(text, file=filename, sep=sep, append=T);
  }
}

#-----------------------------Data pre-processing------------------------------#

### Function to add phenotype data to R/qtl data table 
# Params:
#   pheno.name = name of phenotype column in phenotype spreadsheet 
#   new.col = empty column 
#   col.pos = position of new column in R/qtl table 
add_pheno <- function(pheno.name, new.col, col.pos){
  rqtl <- rqtl %>% add_column(placeholder.name = new.col, .before = col.pos)
  names(rqtl)[names(rqtl) == "placeholder.name"] <- pheno.name
  for (i in 3:nrow(rqtl)){
    rqtl[[pheno.name]][i] <- 
      ifelse(rqtl$mouse_ID[i] %in% pheno$Geno_ID, 
             pheno[[pheno.name]][which(pheno$Geno_ID == rqtl$mouse_ID[i])], NA)
  }
  return(rqtl)
}

#-------------------------Phenotype transformation-----------------------------#

# Truncated Rank Inverse Normal Transformation function 
trint <- function(y, theta=0.01){
  p <- theta + (rank(y)-1)*(1-2*theta)/(length(y)-1)
  qnorm(p)
}





















