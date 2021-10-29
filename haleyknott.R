### This script holds the code for mixed modeling with Haley-Knott regression
### To-do:
###   1. Figure out how to dynamically create random effect terms based on length(randoms).
###   For now, assume there are two (cage and batch) - line 44
###   2. Handle covariates in addition to random effects (also dynamically)
###   3. X-chr approach validated for BC but needs to be validated for F2 (i.e. compare 
###   with results from R/qtl) - line 24

library(lme4) # for modeling random effects 

hk <- function(cross, pheno.col){
  # Calculate genotype probabilities if necessary 
  if (!with(cross$geno[[1]], exists('prob'))){  #calc.genoprob hasn't been run
    cross <- calc.genoprob(cross)
    print("Running calc.genoprob...")
  }
  
  # Calculate dosage (of B allele) 
  if (summary(cross)$type == 'bc'){
    for (c in 1:length(cross$geno)){ # for each chr
      cross$geno[[c]]$dos <- 0*cross$geno[[c]]$prob[,,1] + 1*cross$geno[[c]]$prob[,,2]
    }
  } else if (summary(cross)$type == 'f2'){ ### Needs to be validated  
    for (c in 1:(length(cross$geno)-1)){ # for each autosome
      cross$geno[[c]]$dos <- 0*cross$geno[[c]]$prob[,,1] + 1*cross$geno[[c]]$prob[,,2] + 2*cross$geno[[c]]$prob[,,3]
    }
    cross$geno[[20]]$dos <- 0*cross$geno[[20]]$prob[,,1] + 1*cross$geno[[20]]$prob[,,2] # X-chr; g1 and g2
  }
  
  # Calculate total number of markers 
  total.markers = 0
  for (i in 1:length(cross$geno)){
    total.markers = total.markers + ncol(cross$geno[[i]]$data)
  }
  n = nrow(cross$pheno) # number of mice / observations
  
  # Pre-allocate vectors for results 
  markers <- vector(length = total.markers)
  chrs <- vector(length = total.markers)
  positions <- vector(length = total.markers)
  lods <- vector(length = total.markers)
  
  # Linear model 
  marker.pos = 0
  cage <- cross$pheno[,"Cage"]
  batch <- cross$pheno[,"Batch"]
  for (c in 1:length(cross$geno)){ # for each chr 
    for (k in 1:ncol(cross$geno[[c]]$data)){ # for each marker 
      marker.pos = marker.pos + 1
      #fit <- lmer(cross$pheno[,pheno.col] ~ cross$geno[[c]]$dos[,k] +
      #              (1|cross$pheno[,randoms[1]]) + (1|cross$pheno[,randoms[2]]))
      fit <- lmer(cross$pheno[,pheno.col] ~ cross$geno[[c]]$dos[,k] + (1|batch/cage))
      t <- summary(fit)$coefficients[2,3] # We don't have F, so grab t and square it
      F <- t^2
      df <- 1
      lod <- (n/2)*log10(F*(df/(n-df-1)) + 1) # calculate LOD
      # save data
      marker.name <- colnames(cross$geno[[c]]$data)[k]
      markers[marker.pos] <- marker.name
      chrs[marker.pos] = c
      positions[marker.pos] = cross$geno[[c]]$map[marker.name]
      lods[marker.pos] = lod 
    }
  }
  
  # Create object 
  model.df <- data.frame(chrs, positions, lods) 
  names(model.df) <- c("chr", "pos", "lod") # same column names as scanone object
  rownames(model.df) <- markers 
  class(model.df) <- c("scanone", "data.frame") # so we can use R/qtl functions (see Fig 5.4 in guide)
  
  return(model.df)
}

