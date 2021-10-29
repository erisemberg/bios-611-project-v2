PHONY: purge
PHONY: clean 

purge:
	rm source_data/*
	make clean

clean:
	rm derived_data/*
	rm logs/*
	rm figures/*

derived_data/Rqtl_CC27xC3H_BC.csv derived_data/Geno_CC27xC3H_Y_MT.csv logs/file_processing_notes.md: \
 utils.R \
 source_data/CC27xC3H_backcross_pheno_clean.xlsx \
 source_data/Cr_WB02_miniMUGA-06242021.csv 
	Rscript bc_file_proc.R

# Need to figure out a way to generate a large number of figures at once 
# i.e., run make once for genome scans, once for PxG plots, etc. 
figures/genetic_map.png figures/temp_aac_genome_scan.png figures/temp_aac_hist.png figures/temp_aac_transformed_hist.png figures/temp_histograms.png figures/temp_time-series.png: \
 utils.R \
 derived_data/Rqtl_CC27xC3H_BC.csv 
	Rscript bc_qtl.R