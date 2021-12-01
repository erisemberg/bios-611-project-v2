.PHONY: clean 
clean:
	rm -f report.pdf report.tex report.log report-concordance.tex report.synctex.gz
	rm -rf derived_data/*
	rm -rf logs/*
	rm -rf figures/*

derived_data/Rqtl_CC27xC3H_BC.csv derived_data/Geno_CC27xC3H_Y_MT.csv logs/file_processing_notes.md: utils.R \
 source_data/CC27xC3H_backcross_pheno_clean.xlsx \
 source_data/Cr_WB02_miniMUGA-06242021.csv 
	mkdir -p logs
	mkdir -p figures
	mkdir -p derived_data
	Rscript bc_file_proc.R

## Need to figure out a way to generate a large number of figures at once 
## i.e., run make once for genome scans, once for PxG plots, etc. 
figures/genetic_map.png figures/temp_aac_genome_scan.png figures/temp_aac_hist.png figures/temp_aac_transformed_hist.png figures/temp_histograms.png figures/temp_time-series.png figures/temp_time-series_ex.png: \
 utils.R \
 derived_data/Rqtl_CC27xC3H_BC.csv 
	Rscript bc_qtl.R

report.pdf: report.Rnw
	Rscript -e "options(tinytex.tlmgr.path='/opt/tinytex/bin/x86_64-linux/tlmgr'); library(knitr); knit2pdf(\"report.Rnw\")"