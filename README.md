# Overview
This repository includes all relevant custom R scripts and associated configuration files for the suspension mass cytometry analysis of the J19113 clinical trial in "A phase 2 trial of CXCR4 antagonist AMD3100 and PD1 inhibitor cemiplimab in metastatic pancreatic adenocarcinoma reveals increased tumor-infiltrating T cells but also immunosuppressive macrophages".

# Input Data
All raw FCS files and fully annotated data frames ("backup_output.rds" for myeloid cells and "backup_output.rds" for T cells) are available at 10.5281/zenodo.11371300. The fully annotated data frames can be loaded onto the R script to generate the published figures in the manuscript. In Config, there are metadata, panel, and merged (annotation) files necessary to generate the heatmaps and other plots for the entire clinical dataset (Supplementary Fig. ).

# R Scripts
Scripts stored in Rscripts need to be run separately for myeloid cells and T cells. The R script J19113_Myeloids.R pertains to the myeloid dataset while the R script J19113_T.R pertains to the T cell data set.
