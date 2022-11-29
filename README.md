# SmurfsTrsc

The repository refers to the "Smurfness-based two-phase model of ageing helps deconvolve the ageing transcriptional signature" study, currently avaiable on bioRxiv (doi: https://doi.org/10.1101/2022.11.22.517330).
Author : Flaminia Zane, Michael Rera

## Repository organization

The repository is organized in three main folders. "Code" contains the .R or .Rmd files where the code of the analysis presented in the study is stored; "data" contains the input data used for the analysis, as well as the output results (tables) of the study; 'figures' contains the figures generated by the code and used as main and supplementary figure in the paper. Note the the folder "code" contains an "R_env" subfolder, where all the information about the R version and the packages used is stored. 

## Code 

### RNA-seq analysis

All code used for the RNAseq analysis illustrated in our study is stored here. More specifically:

- **analysis_rnaseq_smurf_non_smurf.Rmd** details the analysis performed on the RNA-seq data to investigate the differences between Smurfs and non-Smurfs. This code generates Fig. 2 and 3, as well as figure S3, S4, S5 and the associated data tables.

- **analysis_rnaseq_nonSmurf.Rmd** details the analysis carried on the non-Smurfs along chronological age. The same framework of the Smurf/non-Smurf analysis is followed. This code generates Fig 4a and 6a, and the associated data tables.

- **analysis_rnaseq_smurf.Rmd** details the analysis carried on the Smurfs along chronological age. The same framework as the previous analysis is followed This code generates Fig 5, and the associated data tables.

- **Figure4_b_c_d.R** as specified by the name, this code details the analysis performed to generate figure 4b, 4c, 4d and the associated data tables.

- **analysis_hetereogeneity_gene_expression.Rmd** details the analysis carried to investigate the effect of chronological age and Smurfness on transcriptional noise (Fig 6b and S12).



- **analysis_age_and_smurfness.Rmd** details the analysis carried to investigate the effect of Smurfness and chronological age on biological pathways (Fig 6c and associated data tables).

### Proteomic and metabolomic analysis

Metabolomic analysis was performed on MetaboAnalyst (https://www.metaboanalyst.ca/, which generated Fig S7). The code used to illustrate the results (Fig. S8, S9 and S10) is stored in **annotate_pathways_genes_metabolites.R**.

Proteomic data were first analyzed using Excel. Results were then manipulated with R (**protein_table_generation.R**) and enrichment analysis was performed on Pantherdb (http://www.pantherdb.org/).


### Longevity experiments analysis
The analysis of the longevity experiments data presented in the study is divided into three files:

- **analysis_longevity_first_screening.Rmd** illustrates the results and analysis of the first longevity screening (Fig S13 and associated data table). It sources the function_get_deaths_effect_pval.R file.

- **analysis_longevity_data_figures.R**  details all the longevity analysis linked to the three longevity genes identified (Fig 7b, S14, S17, S18, S19).

- **analysis_longevity_smurf_counts** illustrates the analysis of the smurf data associated to the presented longevity experiments (Fig 7b, S17).


## Data

The folder contains all the data used as an input for the code, as well as the output generated by the analysis (that can, in some cases, be used as input for another analysis). The folder contains subfolders to facilitate the data retrieving.

- **RNAseq_raw_data_matrix.tabular** is the output of featureCounts, i.e. the non-normalized counts matrix.

- The **annotation** folder contains the annotation data (pathway, longevity genes from GenAge, etc..) that are used in the analysis. Results of the GenAge and Human Ageing Atlas comparison are stored in **human_ageing_atlas_genage_results** folder.

- The **correlation_gene_age_smurfness** folder contains the output (tables) of the **analysis_age_and_smurfness.Rmd** file.

-The **differential_expression** folder contains all the output of the differential gene expression analysis, included GSEA.

- The **longevity** folder contains all the data of the longevity and smurf recording experiments.

- **Metabolites** and **proteomic** folder contains the input files for the metabolomic and proteomic R part of the analysis.

- The **nonSmurf_lm** folder contains the results of the linear regression of gene expression over time for the non-Smurfs samples.



## Figures

The figures folder contains the figures generated by the code, as well as the final figured as found in the main text of the paper (**Final_pictures_main_text**). Individual figures are reported here with the figure number and a short description. The **longevity_plots** folder stores the longevity plots and smurf analysis of the main text, while the supplementary figures can be found in the folder **supplementary**.


<br/>

Please note that all the code is organized so to work according to the organization of the repository. If you wish to run the code without manually modify the file calling (and saving), download the entire repository without changing its organization.


