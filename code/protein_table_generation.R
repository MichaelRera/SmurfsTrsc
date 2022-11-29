setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #setting current working directory
library(org.Dm.eg.db); library(tidyverse)
library(readODS)
library(readxl)


###open the files

###########################young smurfs and ns
ySvsyNS = read_ods("../data/proteomic/youngSmurfs_youngNonSmurf.ods", sheet = 2, col_names = FALSE, strings_as_factors = FALSE )
colnames(ySvsyNS) = c("Accession", "Peptides", "Score", "Anova_pval", "Fold", "Description", "WT_S", "WT_NS")
ySvsyNS$Accession = as.character(ySvsyNS$Accession)
##add the uniprot
uniprot_ySvsyNS = read_excel("../data/proteomic/uniprot_youngSvsyoungNS.xlsx", col_names = TRUE)
colnames(uniprot_ySvsyNS) = c("Accession", "Uniprot" ,colnames(uniprot_ySvsyNS)[-c(1:2)])
##join
compl_young = left_join(ySvsyNS, uniprot_ySvsyNS, by = 'Accession') %>% mutate(new_FC = log2(WT_S/WT_NS), comparison = 'young') 

###########################young smurfs and ns
oSvsoNS = read_ods("../data/proteomic/oldSmurfs_oldNonSmurfs.ods", sheet = 2, col_names = FALSE, strings_as_factors = FALSE )
colnames(oSvsoNS) = c("Accession", "Peptides", "Score", "Anova_pval", "Fold", "Description", "WT_S", "WT_NS")
oSvsoNS$Accession = as.character(oSvsoNS$Accession)
##add the uniprot
uniprot_oSvsoNS = read_excel("../data/proteomic/uniprot_oldSvsoldNS.xlsx", col_names = TRUE)
colnames(uniprot_oSvsoNS) = c("Accession", "Uniprot" ,colnames(uniprot_oSvsoNS)[-c(1:2)])
##join
compl_old = left_join(oSvsoNS, uniprot_oSvsoNS, by = 'Accession') %>% mutate(new_FC = log2(WT_S/WT_NS), comparison = 'old')

#######
compl_young_sup = compl_young %>% select(c("Accession", "Uniprot", "new_FC", "Anova_pval", "Description"))
colnames(compl_young_sup) = c("GI number", "Uniprot", "log2FoldChange", "Anova p-value", "Description")
compl_old_sup =compl_old %>% select(c("Accession", "Uniprot", "new_FC", "Anova_pval", "Description"))
colnames(compl_old_sup) = c("GI number", "Uniprot", "log2FoldChange", "Anova p-value", "Description")


###generation of tables from which significant proteins where manually extracted for enrichment analaysis with Panther
write.table(compl_young_sup, "../data/proteomic/SvsNS_young_all_protein.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep ="\t")
write.table(compl_old_sup, "../data/proteomic/SvsNS_old_all_protein.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep ="\t")






