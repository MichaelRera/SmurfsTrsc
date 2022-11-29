
#code written by flaminiazane
#July 2022

##This code generates Fig 4b,c,d
#It also generates Tabls S4, S5, S6, S7

#############################FIGURE 4b: barplots showing Smurf DEGs and old non-Smurf DEGs overlap###########################

#set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #setting current working directory
##load packages
library(org.Dm.eg.db); library(tidyverse)


###load the required files (DESeq2 results for smurfs and old non-Smurfs)
##Upload all the data needed
res_s = read.table("../data/differential_expression/res_DEGS_SvsNS.txt", header = TRUE)
res_ns = read.table("../data/differential_expression/res_DEGS_NS_40days_20days.txt", header = TRUE)

###################################BP UP WITH AGEING#################
####INFLAMMATION
#load annotation files from flybase
amp = read.table("../data/annotation/FlyBase_IDs_AMP.txt", header = FALSE)# %>% mutate(type = "amp")
toll = read.table("../data/annotation/FlyBase_IDs_Toll.txt", header = FALSE) #%>% mutate(type = "toll")
imd = read.table("../data/annotation/FlyBase_IDs_imd.txt", header = FALSE) #%>% mutate(type = "imd")

infl = rbind(amp, toll, imd) %>% mutate(type = "Toll & Imd pathways")

res_s_infl = res_s %>% filter(rownames(.) %in% infl$V1)
res_ns_infl = res_ns %>% filter(rownames(.) %in% infl$V1)

infl_upset_data = data.frame(rownames(res_s_infl)) %>%  
  cbind(., res_s_infl %>% dplyr::select(c("log2FoldChange", "padj"))) %>% 
  cbind(., res_ns_infl %>% dplyr::select(c("log2FoldChange", "padj"))) %>% 
  .[match(infl %>% distinct %>% filter(V1 %in% rownames(res_ns_infl)) %>% dplyr::select(V1) %>% unlist, rownames(.)),] %>% 
  cbind(., infl %>% distinct %>%  filter(V1 %in% rownames(res_ns_infl)) %>% dplyr::select(type))
colnames(infl_upset_data) = c("gene", "log2FC_s", "padj_s", "log2FC_ns", "padj_ns", "type")
infl_upset_data2 = infl_upset_data %>% replace(is.na(.),1) %>% mutate(all = TRUE,
                                                                      smurf = ifelse(padj_s <= 0.05, TRUE, FALSE),
                                                                      non_smurf = ifelse(padj_ns <= 0.05, TRUE, FALSE))
###STRESS RESPONSE
hsp_all = read.table("../data/annotation/FlyBase_IDs_HSP.txt", header = FALSE) %>% mutate(type = "Chaperones and co-chaperones")
glut_cit = read.table("../data/annotation/FlyBase_IDs_cytGST.txt", header = FALSE) %>% mutate(type = "Cytosolic GSTs")

stress = rbind(hsp_all, glut_cit)


res_s_stress = res_s %>% filter(rownames(.) %in% stress$V1)
res_ns_stress = res_ns %>% filter(rownames(.) %in% stress$V1)

stress_upset_data = data.frame(rownames(res_s_stress)) %>%  
  cbind(., res_s_stress %>% dplyr::select(c("log2FoldChange", "padj"))) %>% 
  cbind(., res_ns_stress %>% dplyr:: select(c("log2FoldChange", "padj"))) %>% 
  .[match(stress %>% filter(V1 %in% rownames(res_ns_stress)) %>% dplyr::select(V1) %>% unlist, rownames(.)),] %>% 
  cbind(., stress %>%  filter(V1 %in% rownames(res_ns_stress)) %>% dplyr::select(type))

colnames(stress_upset_data) = c("gene", "log2FC_s", "padj_s", "log2FC_ns", "padj_ns", "type")
stress_upset_data2 = stress_upset_data %>% replace(is.na(.),1) %>% mutate(all = TRUE,
                                                                          smurf = ifelse(padj_s <= 0.05, TRUE, FALSE),
                                                                          non_smurf = ifelse(padj_ns <= 0.05, TRUE, FALSE))


stress_infl_barplot = rbind(infl_upset_data2, stress_upset_data2) %>% 
  mutate(annotation = ifelse(padj_s > 0.05 & padj_ns > 0.05 , "not DEG", 
                             ifelse(padj_s <= 0.05 & padj_ns <= 0.05, "DEG in Smurf & old non-Smurf",
                                    ifelse(padj_s <= 0.05 & padj_ns > 0.05, "Smurf DEG", "old non-Smurf DEG")))) %>% 
  dplyr::select(c("type", "annotation")) %>% 
  group_by(type, annotation) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count)) %>% as.data.frame() %>% mutate(perc2 = round(.$perc, 2))
stress_infl_barplot$annotation = factor(stress_infl_barplot$annotation, levels = c("old non-Smurf DEG", "DEG in Smurf & old non-Smurf", "Smurf DEG", "not DEG" ))


ggplot(stress_infl_barplot, aes(x = factor(type),  y = perc, fill = annotation, label = scales::percent(perc2))) +
  geom_bar(stat="identity", width = 0.5)+ #, position = position_dodge()) +
  labs(x = "Categories (stress response) ", y = "% detected genes for category", fill = "Analysis result") +
  scale_fill_manual(values = c( "grey73", "lightblue1",  "deepskyblue2", "khaki")) +
  theme_minimal(base_size = 14) +
  #ylim(c(0,1)) +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_text(position = position_dodge(width = 0),    # move to center of bars
            vjust = -0.5,    # nudge above top of bar
            size = 5) + #adjust annotation manually
  theme(axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 14, color = "black"), 
        axis.title.y = element_text(size = 14, color = "black"), legend.position = "top")

ggsave("../figures/Fig4b_ii.pdf", width = 10, height = 6)



#################CATEGORIES DOWNREGULATED WITH AGEING##################################

####ETC GENES

##complex I
complexI = read.table("../data/annotation/FlyBase_IDs_I.txt") %>% mutate(complex = 'Complex I')
##complex II
complexII = read.table("../data/annotation/FlyBase_IDs_II.txt") %>% mutate(complex = 'Complex II')
##complex III
complexIII = read.table("../data/annotation/FlyBase_IDs_III.txt") %>% mutate(complex = 'Complex III')
##complex IV
complexIV =  read.table("../data/annotation/FlyBase_IDs_IV.txt") %>% mutate(complex = 'Complex IV')
##complex V 
complexV = read.table("../data/annotation/FlyBase_IDs_V.txt") %>% mutate(complex = 'Complex V')
#all
mit = rbind(complexI, complexII, complexIII, complexIV, complexV)

res_s_mit = res_s %>% filter(rownames(.) %in% mit$V1)
res_ns_mit = res_ns %>% filter(rownames(.) %in% mit$V1)

mit_upset_data = data.frame(rownames(res_s_mit)) %>%  
  cbind(., res_s_mit %>% dplyr::select(c("log2FoldChange", "padj"))) %>% 
  cbind(., res_ns_mit %>% dplyr::select(c("log2FoldChange", "padj"))) %>% 
  .[match(mit %>% filter(V1 %in% rownames(res_ns_mit)) %>% dplyr::select(V1) %>% unlist, rownames(.)),] %>% 
  cbind(., mit %>% filter(V1 %in% rownames(res_ns_mit)) %>% dplyr::select(complex)) 

colnames(mit_upset_data) = c("gene", "log2FC_s", "padj_s", "log2FC_ns", "padj_ns", "complex")
mit_upset_data2 = mit_upset_data %>%  replace(is.na(.),1) %>% mutate(all = TRUE,
                                                                     smurf = ifelse(padj_s <= 0.05, TRUE, FALSE),
                                                                     non_smurf = ifelse(padj_ns <= 0.05, TRUE, FALSE))

mit_barplot_data = mit_upset_data2 %>% mutate(annotation = ifelse(padj_s > 0.05 & padj_ns > 0.05 , "not DEG", 
                                                                  ifelse(padj_s <= 0.05 & padj_ns <= 0.05, "DEG in Smurf & old non-Smurf",
                                                                         ifelse(padj_s <= 0.05 & padj_ns > 0.05, "Smurf DEG", "old non-Smurf DEG")))) %>% 
  dplyr::select(c("complex", "annotation")) %>% 
  group_by(complex, annotation) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

mit_barplot_data$annotation = factor(mit_barplot_data$annotation, levels = c("not DEG","Smurf DEG", "DEG in Smurf & old non-Smurf", "old non-Smurf DEG"))

####RIBOSOMES
cyt_rib = read.table("../data/annotation/FlyBase_IDs_cytRIB.txt", header = FALSE) %>% mutate(type = "Cytoplasmatic ribosomal proteins")
#get genes for the ribosome biogenesis
require(gage)
gset_go = go.gsets(species = "dme", pkg.name = "org.Dm.eg.db",id.type = "eg" )
# #select the ontology you are interest in. Here we select BP, biological process
go_sets = gset_go$go.sets[gset_go$go.subs$BP]

rib_bio = data.frame(entrez = go_sets$`GO:0042254`, flybase =  mapIds(org.Dm.eg.db,
                                                                      key = go_sets$`GO:0042254`,
                                                                      column = "ENSEMBL", 
                                                                      keytype = "ENTREZID",
                                                                      multiVals = "first")) %>% na.omit() %>% mutate(type = "GO:0042254 ribosome biogenesis")

colnames(rib_bio) = c("entrez", "V1", "type")

rib = rbind(cyt_rib, rib_bio %>% dplyr::select(c("V1", "type")))

res_s_rib = res_s %>% filter(rownames(.) %in% rib$V1)
res_ns_rib = res_ns %>% filter(rownames(.) %in% rib$V1)

rib_upset_data = data.frame(rownames(res_s_rib)) %>%  
  cbind(., res_s_rib %>% dplyr::select(c("log2FoldChange", "padj"))) %>% 
  cbind(., res_ns_rib %>% dplyr::select(c("log2FoldChange", "padj"))) %>% 
  .[match(rib %>% filter(V1 %in% rownames(res_ns_rib)) %>% dplyr::select(V1) %>% unlist, rownames(.)),] %>% 
  cbind(., rib %>% filter(V1 %in% rownames(res_ns_rib)) %>% dplyr::select(type)) 

colnames(rib_upset_data) = c("gene", "log2FC_s", "padj_s", "log2FC_ns", "padj_ns", "type")
rib_upset_data2 = rib_upset_data %>%  replace(is.na(.),1) %>% mutate(all = TRUE,
                                                                     smurf = ifelse(padj_s <= 0.05, TRUE, FALSE),
                                                                     non_smurf = ifelse(padj_ns <= 0.05, TRUE, FALSE))

rib_barplot_data = rib_upset_data2 %>% 
  mutate(annotation = ifelse(padj_s > 0.05 & padj_ns > 0.05 , "not DEG", 
                             ifelse(padj_s <= 0.05 & padj_ns <= 0.05, "DEG in Smurf & old non-Smurf",
                                    ifelse(padj_s <= 0.05 & padj_ns > 0.05, "Smurf DEG", "old non-Smurf DEG")))) %>% 
  dplyr::select(c("type", "annotation")) %>% 
  group_by(type, annotation) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))
rib_barplot_data$annotation = factor(rib_barplot_data$annotation, levels = c("not DEG","Smurf DEG", "DEG in Smurf & old non-Smurf", "old non-Smurf DEG"))


####PROTEOLYSIS GO:0006508
prot = data.frame(entrez = go_sets$`GO:0006508`, flybase =  mapIds(org.Dm.eg.db,
                                                                   key = go_sets$`GO:0006508`,
                                                                   column = "ENSEMBL", 
                                                                   keytype = "ENTREZID",
                                                                   multiVals = "first")) %>%
  na.omit() %>% mutate(type = "GO:0006508 proteolysis") %>% dplyr::select(-entrez)


colnames(prot) = c("V1", "type")


res_s_prot = res_s %>% filter(rownames(.) %in% prot$V1)
res_ns_prot = res_ns %>% filter(rownames(.) %in% prot$V1)

prot_upset_data = data.frame(rownames(res_s_prot)) %>%  
  cbind(., res_s_prot %>% dplyr::select(c("log2FoldChange", "padj"))) %>% 
  cbind(., res_ns_prot %>% dplyr::select(c("log2FoldChange", "padj"))) %>% 
  .[match(prot %>% filter(V1 %in% rownames(res_ns_prot)) %>% dplyr::select(V1) %>% unlist, rownames(.)),] %>% 
  cbind(., prot %>% filter(V1 %in% rownames(res_ns_prot)) %>% dplyr::select(type)) 

colnames(prot_upset_data) = c("gene", "log2FC_s", "padj_s", "log2FC_ns", "padj_ns", "type")
prot_upset_data2 = prot_upset_data %>%  replace(is.na(.),1) %>% mutate(all = TRUE,
                                                                       smurf = ifelse(padj_s <= 0.05, TRUE, FALSE),
                                                                       non_smurf = ifelse(padj_ns <= 0.05, TRUE, FALSE))

prot_barplot_data = prot_upset_data2 %>% 
  mutate(annotation = ifelse(padj_s > 0.05 & padj_ns > 0.05 , "not DEG", 
                             ifelse(padj_s <= 0.05 & padj_ns <= 0.05, "DEG in Smurf & old non-Smurf",
                                    ifelse(padj_s <= 0.05 & padj_ns > 0.05, "Smurf DEG", "old non-Smurf DEG")))) %>% 
  dplyr::select(c("type", "annotation")) %>% 
  group_by(type, annotation) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))
prot_barplot_data$annotation = factor(prot_barplot_data$annotation, levels = c("not DEG","Smurf DEG", "DEG in Smurf & old non-Smurf", "old non-Smurf DEG"))


###OOGENESIS GO BP category

egg = data.frame(entrez = go_sets$`GO:0048477`, flybase =  mapIds(org.Dm.eg.db,
                                                                  key = go_sets$`GO:0048477`,
                                                                  column = "ENSEMBL", 
                                                                  keytype = "ENTREZID",
                                                                  multiVals = "first")) %>%
  na.omit() %>% mutate(type = "GO:0048477 oogenesis") %>% dplyr::select(-entrez)


colnames(egg) = c("V1", "type")


res_s_egg = res_s %>% filter(rownames(.) %in% egg$V1)
res_ns_egg = res_ns %>% filter(rownames(.) %in% egg$V1)

egg_upset_data = data.frame(rownames(res_s_egg)) %>%  
  cbind(., res_s_egg %>% dplyr::select(c("log2FoldChange", "padj"))) %>% 
  cbind(., res_ns_egg %>% dplyr::select(c("log2FoldChange", "padj"))) %>% 
  .[match(egg %>% filter(V1 %in% rownames(res_ns_egg)) %>% dplyr::select(V1) %>% unlist, rownames(.)),] %>% 
  cbind(., egg %>% filter(V1 %in% rownames(res_ns_egg)) %>% dplyr::select(type)) 

colnames(egg_upset_data) = c("gene", "log2FC_s", "padj_s", "log2FC_ns", "padj_ns", "type")
egg_upset_data2 = egg_upset_data %>%  replace(is.na(.),1) %>% mutate(all = TRUE,
                                                                     smurf = ifelse(padj_s <= 0.05, TRUE, FALSE),
                                                                     non_smurf = ifelse(padj_ns <= 0.05, TRUE, FALSE))



egg_barplot_data = egg_upset_data2 %>% 
  mutate(annotation = ifelse(padj_s > 0.05 & padj_ns > 0.05 , "not DEG", 
                             ifelse(padj_s <= 0.05 & padj_ns <= 0.05, "DEG in Smurf & old non-Smurf",
                                    ifelse(padj_s <= 0.05 & padj_ns > 0.05, "Smurf DEG", "old non-Smurf DEG")))) %>% 
  dplyr::select(c("type", "annotation")) %>% 
  group_by(type, annotation) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))
egg_barplot_data$annotation = factor(egg_barplot_data$annotation, levels = c("not DEG","Smurf DEG", "DEG in Smurf & old non-Smurf", "old non-Smurf DEG"))


##put all together to generate the barplot
##rename the columns of mit_barplot to be able to rbind 

colnames(mit_barplot_data) = c("type", "annotation", "count", "perc") 

down_genes = rbind(mit_barplot_data, rib_barplot_data, prot_barplot_data, egg_barplot_data) %>%  as.data.frame() %>% mutate(perc2 = .$perc %>% round(.,2))
down_genes$annotation = factor(down_genes$annotation, levels = c("old non-Smurf DEG", "DEG in Smurf & old non-Smurf", "Smurf DEG", "not DEG" ))

ggplot(down_genes, aes(x = factor(type), y = perc, fill = factor(annotation), label = scales::percent(perc2))) +
  geom_bar(stat="identity", width = 0.7)+ #, position = position_dodge()) +
  labs(x = "Categories", y = "% detected genes for category", fill = "Analysis result") +
  scale_fill_manual(values = c( "grey73", "lightblue1",  "deepskyblue2", "khaki")) +
  theme_minimal(base_size = 14) +
  #ylim(c(0,1)) +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_text(position = position_dodge(width = 0),    # move to center of bars
            vjust = -0.5,    # nudge above top of bar
            size = 5) + #adjust the 
  theme(axis.text.x = element_text(size = 14, color = "black", angle = 70, vjust =  0.6, hjust = 0.5),
        axis.text.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 14, color = "black"), 
        axis.title.y = element_text(size = 14, color = "black"), legend.position = "top")

ggsave("../figures/Fig4b_i.pdf", width = 11, height = 9)



##########INSULINE SIGNALING
###INSULIN SIGNALING
ins_core = read.table("../data/annotation/FlyBase_IDs_insulinCORE.txt", header = FALSE) %>% mutate(type = "Core components")
ins_neg =  read.table("../data/annotation/FlyBase_IDs_insulinNEG.txt", header = FALSE) %>% mutate(type = "Negative regulators")
ins_pos =  read.table("../data/annotation/FlyBase_IDs_insulinPOS.txt", header = FALSE) %>% mutate(type = "Positive regulators")

ins = rbind(ins_core, ins_neg, ins_pos)

res_s_ins = res_s %>% filter(rownames(.) %in% ins$V1)
res_ns_ins = res_ns %>% filter(rownames(.) %in% ins$V1)

ins_upset_data = data.frame(rownames(res_s_ins)) %>%  
  cbind(., res_s_ins %>% dplyr::select(c("log2FoldChange", "padj"))) %>% 
  cbind(., res_ns_ins %>% dplyr::select(c("log2FoldChange", "padj"))) %>% 
  .[match(ins %>% filter(V1 %in% rownames(res_ns_ins)) %>% dplyr::select(V1) %>% unlist, rownames(.)),] %>% 
  cbind(., ins %>% filter(V1 %in% rownames(res_ns_ins)) %>% dplyr::select(type)) 

colnames(ins_upset_data) = c("gene", "log2FC_s", "padj_s", "log2FC_ns", "padj_ns", "type")
ins_upset_data2 = ins_upset_data %>%  replace(is.na(.),1) %>% mutate(all = TRUE,
                                                                     smurf = ifelse(padj_s <= 0.05, TRUE, FALSE),
                                                                     non_smurf = ifelse(padj_ns <= 0.05, TRUE, FALSE))

ins_barplot_data = ins_upset_data2 %>% 
  mutate(annotation = ifelse(padj_s > 0.05 & padj_ns > 0.05 , "not DEG", 
                             ifelse(padj_s <= 0.05 & padj_ns <= 0.05, "DEG in Smurf & old non-Smurf",
                                    ifelse(padj_s <= 0.05 & padj_ns > 0.05, "Smurf DEG", "old non-Smurf DEG")))) %>% 
  dplyr::select(c("type", "annotation")) %>% 
  group_by(type, annotation) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))
ins_barplot_data$annotation = factor(ins_barplot_data$annotation,  levels = c("old non-Smurf DEG", "DEG in Smurf & old non-Smurf", "Smurf DEG", "not DEG" ))


ggplot(ins_barplot_data, aes(x = factor(type), y = perc, fill = factor(annotation), label = scales::percent(perc))) +
  geom_bar(stat="identity", width = 0.5)+ #, position = position_dodge()) +
  labs(x = "Insulin signaling sub-categories ", y = "% detected genes for category", fill = "Analysis result") +
  scale_fill_manual(values = c("deepskyblue2", "khaki")) +
  theme_minimal(base_size = 14) +
  #ylim(c(0,1)) +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_text(position = position_dodge(width = 0),    # move to center of bars
            #vjust = -0.5,    # nudge above top of bar
            size = 5) +
  theme(axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"), 
        axis.title.x = element_text(size = 14, color = "black"), 
        axis.title.y = element_text(size = 14, color = "black"), legend.position = "top")
ggsave("../figures/Fig4b_i.pdf", width = 10, height = 6)



############################################FIGURE 4c: Human Atlas ageing genes################################################## 

#get the ortholog human genes
#file from: https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#Human_Orthologs_.28dmel_human_orthologs_disease_fb_.2A.tsv.gz.29
annotation_file = read_tsv("../data/annotation/dmel_human_orthologs_disease_fb_2022_05.tsv")
colnames(annotation_file) = c("flybase", colnames(annotation_file[-1]))

####first check for the SMURF DEGS
##degs smurfs
smurf_degs = res_s %>% rownames_to_column(var = "flybase") %>% filter(padj <= 0.05)

degs_annotated = left_join(smurf_degs, annotation_file, by = "flybase")

###tables from https://ngdc.cncb.ac.cn/aging/age_related_genes
degs_to_intersect = degs_annotated %>% select(c("flybase", "log2FoldChange", "padj", "Human_gene_symbol"))

human_ageing_table = rbind(read_csv("../data/annotation/HAA_tableExport1.txt", col_names = TRUE),
                           read_csv("../data/annotation/HAA_tableExport2.txt", col_names = TRUE),
                           read_csv("../data/annotation/HAA_tableExport3.txt", col_names = TRUE),
                           read_csv("../data/annotation/HAA_tableExport4.txt", col_names = TRUE),
                           read_csv("../data/annotation/HAA_tableExport5.txt", col_names = TRUE))
colnames(human_ageing_table) = c("Human_gene_symbol", colnames(human_ageing_table)[-1])

intersected = left_join(human_ageing_table, degs_to_intersect, by = 'Human_gene_symbol')
intersected2 = intersected %>% filter(!is.na(flybase))
intersected2$Human_gene_symbol %>% unique %>% length
#there are 134 human genes intersected, corresponding to 121 flybase names 

###Make a table for the Genes overlapping
table_sel = intersected2 %>%  mutate(padj2 = ifelse(padj < 0.001, "***", ifelse(padj < 0.01 & padj > 0.001, "**", "*"))) %>% 
  select(c("Human_gene_symbol",  "flybase", "Gene_Set", "log2FoldChange", "padj2"))

write.table(table_sel, "../data/human_ageing_atlas_genage_results/TABLE_S4_SMURFDEGs_human_ageing_atlas.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

#code to generate the pdf version of the table
require(kableExtra)
table_colnames = c("Human symbol", "Flybase", "Ageing marker", "log2FC (DESeq2)", "FDR (DESeq2)")
kbl_s_human = kbl(table_sel, col.names = table_colnames, row.names = F) %>% kable_classic_2(full_width = F, html_font = "Helvetica")%>% 
  row_spec(row = 0, bold = T, font_size = 14) 
# save_kable(kbl_s_human, "../../tables_and_files/new_tables_atlas/kbl_s_human.html")
# save_kable(kbl_s_human, "../../tables_and_files/new_tables_atlas/kbl_s_human.pdf")


###degs old non-Smurfs
###check the same for the genes that are non-smurfs
ns_degs = res_ns %>% filter(padj <= 0.05) %>% rownames_to_column(var = "flybase")
ns_degs_annotated = left_join(ns_degs, annotation_file, by = "flybase")
ns_degs_to_intersect = ns_degs_annotated %>% select(c("flybase", "log2FoldChange", "padj", "Human_gene_symbol"))
intersected_ns = left_join(human_ageing_table, ns_degs_to_intersect, by = 'Human_gene_symbol')
intersected2_ns = intersected_ns %>% filter(!is.na(flybase))

intersected2_ns$Human_gene_symbol %>% unique %>% length
#25 human genes

###Make a table for the Genes overlapping
table_sel_ns = intersected2_ns %>% mutate(padj2 = ifelse(padj < 0.001, "***", ifelse(padj < 0.01 & padj > 0.001, "**", "*"))) %>% 
  select(c("Human_gene_symbol",  "flybase", "Gene_Set", "log2FoldChange", "padj2"))
write.table(table_sel_ns, "../data/human_ageing_atlas_genage_results/TABLE_S5_SMURFDEGsOldNonSmurfs_human_ageing_atlas.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

#require(kableExtra)
table_colnames = c("Human symbol", "Flybase", "Ageing marker", "log2FC (DESeq2)", "FDR (DESeq2)")
#generate the pdf
kbl_ns_human = kbl(table_sel_ns, col.names = table_colnames, row.names = F) %>% kable_classic_2(full_width = F, html_font = "Helvetica") %>% row_spec(row = 0, bold = T, font_size = 14) 
# save_kable(kbl_ns_human, "../../tables_and_files/new_tables_atlas/kbl_ns_human.html")
# save_kable(kbl_ns_human, "../../tables_and_files/new_tables_atlas/kbl_ns_human.pdf")

####make the general barplots
#for smurfs --> 134 genes out of 500
smurf_pie_chart = data.frame(type = c("Detected as DEGs", "Not detected as DEGs"), values = c(134, 366)) %>% mutate(analysis = "Smurf")
##for non-smurfs -> 25 out of 500
ns_smurf_pie_chart = data.frame(type = c("Detected as DEGs", "Not detected as DEGs"), values = c(25, 475)) %>% mutate(analysis = "non-Smurf")
##plot together
all_overlap = rbind(ns_smurf_pie_chart, smurf_pie_chart)
ggplot(all_overlap, aes(x= analysis, y=values, fill=type)) +
  geom_bar(stat="identity", width=0.5) +
  scale_fill_manual(values = c("deepskyblue2", "khaki")) +
  #coord_polar("y", start=0) +
  coord_flip() +
  geom_text(aes(label=values), position=position_stack(0.5)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"), legend.position = "top") +
  ylab("Genes Annotated in Human Ageing Atlas") + xlab("")
ggsave("../figures/Fig4c_human_ageing_atlas_comparison.pdf", width = 6, height = 5)


############################################FIGURE 4d: GenAge database################################################## 

##open the file, data from GenAge (longevity genes in Drosophila melanogaster)
require(readODS)
gene_longevity = read_ods("../data/annotation/genage_models_export.ods")

#open the results
res_s = res_s %>% rownames_to_column(var = 'flybase') %>% 
  mutate(entrez = mapIds(org.Dm.eg.db, key = .$flybase, keytype = "ENSEMBL", column = 'ENTREZID' )) %>% 
  mutate(symbol = mapIds(org.Dm.eg.db, key = .$flybase, keytype = "ENSEMBL", column = 'SYMBOL' )) 
res_ns = res_ns  %>% rownames_to_column(var = 'flybase') %>% 
  mutate(entrez = mapIds(org.Dm.eg.db, key = .$flybase, keytype = "ENSEMBL", column = 'ENTREZID' )) 


unique_long = gene_longevity$`Entrez Gene ID` %>% unique %>% as.character()
##check the ones that are not detected
unique_long_detected = unique_long[unique_long %in% res_s$entrez]
unique_long[!(unique_long %in% unique_long_detected)]
#16 genes, 8 have a flybase annotation, 6 are miRNA and 2 are human genes artificially expressed in drosophila

res_s %>% filter(flybase %in% c("FBgn0033238", "FBgn0004407", "FBgn0000482", "FBgn0035468", "FBgn0001223", "FBgn0263490", "FBgn0038570", "FBgn0037802"))

res_s = res_s %>% mutate(entrez = ifelse(flybase == "FBgn0000482", "31118", entrez)) %>% #manual annotation
  mutate(entrez = ifelse(flybase == "FBgn0037802", "41254", entrez)) %>% #manual annotation
  mutate(entrez = ifelse(flybase == "FBgn0035468", "38453", entrez)) %>% #manual annotation 
  mutate(entrez = ifelse(flybase == "FBgn0263490", "249462", entrez)) #replacement as the entrez in the GenAge refers to a specific allele of this gene
##I can retrieve 4 our of the 8 genes with flybase annotation, then I loose the 6 miRNA and the other 2 because are expression of human variants in flies (OTC and UCP2)
#check the ones that are not retrived 
unique_long_detected2 = unique_long[unique_long %in% res_s$entrez]

sign_s = res_s %>% filter(padj <= 0.05)
sign_ns = res_ns %>% filter(padj <= 0.05)

###
(unique_long %>% na.omit())  %in% sign_s$entrez %>% sum
(unique_long %>% na.omit())  %in% sign_ns$entrez  %>% sum 

##figure is genereted by using https://eulerr.co/

##supplementary table
################CREATE TABLE FOR PAPER
extract = function(x) {stringr::str_extract_all( x, "(?<=: ).+(?= \")" ) %>% unlist()}
gene_longevity_table = gene_longevity %>% filter(`Entrez Gene ID` %in% sign_s$entrez) %>% mutate(biblio = extract(.$`Bibliographic reference`)) %>% 
  dplyr::select(c("Gene Symbol", "Entrez Gene ID", "Lifespan Effect", "Avg Lifespan Change", "Method", "biblio")) %>% filter(`Gene Symbol` != "hsp70") #to remove the NA in the entrez
colnames(gene_longevity_table) = c("Gene Symbol", "entrez", "Lifespan Effect", "Avg Lifespan Change", "Method", "biblio")
gene_longevity_table$entrez = as.character(gene_longevity_table$entrez)
gene_longevity_table_joined = left_join(gene_longevity_table, sign_s, by = "entrez") %>% 
  dplyr::select(c("Gene Symbol", "log2FoldChange", "Lifespan Effect", "Avg Lifespan Change", "Method", "biblio")) %>% 
  mutate(`Avg Lifespan Change` = ifelse(is.na(`Avg Lifespan Change`), "not reported", round(`Avg Lifespan Change`, 0))) %>% filter(`Lifespan Effect` != "noeffect")

write.table(gene_longevity_table_joined, "../data/human_ageing_atlas_genage_results/TABLE_S6_SMURFDEGs_genAge.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)


table_colnames = c("Symbol", "log2FC", "Effect", "% effect", "Method", "Reference")
kbl_en = kbl(gene_longevity_table_joined, col.names = table_colnames, row.names = F) %>% kable_classic_2(full_width = F, html_font = "Helvetica")%>% 
  row_spec(row = 0, bold = T, font_size = 14) 
#save_kable(kbl_en, "../tables_and_files/geneAGE_Smurf.pdf")
#save_kable(kbl_en, "../tables_and_files/geneAGE")        

###for the non-Smurf
gene_longevity_table_joined_ns= left_join(gene_longevity_table, sign_ns, by = "entrez") %>% 
  dplyr::select(c("Gene Symbol", "log2FoldChange", "Lifespan Effect", "Avg Lifespan Change", "Method", "biblio")) %>% filter(!is.na(log2FoldChange))  %>% 
  mutate(`Avg Lifespan Change` = ifelse(is.na(`Avg Lifespan Change`), "not reported", round(`Avg Lifespan Change`, 0)))

write.table(gene_longevity_table_joined_ns, "../data/human_ageing_atlas_genage_results/TABLE_S6_OldNonSmurfDEGs_genAge.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)


table_colnames = c("Symbol", "log2FC", "Effect", "% effect", "Method", "Reference")
kbl_en_ns = kbl(gene_longevity_table_joined_ns, col.names = table_colnames, row.names = F) %>% kable_classic_2(full_width = F, html_font = "Helvetica")%>% 
  row_spec(row = 0, bold = T, font_size = 14) 
#save_kable(kbl_en_ns, "../tables_and_files/geneAGE_NONSmurf.pdf")