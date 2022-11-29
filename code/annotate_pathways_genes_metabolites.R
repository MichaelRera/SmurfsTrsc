###code written by F.Zane (flaminiazane)
##July 2022
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #setting current working directory

########Annotate the KEGG pathways (dme00020 - TCA cycle, dme00061 - fatty acid biosynthesis) with metabolites and gene expression value
library(org.Dm.eg.db);library(tidyverse);library(pathview);library(gage)

#load the pathways
gset = kegg.gsets(species = "dme", id.type = "entrez", check.new = FALSE)
kegg_gset = gset$kg.sets


###upload the metabolites fold change data (from MetaboAnalyst) - independently of the significance to wilcoxon test
metabolite_data = read.csv("../data/metabolites/fold_change_all.csv")
metabolite_data_noNA = metabolite_data %>% na.omit()

cmp_data = metabolite_data_noNA$log2.FC. 
names(cmp_data) = metabolite_data_noNA$KEGG

###upload significant DEGs in Smurfs
sign_genes = read.table("../data/differential_expression/res_DEGS_SvsNS.txt", header = TRUE) %>% filter(padj <= 0.05) %>% mutate(entrez = mapIds(org.Dm.eg.db,
                                                                                                           key = row.names(.),
                                                                                                           column = "ENTREZID", 
                                                                                                           keytype = "ENSEMBL",
                                                                                                           multiVals = "first"))
gene_data = sign_genes$log2FoldChange
names(gene_data) = sign_genes$entrez

##set the folder to save the pathview output (otherwise it goes on current folder)

setwd("../figures/supplementary/")
pathview(gene.data = gene_data, cpd.data = cmp_data, pathway.id = "dme00020",
         species = 'dme', kegg.native = T, limit = list(gene = c(-0.5, 0.5), cpd = 1)) 
#name of modified file: FigS8_tca_cycle_pathview.pdf
pathview(gene.data = gene_data, cpd.data = cmp_data, pathway.id = "dme00061",
         species = 'dme', kegg.native = T, limit = list(gene = c(-0.5, 0.5), cpd = 1)) 
#name of modified file: FigS10_fatty_acid_biosynthesis.pdf


###significant metabolites are then annotated from the diff_metabolites.csv file, in data/metabolites/diff_metabolites.csv



diff_metabolites = read.csv("../data/metabolites/diff_metabolites.csv")
####plot of lactate metabolite 
#merge the two by X column
merged_df = inner_join(diff_metabolites, metabolite_data, by = "X")
metabolite_data %>% filter(X == "HMDB01311")
diff_metabolites %>% filter(X == "HMDB01311")
#it is merged correctly
merged_df[grep("HMDB01311", merged_df$X, ignore.case = TRUE),]


##plot the levels of metabolites of interest
data = read.csv("../data/metabolites/data_without_S8flies.csv", header = T, sep = ",")
data_norm = data  %>% mutate(across(HMDB00606:HMDB00098, ~ .x/Weight))

#check that I am doing the right thing with the function -> yes its working
# tmp = data.frame(a = c(2,4), b = c(6,2), c = c(8,8))
# tmp %>% mutate(across(b:c, ~ .x/a)) 

data_gather = data_norm %>% gather(., value = 'concentration', key = 'metabolite', -Name, -Pheno, -Weight, -Number)
#merge with the names
names = merged_df %>% select(c("X"))
colnames(names) = c("metabolite")
data_gather2 = merge(data_gather, names, by = "metabolite")


##set the age and the Smurf colors 
smurfs = c("azure4", "deepskyblue2")
days = c("darkgoldenrod1", "chocolate1", "brown")
require(ggsignif)
ggplot(data_gather2 %>% filter(metabolite == "HMDB01311"), aes(Pheno, concentration)) + geom_boxplot(aes(fill = Pheno), notch = TRUE, alpha = 0.5) + 
  geom_jitter(size = 1.5) +  scale_fill_manual(name = "condition", values = smurfs) + theme_light() + 
  ylab("Lactate concentration") + xlab("condition") +
  geom_signif(comparisons=list(c("NS", "S")), annotations="**",
              y_position = max(data_gather2 %>% filter(metabolite == "HMDB01311") %>% select(concentration)) + 1, tip_length = 0, vjust=0.4) 
ggsave("../figures/supplementary/Fig_S9_lactate_boxplot.pdf", width = 6, height = 5)


