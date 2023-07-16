####Infer a linear model on gene expression using smurfness, age and the interaction as variable
####this analysis is mentioned in the main text at the end of the section: "Old Smurfs carry additional related changes"
#set the working directory if needed
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #setting current working directory
#load here libraries that are needed all along
library(org.Dm.eg.db); library(tidyverse)

##the data
##use vst data 
vst_df = read.table("../data/vst_df.txt", header = TRUE)
vst_for_lm = vst_df %>% rownames_to_column("flybase") %>%  gather(., value = 'expression', key = 'sample', -flybase) %>% 
  mutate(time = ifelse(grepl("20",.$sample), "20", 
                       ifelse(grepl("40",.$sample), "40", "30"))) %>% 
  mutate(smurfness = ifelse(grepl("n", .$sample), "NS", "S"))
vst_for_lm$time = as.numeric(vst_for_lm$time)  
vst_for_lm$smurfness = as.factor(vst_for_lm$smurfness)


#the model to infer
#g = k1S + k2A + k3SxA (k1, k2, k3 are the parameters and S is Smurfness and A is age)

###now running on the single genes
# Function to fit the linear model and extract coefficients
fit_model <- function(gene_data) {
  model <- lm(expression ~ smurfness + time + smurfness:time, data = gene_data)
  coefficients <- summary(model)$coefficients
  return(data.frame(gene_name = gene_data$flybase[1], coefficients))
}
# Apply the function to each gene in the data frame
results <- lapply(split(vst_for_lm, vst_for_lm$flybase), fit_model)

# Combine the results into a single data frame
results_df <- do.call(rbind, results)

# Print the results data frame
print(results_df)

###reshape the file
results_df$gene_name = as.character(results_df$gene_name)
results_df = results_df %>% rownames_to_column("info") %>% 
                   mutate(result_type = ifelse(grepl("Intercept", .$info), "intercept", 
                                               ifelse(grepl("smurfnessS:time", .$info), "interaction", 
                                               ifelse(grepl(".time", .$info), "time", "smurfness")))) %>% mutate(symbol = mapIds(org.Dm.eg.db,
                                                                                                                                          keys = .$gene_name,
                                                                                                                                                       keytype="ENSEMBL",
                                                                                                                                                       column="SYMBOL", multiVals = "first")) 
                              
  
###select the smurfness (k1)
smurf_sign = results_df %>% filter(result_type == "smurfness") %>% filter(Pr...t.. < 0.05) #2549
write.table(smurf_sign, 'smurf_sign_lm.txt', sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
###select the time (k2)
time_sign = results_df %>% filter(result_type == "time") %>% filter(Pr...t.. < 0.05) #2027
write.table(time_sign, 'time_sign_lm.txt', sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)

##select the interaction (k3)
interaction_sign = results_df %>% filter(result_type == "interaction") %>% filter(Pr...t.. < 0.05) #1569
write.table(interaction_sign, 'interaction_sign_lm.txt', sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)


##show how they overlap if interested
# venn.diagram(
#   x = list(smurf = smurf_sign$gene_name, time = time_sign$gene_name, interaction = interaction_sign$gene_name),
#   filename = "venn_diagram.pdf",
#   col = "transparent",
#   fill = c("deepskyblue1", "red", "grey"),
#   alpha = 0.5,
#   label.col = "black",
#   cex = 1.5,
#   fontface = "bold",
#   main = "Venn Diagram"
# )

###RUN FGSEA 

#We ran fgsea analysis on the different list of genes to investigate biological processes enriched (if any)
#analysis were run 10 times and categories presents 100% of the times were selected for consideration

####run fgsea on the different categories: start with smurfs
gene_list <- as.numeric(smurf_sign$Estimate)
names(gene_list) <- as.character(smurf_sign$gene_name)
gene_list = gene_list[!is.na(names(gene_list))] %>% .[!duplicated(names(.))] %>% sort(., decreasing = TRUE)

require(clusterProfiler)
#run the analysis
n = 10 #number of times repeating
P = 100 #percentage for selection
pval = 0.05 #threshold to apply
gse_go = vector("list", n)
organism = "org.Dm.eg.db"
library(organism, character.only = TRUE) #library is already loaded
set.seed(123)
for (i in 1:length(gse_go)) {
  gse_go[[i]] <- gseGO(geneList= gene_list,
                       ont ="BP", 
                       keyType = "ENSEMBL", 
                       nPerm = 15000, 
                       minGSSize = 10, 
                       maxGSSize = 600, 
                       pvalueCutoff = 0.05, #change the significance level
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "fdr") # %>% as.data.frame()
}


df_gse = lapply(gse_go, function(x) as.data.frame(x))
#extract the categories of interest
categories = lapply(df_gse, function(x) x$Description)
filtering_vector = tibble(category = unlist(categories)) %>% #unlist the list and making it a tibble
  group_by(category) %>% 
  summarise(count = n()) %>% #group by category and generate a new column counting how many times a category occurs
  mutate(percentage_presence = count/n*100) %>% #adding a column with the persentage 
  #arrange(desc(count)) %>% #adding this only if we want to visualize
  filter(percentage_presence >= P)

df_gse[[1]]$Description %in% filtering_vector$category %>% sum
plot_gse = gse_go[[1]]
emapplot(plot_gse, showCategory = filtering_vector$category, color = "NES", 
         layout = "fr", cex_label_category = 0.1) + scale_color_gradient(low = "blue" , high = "red" )

####genes with time
n = 10 #number of times repeating
P = 100 #percentage for selection
pval = 0.05 #threshold to apply
gse_go2 = vector("list", n)
organism = "org.Dm.eg.db"
library(organism, character.only = TRUE) #library is already loaded
set.seed(123)
for (i in 1:length(gse_go)) {
  gse_go2[[i]] <- gseGO(geneList= gene_list2,
                        ont ="BP", 
                        keyType = "ENSEMBL", 
                        nPerm = 15000, 
                        minGSSize = 10, 
                        maxGSSize = 600, 
                        pvalueCutoff = 0.05, #change the significance level
                        verbose = TRUE, 
                        OrgDb = organism, 
                        pAdjustMethod = "fdr") # %>% as.data.frame()
}


df_gse2 = lapply(gse_go2, function(x) as.data.frame(x))
#extract the categories of interest
categories2 = lapply(df_gse2, function(x) x$Description)
filtering_vector2 = tibble(category = unlist(categories2)) %>% #unlist the list and making it a tibble
  group_by(category) %>% 
  summarise(count = n()) %>% #group by category and generate a new column counting how many times a category occurs
  mutate(percentage_presence = count/n*100) %>% #adding a column with the persentage 
  #arrange(desc(count)) %>% #adding this only if we want to visualize
  filter(percentage_presence >= P)

df_gse2[[1]]$Description %in% filtering_vector2$category %>% sum
plot_gse2 = gse_go2[[1]]
emapplot(plot_gse2, showCategory = filtering_vector2$category, color = "NES", 
         layout = "fr", cex_label_category = 0.1) + scale_color_gradient(low = "blue" , high = "red" )


####do the same for the ones that are interaction only
####run fgsea on the different categories: start with smurfs
gene_list3 <- as.numeric(interaction_sign$Estimate)
names(gene_list3) <- as.character(interaction_sign$gene_name)
gene_list3 = gene_list3[!is.na(names(gene_list3))] %>% .[!duplicated(names(.))] %>% sort(., decreasing = TRUE)

require(clusterProfiler)
#run the analysis
n = 10 #number of times repeating
P = 100 #percentage for selection
pval = 0.05 #threshold to apply
gse_go3 = vector("list", n)
organism = "org.Dm.eg.db"
library(organism, character.only = TRUE) #library is already loaded
set.seed(123)
for (i in 1:length(gse_go3)) {
  gse_go3[[i]] <- gseGO(geneList= gene_list3,
                       ont ="BP", 
                       keyType = "ENSEMBL", 
                       nPerm = 15000, 
                       minGSSize = 10, 
                       maxGSSize = 600, 
                       pvalueCutoff = 0.05, #change the significance level
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "fdr") # %>% as.data.frame()
}


df_gse3 = lapply(gse_go3, function(x) as.data.frame(x))
#extract the categories of interest
categories3 = lapply(df_gse3, function(x) x$Description)
filtering_vector3 = tibble(category = unlist(categories3)) %>% #unlist the list and making it a tibble
  group_by(category) %>% 
  summarise(count = n()) %>% #group by category and generate a new column counting how many times a category occurs
  mutate(percentage_presence = count/n*100) %>% #adding a column with the persentage 
  #arrange(desc(count)) %>% #adding this only if we want to visualize
  filter(percentage_presence >= P)

df_gse3[[1]]$Description %in% filtering_vector2$category %>% sum
plot_gse3 = gse_go3[[1]]
emapplot(plot_gse3, showCategory = filtering_vector3$category, color = "NES", 
         layout = "fr", cex_label_category = 0.1) + scale_color_gradient(low = "blue" , high = "red" )

######We have three "filtering_vectors", representing the categories enriched in the three list of genes 
###check how the filtering vector overlap with venn diagram
# filtering_vector
# filtering_vector2
# filtering_vecto3
# venn.diagram(
#   x = list(smurf = filtering_vector$category, time = filtering_vector2$category, interaction = filtering_vector3$category),
#   filename = "venn_diagram_cat.pdf",
#   col = "transparent",
#   fill = c("deepskyblue1", "red", "grey"),
#   alpha = 0.5,
#   label.col = "black",
#   cex = 1.5,
#   fontface = "bold",
#   main = "Venn Diagram"
# )


#####################HEATMAP TO SHOW THE OVERLAP 
###check the way the categories distribute
smurf_gsea_res = df_gse[[1]] %>% filter(Description %in% filtering_vector$category) %>% select(c("Description", "NES","p.adjust"))
time_gsea_res = df_gse2[[1]] %>% filter(Description %in% filtering_vector2$category) %>% select(c("Description", "NES","p.adjust"))
interaction_gsea_res = df_gse3[[1]] %>% filter(Description %in% filtering_vector3$category) %>% select(c("Description", "NES","p.adjust"))

require(reshape2)
smurf_gsea_res$df_name = "smurf"
time_gsea_res$df_name = 'time'
interaction_gsea_res$df_name = "interaction"
gsea_all <- rbind(smurf_gsea_res, time_gsea_res, interaction_gsea_res) %>% select(-p.adjust)
gsea_wide <- gsea_all %>% spread(df_name, NES) 
gsea_long <- melt(gsea_wide, id.vars = "Description") %>% distinct()

#####correct order in the heatmap
ordering <- gsea_long %>%
  filter(variable == "smurf") %>%
  arrange(value) %>%
  pull(Description)

# Convert Description to a factor with levels ordered by the ordering vector
gsea_long$Description <- factor(gsea_long$Description, levels = ordering)
# Create the plot
ggplot(gsea_long,
       aes(x=variable, y=Description, fill=value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "grey", 
                       midpoint = 0, 
                       name = "NES", na.value = "white") +
  theme_minimal() +
  xlab("Name") +
  ylab("DataFrame") +
  theme(axis.text.y = element_text(size = 8, vjust = 0.5, hjust=1))
ggsave("../figures/extra_figure_heatmap_lm_intercation.pdf")

##############################################################################################
#extra possible checkanalysis
###overlap with the differentially expressed genes - less interesting analysis in this case compared to gsea
smurf_degs = read.table("../data/differential_expression/res_DEGS_SvsNS.txt", header = TRUE) %>% filter(padj < 0.05)
non_smurf_degs = read.table("../data/differential_expression/res_DEGS_NS_40days_20days.txt", header = TRUE) %>% filter(padj < 0.05)
smurf_old_degs = read.table("../data/differential_expression/res_DEGS_S_40days_20days.txt", header = TRUE) %>% filter(padj < 0.05)

library(upsetR)
listInput = list(s_degs = rownames(smurf_degs), ns_degs = rownames(non_smurf_degs), old_s_degs = rownames(smurf_old_degs),
                  lm_s = smurf_sign$gene_name, lm_time = time_sign$gene_name, lm_interaction = interaction_sign$gene_name)

upset(fromList(listInput), nsets = 6, order.by = "freq")
