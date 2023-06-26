#set the working directory if needed
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #setting current working directory
#load here libraries that are needed all along
library(org.Dm.eg.db); library(tidyverse)

##load the files
res_s = read.table("../data/differential_expression/res_DEGS_SvsNS.txt", header = TRUE, stringsAsFactors = FALSE)
res_ns = read.table("../data/differential_expression/res_DEGS_NS_40days_20days.txt", header = TRUE, stringsAsFactors = FALSE)

sign_s = res_s %>% filter(padj < 0.05)

###plot all the genes and in different colors the different cases
corr_plot_file = data.frame(flybase = rownames(res_s), symbol = res_ns$symbol, log2FC_s = res_s$log2FoldChange, log2FC_ns = res_ns$log2FoldChange,
                            padj_s = res_s$padj, padj_ns = res_ns$padj) %>% mutate(sign = ifelse(padj_s < 0.05 & padj_ns < 0.05, "both", ifelse(padj_s < 0.05 & padj_ns > 0.05, "smurf", ifelse(padj_s > 0.05 & padj_ns < 0.05, "non_smurf", "none")))) %>% 
                            mutate(alpha = case_when(sign == "none" ~ 0.05, 
                                                     sign == "non_smurf" ~ 0.4,
                                                     sign == "smurf" ~ 0.4, 
                                                     TRUE ~ 0.8))
  
  
####
ggplot(corr_plot_file %>% na.omit(), aes(log2FC_s, log2FC_ns)) +
  geom_point(aes(color = sign, alpha = alpha)) +
  scale_color_manual(values = c( "none" = "lightgrey", "both" = "pink", "smurf" = "deepskyblue1", "non_smurf" = "red"),  
                                 labels = c("not significant", "significant in S & NS", "significant in S", "significant in NS"),
                                 name = "") +
  geom_smooth(aes(color = sign), method = "lm", se = TRUE, alpha = 0.3) +
  theme_minimal() + guides(alpha ="none") +
  xlab("log2FC S/NS") + ylab("log2FC oldNS/youngNS")
ggsave("../figures/correlation_log2FC_smurf_non_smurf.png", width = 7, height = 5)


# Subset the data to include only the relevant variables and remove NA values
subset_data <- na.omit(corr_plot_file %>% select(log2FC_s, log2FC_ns, sign))

# Calculate correlations by sign level
correlations <- subset_data %>% group_by(sign) %>%
  summarize(correlation = cor(log2FC_s, log2FC_ns, method=c("spearman")), count = n())

# Print the correlations
print(correlations)




# packageurl <- "http://cran.r-project.org/src/contrib/Archive/mgcv/mgcv_1.8-30.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
