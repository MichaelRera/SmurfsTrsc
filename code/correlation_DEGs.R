#set the working directory if needed
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #setting current working directory
#load here libraries that are needed all along
library(org.Dm.eg.db); library(tidyverse)

##load the files
res_s = read.table("../data/differential_expression/res_DEGS_SvsNS.txt", header = TRUE, stringsAsFactors = FALSE)
res_ns = read.table("../data/differential_expression/res_DEGS_NS_40days_20days.txt", header = TRUE, stringsAsFactors = FALSE)

#take only the genes that are differentially expressed in Smurfs
degs_s = res_s %>% filter(padj < 0.05)
degs_s_in_ns = res_ns %>% filter(rownames(.) %in% rownames(degs_s)) %>% mutate(sign_ns = ifelse(padj < 0.05, TRUE, FALSE))
#create the file for the plot
corr_file1 = degs_s_in_ns %>% select(log2FoldChange, symbol, sign_ns) %>% mutate(log2FC_s = degs_s$log2FoldChange)
#plot
p = ggplot(corr_file1 %>% na.omit(), aes(x = log2FoldChange, y = log2FC_s)) + geom_point(aes(col = sign_ns), alpha = 0.5) + 
  geom_smooth(aes(group = sign_ns, color = sign_ns), method = "lm")

# Subset the data where sign_ns is FALSE
subset_data <- corr_file1 %>% na.omit() %>% filter(sign_ns == FALSE)
# Fit linear regression model for the subset
lm_model_false <- lm(log2FC_s ~ log2FoldChange, data = subset_data)
# Extract the coefficients
coefficients_false <- coef(lm_model_false)
# Print the coefficients for sign_ns == FALSE
print(coefficients_false)


# Subset the data where sign_ns is TRUE
subset_data2 <- corr_file1 %>% na.omit() %>% filter(sign_ns == FALSE)
# Fit linear regression model for the subset
lm_model_TRUE <- lm(log2FC_s ~ log2FoldChange, data = subset_data2)
# Extract the coefficients
coefficients_TRUE <- coef(lm_model_TRUE)
# Print the coefficients for sign_ns == TRUE
print(coefficients_TRUE)


# packageurl <- "http://cran.r-project.org/src/contrib/Archive/mgcv/mgcv_1.8-30.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")
