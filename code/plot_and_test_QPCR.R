setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse); library(readxl)

####AFD1

###open from the excel file --> plotting the 2-DDCT and then run wilcoxon test

data_test_adf1 = read_xlsx("../data/qpcr_excel/qPCR_adf1.xlsx", sheet = 7)
colnames(data_test_adf1) = c("RU", "value")
data_plot_adf1 = read_xlsx("../data/qpcr_excel/qPCR_adf1.xlsx", sheet = 8)
colnames(data_plot_adf1) = c("RU", "Av_2DDCT", "STDEV")
data_plot_adf1$RU = factor(data_plot_adf1$RU, levels = as.character(data_plot_adf1$RU) )

##barplot
ggplot(data_plot_adf1, aes(RU, Av_2DDCT)) + 
  geom_bar(stat = 'identity', fill = c("black", "brown4", "firebrick1","coral", "lightpink1"), alpha = 0.8) + 
  geom_errorbar(aes(ymin= Av_2DDCT - STDEV, ymax= Av_2DDCT + STDEV), width=0.1, position=position_dodge(.9)) +
  scale_y_continuous(breaks = round(seq(0, 1.2, by = 0.1),1)) + theme_light() +
  xlab("RU condition") + ylab("Relative expression") +
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
ggsave("../figures/supplementary/Fig_S15_adf1_qpcr_res.pdf", width = 5, height =4)

##test 
wilcox.test(data_test_adf1 %>% filter(RU == "RU0") %>% select(value) %>% unlist, 
              data_test_adf1 %>% filter(RU == "RU10") %>% select(value) %>% unlist) #pval = 0.2

wilcox.test(data_test_adf1 %>% filter(RU == "RU0") %>% select(value) %>% unlist, 
            data_test_adf1 %>% filter(RU == "RU50") %>% select(value) %>% unlist) #pval = 0.7

wilcox.test(data_test_adf1 %>% filter(RU == "RU0") %>% select(value) %>% unlist, 
            data_test_adf1 %>% filter(RU == "RU100") %>% select(value) %>% unlist) #pval = 0.

wilcox.test(data_test_adf1 %>% filter(RU == "RU0") %>% select(value) %>% unlist, 
            data_test_adf1 %>% filter(RU == "RU200") %>% select(value) %>% unlist) #pval =0.4


###CG4360
data_test_cg = read_xlsx("../data/qpcr_excel/qPCR_CG4360.xlsx", sheet = 5)
colnames(data_test_cg) = c("RU", "value")
#data if we want to keep the outlier for the analysis
# data_plot_cg = read_xlsx("../data/qpcr_excel/qPCR_CG4360.xlsx", sheet = 6)
# colnames(data_plot_cg) = c("RU", "Av_2DDCT", "STDEV")
# data_plot_cg$RU = factor(data_plot_cg$RU, levels = as.character(data_plot_cg$RU) )

###Remove the control outlier, probably due to manipulation error
data_test_cg_without_outlier = data_test_cg[-2,]
data_plot_cg= data_test_cg_without_outlier %>% group_by(RU) %>% summarise(mean = mean(value), stdev = sd(value))
data_plot_cg$RU = factor(data_plot_cg$RU, levels = c("RU0", "RU10", "RU50", "RU100", "RU200"))
##barplot
ggplot(data_plot_cg, aes(RU, mean)) + 
  geom_bar(stat = 'identity', fill = c("black", "brown4", "firebrick1","coral", "lightpink1"), alpha = 0.8) + 
  geom_errorbar(aes(ymin= mean - stdev, ymax= mean + stdev), width=0.1, position=position_dodge(.9)) +
  scale_y_continuous(breaks = round(seq(0, 1.5, by = 0.1),1)) + theme_light() +
  xlab("RU condition") + ylab("Relative expression") +
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
ggsave("../figures/supplementary/Fig_S15_cg4360_qpcr_res.pdf", width = 5, height =4)


##test 
wilcox.test(data_test_cg_without_outlier %>% filter(RU == "RU0") %>% select(value) %>% unlist, 
            data_test_cg_without_outlier %>% filter(RU == "RU10") %>% select(value) %>% unlist) #pval = 0.05

wilcox.test(data_test_cg %>% filter(RU == "RU0") %>% select(value) %>% unlist, 
            data_test_cg %>% filter(RU == "RU50") %>% select(value) %>% unlist) #pval = 0.05

wilcox.test(data_test_cg %>% filter(RU == "RU0") %>% select(value) %>% unlist, 
            data_test_cg %>% filter(RU == "RU100") %>% select(value) %>% unlist) #pval = 0.05

wilcox.test(data_test_cg %>% filter(RU == "RU0") %>% select(value) %>% unlist, 
            data_test_cg %>% filter(RU == "RU200") %>% select(value) %>% unlist) #pval =0.05


###please note that, even if running the analysis with the control outlier, the difference in relative expression between RU0 and treated samples is still significant


###############TRL
data_test_trl = read_xlsx("../data/qpcr_excel/qPCR_Trl.xlsx", sheet = 4)
colnames(data_test_trl) = c("RU", "value")
data_plot_trl = read_xlsx("../data/qpcr_excel/qPCR_Trl.xlsx", sheet = 5)
colnames(data_plot_trl) = c("RU", "Av_2DDCT", "STDEV")
data_plot_trl$RU = factor(data_plot_trl$RU, levels = as.character(data_plot_trl$RU) )

##barplot
ggplot(data_plot_trl, aes(RU, Av_2DDCT)) + 
  geom_bar(stat = 'identity', fill = c("black", "brown4", "firebrick1","coral", "lightpink1"), alpha = 0.8) + 
  geom_errorbar(aes(ymin= Av_2DDCT - STDEV, ymax= Av_2DDCT + STDEV), width=0.1, position=position_dodge(.9)) +
  scale_y_continuous(breaks = round(seq(0, 1.5, by = 0.1),1)) + theme_light() +
  xlab("RU condition") + ylab("Relative expression") +
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
ggsave("../figures/supplementary/Fig_S15_trl_qpcr_res.pdf", width = 5, height =4)

##test 
wilcox.test(data_test_trl %>% filter(RU == "RU0") %>% select(value) %>% unlist, 
            data_test_trl %>% filter(RU == "RU10") %>% select(value) %>% unlist) #pval = 0.8

wilcox.test(data_test_trl %>% filter(RU == "RU0") %>% select(value) %>% unlist, 
            data_test_trl %>% filter(RU == "RU50") %>% select(value) %>% unlist) #pval = 0.1

wilcox.test(data_test_trl %>% filter(RU == "RU0") %>% select(value) %>% unlist, 
            data_test_trl %>% filter(RU == "RU100") %>% select(value) %>% unlist) #pval = 0.1

wilcox.test(data_test_trl %>% filter(RU == "RU0") %>% select(value) %>% unlist, 
            data_test_trl %>% filter(RU == "RU200") %>% select(value) %>% unlist) #pval =0.1
