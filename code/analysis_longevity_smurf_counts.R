setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse); library(readxl)




#################TRL ADULT ONLY##########################

##open the smurfs file 
smurf_ao_trl <- read_excel("../data/longevity/smurf_count_data/AdultOnly_trl_FEMALES_SMURFS.xlsm", range = "A9:CN33", col_names = FALSE) %>% 
  select(-c( "...4", "...5")) 
colnames(smurf_ao_trl) <- c("gene", "concentration", "vial", seq(0, ncol(smurf_ao_trl)-3, by = 1)) 
smurf_ao_trl <- smurf_ao_trl %>% select_if(~ !any(is.na(.)))
#open the file
AO_trl <- read_excel("../data/longevity/AdultOnly_Trl_FEMALES.xlsm", range = "A9:CS33", col_names = FALSE) %>% 
  select(-c( "...4", "...5")) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0)))
colnames(AO_trl) <- c("gene", "concentration", "vial", "initial_number", seq(1, ncol(AO_trl)-3, by = 1)) #the initial number is the 0 column for the other
AO_trl <- AO_trl %>% filter(gene %in% smurf_ao_trl$gene) #filter genes that are not needed 

days_sampling_ao <- colnames(smurf_ao_trl)[4:length(colnames(smurf_ao_trl))]

#computing the remaining number of flies at each smurf assay point 

nflies_9 = AO_trl %>% select(matches("init"):matches("8")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies9 = initial_number -(rowsum - initial_number)) %>% select(nflies9)
#you can ignore the warning message here, it is doing the correct thing as I need column 8
nflies_36 = AO_trl %>% select(matches("init"):matches("35")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies36 = initial_number -(rowsum - initial_number)) %>% select(nflies36)
nflies_50 = AO_trl %>% select(matches("init"):matches("49")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies50 = initial_number -(rowsum - initial_number)) %>% select(nflies50)
nflies_57 = AO_trl %>% select(matches("init"):matches("56")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies57 = initial_number -(rowsum - initial_number)) %>% select(nflies57)
nflies_64 = AO_trl %>% select(matches("init"):matches("63")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies64 = initial_number -(rowsum - initial_number)) %>% select(nflies64)
nflies_71 = AO_trl %>% select(matches("init"):matches("70")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies71 = initial_number -(rowsum - initial_number)) %>% select(nflies71)
nflies_78 = AO_trl %>% select(matches("init"):matches("77")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies78 = initial_number -(rowsum - initial_number)) %>% select(nflies78)
nflies_84 = AO_trl %>% select(matches("init"):matches("83")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies84 = initial_number -(rowsum - initial_number)) %>% select(nflies84)



###add those infos to the smurf file
smurfs_ao_modified = smurf_ao_trl %>% cbind(data.frame(nflies_9, nflies_36, nflies_50, nflies_57, nflies_64, nflies_71, nflies_78, nflies_84))
smurf_ao_prop = smurfs_ao_modified %>%
  mutate(p9 = ifelse(nflies9 > 0, smurfs_ao_modified$`9`/nflies9, NA)) %>% 
  mutate(p36 = ifelse(nflies36 > 0, smurfs_ao_modified$`36`/nflies36, NA)) %>% 
  mutate(p50 = ifelse(nflies50 > 0, smurfs_ao_modified$`50`/nflies50, NA)) %>% 
  mutate(p57 = ifelse(nflies57 > 0, smurfs_ao_modified$`57`/nflies57, NA)) %>% 
  mutate(p64 = ifelse(nflies64 > 0, smurfs_ao_modified$`64`/nflies64, NA)) %>% 
  mutate(p71 = ifelse(nflies71 > 0, smurfs_ao_modified$`71`/nflies71, NA)) %>% 
  mutate(p78 = ifelse(nflies78 > 0, smurfs_ao_modified$`78`/nflies78, NA)) %>% 
  mutate(p84 = ifelse(nflies84 > 0, smurfs_ao_modified$`84`/nflies84, NA))


smurf_ao_for_plot = smurf_ao_prop %>%  select(gene, concentration, vial, p9, p36, p50, p57, p64, p71, p78, p84)
colnames(smurf_ao_for_plot) = c("gene", "concentration", "vial", "9", "36", "50", "57", "64", "71", "78", "84")
smurf_ao_for_plot =  smurf_ao_for_plot %>%  gather(., key = 'day', value = "prop_smurf", -gene, -concentration, - vial )
smurf_ao_for_plot$concentration = factor(smurf_ao_for_plot$concentration, levels = c("RU0", "RU10", "RU50", "RU100", "RU200"))


ggplot(smurf_ao_for_plot %>% filter(concentration %in% c("RU0", "RU50")) %>% filter(day !=84), #day 84 removed as too little flies were left and data might not be representative
       aes(as.numeric(day), prop_smurf)) +
  # geom_boxplot(alpha = 0.7) +
  geom_point(aes(col = concentration)) + geom_smooth(aes(col = concentration,fill = concentration), method = "lm") +
  scale_color_manual(values = c("black", "red3")) +
  scale_fill_manual(values = c("grey", "red3")) +
  xlab("days") + ylab("proportion of Smurfs")  +   ylim(0,0.6) + theme_linedraw() +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))

ggsave("../figures/longevity_plots/smurf_counts/Fig7bii_TRL_AD_RU0RU50.pdf", width = 6, height = 4)

###to analyze the effect of time only
lm(prop_smurf ~ as.numeric(day), 
   data = trl_ao_for_plot %>% filter(concentration %in% c("RU0")) %>% 
     filter(day != 84)) %>% summary
lm(prop_smurf ~ as.numeric(day), 
   data = trl_ao_for_plot %>% filter(concentration %in% c("RU50")) %>% 
     filter(day != 84)) %>% summary

###to analyze the effect of the drug treatment
lm(prop_smurf ~ as.numeric(day)  + concentration + concentration*as.numeric(day), 
   data = trl_ao_for_plot %>% filter(concentration %in% c("RU0", "RU50")) %>% 
     filter(day != 84)) %>% summary




#################################ADF1 ADULT ONLY #############################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse); library(readxl)

################## --------  ADULTONLY ------------ ###########

##open the smurfs file 
smurfs_ao_adf1 <- read_excel("../data/longevity/smurf_count_data/AO_adf1_females_smurfs.xlsm", range = "A9:CO33", col_names = FALSE) %>% 
  select(-c( "...4", "...5", "...6")) 
colnames(smurfs_ao_adf1) <- c("gene", "concentration", "vial", seq(0, ncol(smurfs_ao_adf1)-3, by = 1)) 
smurfs_ao_adf1 <- smurfs_ao_adf1 %>% select_if(~!any(is.na(.)))
#open the file
AO_adf1 <- read_excel("../data/longevity/AO_4278ADF1_females.xlsm", range = "A9:CZ33", col_names = FALSE) %>% 
  select(-c("...3", "...4", "...5")) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) 
colnames(AO_adf1) <- c("gene", "concentration", "initial_number", seq(1, ncol(AO_adf1)-3, by = 1)) #the initial number is the 0 column for the other
#AO_val4 <- AO_val4 %>% filter(gene %in% smurfs_ao_val4$gene) #filter genes that are not needed 

days_sampling_ao <- colnames(smurfs_ao_adf1)[4:length(colnames(smurfs_ao_adf1))]

##############################################
#computing the remaining number of flies at each smurf assay point 

nflies_16 = AO_adf1 %>% select(matches("init"):matches("15")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies16 = initial_number -(rowsum - initial_number)) %>% select(nflies16)
nflies_28 = AO_adf1 %>% select(matches("init"):matches("27")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies28 = initial_number -(rowsum - initial_number)) %>% select(nflies28)
nflies_35 = AO_adf1 %>% select(matches("init"):matches("34")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies35 = initial_number -(rowsum - initial_number)) %>% select(nflies35)
nflies_42 = AO_adf1 %>% select(matches("init"):matches("41")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies42 = initial_number -(rowsum - initial_number)) %>% select(nflies42)
nflies_49 = AO_adf1 %>% select(matches("init"):matches("48")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies49 = initial_number -(rowsum - initial_number)) %>% select(nflies49)
nflies_56 = AO_adf1 %>% select(matches("init"):matches("55")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies56 = initial_number -(rowsum - initial_number)) %>% select(nflies56)
nflies_63 = AO_adf1 %>% select(matches("init"):matches("62")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies63 = initial_number -(rowsum - initial_number)) %>% select(nflies63)
nflies_70 = AO_adf1 %>% select(matches("init"):matches("69")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies70 = initial_number -(rowsum - initial_number)) %>% select(nflies70)
nflies_77 = AO_adf1 %>% select(matches("init"):matches("76")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies77 = initial_number -(rowsum - initial_number)) %>% select(nflies77)
nflies_86 = AO_adf1 %>% select(matches("init"):matches("85")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies86 = initial_number -(rowsum - initial_number)) %>% select(nflies86)



###add those infos to the smurf file
smurfs_ao_modified = smurfs_ao_adf1 %>% cbind(data.frame(nflies_16, nflies_28, nflies_35, nflies_42, nflies_49, nflies_56, nflies_63, nflies_70, nflies_77, nflies_86))
smurf_ao_prop = smurfs_ao_modified %>%
  mutate(p16 = ifelse(nflies16 > 0, smurfs_ao_modified$`16`/nflies16, NA)) %>% 
  mutate(p28 = ifelse(nflies28 > 0, smurfs_ao_modified$`28`/nflies28, NA)) %>% 
  mutate(p35 = ifelse(nflies35 > 0, smurfs_ao_modified$`35`/nflies35, NA)) %>% 
  mutate(p42 = ifelse(nflies42 > 0, smurfs_ao_modified$`42`/nflies42, NA)) %>% 
  mutate(p49 = ifelse(nflies49 > 0, smurfs_ao_modified$`49`/nflies49, NA)) %>% 
  mutate(p56 = ifelse(nflies56 > 0, smurfs_ao_modified$`56`/nflies56, NA )) %>% 
  mutate(p63 = ifelse(nflies63 > 0, smurfs_ao_modified$`63`/nflies63, NA)) %>% 
  mutate(p70 = ifelse(nflies70 > 0, smurfs_ao_modified$`70`/nflies70, NA)) %>% 
  mutate(p77 = ifelse(nflies77 > 0, smurfs_ao_modified$`77`/nflies77, NA)) %>% 
  mutate(p86 = ifelse(nflies86 > 0, smurfs_ao_modified$`86`/nflies86, NA)) 




smurf_ao_for_plot = smurf_ao_prop %>%  select(gene, concentration, vial, p16, p28, p35, p42, p49, p56, p63, p70, p77, p86)
colnames(smurf_ao_for_plot) = c("gene", "concentration", "vial", days_sampling_ao)
smurf_ao_for_plot =  smurf_ao_for_plot %>%  gather(., key = 'day', value = "prop_smurf", -gene, -concentration, - vial )
smurf_ao_for_plot$concentration = factor(smurf_ao_for_plot$concentration, levels = c("RU0", "RU10", "RU50", "RU100", "RU200"))

ggplot(smurf_ao_for_plot %>% filter(concentration %in% c("RU0", "RU50")) %>% 
         filter(day != 86), aes(as.numeric(day), prop_smurf)) + 
  scale_fill_manual(values = c("grey", "red3")) +
  #geom_boxplot(alpha = 0.7) +
  geom_point(aes(col = concentration)) +
  scale_color_manual(values = c("black", "red3")) + 
  geom_smooth(aes(color = concentration, fill = concentration), method = "lm") +
  xlab("days") + ylab("proportion of Smurfs")  +   ylim(0, 0.75) + xlim(0,80) + theme_linedraw() + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
ggsave("../figures/longevity_plots/smurf_counts/Fig7bii_ADF1_AD_RU0RU50.pdf", width = 6, height = 4)


###for the slope and time only

lm(prop_smurf ~ as.numeric(day), 
   data = smurf_ao_for_plot %>% filter(concentration %in% c("RU0")) %>% 
     filter(day != 86)) %>% summary
lm(prop_smurf ~ as.numeric(day), 
   data = smurf_ao_for_plot %>% filter(concentration %in% c("RU50")) %>% 
     filter(day != 86)) %>% summary
#for the drug interaction dose
lm(prop_smurf ~ as.numeric(day)  + concentration + concentration*as.numeric(day), 
   data = smurf_ao_for_plot %>% filter(concentration %in% c("RU0", "RU50")) %>% 
     filter(day != 86)) %>% summary


#pvalue = 0.05



##################################CG4360 WHOLE LIFE##############################""""


####WL Smurfs analysis
##open the smurfs file 
smurfs_wl_cg4360 <- read_excel("../data/longevity/smurf_count_data/WholeLife_cg4360_FEMALES_SMURFS.xlsm", range = "A9:CM33", col_names = FALSE) %>% 
  select(-c( "...4", "...5", "...6")) 
colnames(smurfs_wl_cg4360) <- c("gene", "concentration", "vial", seq(0, ncol(smurfs_wl_cg4360)-3, by = 1)) 
smurfs_wl_cg4360 <- smurfs_wl_cg4360 %>% select_if(~!any(is.na(.))) %>% filter(gene == "51813 CG4360")
#open the file
WL_cg4360 <- read_excel("../data/longevity/WholeLife_cg4360_FEMALES.xlsm", range = "A9:DF33", col_names = FALSE) %>% 
  select(-c("...3", "...4", "...5")) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) 
colnames(WL_cg4360) <- c("gene", "concentration", "initial_number", seq(1, ncol(WL_cg4360)-3, by = 1))

#AO_val4 <- AO_val4 %>% filter(gene %in% smurfs_ao_val4$gene) #filter genes that are not needed 

days_sampling_wl <- colnames(smurfs_wl_cg4360)[4:length(colnames(smurfs_wl_cg4360))]
#57 is for setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#computing the remaining number of flies at each smurf assay point 

nflies_21 = WL_cg4360 %>% select(matches("init"):matches("20")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies21 = initial_number -(rowsum - initial_number)) %>% select(nflies21)
nflies_37 = WL_cg4360 %>% select(matches("init"):matches("36")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies37 = initial_number -(rowsum - initial_number)) %>% select(nflies37)
nflies_49 = WL_cg4360 %>% select(matches("init"):matches("48")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies49 = initial_number -(rowsum - initial_number)) %>% select(nflies49)
nflies_56 = WL_cg4360 %>% select(matches("init"):matches("55")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies56 = initial_number -(rowsum - initial_number)) %>% select(nflies56)
nflies_63 = WL_cg4360 %>% select(matches("init"):matches("62")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies63 = initial_number -(rowsum - initial_number)) %>% select(nflies63)
nflies_70 = WL_cg4360 %>% select(matches("init"):matches("69")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies70 = initial_number -(rowsum - initial_number)) %>% select(nflies70)
nflies_77 = WL_cg4360 %>% select(matches("init"):matches("76")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies77 = initial_number -(rowsum - initial_number)) %>% select(nflies77)
nflies_84 = WL_cg4360 %>% select(matches("init"):matches("83")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies84 = initial_number -(rowsum - initial_number)) %>% select(nflies84)

smurfs_wl_modified = smurfs_wl_cg4360 %>% cbind(data.frame(nflies_21, nflies_37, nflies_49, nflies_56, nflies_63, nflies_70, nflies_77, nflies_84))
smurf_wl_prop = smurfs_wl_modified %>%
  mutate(p21 = ifelse(nflies21 > 0, smurfs_wl_modified$`21`/nflies21, NA)) %>% 
  mutate(p37 = ifelse(nflies37 > 0, smurfs_wl_modified$`37`/nflies37, NA)) %>% 
  mutate(p49 = ifelse(nflies49 > 0, smurfs_wl_modified$`49`/nflies49, NA)) %>% 
  mutate(p56 = ifelse(nflies56 > 0, smurfs_wl_modified$`56`/nflies56, NA)) %>% 
  mutate(p63 = ifelse(nflies63 > 0, smurfs_wl_modified$`63`/nflies63, NA)) %>% 
  mutate(p70 = ifelse(nflies70 > 0, smurfs_wl_modified$`70`/nflies70, NA )) %>% 
  mutate(p77 = ifelse(nflies77 > 0, smurfs_wl_modified$`77`/nflies77, NA)) %>% 
  mutate(p84 = ifelse(nflies84 > 0, smurfs_wl_modified$`84`/nflies84, NA)) 

smurf_wl_for_plot = smurf_wl_prop %>%  select(gene, concentration, vial, p21, p37, p49, p56, p63, p70, p77, p84)
colnames(smurf_wl_for_plot) = c("gene", "concentration", "vial", days_sampling_wl)

smurf_wl_plotnflies = smurf_wl_prop %>%  select(gene, concentration, vial, nflies21, nflies37, nflies49, nflies56, nflies63, nflies70, nflies77, nflies84)
colnames(smurf_wl_plotnflies) = c("gene", "concentration", "vial", days_sampling_wl)
smurf_wl_for_plotnflies2 =  smurf_wl_plotnflies %>%  gather(., key = 'day', value = "alive_flies", -gene, -concentration, - vial )

s2 = smurf_wl_for_plotnflies2 %>% select(-c("gene", "vial")) %>% group_by(day,concentration) %>% summarise(count = sum(alive_flies))


###plots
smurf_wl_for_plot =  smurf_wl_for_plot %>%  gather(., key = 'day', value = "prop_smurf", -gene, -concentration, - vial )
smurf_wl_for_plot$concentration = factor(smurf_wl_for_plot$concentration, levels = c("RU0", "RU10", "RU50", "RU100", "RU200"))

ggplot(smurf_wl_for_plot %>% filter(concentration %in% c("RU0", "RU10")) %>%
         filter(day != 84), aes(as.numeric(day), prop_smurf)) + 
  #scale_fill_manual(values = c("grey", "red3")) +
  #geom_boxplot(alpha = 0.7) + #geom_jitter(aes(color = concentration), size = 2, width = 1.5) +  
  geom_point(aes(col = concentration)) + geom_smooth(aes(color = concentration, fill = concentration), method = 'lm') +
  scale_fill_manual(values = c("grey", "red3")) +
  scale_color_manual(values = c("black", "red3")) +
  xlab("days") + ylab("proportion of Smurfs")  +   ylim(0,0.6) + xlim(0,80) + theme_linedraw() +
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
ggsave("../figures/longevity_plots/smurf_counts/Fig7bii_CG4360_WL_RU0RU10.pdf", width = 6, height = 4)

###test independently for the two
lm(prop_smurf ~ as.numeric(day), 
   data = smurf_wl_for_plot %>% 
     filter(day != 84)%>% filter(concentration %in% c("RU0")) ) %>% summary
lm(prop_smurf ~ as.numeric(day), 
   data = smurf_wl_for_plot %>% 
     filter(day != 84)%>% filter(concentration %in% c("RU10")) ) %>% summary



###test for interaction drug dose
lm(prop_smurf ~ as.numeric(day)  + concentration + concentration*as.numeric(day), 
   data = smurf_wl_for_plot %>% 
     filter(day != 84)%>% filter(concentration %in% c("RU0", "RU10")) ) %>% summary

#pvalue = 0.012


#######################GENE NOT DISPLAYING LONGEVITY CHANGES  ---- CG4360 AD THIRD EXPERIMENT ##########################

##open the smurfs file 
smurfs_ao_cg4360 <- read_excel("../data/longevity/smurf_count_data/AdultOnly_cg4360_FEMALES_SMURFS.xlsm", range = "A9:CM33", col_names = FALSE) %>% 
  select(-c( "...4", "...5", "...6")) 
colnames(smurfs_ao_cg4360) <- c("gene", "concentration", "vial", seq(0, ncol(smurfs_ao_cg4360)-3, by = 1)) 
smurfs_ao_cg4360 <- smurfs_ao_cg4360 %>% select_if(~!any(is.na(.))) %>% filter(gene == "51813 CG4360") %>% filter(vial != 30)
#open the file
AO_cg4360_val <- read_excel("../data/longevity/AdultOnly_cg4360_FEMALES_val.xlsm", range = "A9:DB32", col_names = FALSE) %>% 
  select(-c("...3", "...4", "...5")) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) 
colnames(AO_cg4360_val) <- c("gene", "concentration", "initial_number", seq(1, ncol(AO_cg4360_val)-3, by = 1))

#AO_val4 <- AO_val4 %>% filter(gene %in% smurfs_ao_val4$gene) #filter genes that are not needed 

days_sampling_ao <- colnames(smurfs_ao_cg4360)[4:length(colnames(smurfs_ao_cg4360))]
#"21" "37" "49" "56" "63" "70" "77" "84"

#computing the remaining number of flies at each smurf assay point 

nflies_21 = AO_cg4360_val %>% select(matches("init"):matches("20")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies21 = initial_number -(rowsum - initial_number)) %>% select(nflies21)
nflies_37 = AO_cg4360_val %>% select(matches("init"):matches("36")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies37 = initial_number -(rowsum - initial_number)) %>% select(nflies37)
nflies_49 = AO_cg4360_val %>% select(matches("init"):matches("48")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies49 = initial_number -(rowsum - initial_number)) %>% select(nflies49)
nflies_56 = AO_cg4360_val %>% select(matches("init"):matches("55")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies56 = initial_number -(rowsum - initial_number)) %>% select(nflies56)
nflies_63 = AO_cg4360_val %>% select(matches("init"):matches("62")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies63 = initial_number -(rowsum - initial_number)) %>% select(nflies63)
nflies_70 = AO_cg4360_val %>% select(matches("init"):matches("69")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies70 = initial_number -(rowsum - initial_number)) %>% select(nflies70)
nflies_77 = AO_cg4360_val %>% select(matches("init"):matches("76")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies77 = initial_number -(rowsum - initial_number)) %>% select(nflies77)
nflies_84 = AO_cg4360_val %>% select(matches("init"):matches("83")) %>% mutate(rowsum =rowSums(.)) %>% mutate(nflies84 = initial_number -(rowsum - initial_number)) %>% select(nflies84)

#####add those infos to the smurf file
smurfs_ao_modified = smurfs_ao_cg4360 %>% cbind(data.frame(nflies_21, nflies_37, nflies_49, nflies_56, nflies_63, nflies_70, nflies_77, nflies_84))
smurf_ao_prop = smurfs_ao_modified %>%
  mutate(p21 = ifelse(nflies21 > 0, smurfs_ao_modified$`21`/nflies21, NA)) %>% 
  mutate(p37 = ifelse(nflies37 > 0, smurfs_ao_modified$`37`/nflies37, NA)) %>% 
  mutate(p49 = ifelse(nflies49 > 0, smurfs_ao_modified$`49`/nflies49, NA)) %>% 
  mutate(p56 = ifelse(nflies56 > 0, smurfs_ao_modified$`56`/nflies56, NA)) %>% 
  mutate(p63 = ifelse(nflies63 > 0, smurfs_ao_modified$`63`/nflies63, NA)) %>% 
  mutate(p70 = ifelse(nflies70 > 0, smurfs_ao_modified$`70`/nflies70, NA )) %>% 
  mutate(p77 = ifelse(nflies77 > 0, smurfs_ao_modified$`77`/nflies77, NA)) %>% 
  mutate(p84 = ifelse(nflies84 > 0, smurfs_ao_modified$`84`/nflies84, NA)) 


smurf_ao_for_plot = smurf_ao_prop %>%  select(gene, concentration, vial, p21, p37, p49, p56, p63, p70, p77, p84)
colnames(smurf_ao_for_plot) = c("gene", "concentration", "vial", days_sampling_ao)

smurf_ao_for_plot =  smurf_ao_for_plot %>%  gather(., key = 'day', value = "prop_smurf", -gene, -concentration, - vial )
smurf_ao_for_plot$concentration = factor(smurf_ao_for_plot$concentration, levels = c("RU0", "RU10", "RU50", "RU100", "RU200"))

ggplot(smurf_ao_for_plot %>% filter(concentration %in% c("RU0", "RU10")) %>% 
         filter(day != 84), aes(as.numeric(day), prop_smurf)) + 
  scale_fill_manual(values = c("grey", "red3")) +
  #geom_boxplot(alpha = 0.7) + #geom_jitter(aes(color = concentration), size = 2, width = 1.5) +  
  geom_point(aes(col = concentration)) + geom_smooth(aes(color = concentration, fill = concentration), method = 'lm') +
  scale_color_manual(values = c("black", "red3")) +
  xlab("days") + ylab("proportion of Smurfs")  +   ylim(0,0.5) + xlim(0,80) +theme_linedraw() +
  theme(axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14))
ggsave("../figures/supplementary/Fig_S17_smurf_flies_cg4360_AO_ru0_ru10_lm.pdf", width = 6, height = 4)

###test the time dependence for RU0
lm(prop_smurf ~ as.numeric(day) , 
   data = smurf_ao_for_plot %>% filter(concentration %in% c("RU0")) %>% 
     filter(day != 84)) %>% summary 
lm(prop_smurf ~ as.numeric(day) , 
   data = smurf_ao_for_plot %>% filter(concentration %in% c("RU10")) %>% 
     filter(day != 84)) %>% summary 
##test the time dependence for RU10 

###test the drug interaction term
lm(prop_smurf ~ as.numeric(day)+ concentration + concentration*as.numeric(day), 
   data = smurf_ao_for_plot %>% filter(concentration %in% c("RU0", "RU10")) %>% 
     filter(day != 84)) %>% summary #pval = 0.95227
