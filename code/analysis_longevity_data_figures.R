
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse); library(readxl)
library(survival); library(survminer)




####################################################################-------------------FEMALES-------------------------------------------##############################################



######TRL ADULT ONLY (SECOND EXPERIMENT -VALIDATION)

#open the file

AO_trl <- read_excel("../data/longevity/AdultOnly_Trl_FEMALES.xlsm", range = "A9:CS33", col_names = FALSE) %>% 
  select(-c("...3", "...4", "...5")) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) 
colnames(AO_trl) <- c("gene", "concentration", "initial_number", seq(1, ncol(AO_trl)-3, by = 1))
#gene = "41582 Trl KD"
#genes <- unique(AO_trl$gene)

######
#trl41582_ <- AO_trl %>% filter(gene %in% genes[3])
###without line 2, that looks like outlier for the controls
#trl41582_ = AO_val2 %>% filter(gene %in% genes[1]) %>% .[-2,]
trl41582_ = AO_trl %>% select(-gene) %>%  group_by(concentration) %>% 
  summarise(across(.cols = everything(), sum)) #%>% mutate(gene = genes[1])
trl41582_2 <-trl41582_ %>% select(-c(concentration, initial_number)) %>% t() %>% as_tibble() %>% mutate(time = rownames(.))
colnames(trl41582_2) <- c(trl41582_$concentration , "time")
trl41582_3 <- trl41582_2 %>% gather(., key = "concentration", value = "death_count", -time) %>% mutate(gene = genes[3]) %>% filter(death_count != 0) 
trl41582__n_deaths_as_time = as.integer(unlist(mapply(function(i,z) rep(i,z), trl41582_3$time, trl41582_3$death_count)))
trl41582__genes_ = unlist(mapply(function(i,z) rep(i,z), trl41582_3$gene, trl41582_3$death_count))
trl41582__concentration_ = unlist(mapply(function(i,z) rep(i,z), trl41582_3$concentration, trl41582_3$death_count))
trl41582__final = data.frame(time_death = trl41582__n_deaths_as_time, gene = trl41582__genes_, concentration = trl41582__concentration_)
fit_trl41582_ = survfit(formula = Surv(time_death) ~ concentration, data = trl41582__final , conf.type = "log")
p3 = ggsurvplot(fit_trl41582_, linetype = 1, title =  "Trl 41582 KD adult only validation")
p3$plot +  geom_hline(yintercept = 0.5, linetype = 'dashed')

###compute the mean lifespan and T50 (median lifespan) for all the 
mean_lifespan = trl41582__final %>% group_by(concentration) %>% summarise(mean = mean(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))
t50  = trl41582__final %>% group_by(concentration) %>% summarise(mean = median(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))

fit_trl_specific = survfit(formula = Surv(time_death) ~ concentration, data = trl41582__final %>% filter(concentration %in% c("RU0","RU50")), conf.type = "log")


p = ggsurvplot(fit_trl_specific, linetype = 1, palette = c("black", "red"), legend.title = "RU concentration" , 
               legend.labs = c("0 μg/mL", "50 μg/mL"), conf.int = TRUE, pval = TRUE, surv.median.line = "hv", pval.method = TRUE)

p$plot + #geom_hline(yintercept = 0.5, linetype = "dotted") +
  scale_y_continuous(limits = c(0, 1.11), expand = c(0, 0)) + scale_x_continuous(limits = c(0, 110), expand = c(0, 0))  +
  geom_point(data = t50 %>% filter(concentration %in% c("RU0", "RU50")), 
             mapping = aes(x = mean, y = 1.03), col = c("black", "red"), shape = 0, size = 3) +
  geom_point(data = mean_lifespan %>% filter(concentration %in% c("RU0", "RU50")), 
             mapping = aes(x = mean, y = 1.09), col = c("black", "red"), shape = 1, size = 3) +
  annotate("text", y = 1.09, x = 10, label = "Lifespan") +
  annotate("text", y = 1.03, x = 10, label = "T50") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
ggsave("../figures/longevity_plots/Fig7bi_TRL_AD_females_RU0RU50.pdf", width = 5.5, height = 5)




#####ADF1 AD FEMALES 

#open the file
AO_adf1 <- read_excel("../data/longevity/AO_4278ADF1_females.xlsm", range = "A9:CZ33", col_names = FALSE) %>% 
  select(-c("...3", "...4", "...5")) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) 
colnames(AO_adf1) <- c("gene", "concentration", "initial_number", seq(1, ncol(AO_adf1)-3, by = 1))

tmp2 = AO_adf1 %>% select(-gene) %>%  group_by(concentration) %>% summarise(across(.cols = everything(), sum))
tmp3 <- tmp2 %>% select(-c(concentration, initial_number)) %>% t() %>% as_tibble() %>% mutate(time = rownames(.))
colnames(tmp3) <- c(tmp2$concentration , "time")
tmp4 <- tmp3 %>% gather(., key = "concentration", value = "death_count", -time) %>% mutate(gene = "adf1_4278") %>% filter(death_count != 0) 

n_deaths_as_time = as.integer(unlist(mapply(function(i,z) rep(i,z), tmp4$time, tmp4$death_count)))
genes_ = unlist(mapply(function(i,z) rep(i,z), tmp4$gene, tmp4$death_count))
concentration_ = unlist(mapply(function(i,z) rep(i,z), tmp4$concentration, tmp4$death_count))

ao_final = data.frame(time_death =n_deaths_as_time, gene = genes_, concentration = concentration_)

fit_adf1 = survfit(formula = Surv(time_death) ~ concentration, data = ao_final , conf.type = "log")
p = ggsurvplot(fit_adf1, linetype = 1, title = "adf 4278 KD adult only females")
p$plot + geom_hline(yintercept = 0.5, linetype = "dashed")

###trl thesis plot
mean_lifespan = ao_final %>% group_by(concentration) %>% summarise(mean = mean(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))
t50  = trl41582__final %>% group_by(concentration) %>% summarise(mean = median(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))

fit_adf1 = survfit(formula = Surv(time_death) ~ concentration, data = ao_final %>% filter(concentration %in% c("RU0", "RU50")) , conf.type = "log")
p = ggsurvplot(fit_adf1 , linetype = 1, palette = c("black", "red"), legend.title = "RU concentration" , 
               legend.labs = c("0 μg/mL", "50 μg/mL"),conf.int = TRUE, pval = TRUE, surv.median.line = "hv", pval.method = TRUE)

p$plot + #geom_hline(yintercept = 0.5, linetype = "dotted") +
  scale_y_continuous(limits = c(0, 1.11), expand = c(0, 0)) + scale_x_continuous(limits = c(0, 110), expand = c(0, 0))  +
  geom_point(data = t50 %>% filter(concentration %in% c("RU0", "RU50")), 
             mapping = aes(x = mean, y = 1.03), col = c("black", "red"), shape = 0, size = 3) +
  geom_point(data = mean_lifespan %>% filter(concentration %in% c("RU0", "RU50")), 
             mapping = aes(x = mean, y = 1.09), col = c("black", "red"), shape = 1, size = 3) +
  annotate("text", y = 1.09, x = 10, label = "Lifespan") +
  annotate("text", y = 1.03, x = 10, label = "T50") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
ggsave("../figures/longevity_plots/Fig7bi_ADF1_AD_females_RU0RU50.pdf", width = 5.5, height = 5)



################CG4360 WL
#open the file
WL_cg4360 <- read_excel("../data/longevity/WholeLife_cg4360_FEMALES.xlsm", range = "A9:DF33", col_names = FALSE) %>% 
  select(-c("...3", "...4", "...5")) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) 
colnames(WL_cg4360) <- c("gene", "concentration", "initial_number", seq(1, ncol(WL_cg4360)-3, by = 1))

#cg3260
#name = c("51813 CG4360")
#tmp = WL_cg4360 %>% filter(gene == name)
# get the sum, merging all the vials
# for doing this use the function across inside summarise, as summarise_each is now deprecated 
tmp2 = WL_cg4360 %>% select(-gene) %>%  group_by(concentration) %>% summarise(across(.cols = everything(), sum))

tmp3 <- tmp2 %>% select(-c(concentration, initial_number)) %>% t() %>% as_tibble() %>% mutate(time = rownames(.))
colnames(tmp3) <- c(tmp2$concentration , "time")
tmp4 <- tmp3 %>% gather(., key = "concentration", value = "death_count", -time) %>% mutate(gene = "51813 CG4360") %>% filter(death_count != 0) 

n_deaths_as_time = as.integer(unlist(mapply(function(i,z) rep(i,z), tmp4$time, tmp4$death_count)))
genes_ = unlist(mapply(function(i,z) rep(i,z), tmp4$gene, tmp4$death_count))
concentration_ = unlist(mapply(function(i,z) rep(i,z), tmp4$concentration, tmp4$death_count))

wl_final = data.frame(time_death =n_deaths_as_time, gene = genes_, concentration = concentration_)

fit_cg4360 = survfit(formula = Surv(time_death) ~ concentration, data = wl_final , conf.type = "log")
#general plot
p = ggsurvplot(fit_cg4360, linetype = 1, title = "cg4360 KD whole females")
p$plot + geom_hline(yintercept = 0.5, linetype = "dashed")
#cgWL_LR = lapply(concentrations, logrank_auto, df =  wl_final)

##mean lifespan and t50
mean_lifespan = wl_final %>% group_by(concentration) %>% summarise(mean = mean(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))
t50  = wl_final %>% group_by(concentration) %>% summarise(mean = median(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))


###plot for the thesis
fit_cg_specific = survfit(formula = Surv(time_death) ~ concentration, 
                          data = wl_final %>% filter(concentration %in% c("RU0", "RU10")) , conf.type = "log")
p = ggsurvplot(fit_cg_specific , linetype = 1, palette = c("black", "red"), legend.title = "RU concentration" , 
               legend.labs = c("0 μg/mL", "10 μg/mL"), conf.int = TRUE, pval = TRUE, surv.median.line = "hv", pval.method = TRUE)

p$plot + #geom_hline(yintercept = 0.5, linetype = "dotted") +
  scale_y_continuous(limits = c(0, 1.11), expand = c(0, 0)) + scale_x_continuous(limits = c(0, 110), expand = c(0, 0))  +
  geom_point(data = t50 %>% filter(concentration %in% c("RU0", "RU10")), 
             mapping = aes(x = mean, y = 1.03), col = c("black", "red"), shape = 0, size = 3) +
  geom_point(data = mean_lifespan %>% filter(concentration %in% c("RU0", "RU10")), 
             mapping = aes(x = mean, y = 1.09), col = c("black", "red"), shape = 1, size = 3) +
  annotate("text", y = 1.09, x = 10, label = "Lifespan") +
  annotate("text", y = 1.03, x = 10, label = "T50") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
ggsave("../figures/longevity_plots/Fig7bi_cg4360__WL_ru0_ru10.pdf", width = 5.5, height = 5)



####################################################################################CONFIRMATION OF CG4360 LONGEVITY GENE - WL experiment showing lifespan extension - 
WL_cg4360_val <- read_excel("../data/longevity/WholeLife_cg4360_validation_females.xlsm", range = "A9:DD32", col_names = FALSE) %>% 
  select(-c("...3", "...4", "...5")) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) 
colnames(WL_cg4360_val) <- c("gene", "concentration", "initial_number", seq(1, ncol(WL_cg4360_val)-3, by = 1))

#cg3260
#name = c("51813 CG4360 KD")
#tmp = WL_cg4360_val %>% filter(gene == name)
# get the sum, merging all the vials
# for doing this use the function across inside summarise, as summarise_each is now deprecated 
tmp2 = WL_cg4360_val %>% select(-gene) %>%  group_by(concentration) %>% summarise(across(.cols = everything(), sum))

tmp3 <- tmp2 %>% select(-c(concentration, initial_number)) %>% t() %>% as_tibble() %>% mutate(time = rownames(.)) #you can ignore the warning message
colnames(tmp3) <- c(tmp2$concentration , "time")
tmp4 <- tmp3 %>% gather(., key = "concentration", value = "death_count", -time) %>% mutate(gene = "51813 CG4360 KD") %>% filter(death_count != 0) 

n_deaths_as_time = as.integer(unlist(mapply(function(i,z) rep(i,z), tmp4$time, tmp4$death_count)))
genes_ = unlist(mapply(function(i,z) rep(i,z), tmp4$gene, tmp4$death_count))
concentration_ = unlist(mapply(function(i,z) rep(i,z), tmp4$concentration, tmp4$death_count))

wl_final = data.frame(time_death =n_deaths_as_time, gene = genes_, concentration = concentration_)

fit_cg4360 = survfit(formula = Surv(time_death) ~ concentration, data = wl_final , conf.type = "log")
#general plot
p = ggsurvplot(fit_cg4360, linetype = 1, title = "cg4360 KD whole females")
p$plot + geom_hline(yintercept = 0.5, linetype = "dashed")
#cgWL_LR = lapply(concentrations, logrank_auto, df =  wl_final)

##mean lifespan and t50
mean_lifespan = wl_final %>% group_by(concentration) %>% summarise(mean = mean(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))
t50  = wl_final %>% group_by(concentration) %>% summarise(mean = median(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))


###plot for the thesis
fit_cg_specific = survfit(formula = Surv(time_death) ~ concentration, 
                          data = wl_final %>% filter(concentration %in% c("RU0", "RU10")) , conf.type = "log")
p = ggsurvplot(fit_cg_specific , linetype = 1, palette = c("black", "red"), legend.title = "RU concentration" , 
               legend.labs = c("0 μg/mL", "10 μg/mL"), conf.int = TRUE, pval = TRUE, surv.median.line = "hv", pval.method = TRUE)

p$plot + #geom_hline(yintercept = 0.5, linetype = "dotted") +
  scale_y_continuous(limits = c(0, 1.11), expand = c(0, 0)) + scale_x_continuous(limits = c(0, 110), expand = c(0, 0))  +
  geom_point(data = t50 %>% filter(concentration %in% c("RU0", "RU10")), 
             mapping = aes(x = mean, y = 1.03), col = c("black", "red"), shape = 0, size = 3) +
  geom_point(data = mean_lifespan %>% filter(concentration %in% c("RU0", "RU10")), 
             mapping = aes(x = mean, y = 1.09), col = c("black", "red"), shape = 1, size = 3) +
  annotate("text", y = 1.09, x = 10, label = "Lifespan") +
  annotate("text", y = 1.03, x = 10, label = "T50") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"),axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
ggsave("../figures/supplementary/Fig_S14_cg4360_WL_ru0_ru10_VALIDATIONSUPP.pdf", width = 5.5, height = 5)

#######################GENE NOT DISPLAYING LONGEVITY CHANGES  ---- CG4360 AD THIRD EXPERIMENT ##########################
AO_cg4360_val <- read_excel("../data/longevity/AdultOnly_cg4360_FEMALES_val.xlsm", range = "A9:DB32", col_names = FALSE) %>% 
  select(-c("...3", "...4", "...5")) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) 
colnames(AO_cg4360_val) <- c("gene", "concentration", "initial_number", seq(1, ncol(AO_cg4360_val)-3, by = 1))


# name = c("51813 CG4360")
# tmp = AO_cg4360_val %>% filter(gene == name)
# get the sum, merging all the vials
# for doing this use the function across inside summarise, as summarise_each is now deprecated 
tmp2 = AO_cg4360_val %>% select(-gene) %>%  group_by(concentration) %>% summarise(across(.cols = everything(), sum))

tmp3 <- tmp2 %>% select(-c(concentration, initial_number)) %>% t() %>% as_tibble() %>% mutate(time = rownames(.))
colnames(tmp3) <- c(tmp2$concentration , "time")
tmp4 <- tmp3 %>% gather(., key = "concentration", value = "death_count", -time) %>% mutate(gene = "51813 CG4360") %>% filter(death_count != 0) 

n_deaths_as_time = as.integer(unlist(mapply(function(i,z) rep(i,z), tmp4$time, tmp4$death_count)))
genes_ = unlist(mapply(function(i,z) rep(i,z), tmp4$gene, tmp4$death_count))
concentration_ = unlist(mapply(function(i,z) rep(i,z), tmp4$concentration, tmp4$death_count))

ao_final = data.frame(time_death =n_deaths_as_time, gene = genes_, concentration = concentration_)

fit_cg4360 = survfit(formula = Surv(time_death) ~ concentration, data = ao_final , conf.type = "log")
p = ggsurvplot(fit_cg4360, linetype = 1, title = "cg4360 KD adult only females")
p$plot + geom_hline(yintercept = 0.5, linetype = "dashed")
#ggsave("./cg4360_kd_AO.pdf", width = 10, height = 6)

##mean lifespan and t50
mean_lifespan = ao_final %>% group_by(concentration) %>% summarise(mean = mean(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))
t50  = ao_final %>% group_by(concentration) %>% summarise(mean = median(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))


###plot for the thesis
fit_cg_specific = survfit(formula = Surv(time_death) ~ concentration, 
                          data = ao_final %>% filter(concentration %in% c("RU0", "RU10")) , conf.type = "log")
p = ggsurvplot(fit_cg_specific , linetype = 1, palette = c("black", "red"), legend.title = "RU concentration" , 
               legend.labs = c("0 μg/mL", "10 μg/mL"),conf.int = TRUE, pval = TRUE, surv.median.line = "hv", pval.method = TRUE)

p$plot + #geom_hline(yintercept = 0.5, linetype = "dotted") +
  scale_y_continuous(limits = c(0, 1.11), expand = c(0, 0)) + scale_x_continuous(limits = c(0, 110), expand = c(0, 0))  +
  geom_point(data = t50 %>% filter(concentration %in% c("RU0", "RU50")), 
             mapping = aes(x = mean, y = 1.03), col = c("black", "red"), shape = 0, size = 3) +
  geom_point(data = mean_lifespan %>% filter(concentration %in% c("RU0", "RU50")), 
             mapping = aes(x = mean, y = 1.09), col = c("black", "red"), shape = 1, size = 3) +
  annotate("text", y = 1.09, x = 10, label = "Lifespan") +
  annotate("text", y = 1.03, x = 10, label = "T50") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
ggsave("../figures/supplementary/Fig_S17_cg4360_AO_ru0_ru10.pdf", width = 5.5, height = 5)







##########################AEF1, SUPPLEMENTARY FIGURE###################################"
##open file and rename colnames

AO_aef1 <- read_excel("../data/longevity/Adult_only_aef1.xlsm", range = "A9:DF33", col_names = FALSE) %>% 
  select(-c("...3", "...4", "...5")) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) 
  select(-c("...3", "...4")) 
colnames(AO_aef1) <- c("gene", "concentration", "initial_number", seq(1, ncol(AO_aef1)-2, by = 1))

#name = c("80390 Aef1 KD")
#tmp = AO_aef1 %>% filter(gene == name)
tmp2 = AO_aef1 %>% select(-gene) %>%  group_by(concentration) %>% summarise(across(.cols = everything(), sum))
tmp3 <- tmp2 %>% select(-c(concentration, initial_number)) %>% t() %>% as_tibble() %>% mutate(time = rownames(.))
colnames(tmp3) <- c(tmp2$concentration , "time")
tmp4 <- tmp3 %>% gather(., key = "concentration", value = "death_count", -time) %>% mutate(gene = name) %>% filter(death_count != 0) 

n_deaths_as_time = as.integer(unlist(mapply(function(i,z) rep(i,z), tmp4$time, tmp4$death_count)))
genes_ = unlist(mapply(function(i,z) rep(i,z), tmp4$gene, tmp4$death_count))
concentration_ = unlist(mapply(function(i,z) rep(i,z), tmp4$concentration, tmp4$death_count))

ao_final = data.frame(time_death =n_deaths_as_time, gene = genes_, concentration = concentration_)

##mean lifespan and t50
mean_lifespan = ao_final %>% group_by(concentration) %>% summarise(mean = mean(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))
t50  = ao_final %>% group_by(concentration) %>% summarise(mean = median(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))


aef1_survfit = survfit(formula = Surv(time_death) ~ concentration, data = ao_final , conf.type = "log")
p = ggsurvplot(aef1_survfit, data = ao_final, conf.int = TRUE, surv.median.line = "hv")
p$plot + 
  # theme(legend.text=element_text(size=6)) +
  # ylim(0, 1.10) +
  #geom_point(data = lifespan_df_ad1$`80390 Aef1 KD`, mapping = aes(x = mean, y = 1.05, color = concentration), shape = 4, size = 3) +
  # geom_point(data = t50_dfwl, mapping = aes(x = t50, y = 1.10, color = concentration2), shape = 6, size = 3) +
  #annotate("text", y = 1.05, x = 10, label = "Lifespan") +
  # annotate("text", y = 1.10, x = 10, label = "T50") +
  #geom_hline(yintercept = 0.5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, max(ao_final$time_death)), expand = c(0, 0))  +
  scale_color_manual(values = c("black", "brown4", "coral", "lightpink1", "firebrick1")) +
  scale_fill_manual(values = c("black", "brown4", "coral", "lightpink1", "firebrick1")) +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
ggsave("../figures/supplementary/Fig_S19_aef1_graph2.pdf", width = 5.5, height = 5)

#logrank for Aef1
##compared to control
get_effect_significance = function(df){
  df$concentration = as.character(df$concentration) %>% as.factor() #step needed as otherwise the factors are the ones before the splitting
  pval = vector("integer", length(levels(df$concentration)))
  for (i in 2:length(levels(df$concentration))) {
    df2 = filter(df, concentration %in% c("RU0", levels(df$concentration)[i]))
    test <- survdiff(Surv(time_death) ~ concentration, data = df2)
    pval[i] <- pchisq(test$chisq, df = 1, lower = FALSE)
  }
  pval = pval[-1]
  names(pval) = levels(df$concentration)[-1]
  return(pval)  
}
get_effect_significance(ao_final)
#RU10 RU50
survdiff(Surv(time_death) ~ concentration, data = ao_final %>% filter(concentration %in% c("RU10", "RU50")))
#RU50 RU100
survdiff(Surv(time_death) ~ concentration, data = ao_final %>% filter(concentration %in% c("RU50", "RU100"))) #no difference
#RU100 RU200
survdiff(Surv(time_death) ~ concentration, data = ao_final %>% filter(concentration %in% c("RU100", "RU200"))) #no difference



##################MALES PLOT

#######################################TRL ADULT ONLY
#open the file

#open the file
AO_trl_males <- read_excel("../data/longevity/AdultOnly_trl_MALES.xlsm", range = "A9:CZ58", col_names = FALSE) %>% 
  select(-c("...3", "...4", "...5")) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) 
colnames(AO_trl_males) <- c("gene", "concentration", "initial_number", seq(1, ncol(AO_trl_males)-3, by = 1))

#AO_trl_males
#genes <- unique(AO_trl_males$gene)

##Trl knockdown validation

#trl41582_ <- AO_trl_males %>% filter(gene %in% genes[2])
###without line 2, that looks like outlier for the controls
#trl41582_ = AO_val2 %>% filter(gene %in% genes[1]) %>% .[-2,]
trl41582_ = AO_trl_males %>% select(-gene) %>%  group_by(concentration) %>% 
  summarise(across(.cols = everything(), sum)) #%>% mutate(gene = genes[1])
trl41582_2 <-trl41582_ %>% select(-c(concentration, initial_number)) %>% t() %>% as_tibble() %>% mutate(time = rownames(.))
colnames(trl41582_2) <- c(trl41582_$concentration , "time")
trl41582_3 <- trl41582_2 %>% gather(., key = "concentration", value = "death_count", -time) %>% mutate(gene = genes[2]) %>% filter(death_count != 0) 
trl41582__n_deaths_as_time = as.integer(unlist(mapply(function(i,z) rep(i,z), trl41582_3$time, trl41582_3$death_count)))
trl41582__genes_ = unlist(mapply(function(i,z) rep(i,z), trl41582_3$gene, trl41582_3$death_count))
trl41582__concentration_ = unlist(mapply(function(i,z) rep(i,z), trl41582_3$concentration, trl41582_3$death_count))
trl41582__final = data.frame(time_death = trl41582__n_deaths_as_time, gene = trl41582__genes_, concentration = trl41582__concentration_)
fit_trl41582_ = survfit(formula = Surv(time_death) ~ concentration, data = trl41582__final %>% filter(concentration %in% c("RU0", "RU50")), conf.type = "log")
p3 = ggsurvplot(fit_trl41582_, linetype = 1, title =  "Trl 41582 KD adult only validation")
p3$plot +  geom_hline(yintercept = 0.5, linetype = 'dashed')


##plot for the thesis
p = ggsurvplot(fit_trl41582_, linetype = 1, legend.labs = c("0 μg/mL", "50 μg/mL"), conf.int = TRUE, pval = TRUE, surv.median.line = "hv" )

p$plot + 
  # theme(legend.text=element_text(size=6)) +
  # ylim(0, 1.10) +
  #geom_point(data = lifespan_df_ad1$`80390 Aef1 KD`, mapping = aes(x = mean, y = 1.05, color = concentration), shape = 4, size = 3) +
  # geom_point(data = t50_dfwl, mapping = aes(x = t50, y = 1.10, color = concentration2), shape = 6, size = 3) +
  #annotate("text", y = 1.05, x = 10, label = "Lifespan") +
  # annotate("text", y = 1.10, x = 10, label = "T50") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, max(trl41582__final$time_death)), expand = c(0, 0))  +
  scale_color_manual(values = c("black", "red")) +
  scale_fill_manual(values = c("black",  "red")) +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
ggsave("../figures/supplementary/Fig_S18_trl_kD_males_AO.pdf", width = 5.5, height = 5)

##info for the lifespan
trl_AO_mean_lifespan <- trl41582__final %>% group_by(concentration) %>% summarise(mean_lifespan = mean(time_death))
trl_AO_t50 <- trl41582__final %>% group_by(concentration) %>% summarise(t50 = median(time_death))


get_effect_significance(trl41582__final)


#######################################CG4360 WHOLE LIFE
#open the file

WL_M_cg4360 <- read_excel("../data/longevity/WholeLife_cg4360_MALES.xlsm", range = "A9:DB24", col_names = FALSE) %>% 
  select(-c("...3", "...4", "...5")) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) 
colnames(WL_M_cg4360) <- c("gene", "concentration", "initial_number", seq(1, ncol(WL_M_cg4360)-3, by = 1))

#genes <- unique(WL_M_cg4360$gene)

##CG4360 KD
#cg4360KD <- WL_M_cg4360 %>% filter(gene %in% genes[2])
cg4360KD = WL_M_cg4360 %>% select(-gene) %>%  group_by(concentration) %>% 
  summarise(across(.cols = everything(), sum)) #%>% mutate(gene = genes[1])
cg4360KD2 <-cg4360KD %>% select(-c(concentration, initial_number)) %>% t() %>% as_tibble() %>% mutate(time = rownames(.))
colnames(cg4360KD2) <- c(cg4360KD$concentration , "time")
cg4360KD3 <- cg4360KD2 %>% gather(., key = "concentration", value = "death_count", -time) %>% mutate(gene = genes[2]) %>% filter(death_count != 0) 
cg4360KD_n_deaths_as_time = as.integer(unlist(mapply(function(i,z) rep(i,z), cg4360KD3$time, cg4360KD3$death_count)))
cg4360KD_genes_ = unlist(mapply(function(i,z) rep(i,z), cg4360KD3$gene, cg4360KD3$death_count))
cg4360KD_concentration_ = unlist(mapply(function(i,z) rep(i,z), cg4360KD3$concentration, cg4360KD3$death_count))
cg4360KD_final = data.frame(time_death = cg4360KD_n_deaths_as_time, gene = cg4360KD_genes_, concentration = cg4360KD_concentration_)
fit_cg4360KD = survfit(formula = Surv(time_death) ~ concentration, data = cg4360KD_final %>% filter(concentration %in% c("RU0", "RU10")) , conf.type = "log")
#ggsurvplot(fit_cg4360KD, linetype = 1, title = "cg4360 KD whole life males")
#ggsave("./plots/cg4360KD_WL_val2_males.pdf", width = 10, height = 6)


##plot for thesis
p = ggsurvplot(fit_cg4360KD, linetype = 1, legend.labs = c("0 μg/mL", "10 μg/mL"), conf.int = TRUE, pval = TRUE, surv.median.line = "hv" )

p$plot + 
  # theme(legend.text=element_text(size=6)) +
  # ylim(0, 1.10) +
  #geom_point(data = lifespan_df_ad1$`80390 Aef1 KD`, mapping = aes(x = mean, y = 1.05, color = concentration), shape = 4, size = 3) +
  # geom_point(data = t50_dfwl, mapping = aes(x = t50, y = 1.10, color = concentration2), shape = 6, size = 3) +
  #annotate("text", y = 1.05, x = 10, label = "Lifespan") +
  # annotate("text", y = 1.10, x = 10, label = "T50") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, max(cg4360KD_final$time_death)), expand = c(0, 0))  +
  scale_color_manual(values = c("black", "red")) +
  scale_fill_manual(values = c("black", "red")) +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
ggsave("../figures/supplementary/Fig_S18_cg4360_kD_males_WL.pdf", width = 5.5, height = 5)

cg_wl_mean_lifespan <- cg4360KD_final %>% group_by(concentration) %>% summarise(mean_lifespan = mean(time_death))
cg_wl_t50 <- cg4360KD_final %>% group_by(concentration) %>% summarise(t50 = median(time_death))

get_effect_significance(cg4360KD_final)


#######################################adf1 WHOLE LIFE
#open the file

AO_M_adf1 <- read_excel("../data/longevity/AO_males_daGSadf1.xlsm", range = "A9:CY33", col_names = FALSE) %>% 
  select(-c("...3", "...4", "...5")) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) 
colnames(AO_M_adf1) <- c("gene", "concentration", "initial_number", seq(1, ncol(AO_M_adf1)-3, by = 1))

#genes <- unique(AO_M_adf1$gene)

##adf1 KD
#adf1KD <- AO_M_adf1 %>% filter(gene %in% genes[1])
adf1KD = AO_M_adf1 %>% select(-gene) %>%  group_by(concentration) %>% 
  summarise(across(.cols = everything(), sum)) #%>% mutate(gene = genes[1])
adf1KD2 <-adf1KD %>% select(-c(concentration, initial_number)) %>% t() %>% as_tibble() %>% mutate(time = rownames(.))
colnames(adf1KD2) <- c(adf1KD$concentration , "time")
adf1KD3 <- adf1KD2 %>% gather(., key = "concentration", value = "death_count", -time) %>% mutate(gene = "daGSadf1 males") %>% filter(death_count != 0) 
adf1KD_n_deaths_as_time = as.integer(unlist(mapply(function(i,z) rep(i,z), adf1KD3$time, adf1KD3$death_count)))
adf1KD_genes_ = unlist(mapply(function(i,z) rep(i,z), adf1KD3$gene, adf1KD3$death_count))
adf1KD_concentration_ = unlist(mapply(function(i,z) rep(i,z), adf1KD3$concentration, adf1KD3$death_count))
adf1KD_final = data.frame(time_death = adf1KD_n_deaths_as_time, gene = adf1KD_genes_, concentration = adf1KD_concentration_)
fit_adf1KD = survfit(formula = Surv(time_death) ~ concentration, data = adf1KD_final %>% filter(concentration %in% c("RU0", "RU50")) , conf.type = "log")
ggsurvplot(fit_adf1KD, linetype = 1, title = "adf1 KD whole life males")
#ggsave("./plots/adf1KD_WL_val2_males.pdf", width = 10, height = 6)

p = ggsurvplot(fit_adf1KD, linetype = 1,  legend.labs = c("0 μg/mL", "50 μg/mL"), conf.int = TRUE, pval = TRUE, surv.median.line = "hv" )

p$plot + 
  # theme(legend.text=element_text(size=6)) +
  # ylim(0, 1.10) +
  #geom_point(data = lifespan_df_ad1$`80390 Aef1 KD`, mapping = aes(x = mean, y = 1.05, color = concentration), shape = 4, size = 3) +
  # geom_point(data = t50_dfwl, mapping = aes(x = t50, y = 1.10, color = concentration2), shape = 6, size = 3) +
  #annotate("text", y = 1.05, x = 10, label = "Lifespan") +
  # annotate("text", y = 1.10, x = 10, label = "T50") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0, max(trl41582__final$time_death)), expand = c(0, 0))  +
  scale_color_manual(values = c("black", "red")) +
  scale_fill_manual(values = c("black", "red")) +
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
ggsave("../figures/supplementary/Fig_S18_adf1_kD_males_WL.pdf", width = 5.5, height = 5)

adf1_wl_mean_lifespan <- adf1KD_final %>% group_by(concentration) %>% summarise(mean_lifespan = mean(time_death))
adf1_wl_t50 <- adf1KD_final %>% group_by(concentration) %>% summarise(t50 = median(time_death))

get_effect_significance(adf1KD_final)



####control experiments, daGS sh-w, adult only and whole life. Fig S16
##code adapted started from the analysis of Sofia Sosa Marmol

AO_days <-read_excel("../data/longevity/AO_only_dags_sh.xlsx",  range = "A7:BA31", col_names = FALSE) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) 
colnames(AO_days) <- c("crosses", "concentration", "name", "id_counts", "real_counts", "2",	"4","6", "8","10","12","14","15","17","19","21","23","26","28",
                       "30","33","35","37","40","42","44","46","48","50","52","54","56","58","61","62","63","65","68","69","71","73","76",
                       "77", "79", "82", "85", "86", "89","91","93","96","97","100")
###sh white
sh_white = AO_days %>% select(-c(crosses, name, id_counts)) %>%  group_by(concentration) %>% 
  summarise(across(.cols = everything(), sum)) #%>% mutate(gene = genes[1])
sh_white2 <-sh_white %>% select(-c(concentration, real_counts)) %>% t() %>% as_tibble() %>% mutate(time = rownames(.))
colnames(sh_white2) <- c(sh_white$concentration , "time")
sh_white3 <- sh_white2 %>% gather(., key = "concentration", value = "death_count", -time) %>% mutate(gene = AO_days$crosses[1]) %>% filter(death_count != 0) 
sh_white_n_deaths_as_time = as.integer(unlist(mapply(function(i,z) rep(i,z), sh_white3$time, sh_white3$death_count)))
sh_white_genes_ = unlist(mapply(function(i,z) rep(i,z), sh_white3$gene, sh_white3$death_count))
sh_white_concentration_ = unlist(mapply(function(i,z) rep(i,z), sh_white3$concentration, sh_white3$death_count))
sh_white_final = data.frame(time_death = sh_white_n_deaths_as_time, gene = sh_white_genes_, concentration = sh_white_concentration_)
fit_sh_white = survfit(formula = Surv(time_death) ~ concentration, data = sh_white_final %>% filter(concentration %in% c("RU0", "RU200")) , conf.type = "log")
ggsurvplot(fit_sh_white, linetype = 1, title = "sh_dags")

###compute the mean lifespan and T50 (median lifespan) for all the 
mean_lifespan = sh_white_final %>% group_by(concentration) %>% summarise(mean = mean(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))
t50  = sh_white_final %>% group_by(concentration) %>% summarise(mean = median(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))

fit_shwhite_specific = survfit(formula = Surv(time_death) ~ concentration, data = sh_white_final %>% filter(concentration %in% c("RU0","RU200")), conf.type = "log")

p = ggsurvplot(fit_shwhite_specific, linetype = 1, palette = c("black", "red"), legend.title = "RU concentration" , 
               legend.labs = c("0 μg/mL", "200 μg/mL"), conf.int = TRUE, pval = TRUE, surv.median.line = "hv", pval.method = TRUE)

p$plot + #geom_hline(yintercept = 0.5, linetype = "dotted") +
  scale_y_continuous(limits = c(0, 1.11), expand = c(0, 0)) + scale_x_continuous(limits = c(0, 110), expand = c(0, 0))  +
  geom_point(data = t50 %>% filter(concentration %in% c("RU0", "RU200")), 
             mapping = aes(x = mean, y = 1.03), col = c("black", "red"), shape = 0, size = 3) +
  geom_point(data = mean_lifespan %>% filter(concentration %in% c("RU0", "RU200")), 
             mapping = aes(x = mean, y = 1.09), col = c("black", "red"), shape = 1, size = 3) +
  annotate("text", y = 1.09, x = 10, label = "Lifespan") +
  annotate("text", y = 1.03, x = 10, label = "T50") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
ggsave("../figures/supplementary/Fig_S16_daGS_shwhite_RU0RU200_AO.pdf", width = 5.5, height = 5)

survdiff(Surv(time_death) ~ concentration, data = sh_white_final %>% filter(concentration %in% c("RU0", "RU200"))) 

####whole life

WL_days <-read_excel("../data/longevity/WL_only_dags_sh.xlsx",  range = "A7:BA31", col_names = FALSE) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) 
colnames(WL_days) <- c("crosses", "concentration", "name", "id_counts", "real_counts", "1","3","5","7","9","11","12","14","16","18","20","23",
                       "25","27","30","32","34","37","39","41","43","45","47","49","51","53","55","58","59","60","62","65","66","68","70","73","74", "76", "79", "82", "83", "86",
                       "88", "90", "93", "94", "97", "100", "103", "106", "109", "112", "114")
###analysis
sh_white = WL_days %>% select(-c(crosses, name, id_counts)) %>%  group_by(concentration) %>% 
  summarise(across(.cols = everything(), sum)) #%>% mutate(gene = genes[1])
sh_white2 <-sh_white %>% select(-c(concentration, real_counts)) %>% t() %>% as_tibble() %>% mutate(time = rownames(.))
colnames(sh_white2) <- c(sh_white$concentration , "time")
sh_white3 <- sh_white2 %>% gather(., key = "concentration", value = "death_count", -time) %>% mutate(gene = AO_days$crosses[1]) %>% filter(death_count != 0) 
sh_white_n_deaths_as_time = as.integer(unlist(mapply(function(i,z) rep(i,z), sh_white3$time, sh_white3$death_count)))
sh_white_genes_ = unlist(mapply(function(i,z) rep(i,z), sh_white3$gene, sh_white3$death_count))
sh_white_concentration_ = unlist(mapply(function(i,z) rep(i,z), sh_white3$concentration, sh_white3$death_count))
sh_white_final = data.frame(time_death = sh_white_n_deaths_as_time, gene = sh_white_genes_, concentration = sh_white_concentration_)
fit_sh_white = survfit(formula = Surv(time_death) ~ concentration, data = sh_white_final %>% filter(concentration %in% c("RU0", "RU200")) , conf.type = "log")
ggsurvplot(fit_sh_white, linetype = 1, title = "sh_dags")

###compute the mean lifespan and T50 (median lifespan) for all the 
mean_lifespan = sh_white_final %>% group_by(concentration) %>% summarise(mean = mean(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))
t50  = sh_white_final %>% group_by(concentration) %>% summarise(mean = median(time_death)) %>% 
  mutate(concentration2 = paste0("concentration=", concentration))

fit_shwhite_specific = survfit(formula = Surv(time_death) ~ concentration, data = sh_white_final %>% filter(concentration %in% c("RU0","RU200")), conf.type = "log")

p = ggsurvplot(fit_shwhite_specific, linetype = 1, palette = c("black", "red"), legend.title = "RU concentration" , 
               legend.labs = c("0 μg/mL", "200 μg/mL"), conf.int = TRUE, pval = TRUE, surv.median.line = "hv", pval.method = TRUE)

p$plot + #geom_hline(yintercept = 0.5, linetype = "dotted") +
  scale_y_continuous(limits = c(0, 1.11), expand = c(0, 0)) + scale_x_continuous(limits = c(0, 110), expand = c(0, 0))  +
  geom_point(data = t50 %>% filter(concentration %in% c("RU0", "RU200")), 
             mapping = aes(x = mean, y = 1.03), col = c("black", "red"), shape = 0, size = 3) +
  geom_point(data = mean_lifespan %>% filter(concentration %in% c("RU0", "RU200")), 
             mapping = aes(x = mean, y = 1.09), col = c("black", "red"), shape = 1, size = 3) +
  annotate("text", y = 1.09, x = 10, label = "Lifespan") +
  annotate("text", y = 1.03, x = 10, label = "T50") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"), axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14))
ggsave("../figures/supplementary/Fig_S16_daGS_shwhite_RU0RU200_WL.pdf", width = 5.5, height = 5)




