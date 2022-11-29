#function to get deaths count per  day from the number of alive flies
gettingDeathsCounts = function(input) {
  output = vector("integer", length(input) -1) 
  for (pos in 1:length(input)) {
    initial= input[pos]
    alive = input[pos + 1]
    if (is.na(alive)) {   #will need to remove the last one at the end 
      alive = 0
    }
    if (initial - alive > 0) {
      output[pos] = initial - alive 
    } else {  output[pos] = 0 }
  }
  #return(output)
  return(output = c(0, output[1: length(output)-1])) 
  #[1: length(output)-1]
}

#function to get effect on longevity on the different concentrations, where df is single gene df as the ones stored in the lifespan_df
#this function is supposed to be lapplied on such list 
#first element has to be mean lifespan of RU0 and RU0 on concentration
get_longevity_effect = function(df_lifespan) {
  effect = vector("integer", nrow(df_lifespan))
  for (i in 2:length(df_lifespan$mean)){
    RU0 = df_lifespan$mean[1]
    effect[i] = (df_lifespan$mean[i] - RU0)/RU0 #to get the effect
  }
  effect = effect[-1] #remove the first value which is 0 as it corresponds to RU0
  names(effect) = df_lifespan$concentration[2:length(df_lifespan$concentration)]
  return(effect)
}

#function to get significance on longevity effect on the different concentrations, where df is single gene df as the ones stored in the adultonly_list
#the function is supposed to be applied on such a list 
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
