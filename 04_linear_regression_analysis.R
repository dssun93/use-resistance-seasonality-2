#This script performs linear regression analyses to assess the contribution of seasonal use and demographic sampling to seasonality in resistance

#Load libraries
library(tidyverse)
library(magrittr)

# ##################################################
# Inputs
# ##################################################

#Load raw resistance data
data.SA = read_csv("raw_data/Saureus_resistance_data.csv")
data.EC = read_csv("raw_data/Ecoli_resistance_data.csv")
data.KP = read_csv("raw_data/Kpneumoniae_resistance_data.csv")

#Combine resistance data, filtered to just bug/drugs with seasonality 
data.res = bind_rows(
  data.SA %>%
    filter(drug_code %in% c("CIP", "ERY", "NIT", "OXA")),
  data.EC %>%
    filter(drug_code %in% c("AMP", "CIP", "NIT")),
  data.KP %>%
    filter(drug_code %in% c("CIP", "NIT"))
) %>%
  mutate(site_of_infection = case_when(
    site_of_infection == "abscess_or_fluid_nos" ~ "AB",
    site_of_infection == "blood" ~ "BL",
    site_of_infection == "respiratory_tract" ~ "RT",
    site_of_infection == "skin_softtissue" ~ "SST",
    site_of_infection == "urinary_tract" ~ "UT",
  ))

#Load use and MIC weekly seasonal deviates
use.deviates = read_csv("tables/01_use_weeklyMean_seasonal_deviates.csv")
res.deviates = read_csv("tables/01_resistance_weeklyMean_seasonal_deviates.csv")

# ##################################################
# Make input table for linear regressions
# ##################################################

#Calculate weekly proportion of isolates from each demographic group 
#(input as covariates for linear regression)
weekly_demographic_proportions = data.res %>%
  group_by(organism, drug_code, week) %>%
  summarize(age = mean(age)) %>%
  ungroup() %>%
  
  left_join(
    data.res %>%
      group_by(organism, drug_code, week) %>%
      count(sex) %>%
      ungroup() %>%
      spread(sex, n),
    by = c("organism", "drug_code", "week")
  ) %>%
  
  left_join(
    data.res %>%
      group_by(organism, drug_code, week) %>%
      count(site_of_infection) %>%
      ungroup() %>%
      spread(site_of_infection, n),
    by = c("organism", "drug_code", "week")
  ) %>%
  
  left_join(
    data.res %>%
      count(organism, drug_code, week) %>%
      rename(total = n),
    by = c("organism", "drug_code", "week")
    
  ) %>%
  
  mutate_at(vars(-organism, -drug_code, -week, -age, -total), ~ . / total) %>%
  select(-total)

#Make regression input table, with lags of 0-12 weeks between use and resistance 
data.toModel = left_join(
  res.deviates,
  
  bind_rows(
    use.deviates %>%
      select(drug_class, t_use = t, mean_deviate) %>%
      spread(drug_class, mean_deviate) %>%
      mutate(t = t_use) %>%
      mutate(lag = 0),
    
    use.deviates %>%
      select(drug_class, t_use = t, mean_deviate) %>%
      spread(drug_class, mean_deviate) %>%
      mutate(t = t_use + 4) %>%
      mutate(t = ifelse(t > 52, t - 52, t)) %>%
      mutate(lag = 4),
    
    use.deviates %>%
      select(drug_class, t_use = t, mean_deviate) %>%
      spread(drug_class, mean_deviate) %>%
      mutate(t = t_use + 8) %>%
      mutate(t = ifelse(t > 52, t - 52, t)) %>%
      mutate(lag = 8),
    
    use.deviates %>%
      select(drug_class, t_use = t, mean_deviate) %>%
      spread(drug_class, mean_deviate) %>%
      mutate(t = t_use + 12) %>%
      mutate(t = ifelse(t > 52, t - 52, t)) %>%
      mutate(lag = 12)
  ),
  
  by = c("t")
  
) %>%
  left_join(weekly_demographic_proportions %>% rename(t = week), by = c("organism", "drug_code", "t")) %>%
  mutate(org_drug = paste(organism, drug_code, sep = " / ")) %>%
  mutate(Pen_Mac = Penicillins + Macrolides) %>%
  select(org_drug, t, t_use, lag, y = mean_deviate, t, Macrolides, Nitrofurans, Penicillins, Quinolones, Tetracyclines, Pen_Mac, age, "F", M, AB, BL, RT, SST, UT)

# ##################################################
# Perform linear regressions part 1
# non-negative regression (TO DO)
# ##################################################


# TO ADD


# ##################################################
# Perform linear regressions part 2
# ##################################################

#Perform linear regression on following model and make VIF table
#y ~ Pen + Mac + Age + Sex + InfSite

#Run model and compute VIF
model.vifs = data.toModel %>%
  nest(-org_drug, -lag) %>%
  mutate(model = map(data, function(d){
    data = d %>% select(y, c("Penicillins", "Macrolides", "age", "M", "AB", "BL", "SST", "UT")) 
    lm(as.formula(y ~ .), data = data)
  })) %>%
  mutate(vif = map(model, ~car::vif(.))) %>%
  mutate(vif_table = map(
    vif,
    function(x) {
      x %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        set_colnames(c("covariate", "VIF"))
    }
  )) 

#Make full VIF table with models for all bug/drugs and lags
vif_table_full = model.vifs %>%
  unnest(vif_table) %>%
  select(org_drug, lag, covariate, VIF) %>%
  mutate(covariate = ifelse(covariate %in% c("Penicillins", "Macrolides"), paste(covariate, "use"), covariate)) %>%
  mutate(covariate = factor(covariate, levels = c("Penicillins use", "Macrolides use", "age", "M", "BL", "UT", "SST", "AB"))) %>%
  arrange(org_drug, lag, covariate) 

#Calculate average VIF for each covariate across all models
vif_table_avg = vif_table_full %>%
  group_by(covariate) %>%
  summarise(mean_VIF = mean(VIF)) %>%
  ungroup() %>%
  mutate(mean_VIF = round(mean_VIF, 1))

# ##################################################
# Perform linear regressions part 3
# ##################################################

#Perform linear regressions on the following model and calculate relative importance of covariates
#y ~ (Pen+Mac) + Age + Sex + InfSite

#Run linear regressions and calculate relative importance (lmg) on use and demographics covariates grouped
relImp.grouped = data.toModel %>%
  nest(-org_drug, -lag) %>%
  mutate(model = map(data, function(d){
    data = d %>% select(y, c("Pen_Mac", "age", "M", "AB", "BL", "SST", "UT")) 
    lm(as.formula(y ~ .), data = data)
  })) %>%
  #Get R_squared
  mutate(R_squared = map_dbl(model, ~ as.double(summary(.)$r.squared))) %>%
  #Calculate relative importance (lmg) of covariates
  mutate(relImp = map(
    model,
    ~ relaimpo::calc.relimp(
      .,
      type = c("lmg"),
      group = list(c("Pen_Mac"), c("age", "M", "AB", "BL", "SST", "UT")),
      groupnames = c("Use", "Demographics"),
      rela=FALSE
    )
  )) %>%
  #Make summary table
  mutate(relImp_summary = map(relImp, function(r) {
    lmg = r$lmg %>% 
      as.data.frame() %>%
      rownames_to_column("covariate") %>%
      set_colnames(c("covariate", "lmg"))
    return(lmg)
  })) 

relImp.grouped.table = relImp.grouped %>%
  select(org_drug, lag, R_squared, relImp_summary) %>%
  unnest(relImp_summary) %>%
  mutate(R_squared = R_squared*100) %>%
  mutate(lmg = lmg*100) %>%
  mutate(covariate = ifelse(covariate == "Pen_Mac", "Use", covariate))

#Run linear regressions and calculate RELATIVE lmg of each individual covariate
relImp.separate.relative = data.toModel %>%
  nest(-org_drug, -lag) %>%
  mutate(model = map(data, function(d){
    data = d %>% select(y, c("Pen_Mac", "age", "M", "AB", "BL", "SST", "UT")) 
    lm(as.formula(y ~ .), data = data)
  })) %>%
  #Get R_squared
  mutate(R_squared = map_dbl(model, ~ as.double(summary(.)$r.squared) )) %>%
  #Calculate relative importance (lmg) of covariates
  mutate(relImp = map(
    model,
    ~ relaimpo::calc.relimp(
      .,
      type = c("lmg"),
      group = c("AB", "BL", "SST", "UT"),
      groupnames = c("Site of infection"),
      rela=TRUE
    )
  )) %>%
  #Make summary table
  mutate(relImp_summary = map(relImp, function(r) {
    lmg = r$lmg %>% 
      as.data.frame() %>%
      rownames_to_column("covariate") %>%
      set_colnames(c("covariate", "rel_lmg"))
    return(lmg)
  })) 

relImp.separate.relative.table = relImp.separate.relative %>%
  select(org_drug, lag, R_squared, relImp_summary) %>%
  unnest(relImp_summary) %>%
  mutate(covariate = case_when(covariate == "M" ~ "Sex", covariate == "age" ~ "Age", TRUE ~ covariate)) %>%
  mutate(R_squared = R_squared*100) %>%
  mutate(rel_lmg = rel_lmg*100) %>%
  mutate(covariate = ifelse(covariate == "Pen_Mac", "Pen + Mac use", covariate)) 


# ##################################################
# Save outputs
# ##################################################

#Save supplemental tables
#TO DO: [placeholder of nnls table]
write_csv(vif_table_avg, "figures/Table_S3.csv")

#Save relative importance tables
write_csv(relImp.grouped.table, "tables/04_relative_importance_grouped_v4.2.csv")
write_csv(relImp.separate.relative.table, "tables/04_relative_importance_separate_v4.2.csv")
