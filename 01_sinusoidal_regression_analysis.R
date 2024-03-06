#This script runs use and resistance sinusoidal regressions (sinusoid + linear) on:
#1) full antibiotic use data
#2) full antibiotic resistance data 
#3) CIP, NIT resistance data stratified by R/S to penicillin and macrolide antibiotics

#Load libraries
library(tidyverse)
library(magrittr)
library(lubridate)
library(plotrix) #for std.error

# ##################################################
# Inputs
# ##################################################

#Load raw use data
data.use = read_csv("raw_data/use_data.csv")

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
)

# ##################################################
# Functions
# ##################################################

#Function to create a model matrix: where rows represent data rows and columns represent
#each year or clinic/year. The model matrix consists of 0s and 1s, where a 1 means that this 
#data row is from the corresponding year or clinic/year. The model matrix is used in the 
#nls regression to allow for a separate slope and intercept term to be estimated for each
#year or clinic/year.
make_model_matrix_func = function(data, on) {
  #get number of clinical per year combinations
  n_cys = length(unique(data[[on]]))
  
  #create model matrix - same number of rows as data rows, 0/1s, 1 if this row in clinic year
  model_matrix = model.matrix(formula(paste("~ 0 +", on)), data) %>%
    set_colnames(str_remove(colnames(.), on))
  
  #check that num rows in model matrix same as data
  stopifnot(dim(model_matrix) == c(nrow(data), n_cys))
  
  return(model_matrix)
}

#Function to run the nls regression: sinusoidal(t) + linear(t)
run_regression_func = function(data, model_matrix, omega, A_init, on) {
  
  #get number of clinical per year combinations
  n_cys = length(unique(data[[on]]))
  
  #get first guess for clinic per year intercepts: use mean values
  start_intercepts = data %>%
    group_by(!!sym(on)) %>%
    summarize(y = mean(y)) %>%
    arrange(!!sym(on)) %T>%
    #check that order matches model matrix columns
    {stopifnot(all(.[[on]] == colnames(model_matrix)))} %>%
    pull(y)
  
  #get first guess for clinic per year slopes: set as 0
  start_slopes = rep(0, n_cys)
  
  #fit the model
  model = nls(y ~ amplitude * cos(omega * (t - phase)) + drop(model_matrix %*% slope) * t + drop(model_matrix %*% intercept),
              start = list(amplitude = A_init, phase = 0, slope = start_slopes, intercept = start_intercepts),
              data = data) 
  
  return(model)
}


#Function to make table of regression output params
get_model_summary_func = function(model, model_matrix, on) {
  
  names = colnames(model_matrix)
  
  model_values1 = summary(model)$coefficients %>%
    as_tibble(rownames = "term") %>%
    set_colnames(c("term","estimate","std.error","statistic","p.value")) 
  
  model_values2 = confint.default(model) %>%
    set_colnames(c("ci.lower", "ci.upper")) %>%
    as_tibble(rownames = "term")
  
  model_values = left_join(model_values1, model_values2, by="term") %>%
    #add hospital/year to intercept/slope terms
    mutate(!!on := c("","",names,names)) %>%
    mutate(term = case_when(str_detect(term, "slope") ~ "slope",
                            str_detect(term, "intercept") ~ "intercept",
                            TRUE ~ term))
  
  return(model_values)
}


#Function to make seasonal deviates table 
make_deviates_table_func = function(data, model_summary, on) {
  
  #make table of slopes and intercepts
  slopes_intercepts = model_summary %>%
    filter(term %in% c("slope", "intercept")) %>%
    select(term, estimate, !!sym(on)) %>%
    spread(term,estimate)
  
  #make table of seasonal deviates
  deviates = data %>%
    #add slopes and intercepts
    left_join(slopes_intercepts, by = c(on)) %>%
    #calculate detrended seasonal deviate from slope and intercept
    mutate(deviate = y - (slope*t + intercept)) 
  
  return(deviates)
}

#Function to convert negative amplitudes to positive and phases to range from 0 months to model period. 
#A sinusoid with amplitude -A is the same as a sinusoid with amplitude A and phase + 0.5*period.  
#Similarly, a sinusoid with phase p is the same a sinuoid with phase p-period or p+period. 
convert_a_phases_func = function(a_estimate, a_ci.lower, a_ci.upper, phase_estimate, phase_ci.lower, phase_ci.upper, period) {
  
  #if the amplitude is negative, add (-2*a) to a, a_ci.lower, a_ci.upper; add 1/2*period to phase, phase_ci.lower, phase_ci.upper
  if (a_estimate < 0) {
    a_ci.lower = a_ci.lower + -2*a_estimate
    a_ci.upper = a_ci.upper + -2*a_estimate
    a_estimate = a_estimate + -2*a_estimate
    
    phase_estimate = phase_estimate + period/2
    phase_ci.lower = phase_ci.lower + period/2
    phase_ci.upper = phase_ci.upper + period/2
  }
  
  #while phase is not between 0-period, add or subtract period 
  while (!(phase_estimate >= 0 & phase_estimate <=period)) {
    if(phase_estimate < 0) {
      phase_estimate = phase_estimate + period
      phase_ci.lower = phase_ci.lower + period
      phase_ci.upper = phase_ci.upper + period
    }
    
    if(phase_estimate > period) {
      phase_estimate = phase_estimate - period
      phase_ci.lower = phase_ci.lower - period
      phase_ci.upper = phase_ci.upper - period
    }
  }
  
  dat = data.frame(amplitude_estimate = a_estimate, amplitude_ci.lower = a_ci.lower,
                   amplitude_ci.upper = a_ci.upper, phase_estimate = phase_estimate,
                   phase_ci.lower = phase_ci.lower, phase_ci.upper = phase_ci.upper)
  return(dat)
}

# ##################################################
# Run regressions on use data
# ##################################################

#Run regressions
results.use = data.use %>%
  mutate(year = as.character(year)) %>%
  mutate(y = claims_per_1000ppl) %>%
  mutate(t = week) %>%
  nest(-drug_class) %>%
  left_join(crossing(drug_class = c("Macrolides", "Nitrofurans", "Penicillins", "Quinolones", "Tetracyclines"), period = c(26, 52)), by = c("drug_class")) %>%
  mutate(omega = 2*pi/period) %>%
  mutate(model_matrix = map(data, ~ make_model_matrix_func(., "year"))) %>%
  mutate(model = pmap(.l = list(data = data, model_matrix = model_matrix, omega = omega, A_init=0.1, on="year"), .f = run_regression_func)) %>%
  mutate(AIC = map_dbl(model, ~ AIC(.))) %>%
  mutate(model_summary = map2(model, model_matrix, ~ get_model_summary_func(.x, .y, "year"))) %>%
  mutate(deviates = pmap(.l = list(data = data, model_summary = model_summary, on = "year"), .f = make_deviates_table_func)) 

#Edit model parameters so that all amplitudes are positive
#and phases are between 0 and 12 (or 6, depending on the period) months
#Make tables of raw model parameters
use.params.raw = results.use %>%
  select(drug_class, period, omega, AIC, model_summary) %>%
  unnest(model_summary)

#Make spread table of just amplitude, phase, and period parameters, then edit amplitudes and phases
use.params.edit = use.params.raw %>%
  filter(term %in% c("amplitude", "phase")) %>%
  select(-year) %>%
  gather(variable, value, -(c("drug_class", "term", "period", "omega", "AIC"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  mutate(estimates_edit = pmap(.l = list(a_estimate=amplitude_estimate, a_ci.lower=amplitude_ci.lower, a_ci.upper=amplitude_ci.upper,
                                         phase_estimate=phase_estimate, phase_ci.lower=phase_ci.lower, phase_ci.upper=phase_ci.upper,
                                         period=period),
                               .f = convert_a_phases_func)) %>%
  select(-amplitude_estimate, -amplitude_ci.lower, -amplitude_ci.upper, -phase_estimate, - phase_ci.lower, -phase_ci.upper) %>%
  unnest(estimates_edit) %>%
  gather(variable, value, -(c("drug_class", "period", "omega", "AIC"))) %>%
  separate(variable, c("term", "temp"), "_") %>%
  spread(temp, value)

#Combine raw and edited model params tables
use.params.full = use.params.raw %>%
  filter(!(term %in% c("amplitude", "phase"))) %>%
  bind_rows(use.params.edit) %>%
  select(drug_class, period, omega, AIC, term, year, estimate, ci.lower, ci.upper, std.error, statistic, p.value) %>%
  mutate(term = factor(term, levels = c("amplitude", "phase", "slope", "intercept"))) %>%
  #apply Benjamini-Hochberg corrections to p-values
  group_by(term) %>%
  mutate(p.value.BH = p.adjust(p.value, method="BH")) %>%
  ungroup() %>%
  arrange(drug_class, period, term)

#Make mean weekly seasonal deviates table 
use.mean.deviates = results.use %>%
  #filter to just model with lower AIC for each org/drug
  group_by(drug_class) %>%
  mutate(rank = dense_rank(AIC)) %>%
  ungroup() %>%
  filter(rank == 1) %>%
  select(-AIC, -rank) %>%
  unnest(deviates) %>%
  group_by(drug_class, period, t) %>%
  summarize(mean_deviate = mean(deviate), sem = std.error(deviate)) %>%
  ungroup() %>%
  arrange(drug_class, t)

# ##################################################
# Run regressions on resistance data
# ##################################################

#Run regressions
results.res = data.res %>%
  mutate(hos_year = paste(hospital, as.character(year), sep = "_")) %>%
  mutate(y = log2(as.double(MIC))) %>%
  mutate(t = week) %>%
  nest(-organism, -drug_code, -drug_name, -drug_class) %>%
  left_join(
    crossing(drug_code = c("AMP", "CIP", "ERY", "NIT", "OXA"), period = c(26, 52)),
    by = c("drug_code")
  ) %>%
  mutate(omega = 2*pi/period) %>%
  mutate(model_matrix = map(data, ~ make_model_matrix_func(., "hos_year"))) %>%
  mutate(model = pmap(.l = list(data = data, model_matrix = model_matrix, omega = omega, A_init = 0.01, on = "hos_year"), .f = run_regression_func)) %>%
  mutate(AIC = map_dbl(model, ~ AIC(.))) %>%
  mutate(model_summary = map2(model, model_matrix, ~ get_model_summary_func(.x, .y, "hos_year"))) %>%
  mutate(deviates = pmap(.l = list(data = data, model = model_summary, on = "hos_year"), .f = make_deviates_table_func)) 

#Edit model parameters so that all sinusoid amplitudes are positive
#and phases are between 0 and 52 (or 26, depending on the period) weeks
#Make tables of raw model parameters
res.params.raw = results.res %>%
  select(organism, drug_code, drug_name, drug_class, period, omega, AIC, model_summary) %>% 
  unnest(model_summary)

#Make spread table of just amplitude, phase, and period parameters, then edit amplitudes and phases
res.params.edit = res.params.raw %>%
  filter(term %in% c("amplitude", "phase")) %>%
  select(-hos_year) %>%
  gather(variable, value, -(c("organism", "drug_code", "drug_name", "drug_class", "term", "omega", "period", "AIC"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  mutate(estimates_edit = pmap(.l = list(a_estimate=amplitude_estimate, a_ci.lower=amplitude_ci.lower, a_ci.upper=amplitude_ci.upper,
                                         phase_estimate=phase_estimate, phase_ci.lower=phase_ci.lower, phase_ci.upper=phase_ci.upper,
                                         period=period),
                               .f = convert_a_phases_func)) %>%
  select(-amplitude_estimate, -amplitude_ci.lower, -amplitude_ci.upper, -phase_estimate, - phase_ci.lower, -phase_ci.upper) %>%
  unnest(estimates_edit) %>%
  gather(variable, value, -(c("organism", "drug_code", "drug_name", "drug_class", "omega", "period", "AIC"))) %>%
  separate(variable, c("term", "temp"), "_") %>%
  spread(temp, value)

#Combine raw and edited model params tables
res.params.full = res.params.raw %>%
  filter(!(term %in% c("amplitude", "phase"))) %>%
  bind_rows(res.params.edit) %>%
  select(organism, drug_code, drug_name, drug_class, period, omega, AIC, term, hos_year, estimate, ci.lower, ci.upper, std.error, statistic, p.value) %>%
  mutate(term = factor(term, levels = c("amplitude", "phase", "slope", "intercept"))) %>%
  #apply Benjamini-Hochberg corrections
  group_by(term) %>%
  mutate(p.value.BH = p.adjust(p.value, method="BH")) %>%
  ungroup() %>%
  arrange(organism, drug_code, period, term)

#Make full seasonal deviates table
res.deviates = bind_rows(results.res) %>%
  #filter to just model with lower AIC for each org/drug
  group_by(organism, drug_code, drug_name, drug_class) %>%
  mutate(rank = dense_rank(AIC)) %>%
  ungroup() %>%
  filter(rank == 1) %>%
  select(-AIC, -rank) %>%
  unnest(deviates) %>%
  select(isolate_ID, t, hospital, age, sex, site_of_infection, organism, drug_code, drug_name, drug_class, period, seasonal_deviate = deviate) %>%
  arrange(organism, drug_code, t)

#Make mean weekly seasonal deviates table
res.mean.deviates = res.deviates %>%
  group_by(organism, drug_code, drug_name, drug_class, period, t) %>%
  summarise(mean_deviate = mean(seasonal_deviate), sem = std.error(seasonal_deviate)) %>%
  ungroup() %>%
  arrange(organism, drug_code, t)

# ##################################################
# Run regression on stratified antibiotic resistance data
# for CIP and NIT resistance stratified by R/S to
# penicillins and macrolides
# ##################################################

#Make stratified dataset
data.stratified = bind_rows(data.SA, data.EC, data.KP) %>%
  mutate(hos_year = paste(hospital, as.character(year), sep = "_")) %>%
  mutate(y = log2(as.double(MIC))) %>%
  mutate(t = week) %>%
  filter(drug_code %in% c("NIT", "CIP")) %>%
  select(isolate_ID, organism, t = week, hos_year, drugA = drug_code, y) %>%
  left_join(
    bind_rows(data.SA, data.EC, data.KP) %>%
      filter(drug_class %in% c("Penicillins", "Macrolides")) %>%
      select(organism, isolate_ID, drugB = drug_code, drugB_pheno = phenotype),
    by = c("organism", "isolate_ID")
  ) %>%
  filter(!is.na(drugB)) %>%
  filter(!is.na(drugB_pheno)) 

#Run regressions
results.stratified = data.stratified %>%
  nest(-organism, -drugA, -drugB, -drugB_pheno) %>%
  mutate(omega = 2*pi/52, period = 52) %>%
  mutate(model_matrix = map(data, ~ make_model_matrix_func(., "hos_year"))) %>%
  mutate(model = pmap(.l = list(data = data, model_matrix = model_matrix, omega = omega, A_init = 0.01, on = "hos_year"), .f = run_regression_func)) %>%
  mutate(AIC = map_dbl(model, ~ AIC(.))) %>%
  mutate(model_summary = map2(model, model_matrix, ~ get_model_summary_func(.x, .y, "hos_year"))) %>%
  mutate(deviates = pmap(.l = list(data = data, model = model_summary, on = "hos_year"), .f = make_deviates_table_func)) 

#Edit model parameters so that all sinusoid amplitudes are positive
#and phases are between 0 and 52 (or 26, depending on the period) weeks
#Make tables of raw model parameters
res.strat.params.raw = results.stratified %>%
  select(organism, drugA, drugB, drugB_pheno, period, omega, AIC, model_summary) %>% 
  unnest(model_summary)

#Make spread table of just amplitude, phase, and period parameters, then edit amplitudes and phases
res.strat.params.edit = res.strat.params.raw %>%
  filter(term %in% c("amplitude", "phase")) %>%
  select(-hos_year) %>%
  gather(variable, value, -(c("organism", "drugA", "drugB", "drugB_pheno", "term", "omega", "period", "AIC"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  mutate(estimates_edit = pmap(.l = list(a_estimate=amplitude_estimate, a_ci.lower=amplitude_ci.lower, a_ci.upper=amplitude_ci.upper,
                                         phase_estimate=phase_estimate, phase_ci.lower=phase_ci.lower, phase_ci.upper=phase_ci.upper,
                                         period=period),
                               .f = convert_a_phases_func)) %>%
  select(-amplitude_estimate, -amplitude_ci.lower, -amplitude_ci.upper, -phase_estimate, - phase_ci.lower, -phase_ci.upper) %>%
  unnest(estimates_edit) %>%
  gather(variable, value, -(c("organism", "drugA", "drugB", "drugB_pheno", "omega", "period", "AIC"))) %>%
  separate(variable, c("term", "temp"), "_") %>%
  spread(temp, value)

#Combine raw and edited model params tables
res.strat.params.full = res.strat.params.raw %>%
  filter(!(term %in% c("amplitude", "phase"))) %>%
  bind_rows(res.strat.params.edit) %>%
  select(organism, drugA, drugB, drugB_pheno, period, omega, AIC, term, hos_year, estimate, ci.lower, ci.upper, std.error, statistic, p.value) %>%
  mutate(term = factor(term, levels = c("amplitude", "phase", "slope", "intercept"))) %>%
  #apply Benjamini-Hochberg corrections
  group_by(term) %>%
  mutate(p.value.BH = p.adjust(p.value, method="BH")) %>%
  ungroup() %>%
  arrange(organism, drugA, drugB, drugB_pheno, term)

#Make mean weekly seasonal deviates table
res.strat.mean.deviates = results.stratified %>%
  unnest(deviates) %>%
  group_by(organism, drugA, drugB, drugB_pheno, period, t) %>%
  summarise(mean_deviate = mean(deviate), sem = std.error(deviate)) %>%
  ungroup() %>%
  arrange(organism, drugA, drugB, drugB_pheno, t)


# ##################################################
# Save outputs
# ##################################################

#Save regression model values
write_csv(use.params.full, "tables/01_use_model_values.csv")
write_csv(res.params.full, "tables/01_resistance_model_values.csv")
write_csv(res.strat.params.full, "tables/01_resistance_stratified_model_values.csv")

#Save resistance seasonal deviates
write_csv(res.deviates, "tables/01_resistance_seasonal_deviates_full.csv")

#Save weekly mean seasonal deviates
write_csv(use.mean.deviates, "tables/01_use_weeklyMean_seasonal_deviates.csv")
write_csv(res.mean.deviates, "tables/01_resistance_weeklyMean_seasonal_deviates.csv")
write_csv(res.strat.mean.deviates, "tables/01_resistance_stratified_weeklyMean_seasonal_deviates.csv")
