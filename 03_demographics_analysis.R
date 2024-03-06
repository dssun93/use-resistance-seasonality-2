#This script runs analyses to explore the impact of demographic and
#site of infection sampling on seasonality of resistance

#Load libraries
library(tidyverse)
library(magrittr)
library(ggplot2)
library(ggpubr)

# ##################################################
# Inputs
# ##################################################

#Load full resistance seasonal deviates
res.deviates = read_csv("tables/01_resistance_seasonal_deviates_full.csv")

#Add age groups to resistance deviates
res.deviates.edit = res.deviates %>%
  mutate(age_group = case_when(
    age <= 19 ~ "00-19",
    age > 19 & age <= 39 ~ "20-39",
    age > 39 & age <= 64 ~ "40-64",
    age > 64 ~ "65+"
  ))

# ##################################################
# Calculate actual and predicted 
# weekly mean MIC seasonal deviates
# ##################################################

predicted_actual_deviates = left_join(
  #calculate proportions of each sub-group by week
  res.deviates.edit %>%
    count(organism, t, sex, age_group, site_of_infection) %>%
    left_join(res.deviates.edit %>% count(organism, t) %>% rename(total = n), by=c("organism", "t")) %>%
    mutate(proportion = n/total),
  #calculate mean deviates by type
  res.deviates.edit %>%
    group_by(organism, drug_code, drug_class, sex, age_group, site_of_infection) %>%
    summarize(mean_dev = mean(seasonal_deviate), sem = std.error(seasonal_deviate)) %>%
    ungroup(),
  by = c("organism", "sex", "age_group", "site_of_infection")
) %>%
  #Calculate predicted MIC as weighted sum of mean MIC by type and proportion in each week
  mutate(predicted = proportion*mean_dev) %>%
  group_by(t, organism, drug_code) %>%
  summarize(y = sum(predicted)) %>%
  ungroup() %>%
  mutate(type = "Predicted mean MIC seasonal deviate") %>%
  
  #Calculate actual mean MIC per week
  bind_rows(
    res.deviates.edit %>%
      group_by(organism, drug_code, t) %>%
      summarize(y = mean(seasonal_deviate)) %>%
      ungroup() %>%
      mutate(type = "Actual mean MIC seasonal deviate")
  ) %>%
  select(organism, drug_code, week = t, type, mean_deviate = y) %>%
  arrange(organism, drug_code, week, type)

# ##################################################
# Fit actual and predicted MIC seasonal deviates 
# to a sinusoid model
# ##################################################

#Fit to sinuosoid
predicted_actual_model = predicted_actual_deviates %>%
  nest(-organism, -drug_code, -type) %>%
  mutate(omega = ifelse(organism == "E. coli" & drug_code %in% c("AMC", "AMP"), 2*pi/26, 2*pi/52)) %>%
  mutate(model = map2(data, omega, ~ nls(mean_deviate ~ amplitude * cos(.y * (week - phase)),
                                         start = list(amplitude = 0.01, phase = 0),
                                         data = .x) 
                      
  )) %>%
  mutate(model_summary = map(model, function(m){
    model_values1 = summary(m)$coefficients %>%
      as_tibble(rownames = "term") %>%
      set_colnames(c("term","estimate","std.error","statistic","p.value")) 
    
    model_values2 = confint.default(m) %>%
      set_colnames(c("ci.lower", "ci.upper")) %>%
      as_tibble(rownames = "term")
    
    return(left_join(model_values1, model_values2, by="term"))
  })) %>%
  unnest(model_summary) %>%
  select(-data, -model)

#Edit model values so that all amplitudes are positive, all phases between 0-52
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
  
  #while phase is not between 1-period, add or subtract period 
  while (!(phase_estimate >= 0 & phase_estimate <= period)) {
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

#Make spread table of just amplitude, phase, and period parameters, then edit amplitudes and phases
model_values.edit = predicted_actual_model %>%
  gather(variable, value, -(c("organism", "drug_code", "type", "term", "omega"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  mutate(period = 2*pi/omega) %>%
  mutate(estimates_edit = pmap(.l = list(a_estimate=amplitude_estimate, a_ci.lower=amplitude_ci.lower, a_ci.upper=amplitude_ci.upper,
                                         phase_estimate=phase_estimate, phase_ci.lower=phase_ci.lower, phase_ci.upper=phase_ci.upper,
                                         period=period),
                               .f = convert_a_phases_func)) %>%
  select(-amplitude_estimate, -amplitude_ci.lower, -amplitude_ci.upper, -phase_estimate, - phase_ci.lower, -phase_ci.upper) %>%
  unnest(estimates_edit) %>%
  gather(variable, value, -(c("organism", "drug_code", "type", "omega", "period"))) %>%
  separate(variable, c("term", "temp"), "_") %>%
  spread(temp, value) %>%
  #apply BH correction of p values
  group_by(term) %>%
  mutate(p.value.BH = p.adjust(p.value, method="BH")) %>%
  ungroup() %>%
  select(organism, drug_code, type, omega, term, estimate, ci.lower, ci.upper, std.error, statistic, p.value, p.value.BH)


# ##################################################
# Save outputs
# ##################################################

write_csv(predicted_actual_deviates, "tables/03_demographics_predicted_actual_weeklyMean_deviates.csv")
write_csv(model_values.edit, "tables/03_demographics_predicted_actual_model_values.csv")