#This script makes the figures and tables included in the publication

#Figure 1: schematic
#Figure 2: use-resistance seasonality
#Figure S1: amp resistance E. coli seasonality
#Figure 3: phi coefficients
#Figure 4, S2: stratified analysis
#Figure 5, S3, S4: demographics MICs and proportions
#Figure 6: actual vs predicted demographic analysis
#Figure 7: relative importance


#Load libraries
library(tidyverse)
library(magrittr)
library(DescTools)
library(plotrix)
library(cowplot)
library(ggpubr)

# ##################################################
# Inputs
# ##################################################

#Inputs for Figure 2, S1
#Load sinusoidal model values and deviates (for fig 1)
data.f2.model.use = read_csv("tables/01_use_model_values.csv")
data.f2.model.res = read_csv("tables/01_resistance_model_values.csv")
data.f2.deviates.res = read_csv("tables/01_resistance_weeklyMean_seasonal_deviates.csv")

#Filter model parameter tables to the best-fitting model (by AIC)
#between a 26 and 52 week period for each drug class or bug/drug
data.f2.model.use.fil = data.f2.model.use %>%
  group_by(drug_class) %>%
  mutate(rank = dense_rank(AIC)) %>%
  ungroup() %>%
  filter(rank == 1) %>%
  select(-AIC, -rank)

data.f2.model.res.fil = data.f2.model.res %>%
  group_by(organism, drug_name) %>%
  mutate(rank = dense_rank(AIC)) %>%
  ungroup() %>%
  mutate(rank = ifelse(organism == "E. coli" & drug_code == "AMP" & period == 52, 1, rank)) %>%
  mutate(rank = ifelse(organism == "E. coli" & drug_code == "AMP" & period == 26, 2, rank)) %>%
  filter(rank == 1) %>%
  select(-AIC, -rank)

#Inputs for Figure 3
#Load co-resistance phi coefficient table
data.f3.phi = read_csv("tables/02_phi_coefficients.csv")

#Inputs for Figure 4, S2
#Load stratified regression model values and deviates tables
data.f4.strat.model.values = read_csv("tables/01_resistance_stratified_model_values.csv")
data.f4.strat.deviates = read_csv("tables/01_resistance_stratified_weeklyMean_seasonal_deviates.csv")

#Inputs for Figure 5, S3, S4
#Load raw resistance data (for fig 5)
data.SA = read_csv("raw_data/Saureus_resistance_data.csv")
data.EC = read_csv("raw_data/Ecoli_resistance_data.csv")
data.KP = read_csv("raw_data/Kpneumoniae_resistance_data.csv")

#Make combined data table of raw data
data.res = bind_rows(
  data.SA %>%
    filter(drug_code %in% c("CIP", "ERY", "NIT", "OXA")),
  data.EC %>%
    filter(drug_code %in% c("AMP", "CIP", "NIT")),
  data.KP %>%
    filter(drug_code %in% c("CIP", "NIT"))
) %>%
  mutate(y = log2(MIC)) %>%
  mutate(age_group = case_when(
    age <= 19 ~ "00-19",
    age > 19 & age <= 39 ~ "20-39",
    age > 39 & age <= 64 ~ "40-64",
    age > 64 ~ "65+"
  )) %>%
  mutate(site_of_infection = case_when(
    site_of_infection == "abscess_or_fluid_nos" ~ "AB",
    site_of_infection == "blood" ~ "BL",
    site_of_infection == "respiratory_tract" ~ "RT",
    site_of_infection == "skin_softtissue" ~ "SST",
    site_of_infection == "urinary_tract" ~ "UT",
  ))

#Inputs for Figure 6
#Load output of demographic analysis
data.f6.predicted.actual = read_csv("tables/03_demographics_predicted_actual_weeklyMean_deviates.csv")
data.f6.model_values = read_csv("tables/03_demographics_predicted_actual_model_values.csv")

#Inputs for Figure 7
#Load relative importance tables
data.f7.relImp.grp = read_csv("tables/04_relative_importance_grouped_v4.2.csv")
data.f7.relImp.sep = read_csv("tables/04_relative_importance_separate_v4.2.csv")


# ##################################################
# User defined plotting inputs
# ##################################################

bug_drug_labels = c("S. aureus / CIP" = expression(paste(italic("S. aureus"), " / CIP")),
                    "S. aureus / ERY" = expression(paste(italic("S. aureus"), " / ERY")),
                    "S. aureus / NIT" = expression(paste(italic("S. aureus"), " / NIT")),
                    "S. aureus / OXA" = expression(paste(italic("S. aureus"), " / OXA")),
                    "E. coli / AMP" = expression(paste(italic("E. coli"), " / AMP")),
                    "E. coli / CIP" = expression(paste(italic("E. coli"), " / CIP")),
                    "E. coli / NIT" = expression(paste(italic("E. coli"), " / NIT")),
                    "K. pneumoniae / CIP" = expression(paste(italic("K. pneumoniae"), " / CIP")),
                    "K. pneumoniae / NIT" = expression(paste(italic("K. pneumoniae"), " / NIT"))
)

bug_drug_order = c("S. aureus / ERY", "S. aureus / OXA", "E. coli / AMP",
                   "S. aureus / NIT", "E. coli / NIT", "K. pneumoniae / NIT",
                   "S. aureus / CIP", "E. coli / CIP", "K. pneumoniae / CIP")  

bug_drug_order_vert = c("S. aureus / ERY", "S. aureus / NIT", "S. aureus / CIP",
                        "S. aureus / OXA", "E. coli / NIT", "E. coli / CIP",
                        "E. coli / AMP", "K. pneumoniae / NIT", "K. pneumoniae / CIP")

#Base colors
color_sampling = "#1ca084"
color_use = "#854696"

#Cosine function
cos_func = function(t, amplitude, phase, omega, intercept) {
  amplitude * cos(omega *(t - phase)) + intercept
}


# ##################################################
# Make Figure 2
# ##################################################

#Define colors to be used across figures 2, S1
colors = setNames( c("#220050", "#0091a8", "#b30059", "#359023"), 
                   c("Macrolides", "Penicillins", "Nitrofurans",  "Quinolones") )


#Function to plot use and resistance models together
plot_use_resistance_func = function(org_drug, class, u_a, u_p, u_o, u_low, u_up, r_a, r_p, r_o, r_low, r_up, ratio) {
  col = as.character(colors[class])
  
  regressions = data.frame(t=seq(1,52,0.1)) %>%
    mutate(r_actual = map_dbl(t, ~cos_func(., r_a, r_p, r_o, 0))) %>%
    mutate(u_actual = map_dbl(t, ~cos_func(., u_a/ratio, u_p, u_o, 0))) %>%
    gather(type, value, -t) %>%
    mutate(leg = case_when(type == "r_actual" ~ "Resistance", type == "u_actual" ~ paste(class, "use"))) %>%
    mutate(leg = factor(leg, levels = c("Resistance", paste(class, "use"))))
  
  ci = data.frame(t=seq(1,52,0.1)) %>%
    mutate(r_lower = map_dbl(t, ~cos_func(., r_low, r_p, r_o, 0))) %>%
    mutate(r_upper = map_dbl(t, ~cos_func(., r_up, r_p, r_o, 0))) %>%
    mutate(u_lower = map_dbl(t, ~cos_func(., u_low, u_p, u_o, 0))) %>%
    mutate(u_upper = map_dbl(t, ~cos_func(., u_up, u_p, u_o, 0)))
  
  p = ggplot(data = regressions) +
    geom_line(aes(x = t, y = value, color = leg, linetype = leg), size = 1) +
    geom_ribbon(data = ci, aes(x = t, ymin = r_lower, ymax = r_upper), fill = col, alpha = 0.3) +
    geom_ribbon(data = ci, aes(x = t, ymin = u_upper/ratio, ymax = u_lower/ratio), fill = col, alpha = 0.3) +
    scale_color_manual(values = c(col, col)) +
    scale_y_continuous(sec.axis = sec_axis(~. * ratio), limits = c(-0.12, 0.12)) +
    ggtitle(bug_drug_labels[org_drug]) +
    xlab("Week") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
          axis.text = element_text(size = 8),
          axis.title.x = element_text(size = 9),
          axis.title.y = element_blank()
    ) 
  
  return(p)
} 

#Make plotting dataframe for figure 2
f2_data = data.f2.model.res.fil %>%
  filter(term %in% c("amplitude", "phase")) %>%
  
  select(organism, drug_code, drug_class, omega, term, estimate, ci.lower, ci.upper) %>%
  gather(variable, value, -(c("organism", "drug_code", "drug_class", "term", "omega"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  mutate(org_drug = paste(organism, "/", drug_code)) %>%
  select(org_drug, drug_class, res_amplitude = amplitude_estimate, res_phase = phase_estimate,
         res_omega = omega, res_upper = amplitude_ci.upper, res_lower = amplitude_ci.lower) %>%
  #add use model
  left_join(
    data.f2.model.use.fil %>%
      filter(term %in% c("amplitude", "phase")) %>%
      select(drug_class, omega, term, estimate, ci.lower, ci.upper) %>%
      gather(variable, value, -(c("drug_class", "term", "omega"))) %>%
      unite(temp, term, variable) %>%
      spread(temp, value) %>%
      select(drug_class, use_amplitude = amplitude_estimate, use_phase = phase_estimate, use_omega = omega,
             use_upper = amplitude_ci.upper, use_lower = amplitude_ci.lower),
    by = c("drug_class")
  ) 

#Make plots
use_res_plots = f2_data %>%
  #get use-resistance scaling ratio
  mutate(u.r_ratio = abs(use_amplitude/res_amplitude)) %>%
  group_by(drug_class) %>%
  mutate(ratio = min(u.r_ratio)) %>%
  ungroup() %>%
  #make plots
  mutate(plot = pmap(.l = list(org_drug = org_drug, class = drug_class, u_a = use_amplitude,
                               u_p = use_phase, u_o = use_omega, u_low = use_lower, u_up = use_upper,
                               r_a = res_amplitude, r_p = res_phase, r_o = res_omega, r_low = res_lower,
                               r_up = res_upper, ratio = ratio),
                     .f = plot_use_resistance_func)) %>%
  #reorder plots
  mutate(org_drug = factor(org_drug, levels = bug_drug_order_vert)) %>%
  arrange(org_drug) %>%
  pull(plot)

#Make figure
f2 = do.call(ggarrange, c(use_res_plots, nrow = 3, ncol = 3, common.legend = T, legend = "none", align = "v")) %>%
  annotate_figure(left = text_grob(expression("Seasonal deviates in resistance ("*log["2"]*"(MIC))"), size = 9, rot = 90)) %>%
  annotate_figure(right = text_grob("Seasonal deviates in use (claims/1000 people)", size = 9, rot = 270)) 

#Make figure data table
f2_fig_data = f2_data %>%
  mutate_at(vars(-org_drug, -drug_class), ~ signif(., 2))

# ##################################################
# Make Figure S1
# ##################################################

# [TO DO: Placeholder for making S1 figure (potentially add deviates points back?)]


# ##################################################
# Make Figure 3
# ##################################################

plot_fig3_func = function(dat, title) {
  
  dat %>%
    ggplot(aes(x = drugA, y = reorder(drugB, desc(drugB)), fill = phi)) +
    geom_tile(color = "black") +
    geom_text(aes(label = txt), size=3) +
    scale_fill_gradient(low = "#f2ecf4", high = "#854696", limits=c(-.1, 1)) +
    ggtitle(title) +
    labs(fill = "phi\ncoefficient") +
    theme_classic() +
    theme(plot.title = element_text(size = 11, face = "bold.italic"),
          axis.text = element_text(size = 10),
          axis.title = element_blank(),
          legend.title = element_text(size = 10),
          legend.position = "right")
}

#Transform phi data for plotting
data.f3.phi.edit = data.f3.phi %>%
  mutate_at(c("drugA", "drugB"), ~ case_when(
    . %in% c("AMC", "AMP", "OXA", "PEN") ~ paste0(., "\n(Penicillin)"),
    . == "ERY" ~ paste0(., "\n(Macrolide)"),
    . == "CIP" ~ paste0(., "\n(Quinolone)"),
    . == "NIT" ~ paste0(., "\n(Nitrofuran)")
  )) %>%
  mutate_at(vars(drugA, drugB), ~ factor(., levels = c("AMC\n(Penicillin)", "AMP\n(Penicillin)", "ERY\n(Macrolide)", "OXA\n(Penicillin)", "PEN\n(Penicillin)", "CIP\n(Quinolone)", "NIT\n(Nitrofuran)"))) %>%
  mutate_at(c("phi", "ci.lower", "ci.upper"), ~ ifelse(. < 0.1, sprintf("%.3f", signif(.,2)), sprintf("%.2f", signif(.,2)))) %>%
  mutate(txt = paste0(phi, "\n(",  ci.lower, ", ", ci.upper, ")" )) %>%
  mutate_at(c("phi", "ci.lower", "ci.upper"), ~ as.double(.))

#Make [hi plots for each organism
f3a = plot_fig3_func(data.f3.phi.edit %>% filter(organism == "E. coli"), "E. coli")
f3b = plot_fig3_func(data.f3.phi.edit %>% filter(organism == "K. pneumoniae"), "K. pneumoniae")
f3c = plot_fig3_func(data.f3.phi.edit %>% filter(organism == "S. aureus"), "S. aureus")

#Combine plots
f3 = ggarrange(
  ggarrange(f3a+theme(legend.position = "none"), f3b+theme(legend.position = "none"), nrow = 1, ncol = 2, labels = c("A.", "B."), widths = c(4,3)),
  ggarrange(f3c, labels = c("C.")),
  nrow = 2, ncol = 1, heights = c(4,  5)
)

#Make fig data table
f3_fig_data = data.f3.phi %>%
  mutate_at(c("phi", "ci.lower", "ci.upper"), ~ ifelse(. < 0.1, sprintf("%.3f", signif(.,2)), sprintf("%.2f", signif(.,2))))
            
# ##################################################
# Make Figures 4b, S2b
# ##################################################

#Sinusoid plotting function
plot_fig_4b = function(organism, drugB, R_deviates, S_deviates, R_o, R_a, R_p, R_upper, R_lower, R_sig, S_o, S_a, S_p, S_upper, S_lower, S_sig) {
  
  colors = setNames( c(color_use, "black"), c("R", "S"))
  lines = setNames(c("solid", "dotted"), c("FDR < 0.05", "n.s."))
  shapes = setNames(c(5, 0), c("R", "S"))
  
  regressions = data.frame(t=seq(1,52,0.1)) %>%
    mutate(R_actual = map_dbl(t, ~cos_func(., R_a, R_p, R_o, 0))) %>%
    mutate(S_actual = map_dbl(t, ~cos_func(., S_a, S_p, S_o, 0))) %>%
    gather(type, value, -t) %>%
    mutate(`Drug B\nphenotype` = case_when(type == "R_actual" ~ "R", type == "S_actual" ~ "S"))
  
  ci = data.frame(t=seq(1,52,0.1)) %>%
    mutate(R_lower = map_dbl(t, ~cos_func(., R_lower, R_p, R_o, 0))) %>%
    mutate(R_upper = map_dbl(t, ~cos_func(., R_upper, R_p, R_o, 0))) %>%
    mutate(S_lower = map_dbl(t, ~cos_func(., S_lower, S_p, S_o, 0))) %>%
    mutate(S_upper = map_dbl(t, ~cos_func(., S_upper, S_p, S_o, 0)))
  
  deviates = bind_rows(
    R_deviates %>% mutate(`Drug B\nphenotype` = "R"),
    S_deviates %>% mutate(`Drug B\nphenotype` = "S")
  )
  
  p = ggplot(data = regressions) +
    geom_point(data = deviates, aes(x = t, y = mean_deviate, color = `Drug B\nphenotype`, shape = `Drug B\nphenotype`), alpha = 0.4) +
    geom_errorbar(data = deviates, aes(x = t, ymin = mean_deviate - sem, ymax = mean_deviate + sem, color = `Drug B\nphenotype`), width = 0.5, alpha = 0.4) +
    geom_line(aes(x = t, y = value, color = `Drug B\nphenotype`), size = 1) +
    geom_ribbon(data = ci, aes(x = t, ymin = R_lower, ymax = R_upper), fill = color_use, alpha = 0.3) +
    geom_ribbon(data = ci, aes(x = t, ymin = S_lower, ymax = S_upper), fill = "black", alpha = 0.3) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    scale_linetype_manual(values = lines) +
    ggtitle(expr(paste(italic(!!organism), "/Drug B:", " ", !!drugB, !!R_sig, !!S_sig))) +
    xlab("Week") +
    labs(color = "Drug B\nphenotype") +
    theme_minimal() +
    theme(legend.position = "right",
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 9),
          plot.title = element_text(size = 9, hjust = 0.5),
          axis.text = element_text(size = 9),
          axis.title.x = element_text(size = 9),
          axis.title.y = element_blank()
    ) 
  
  return(p)
}

#Make table for plotting
model.params.spread = data.f4.strat.model.values %>%
  filter(term %in% c("amplitude", "phase")) %>%
  select(organism, drugA, drugB, drugB_pheno, omega, term, estimate, ci.lower, ci.upper, p.value.BH) %>%
  gather(variable, value, -(c("organism", "drugA", "drugB", "drugB_pheno", "term", "omega"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  mutate(amplitude_sig = ifelse(amplitude_p.value.BH < 0.05, T, F))

f4_plots = left_join(
  model.params.spread %>%
    filter(drugB_pheno == "NS") %>%
    select(organism, drugA, drugB, R_o = omega, R_a = amplitude_estimate, R_p = phase_estimate,
           R_upper = amplitude_ci.upper, R_lower = amplitude_ci.lower, R_sig = amplitude_sig),
  model.params.spread %>%
    filter(drugB_pheno == "S") %>%
    select(organism, drugA, drugB, S_o = omega, S_a = amplitude_estimate, S_p = phase_estimate,
           S_upper = amplitude_ci.upper, S_lower = amplitude_ci.lower, S_sig = amplitude_sig),
  by = c("organism", "drugA", "drugB")
) %>%
  mutate(R_sig = ifelse(R_sig == T, "*", "")) %>%
  mutate(S_sig = ifelse(S_sig == T, "+", "")) %>%
  #add seasonal deviates
  left_join(
    data.f4.strat.deviates %>%
      # filter(t %in% seq(1, 52, 4)) %>%
      nest(-organism, -drugA, -drugB, -drugB_pheno) %>%
      spread(drugB_pheno, data) %>%
      rename(R_deviates = NS, S_deviates = S),
    by = c("organism", "drugA", "drugB")
  ) %>%
  #Make plots
  mutate(plot = pmap(.l = list(organism = organism, drugB = drugB, R_deviates = R_deviates, S_deviates = S_deviates,
                               R_o = R_o, R_a = R_a, R_p = R_p, R_upper = R_upper, R_lower = R_lower, R_sig = R_sig,
                               S_o = S_o, S_a = S_a, S_p = S_p, S_upper = S_upper, S_lower = S_lower, S_sig = S_sig),
                     .f = plot_fig_4b)) 

#Make figures
f4B_plots = f4_plots %>%
  filter(drugA == "CIP") %>%
  arrange(organism, drugB) %>%
  pull(plot)

fS3B_plots = f4_plots %>%
  filter(drugA == "NIT") %>%
  arrange(organism, drugB) %>%
  pull(plot)

f4B = do.call(ggarrange, c(f4B_plots, nrow = 2, ncol = 3, common.legend = T, legend = "bottom")) %>%
  annotate_figure(left = text_grob(expression("Seasonal deviates in CIP resistance ("*log[2]~"(MIC))"), size = 9, rot = 90))

fS1B = do.call(ggarrange, c(fS3B_plots, nrow = 2, ncol = 3, common.legend = T, legend = "bottom")) %>%
  annotate_figure(left = text_grob(expression("Seasonal deviates in NIT resistance ("*log[2]~"(MIC))"), size = 9, rot = 90))


#Make figure data tables
f4_fig_data_model = model.params.spread %>%
  filter(drugA == "CIP") %>%
  select(-amplitude_sig) %>%
  mutate_at(vars(-organism, -drugA, -drugB, -drugB_pheno), ~ signif(., 2))

f4_fig_data_dev = data.f4.strat.deviates %>%
  filter(drugA == "CIP") %>%
  mutate_at(vars(mean_deviate, sem), ~ signif(., 2))

fS2_fig_data_model = model.params.spread %>%
  filter(drugA == "NIT") %>%
  select(-amplitude_sig) %>%
  mutate_at(vars(-organism, -drugA, -drugB, -drugB_pheno), ~ signif(., 2))

fS2_fig_data_dev = data.f4.strat.deviates %>%
  filter(drugA == "NIT") %>%
  mutate_at(vars(mean_deviate, sem), ~ signif(., 2))


# ##################################################
# Make Figures 5, S3, S4
# ##################################################

#Perform Kruskal-Wallis/Mann Whitney tests to compare MICs between groups
data.f5.pvalues = data.res %>%
  nest(-organism, -drug_code) %>%
  mutate(KW_age = map(data, ~kruskal.test(y ~ age_group, data = .))) %>%
  mutate(KW_site = map(data, ~kruskal.test(y ~ site_of_infection, data = .))) %>%
  mutate(Wilcox_sex = map(data, ~wilcox.test(y ~ sex, data = .))) %>%
  mutate(`Age group` = map_dbl(KW_age, ~ .$p.value)) %>%
  mutate(`Site of infection` = map_dbl(KW_site, ~ .$p.value)) %>%
  mutate(Sex = map_dbl(Wilcox_sex, ~.$p.value)) %>%
  select(organism, drug_code, `Age group`, Sex, `Site of infection`) %>%
  gather(demographic_type, p.value, c("Age group", "Sex", "Site of infection")) %>%
  mutate(p.value.BH = p.adjust(p.value, method="BH")) %>%
  mutate(significance = case_when(p.value.BH < 0.05 & p.value.BH >= 0.01 ~ "*",
                                  p.value.BH < 0.01 & p.value.BH >= 0.001 ~ "**",
                                  p.value.BH < 0.001 ~ "***",
                                  TRUE ~ "n.s."
  ))

#Calculate average MICs by demographic group, add p values
data.f5.meanMIC = lapply(c("sex", "age_group", "site_of_infection"), function(x){
  data.res %>%
    group_by(organism, drug_code, !!sym(x)) %>%
    summarise(mean_MIC = mean(y), sem = std.error(y)) %>%
    ungroup() %>%
    mutate(demographic_type = StrCap(str_replace_all(x, "_", " "))) %>%
    rename(demographic_group = !!sym(x))
}) %>%
  bind_rows() %>%
  left_join(data.f5.pvalues, by = c("organism", "drug_code", "demographic_type")) %>%
  mutate(drug_code = case_when(
    drug_code %in% c("AMC", "AMP", "OXA", "PEN") ~ paste0(drug_code, "\n(Penicillin)"),
    drug_code == "ERY" ~ paste0(drug_code, "\n(Macrolide)"),
    drug_code == "CIP" ~ paste0(drug_code, "\n(Quinolone)"),
    drug_code == "NIT" ~ paste0(drug_code, "\n(Nitrofuran)")
  )) %>%
  mutate(drug_code = factor(drug_code, levels = c("AMC\n(Penicillin)", "AMP\n(Penicillin)", "ERY\n(Macrolide)", "OXA\n(Penicillin)", "PEN\n(Penicillin)", "CIP\n(Quinolone)", "NIT\n(Nitrofuran)"))) %>%
  select(organism, drug_code, demographic_type, p.value, p.value.BH, significance, demographic_group, mean_MIC, sem) %>%
  arrange(organism, drug_code, demographic_type, demographic_group)

#Calculate monthly proportion of isolates from each demographic group
data.f5.proportions = lapply(c("sex", "age_group", "site_of_infection"), function(x){
  
  data.res %>%
    count(isolate_ID, organism, month, !!sym(x)) %>%
    select(-n) %>%
    count(organism, month, !!sym(x)) %>%
    mutate(demographic_type = StrCap(str_replace_all(x, "_", " "))) %>%
    rename(demographic_group = !!sym(x))
  
}) %>%
  bind_rows() %>%
  left_join(
    data.res %>%
      count(isolate_ID, organism, month) %>%
      select(-n) %>%
      count(organism, month) %>%
      rename(total = n),
    by = c("month", "organism")
  ) %>%
  mutate(proportion = n/total) %>%
  select(organism, demographic_type, demographic_group, month, proportion) %>%
  arrange(organism, demographic_type, demographic_group, month)


#Make figures
demographic_colors = setNames(c("#A4D9CD", "#76C6B5", "#49B39C", "#1ca084", "#168069", "#76C6B5", "#1ca084", "#A4D9CD", "#76C6B5", "#49B39C", "#1ca084"), 
                              c("SST", "AB", "BL", "RT", "UT", "F", "M", "00-19", "20-39", "40-64", "65+"))

plot_fig5_func = function(dat1, dat2, y_range, height_ratio) {
  
  #Make blank data (used to manually adjust y-axis limits)
  blank_data = dat1 %>%
    group_by(organism, drug_code, demographic_type) %>%
    summarise(max_val = max(mean_MIC)) %>%
    ungroup() %>%
    mutate(x = 0) %>%
    mutate(y = max_val + 0.05) %>%
    select(drug_code, demographic_type, x, y)
  
  #1A. Make mean MIC plot
  p1_dat = dat1 %>%
    mutate(demographic_group = factor(demographic_group, c("00-19", "20-39", "F", "M", "40-64", "65+", "SST", "AB", "BL", "RT", "UT"))) %>%
    mutate(demographic_type = factor(demographic_type, c("Age group", "Sex", "Site of infection")))
  
  p1 = ggplot(data = p1_dat) +
    facet_grid(drug_code ~ demographic_type, scales = "free", space = "free_x") +
    geom_point(aes(x = demographic_group, y = mean_MIC, color = demographic_group), size = 2) + 
    geom_errorbar(aes(x = demographic_group, y = mean_MIC, ymin = mean_MIC - sem, ymax = mean_MIC + sem, color = demographic_group,), width = 0.4) +
    geom_blank(data = blank_data, aes(x = x, y = y)) +
    geom_label(aes(label = significance), x = -Inf, y = Inf, hjust = "inward", vjust = "inward", color = "black") +
    scale_color_manual(values = demographic_colors) + 
    ylab(expr("Mean  "~log["2"]*"(MIC)")) +
    theme_light() + 
    theme(strip.text = element_text(size = 10, color = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 10),
          axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 9),
          strip.background = element_blank(),
          legend.position = "none"
    )
  
  #1B. Make monthly proportions plot
  p2_dat = dat2 %>%
    mutate(demographic_type = factor(demographic_type, c("Age group", "Sex", "Site of infection"))) %>%
    mutate(lab = paste0("  ", str_pad(demographic_group, 5, side=c("right"), pad=" "))) %>%
    mutate(demographic_group = factor(demographic_group, c("00-19", "20-39", "40-64", "65+", "F", "M", "SST", "AB", "BL", "RT", "UT"))) 
  
  p2 = ggplot(data = p2_dat, aes(x = month, y = proportion, color = demographic_group, group = demographic_group)) +
    facet_wrap(~ demographic_type) +
    geom_point(size = 1.5) +
    geom_line() +
    directlabels::geom_dl(aes(label = lab), method = "last.qp") + 
    scale_color_manual(values = demographic_colors) + 
    ylab("Proportion") +
    xlab("Month") +
    ylim(y_range) +
    xlim(c(1, 15)) +
    theme_light() + 
    theme(strip.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 9),
          legend.position = "none",
          strip.background = element_blank()
    )
  
  ggarrange(p1, p2 + theme(legend.position = "none"), nrow = 2, ncol = 1, heights = height_ratio, labels = c("A.", "B."))
  
  
}

f5 = plot_fig5_func(data.f5.meanMIC %>% filter(organism == "S. aureus"), data.f5.proportions %>% filter(organism == "S. aureus"), c(0,0.6), c(3,1.5))
fS2 = plot_fig5_func(data.f5.meanMIC %>% filter(organism == "E. coli"), data.f5.proportions %>% filter(organism == "E. coli"), c(0,1), c(2.5, 1.5))
fS3 = plot_fig5_func(data.f5.meanMIC %>% filter(organism == "K. pneumoniae"), data.f5.proportions %>% filter(organism == "K. pneumoniae"), c(0,1), c(2, 1.5))


#Make figure data tables
f5_data_meanMIC = data.f5.meanMIC %>% filter(organism == "S. aureus") %>% mutate_at(vars(p.value, p.value.BH, mean_MIC, sem), ~signif(., 2))
fS2_data_meanMIC = data.f5.meanMIC %>% filter(organism == "E. coli") %>% mutate_at(vars(p.value, p.value.BH, mean_MIC, sem), ~signif(., 2))
fS3_data_meanMIC = data.f5.meanMIC %>% filter(organism == "K. pneumoniae") %>% mutate_at(vars(p.value, p.value.BH, mean_MIC, sem), ~signif(., 2))

f5_data_proportions = data.f5.proportions %>% filter(organism == "S. aureus") %>% mutate_at(vars(proportion), ~ signif(.,2))
fS2_data_proportions = data.f5.proportions %>% filter(organism == "E. coli") %>% mutate_at(vars(proportion), ~ signif(.,2))
fS3_data_proportions = data.f5.proportions %>% filter(organism == "K. pneumoniae") %>% mutate_at(vars(proportion), ~ signif(.,2))

# ##################################################
# Make Figure 6
# ##################################################

#Sinusoid plotting function
make_sinusoid_plot_func = function(org, drug, Ac_points, Pr_points, Ac_o, Ac_a, Ac_p, Ac_upper, Ac_lower, Ac_sig, Pr_o, Pr_a, Pr_p, Pr_upper, Pr_lower, Pr_sig) {
  
  colors = setNames( c("black", color_sampling), c("Actual MIC\nseasonal\ndeviates", "Predicted MIC\nseasonal\ndeviates"))
  shapes = setNames(c(0, 5), c("Actual MIC\nseasonal\ndeviates", "Predicted MIC\nseasonal\ndeviates"))
  
  regressions = data.frame(t=seq(1,52,0.1)) %>%
    mutate(Ac_model = map_dbl(t, ~cos_func(., Ac_a, Ac_p, Ac_o, 0))) %>%
    mutate(Pr_model = map_dbl(t, ~cos_func(., Pr_a, Pr_p, Pr_o, 0))) %>%
    gather(type, value, -t) %>%
    mutate(`Data` = case_when(type == "Ac_model" ~ "Actual MIC\nseasonal\ndeviates", type == "Pr_model" ~ "Predicted MIC\nseasonal\ndeviates")) 
  
  ci = data.frame(t=seq(1,52,0.1)) %>%
    mutate(Ac_lower = map_dbl(t, ~cos_func(., Ac_lower, Ac_p, Ac_o, 0))) %>%
    mutate(Ac_upper = map_dbl(t, ~cos_func(., Ac_upper, Ac_p, Ac_o, 0))) %>%
    mutate(Pr_lower = map_dbl(t, ~cos_func(., Pr_lower, Pr_p, Pr_o, 0))) %>%
    mutate(Pr_upper = map_dbl(t, ~cos_func(., Pr_upper, Pr_p, Pr_o, 0)))
  
  points = bind_rows(
    Ac_points %>% mutate(`Data` = "Actual MIC\nseasonal\ndeviates"),
    Pr_points %>% mutate(`Data` = "Predicted MIC\nseasonal\ndeviates")
  )
  
  p = ggplot(data = regressions) +
    geom_point(data = points, aes(x = t, y = mean_deviate, color = `Data`, shape = `Data`)) +
    geom_line(aes(x = t, y = value, color = `Data`), size = 0.7) +
    geom_ribbon(data = ci, aes(x = t, ymin = Pr_lower, ymax = Pr_upper), fill = color_sampling, alpha = 0.3) +
    geom_ribbon(data = ci, aes(x = t, ymin = Ac_lower, ymax = Ac_upper), fill = "black", alpha = 0.3) +
    scale_color_manual(values = colors) +
    scale_shape_manual(values = shapes) +
    ggtitle(expr(italic(!!org)~"/"~(!!drug))) +
    guides(color = guide_legend(nrow=2,byrow=TRUE)) +
    xlab("week") +
    theme_minimal() +
    theme(legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          plot.title = element_text(size = 10, hjust = 0.5),
          axis.text = element_text(size = 9),
          axis.title.x = element_text(size = 9),
          axis.title.y = element_blank(),
          legend.position = "right"
    )
  
  return(p)
}

#Make figure 6A
model_values_spread = data.f6.model_values %>%
  select(organism, drug_code, type, omega, term, estimate, ci.lower, ci.upper, p.value.BH) %>%
  gather(variable, value, -(c("organism", "drug_code", "type", "term", "omega"))) %>%
  unite(temp, term, variable) %>%
  spread(temp, value) %>%
  mutate(amplitude_sig = ifelse(amplitude_p.value.BH < 0.05, "FDR < 0.05", "n.s."))

f6a_plots_full = left_join(
  model_values_spread %>%
    filter(type == "Actual mean MIC seasonal deviate") %>%
    select(organism, drug_code, Ac_o = omega, Ac_a = amplitude_estimate, Ac_p = phase_estimate,
           Ac_upper = amplitude_ci.upper, Ac_lower = amplitude_ci.lower, Ac_sig = amplitude_sig),
  model_values_spread %>%
    filter(type == "Predicted mean MIC seasonal deviate") %>%
    select(organism, drug_code, Pr_o = omega, Pr_a = amplitude_estimate, Pr_p = phase_estimate,
           Pr_upper = amplitude_ci.upper, Pr_lower = amplitude_ci.lower, Pr_sig = amplitude_sig),
  by = c("organism", "drug_code")
) %>%
  #add seasonal deviates
  left_join(
    data.f6.predicted.actual %>%
      rename(t = week) %>%
      nest(-organism, -drug_code, -type) %>%
      spread(type, data) %>%
      rename(Ac_points = `Actual mean MIC seasonal deviate`, Pr_points = `Predicted mean MIC seasonal deviate`),
    by = c("organism", "drug_code")
  ) %>%
  #Make plots
  mutate(plot = pmap(.l = list(org = organism, drug = drug_code, Ac_points = Ac_points, Pr_points = Pr_points,
                               Ac_o = Ac_o, Ac_a = Ac_a, Ac_p = Ac_p, Ac_upper = Ac_upper, Ac_lower = Ac_lower, Ac_sig = Ac_sig,
                               Pr_o = Pr_o, Pr_a = Pr_a, Pr_p = Pr_p, Pr_upper = Pr_upper, Pr_lower = Pr_lower, Pr_sig = Pr_sig),
                     .f = make_sinusoid_plot_func)) %>%
  mutate(plot2 = map(plot, ~ . + theme(legend.position = "none"))) 

#get f2a plots with out legend
f6a_plots = f6a_plots_full %>%
  mutate(org_drug = paste(organism, "/", drug_code)) %>%
  mutate(org_drug = factor(org_drug, levels = bug_drug_order_vert)) %>%
  arrange(org_drug) %>%
  pull(plot)

#get legend for f6a
f6a_leg = get_legend(pull(f6a_plots_full, plot)[[1]])


#Make figure 6B
f6b_toplot = data.f6.model_values %>%
  filter(term == "amplitude") %>%
  mutate(Data = case_when(type == "Actual mean MIC seasonal deviate" ~ "Actual MIC\nseasonal deviates",
                          type == "Predicted mean MIC seasonal deviate" ~ "Predicted MIC\nseasonal deviates")) %>%
  mutate(`Amplitude significance` = ifelse(p.value.BH < 0.05, "FDR < 0.05", "n.s.")) %>%
  mutate(org_drug = paste(organism, "/", drug_code)) %>%
  mutate(org_drug = factor(org_drug, levels = bug_drug_order))

f6b_plot = ggplot(aes(x = org_drug, y = estimate, color = Data), data = f6b_toplot) +
  geom_point(aes(shape = `Amplitude significance`), size = 2.5) +
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper), width = 0.3) +
  scale_x_discrete(labels = bug_drug_labels) +
  scale_color_manual(values = c("black", color_sampling)) +
  scale_shape_manual(values = c(16, 1)) +
  
  ylab(expression(atop("Amplitude of seasonality", "("*log[2]*"(MIC) deviates)"))) +
  guides(color = guide_legend(nrow=2,byrow=TRUE), shape = guide_legend(nrow=2, byrow=TRUE)) +
  theme_minimal() +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        axis.text = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 9),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.justification = "center"
  )

#Make combined figure
f6a = do.call(ggarrange, c(f6a_plots, nrow = 3, ncol = 3, common.legend = T, legend = "right")) %>%
  annotate_figure(left = text_grob(expression("Weekly mean seasonal deviate ("*log[2]*"(MIC))"), size = 9, rot = 90)) 

f6b = ggarrange(f6b_plot + theme(legend.position = "none"), get_legend(f6b_plot), NULL, nrow = 3, ncol = 1, heights = c(3,1.25,0.75))

f6 = ggarrange(f6a, f6b_plot, ncol = 1, nrow = 2, heights = c(1.75,1))

#Make figure data
f6_data_model = model_values_spread %>%
  select(-amplitude_sig) %>%
  mutate_at(vars(-organism, -drug_code, -type), ~signif(., 2))

f6_data_deviates = data.f6.predicted.actual %>%
  mutate_at(vars(mean_deviate), ~signif(., 2))

# ##################################################
# Make Figure 7
# ##################################################

relImp_color_list = setNames(c("#1ca084", "#60bca8", "#a4d9cd", "#854696"),
                             c("Age", "Sex", "Site of infection", "Pen + Mac use"))

#Function to make single relative importance plot
plot_relImp_func = function(dat1, dat2, plot_order, title) {
  
  data1 = dat1 %>%
    mutate(org_drug = factor(org_drug, levels = plot_order)) 
  
  data2 = dat2 %>%
    mutate(org_drug = factor(org_drug, levels = plot_order)) %>%
    mutate(covariate = factor(covariate, levels = c("Age", "Sex", "Site of infection", "Pen + Mac use")))
  
  relImp_summed = data2 %>%
    mutate(covariate = case_when(covariate == "Pen + Mac Use" ~ "Use",
                                 covariate %in% c("Age", "Sex", "Site of infection") ~ "Demographics")) %>%
    group_by(org_drug, R_squared, covariate) %>%
    summarize(rel_lmg = sum(rel_lmg)) %>%
    ungroup() %>%
    mutate(org_drug = factor(org_drug, levels = plot_order))
  
  p1 = ggplot(data = data1, aes(x = org_drug, y = lmg)) +
    geom_bar(stat = "identity", aes(fill = covariate), size = 0.8) + 
    scale_fill_manual(values = c("grey90", "grey50")) +
    ylim(0, 80) +
    labs(x = "", y = expression(R^2~"(%)"), fill = "Grouped covariate") +
    theme_minimal() +
    theme(title = element_text(size = 10),
          axis.title.y = element_text(size = 9),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 9),
          axis.text.x = element_blank(),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          plot.margin = unit(c(0,0,0,0), "cm"))
  
  p2 = ggplot() +
    geom_bar(data = data2, aes(x = org_drug, y = rel_lmg, fill = covariate), stat = "identity") +
    geom_bar(data = relImp_summed, aes(x = org_drug, y = rel_lmg, group = covariate), stat = "identity", color = "black", alpha = 0, size = 0.8) +
    scale_fill_manual(values = relImp_color_list) +
    scale_x_discrete(labels = bug_drug_labels) +
    labs(x = "", y = "Relative lmg (%)", fill = "Covariate") +
    theme_minimal() +
    theme(axis.title.y = element_text(size = 9),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 9, angle = 45, hjust = 1, vjust = 1),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8))
  
  plot = ggarrange(p1 + theme(legend.position = "none"), p2 + theme(legend.position = "none"), nrow = 2, ncol = 1, heights = c(1, 1.5), align = "v") %>%
    annotate_figure(top = text_grob(title, size = 11))
  
  leg1 = get_legend(p1)
  leg2 = get_legend(p2)
  
  return(list(plot = plot + theme(legend.position = "none"), leg1 = leg1, leg2 = leg2))
  
}

#Make Figure 7
#set plot order
plot_order = data.f7.relImp.grp %>%
  filter(lag == 0) %>%
  filter(covariate == "Use") %>%
  arrange(desc(R_squared)) %>%
  pull(org_drug)

#make relimp plot for each lag
f7_p0 = plot_relImp_func(data.f7.relImp.grp %>% filter(lag == 0), data.f7.relImp.sep %>% filter(lag == 0), plot_order, "0 weeks lag")
f7_p4 = plot_relImp_func(data.f7.relImp.grp %>% filter(lag == 4), data.f7.relImp.sep %>% filter(lag == 4), plot_order, "4 weeks lag")
f7_p8 = plot_relImp_func(data.f7.relImp.grp %>% filter(lag == 8), data.f7.relImp.sep %>% filter(lag == 8), plot_order, "8 weeks lag")
f7_p12 = plot_relImp_func(data.f7.relImp.grp %>% filter(lag == 12), data.f7.relImp.sep %>% filter(lag == 12), plot_order, "12 weeks lag")

#make full figure
f7_leg = ggarrange(NULL, f7_p0$leg1, f7_p0$leg2, NULL, ncol = 1, nrow = 4, heights = c(2, 1, 1, 2))
f7_plots = ggarrange(f7_p0$plot, f7_p4$plot, f7_p8$plot, f7_p12$plot, ncol = 2, nrow = 2, labels = c("A.", "B.", "C.", "D."))

f7 = ggarrange(f7_plots, f7_leg, nrow = 1, ncol = 2, widths = c(4, 1))

# ##################################################
# Make Tables S1, S2
# ##################################################

#Make use AIC table comparing 26 and 52-week period models
tableS1 = data.f2.model.use %>%
  filter(term == "amplitude") %>%
  select(drug_class, period, AIC) %>%
  group_by(drug_class) %>%
  mutate(min = min(AIC)) %>%
  mutate(diff = AIC-min) %>%
  mutate_at(c("AIC", "diff"), ~ as.character(round(., 1))) %>%
  mutate(AIC_2 = paste0(AIC, " (+", diff, ")")) %>%
  select(drug_class, period, AIC_2) %>%
  spread(period, AIC_2) %>%
  rename(`26-week period` = `26`, `52-week period` = `52`)

#Make resistance AIC table comparing 26 and 52-week period models
tableS2 = data.f2.model.res %>%
  filter(term == "amplitude") %>%
  select(organism, drug_name, period, AIC) %>%
  group_by(organism, drug_name) %>%
  mutate(min = min(AIC)) %>%
  ungroup() %>%
  mutate(diff = AIC-min) %>%
  mutate_at(c("AIC", "diff"), ~ as.character(round(., 1))) %>%
  mutate(AIC_2 = paste0(AIC, " (+", diff, ")")) %>%
  select(organism, drug_name, period, AIC_2) %>%
  spread(period, AIC_2) %>%
  rename(`26-week period` = `26`, `52-week period` = `52`)

# ##################################################
# Save outputs
# ##################################################

#Save figures 
ggsave(f2, filename = "figures/Figure2.png", height = 5, width = 6.5)
#[Placeholder for fS1]
ggsave(f3, filename = "figures/Figure3.png", height = 6.5, width = 6.5)

ggsave(f4B, filename = "figures/Figure4B.png", height = 4, width = 6.5)
ggsave(fS1B, filename = "figures/FigureS1B.png", height = 4, width = 6.5)

ggsave(f5, filename = "figures/Figure5.png", height = 8, width = 6.5)
ggsave(fS2, filename = "figures/FigureS2.png", height = 7, width = 6.5)
ggsave(fS3, filename = "figures/FigureS3.png", height = 6, width = 6.5)

ggsave(f6, filename = "figures/Figure6.png", height = 7, width = 6.5)

ggsave(f7, filename = "figures/Figure7.png", height = 7, width = 6.5)

#Save tables
write_csv(tableS1, "figures/TableS1.csv")
write_csv(tableS2, "figures/TableS12.csv")


#Save figure data
write_csv(f2_fig_data, "figure_data/Figure2_data.csv")
#[Placeholder for fS1]
write_csv(f3_fig_data, "figure_data/Figure3_data.csv")

write_csv(f4_fig_data_model, "figure_data/Figure4_model_data.csv")
write_csv(f4_fig_data_dev, "figure_data/Figure4_deviates_data.csv")
write_csv(fS2_fig_data_model, "figure_data/FigureS2_model_data.csv")
write_csv(fS2_fig_data_dev, "figure_data/FigureS2_deviates_data.csv")

write_csv(f5_data_meanMIC, "figure_data/Figure5_meanMIC_data.csv")
write_csv(f5_data_proportions, "figure_data/Figure5_proportions_data.csv")
write_csv(fS2_data_meanMIC, "figure_data/FigureS2_meanMIC_data.csv")
write_csv(fS2_data_proportions, "figure_data/FigureS2_proportions_data.csv")
write_csv(fS3_data_meanMIC, "figure_data/FigureS3_meanMIC_data.csv")
write_csv(fS3_data_proportions, "figure_data/FigureS3_proportions_data.csv")

write_csv(f6_data_model, "figure_data/Figure6_model_data.csv")
write_csv(f6_data_deviates, "figure_data/Figure6_deviates_data.csv")



