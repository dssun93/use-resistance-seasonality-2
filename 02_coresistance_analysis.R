#This script calculates the phi coefficient of co-resistance between antibiotic pairs

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

# ##################################################
# Calculate phi coefficients of co-resistance
# ##################################################

#Make spread antibiogram data tables, include all data even if missing 
drugs.SA = c("CIP", "ERY", "NIT", "OXA", "PEN")
drugs.EC = c("AMC", "AMP", "CIP", "NIT")
drugs.KP = c("AMC", "CIP", "NIT")

data.Atype.SA = data.SA %>%
  select(organism, isolate_ID, drug_code, phenotype) %>%
  spread(drug_code, phenotype) %>%
  unite("A_type", all_of(drugs.SA), remove = FALSE)

data.Atype.EC = data.EC %>%
  select(organism, isolate_ID, drug_code, phenotype) %>%
  spread(drug_code, phenotype) %>%
  unite("A_type", all_of(drugs.EC), remove = FALSE)

data.Atype.KP = data.KP %>%
  select(organism, isolate_ID, drug_code, phenotype) %>%
  spread(drug_code, phenotype) %>%
  unite("A_type", all_of(drugs.KP), remove = FALSE)

#Calculate phi coefficient of resistance between abx pairs
calculate_phi_func = function(dat, drug_list, org) {
  
  lapply(drug_list, function(drA) {
    
    dat %>%
      #gather by drug code and pheno
      gather(drugB, phenoB, drug_list[drug_list != drA]) %>%
      rename(phenoA = !!sym(drA)) %>%
      #filter out pairs with NA as phenoA or phenoB
      filter(!is.na(phenoA)) %>%
      filter(!is.na(phenoB)) %>%
      #count pheno combinations
      count(drugB, phenoA, phenoB) %>%
      mutate(phenoA = ifelse(phenoA == "S", "a", "A")) %>%
      mutate(phenoB = ifelse(phenoB == "S", "b", "B")) %>%
      unite(type, phenoA, phenoB, sep="") %>%
      #spread by pheno combinaion
      spread(type,n) %>%
      mutate(drugA = drA) %>%
      #calculate phi coefficient and 95% CI
      mutate_at(c("ab", "aB", "Ab", "AB"), ~ as.double(.)) %>%
      mutate(calc_phi = pmap(
        .l = list(A = ab, B = Ab, C = aB, D = AB),
        .f = function(A, B, C, D){statpsych::ci.phi(0.05, A, B, C, D)}
      )) %>%
      mutate(phi = map_dbl(calc_phi, ~ .[[1]])) %>%
      mutate(ci.lower = map_dbl(calc_phi, ~ .[[3]])) %>%
      mutate(ci.upper = map_dbl(calc_phi, ~ .[[4]]))
    
  }) %>%
    bind_rows() %>%
    mutate(organism = org) 
}


phi.SA = calculate_phi_func(data.Atype.SA, drugs.SA, "S. aureus") 
phi.EC = calculate_phi_func(data.Atype.EC, drugs.EC, "E. coli") 
phi.KP = calculate_phi_func(data.Atype.KP, drugs.KP, "K. pneumoniae") 

#Make full odds ratio table
phi_table = bind_rows(phi.SA, phi.EC, phi.KP) %>%
  #remove duplicate pairs of drugs
  mutate_at(vars(drugA, drugB), ~ factor(., levels = c("AMC", "AMP",  "ERY", "OXA", "PEN", "CIP", "NIT"))) %>%
  mutate(A_rank = dense_rank(drugA)) %>%
  mutate(B_rank = dense_rank(drugB)) %>%
  filter(A_rank < B_rank) %>%
  select(organism, drugA, drugB, phi, ci.lower, ci.upper) %>%
  arrange(organism, drugA, drugB)


# ##################################################
# Save output
# ##################################################

write_csv(phi_table, "tables/02_phi_coefficients.csv")