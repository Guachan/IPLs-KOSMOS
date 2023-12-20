library(tidyverse)
library(modelr)
library(stringr)
library(ggrepel)

# ---------------------------------------------------------------
# ----------------------- Load Data ----------------------------------
# ----------------------------------------------------------------

# read lipid/nutrient/phytoplankton data
# Inputs:
#   1: path -> location of "lipid_abundances_012021.csv" file.
load_data = function(path){
  raw_lipid = read.csv(paste0(path, 'lipid_abundances_012021.csv')) %>%
    select(-X)
  # pivot on lipid class
  lipid_class = raw_lipid %>%
    select(-'Compound') %>%
    group_by(Depth, Day, Mesocosm, Class) %>%
    summarize(total_L = sum(total_L), .groups='drop_last') %>%
    pivot_wider(names_from='Class', values_from='total_L')
  # pivot on lipid compound
  lipid_compound = raw_lipid %>%
    select(-'Class') %>%
    pivot_wider(names_from='Compound', values_from='total_L')
  phyto = read.csv(paste0(path, 'phytoplankton_abundances.csv'))
  water = read.csv(paste0(path, 'water chemistry.csv'))
  # join data frames
  data = left_join(water, phyto, by=c('Depth', 'Day', 'Mesocosm')) %>%
    left_join(lipid_class, by=c('Depth', 'Day', 'Mesocosm')) %>%
    left_join(lipid_compound, by=c('Depth', 'Day', 'Mesocosm'))
  # add features
  data[is.na(data)] = 0
  data$Depth = as.factor(data$Depth)
  data$Mesocosm = as.factor(data$Mesocosm)
  # add n-to-p ratio
  data = data %>%
    mutate(NP_Ratio=(NO3+NO2+NH4)/PO4)
  return(data)
}

# read in data
path = "C:/Users/canta/MLRs and GAMs -- KOSMOS/"
data = load_data(path)

# library(openxlsx)
# path = "C:/Users/Johnny Tamanaha/Desktop/"
# write.xlsx(data, paste0(path, 'phyto_data.xlsx'), sheetName='Sheet1')

# --------------------------------------------------------------
# ------------- Important Lists ------------------------------------
# ----------------------------------------------------------------

# nutrients
nutrients = c(colnames(data)[5:13], 'NP_Ratio')
# phytoplankton
phytos = c(colnames(data)[14:18], "Pelago")
# lipid compounds
lipids_compound = colnames(data)[29:193]
# lipid classes
lipids_class = colnames(data)[21:28]


# -------------------------------------------------------------
# ------------------ MLR Analysis ------------------------------
# --------------------------------------------------------------

# create mlr models for each lipid-nutrient pair
# Inputs:
#   1: data -> output of load_data function.
#   2: lipids -> list of lipid strings (classes or compounds) to use as responses.
#   3: nutrient -> list nutrient strings to use as predictors. 
mlr_analysis_bio = function(data, lipids, phytos){
  # data frame will have one row for every lipid nutrient pair
  l = length(lipids)*length(phytos)
  # information recorded:
  #   1: lipid
  #   2: nutrient
  #   3: nutrient predictor coefficient
  #   4: nutrient predictor p-value
  #   5: r-squared of full model
  #   6: r-squared of model excluding nutrient predictor
  #   7: difference in r-squared of previous models
  #   8: number of zeros in response
  ls = ns = coef = pvalue = r_full = r_null = r_diff = zeros = rep(0, l)
  # scale nutrients, phytoplankton, and lipids (min-max)
  nutrients=c("Light", "NH4", "NO2", "NO3", "O2", "pH", "PO4", "SiO2", "Temp", "NP_Ratio")
  scale_variables = c(phytos, lipids)
  for(variable in scale_variables){
    data[variable]=(scale(data[variable]))
  }
  # fit each lipid-nutrient pair model
  for (i in 1:length(lipids)) {
    for (j in 1:length(phytos)) {
      index = (i-1)*length(phytos) + j
      # record lipid, nutrient, and number of zeros
      ls[index] = lipids[i]
      ns[index] = phytos[j]
      zeros[index] = nrow(filter(data, get(lipids[i])==0))
      # ignore models that cannot fit
      tryCatch({
        model=lm(get(lipids[i])~get(phytos[j]), data=data)
        model_null=lm(get(lipids[i])~get(phytos[j]), data=data)
        # recored r-squared, coef, and pvalue
        r_diff[index] = summary(model)$r.squared - summary(model_null)$r.squared
        r_full[index] = summary(model)$r.squared
        r_null[index] = summary(model_null)$r.squared
        coef[index] = coef(summary(model))[2, 1]
        pvalue[index] = coef(summary(model))[2, 4]
      }, error=function(e){cat('ERROR:', conditionMessage(e), '\n')})
    }
  }
  # populate output data frame
  df = data.frame('Lipid'=ls,
                  'Phyto'=ns,
                  'Coef'=coef,
                  'p_value'=pvalue,
                  'r_sq_diff'=r_diff,
                  'r_sq_full'=r_full,
                  'r_sq_null'=r_null,
                  'zeros'=zeros)
  # delete rows where model didn't fit
  df = df[df$Coef!=0,]
  # add significance flag
  df = mutate(df, sig95=ifelse(p_value<0.05, T, F))
  # return data frame in order of r-squared diff
  return(arrange(df, desc(r_sq_diff)))
}

# mlr analysis on depth 1
mlr_d <- mlr_analysis_bio(data, lipids_compound, phytos)
data_d1 = filter(data, Depth==1)
mlr_d1 = mlr_analysis_bio(data_d1, lipids_compound, phytos)


# mlr analysis on depth 2
data_d2 = filter(data, Depth==2)
mlr_d2 = mlr_analysis_bio(data_d2, lipids_compound, phytos)

#mlr analysis on phases 1and 2
data_p12 <- filter(data, Day <21)
mlr_p12 <- mlr_analysis_bio(data_p12, lipids_compound, phytos)

data_p34 <- filter(data, Day >21)
mlr_p34 <- mlr_analysis_bio(data_p34, lipids_compound, phytos)

# mlr phases and depth
data_p12_d1 <- filter(data, Day <21 & Depth ==1)
mlr_p12_d1 <- mlr_analysis_bio(data_p12_d1, lipids_compound, phytos)

data_p34_d1 <- filter(data, Day >21 & Depth ==1)
mlr_p34_d1 <- mlr_analysis_bio(data_p34_d1, lipids_compound, phytos)

# --------------------------------------------------------------------------
# -------- MLR Additional Analysis (classes, carbon atoms, double bonds) --------------
# ---------------------------------------------------------------------------

# augment the mlr_analysis output data frame for compounds
# Inputs:
#   1: mlr_data -> output of mlr_analysis function using lipid compounds.
#   2: path ->location of "lipid_abundances_012021.csv" file.
# Added Columns:
#   1: Class -> class associated to compound.
#   2: carbon_atoms
#   3: double_bonds
#   4: ca_db -> carbon_atoms:double_bonds
mlr_additional = function(mlr_data, path){
  # create a key between lipid classes and compounds
  class_compound_keys = read.csv(paste0(path, 'lipid_abundances_012021.csv')) %>%
    select(Compound, Class) %>%
    unique() 
  # add classes column
  mlr_data = left_join(mlr_data, class_compound_keys, c('Lipid'='Compound'))
  # filter out the "other" lipid class
  mlr_data = filter(mlr_data, Class!='Other')
  # extract carbon atom and double bond info from compound names
  ca_db = class_compound_keys %>%
    mutate(carbon_atoms=as.integer(str_extract(Compound, '\\d+(?=:)')),
           double_bonds=as.integer(str_extract(Compound, '(?<=:)\\d+'))) %>%
    select(-Class) %>%
    mutate(carbon_atoms = as.character(carbon_atoms)) %>%
    mutate(double_bonds = as.character(double_bonds))
  # add carbon atom and double bond info to data frame
  mlr_data = left_join(mlr_data, ca_db, by=c('Lipid'='Compound')) %>%
    distinct() %>%
    mutate(ca_db = paste0(as.character(carbon_atoms), ':', as.character(double_bonds)))
  return(mlr_data)
}

# augmented mlr analysis on depth 1

path = "C:/Users/canta/MLRs and GAMs -- KOSMOS/"

mlr_aug_d = mlr_additional(mlr_d, path) %>% 
  mutate(Class =
           case_when(grepl("MGDG", Lipid) ~ "MG",
                     grepl("PE.", Lipid) ~ "PE",
                     grepl("PG.", Lipid) ~ "PG",
                     grepl("PC.", Lipid) ~ "PC",
                     grepl("DGTS", Lipid) ~ "BL",
                     grepl("DGTA", Lipid) ~ "BL",
                     grepl("SQ", Lipid) ~ "SQ",
                     grepl("DGDG", Lipid) ~ "DG",
                     grepl("Gly.Cer", Lipid) ~ "Other",
                     grepl("PI.", Lipid) ~ "Other",
                     grepl("DGCC", Lipid) ~ "BL"))

mlr_aug_d1 = mlr_additional(mlr_d1, path) %>% 
  mutate(Class =
           case_when(grepl("MGDG", Lipid) ~ "MG",
                     grepl("PE.", Lipid) ~ "PE",
                     grepl("PG.", Lipid) ~ "PG",
                     grepl("PC.", Lipid) ~ "PC",
                     grepl("DGTS", Lipid) ~ "BL",
                     grepl("DGTA", Lipid) ~ "BL",
                     grepl("SQ", Lipid) ~ "SQ",
                     grepl("DGDG", Lipid) ~ "DG",
                     grepl("Gly.Cer", Lipid) ~ "Other",
                     grepl("PI.", Lipid) ~ "Other",
                     grepl("DGCC", Lipid) ~ "BL"))

# augmented mlr analysis on depth 2
path = "C:/Users/canta/MLRs and GAMs -- KOSMOS/"
mlr_aug_d2 = mlr_additional(mlr_d2, path) %>% 
  mutate(Class =
           case_when(grepl("MGDG", Lipid) ~ "MG",
                     grepl("PE.", Lipid) ~ "PE",
                     grepl("PG.", Lipid) ~ "PG",
                     grepl("PC.", Lipid) ~ "PC",
                     grepl("DGTS", Lipid) ~ "BL",
                     grepl("DGTA", Lipid) ~ "BL",
                     grepl("SQ", Lipid) ~ "SQ",
                     grepl("DGDG", Lipid) ~ "DG",
                     grepl("Gly.Cer", Lipid) ~ "Other",
                     grepl("PI.", Lipid) ~ "Other",
                     grepl("DGCC", Lipid) ~ "BL"))


mlr_aug_p12 = mlr_additional(mlr_p12, path) %>% 
  mutate(Class =
           case_when(grepl("MGDG", Lipid) ~ "MG",
                     grepl("PE.", Lipid) ~ "PE",
                     grepl("PG.", Lipid) ~ "PG",
                     grepl("PC.", Lipid) ~ "PC",
                     grepl("DGTS", Lipid) ~ "BL",
                     grepl("DGTA", Lipid) ~ "BL",
                     grepl("SQ", Lipid) ~ "SQ",
                     grepl("DGDG", Lipid) ~ "DG",
                     grepl("Gly.Cer", Lipid) ~ "Other",
                     grepl("PI.", Lipid) ~ "Other",
                     grepl("DGCC", Lipid) ~ "BL"))

mlr_aug_p34 = mlr_additional(mlr_p34, path) %>% 
  mutate(Class =
           case_when(grepl("MGDG", Lipid) ~ "MG",
                     grepl("PE.", Lipid) ~ "PE",
                     grepl("PG.", Lipid) ~ "PG",
                     grepl("PC.", Lipid) ~ "PC",
                     grepl("DGTS", Lipid) ~ "BL",
                     grepl("DGTA", Lipid) ~ "BL",
                     grepl("SQ", Lipid) ~ "SQ",
                     grepl("DGDG", Lipid) ~ "DG",
                     grepl("Gly.Cer", Lipid) ~ "Other",
                     grepl("PI.", Lipid) ~ "Other",
                     grepl("DGCC", Lipid) ~ "BL"))


mlr_aug_p12_d1 = mlr_additional(mlr_p12_d1, path) %>% 
  mutate(Class =
           case_when(grepl("MGDG", Lipid) ~ "MG",
                     grepl("PE.", Lipid) ~ "PE",
                     grepl("PG.", Lipid) ~ "PG",
                     grepl("PC.", Lipid) ~ "PC",
                     grepl("DGTS", Lipid) ~ "BL",
                     grepl("DGTA", Lipid) ~ "BL",
                     grepl("SQ", Lipid) ~ "SQ",
                     grepl("DGDG", Lipid) ~ "DG",
                     grepl("Gly.Cer", Lipid) ~ "Other",
                     grepl("PI.", Lipid) ~ "Other",
                     grepl("DGCC", Lipid) ~ "BL"))

mlr_aug_p34_d1 = mlr_additional(mlr_p34_d1, path) %>% 
  mutate(Class =
           case_when(grepl("MGDG", Lipid) ~ "MG",
                     grepl("PE.", Lipid) ~ "PE",
                     grepl("PG.", Lipid) ~ "PG",
                     grepl("PC.", Lipid) ~ "PC",
                     grepl("DGTS", Lipid) ~ "BL",
                     grepl("DGTA", Lipid) ~ "BL",
                     grepl("SQ", Lipid) ~ "SQ",
                     grepl("DGDG", Lipid) ~ "DG",
                     grepl("Gly.Cer", Lipid) ~ "Other",
                     grepl("PI.", Lipid) ~ "Other",
                     grepl("DGCC", Lipid) ~ "BL"))
# -----------------------------------------------------------------------
# ------------------ MLR Plotting ----------------------------------------
# ------------------------------------------------------------------------

rel_ab <- read.csv('lipid_rel_abundances.csv') %>% 
  filter(percentage >.99) %>% 
  filter(Mesocosm <9) %>% 
  select(Compound) %>% 
  unique() %>% 
  rename(Lipid = Compound)

mlr_aug_d1 <- mlr_aug_d1 %>% 
  mutate(Depth = 1)
mlr_aug_d2 <- mlr_aug_d2 %>% 
  mutate(Depth =2)
combined_plot <- bind_rows(mlr_aug_d1, mlr_aug_d2)
test_rel_ab <- merge(combined_plot, rel_ab, by = "Lipid")
integrated <- merge(mlr_aug_d, rel_ab, by = "Lipid")

mlr_aug_p12 <- mlr_aug_p12 %>% 
  mutate(phase =12)
mlr_aug_p34 <- mlr_aug_p34 %>% 
  mutate(phase =34)
combined_phases <- bind_rows(mlr_aug_p12, mlr_aug_p34)
phases_rel_ab <- merge(combined_phases, rel_ab, by = "Lipid")

mlr_aug_p12_d1 <- mlr_aug_p12_d1 %>% 
  mutate(phase =12)
mlr_aug_p34_d1 <- mlr_aug_p34_d1 %>% 
  mutate(phase =34)
combined_phases_d1 <- bind_rows(mlr_aug_p12_d1, mlr_aug_p34_d1)
phases_rel_ab_d1 <- merge(combined_phases_d1, rel_ab, by = "Lipid")

############ Seb C Plots


ggplot(filter(mlr_aug_d1, p_value <.05), aes(Coef, Nutrient, color = Class))+
  geom_point(size = 2)+
  scale_color_manual(values = cbPalette)+
  geom_vline(xintercept=0)
ggsave("mlrs_all_classes_june2021.png", width = 7, height = 6)


ggplot(filter(mlr_aug_d, sig95 == T))+
  geom_point(aes(Coef, Nutrient, color = Class), size = 2.5)+
  geom_text_repel(aes(Coef, Nutrient, label = ca_db), size = 3, segment.size = 0.25, max.overlaps = 100)+
  scale_color_manual(values = cbPalette)+
  geom_vline(xintercept=0)+
  facet_grid(~Class)+
  theme(text = element_text(size = 18))+
  ggsave("GLMs_core_scale_noprym_Aug2021.png", width = 14, height = 6)


ggplot(filter(integrated, sig95 == T))+
  geom_point(alpha = .4, aes(Coef, Phyto, color = Class, size = r_sq_full))+
  geom_text_repel(aes(Coef, Phyto, label = ca_db), size = 3, segment.size = 0.25, max.overlaps = 100)+
  scale_color_manual("IPL Class", values = cbPalette)+
  geom_vline(xintercept=0)+
  facet_grid(~Class)+
  scale_x_continuous(breaks = seq(-1,1, 1), limits = c(-1.1,1.1))+
  # scale_y_discrete(labels = c("Light", expression(NH[4]), expression(NO[2]), expression(NO[3]), "N:P", expression(O[2]), "pH", expression(PO[4]), expression(SiO[2]), "Temp"))+
  xlab("Regression Coeficient")+ ylab("Physiochemical Variable")+
  theme(text = element_text(size = 18))+
  scale_size(name = expression(R^2), breaks = seq(.1,.35,.05), labels = seq(.1,.35,.05), range = c(2,7))
ggsave("mlr_bio_integrated.png", width = 12, height = 10)

ggplot(filter(combined_phases_d1, sig95 == T))+
  geom_point(alpha = .4, aes(Coef, Phyto, color = Class, size = r_sq_full))+
  geom_text_repel(aes(Coef, Phyto, label = ca_db), size = 3, segment.size = 0.25, max.overlaps = 100)+
  scale_color_manual("IPL Class", values = cbPalette)+
  geom_vline(xintercept=0)+
  facet_grid(phase~Class)+
  scale_x_continuous(breaks = seq(-1,1, 1), limits = c(-1,1))+
  # scale_y_discrete(labels = c("Light", expression(NH[4]), expression(NO[2]), expression(NO[3]), "N:P", expression(O[2]), "pH", expression(PO[4]), expression(SiO[2]), "Temp"))+
  xlab("Regression Coeficient")+ ylab("Physiochemical Variable")+
  theme(text = element_text(size = 18))+
  scale_size(name = expression(R^2), breaks = seq(.1,.35,.05), labels = seq(.1,.35,.05), range = c(2,7))
ggsave("mlr_bio_integrated_byphase_sigonly_bydepth.png", width = 12, height = 10)

ggplot(filter(phases_rel_ab, sig95 == T))+
  geom_point(alpha = .4, aes(Coef, Phyto, color = Class, size = r_sq_full))+
  geom_text_repel(aes(Coef, Phyto, label = ca_db), size = 3, segment.size = 0.25, max.overlaps = 100)+
  scale_color_manual("IPL Class", values = cbPalette)+
  geom_vline(xintercept=0)+
  facet_grid(phase~Class)+
  scale_x_continuous(breaks = seq(-1,1, 1), limits = c(-1.1,1.1))+
  # scale_y_discrete(labels = c("Light", expression(NH[4]), expression(NO[2]), expression(NO[3]), "N:P", expression(O[2]), "pH", expression(PO[4]), expression(SiO[2]), "Temp"))+
  xlab("Regression Coeficient")+ ylab("Physiochemical Variable")+
  theme(text = element_text(size = 18))+
  scale_size(name = expression(R^2), breaks = seq(.1,.35,.05), labels = seq(.1,.35,.05), range = c(2,7))
ggsave("mlr_bio_integrated_byphase_sigonly.png", width = 12, height = 10)

ggplot(filter(mlr_aug_d2, sig95 == T))+
  geom_point(aes(Coef, Nutrient, color = Class), size = 2.5)+
  geom_text_repel(aes(Coef, Nutrient, label = ca_db), size = 3, segment.size = 0.25, max.overlaps = 100)+
  scale_color_manual(values = cbPalette)+
  geom_vline(xintercept=0)+
  facet_grid(~Class)+
  theme(text = element_text(size = 18)
  )+
  ggsave("GLMs_deep_compounds_core_june2021.png", width = 14, height = 6)


################## Combined mlr plots



ggplot(filter(combined_plot, sig95 == T))+
  geom_point(aes(Coef, Phyto, color = Class, size = r_sq_full))+
  geom_text_repel(aes(Coef, Phyto, label = ca_db), size = 3, segment.size = 0.25, max.overlaps = 100)+
  scale_color_manual("IPL Class", values = cbPalette)+
  # geom_vline(xintercept=0)+
  facet_grid(Depth~Class)+
  scale_x_continuous(breaks = seq(-.5,.5, 1), limits = c(-0,1))+
  # scale_y_discrete(labels = c("Light", expression(NH[4]), expression(NO[2]), expression(NO[3]), "N:P", expression(O[2]), "pH", expression(PO[4]), expression(SiO[2]), "Temp"))+
  xlab("Regression Coeficient")+ ylab("Physiochemical Variable")+
  theme(text = element_text(size = 18))+
  scale_size(name = expression(R^2), breaks = seq(.1,.35,.05), labels = seq(.1,.35,.05), range = c(2,7))
  ggsave("biomlrs_combined_plot_sigonly_Sept2021.png", width = 12, height = 10)


ggplot(filter(test_rel_ab, sig95 == T))+
  geom_point(alpha = .4, aes(Coef, Nutrient, color = Class, size = r_sq_diff))+
  geom_text_repel(aes(Coef, Nutrient, label = ca_db), size = 3, segment.size = 0.25, max.overlaps = 100)+
  scale_color_manual("IPL Class", values = cbPalette)+
  geom_vline(xintercept=0)+
  facet_grid(Depth~Class)+
  scale_x_continuous(breaks = seq(-1,1, 1), limits = c(-1.1,1.1))+
  scale_y_discrete(labels = c("Light", expression(NH[4]), expression(NO[2]), expression(NO[3]), "N:P", expression(O[2]), "pH", expression(PO[4]), expression(SiO[2]), "Temp"))+
  xlab("Regression Coeficient")+ ylab("Physiochemical Variable")+
  theme(text = element_text(size = 18))+
  scale_size(name = expression(R^2), breaks = seq(.1,.35,.05), labels = seq(.1,.35,.05), range = c(2,7))
  ggsave("mlrs_combined_plot_sigonly_Sept2021.png", width = 12, height = 10)
  
 
  
  ###############
  ggplot(filter(combined_plot, sig95 == T), aes(Coef, Nutrient))+
    geom_point(aes(color = r_sq_diff), size = 4)+
    scale_color_gradient2()+
    geom_text_repel(aes(Coef, Nutrient, label = ca_db), size = 3, segment.size = 0.25, max.overlaps = 100)+
    geom_vline(xintercept=0)+
    facet_grid(Depth~Class)+
    scale_x_continuous(breaks = seq(-.5,.5, 1), limits = c(-1,1))+
    scale_y_discrete(labels = c("Light", expression(NH[4]), expression(NO[2]), expression(NO[3]), "N:P", expression(O[2]), "pH", expression(PO[4]), expression(SiO[2]), "Temp"))+
    xlab("Regression Coeficient")+ ylab("Physiochemical Variable")+
    theme(text = element_text(size = 18))+
    scale_size(name = expression(R^2), breaks = seq(.1,.35,.05), labels = seq(.1,.35,.05), range = c(1,4))
  ggsave("test.png", width = 12, height = 10)
