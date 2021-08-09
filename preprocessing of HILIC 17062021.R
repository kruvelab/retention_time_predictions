library(tidyverse)
library(htmlwidgets)
library(compound)
library(rcdk)
library(OrgMassSpecR)
library(webchem)
library(rJava)
library(rcdklibs)
library(rcdk)
library(grid)
library(OrgMassSpecR)
library(enviPat)

# Make a toy dataset

#ToDo
#filter out M+Na and other ions
#fragmendid - OK?
#isotop correction - OK
#ret aegade analyys, ega pole hälbijaid
#compound structure to the pag
#parem sortimine x-teljel

#function to find molecular mass from SMILES code
fn_molecular_mass <- function(smiles){
  #convert SMILES to chemical formula
  molecule <- parse.smiles(smiles)[[1]]
  formula <- get.mol2formula(molecule,charge=0)
  formula <- formula@string
  #calcuate molecular weight
  MW <- MolecularWeight(formula = ListFormula(formula))
  return(MW)
}

#function to find isotope distribution parameter from SMILES code
fn_isotope_distribution <- function(smiles){
  #convert SMILES to chemical formula
  molecule <- parse.smiles(smiles)[[1]]
  formula <- get.mol2formula(molecule,charge=0)
  formula <- formula@string
  
  # Chemical formula to isotope distribution
  data(isotopes)
  pattern<-isopattern(isotopes,
                      formula,
                      threshold=0.1,
                      plotit=FALSE,
                      charge=FALSE,
                      algo=1)
  isotopes <- as.data.frame(pattern[[1]])
  isotope_dist <- as.numeric(sum(isotopes$abundance))
  return(isotope_dist)
}

linear_regression <- function(y, x) {
  if(length(y) > 3) {
    for (i in length(y):3){
      y = y[1:i]
      x = x[1:i]
      slope = summary(lm(y ~ x))$coefficients[2]
      intercept = summary(lm(y ~ x))$coefficients[1]
      residuals = (y - (slope*x +intercept))/y*100
      regression_parameters <- list("slope" = slope, "intercept" = intercept)
      if (max(abs(residuals)) < 10) {
        return(regression_parameters)
        break
      }
    }
    return(regression_parameters)
  }
  else {
    slope = summary(lm(y ~ x))$coefficients[2]
    intercept = summary(lm(y ~ x))$coefficients[1]
    regression_parameters <- list("slope" = slope, "intercept" = intercept)
    return(regression_parameters)
    
  }
}


metadata <- read_delim('Solutions.csv',
                       delim = ',',
                       col_names = TRUE)
metadata <- metadata %>% select(-SMILES)
#NEEDS CHECKING:
#Chloramphenicol pH 10, 8 (formic) and 0.1% formic, ret time of M+H and M-H2O does not match for MeCH runs
#Tri-isobutyl phosphate at acetic acid with MeOH, fragment and M+H ret times do not match
#I filtered these out at the moment
#Why is there only one 0.1% formic acid?

#Chloramphenicol - water loss can be real


dataset <- dataset %>%
  #filter out ions which correspond to different ionization mechanism
  filter(Ion != "M+Na" & Ion != "M?" & Ion != "2M+Na") %>%
  #check if all of the ions have the same ret time and filter out the ones where dif is more then 0.4 min
  group_by(File_name, Compound_name, SMILES, Column, Organic_modifier, pH, Buffer, Conc) %>%
  mutate(RT_res = ret_time - median(ret_time)) %>%
  ungroup() %>%
  filter(abs(RT_res) < 0.2) %>%
  dplyr::select(-RT_res) %>%
  #check if all of the peaks for one compound accross concentrations have same ret time and filter out the ones where it is not true
  group_by(File_name, Compound_name, SMILES, Column, Organic_modifier, pH, Buffer) %>%
  mutate(RT_res = abs(ret_time - median(ret_time))) %>%
  ungroup() %>%
  filter(abs(RT_res) < 0.2) %>%
  dplyr::select(-RT_res) %>%
  #combine the peak areas of different ions belonging to the same compound
  group_by(File_name, Compound_name, SMILES, Column, Organic_modifier, pH, Buffer, Conc) %>%
  summarise(Peak_Area = sum(Peak_Area),
            ret_time = mean(ret_time),
            `S/N` = max(`S/N`)) %>%
  ungroup()
dataset <- dataset %>%
  #check if there are any compounds with less then 3 data points and if so, remove these
  group_by(File_name, Compound_name, SMILES, Column, Organic_modifier, pH, Buffer) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  filter(count > 2) %>%
  #add MW and isotope correction
  group_by(SMILES) %>%
  mutate(MW = fn_molecular_mass(SMILES)) %>%
  ungroup() %>%
  group_by(SMILES) %>%
  mutate(IC = fn_isotope_distribution(SMILES)) %>%
  ungroup() %>%
  #add the concentrations
  left_join(metadata) %>%
  mutate(conc_M = as.numeric(mg) /1000 / MW / (solv_g / 0.76 / 1000) *pip_ul / 10 *Conc / 20,
         Peak_Area_IC = IC*Peak_Area)
dataset <- dataset %>%
  #add MW filter 100 because 3quads are known to have very low sensitivity on ions wit m/z below 100
  filter(MW > 100) %>%
  dplyr::select(File_name, Compound_name, Column, Organic_modifier, pH, Buffer, SMILES, ret_time, Peak_Area_IC, Peak_Area, `S/N`, conc_M, MW) %>%
  #remove everything where some data is missing
  na.omit() %>%
  #start calculating the slopes
  group_by(Compound_name, Column, Organic_modifier, pH, Buffer) %>%
  mutate(slope = linear_regression(Peak_Area_IC, conc_M)$slope,
         intercept = linear_regression(Peak_Area_IC, conc_M)$intercept) %>%
  ungroup() %>%
  #lets add the residuals as well
  mutate(residual = (Peak_Area_IC - (slope*conc_M + intercept)) / Peak_Area_IC *100) 

dataset <- dataset %>%
  filter(abs(residual) < 20) %>%
  group_by(File_name, Compound_name, SMILES, Column, Organic_modifier, pH, Buffer) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  #only compounds where there are more then 3 points on the graph are left into the dataset
  filter(count > 3) %>%
  dplyr::select(-count)

dataset <- dataset %>%
  group_by(File_name, Compound_name, SMILES, Column, Organic_modifier, pH, Buffer) %>%
  mutate(max_c = max(conc_M),
         min_c = min(conc_M),
         max_peak = max(Peak_Area_IC),
         min_peak = min(Peak_Area_IC)) %>%
  ungroup() %>%
  mutate(Peak_area_max_min_ratio = max_peak / min_peak,
         c_max_min_ratio = max_c / min_c,
         ratio = case_when(
           Peak_area_max_min_ratio > c_max_min_ratio ~ Peak_area_max_min_ratio / c_max_min_ratio,
           TRUE ~ c_max_min_ratio / Peak_area_max_min_ratio)
  ) %>%
  ungroup() %>%
  filter(ratio < 2) %>%
  dplyr::select(File_name, Compound_name, SMILES, Column, Organic_modifier, pH, Buffer, Peak_Area_IC, slope, 'S/N', conc_M, ret_time)

write_delim(dataset, "Filter_RP_2010.csv", delim = ",")
