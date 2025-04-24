# Liquid Chormatography Retention Time Predictions
This repository contains code for developing liquid chromatography (LC) retetnion time (RT) prediction models. 
Firstly, folder **MultiConditionRT** contains data and different models for predicting RT depending on the mobile phases by Souihi et al. [https://doi.org/10.1016/j.chroma.2022.462867].
Secondly, folder **RT_structural_isoemrs** contains materials for Kruve et al. "Retention time prediction for structural isomers in reversed-phase liquid chromatography" (in press): \
(1) *data.csv* contains cleaned data for model training; \
(2) *standardization_Mondred_description.jpyter* contains **python** code for standardization of SMILES and calculation of Mordred descriptors from mol representation; \
(3) *modelling.R* contains code used for splitting the data for training, uncertainty quantification with Mondrian conformal predictive systems, and testing; \
(4) *functions.R* contains all functions regarding data cleaning prior to model training, model training, and uncertainty quantification.
