library(caret)
library(caTools)
library(cluster)
#library(factoextra)
#library(ggbeeswarm)
#library(OrgMassSpecR)
#library(plotly)
# library(rcdk)
library(readxl)
#library(rjson)
library(tidyverse)
#library(umap)


# data cleaning functions ----
remove_correlated_nearZero = function(data = tibble(),
                                      cutoff = numeric()) {
  data = data %>%
    mutate(Lipinski = as.numeric(Lipinski),
           GhoseFilter = as.numeric(GhoseFilter))
  
  data = data %>%
    relocate(c(inchikey14, standardized_inchikey, isomeric, SMILES_enumerated, MF, RT), .after = last_col())
  
  data = data %>%
    select(-nearZeroVar(data %>% 
                          select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT)))
  
  correlation_matrix = cor(data %>% 
                             select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT))
  
  highly_correlated = findCorrelation(correlation_matrix,
                                      cutoff = cutoff)
  
  data = data %>%
    select(-all_of(highly_correlated)) 
  
  data = data %>%
    select(inchikey14, standardized_inchikey, isomeric, SMILES_enumerated, RT, MF, colnames(data))
  
  return(data)
  
}

# modelling functions ----

# training_data = trainTransformed
# test_data = testTransformed
# candidate_structures = allTransformed

model_smart_combined = function(training_data = tibble(),
                                test_data = tibble(),
                                cutoff_correlated = 0.7,
                                cross_validation_folds = 2,
                                candidate_structures = tibble()) {
  
  training_data = remove_correlated_nearZero(training_data,
                                             cutoff_correlated)
  
  #train and apply simple model
  model_simple = simple_model(training_data = training_data %>%
                                select(-c(inchikey14, standardized_inchikey, isomeric, SMILES_enumerated, MF)))
  
  training_data = training_data %>%
    mutate(RT_simple_pred = predict(model_simple,
                                    newdata = training_data),
           RT_delta = RT - RT_simple_pred)
  
  folds = groupKFold(training_data$MF, 
                     k = cross_validation_folds) 
  
  fitControl = trainControl(
    method = "boot",
    index = folds)
  
  model_rf = train(RT ~ ., 
                   data = training_data %>%
                     select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT_delta, -RT_simple_pred), 
                   method = "rf", 
                   trControl = fitControl)
  
  model_svm = train(RT ~ ., 
                    data = training_data %>%
                      select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT_delta, -RT_simple_pred), 
                    method = "svmLinear", 
                    trControl = fitControl)
  
  model_knn = train(RT ~ ., 
                    data = training_data %>%
                      select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT_delta, -RT_simple_pred), 
                    method = "knn", 
                    trControl = fitControl)
  
  
  model_xgbTREE = train(RT ~ ., 
                        data = training_data %>%
                          select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT_delta, -RT_simple_pred), 
                        method = "xgbTree", 
                        trControl = fitControl)
  
  model_rf_on_simple = train(RT_delta ~ ., 
                             data = training_data %>%
                               select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT, -RT_simple_pred),
                             method = "rf", 
                             trControl = fitControl)
  
  model_svm_on_simple = train(RT_delta ~ ., 
                              data = training_data %>%
                                select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT, -RT_simple_pred),
                              method = "svmLinear", 
                              trControl = fitControl)
  
  model_knn_on_simple = train(RT_delta ~ ., 
                              data = training_data %>%
                                select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT, -RT_simple_pred),
                              method = "knn", 
                              trControl = fitControl)
  
  model_xgbTREE_on_simple = train(RT_delta ~ ., 
                                  data = training_data %>%
                                    select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT, -RT_simple_pred),
                                  method = "xgbTree", 
                                  trControl = fitControl)
  
  training_data = training_data %>%
    mutate(RT_simple_pred = predict(model_simple, 
                                    newdata = training_data), 
           RT_delta_pred_rf = predict(model_rf_on_simple,
                                      newdata = training_data),
           RT_smart_pred_rf = RT_simple_pred + RT_delta_pred_rf,
           RT_rf_pred = predict(model_rf,
                                newdata = training_data),
           RT_delta_pred_svm = predict(model_svm_on_simple,
                                       newdata = training_data),
           RT_smart_pred_svm = RT_simple_pred + RT_delta_pred_svm,
           RT_svm_pred = predict(model_svm,
                                 newdata = training_data),
           RT_delta_pred_knn = predict(model_knn_on_simple,
                                       newdata = training_data),
           RT_smart_pred_knn = RT_simple_pred + RT_delta_pred_knn,
           RT_knn_pred = predict(model_knn,
                                 newdata = training_data),
           RT_delta_pred_xgbTREE = predict(model_xgbTREE_on_simple,
                                           newdata = training_data),
           RT_smart_pred_xgbTREE = RT_simple_pred + RT_delta_pred_xgbTREE,
           RT_xgbTREE_pred = predict(model_xgbTREE,
                                     newdata = training_data)) %>%
    select(inchikey14, RT, standardized_inchikey, isomeric, SMILES_enumerated, RT_simple_pred, RT_rf_pred, RT_svm_pred, RT_knn_pred, RT_xgbTREE_pred, RT_smart_pred_rf, RT_smart_pred_svm, RT_smart_pred_knn, RT_smart_pred_xgbTREE)
  
  test_data = test_data %>%
    mutate(Lipinski = as.numeric(Lipinski),
           GhoseFilter = as.numeric(GhoseFilter)) 
  
  test_data = test_data %>%
    mutate(RT_simple_pred = predict(model_simple, 
                                    newdata = test_data),
           RT_delta = RT - RT_simple_pred, 
           RT_delta_pred_rf = predict(model_rf_on_simple,
                                      newdata = test_data),
           RT_smart_pred_rf = RT_simple_pred + RT_delta_pred_rf,
           RT_rf_pred = predict(model_rf,
                                newdata = test_data),
           RT_delta_pred_svm = predict(model_svm_on_simple,
                                       newdata = test_data),
           RT_smart_pred_svm = RT_simple_pred + RT_delta_pred_svm,
           RT_svm_pred = predict(model_svm,
                                 newdata = test_data),
           RT_delta_pred_knn = predict(model_knn_on_simple,
                                       newdata = test_data),
           RT_smart_pred_knn = RT_simple_pred + RT_delta_pred_knn,
           RT_knn_pred = predict(model_knn,
                                 newdata = test_data),
           RT_delta_pred_xgbTREE = predict(model_xgbTREE_on_simple,
                                           newdata = test_data),
           RT_smart_pred_xgbTREE = RT_simple_pred + RT_delta_pred_xgbTREE,
           RT_xgbTREE_pred = predict(model_xgbTREE,
                                     newdata = test_data)) %>%
    select(inchikey14, standardized_inchikey, RT, isomeric, SMILES_enumerated, RT_simple_pred, RT_rf_pred, RT_svm_pred, RT_knn_pred, RT_xgbTREE_pred, RT_smart_pred_rf, RT_smart_pred_svm, RT_smart_pred_knn, RT_smart_pred_xgbTREE)
  
  
  candidate_structures = candidate_structures %>%
    mutate(Lipinski = as.numeric(Lipinski),
           GhoseFilter = as.numeric(GhoseFilter)) 
  
  candidate_structures = candidate_structures %>%
    mutate(RT_simple_pred = predict(model_simple, 
                                    newdata = candidate_structures), 
           RT_delta_pred_rf = predict(model_rf_on_simple,
                                      newdata = candidate_structures),
           RT_smart_pred_rf = RT_simple_pred + RT_delta_pred_rf,
           RT_rf_pred = predict(model_rf,
                                newdata = candidate_structures),
           RT_delta_pred_svm = predict(model_svm_on_simple,
                                       newdata = candidate_structures),
           RT_smart_pred_svm = RT_simple_pred + RT_delta_pred_svm,
           RT_svm_pred = predict(model_svm,
                                 newdata = candidate_structures),
           RT_delta_pred_knn = predict(model_knn_on_simple,
                                       newdata = candidate_structures),
           RT_smart_pred_knn = RT_simple_pred + RT_delta_pred_knn,
           RT_knn_pred = predict(model_knn,
                                 newdata = candidate_structures),
           RT_delta_pred_xgbTREE = predict(model_xgbTREE_on_simple,
                                           newdata = candidate_structures),
           RT_smart_pred_xgbTREE = RT_simple_pred + RT_delta_pred_xgbTREE,
           RT_xgbTREE_pred = predict(model_xgbTREE,
                                     newdata = candidate_structures)) %>%
    select(type, standardized_SMILES, RT_simple_pred, RT_rf_pred, RT_svm_pred, RT_knn_pred, RT_xgbTREE_pred, RT_smart_pred_rf, RT_smart_pred_svm, RT_smart_pred_knn, RT_smart_pred_xgbTREE)
  
  
  all = list("training_data" = training_data, 
             "test_data" = test_data, 
             "candidate_structures" = candidate_structures)
  
  return(all)
}

simple_model = function(training_data = tibble()) {
  data = training_data %>% 
    select(RT, everything()) %>%
    as.matrix() %>%
    cor() %>%
    as.data.frame() %>%
    rownames_to_column(var = 'var') %>%
    select(var, RT) %>%
    filter(var != "RT") %>%
    filter(RT == max(RT)) 
  
  correlated_feature = data %>%
    spread(key = var, value = RT) 
  
  training_data = training_data %>%
    select(RT, colnames(correlated_feature))
  
  names(training_data)[2] = colnames(correlated_feature)
  #names(training_data)[2] ="feature"
  
  model = lm(RT ~ .,
             data = training_data)
  
  return(model)
}


model_final = function(training_data = tibble(),
                       test_data = tibble(),
                       cutoff_correlated = 0.7,
                       cross_validation_folds = 2,
                       candidate_structures = tibble()) {
  
  training_data = remove_correlated_nearZero(training_data,
                                             cutoff_correlated)
  
  #train and apply simple model
  model_simple = simple_model(training_data = training_data %>%
                                select(-c(inchikey14, standardized_inchikey, isomeric, SMILES_enumerated, MF)))
  
  training_data = training_data %>%
    mutate(RT_simple_pred = predict(model_simple,
                                    newdata = training_data),
           RT_delta = RT - RT_simple_pred)
  
  folds = groupKFold(training_data$MF, 
                     k = cross_validation_folds) 
  
  fitControl = trainControl(
    method = "boot",
    index = folds)
  
  model_rf = train(RT ~ ., 
                   data = training_data %>%
                     select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT_delta, -RT_simple_pred), 
                   method = "rf", 
                   trControl = fitControl)
  
  model_xgbTREE = train(RT ~ ., 
                        data = training_data %>%
                          select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT_delta, -RT_simple_pred), 
                        method = "xgbTree", 
                        trControl = fitControl)
  
  model_rf_on_simple = train(RT_delta ~ ., 
                             data = training_data %>%
                               select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT, -RT_simple_pred),
                             method = "rf", 
                             trControl = fitControl)
  
  
  model_xgbTREE_on_simple = train(RT_delta ~ ., 
                                  data = training_data %>%
                                    select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT, -RT_simple_pred),
                                  method = "xgbTree", 
                                  trControl = fitControl)
  
  training_data = training_data %>%
    mutate(RT_simple_pred = predict(model_simple, 
                                    newdata = training_data), 
           RT_delta_pred_rf = predict(model_rf_on_simple,
                                      newdata = training_data),
           RT_smart_pred_rf = RT_simple_pred + RT_delta_pred_rf,
           RT_rf_pred = predict(model_rf,
                                newdata = training_data),
           RT_delta_pred_xgbTREE = predict(model_xgbTREE_on_simple,
                                           newdata = training_data),
           RT_smart_pred_xgbTREE = RT_simple_pred + RT_delta_pred_xgbTREE,
           RT_xgbTREE_pred = predict(model_xgbTREE,
                                     newdata = training_data)) %>%
    select(inchikey14, RT, standardized_inchikey, isomeric, SMILES_enumerated, RT_simple_pred, RT_rf_pred, RT_xgbTREE_pred, RT_smart_pred_rf, RT_smart_pred_xgbTREE)
  
  test_data = test_data %>%
    mutate(Lipinski = as.numeric(Lipinski),
           GhoseFilter = as.numeric(GhoseFilter)) 
  
  test_data = test_data %>%
    mutate(RT_simple_pred = predict(model_simple, 
                                    newdata = test_data),
           RT_delta = RT - RT_simple_pred, 
           RT_delta_pred_rf = predict(model_rf_on_simple,
                                      newdata = test_data),
           RT_smart_pred_rf = RT_simple_pred + RT_delta_pred_rf,
           RT_rf_pred = predict(model_rf,
                                newdata = test_data),
           RT_delta_pred_xgbTREE = predict(model_xgbTREE_on_simple,
                                           newdata = test_data),
           RT_smart_pred_xgbTREE = RT_simple_pred + RT_delta_pred_xgbTREE,
           RT_xgbTREE_pred = predict(model_xgbTREE,
                                     newdata = test_data)) %>%
    select(inchikey14, RT, standardized_inchikey, isomeric, SMILES_enumerated, RT_simple_pred, RT_rf_pred, RT_xgbTREE_pred, RT_smart_pred_rf, RT_smart_pred_xgbTREE)
  
  
  candidate_structures = candidate_structures %>%
    mutate(Lipinski = as.numeric(Lipinski),
           GhoseFilter = as.numeric(GhoseFilter)) 
  
  candidate_structures = candidate_structures %>%
    mutate(RT_simple_pred = predict(model_simple, 
                                    newdata = candidate_structures), 
           RT_delta_pred_rf = predict(model_rf_on_simple,
                                      newdata = candidate_structures),
           RT_smart_pred_rf = RT_simple_pred + RT_delta_pred_rf,
           RT_rf_pred = predict(model_rf,
                                newdata = candidate_structures),
           RT_delta_pred_xgbTREE = predict(model_xgbTREE_on_simple,
                                           newdata = candidate_structures),
           RT_smart_pred_xgbTREE = RT_simple_pred + RT_delta_pred_xgbTREE,
           RT_xgbTREE_pred = predict(model_xgbTREE,
                                     newdata = candidate_structures)) %>%
    select(id, correct, RT, standardized_inchikey, standardized_SMILES, RT_simple_pred, RT_rf_pred, RT_xgbTREE_pred, RT_smart_pred_rf, RT_smart_pred_xgbTREE)
  
  
  all = list("training_data" = training_data, 
             "test_data" = test_data, 
             "candidate_structures" = candidate_structures)
  
  return(all)
}


# cf for test set and pubchem candidates -----

model_mondrean_cf = function(training_data = tibble(),
                             calibration_data = tibble(),
                             test_data = tibble(),
                             cutoff_correlated = 0.7,
                             cross_validation_folds = 2,
                             candidate_structures = tibble()) {
  
  training_data = remove_correlated_nearZero(training_data,
                                             cutoff_correlated)
  
  #train and apply simple model
  model_simple = simple_model(training_data = training_data %>%
                                select(-c(inchikey14, standardized_inchikey, isomeric, SMILES_enumerated, MF)))
  
  training_data = training_data %>%
    mutate(RT_simple_pred = predict(model_simple,
                                    newdata = training_data),
           RT_delta = RT - RT_simple_pred)
  
  folds = groupKFold(training_data$MF, 
                     k = cross_validation_folds) 
  
  fitControl = trainControl(
    method = "boot",
    index = folds)
  
  model_rf = train(RT ~ ., 
                   data = training_data %>%
                     select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT_delta, -RT_simple_pred), 
                   method = "rf", 
                   trControl = fitControl)
  
  model_xgbTREE = train(RT ~ ., 
                        data = training_data %>%
                          select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT_delta, -RT_simple_pred), 
                        method = "xgbTree", 
                        trControl = fitControl)
  
  model_rf_on_simple = train(RT_delta ~ ., 
                             data = training_data %>%
                               select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT, -RT_simple_pred),
                             method = "rf", 
                             trControl = fitControl)
  
  
  model_xgbTREE_on_simple = train(RT_delta ~ ., 
                                  data = training_data %>%
                                    select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT, -RT_simple_pred),
                                  method = "xgbTree", 
                                  trControl = fitControl)
  
  training_data = training_data %>%
    mutate(RT_simple_pred = predict(model_simple, 
                                    newdata = training_data), 
           RT_delta_pred_rf = predict(model_rf_on_simple,
                                      newdata = training_data),
           RT_smart_pred_rf = RT_simple_pred + RT_delta_pred_rf,
           RT_rf_pred = predict(model_rf,
                                newdata = training_data),
           RT_delta_pred_xgbTREE = predict(model_xgbTREE_on_simple,
                                           newdata = training_data),
           RT_smart_pred_xgbTREE = RT_simple_pred + RT_delta_pred_xgbTREE,
           RT_xgbTREE_pred = predict(model_xgbTREE,
                                     newdata = training_data)) %>%
    select(inchikey14, RT, standardized_inchikey, isomeric, SMILES_enumerated, RT_simple_pred, RT_rf_pred, RT_xgbTREE_pred, RT_smart_pred_rf, RT_smart_pred_xgbTREE)
  
  
  calibration_data = calibration_data %>%
    mutate(Lipinski = as.numeric(Lipinski),
           GhoseFilter = as.numeric(GhoseFilter)) 
  
  calibration_data = calibration_data %>%
    mutate(RT_simple_pred = predict(model_simple, 
                                    newdata = calibration_data),
           RT_delta = RT - RT_simple_pred, 
           RT_delta_pred_rf = predict(model_rf_on_simple,
                                      newdata = calibration_data),
           RT_smart_pred_rf = RT_simple_pred + RT_delta_pred_rf,
           RT_rf_pred = predict(model_rf,
                                newdata = calibration_data),
           RT_delta_pred_xgbTREE = predict(model_xgbTREE_on_simple,
                                           newdata = calibration_data),
           RT_smart_pred_xgbTREE = RT_simple_pred + RT_delta_pred_xgbTREE,
           RT_xgbTREE_pred = predict(model_xgbTREE,
                                     newdata = calibration_data)) %>%
    select(inchikey14, RT, standardized_inchikey, isomeric, SMILES_enumerated, RT_simple_pred, RT_rf_pred, RT_xgbTREE_pred, RT_smart_pred_rf, RT_smart_pred_xgbTREE)
  
  mondrean_summary = mondrean_cf_one(calibration_data = calibration_data)
  
  test_data = test_data %>%
    mutate(Lipinski = as.numeric(Lipinski),
           GhoseFilter = as.numeric(GhoseFilter)) 
  
  test_data = test_data %>%
    mutate(RT_simple_pred = predict(model_simple, 
                                    newdata = test_data),
           RT_delta = RT - RT_simple_pred, 
           RT_delta_pred_rf = predict(model_rf_on_simple,
                                      newdata = test_data),
           RT_smart_pred_rf = RT_simple_pred + RT_delta_pred_rf,
           RT_rf_pred = predict(model_rf,
                                newdata = test_data),
           RT_delta_pred_xgbTREE = predict(model_xgbTREE_on_simple,
                                           newdata = test_data),
           RT_smart_pred_xgbTREE = RT_simple_pred + RT_delta_pred_xgbTREE,
           RT_xgbTREE_pred = predict(model_xgbTREE,
                                     newdata = test_data)) %>%
    select(inchikey14, RT, standardized_inchikey, isomeric, SMILES_enumerated, RT_simple_pred, RT_rf_pred, RT_xgbTREE_pred, RT_smart_pred_rf, RT_smart_pred_xgbTREE)
  
  test_data = mondrean_apply_one(mondrean_summary = mondrean_summary,
                                 new_data = test_data)
  
  candidate_structures = candidate_structures %>%
    mutate(Lipinski = as.numeric(Lipinski),
           GhoseFilter = as.numeric(GhoseFilter)) 
  
  candidate_structures = candidate_structures %>%
    mutate(RT_simple_pred = predict(model_simple, 
                                    newdata = candidate_structures), 
           RT_delta_pred_rf = predict(model_rf_on_simple,
                                      newdata = candidate_structures),
           RT_smart_pred_rf = RT_simple_pred + RT_delta_pred_rf,
           RT_rf_pred = predict(model_rf,
                                newdata = candidate_structures),
           RT_delta_pred_xgbTREE = predict(model_xgbTREE_on_simple,
                                           newdata = candidate_structures),
           RT_smart_pred_xgbTREE = RT_simple_pred + RT_delta_pred_xgbTREE,
           RT_xgbTREE_pred = predict(model_xgbTREE,
                                     newdata = candidate_structures)) %>%
    select(MF, RT, standardized_inchikey, standardized_SMILES, correct, RT_simple_pred, RT_rf_pred, RT_xgbTREE_pred, RT_smart_pred_rf, RT_smart_pred_xgbTREE)
  
  candidate_structures = mondrean_apply_one(mondrean_summary = mondrean_summary,
                                            new_data = candidate_structures)
  
  all = list("training_data" = training_data, 
             "test_data" = test_data, 
             "candidate_structures" = candidate_structures,
             "mondrean_summary" = mondrean_summary)
  
  return(all)
}


mondrean_cf = function(calibration_data = tibble(),
                       p_value = 0.95) {
  
  calibration_data = calibration_data %>%
    gather(key = "model", value = "RT_pred", 6:10) %>%
    unique()
  
  scatter_plot = ggplot(data = calibration_data, 
                        mapping = aes(x=RT, 
                                      y=RT_pred)) +
    geom_point() +
    facet_wrap(~model)
  
  print(scatter_plot)
  
  calibration_data = calibration_data %>%
    group_by(model) %>%
    mutate(q33 = quantile(RT_pred, probs = 0.33),
           q67 = quantile(RT_pred, probs = 0.67)) %>%
    ungroup()
  
  calibration_data = calibration_data %>%
    mutate(taxonomy = case_when(
      RT_pred < q33 ~ "set_low",
      RT_pred < q67 ~ "set_medium",
      TRUE ~ "set_high")) 
  
  calibration_data = calibration_data %>%
    mutate(delta_RT = RT - RT_pred)
  
  density_plot = ggplot(data = calibration_data, 
                        mapping = aes(x=delta_RT, 
                                      color=taxonomy)) +
    geom_density() +
    facet_wrap(~model)
  
  print(density_plot)
  
  summary = calibration_data %>%
    group_by(model, taxonomy, q33, q67) %>%
    summarize(delta_RT_low = quantile(delta_RT, probs = (1-p_value)/2),
              delta_RT_high = quantile(delta_RT, probs = p_value + (1-p_value)/2)) %>%
    ungroup()
  
  return(summary)
  
}

mondrean_apply = function(mondrean_summary = tibble(),
                          new_data = tibble()) {
  
  new_data = new_data %>%
    gather(key = "model", value = "RT_pred", 6:10) %>%
    unique()
  
  new_data = new_data %>%
    left_join(mondrean_summary %>%
                select(model, q33, q67) %>%
                unique())
  
  new_data = new_data %>%
    mutate(taxonomy = case_when(
      RT_pred < q33 ~ "set_low",
      RT_pred < q67 ~ "set_medium",
      TRUE ~ "set_high")) 
  
  new_data = new_data %>% 
    left_join(mondrean_summary %>%
                select(model, taxonomy, delta_RT_low, delta_RT_high)) %>%
    mutate(RT_pred_low = RT_pred + delta_RT_low,
           RT_pred_high = RT_pred + delta_RT_high)
  
  scatter_plot = ggplot(data = new_data) +
    geom_abline(intercept = 0, slope = 1) +
    geom_segment(mapping = aes(x=RT,
                               y=RT_pred_low,
                               xend=RT,
                               yend=RT_pred_high,
                               color = standardized_inchikey)) +
    facet_wrap(~model) +
    theme(legend.position = "none")
  
  print(scatter_plot)
  
  return(new_data)
  
}

#cf and monte carlo sampling for NTS samples -----

model_mondrean_cf_mc = function(training_data = tibble(),
                                calibration_data = tibble(),
                                cutoff_correlated = 0.7,
                                cross_validation_folds = 2,
                                candidate_structures = tibble()) {
  
  training_data = remove_correlated_nearZero(training_data,
                                             cutoff_correlated)
  
  #train and apply simple model
  model_simple = simple_model(training_data = training_data %>%
                                select(-c(inchikey14, standardized_inchikey, isomeric, SMILES_enumerated, MF)))
  
  training_data = training_data %>%
    mutate(RT_simple_pred = predict(model_simple,
                                    newdata = training_data),
           RT_delta = RT - RT_simple_pred)
  
  folds = groupKFold(training_data$MF, 
                     k = cross_validation_folds) 
  
  fitControl = trainControl(
    method = "boot",
    index = folds)
  
  model_rf = train(RT ~ ., 
                   data = training_data %>%
                     select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT_delta, -RT_simple_pred), 
                   method = "rf", 
                   trControl = fitControl)
  
  model_xgbTREE = train(RT ~ ., 
                        data = training_data %>%
                          select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT_delta, -RT_simple_pred), 
                        method = "xgbTree", 
                        trControl = fitControl)
  
  model_rf_on_simple = train(RT_delta ~ ., 
                             data = training_data %>%
                               select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT, -RT_simple_pred),
                             method = "rf", 
                             trControl = fitControl)
  
  
  model_xgbTREE_on_simple = train(RT_delta ~ ., 
                                  data = training_data %>%
                                    select(-inchikey14, -standardized_inchikey, -isomeric, -SMILES_enumerated, -MF, -RT, -RT_simple_pred),
                                  method = "xgbTree", 
                                  trControl = fitControl)
  
  training_data = training_data %>%
    mutate(RT_simple_pred = predict(model_simple, 
                                    newdata = training_data), 
           RT_delta_pred_rf = predict(model_rf_on_simple,
                                      newdata = training_data),
           RT_smart_pred_rf = RT_simple_pred + RT_delta_pred_rf,
           RT_rf_pred = predict(model_rf,
                                newdata = training_data),
           RT_delta_pred_xgbTREE = predict(model_xgbTREE_on_simple,
                                           newdata = training_data),
           RT_smart_pred_xgbTREE = RT_simple_pred + RT_delta_pred_xgbTREE,
           RT_xgbTREE_pred = predict(model_xgbTREE,
                                     newdata = training_data)) %>%
    select(inchikey14, RT, standardized_inchikey, isomeric, SMILES_enumerated, RT_simple_pred, RT_rf_pred, RT_xgbTREE_pred, RT_smart_pred_rf, RT_smart_pred_xgbTREE)
  
  
  calibration_data = calibration_data %>%
    mutate(Lipinski = as.numeric(Lipinski),
           GhoseFilter = as.numeric(GhoseFilter)) 
  
  calibration_data = calibration_data %>%
    mutate(RT_simple_pred = predict(model_simple, 
                                    newdata = calibration_data),
           RT_delta = RT - RT_simple_pred, 
           RT_delta_pred_rf = predict(model_rf_on_simple,
                                      newdata = calibration_data),
           RT_smart_pred_rf = RT_simple_pred + RT_delta_pred_rf,
           RT_rf_pred = predict(model_rf,
                                newdata = calibration_data),
           RT_delta_pred_xgbTREE = predict(model_xgbTREE_on_simple,
                                           newdata = calibration_data),
           RT_smart_pred_xgbTREE = RT_simple_pred + RT_delta_pred_xgbTREE,
           RT_xgbTREE_pred = predict(model_xgbTREE,
                                     newdata = calibration_data)) %>%
    select(inchikey14, RT, standardized_inchikey, isomeric, SMILES_enumerated, RT_simple_pred, RT_rf_pred, RT_xgbTREE_pred, RT_smart_pred_rf, RT_smart_pred_xgbTREE)
  
  mondrean_summary = mondrean_cf_mc(calibration_data = calibration_data)
  
  candidate_structures = candidate_structures %>%
    mutate(Lipinski = as.numeric(Lipinski),
           GhoseFilter = as.numeric(GhoseFilter)) 
  
  candidate_structures = candidate_structures %>%
    mutate(RT_simple_pred = predict(model_simple, 
                                    newdata = candidate_structures), 
           RT_delta_pred_rf = predict(model_rf_on_simple,
                                      newdata = candidate_structures),
           RT_smart_pred_rf = RT_simple_pred + RT_delta_pred_rf,
           RT_rf_pred = predict(model_rf,
                                newdata = candidate_structures),
           RT_delta_pred_xgbTREE = predict(model_xgbTREE_on_simple,
                                           newdata = candidate_structures),
           RT_smart_pred_xgbTREE = RT_simple_pred + RT_delta_pred_xgbTREE,
           RT_xgbTREE_pred = predict(model_xgbTREE,
                                     newdata = candidate_structures)) %>%
    select(id, sample, correct, MF, standardized_SMILES, standardized_inchikey, candidate_type, RT, RT_simple_pred, RT_rf_pred, RT_xgbTREE_pred, RT_smart_pred_rf, RT_smart_pred_xgbTREE)
  
  candidate_structures = mondrean_mc_apply(mondrean_summary = mondrean_summary,
                                           new_data = candidate_structures)
  
  all = list("training_data" = training_data, 
             "candidate_structures" = candidate_structures,
             "mondrean_summary" = mondrean_summary)
  
  return(all)
}

mondrean_cf_mc = function(calibration_data = tibble(),
                          p_value = 0.95) {
  
  calibration_data = calibration_data %>%
    gather(key = "model", value = "RT_pred", 6:10) %>%
    unique()
  
  scatter_plot = ggplot(data = calibration_data, 
                        mapping = aes(x=RT, 
                                      y=RT_pred)) +
    geom_point() +
    facet_wrap(~model)
  
  print(scatter_plot)
  
  calibration_data = calibration_data %>%
    group_by(model) %>%
    mutate(q33 = quantile(RT_pred, probs = 0.33),
           q67 = quantile(RT_pred, probs = 0.67)) %>%
    ungroup()
  
  calibration_data = calibration_data %>%
    mutate(taxonomy = case_when(
      RT_pred < q33 ~ "set_low",
      RT_pred < q67 ~ "set_medium",
      TRUE ~ "set_high")) 
  
  calibration_data = calibration_data %>%
    mutate(delta_RT = RT - RT_pred) %>%
    select(q33, q67, taxonomy, model, delta_RT)
  
  density_plot = ggplot(data = calibration_data, 
                        mapping = aes(x=delta_RT, 
                                      color=taxonomy)) +
    geom_density() +
    facet_wrap(~model)
  
  print(density_plot)
  
  return(calibration_data)
  
}

mondrean_mc_apply = function(mondrean_summary = tibble(),
                             new_data = tibble()) {
  
  new_data = new_data %>%
    gather(key = "model", value = "RT_pred", 9:13) %>%
    unique()
  
  new_data = new_data %>%
    left_join(mondrean_summary %>%
                select(model, q33, q67) %>%
                unique())
  
  new_data = new_data %>%
    mutate(taxonomy = case_when(
      RT_pred < q33 ~ "set_low",
      RT_pred < q67 ~ "set_medium",
      TRUE ~ "set_high")) 
  
  new_data = new_data %>% 
    left_join(mondrean_summary) %>%
    mutate(RT_pred_mc = RT_pred + delta_RT)
  
  mc_data = monte_carlo(data = new_data)
  
  scatter_plot = ggplot(data = mc_data) +
    geom_abline(intercept = 0, slope = 1) +
    geom_segment(mapping = aes(x=RT,
                               y=RT_pred_low,
                               xend=RT,
                               yend=RT_pred_high,
                               color = standardized_inchikey)) +
    facet_wrap(~model) +
    theme(legend.position = "none")
  
  print(scatter_plot)
  
  return(mc_data)
  
}

monte_carlo = function(data = tibble(),
                       p_value = 0.95) {
  for (sample_this in levels(as.factor(data$sample))) {
    data_this = data %>%
      filter(sample == sample_this)
    for (id_this in levels(as.factor(data_this$id))) {
      print(id_this)
      for (model_this in levels(as.factor(data_this$model))) {
        print(model_this)
        calibrant = data_this %>%
          filter((id != id_this) & correct & (model == model_this)) 
        candidates_this = data %>%
          filter((id == id_this) & (model == model_this))
        candidates_sum = tibble()
        for (i in 1:1000) {
          set.seed(i)
          print(i)
          calibrants_this = calibrant %>%
            group_by(standardized_SMILES) %>%
            mutate(sampled = sample.split(standardized_SMILES, SplitRatio = 1)) %>%
            ungroup() %>%
            filter(sampled) %>%
            select(-sampled)
          
          candidates_i = candidates_this %>%
            group_by(standardized_SMILES) %>%
            mutate(sampled = sample.split(standardized_SMILES, SplitRatio = 1)) %>%
            ungroup() %>%
            filter(sampled) %>%
            select(-sampled)
          
          gam_this = mgcv::gam(RT ~ s(RT_pred_mc, bs = "cr", k = 6),
                               data = calibrants_this)
          
          RT_pred_mc_i = mgcv::predict.gam(gam_this,
                                           newdata = candidates_i)
          
          RT_pred_mc_i = as_tibble(RT_pred_mc_i)
          
          candidates_i = candidates_i %>%
            bind_cols(RT_pred_mc_i)
          
          candidates_sum = candidates_sum %>%
            bind_rows(candidates_i %>%
                        mutate(i = i))
        }
        # candidates_sum = candidates_sum %>%
        #   group_by(standardized_SMILES) %>%
        #   summarize(RT_pred_low = quantile(RT_pred_mc_i, probs = (1-p_value)/2),
        #             RT_pred_high = quantile(RT_pred_mc_i, probs = p_value + (1-p_value)/2)) %>%
        #   ungroup()
        
        candidates_this = candidates_this %>%
          select(-delta_RT, -RT_pred_mc) %>%
          unique() %>%
          left_join(candidates_sum)
        
        write_delim(candidates_this,
                    paste("~/Library/CloudStorage/OneDrive-Kruvelab/Kruvelab/Research/artiklite kirjutamine/isomer_paper/code_2024-10-05/modelling/2025-02-28/data/sirius_pubchem_all/_", sample_this, id_this, model_this, ".csv", sep = " "))
        
      }
    }
  }
  
}