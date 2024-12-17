library(caret)
library(caTools)
library(cluster)
library(factoextra)
library(readxl)
library(tidyverse)
library(umap)


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
  
  model = lm(RT ~ .,
             data = training_data)
  
  return(model)
}


