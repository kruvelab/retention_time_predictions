
library(tidyverse)
library(rcdk)
library(caret)
library(caTools)
library(Metrics)
setwd("C:/Users/amso1873/OneDrive - Kruvelab/Amina/R code")

#importing the data 
dataset <-  read_delim('descs200502.csv',
                       delim = ",",
                       col_names = TRUE,
                       trim_ws = TRUE) %>%
  drop_na()

dataset2 <-  read_delim('Combined dataset.csv',
                        delim = ",",
                        col_names = TRUE,
                        trim_ws = TRUE) %>%
  drop_na()
SMILES_unique <-  read_delim('SMILES_unique.csv',
                             delim = ",",
                             col_names = TRUE,
                             trim_ws = TRUE) %>%
  drop_na()

dataset2 <- dataset2 %>% 
  left_join(dataset)
#Encoding categorical data
dataset2$Organic_modifier=factor(dataset2$Organic_modifier,
                                 levels = c('Acetonitrile','Methanol'),
                                 labels=c(1,2))
dataset2$Buffer=factor(dataset2$Buffer, levels = c('Formic acid', 'Ammonia', 'TFA', 'Acetic acid','Acetate', 'Formiate','Bicarbonate'),
                       labels=c(1,2,3,4,5,6,7))
dataset2$Column=factor(dataset2$Column, levels = c('C18-Mixed','C18-RP', 'HILIC','Biphenyl'),
                       labels=c(1,2,3,4))

# this will be removed, we will keep only M+H
dataset <- dataset2 
#Replacing the NAs with 0
dataset[is.na(dataset)] <- 0

dataset3 <- dataset %>%
  dplyr::select(-c(Column, Organic_modifier, Buffer,  bruto_formula, Compound_name, RT))


dataset3 <- dataset3 %>% 
  dplyr::select(-SMILES)
dataset3 <- dataset3 %>%
  dplyr::select(-nearZeroVar(dataset3)) 
#Calculating  the correlation between features
correlationMatrix <- cor(dataset3, use = "complete.obs")
# Finding the  highly corrected features (ideally >0.8)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.8)

dataset3 <- dataset3 %>%
  dplyr::select(-highlyCorrelated) %>%
  dplyr::mutate(Column= dataset$Column, Organic_modifier= dataset$Organic_modifier, Buffer=dataset$Buffer)

dataset3 <- dataset3 %>%
  dplyr::mutate(RT= dataset$RT, Compound_name=dataset$Compound_name, SMILES = dataset$SMILES)
dataset <- dataset3
#------ alpha calculations------------- 
dataset4 <-  read_delim('pka and logP for MM.csv',
                        delim = ",",
                        col_names = TRUE,
                        trim_ws = TRUE) %>%
  drop_na()
fitControl <- trainControl(method = "repeatedcv", number = 5, repeats = 5) 
lm_pred_logDan <- train(`logD A-` ~`logP AH`, 
                                  data = dataset4, 
                                  method = "lm", 
                                  trControl = fitControl)
print(lm_pred_logDan)
summary(lm_pred_logDan)
#R2=86,59

lm_pred_logDcat <- train(`logD AH2+` ~`logP AH`, 
                                   data = dataset4, 
                                   method = "lm", 
                                   trControl = fitControl)
print(lm_pred_logDcat)
summary(lm_pred_logDcat)
#R2=76,73


dataset4 <-  read_delim('pka and logP for MM.csv',
                        delim = ",",
                        col_names = TRUE,
                        trim_ws = TRUE)

logDA <- dataset4[4]
for(i in 1:78) { if (is.na(logDA[i,1])) {
  logDA[i,1] <-predict(lm_pred_logDan,dataset4[i,5])  }
}
logDAH2 <- dataset4[6]
for(i in 1:78) { if (is.na(logDAH2[i,1])) {
  logDAH2[i,1] <-predict(lm_pred_logDcat,dataset4[i,5])  }
}
pkaca <- dataset4[3]
pkaca[is.na(pkaca)] <--5
pkaan <- dataset4[2]
pkaan[is.na(pkaan)] <-20
logP <- dataset4[5]

dataset4 <- dataset4[1]
dataset4 <- data.frame(dataset4, pkaca, pkaan, logDA, logP, logDAH2 )

dataset <- dataset %>%
  left_join(dataset4)

dataset <- data.frame(dataset, alphaAH2= 1/(10^(dataset$pH - dataset$pKa.base)+1))
dataset <- data.frame(dataset, alphaA_H= (10^(dataset$pH-dataset$pKa.acid))/(10^(dataset$pH-dataset$pKa.acid )+1))
dataset <- data.frame(dataset, alphaAH= 1-(dataset$alphaA_H +dataset$alphaAH2))

dataset <- data.frame(dataset, DA= 10^dataset$logD.A.)
dataset <- data.frame(dataset, P=10^dataset$logP.AH)
dataset <- data.frame(dataset, DAH2 = 10^dataset$logD.AH2. )

dataset <- data.frame(dataset, P_sum = dataset$alphaA_H*dataset$DA +dataset$alphaAH*dataset$P + dataset$alphaAH2*dataset$DAH2)
dataset <- data.frame(dataset, log_P_sum = log10(dataset$P_sum))
n <- as.numeric(which(colnames(dataset)=="log_P_sum"))
y <- as.numeric(which(colnames(dataset)=="logP.AH"))
for(i in 1:2961) { if (is.nan(dataset[i,n])) {
  dataset[i,n] <-dataset[i,y]  }
}
dataset <- dataset %>% dplyr::select(-c(pKa.base, pKa.acid, logD.A., logP.AH, logD.AH2., DA, P, DAH2, P_sum, alphaAH2, alphaA_H, alphaAH))
Combined <- dataset
#------------- MM model -----------------
dataset <- Combined %>% filter(Column==1)
dataset <- dataset %>% dplyr::select(-Column)
dataset <- dataset %>%   left_join(SMILES_unique) %>%  
  dplyr::select(SPLIT, everything()) 
for(i in 1:1153) { if (is.na(dataset[i,1])) {
  dataset[i,1] <-TRUE  }
}
#trainingset 
Mixed_training_set <- dataset %>%   filter(SPLIT == TRUE) %>%  
  dplyr::select(-SPLIT) 
#940 datapoints
#testset 
Mixed_test_set <- dataset %>%   filter(SPLIT == FALSE) %>%
  dplyr::select(-SPLIT) 
#213 datapoints 

#-------- running all the models --------
set.seed(100)
folds <- groupKFold(Mixed_training_set$SMILES, k = 20)

fitControl <- trainControl(## 10-fold CV
  method = "boot",
  index = folds)
#-----------------------Training and cross-validation -----------------------------
#-------------------------Independent Component Regression----------------------
#install.packages("fastICA")
library(fastICA)
set.seed(825)
regressor_icr <- train(RT ~ ., 
                       data = Mixed_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "icr", 
                       trControl = fitControl,
                       verbose=FALSE)

print(regressor_icr)
summary(regressor_icr)
regressor_icr$finalModel

#-------------------------Linear Regression with Backwards Selection----------------------
#install.packages("leaps")
library(leaps)
set.seed(825)
regressor_LRBS <- train(RT ~ ., 
                        data = Mixed_training_set%>% dplyr::select(-c(SMILES, Compound_name)), 
                        method = "leapBackward", 
                        trControl = fitControl,
                        verbose = FALSE)

print(regressor_LRBS)
summary(regressor_LRBS)
regressor_LRBS$finalModel
#-------------------------Partial least squares ----------------------
#install.packages("pls")
library(pls)
set.seed(825)
regressor_pls <- train(RT ~ ., 
                       data = Mixed_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "pls", 
                       trControl = fitControl,
                       verbose = FALSE)

print(regressor_pls)
summary(regressor_pls)
regressor_pls$finalModel

#-------------------------Multivariate Adaptive Regression Spline-------------------------
#install.packages("earth")
library(earth)
set.seed(825)
regressor_MARS <- train(RT ~ ., 
                        data = Mixed_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                        method = "earth", 
                        trControl = fitControl)

print(regressor_MARS)
summary(regressor_MARS)

#-------------------------KNN-------------------------
#install.packages("kknn")
library(kknn)
set.seed(825)
regressor_knn <- train(RT ~ ., 
                       data = Mixed_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "kknn", 
                       trControl = fitControl)

print(regressor_knn)
summary(regressor_knn)
#-------------------------RF-------------------------
set.seed(825)
regressor_RF <- train(RT ~ ., 
                      data = Mixed_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                      method = "ranger", 
                      trControl = fitControl, importance='impurity')

print(regressor_RF)
summary(regressor_RF)
varImp(regressor_RF)
regressor_RF$finalModel
saveRDS(regressor_RF,  "regressorRF for mixed mode 16082021.rds")
#-------------------------Validation--------------------------------
#Predicting for the test set
pred_icr_test = predict(regressor_icr, newdata = Mixed_test_set)
pred_LRBS_test = predict(regressor_LRBS, newdata = Mixed_test_set)
pred_pls_test = predict(regressor_pls, newdata = Mixed_test_set)
pred_MARS_test = predict(regressor_MARS, newdata = Mixed_test_set) 
pred_knn_test = predict(regressor_knn, newdata = Mixed_test_set) 
pred_RF_test = predict(regressor_RF, newdata = Mixed_test_set) 

#Predicting tR for the training set 
pred_icr_training = predict(regressor_icr, newdata = Mixed_training_set)
pred_LRBS_training = predict(regressor_LRBS, newdata = Mixed_training_set)
pred_pls_training = predict(regressor_pls, newdata = Mixed_training_set)
pred_MARS_training= predict(regressor_MARS, newdata = Mixed_training_set) 
pred_knn_training= predict(regressor_knn, newdata = Mixed_training_set) 
pred_RF_training = predict(regressor_RF, newdata = Mixed_training_set) 

Mixed_test_set<- Mixed_test_set %>%
  mutate(pred_icr_test, pred_LRBS_test, pred_pls_test,  pred_MARS_test,  pred_knn_test, pred_RF_test)

Mixed_training_set<- Mixed_training_set %>%
  mutate(pred_icr_training, pred_LRBS_training, pred_pls_training,  pred_MARS_training,  pred_knn_training, pred_RF_training)


#Calculation the relative error 
#Predicting tR for the training set 
ICR <- c()
LRBS <- c()
pls <- c()
MARS <- c()
knn <- c()
RF <- c()
for( i in 1:213) {
  ICR[i] = Mixed_test_set$RT[i]-Mixed_test_set$pred_icr_test[i]
  LRBS[i] = Mixed_test_set$RT[i]-Mixed_test_set$pred_LRBS_test[i]
  pls[i] = Mixed_test_set$RT[i]-Mixed_test_set$pred_pls_test[i]
  MARS[i] = Mixed_test_set$RT[i]-Mixed_test_set$pred_MARS_test[i]
  knn[i] = Mixed_test_set$RT[i]-Mixed_test_set$pred_knn_test[i]
  RF[i] =  Mixed_test_set$RT[i]-Mixed_test_set$pred_RF_test[i]
}
Mixed_test_set<- Mixed_test_set %>%
  mutate(ICR , LRBS, pls ,MARS ,knn, RF )

for( i in 1:940) {
  ICR[i] = Mixed_training_set$RT[i]-Mixed_training_set$pred_icr_training[i]
  LRBS[i] = Mixed_training_set$RT[i]-Mixed_training_set$pred_LRBS_training[i]
  pls[i] =  Mixed_training_set$RT[i]-Mixed_training_set$pred_pls_training[i]
  MARS[i] = Mixed_training_set$RT[i]-Mixed_training_set$pred_MARS_training[i]
  knn[i] = Mixed_training_set$RT[i]-Mixed_training_set$pred_knn_training[i]
  RF[i] = Mixed_training_set$RT[i]-Mixed_training_set$pred_RF_training[i]
}
Mixed_training_set<- Mixed_training_set %>%
  mutate(ICR , LRBS, pls ,MARS ,knn, RF )


#calculating the RMSE for the test set 
RSME_icr_testset <- RMSE(Mixed_test_set$pred_icr_test, Mixed_test_set$RT)
RSME_LRBS_testset <- RMSE(Mixed_test_set$pred_LRBS_test, Mixed_test_set$RT)
RSME_pls_testset <- RMSE(Mixed_test_set$pred_pls_test, Mixed_test_set$RT)
RSME_MARS_testset <- RMSE(Mixed_test_set$pred_MARS_test, Mixed_test_set$RT)
RSME_knn_testset <- RMSE(Mixed_test_set$pred_knn_test, Mixed_test_set$RT)
RSME_RF_testset <- RMSE(Mixed_test_set$pred_RF_test, Mixed_test_set$RT)

RSME_testset <- c(RSME_icr_testset, RSME_LRBS_testset, RSME_pls_testset, RSME_MARS_testset,  RSME_knn_testset, RSME_RF_testset)

#calculating the RMSE for the training set 
RSME_icr_trainingset <- RMSE(Mixed_training_set$pred_icr_training, Mixed_training_set$RT)
RSME_LRBS_trainingset <- RMSE(Mixed_training_set$pred_LRBS_training, Mixed_training_set$RT)
RSME_pls_trainingset<- RMSE(Mixed_training_set$pred_pls_training, Mixed_training_set$RT)
RSME_MARS_trainingset <- RMSE(Mixed_training_set$pred_MARS_training, Mixed_training_set$RT)
RSME_knn_trainingset <- RMSE(Mixed_training_set$pred_knn_training, Mixed_training_set$RT)
RSME_RF_trainingset <- RMSE(Mixed_training_set$pred_RF_training, Mixed_training_set$RT)

RSME_trainingset <- c(RSME_icr_trainingset, RSME_LRBS_trainingset, RSME_pls_trainingset,  RSME_MARS_trainingset,  RSME_knn_trainingset, RSME_RF_trainingset)

Models <- c("ICR", "LRBS", "PLS",  "´MARS",  "KNN", "RF")
library(extrafont)
font_import()
font <- choose_font("Arial")
fontsize <- 12
basecolor <- "black"
my_theme <- theme(plot.background = element_blank(), 
                  panel.background = element_blank(), axis.line=element_line(  size=0.5), 
                  legend.background = element_blank(),
                  legend.title = element_text(),
                  legend.position = c(0.9,0.3), 
                  aspect.ratio = 1, text = element_text(family = font, size=fontsize),
                  plot.title = element_text(color=basecolor, size=14, face="bold"),
                  legend.text=element_text(family = font, size=fontsize, color=basecolor),
                  axis.text=element_text(family = font, size=fontsize, color=basecolor), 
                  panel.grid.major =  element_blank(),
                  panel.grid.minor=element_blank(), 
                  panel.border = element_blank())

p1 <-ggplot() +
  geom_point(mapping = aes(x = Models, y =RSME_trainingset), color =  "#3d5a80", size = 2, alpha = 1) +
  geom_point(mapping = aes(x = Models, y =RSME_testset),color ="#ee6c4d",  size = 2, alpha = 1) +
  ylab("RMSE (min)")+
  theme_bw()

p1+my_theme
 
p2 <- ggplot() +
  geom_point(mapping = aes(x =Mixed_training_set$RT, y =Mixed_training_set$pred_RF_training), color =  "#3d5a80", size = 2, alpha = 1) +
  geom_point(mapping = aes(x = Mixed_test_set$RT, y =Mixed_test_set$pred_RF_test), color ="#ee6c4d",  size = 2, alpha = 1) +
  xlab("Measured"~italic("t")["R"] ~"(min)")+
  ylab("Predicted"~ italic("t")["R"] ~"(min)")+
  theme_bw()
p2+my_theme

print(RSME_RF_testset)
#1.79
print(RSME_RF_trainingset)
#0,51

Mixed_RF_training <- Mixed_training_set %>% dplyr::select(-c(ICR,LRBS, pls, MARS, knn))
names(Mixed_RF_training)[names(Mixed_RF_training) == "RF"] <- "Error"
Mixed_RF_training <- Mixed_RF_training %>% mutate(Model="RF")
ICR <- Mixed_training_set %>% dplyr::select(-c(RF,LRBS, pls, MARS, knn))
names(ICR)[names(ICR) == "ICR"] <- "Error"
ICR <- ICR %>% mutate(Model="ICR")
LRBS <- Mixed_training_set %>% dplyr::select(-c(RF,ICR, pls, MARS, knn))
names(LRBS)[names(LRBS) == "LRBS"] <- "Error"
LRBS <- LRBS %>% mutate(Model="LRBS")
pls <- Mixed_training_set %>% dplyr::select(-c(RF,ICR,LRBS, MARS, knn))
names(pls)[names(pls) == "pls"] <- "Error"
pls <- pls %>% mutate(Model="PLS")
MARS <- Mixed_training_set %>% dplyr::select(-c(RF,ICR,LRBS, pls, knn))
names(MARS)[names(MARS) == "MARS"] <- "Error"
MARS <- MARS %>% mutate(Model="MARS")
knn <- Mixed_training_set %>% dplyr::select(-c(RF,ICR,LRBS, pls,MARS))
names(knn)[names(knn) == "knn"] <- "Error"
knn <- knn %>% mutate(Model="KNN")

Mixed_RF_training <- Mixed_RF_training %>%
  bind_rows(ICR) %>% bind_rows(LRBS) %>% 
  bind_rows(pls) %>% bind_rows(MARS) %>% bind_rows(knn)

p3 <- ggplot(data = Mixed_RF_training,
              aes(x = Model, y = Error)) +
  geom_boxplot()+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               binwidth = 0.04,
               #fill = "grey",
               alpha = 0.5
  )+ ylab("Error (min)")+xlab("Models")
p3+my_theme

Mixed_RF_test <- Mixed_test_set %>% dplyr::select(-c(ICR,LRBS, pls, MARS, knn))
names(Mixed_RF_test)[names(Mixed_RF_test) == "RF"] <- "Error"
Mixed_RF_test <- Mixed_RF_test %>% mutate(Model="RF")
ICR <- Mixed_test_set %>% dplyr::select(-c(RF,LRBS, pls, MARS, knn))
names(ICR)[names(ICR) == "ICR"] <- "Error"
ICR <- ICR %>% mutate(Model="ICR")
LRBS <- Mixed_test_set %>% dplyr::select(-c(RF,ICR, pls, MARS, knn))
names(LRBS)[names(LRBS) == "LRBS"] <- "Error"
LRBS <- LRBS %>% mutate(Model="LRBS")
pls <- Mixed_test_set %>% dplyr::select(-c(RF,ICR,LRBS, MARS, knn))
names(pls)[names(pls) == "pls"] <- "Error"
pls <- pls %>% mutate(Model="PLS")
MARS <- Mixed_test_set %>% dplyr::select(-c(RF,ICR,LRBS, pls, knn))
names(MARS)[names(MARS) == "MARS"] <- "Error"
MARS <- MARS %>% mutate(Model="MARS")
knn <- Mixed_test_set %>% dplyr::select(-c(RF,ICR,LRBS, pls,MARS))
names(knn)[names(knn) == "knn"] <- "Error"
knn <- knn %>% mutate(Model="KNN")

Mixed_RF_test<- Mixed_RF_test %>%
  bind_rows(ICR) %>% bind_rows(LRBS) %>% 
  bind_rows(pls) %>% bind_rows(MARS) %>% bind_rows(knn)

p4<- ggplot(data = Mixed_RF_test,
             aes(x = Model, y = Error)) +
  geom_boxplot()+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               binwidth = 0.04, #give it different values
               #fill = "grey",
               alpha = 0.5
  )+ ylim(-10,10) +ylab("Error (min)")+xlab("Models")
p4+my_theme


library(readxl)
varImp(regressor_RF)
Important_factors_in_RF_mixed_mode <- read_excel("Important factors in RF.xlsx", 
                                      sheet = "Mixed mode ")

p5 <- ggplot(data=Important_factors_in_RF_mixed_mode)+
  geom_col(mapping=aes(x=reorder(Feature, Importance), y= Importance))+
  xlab("Variables")+
  ylab("Importance in the random forest model")+
  theme_bw()

p5 + my_theme
#-------------------------LM-------------------------
set.seed(825)
regressor_lm <- train(RT ~ log_P_sum,
                      data = Mixed_training_set,
                      method = "lm",
                      trControl = fitControl)

print(regressor_lm)
summary(regressor_lm)

dataset <- dataset %>%
  mutate(tR_pred = predict(regressor_lm, newdata = dataset))
dataset <- dataset %>% left_join(SMILES_unique)
for (i in 1:1153 ) {
  if (dataset$SPLIT[i]=="TRUE") {dataset$SPLIT[i]<-"training set"} 
  else dataset$SPLIT[i]<-"test set"}
training_set <- dataset %>% filter(SPLIT=="training set")
test_set <- dataset %>% filter(SPLIT=="test set")
RSME_lm_trainingset <- RMSE( training_set$tR_pred, training_set$RT)
#2.59min
RSME_lm_testset <- RMSE(test_set$tR_pred, test_set$RT)
#2.36
p6 <- ggplot(data=dataset)+
  geom_point(mapping=aes(x=RT, y= tR_pred, color = SPLIT),size=3, alpha=1)+
  scale_color_manual(values=c("#ee6c4d","#3d5a80")) +
  xlab("Measured" ~italic("t")["R"]~ "(min)")+
  ylab("Predicted " ~italic("t")["R"]~ "(min)")+
  labs("Predicted VS Actual RT")
p6+my_theme
F <- (RSME_lm_testset*RSME_lm_testset)/(RSME_RF_testset*RSME_RF_testset)
pf(F, 211, 211, lower.tail = FALSE)*100
# 0.0036272%

#---------------RP model------
dataset <- Combined %>% filter(Column==2)
dataset <- dataset %>% dplyr::select(-Column)
dataset <- dataset %>%   left_join(SMILES_unique) %>%  
  dplyr::select(SPLIT, everything()) 
for(i in 1:839) { if (is.na(dataset[i,1])) {
  dataset[i,1] <-TRUE  }
}
#trainingset 
RP_training_set <- dataset %>%   filter(SPLIT == TRUE) %>%  
  dplyr::select(-SPLIT)
#697 datapoints
#testset 
RP_test_set <- dataset %>%   filter(SPLIT == FALSE) %>%
  dplyr::select(-SPLIT) 
#142 datapoints
set.seed(100)
folds <- groupKFold(RP_training_set$SMILES, k = 20)
 
fitControl <- trainControl(## 10-fold CV
  method = "boot",
  index = folds)

#-------------------------LM-------------------------
set.seed(825)
regressor_lm <- train(RT ~ log_P_sum, 
                      data = RP_training_set, 
                      method = "lm", 
                      trControl = fitControl)

print(regressor_lm)
summary(regressor_lm)
dataset <- dataset %>% 
  mutate(tR_pred = predict(regressor_lm, newdata = dataset))

for (i in 1:839 ) {
  if (dataset$SPLIT[i]=="TRUE") {dataset$SPLIT[i]<-"training set"} 
  else dataset$SPLIT[i]<-"test set"}
training_set <- dataset %>% filter(SPLIT=="training set")
test_set <- dataset %>% filter(SPLIT=="test set")
RSME_lm_trainingset <- RMSE( training_set$tR_pred, training_set$RT)
#2.48min
RSME_lm_testset <- RMSE(test_set$tR_pred, test_set$RT)
#2.13
p7 <- ggplot(data=dataset)+
  geom_point(mapping=aes(x=RT, y= tR_pred, color=SPLIT), size=3, alpha=1, show.legend = FALSE)+
  scale_color_manual(values=c("#ee6c4d","#3d5a80")) +
  xlab("Measured"~ italic("t")["R"]~ "(min)")+
  ylab("Predicted"~ italic("t")["R"]~ "(min)")
p7+my_theme

#-----------------------Training and cross-validation -----------------------------
#install.packages("fastICA")
library(fastICA)
set.seed(825)
regressor_lm <- train(RT ~ ., 
                      data = RP_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                      method = "lm", 
                      trControl = fitControl,
                      verbose=FALSE)

print(regressor_lm)
summary(regressor_lm)
regressor_lm$finalModel
#-------------------------Independent Component Regression----------------------
#install.packages("fastICA")
library(fastICA)
set.seed(825)
regressor_icr <- train(RT ~ ., 
                       data = RP_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "icr", 
                       trControl = fitControl,
                       verbose=FALSE)

print(regressor_icr)
summary(regressor_icr)
regressor_icr$finalModel

#-------------------------Linear Regression with Backwards Selection----------------------
#install.packages("leaps")
library(leaps)
set.seed(825)
regressor_LRBS <- train(RT ~ ., 
                        data = RP_training_set%>% dplyr::select(-c(SMILES, Compound_name)), 
                        method = "leapBackward", 
                        trControl = fitControl,
                        verbose = FALSE)

print(regressor_LRBS)
summary(regressor_LRBS)
regressor_LRBS$finalModel
#-------------------------Partial least squares ----------------------
#install.packages("pls")
library(pls)
set.seed(825)
regressor_pls <- train(RT ~ ., 
                       data = RP_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "pls", 
                       trControl = fitControl,
                       verbose = FALSE)

print(regressor_pls)
summary(regressor_pls)
regressor_pls$finalModel


#-------------------------Multivariate Adaptive Regression Spline-------------------------
#install.packages("earth")
library(earth)
set.seed(825)
regressor_MARS <- train(RT ~ ., 
                        data = RP_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                        method = "earth", 
                        trControl = fitControl)

print(regressor_MARS)
summary(regressor_MARS)

#-------------------------KNN-------------------------
#install.packages("kknn")
library(kknn)
set.seed(825)
regressor_knn <- train(RT ~ ., 
                       data = RP_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "kknn", 
                       trControl = fitControl)

print(regressor_knn)
summary(regressor_knn)
#-------------------------RF-------------------------
set.seed(825)
regressor_RF <- train(RT ~ ., 
                      data = RP_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                      method = "ranger", 
                      trControl = fitControl, importance='impurity')

print(regressor_RF)
summary(regressor_RF)
varImp(regressor_RF) 
saveRDS(regressor_RF,  "regressorRF for RP 16082021.rds")


#-------------------------Validation--------------------------------
#Predicting for the test set
pred_icr_test = predict(regressor_icr, newdata = RP_test_set)
pred_LRBS_test = predict(regressor_LRBS, newdata = RP_test_set)
pred_pls_test = predict(regressor_pls, newdata = RP_test_set)
pred_MARS_test = predict(regressor_MARS, newdata = RP_test_set) 
pred_knn_test = predict(regressor_knn, newdata = RP_test_set) 
pred_RF_test = predict(regressor_RF, newdata = RP_test_set) 

#Predicting tR for the training set 
pred_icr_training = predict(regressor_icr, newdata = RP_training_set)
pred_LRBS_training = predict(regressor_LRBS, newdata = RP_training_set)
pred_pls_training = predict(regressor_pls, newdata = RP_training_set)
pred_MARS_training= predict(regressor_MARS, newdata = RP_training_set) 
pred_knn_training= predict(regressor_knn, newdata = RP_training_set) 
pred_RF_training = predict(regressor_RF, newdata = RP_training_set) 

RP_test_set<- RP_test_set %>%
  mutate(pred_icr_test, pred_LRBS_test, pred_pls_test,  pred_MARS_test,  pred_knn_test, pred_RF_test)

RP_training_set<- RP_training_set %>%
  mutate(pred_icr_training, pred_LRBS_training, pred_pls_training,  pred_MARS_training,  pred_knn_training, pred_RF_training)

#calculating the RMSE for the test set 
RSME_icr_testset <- RMSE(RP_test_set$pred_icr_test, RP_test_set$RT)
RSME_LRBS_testset <- RMSE(RP_test_set$pred_LRBS_test, RP_test_set$RT)
RSME_pls_testset <- RMSE(RP_test_set$pred_pls_test, RP_test_set$RT)
RSME_MARS_testset <- RMSE(RP_test_set$pred_MARS_test, RP_test_set$RT)
RSME_knn_testset <- RMSE(RP_test_set$pred_knn_test, RP_test_set$RT)
RSME_RF_testset <- RMSE(RP_test_set$pred_RF_test, RP_test_set$RT)

RSME_testset <- c(RSME_icr_testset, RSME_LRBS_testset, RSME_pls_testset, RSME_MARS_testset,  RSME_knn_testset, RSME_RF_testset)

#calculating the RMSE for the training set 
RSME_icr_trainingset <- RMSE(RP_training_set$pred_icr_training, RP_training_set$RT)
RSME_LRBS_trainingset <- RMSE(RP_training_set$pred_LRBS_training, RP_training_set$RT)
RSME_pls_trainingset<- RMSE(RP_training_set$pred_pls_training, RP_training_set$RT)
RSME_MARS_trainingset <- RMSE(RP_training_set$pred_MARS_training, RP_training_set$RT)
RSME_knn_trainingset <- RMSE(RP_training_set$pred_knn_training, RP_training_set$RT)
RSME_RF_trainingset <- RMSE(RP_training_set$pred_RF_training, RP_training_set$RT)

RSME_trainingset <- c(RSME_icr_trainingset, RSME_LRBS_trainingset, RSME_pls_trainingset,  RSME_MARS_trainingset,  RSME_knn_trainingset, RSME_RF_trainingset)

Models <- c("ICR", "LRBS", "PLS",  "´MARS",  "KNN", "RF")


p8<- ggplot() +
  geom_point(mapping = aes(x = Models, y =RSME_trainingset), color = "#3d5a80", size = 2, alpha = 1) +
  geom_point(mapping = aes(x = Models, y =RSME_testset),color = "#ee6c4d",  size = 2, alpha = 1) +
  xlab("Models")+
  ylab("RMSE (min)")+
  theme_bw()
p8+my_theme

p9 <- ggplot() +
  geom_point(mapping = aes(x = RP_training_set$RT, y =RP_training_set$pred_RF_training), color = "#3d5a80", size = 2, alpha = 1) +
  geom_point(mapping = aes(x = RP_test_set$RT, y =RP_test_set$pred_RF_test), color = "#ee6c4d", size = 2, alpha = 1) +
  xlab("Measured " ~ italic("t")["R"]~" (min)")+
  ylab("Predicted" ~italic("t")["R"]~" (min)")+
  theme_bw()
p9 + my_theme

print(RSME_RF_testset)
#1,53
print(RSME_RF_trainingset)
#0,51


F <- (RSME_lm_testset*RSME_lm_testset)/(RSME_RF_testset*RSME_RF_testset)
pf(F, 141,141, lower.tail = FALSE)*100
# 0.005822968%

#Calculation the relative error 
#Predicting tR for the training set 
ICR <- c()
LRBS <- c()
pls <- c()
MARS <- c()
knn <- c()
RF <- c()
for( i in 1:142) {
  ICR[i] = RP_test_set$RT[i]-RP_test_set$pred_icr_test[i]
  LRBS[i] = RP_test_set$RT[i]-RP_test_set$pred_LRBS_test[i]
  pls[i] = RP_test_set$RT[i]-RP_test_set$pred_pls_test[i]
  MARS[i] = RP_test_set$RT[i]-RP_test_set$pred_MARS_test[i]
  knn[i] = RP_test_set$RT[i]-RP_test_set$pred_knn_test[i]
  RF[i] =  RP_test_set$RT[i]-RP_test_set$pred_RF_test[i]
}
RP_test_set<- RP_test_set %>%
  mutate(ICR , LRBS, pls ,MARS ,knn, RF )

for( i in 1:697) {
  ICR[i] = RP_training_set$RT[i]-RP_training_set$pred_icr_training[i]
  LRBS[i] = RP_training_set$RT[i]-RP_training_set$pred_LRBS_training[i]
  pls[i] =  RP_training_set$RT[i]-RP_training_set$pred_pls_training[i]
  MARS[i] = RP_training_set$RT[i]-RP_training_set$pred_MARS_training[i]
  knn[i] = RP_training_set$RT[i]-RP_training_set$pred_knn_training[i]
  RF[i] = RP_training_set$RT[i]-RP_training_set$pred_RF_training[i]
}
RP_training_set<- RP_training_set %>%
  mutate(ICR , LRBS, pls ,MARS ,knn, RF )

RP_RF_training <- RP_training_set %>% dplyr::select(-c(ICR,LRBS, pls, MARS, knn))
names(RP_RF_training)[names(RP_RF_training) == "RF"] <- "Error"
RP_RF_training <- RP_RF_training %>% mutate(Model="RF")
ICR <- RP_training_set %>% dplyr::select(-c(RF,LRBS, pls, MARS, knn))
names(ICR)[names(ICR) == "ICR"] <- "Error"
ICR <- ICR %>% mutate(Model="ICR")
LRBS <- RP_training_set %>% dplyr::select(-c(RF,ICR, pls, MARS, knn))
names(LRBS)[names(LRBS) == "LRBS"] <- "Error"
LRBS <- LRBS %>% mutate(Model="LRBS")
pls <- RP_training_set %>% dplyr::select(-c(RF,ICR,LRBS, MARS, knn))
names(pls)[names(pls) == "pls"] <- "Error"
pls <- pls %>% mutate(Model="PLS")
MARS <- RP_training_set %>% dplyr::select(-c(RF,ICR,LRBS, pls, knn))
names(MARS)[names(MARS) == "MARS"] <- "Error"
MARS <- MARS %>% mutate(Model="MARS")
knn <- RP_training_set %>% dplyr::select(-c(RF,ICR,LRBS, pls,MARS))
names(knn)[names(knn) == "knn"] <- "Error"
knn <- knn %>% mutate(Model="KNN")

RP_RF_training <- RP_RF_training %>%
  bind_rows(ICR) %>% bind_rows(LRBS) %>% 
  bind_rows(pls) %>% bind_rows(MARS) %>% bind_rows(knn)

p10 <- ggplot(data = RP_RF_training,
             aes(x = Model, y = Error)) +
  geom_boxplot()+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               binwidth = 0.04,
               #fill = "grey",
               alpha = 0.5
  )+ ylab("Error (min)")+xlab("Models")
p10 +my_theme
RP_RF_test <- RP_test_set %>% dplyr::select(-c(ICR,LRBS, pls, MARS, knn))
names(RP_RF_test)[names(RP_RF_test) == "RF"] <- "Error"
RP_RF_test <- RP_RF_test %>% mutate(Model="RF")
ICR <- RP_test_set %>% dplyr::select(-c(RF,LRBS, pls, MARS, knn))
names(ICR)[names(ICR) == "ICR"] <- "Error"
ICR <- ICR %>% mutate(Model="ICR")
LRBS <- RP_test_set %>% dplyr::select(-c(RF,ICR, pls, MARS, knn))
names(LRBS)[names(LRBS) == "LRBS"] <- "Error"
LRBS <- LRBS %>% mutate(Model="LRBS")
pls <- RP_test_set %>% dplyr::select(-c(RF,ICR,LRBS, MARS, knn))
names(pls)[names(pls) == "pls"] <- "Error"
pls <- pls %>% mutate(Model="PLS")
MARS <- RP_test_set %>% dplyr::select(-c(RF,ICR,LRBS, pls, knn))
names(MARS)[names(MARS) == "MARS"] <- "Error"
MARS <- MARS %>% mutate(Model="MARS")
knn <- RP_test_set %>% dplyr::select(-c(RF,ICR,LRBS, pls,MARS))
names(knn)[names(knn) == "knn"] <- "Error"
knn <- knn %>% mutate(Model="KNN")

RP_RF_test<- RP_RF_test %>%
  bind_rows(ICR) %>% bind_rows(LRBS) %>% 
  bind_rows(pls) %>% bind_rows(MARS) %>% bind_rows(knn)

p11 <- ggplot(data = RP_RF_test,
             aes(x = Model, y = Error)) +
  geom_boxplot()+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               binwidth = 0.04, #give it different values
               #fill = "grey",
               alpha = 0.5
  )+ ylim(-10,10) +ylab("Error (min)")+xlab("Models")
p11+my_theme
varImp(regressor_RF)
Important_factors_in_RF_RP <- read_excel("Important factors in RF.xlsx", 
                                                 sheet = "RP") 

p12 <- ggplot(data=Important_factors_in_RF_RP)+
  geom_col(mapping=aes(x=reorder(Feature, Importance), y= Importance))+
  xlab("Variables")+
  ylab("Importance in the random forest model")+
  theme_bw()

p12 + my_theme

#---------------HILIC model------
dataset <- Combined %>% filter(Column==3)
dataset <- dataset %>% dplyr::select(-Column)
dataset <- dataset %>%   left_join(SMILES_unique) %>%  
  dplyr::select(SPLIT, everything()) 
dataset <-  dataset %>%
  mutate(NH4 = case_when(
    pH == 3.0 ~ 1,
    pH == 5.0 ~ 1,
    pH == 6.7 ~1,
    TRUE ~ 0))

for(i in 1:836) { if (is.na(dataset[i,1])) {
  dataset[i,1] <-TRUE  }
}
#trainingset 
HILIC_training_set <- dataset %>%   filter(SPLIT == TRUE) %>%  
  dplyr::select(-SPLIT)
#685 datapoints
#testset 
HILIC_test_set <- dataset %>%   filter(SPLIT == FALSE) %>%
  dplyr::select(-SPLIT) 
#151 datapoints
set.seed(100)
folds <- groupKFold(HILIC_training_set$SMILES, k = 20)

fitControl <- trainControl(## 10-fold CV
  method = "boot",
  index = folds)

#-----------------------Training and cross-validation -----------------------------
#install.packages("fastICA")
library(fastICA)
set.seed(825)
regressor_lm <- train(RT ~ ., 
                      data = HILIC_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                      method = "lm", 
                      trControl = fitControl,
                      verbose=FALSE)

print(regressor_lm)
summary(regressor_lm)
regressor_lm$finalModel
#-------------------------Independent Component Regression----------------------
#install.packages("fastICA")
library(fastICA)
set.seed(825)
regressor_icr <- train(RT ~ ., 
                       data = HILIC_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "icr", 
                       trControl = fitControl,
                       verbose=FALSE)

print(regressor_icr)
summary(regressor_icr)
regressor_icr$finalModel

#-------------------------Linear Regression with Backwards Selection----------------------
#install.packages("leaps")
library(leaps)
set.seed(825)
regressor_LRBS <- train(RT ~ ., 
                        data = HILIC_training_set%>% dplyr::select(-c(SMILES, Compound_name)), 
                        method = "leapBackward", 
                        trControl = fitControl,
                        verbose = FALSE)

print(regressor_LRBS)
summary(regressor_LRBS)
regressor_LRBS$finalModel
#-------------------------Partial least squares ----------------------
#install.packages("pls")
library(pls)
set.seed(825)
regressor_pls <- train(RT ~ ., 
                       data = HILIC_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "pls", 
                       trControl = fitControl,
                       verbose = FALSE)

print(regressor_pls)
summary(regressor_pls)
regressor_pls$finalModel


#-------------------------Multivariate Adaptive Regression Spline-------------------------
#install.packages("earth")
library(earth)
set.seed(825)
regressor_MARS <- train(RT ~ ., 
                        data = HILIC_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                        method = "earth", 
                        trControl = fitControl)

print(regressor_MARS)
summary(regressor_MARS)

#-------------------------KNN-------------------------
#install.packages("kknn")
library(kknn)
set.seed(825)
regressor_knn <- train(RT ~ ., 
                       data = HILIC_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "kknn", 
                       trControl = fitControl)

print(regressor_knn)
summary(regressor_knn)
#-------------------------RF-------------------------
set.seed(825)
regressor_RF <- train(RT ~ ., 
                      data = HILIC_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                      method = "ranger", 
                      trControl = fitControl, importance='impurity')

print(regressor_RF)
summary(regressor_RF)
varImp(regressor_RF) 
saveRDS(regressor_RF,  "regressorRF for HILIC 16082021.rds")


#-------------------------Validation--------------------------------
#Predicting for the test set
pred_icr_test = predict(regressor_icr, newdata = HILIC_test_set)
pred_LRBS_test = predict(regressor_LRBS, newdata = HILIC_test_set)
pred_pls_test = predict(regressor_pls, newdata = HILIC_test_set)
pred_MARS_test = predict(regressor_MARS, newdata = HILIC_test_set) 
pred_knn_test = predict(regressor_knn, newdata = HILIC_test_set) 
pred_RF_test = predict(regressor_RF, newdata = HILIC_test_set) 

#Predicting tR for the training set 
pred_icr_training = predict(regressor_icr, newdata = HILIC_training_set)
pred_LRBS_training = predict(regressor_LRBS, newdata = HILIC_training_set)
pred_pls_training = predict(regressor_pls, newdata = HILIC_training_set)
pred_MARS_training= predict(regressor_MARS, newdata = HILIC_training_set) 
pred_knn_training= predict(regressor_knn, newdata = HILIC_training_set) 
pred_RF_training = predict(regressor_RF, newdata = HILIC_training_set) 

HILIC_test_set<- HILIC_test_set %>%
  mutate(pred_icr_test, pred_LRBS_test, pred_pls_test,  pred_MARS_test,  pred_knn_test, pred_RF_test)

HILIC_training_set<- HILIC_training_set %>%
  mutate(pred_icr_training, pred_LRBS_training, pred_pls_training,  pred_MARS_training,  pred_knn_training, pred_RF_training)

#calculating the RMSE for the test set 
RSME_icr_testset <- RMSE(HILIC_test_set$pred_icr_test, HILIC_test_set$RT)
RSME_LRBS_testset <- RMSE(HILIC_test_set$pred_LRBS_test, HILIC_test_set$RT)
RSME_pls_testset <- RMSE(HILIC_test_set$pred_pls_test, HILIC_test_set$RT)
RSME_MARS_testset <- RMSE(HILIC_test_set$pred_MARS_test, HILIC_test_set$RT)
RSME_knn_testset <- RMSE(HILIC_test_set$pred_knn_test, HILIC_test_set$RT)
RSME_RF_testset <- RMSE(HILIC_test_set$pred_RF_test, HILIC_test_set$RT)

RSME_testset <- c(RSME_icr_testset, RSME_LRBS_testset, RSME_pls_testset, RSME_MARS_testset,  RSME_knn_testset, RSME_RF_testset)

#calculating the RMSE for the training set 
RSME_icr_trainingset <- RMSE(HILIC_training_set$pred_icr_training, HILIC_training_set$RT)
RSME_LRBS_trainingset <- RMSE(HILIC_training_set$pred_LRBS_training, HILIC_training_set$RT)
RSME_pls_trainingset<- RMSE(HILIC_training_set$pred_pls_training, HILIC_training_set$RT)
RSME_MARS_trainingset <- RMSE(HILIC_training_set$pred_MARS_training, HILIC_training_set$RT)
RSME_knn_trainingset <- RMSE(HILIC_training_set$pred_knn_training, HILIC_training_set$RT)
RSME_RF_trainingset <- RMSE(HILIC_training_set$pred_RF_training, HILIC_training_set$RT)

RSME_trainingset <- c(RSME_icr_trainingset, RSME_LRBS_trainingset, RSME_pls_trainingset,  RSME_MARS_trainingset,  RSME_knn_trainingset, RSME_RF_trainingset)

Models <- c("ICR", "LRBS", "PLS",  "´MARS",  "KNN", "RF")


p13 <- ggplot() +
  geom_point(mapping = aes(x = HILIC_training_set$RT, y =HILIC_training_set$pred_RF_training), color = "#3d5a80", size = 2, alpha = 1) +
  geom_point(mapping = aes(x = HILIC_test_set$RT, y =HILIC_test_set$pred_RF_test), color = "#ee6c4d", size = 2, alpha = 1) +
  xlab("Measured " ~italic("t")["R"]~" (min)")+
  ylab("Predicted" ~italic("t")["R"]~" (min)")+
  theme_bw()
p13 + my_theme

print(RSME_RF_testset)
#1,93
print(RSME_RF_trainingset)
#0,43
#Calculation the relative error 
#Predicting tR for the training set 
ICR <- c()
LRBS <- c()
pls <- c()
MARS <- c()
knn <- c()
RF <- c()
for( i in 1:151) {
  ICR[i] = HILIC_test_set$RT[i]-HILIC_test_set$pred_icr_test[i]
  LRBS[i] = HILIC_test_set$RT[i]-HILIC_test_set$pred_LRBS_test[i]
  pls[i] = HILIC_test_set$RT[i]-HILIC_test_set$pred_pls_test[i]
  MARS[i] = HILIC_test_set$RT[i]-HILIC_test_set$pred_MARS_test[i]
  knn[i] = HILIC_test_set$RT[i]-HILIC_test_set$pred_knn_test[i]
  RF[i] =  HILIC_test_set$RT[i]-HILIC_test_set$pred_RF_test[i]
}
HILIC_test_set<- HILIC_test_set %>%
  mutate(ICR , LRBS, pls ,MARS ,knn, RF )

for( i in 1:685) {
  ICR[i] = HILIC_training_set$RT[i]-HILIC_training_set$pred_icr_training[i]
  LRBS[i] = HILIC_training_set$RT[i]-HILIC_training_set$pred_LRBS_training[i]
  pls[i] =  HILIC_training_set$RT[i]-HILIC_training_set$pred_pls_training[i]
  MARS[i] = HILIC_training_set$RT[i]-HILIC_training_set$pred_MARS_training[i]
  knn[i] = HILIC_training_set$RT[i]-HILIC_training_set$pred_knn_training[i]
  RF[i] = HILIC_training_set$RT[i]-HILIC_training_set$pred_RF_training[i]
}
HILIC_training_set<- HILIC_training_set %>%
  mutate(ICR , LRBS, pls ,MARS ,knn, RF )

HILIC_RF_training <- HILIC_training_set %>% dplyr::select(-c(ICR,LRBS, pls, MARS, knn))
names(HILIC_RF_training)[names(HILIC_RF_training) == "RF"] <- "Error"
HILIC_RF_training <- HILIC_RF_training %>% mutate(Model="RF")
ICR <- HILIC_training_set %>% dplyr::select(-c(RF,LRBS, pls, MARS, knn))
names(ICR)[names(ICR) == "ICR"] <- "Error"
ICR <- ICR %>% mutate(Model="ICR")
LRBS <- HILIC_training_set %>% dplyr::select(-c(RF,ICR, pls, MARS, knn))
names(LRBS)[names(LRBS) == "LRBS"] <- "Error"
LRBS <- LRBS %>% mutate(Model="LRBS")
pls <- HILIC_training_set %>% dplyr::select(-c(RF,ICR,LRBS, MARS, knn))
names(pls)[names(pls) == "pls"] <- "Error"
pls <- pls %>% mutate(Model="PLS")
MARS <- HILIC_training_set %>% dplyr::select(-c(RF,ICR,LRBS, pls, knn))
names(MARS)[names(MARS) == "MARS"] <- "Error"
MARS <- MARS %>% mutate(Model="MARS")
knn <- HILIC_training_set %>% dplyr::select(-c(RF,ICR,LRBS, pls,MARS))
names(knn)[names(knn) == "knn"] <- "Error"
knn <- knn %>% mutate(Model="KNN")

HILIC_RF_training <- HILIC_RF_training %>%
  bind_rows(ICR) %>% bind_rows(LRBS) %>% 
  bind_rows(pls) %>% bind_rows(MARS) %>% bind_rows(knn)

p14 <- ggplot(data = HILIC_RF_training,
              aes(x = Model, y = Error)) +
  geom_boxplot()+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               binwidth = 0.04,
               #fill = "grey",
               alpha = 0.5
  )+ ylab("Error (min)")+xlab("Models")+ylim(-10,10)
p14 +my_theme
HILIC_RF_test <- HILIC_test_set %>% dplyr::select(-c(ICR,LRBS, pls, MARS, knn))
names(HILIC_RF_test)[names(HILIC_RF_test) == "RF"] <- "Error"
HILIC_RF_test <- HILIC_RF_test %>% mutate(Model="RF")
ICR <- HILIC_test_set %>% dplyr::select(-c(RF,LRBS, pls, MARS, knn))
names(ICR)[names(ICR) == "ICR"] <- "Error"
ICR <- ICR %>% mutate(Model="ICR")
LRBS <- HILIC_test_set %>% dplyr::select(-c(RF,ICR, pls, MARS, knn))
names(LRBS)[names(LRBS) == "LRBS"] <- "Error"
LRBS <- LRBS %>% mutate(Model="LRBS")
pls <- HILIC_test_set %>% dplyr::select(-c(RF,ICR,LRBS, MARS, knn))
names(pls)[names(pls) == "pls"] <- "Error"
pls <- pls %>% mutate(Model="PLS")
MARS <- HILIC_test_set %>% dplyr::select(-c(RF,ICR,LRBS, pls, knn))
names(MARS)[names(MARS) == "MARS"] <- "Error"
MARS <- MARS %>% mutate(Model="MARS")
knn <- HILIC_test_set %>% dplyr::select(-c(RF,ICR,LRBS, pls,MARS))
names(knn)[names(knn) == "knn"] <- "Error"
knn <- knn %>% mutate(Model="KNN")

HILIC_RF_test<- HILIC_RF_test %>%
  bind_rows(ICR) %>% bind_rows(LRBS) %>% 
  bind_rows(pls) %>% bind_rows(MARS) %>% bind_rows(knn)

p15 <- ggplot(data = HILIC_RF_test,
              aes(x = Model, y = Error)) +
  geom_boxplot()+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               binwidth = 0.04, #give it different values
               #fill = "grey",
               alpha = 0.5
  )+ ylim(-10,10) +ylab("Error (min)")+xlab("Models")
p15+my_theme
varImp(regressor_RF)
Important_factors_in_RF_HILIC <- read_excel("Important factors in RF.xlsx", 
                                         sheet = "HILIC")

p16 <- ggplot(data=Important_factors_in_RF_HILIC)+
  geom_col(mapping=aes(x=reorder(Feature, Importance), y= Importance))+
  xlab("Variables")+
  ylab("Importance in the random forest model")+
  theme_bw()

p16 + my_theme
 





#---------------Biphenyl model------
dataset <- Combined %>% filter(Column==4)
dataset <- dataset %>% dplyr::select(-Column)
dataset <- dataset %>%   left_join(SMILES_unique) %>%  
  dplyr::select(SPLIT, everything()) 

for(i in 1:133) { if (is.na(dataset[i,1])) {
  dataset[i,1] <-TRUE  }
}
#trainingset 
Biphenyl_training_set <- dataset %>%   filter(SPLIT == TRUE) %>%  
  dplyr::select(-SPLIT)
#99 datapoints
#testset 
Biphenyl_test_set <- dataset %>%   filter(SPLIT == FALSE) %>%
  dplyr::select(-SPLIT) 
#33 datapoints
set.seed(100)
folds <- groupKFold(Biphenyl_training_set$SMILES, k = 20)

fitControl <- trainControl(## 10-fold CV
  method = "boot",
  index = folds)

#-----------------------Training and cross-validation -----------------------------
#install.packages("fastICA")
library(fastICA)
set.seed(825)
regressor_lm <- train(RT ~ ., 
                      data = Biphenyl_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                      method = "lm", 
                      trControl = fitControl,
                      verbose=FALSE)

print(regressor_lm)
summary(regressor_lm)
regressor_lm$finalModel
#-------------------------Independent Component Regression----------------------
#install.packages("fastICA")
library(fastICA)
set.seed(825)
regressor_icr <- train(RT ~ ., 
                       data = Biphenyl_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "icr", 
                       trControl = fitControl,
                       verbose=FALSE)

print(regressor_icr)
summary(regressor_icr)
regressor_icr$finalModel

#-------------------------Linear Regression with Backwards Selection----------------------
#install.packages("leaps")
library(leaps)
set.seed(825)
regressor_LRBS <- train(RT ~ ., 
                        data = Biphenyl_training_set%>% dplyr::select(-c(SMILES, Compound_name)), 
                        method = "leapBackward", 
                        trControl = fitControl,
                        verbose = FALSE)

print(regressor_LRBS)
summary(regressor_LRBS)
regressor_LRBS$finalModel
#-------------------------Partial least squares ----------------------
#install.packages("pls")
library(pls)
set.seed(825)
regressor_pls <- train(RT ~ ., 
                       data = Biphenyl_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "pls", 
                       trControl = fitControl,
                       verbose = FALSE)

print(regressor_pls)
summary(regressor_pls)
regressor_pls$finalModel


#-------------------------Multivariate Adaptive Regression Spline-------------------------
#install.packages("earth")
library(earth)
set.seed(825)
regressor_MARS <- train(RT ~ ., 
                        data = Biphenyl_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                        method = "earth", 
                        trControl = fitControl)

print(regressor_MARS)
summary(regressor_MARS)

#-------------------------KNN-------------------------
#install.packages("kknn")
library(kknn)
set.seed(825)
regressor_knn <- train(RT ~ ., 
                       data = Biphenyl_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "kknn", 
                       trControl = fitControl)

print(regressor_knn)
summary(regressor_knn)
#-------------------------RF-------------------------
set.seed(825)
regressor_RF <- train(RT ~ ., 
                      data = Biphenyl_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                      method = "ranger", 
                      trControl = fitControl, importance='impurity')

print(regressor_RF)
summary(regressor_RF)
varImp(regressor_RF) 
saveRDS(regressor_RF,  "regressorRF for Biphenyl 16082021.rds")


#-------------------------Validation--------------------------------
#Predicting for the test set
pred_icr_test = predict(regressor_icr, newdata = Biphenyl_test_set)
pred_LRBS_test = predict(regressor_LRBS, newdata = Biphenyl_test_set)
pred_pls_test = predict(regressor_pls, newdata = Biphenyl_test_set)
pred_MARS_test = predict(regressor_MARS, newdata = Biphenyl_test_set) 
pred_knn_test = predict(regressor_knn, newdata = Biphenyl_test_set) 
pred_RF_test = predict(regressor_RF, newdata = Biphenyl_test_set) 

#Predicting tR for the training set 
pred_icr_training = predict(regressor_icr, newdata = Biphenyl_training_set)
pred_LRBS_training = predict(regressor_LRBS, newdata = Biphenyl_training_set)
pred_pls_training = predict(regressor_pls, newdata = Biphenyl_training_set)
pred_MARS_training= predict(regressor_MARS, newdata = Biphenyl_training_set) 
pred_knn_training= predict(regressor_knn, newdata = Biphenyl_training_set) 
pred_RF_training = predict(regressor_RF, newdata = Biphenyl_training_set) 

Biphenyl_test_set<- Biphenyl_test_set %>%
  mutate(pred_icr_test, pred_LRBS_test, pred_pls_test,  pred_MARS_test,  pred_knn_test, pred_RF_test)

Biphenyl_training_set<- Biphenyl_training_set %>%
  mutate(pred_icr_training, pred_LRBS_training, pred_pls_training,  pred_MARS_training,  pred_knn_training, pred_RF_training)

#calculating the RMSE for the test set 
RSME_icr_testset <- RMSE(Biphenyl_test_set$pred_icr_test, Biphenyl_test_set$RT)
RSME_LRBS_testset <- RMSE(Biphenyl_test_set$pred_LRBS_test, Biphenyl_test_set$RT)
RSME_pls_testset <- RMSE(Biphenyl_test_set$pred_pls_test, Biphenyl_test_set$RT)
RSME_MARS_testset <- RMSE(Biphenyl_test_set$pred_MARS_test, Biphenyl_test_set$RT)
RSME_knn_testset <- RMSE(Biphenyl_test_set$pred_knn_test, Biphenyl_test_set$RT)
RSME_RF_testset <- RMSE(Biphenyl_test_set$pred_RF_test, Biphenyl_test_set$RT)

RSME_testset <- c(RSME_icr_testset, RSME_LRBS_testset, RSME_pls_testset, RSME_MARS_testset,  RSME_knn_testset, RSME_RF_testset)

#calculating the RMSE for the training set 
RSME_icr_trainingset <- RMSE(Biphenyl_training_set$pred_icr_training, Biphenyl_training_set$RT)
RSME_LRBS_trainingset <- RMSE(Biphenyl_training_set$pred_LRBS_training, Biphenyl_training_set$RT)
RSME_pls_trainingset<- RMSE(Biphenyl_training_set$pred_pls_training, Biphenyl_training_set$RT)
RSME_MARS_trainingset <- RMSE(Biphenyl_training_set$pred_MARS_training, Biphenyl_training_set$RT)
RSME_knn_trainingset <- RMSE(Biphenyl_training_set$pred_knn_training, Biphenyl_training_set$RT)
RSME_RF_trainingset <- RMSE(Biphenyl_training_set$pred_RF_training, Biphenyl_training_set$RT)

RSME_trainingset <- c(RSME_icr_trainingset, RSME_LRBS_trainingset, RSME_pls_trainingset,  RSME_MARS_trainingset,  RSME_knn_trainingset, RSME_RF_trainingset)

Models <- c("ICR", "LRBS", "PLS",  "´MARS",  "KNN", "RF")


p17 <- ggplot() +
  geom_point(mapping = aes(x = Biphenyl_training_set$RT, y =Biphenyl_training_set$pred_RF_training), color = "#3d5a80", size = 2, alpha = 1) +
  geom_point(mapping = aes(x = Biphenyl_test_set$RT, y =Biphenyl_test_set$pred_RF_test), color = "#ee6c4d", size = 2, alpha = 1) +
  xlab("Measured " ~italic("t")["R"]~" (min)")+
  ylab("Predicted" ~italic("t")["R"]~" (min)")+
  theme_bw()
p17 + my_theme

print(RSME_RF_testset)
#1.56
print(RSME_RF_trainingset)
#0.88
#Calculation the relative error 
#Predicting tR for the training set 
ICR <- c()
LRBS <- c()
pls <- c()
MARS <- c()
knn <- c()
RF <- c()
for( i in 1:23) {
  ICR[i] = Biphenyl_test_set$RT[i]-Biphenyl_test_set$pred_icr_test[i]
  LRBS[i] = Biphenyl_test_set$RT[i]-Biphenyl_test_set$pred_LRBS_test[i]
  pls[i] = Biphenyl_test_set$RT[i]-Biphenyl_test_set$pred_pls_test[i]
  MARS[i] = Biphenyl_test_set$RT[i]-Biphenyl_test_set$pred_MARS_test[i]
  knn[i] = Biphenyl_test_set$RT[i]-Biphenyl_test_set$pred_knn_test[i]
  RF[i] =  Biphenyl_test_set$RT[i]-Biphenyl_test_set$pred_RF_test[i]
}
Biphenyl_test_set<- Biphenyl_test_set %>%
  mutate(ICR , LRBS, pls ,MARS ,knn, RF )

for( i in 1:110) {
  ICR[i] = Biphenyl_training_set$RT[i]-Biphenyl_training_set$pred_icr_training[i]
  LRBS[i] = Biphenyl_training_set$RT[i]-Biphenyl_training_set$pred_LRBS_training[i]
  pls[i] =  Biphenyl_training_set$RT[i]-Biphenyl_training_set$pred_pls_training[i]
  MARS[i] = Biphenyl_training_set$RT[i]-Biphenyl_training_set$pred_MARS_training[i]
  knn[i] = Biphenyl_training_set$RT[i]-Biphenyl_training_set$pred_knn_training[i]
  RF[i] = Biphenyl_training_set$RT[i]-Biphenyl_training_set$pred_RF_training[i]
}
Biphenyl_training_set<- Biphenyl_training_set %>%
  mutate(ICR , LRBS, pls ,MARS ,knn, RF )

Biphenyl_RF_training <- Biphenyl_training_set %>% dplyr::select(-c(ICR,LRBS, pls, MARS, knn))
names(Biphenyl_RF_training)[names(Biphenyl_RF_training) == "RF"] <- "Error"
Biphenyl_RF_training <- Biphenyl_RF_training %>% mutate(Model="RF")
ICR <- Biphenyl_training_set %>% dplyr::select(-c(RF,LRBS, pls, MARS, knn))
names(ICR)[names(ICR) == "ICR"] <- "Error"
ICR <- ICR %>% mutate(Model="ICR")
LRBS <- Biphenyl_training_set %>% dplyr::select(-c(RF,ICR, pls, MARS, knn))
names(LRBS)[names(LRBS) == "LRBS"] <- "Error"
LRBS <- LRBS %>% mutate(Model="LRBS")
pls <- Biphenyl_training_set %>% dplyr::select(-c(RF,ICR,LRBS, MARS, knn))
names(pls)[names(pls) == "pls"] <- "Error"
pls <- pls %>% mutate(Model="PLS")
MARS <- Biphenyl_training_set %>% dplyr::select(-c(RF,ICR,LRBS, pls, knn))
names(MARS)[names(MARS) == "MARS"] <- "Error"
MARS <- MARS %>% mutate(Model="MARS")
knn <- Biphenyl_training_set %>% dplyr::select(-c(RF,ICR,LRBS, pls,MARS))
names(knn)[names(knn) == "knn"] <- "Error"
knn <- knn %>% mutate(Model="KNN")

Biphenyl_RF_training <- Biphenyl_RF_training %>%
  bind_rows(ICR) %>% bind_rows(LRBS) %>% 
  bind_rows(pls) %>% bind_rows(MARS) %>% bind_rows(knn)

p18 <- ggplot(data = Biphenyl_RF_training,
              aes(x = Model, y = Error)) +
  geom_boxplot()+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               binwidth = 0.04,
               #fill = "grey",
               alpha = 0.5
  )+ ylab("Error (min)")+xlab("Models")
p18 +my_theme
Biphenyl_RF_test <- Biphenyl_test_set %>% dplyr::select(-c(ICR,LRBS, pls, MARS, knn))
names(Biphenyl_RF_test)[names(Biphenyl_RF_test) == "RF"] <- "Error"
Biphenyl_RF_test <- Biphenyl_RF_test %>% mutate(Model="RF")
ICR <- Biphenyl_test_set %>% dplyr::select(-c(RF,LRBS, pls, MARS, knn))
names(ICR)[names(ICR) == "ICR"] <- "Error"
ICR <- ICR %>% mutate(Model="ICR")
LRBS <- Biphenyl_test_set %>% dplyr::select(-c(RF,ICR, pls, MARS, knn))
names(LRBS)[names(LRBS) == "LRBS"] <- "Error"
LRBS <- LRBS %>% mutate(Model="LRBS")
pls <- Biphenyl_test_set %>% dplyr::select(-c(RF,ICR,LRBS, MARS, knn))
names(pls)[names(pls) == "pls"] <- "Error"
pls <- pls %>% mutate(Model="PLS")
MARS <- Biphenyl_test_set %>% dplyr::select(-c(RF,ICR,LRBS, pls, knn))
names(MARS)[names(MARS) == "MARS"] <- "Error"
MARS <- MARS %>% mutate(Model="MARS")
knn <- Biphenyl_test_set %>% dplyr::select(-c(RF,ICR,LRBS, pls,MARS))
names(knn)[names(knn) == "knn"] <- "Error"
knn <- knn %>% mutate(Model="KNN")

Biphenyl_RF_test<- Biphenyl_RF_test %>%
  bind_rows(ICR) %>% bind_rows(LRBS) %>% 
  bind_rows(pls) %>% bind_rows(MARS) %>% bind_rows(knn)

p19 <- ggplot(data = Biphenyl_RF_test,
              aes(x = Model, y = Error)) +
  geom_boxplot()+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               binwidth = 0.04, #give it different values
               #fill = "grey",
               alpha = 0.5
  )+ ylim(-10,10) +ylab("Error (min)")+xlab("Models")
p19+my_theme
Important_factors_in_RF_Biphenyl <- read_excel("Important factors in RF.xlsx", 
                                            sheet = "biphenyl")
varImp(regressor_RF)
p20 <- ggplot(data=Important_factors_in_RF_Biphenyl)+
  geom_col(mapping=aes(x=reorder(Feature, Importance), y= Importance))+
  xlab("Variables")+
  ylab("Importance in the random forest model")+
  theme_bw()

p20 + my_theme






#---------------General model------
dataset <- Combined %>% filter(Column!=3)
dataset <- dataset %>%   left_join(SMILES_unique) %>%  
  dplyr::select(SPLIT, everything()) 

for(i in 1:2125) { if (is.na(dataset[i,1])) {
  dataset[i,1] <-TRUE  }
}
#trainingset 
General_training_set <- dataset %>%   filter(SPLIT == TRUE) %>%  
  dplyr::select(-SPLIT)
#378 datapoints
#testset 
General_test_set <- dataset %>%   filter(SPLIT == FALSE) %>%
  dplyr::select(-SPLIT) 
#1747 datapoints
set.seed(100)
folds <- groupKFold(General_training_set$SMILES, k = 20)

fitControl <- trainControl(## 10-fold CV
  method = "boot",
  index = folds)

#-----------------------Training and cross-validation -----------------------------
#install.packages("fastICA")
library(fastICA)
set.seed(825)
regressor_lm <- train(RT ~ ., 
                      data = General_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                      method = "lm", 
                      trControl = fitControl,
                      verbose=FALSE)

print(regressor_lm)
summary(regressor_lm)
regressor_lm$finalModel
#-------------------------Independent Component Regression----------------------
#install.packages("fastICA")
library(fastICA)
set.seed(825)
regressor_icr <- train(RT ~ ., 
                       data = General_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "icr", 
                       trControl = fitControl,
                       verbose=FALSE)

print(regressor_icr)
summary(regressor_icr)
regressor_icr$finalModel

#-------------------------Linear Regression with Backwards Selection----------------------
#install.packages("leaps")
library(leaps)
set.seed(825)
regressor_LRBS <- train(RT ~ ., 
                        data = General_training_set%>% dplyr::select(-c(SMILES, Compound_name)), 
                        method = "leapBackward", 
                        trControl = fitControl,
                        verbose = FALSE)

print(regressor_LRBS)
summary(regressor_LRBS)
regressor_LRBS$finalModel
#-------------------------Partial least squares ----------------------
#install.packages("pls")
library(pls)
set.seed(825)
regressor_pls <- train(RT ~ ., 
                       data = General_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "pls", 
                       trControl = fitControl,
                       verbose = FALSE)

print(regressor_pls)
summary(regressor_pls)
regressor_pls$finalModel


#-------------------------Multivariate Adaptive Regression Spline-------------------------
#install.packages("earth")
library(earth)
set.seed(825)
regressor_MARS <- train(RT ~ ., 
                        data = General_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                        method = "earth", 
                        trControl = fitControl)

print(regressor_MARS)
summary(regressor_MARS)

#-------------------------KNN-------------------------
#install.packages("kknn")
library(kknn)
set.seed(825)
regressor_knn <- train(RT ~ ., 
                       data = General_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "kknn", 
                       trControl = fitControl)

print(regressor_knn)
summary(regressor_knn)
#-------------------------RF-------------------------
set.seed(825)
regressor_RF <- train(RT ~ ., 
                      data = General_training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                      method = "ranger", 
                      trControl = fitControl, importance='impurity')

print(regressor_RF)
summary(regressor_RF)
saveRDS(regressor_RF,  "regressorRF for General 16082021.rds")


#-------------------------Validation--------------------------------
#Predicting for the test set
pred_icr_test = predict(regressor_icr, newdata = General_test_set)
pred_LRBS_test = predict(regressor_LRBS, newdata = General_test_set)
pred_pls_test = predict(regressor_pls, newdata = General_test_set)
pred_MARS_test = predict(regressor_MARS, newdata = General_test_set) 
pred_knn_test = predict(regressor_knn, newdata = General_test_set) 
pred_RF_test = predict(regressor_RF, newdata = General_test_set) 

#Predicting tR for the training set 
pred_icr_training = predict(regressor_icr, newdata = General_training_set)
pred_LRBS_training = predict(regressor_LRBS, newdata = General_training_set)
pred_pls_training = predict(regressor_pls, newdata = General_training_set)
pred_MARS_training= predict(regressor_MARS, newdata = General_training_set) 
pred_knn_training= predict(regressor_knn, newdata = General_training_set) 
pred_RF_training = predict(regressor_RF, newdata = General_training_set) 

General_test_set<- General_test_set %>%
  mutate(pred_icr_test, pred_LRBS_test, pred_pls_test,  pred_MARS_test,  pred_knn_test, pred_RF_test)

General_training_set<- General_training_set %>%
  mutate(pred_icr_training, pred_LRBS_training, pred_pls_training,  pred_MARS_training,  pred_knn_training, pred_RF_training)

#calculating the RMSE for the test set 
RSME_icr_testset <- RMSE(General_test_set$pred_icr_test, General_test_set$RT)
RSME_LRBS_testset <- RMSE(General_test_set$pred_LRBS_test, General_test_set$RT)
RSME_pls_testset <- RMSE(General_test_set$pred_pls_test, General_test_set$RT)
RSME_MARS_testset <- RMSE(General_test_set$pred_MARS_test, General_test_set$RT)
RSME_knn_testset <- RMSE(General_test_set$pred_knn_test, General_test_set$RT)
RSME_RF_testset <- RMSE(General_test_set$pred_RF_test, General_test_set$RT)

RSME_testset <- c(RSME_icr_testset, RSME_LRBS_testset, RSME_pls_testset, RSME_MARS_testset,  RSME_knn_testset, RSME_RF_testset)

#calculating the RMSE for the training set 
RSME_icr_trainingset <- RMSE(General_training_set$pred_icr_training, General_training_set$RT)
RSME_LRBS_trainingset <- RMSE(General_training_set$pred_LRBS_training, General_training_set$RT)
RSME_pls_trainingset<- RMSE(General_training_set$pred_pls_training, General_training_set$RT)
RSME_MARS_trainingset <- RMSE(General_training_set$pred_MARS_training, General_training_set$RT)
RSME_knn_trainingset <- RMSE(General_training_set$pred_knn_training, General_training_set$RT)
RSME_RF_trainingset <- RMSE(General_training_set$pred_RF_training, General_training_set$RT)

RSME_trainingset <- c(RSME_icr_trainingset, RSME_LRBS_trainingset, RSME_pls_trainingset,  RSME_MARS_trainingset,  RSME_knn_trainingset, RSME_RF_trainingset)

Models <- c("ICR", "LRBS", "PLS",  "´MARS",  "KNN", "RF")


p21 <- ggplot() +
  geom_point(mapping = aes(x = General_training_set$RT, y =General_training_set$pred_RF_training), color = "#3d5a80", size = 2, alpha = 1) +
  geom_point(mapping = aes(x = General_test_set$RT, y =General_test_set$pred_RF_test), color = "#ee6c4d", size = 2, alpha = 1) +
  xlab("Measured " ~italic("t")["R"]~" (min)")+
  ylab("Predicted" ~italic("t")["R"]~" (min)")+
  theme_bw()
p21 + my_theme

print(RSME_RF_testset)
#1.66
print(RSME_RF_trainingset)
#0,47
#Calculation the relative error 
#Predicting tR for the training set 
ICR <- c()
LRBS <- c()
pls <- c()
MARS <- c()
knn <- c()
RF <- c()
for( i in 1:378) {
  ICR[i] = General_test_set$RT[i]-General_test_set$pred_icr_test[i]
  LRBS[i] = General_test_set$RT[i]-General_test_set$pred_LRBS_test[i]
  pls[i] = General_test_set$RT[i]-General_test_set$pred_pls_test[i]
  MARS[i] = General_test_set$RT[i]-General_test_set$pred_MARS_test[i]
  knn[i] = General_test_set$RT[i]-General_test_set$pred_knn_test[i]
  RF[i] =  General_test_set$RT[i]-General_test_set$pred_RF_test[i]
}
General_test_set<- General_test_set %>%
  mutate(ICR , LRBS, pls ,MARS ,knn, RF )

for( i in 1:1747) {
  ICR[i] = General_training_set$RT[i]-General_training_set$pred_icr_training[i]
  LRBS[i] = General_training_set$RT[i]-General_training_set$pred_LRBS_training[i]
  pls[i] =  General_training_set$RT[i]-General_training_set$pred_pls_training[i]
  MARS[i] = General_training_set$RT[i]-General_training_set$pred_MARS_training[i]
  knn[i] = General_training_set$RT[i]-General_training_set$pred_knn_training[i]
  RF[i] = General_training_set$RT[i]-General_training_set$pred_RF_training[i]
}
General_training_set<- General_training_set %>%
  mutate(ICR , LRBS, pls ,MARS ,knn, RF )

General_RF_training <- General_training_set %>% dplyr::select(-c(ICR,LRBS, pls, MARS, knn))
names(General_RF_training)[names(General_RF_training) == "RF"] <- "Error"
General_RF_training <- General_RF_training %>% mutate(Model="RF")
ICR <- General_training_set %>% dplyr::select(-c(RF,LRBS, pls, MARS, knn))
names(ICR)[names(ICR) == "ICR"] <- "Error"
ICR <- ICR %>% mutate(Model="ICR")
LRBS <- General_training_set %>% dplyr::select(-c(RF,ICR, pls, MARS, knn))
names(LRBS)[names(LRBS) == "LRBS"] <- "Error"
LRBS <- LRBS %>% mutate(Model="LRBS")
pls <- General_training_set %>% dplyr::select(-c(RF,ICR,LRBS, MARS, knn))
names(pls)[names(pls) == "pls"] <- "Error"
pls <- pls %>% mutate(Model="PLS")
MARS <- General_training_set %>% dplyr::select(-c(RF,ICR,LRBS, pls, knn))
names(MARS)[names(MARS) == "MARS"] <- "Error"
MARS <- MARS %>% mutate(Model="MARS")
knn <- General_training_set %>% dplyr::select(-c(RF,ICR,LRBS, pls,MARS))
names(knn)[names(knn) == "knn"] <- "Error"
knn <- knn %>% mutate(Model="KNN")

General_RF_training <- General_RF_training %>%
  bind_rows(ICR) %>% bind_rows(LRBS) %>% 
  bind_rows(pls) %>% bind_rows(MARS) %>% bind_rows(knn)

p22 <- ggplot(data = General_RF_training,
              aes(x = Model, y = Error)) +
  geom_boxplot()+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               binwidth = 0.04,
               #fill = "grey",
               alpha = 0.5
  )+ ylab("Error (min)")+xlab("Models")
p22 +my_theme
General_RF_test <- General_test_set %>% dplyr::select(-c(ICR,LRBS, pls, MARS, knn))
names(General_RF_test)[names(General_RF_test) == "RF"] <- "Error"
General_RF_test <- General_RF_test %>% mutate(Model="RF")
ICR <- General_test_set %>% dplyr::select(-c(RF,LRBS, pls, MARS, knn))
names(ICR)[names(ICR) == "ICR"] <- "Error"
ICR <- ICR %>% mutate(Model="ICR")
LRBS <- General_test_set %>% dplyr::select(-c(RF,ICR, pls, MARS, knn))
names(LRBS)[names(LRBS) == "LRBS"] <- "Error"
LRBS <- LRBS %>% mutate(Model="LRBS")
pls <- General_test_set %>% dplyr::select(-c(RF,ICR,LRBS, MARS, knn))
names(pls)[names(pls) == "pls"] <- "Error"
pls <- pls %>% mutate(Model="PLS")
MARS <- General_test_set %>% dplyr::select(-c(RF,ICR,LRBS, pls, knn))
names(MARS)[names(MARS) == "MARS"] <- "Error"
MARS <- MARS %>% mutate(Model="MARS")
knn <- General_test_set %>% dplyr::select(-c(RF,ICR,LRBS, pls,MARS))
names(knn)[names(knn) == "knn"] <- "Error"
knn <- knn %>% mutate(Model="KNN")

General_RF_test<- General_RF_test %>%
  bind_rows(ICR) %>% bind_rows(LRBS) %>% 
  bind_rows(pls) %>% bind_rows(MARS) %>% bind_rows(knn)

p23 <- ggplot(data = General_RF_test,
              aes(x = Model, y = Error)) +
  geom_boxplot()+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               binwidth = 0.04, #give it different values
               #fill = "grey",
               alpha = 0.5
  )+ ylim(-10,10) +ylab("Error (min)")+xlab("Models")
p23+my_theme
varImp(regressor_RF) 
Important_factors_in_RF_General <- read_excel("Important factors in RF.xlsx", 
                                               sheet = "General")

p24 <- ggplot(data=Important_factors_in_RF_General)+
  geom_col(mapping=aes(x=reorder(Feature, Importance), y= Importance))+
  xlab("Variables")+
  ylab("Importance in the random forest model")+
  theme_bw()

p24 + my_theme
library(cowplot)
plot_grid(p2+my_theme, p5+my_theme1, p9+my_theme, p12+my_theme1, p13+my_theme, p16+my_theme1, p17+my_theme, p20+my_theme1, p21+my_theme, p24+my_theme1,
          labels =c("a) Predictions for mixed mode", "b) Important variables  for mixed mode ",
                    "c) Predictions for reversed-C18", "d) Important variables for reversed-C18",
                    "e) Predictions for HILIC", "f) Important variables for HILIC", 
                    "g) Predictions for biphenyl", "h) Important variables for biphenyl", 
                    "i) Predictions for the combined dataset", "j) Important variables for combined dataset"),
          align = "h",
          hjust = 0.1, vjust = 1.5,
          label_x = 0.36,
          label_y = 0.999,
          nrow = 5,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("Figure 1 in manuscript -style 1.svg",
       width = 35,
       height = 90/668*380,
       units = "cm")
my_theme1 <- theme(plot.background = element_blank(), 
                  panel.background = element_blank(), axis.line=element_line(  size=0.5), 
                  legend.background = element_blank(),
                  legend.title = element_text(),
                  legend.position = c(0.9,0.3), 
                  axis.text.x = element_text(angle =10, vjust = 0.5, hjust=0.4),
                  aspect.ratio = 1, text = element_text(family = font, size=fontsize),
                  plot.title = element_text(color=basecolor, size=14, face="bold"),
                  legend.text=element_text(family = font, size=fontsize, color=basecolor),
                  axis.text=element_text(family = font, size=fontsize, color=basecolor), 
                  panel.grid.major =  element_blank(),
                  panel.grid.minor=element_blank(), 
                  panel.border = element_blank())
my_theme2 <- theme(plot.background = element_blank(), 
                   panel.background = element_blank(), axis.line=element_line(  size=0.5), 
                   legend.background = element_blank(),
                   legend.title = element_text(),
                   legend.position = c(0.9,0.3), 
                   axis.text.x = element_text(angle =45, vjust = 0.5, hjust=0.5),
                   aspect.ratio = 1, text = element_text(family = font, size=fontsize),
                   plot.title = element_text(color=basecolor, size=14, face="bold"),
                   legend.text=element_text(family = font, size=fontsize, color=basecolor),
                   axis.text=element_text(family = font, size=fontsize, color=basecolor), 
                   panel.grid.major =  element_blank(),
                   panel.grid.minor=element_blank(), 
                   panel.border = element_blank())
plot_grid(p2+my_theme,  p9+my_theme , p13+my_theme,  p17+my_theme, p21+my_theme,
          labels =c("a) Mixed mode", "b) Reversed phase","c) HILIC column",
                    "d) Biphenyl column","e) Combined dataset"),
          align = "h",
          hjust =  0.1, vjust = 1.5 ,
          label_x = 0.3,
          label_y = 0.999,
          nrow = 2,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("Figure 1 in manuscript -style 2.svg", height = 9, width=16)
plot_grid(p5+my_theme1,  p12+my_theme1 , p16+my_theme1,  p20+my_theme1, p24+my_theme1,
          labels =c("a) Mixed mode", "b) Reversed phase","c) HILIC column",
                    "d) Biphenyl column","e) Combined dataset"),
          align = "h",
          hjust =  0.1, vjust = 1.5 ,
          label_x = 0.3,
          label_y = 0.999,
          nrow = 2,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("Figure S7 -style 2.svg", height = 9, width=16)


plot_grid(p3+my_theme,  p10+my_theme , p14+my_theme,  p18+my_theme, p22+my_theme, p4+my_theme ,p11+my_theme , p15+my_theme , p19+my_theme ,p23+my_theme ,
          labels =c("a) Taining set mixed mode", "b) Taining set reversed-C18","c) Training set HILIC",
                    "d) Training set biphenyl","e) Training set combined dataset","f) Test set mixed mode", "g) Test set reversed-C18",
                    "h) Test set HILIC", 
                    "i) Test set biphenyl", 
                    "j) Test set combined dataset"),
          align = "h",
          hjust =  0.1, vjust = 1.5 ,
          label_x = 0.3,
          label_y = 0.999,
          nrow = 2,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("Figure 2 in manuscript -style 2.svg", height = 8, width=20)

plot_grid(p6+my_theme,  p7+my_theme , 
          labels =c("a) Reversed phase ", "b) Mixed mode "),
          align = "h",
          hjust =  0.1, vjust = 1.5 ,
          label_x = 0.3,
          label_y = 0.999,
          nrow = 1,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("Linear regression for MM and RP.svg", height = 5, width=10)
#--------------- calculating the RMSE for each LC condition -----------------
#-------RP------
training2ACN <- RP_training_set %>% filter(pH == 2.7 & Organic_modifier==1)
test2ACN <- RP_test_set %>% filter(pH == 2.7 & Organic_modifier==1)
RMSE_train2ACN <- RMSE(training2ACN$pred_RF_training, training2ACN$RT)
RMSE_test2ACN <- RMSE(test2ACN$pred_RF_test, test2ACN$RT)


training2MEOH <- RP_training_set %>% filter(pH == 2.7 & Organic_modifier==2)
test2MEOH<- RP_test_set %>% filter(pH == 2.7 & Organic_modifier==2)
RMSE_train2MEOH <- RMSE(training2MEOH$pred_RF_training, training2MEOH$RT)
RMSE_test2MEOH <- RMSE(test2MEOH$pred_RF_test, test2MEOH$RT)


training5ACN <- RP_training_set %>% filter(pH == 5 & Organic_modifier==1)
test5ACN <- RP_test_set %>% filter(pH == 5 & Organic_modifier==1)
RMSE_train5ACN <- RMSE(training5ACN$pred_RF_training, training5ACN$RT)
RMSE_test5ACN <- RMSE(test5ACN$pred_RF_test, test5ACN$RT)


training5MEOH <- RP_training_set %>% filter(pH == 5 & Organic_modifier==2)
test5MEOH<- RP_test_set %>% filter(pH == 5 & Organic_modifier==2)
RMSE_train5MEOH <- RMSE(training5MEOH$pred_RF_training, training5MEOH$RT)
RMSE_test5MEOH <- RMSE(test5MEOH$pred_RF_test, test5MEOH$RT)

training8ACN <- RP_training_set %>% filter(pH == 8 & Organic_modifier==1)
test8ACN <- RP_test_set %>% filter(pH == 8 & Organic_modifier==1)
RMSE_train8ACN <- RMSE(training8ACN$pred_RF_training, training8ACN$RT)
RMSE_test8ACN <- RMSE(test8ACN$pred_RF_test, test8ACN$RT)


training8MEOH <- RP_training_set %>% filter(pH == 8 & Organic_modifier==2)
test8MEOH<- RP_test_set %>% filter(pH == 8 & Organic_modifier==2)
RMSE_train8MEOH <- RMSE(training8MEOH$pred_RF_training, training8MEOH$RT)
RMSE_test8MEOH <- RMSE(test8MEOH$pred_RF_test, test8MEOH$RT)


training10ACN <- RP_training_set %>% filter(pH == 10 & Organic_modifier==1)
test10ACN <- RP_test_set %>% filter(pH == 10 & Organic_modifier==1)
RMSE_train10ACN <- RMSE(training10ACN$pred_RF_training, training10ACN$RT)
RMSE_test10ACN <- RMSE(test10ACN$pred_RF_test, test10ACN$RT)


training10MEOH <- RP_training_set %>% filter(pH == 10 & Organic_modifier==2)
test10MEOH<- RP_test_set %>% filter(pH == 10 & Organic_modifier==2)
RMSE_train10MEOH <- RMSE(training10MEOH$pred_RF_training, training10MEOH$RT)
RMSE_test10MEOH <- RMSE(test10MEOH$pred_RF_test, test10MEOH$RT)


RMSE <- c("2.7ACN", "2.7MEOH", "5ACN",  "5MEOH",  "8ACN", "8MEOH", "10ACN", "10MEOH")
RMSE_training_RP <- c(RMSE_train2ACN, RMSE_train2MEOH, RMSE_train5ACN, RMSE_train5MEOH, RMSE_train8ACN, RMSE_train8MEOH, RMSE_train10ACN, RMSE_train10MEOH)
RMSE_test_RP <- c(RMSE_test2ACN, RMSE_test2MEOH, RMSE_test5ACN, RMSE_test5MEOH, RMSE_test8ACN, RMSE_test8MEOH, RMSE_test10ACN, RMSE_test10MEOH)

p <-ggplot() +
  geom_point(mapping = aes(x = RMSE, y =RMSE_training_RP), color = "blue", size = 2, alpha = 0.25) +
  geom_point(mapping = aes(x = RMSE, y =RMSE_test_RP),color = "red",  size = 2, alpha = 0.25) +
  ylab("RMSE Values")+
  theme_bw()
p+ my_theme

#------General model -----
training2ACN <- General_training_set %>% filter(pH == 2.7 & Organic_modifier==1)
test2ACN <- General_test_set %>% filter(pH == 2.7 & Organic_modifier==1)
RMSE_train2ACN <- RMSE(training2ACN$pred_RF_training, training2ACN$RT)
RMSE_test2ACN <- RMSE(test2ACN$pred_RF_test, test2ACN$RT)


training2MEOH <- General_training_set %>% filter(pH == 2.7 & Organic_modifier==2)
test2MEOH<- General_test_set %>% filter(pH == 2.7 & Organic_modifier==2)
RMSE_train2MEOH <- RMSE(training2MEOH$pred_RF_training, training2MEOH$RT)
RMSE_test2MEOH <- RMSE(test2MEOH$pred_RF_test, test2MEOH$RT)


training5ACN <- General_training_set %>% filter(pH == 5 & Organic_modifier==1)
test5ACN <- General_test_set %>% filter(pH == 5 & Organic_modifier==1)
RMSE_train5ACN <- RMSE(training5ACN$pred_RF_training, training5ACN$RT)
RMSE_test5ACN <- RMSE(test5ACN$pred_RF_test, test5ACN$RT)


training5MEOH <- General_training_set %>% filter(pH == 5 & Organic_modifier==2)
test5MEOH<- General_test_set %>% filter(pH == 5 & Organic_modifier==2)
RMSE_train5MEOH <- RMSE(training5MEOH$pred_RF_training, training5MEOH$RT)
RMSE_test5MEOH <- RMSE(test5MEOH$pred_RF_test, test5MEOH$RT)

training8ACN <- General_training_set %>% filter(pH == 8 & Organic_modifier==1)
test8ACN <- General_test_set %>% filter(pH == 8 & Organic_modifier==1)
RMSE_train8ACN <- RMSE(training8ACN$pred_RF_training, training8ACN$RT)
RMSE_test8ACN <- RMSE(test8ACN$pred_RF_test, test8ACN$RT)


training8MEOH <- General_training_set %>% filter(pH == 8 & Organic_modifier==2)
test8MEOH<- General_test_set %>% filter(pH == 8 & Organic_modifier==2)
RMSE_train8MEOH <- RMSE(training8MEOH$pred_RF_training, training8MEOH$RT)
RMSE_test8MEOH <- RMSE(test8MEOH$pred_RF_test, test8MEOH$RT)


training10ACN <- General_training_set %>% filter(pH == 10 & Organic_modifier==1)
test10ACN <- General_test_set %>% filter(pH == 10 & Organic_modifier==1)
RMSE_train10ACN <- RMSE(training10ACN$pred_RF_training, training10ACN$RT)
RMSE_test10ACN <- RMSE(test10ACN$pred_RF_test, test10ACN$RT)


training10MEOH <- General_training_set %>% filter(pH == 10 & Organic_modifier==2)
test10MEOH<- General_test_set %>% filter(pH == 10 & Organic_modifier==2)
RMSE_train10MEOH <- RMSE(training10MEOH$pred_RF_training, training10MEOH$RT)
RMSE_test10MEOH <- RMSE(test10MEOH$pred_RF_test, test10MEOH$RT)


RMSE <- c("2.7ACN", "2.7MEOH", "5ACN",  "5MEOH",  "8ACN", "8MEOH", "10ACN", "10MEOH")
RMSE_training_general <- c(RMSE_train2ACN, RMSE_train2MEOH, RMSE_train5ACN, RMSE_train5MEOH, RMSE_train8ACN, RMSE_train8MEOH, RMSE_train10ACN, RMSE_train10MEOH)
RMSE_test_general <- c(RMSE_test2ACN, RMSE_test2MEOH, RMSE_test5ACN, RMSE_test5MEOH, RMSE_test8ACN, RMSE_test8MEOH, RMSE_test10ACN, RMSE_test10MEOH)

p <-ggplot() +
  geom_point(mapping = aes(x = RMSE, y =RMSE_training_general), color = "blue", size = 2, alpha = 0.25) +
  geom_point(mapping = aes(x = RMSE, y =RMSE_test_general),color = "red",  size = 2, alpha = 0.25) +
  ylab("RMSE Values")+
  theme_bw()
p+ my_theme
