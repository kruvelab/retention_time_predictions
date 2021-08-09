#-------------------------- Libraries -----------------------------------------------
library(caret) #machine learinign workflow
library(leaps) #library for stepwise regression
library(MASS) #contains some important linear regression tools
library(caTools) #sample split is from this package
library(tidyverse) #helps us to write concise code
library(rcdk)  
library(Metrics)
#------------------------importing the data---------------------------------------------- 
dataset <-  read_delim('descs200502.csv',
                       delim = ",",
                       col_names = TRUE,
                       trim_ws = TRUE) %>%
  drop_na()
dataset2 <-  read_delim('Filter_RP_2010.csv',
                        delim = ",",
                        col_names = TRUE,
                        trim_ws = TRUE) %>%
  drop_na()
SMILES_unique <-  read_delim('SMILES_unique.csv',
                             delim = ",",
                             col_names = TRUE,
                             trim_ws = TRUE) %>%
  drop_na()
dataset3 <- dataset2 %>%
  group_by(Compound_name, Column, Organic_modifier, pH, Buffer, SMILES ) %>%
  summarise(slope = mean(slope),
            RT= mean(ret_time))

dataset2 <- dataset3 %>% 
  left_join(dataset)          

#Encoding categorical data
dataset2$Organic_modifier=factor(dataset2$Organic_modifier,
                                 levels = c('Acetonitrile','Methanol'),
                                 labels=c(1,2))
dataset2$Buffer=factor(dataset2$Buffer, levels = c('Formic acid', 'Ammonia', 'TFA', 'Acetic acid','Acetate', 'Formiate','Bicarbonate'),
                       labels=c(1,2,3,4,5,6,7))

# this will be removed, we will keep only M+H
#First model the RT, and use the ionization effenciencies model. 
dataset <- dataset2 %>%
  dplyr::select(-c(SMILES,slope, Column, Compound_name, RT))
#Removing the non-numeric variables 
dataset <- dataset[3:1222]

dataset3 <- dataset %>%
  dplyr::select(-c(Organic_modifier, Buffer,  bruto_formula))

dataset3 <- dataset3 %>%
  dplyr::select(-nearZeroVar(dataset3)) 

#Calculating  the correlation between features
correlationMatrix <- cor(dataset3, use = "complete.obs")
# Finding the  highly corrected features (ideally >0.8)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.8)

dataset3 <- dataset3 %>%
  dplyr::select(-highlyCorrelated) %>%
  dplyr::mutate(Organic_modifier= dataset$Organic_modifier, Buffer=dataset$Buffer)

dataset <- dataset3

dataset <- dataset %>%
  dplyr::mutate(RT= dataset2$RT, Compound_name=dataset2$Compound_name, SMILES = dataset2$SMILES)

#------ alpha calculations------------- 
dataset4 <-  read_delim('pka and logP for MM.csv',
                        delim = ",",
                        col_names = TRUE,
                        trim_ws = TRUE) %>%
  drop_na()
fitControl <- trainControl(method = "repeatedcv", number = 5, repeats = 5) 
regressor_lm_pred_logDan <- train(`logD A-` ~`logP AH`, 
                                  data = dataset4, 
                                  method = "lm", 
                                  trControl = fitControl)
print(regressor_lm_pred_logDan)
summary(regressor_lm_pred_logDan)
#R2=86,59

regressor_lm_pred_logDcat <- train(`logD AH2+` ~`logP AH`, 
                                   data = dataset4, 
                                   method = "lm", 
                                   trControl = fitControl)
print(regressor_lm_pred_logDcat)
summary(regressor_lm_pred_logDcat)
#R2=76,73


dataset4 <-  read_delim('pka and logP for MM.csv',
                        delim = ",",
                        col_names = TRUE,
                        trim_ws = TRUE)
logDA <- dataset4[7]
for(i in 1:77) { if (is.na(logDA[i,1])) {
  logDA[i,1] <-predict(regressor_lm_pred_logDan,dataset4[i,8])  }
}
logDAH2 <- dataset4[9]
for(i in 1:77) { if (is.na(logDAH2[i,1])) {
  logDAH2[i,1] <-predict(regressor_lm_pred_logDcat,dataset4[i,8])  }
}
pkaca <- dataset4[6]
pkaca[is.na(pkaca)] <--5
pkaan <- dataset4[5]
pkaan[is.na(pkaan)] <-20
logP <- dataset4[8]

dataset4 <- dataset4[1]
dataset4 <- data.frame(dataset4, pkaca, pkaan, logDA, logP, logDAH2 )


dataset <- dataset %>% 
  left_join(dataset4)

dataset <- dataset %>%  
  drop_na()

dataset <- data.frame(dataset, alphaAH2= 1/(10^(dataset$pH - dataset$pKa.base)+1))
dataset <- data.frame(dataset, alphaA_H= (10^(dataset$pH-dataset$pKa.acid))/(10^(dataset$pH-dataset$pKa.acid )+1))
dataset <- data.frame(dataset, alphaAH= 1-(dataset$alphaA_H +dataset$alphaAH2))

dataset <- data.frame(dataset, DA= 10^dataset$logD.A.)
dataset <- data.frame(dataset, P=10^dataset$logP.AH)
dataset <- data.frame(dataset, DAH2 = 10^dataset$logD.AH2. )

dataset <- data.frame(dataset, P_sum = dataset$alphaA_H*dataset$DA +dataset$alphaAH*dataset$P + dataset$alphaAH2*dataset$DAH2)
dataset <- data.frame(dataset, log_P_sum = log10(dataset$P_sum))

for(i in 1:807) { if (is.nan(dataset[i,213])) {
  dataset[i,213] <-dataset[i,205]  }
}

dataset <- dataset %>% dplyr::select(-c(alphaAH2, alphaA_H, alphaAH, DA, P, DAH2, P_sum, pKa.base, pKa.acid, logD.A., logP.AH,logD.AH2.))

dataset <- dataset %>%   left_join(SMILES_unique) %>%  
  dplyr::select(SPLIT, everything()) 
for(i in 1:807) { if (is.na(dataset[i,1])) {
  dataset[i,1] <-TRUE  }
}
#trainingset 
training_set <- dataset %>%   filter(SPLIT == TRUE) %>%  
  dplyr::select(-SPLIT)  
#testset 
test_set <- dataset %>%   filter(SPLIT == FALSE) %>%
  dplyr::select(-SPLIT) 

folds <- groupKFold(training_set$SMILES, k = 20)
# training_set <- training_set %>%
#   dplyr::select(-c(SMILES, Compound_name))
# test_set <- test_set %>%
#   dplyr::select(-c(SMILES, Compound_name))
# 
# #Trainign the model
fitControl <- trainControl(## 10-fold CV
  method = "boot",
  index = folds)
#setting the cross-validation parameters
# fitControl <- trainControl(method = "repeatedcv", number = 5, repeats = 5) 

#-------------------------LM-------------------------
set.seed(825)
regressor_lm <- train(RT ~ log_P_sum, 
                      data = training_set, 
                      method = "lm", 
                      trControl = fitControl)

print(regressor_lm)
summary(regressor_lm)
library(extrafont)
font_import()
font <- choose_font("Arial")
fontsize <- 12
basecolor <- "black"
my_theme <- theme(plot.background = element_blank(), 
                  panel.background = element_blank(), axis.line=element_line(  size=0.5), 
                  legend.background = element_blank(),
                  legend.title = element_text(),
                  legend.position = c(1,0.3), 
                  aspect.ratio = 1, text = element_text(family = font, size=fontsize),
                  plot.title = element_text(color=basecolor, size=14, face="bold"),
                  legend.text=element_text(family = font, size=fontsize, color=basecolor),
                  axis.text=element_text(family = font, size=fontsize, color=basecolor), 
                  panel.grid.major =  element_blank(),
                  panel.grid.minor=element_blank(), 
                  panel.border = element_blank())

dataset_lm <- dataset %>% 
  mutate(tR_pred = predict(regressor_lm, newdata = dataset))
dataset_lm <- dataset_lm %>%
  mutate(SPLIT)
p1 <- ggplot(data=dataset_lm)+
  geom_point(mapping=aes(x=RT, y= tR_pred, color=SPLIT), size=3, alpha=1)+
  scale_color_manual(values=c("#ee6c4d","#3d5a80")) +
  xlab("Measured"~ t[R]~ "(min)")+
  ylab("Predicted"~ t[R]~ "(min)")
p1 +my_theme
# ggsave("RP linear regression with calculated logP_sum.svg", 
#        width=10,
#        height=10, 
#        units="cm")


#-----------------------Training and cross-validation -----------------------------
#install.packages("fastICA")
library(fastICA)
set.seed(825)
regressor_lm <- train(RT ~ ., 
                       data = training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
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
                       data = training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
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
                        data = training_set%>% dplyr::select(-c(SMILES, Compound_name)), 
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
                       data = training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
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
                        data = training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                        method = "earth", 
                        trControl = fitControl)

print(regressor_MARS)
summary(regressor_MARS)

#-------------------------KNN-------------------------
#install.packages("kknn")
library(kknn)
set.seed(825)
regressor_knn <- train(RT ~ ., 
                       data = training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                       method = "kknn", 
                       trControl = fitControl)

print(regressor_knn)
summary(regressor_knn)
#-------------------------RF-------------------------
set.seed(825)
regressor_RF <- train(RT ~ ., 
                      data = training_set %>% dplyr::select(-c(SMILES, Compound_name)), 
                      method = "ranger", 
                      trControl = fitControl, importance='impurity')

print(regressor_RF)
summary(regressor_RF)
varImp(regressor_RF) 

# #-------------------------RF with xgbtree-------------------------
# set.seed(825)
# regressor_RF2 <- train(RT ~ ., 
#                       data = training_set  %>% dplyr::select(-c(SMILES, Compound_name)), 
#                       method = "RRF", 
#                       trControl = fitControl, importance='impurity')
# 
# print(regressor_RF2)
# summary(regressor_RF2)
# varImp(regressor_RF2) 
# #-------------------------RF with xgbtree-------------------------
# set.seed(825)
# regressor_RF3 <- train(RT ~ log_P_sum, 
#                        data = training_set  %>% dplyr::select(-c(SMILES, Compound_name)), 
#                        method = "RRF", 
#                        trControl = fitControl)
# 
# print(regressor_RF3)
# summary(regressor_RF3)
# #-------------------------RF with xgbtree-------------------------
# set.seed(825)
# regressor_RF4 <- train(RT ~ log_P_sum, 
#                        data = training_set  %>% dplyr::select(-c(SMILES, Compound_name)), 
#                        method = "xgbTree", 
#                        trControl = fitControl)
# 
# print(regressor_RF4)
# summary(regressor_RF4)

#-------------------------Validation--------------------------------
#Predicting for the test set
pred_icr_test = predict(regressor_icr, newdata = test_set)
pred_LRBS_test = predict(regressor_LRBS, newdata = test_set)
pred_pls_test = predict(regressor_pls, newdata = test_set)
pred_MARS_test = predict(regressor_MARS, newdata = test_set) 
pred_knn_test = predict(regressor_knn, newdata = test_set) 
pred_RF_test = predict(regressor_RF, newdata = test_set) 
# pred_RF2_test = predict(regressor_RF2, newdata = test_set) 
# pred_RF3_test = predict(regressor_RF3, newdata = test_set) 
# pred_RF4_test = predict(regressor_RF4, newdata = test_set) 

#Predicting tR for the training set 
pred_icr_training = predict(regressor_icr, newdata = training_set)
pred_LRBS_training = predict(regressor_LRBS, newdata = training_set)
pred_pls_training = predict(regressor_pls, newdata = training_set)
pred_MARS_training= predict(regressor_MARS, newdata = training_set) 
pred_knn_training= predict(regressor_knn, newdata = training_set) 
pred_RF_training = predict(regressor_RF, newdata = training_set) 
# pred_RF2_training = predict(regressor_RF2, newdata = training_set)
# pred_RF3_training = predict(regressor_RF3, newdata = training_set) 
# pred_RF4_training = predict(regressor_RF4, newdata = training_set) 

test_set<- test_set %>%
  mutate(pred_icr_test, pred_LRBS_test, pred_pls_test,  pred_MARS_test,  pred_knn_test, pred_RF_test)

training_set<- training_set %>%
  mutate(pred_icr_training, pred_LRBS_training, pred_pls_training,  pred_MARS_training,  pred_knn_training, pred_RF_training)

#calculating the RMSE for the test set 
RSME_icr_testset <- RMSE(test_set$pred_icr_test, test_set$RT)
RSME_LRBS_testset <- RMSE(test_set$pred_LRBS_test, test_set$RT)
RSME_pls_testset <- RMSE(test_set$pred_pls_test, test_set$RT)
RSME_MARS_testset <- RMSE(test_set$pred_MARS_test, test_set$RT)
RSME_knn_testset <- RMSE(test_set$pred_knn_test, test_set$RT)
RSME_RF_testset <- RMSE(test_set$pred_RF_test, test_set$RT)
# RSME_RF2_testset <- RMSE(test_set$pred_RF2_test, test_set$RT)
# RSME_RF3_testset <- RMSE(test_set$pred_RF3_test, test_set$RT)
# RSME_RF4_testset <- RMSE(test_set$pred_RF4_test, test_set$RT)

RSME_testset <- c(RSME_icr_testset, RSME_LRBS_testset, RSME_pls_testset, RSME_MARS_testset,  RSME_knn_testset, RSME_RF_testset)

#calculating the RMSE for the training set 
RSME_icr_trainingset <- RMSE(training_set$pred_icr_training, training_set$RT)
RSME_LRBS_trainingset <- RMSE(training_set$pred_LRBS_training, training_set$RT)
RSME_pls_trainingset<- RMSE(training_set$pred_pls_training, training_set$RT)
RSME_MARS_trainingset <- RMSE(training_set$pred_MARS_training, training_set$RT)
RSME_knn_trainingset <- RMSE(training_set$pred_knn_training, training_set$RT)
RSME_RF_trainingset <- RMSE(training_set$pred_RF_training, training_set$RT)
# RSME_RF2_trainingset <- RMSE(training_set$pred_RF2_training, training_set$RT)
# RSME_RF3_trainingset <- RMSE(training_set$pred_RF3_training, training_set$RT)
# RSME_RF4_trainingset <- RMSE(training_set$pred_RF4_training, training_set$RT)

RSME_trainingset <- c(RSME_icr_trainingset, RSME_LRBS_trainingset, RSME_pls_trainingset,  RSME_MARS_trainingset,  RSME_knn_trainingset, RSME_RF_trainingset)

Models <- c("ICR", "LRBS", "PLS",  "´MARS",  "KNN", "RF")


p2<- ggplot() +
  geom_point(mapping = aes(x = Models, y =RSME_trainingset), color = "#3d5a80", size = 2, alpha = 1) +
  geom_point(mapping = aes(x = Models, y =RSME_testset),color = "#ee6c4d",  size = 2, alpha = 1) +
  xlab("Models")+
  ylab("RMSE (min)")+
  theme_bw()
p2+my_theme
ggsave("RP with calculated logP_sum with xgbtree.svg", 
       width=10,
       height=10, 
       units="cm")

p3 <- ggplot() +
  geom_point(mapping = aes(x = training_set$RT, y =training_set$pred_RF_training), color = "#3d5a80", size = 2, alpha = 1) +
  geom_point(mapping = aes(x = test_set$RT, y =test_set$pred_RF_test), color = "#ee6c4d", size = 2, alpha = 1) +
  xlab("Measured " ~t[R]~" (min)")+
  ylab("Predicted" ~t[R]~" (min)")+ 
  theme_bw()
p3+my_theme
ggsave("RP MODE RF with calculated logP_sum.svg", 
       width=10,
       height=10, 
       units="cm")
print(RSME_RF_testset)
#1,83
print(RSME_RF_trainingset)
#0,51
#Calculation the relative error 
#Predicting tR for the training set 
ICR <- c()
LRBS <- c()
pls <- c()
MARS <- c()
knn <- c()
RF <- c()
for( i in 1:142) {
  ICR[i] = test_set$RT[i]-test_set$pred_icr_test[i]
  LRBS[i] = test_set$RT[i]-test_set$pred_LRBS_test[i]
  pls[i] = test_set$RT[i]-test_set$pred_pls_test[i]
  MARS[i] = test_set$RT[i]-test_set$pred_MARS_test[i]
  knn[i] = test_set$RT[i]-test_set$pred_knn_test[i]
  RF[i] =  test_set$RT[i]-test_set$pred_RF_test[i]
}
test_set<- test_set %>%
  mutate(ICR , LRBS, pls ,MARS ,knn, RF )

for( i in 1:665) {
  ICR[i] = training_set$RT[i]-training_set$pred_icr_training[i]
  LRBS[i] = training_set$RT[i]-training_set$pred_LRBS_training[i]
  pls[i] =  training_set$RT[i]-training_set$pred_pls_training[i]
  MARS[i] = training_set$RT[i]-training_set$pred_MARS_training[i]
  knn[i] = training_set$RT[i]-training_set$pred_knn_training[i]
  RF[i] = training_set$RT[i]-training_set$pred_RF_training[i]
}
training_set<- training_set %>%
  mutate(ICR , LRBS, pls ,MARS ,knn, RF )
# p <- ggplot(test_set) + 
#   geom_boxplot(aes(x="ICR", y=ICR))+ 
#   geom_boxplot(aes(x="LRBS", y=LRBS))+
#   geom_boxplot(aes(x="PLS", y=pls))+
#   geom_boxplot(aes(x="MARS", y=MARS))+
#   geom_boxplot(aes(x="KNN", y=knn))+
#   geom_boxplot(aes(x="RF", y=RF))+
#   ylim(-10,10)+
#   xlab("Models")+
#   ylab("Error (min)")
# 
# p+my_theme
# 
# ggsave("box plot for test set - RP mode .svg",
#        width=10,
#        height=10,
#        units="cm")
# p <- ggplot(training_set) + 
#   geom_boxplot(aes(x="ICR", y=ICR))+ 
#   geom_boxplot(aes(x="LRBS", y=LRBS))+
#   geom_boxplot(aes(x="PLS", y=pls))+
#   geom_boxplot(aes(x="MARS", y=MARS))+
#   geom_boxplot(aes(x="KNN", y=knn))+
#   geom_boxplot(aes(x="RF", y=RF))+
#   xlab("Models")+
#   ylab("Error (min)")
# 
# p+my_theme
# ggsave("box plot for training set - RP .svg",
#        width=10,
#        height=10,
#        units="cm")
RF_training <- training_set %>% dplyr::select(-c(ICR,LRBS, pls, MARS, knn))
names(RF_training)[names(RF_training) == "RF"] <- "Error"
RF_training <- RF_training %>% mutate(Model="RF")
ICR <- training_set %>% dplyr::select(-c(RF,LRBS, pls, MARS, knn))
names(ICR)[names(ICR) == "ICR"] <- "Error"
ICR <- ICR %>% mutate(Model="ICR")
LRBS <- training_set %>% dplyr::select(-c(RF,ICR, pls, MARS, knn))
names(LRBS)[names(LRBS) == "LRBS"] <- "Error"
LRBS <- LRBS %>% mutate(Model="LRBS")
pls <- training_set %>% dplyr::select(-c(RF,ICR,LRBS, MARS, knn))
names(pls)[names(pls) == "pls"] <- "Error"
pls <- pls %>% mutate(Model="PLS")
MARS <- training_set %>% dplyr::select(-c(RF,ICR,LRBS, pls, knn))
names(MARS)[names(MARS) == "MARS"] <- "Error"
MARS <- MARS %>% mutate(Model="MARS")
knn <- training_set %>% dplyr::select(-c(RF,ICR,LRBS, pls,MARS))
names(knn)[names(knn) == "knn"] <- "Error"
knn <- knn %>% mutate(Model="KNN")

RF_training <- RF_training %>%
  bind_rows(ICR) %>% bind_rows(LRBS) %>% 
  bind_rows(pls) %>% bind_rows(MARS) %>% bind_rows(knn)

p4 <- ggplot(data = RF_training,
            aes(x = Model, y = Error)) +
  geom_boxplot()+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               binwidth = 0.04,
               #fill = "grey",
               alpha = 0.5
  )+ ylab("Error (min)")+xlab("Models")
p4+my_theme
# ggsave("dotplot x boxplot for training set - RP .svg",
#        width=10,
#        height=10,
#        units="cm")
RF_test <- test_set %>% dplyr::select(-c(ICR,LRBS, pls, MARS, knn))
names(RF_test)[names(RF_test) == "RF"] <- "Error"
RF_test <- RF_test %>% mutate(Model="RF")
ICR <- test_set %>% dplyr::select(-c(RF,LRBS, pls, MARS, knn))
names(ICR)[names(ICR) == "ICR"] <- "Error"
ICR <- ICR %>% mutate(Model="ICR")
LRBS <- test_set %>% dplyr::select(-c(RF,ICR, pls, MARS, knn))
names(LRBS)[names(LRBS) == "LRBS"] <- "Error"
LRBS <- LRBS %>% mutate(Model="LRBS")
pls <- test_set %>% dplyr::select(-c(RF,ICR,LRBS, MARS, knn))
names(pls)[names(pls) == "pls"] <- "Error"
pls <- pls %>% mutate(Model="PLS")
MARS <- test_set %>% dplyr::select(-c(RF,ICR,LRBS, pls, knn))
names(MARS)[names(MARS) == "MARS"] <- "Error"
MARS <- MARS %>% mutate(Model="MARS")
knn <- test_set %>% dplyr::select(-c(RF,ICR,LRBS, pls,MARS))
names(knn)[names(knn) == "knn"] <- "Error"
knn <- knn %>% mutate(Model="KNN")

RF_test<- RF_test %>%
  bind_rows(ICR) %>% bind_rows(LRBS) %>% 
  bind_rows(pls) %>% bind_rows(MARS) %>% bind_rows(knn)

p5 <- ggplot(data = RF_test,
            aes(x = Model, y = Error)) +
  geom_boxplot()+
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               binwidth = 0.04, #give it different values 
               #fill = "grey",
               alpha = 0.5
  )+ ylim(-10,10) +ylab("Error (min)")+xlab("Models")
p5+my_theme

ggsave("dotplot x boxplot for test set - RP.svg",
       width=10,
       height=10,
       units="cm")

pred_RT = predict(regressor_RF, newdata = dataset)

dataset <- dataset %>% mutate(pred_RT)
dataset  <- data.frame(dataset, error=dataset$RT-dataset$pred_RT)
mean <- mean(dataset$error)
median <- median(dataset$error)
dataset  <- data.frame(dataset, relative_error=abs(dataset$RT-dataset$pred_RT)/dataset$RT)
relative_mean <- mean(dataset$relative_error)
relative_median <- median(dataset$relative_error)

varImp(regressor_RF)
my_theme1 <- theme(panel.grid.major =  element_blank(),
                   panel.grid.minor=element_blank())
#get the excel file which contains the important features

p6 <- ggplot(data=Important_factors_in_RF)+
  geom_col(mapping=aes(x=Feature, y= Importance))+
  xlab("Variables")+
  ylab("Importance in the RF model")+
  theme_bw()

p6 + my_theme
ggsave("feature importance RF in RP mode.svg",
       width=10,
       height=10,
       units="cm")
saveRDS(regressor_RF,  "regressorRF for RP16062021.rds")
library(cowplot)
plot_grid(p3+my_theme, p6+my_theme,
          labels =c("a) Predictions for reversed-C18", "b) Important variables in RF model"),
          align = "h",
          hjust = 0.1, vjust = 0.8,
          label_x = 0.3,
          label_y = 0.99,
          nrow = 1,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("prediction of RF in RF.svg",
       width = 30,
       height = 22.5/668*380,
       units = "cm")

