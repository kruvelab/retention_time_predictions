#-------------------------- Libraries -----------------------------------------------
library(caret) #machine learinign workflow
library(leaps) #library for stepwise regression
library(MASS) #contains some important linear regression tools
library(caTools) #sample split is from this package
library(tidyverse) #helps us to write concise code
library(rcdk)
library(Metrics)
#------------------------importing the data---------------------------------------------- 
#run the RP code first 
dataset <-  read_delim('descs200502.csv',
                       delim = ",",
                       col_names = TRUE,
                       trim_ws = TRUE) %>%
  drop_na()
dataset2 <-  read_delim('Filtered Mixed data1510.csv',
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
#Replacing the NAs with 0
dataset[is.na(dataset)] <- 0

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
pkaan[is.na(pkaan)] <--5
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
dataset <- dataset %>%
  drop_na()


p <- ggplot() +
  geom_point(mapping = aes(x = dataset$logP.AH, y =dataset$log_P_sum), color = "blue", size = 2, alpha = 0.25) +
  xlab("LogP_AH ")+
  ylab("logP_sum")+
  theme_bw()
p

# ggsave("logP_AH VS logP_sum.svg",
#        width=10,
#        height=10,
#        units="cm")

dataset <- dataset %>%   left_join(SMILES_unique) %>%  
  dplyr::select(SPLIT, everything()) %>%   na.omit() 
#trainingset 
training_set <- dataset %>%   filter(SPLIT == TRUE) %>%  
  dplyr::select(-SPLIT)  
#testset 
test_set <- dataset %>%   filter(SPLIT == FALSE) %>%
  dplyr::select(-SPLIT) 


Comp_training <- training_set$Compound_name
Comp_test <- test_set$Compound_name
intersect(Comp_training, Comp_test)

training_set <- training_set %>%
  dplyr::select(-c(SMILES, Compound_name))
test_set <- test_set %>%
  dplyr::select(-c(SMILES, Compound_name))



#setting the cross-validation parameters
fitControl <- trainControl(method = "repeatedcv", number = 5, repeats = 5) 

#-------------------------LM-------------------------
set.seed(825)
regressor_lm <- train(RT ~ log_P_sum, 
                      data = training_set, 
                      method = "lm", 
                      trControl = fitControl)

print(regressor_lm)
summary(regressor_lm)

dataset <- dataset %>% 
  mutate(tR_pred = predict(regressor_lm, newdata = dataset))
dataset <- dataset %>%
  mutate(SPLIT)
p <- ggplot(data=dataset)+
  geom_point(mapping=aes(x=RT, y= tR_pred, color = SPLIT), size=3, alpha=3/4)+
  scale_color_manual(values=c("#d90429", "#023e8a")) +
  xlab("Actual" ~t[R]~ "(min)")+
  ylab("Predicted " ~t[R]~ "(min)")+ 
  labs("Predicted VS Actual RT")
p+my_theme
# ggsave("linear regression with calculated logP_sum.svg",
#        width=10,
#        height=10,
#        units="cm")
#-------- running all the models --------
LogP_sum <- dataset$log_P_sum
dataset <- dataset[1:206]
dataset <- dataset %>% select(-SPLIT)
dataset <- dataset %>% mutate(LogP_sum)
library(readxl)
dataset2 <- read_excel("mass values.xlsx")
dataset <- dataset %>% left_join(dataset2)


#join the split to the data 
dataset <- dataset %>%   left_join(SMILES_unique) %>%  
  dplyr::select(SPLIT, everything()) %>%   na.omit() 
#trainingset 
training_set <- dataset %>%   filter(SPLIT == TRUE) %>%  
  dplyr::select(-SPLIT)  
#testset 
test_set <- dataset %>%   filter(SPLIT == FALSE) %>%
  dplyr::select(-SPLIT) 


Comp_training <- training_set$Compound_name
Comp_test <- test_set$Compound_name
intersect(Comp_training, Comp_test)
folds <- groupKFold(training_set$SMILES, k = 20)
training_set <- training_set %>%
  dplyr::select(-c(SMILES, Compound_name))
test_set <- test_set %>%
  dplyr::select(-c(SMILES, Compound_name))

#setting the cross-validation parameters
# fitControl <- trainControl(method = "repeatedcv", number = 5, repeats = 5) 
# 
# #Trainign the model
fitControl <- trainControl(## 10-fold CV
  method = "boot",
  index = folds)
#-----------------------Training and cross-validation -----------------------------
#-------------------------Independent Component Regression----------------------
#install.packages("fastICA")
library(fastICA)
set.seed(825)
regressor_icr <- train(RT ~ ., 
                       data = training_set, 
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
                        data = training_set, 
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
                       data = training_set, 
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
                        data = training_set, 
                        method = "earth", 
                        trControl = fitControl)

print(regressor_MARS)
summary(regressor_MARS)

#-------------------------KNN-------------------------
#install.packages("kknn")
library(kknn)
set.seed(825)
regressor_knn <- train(RT ~ ., 
                       data = training_set, 
                       method = "kknn", 
                       trControl = fitControl)

print(regressor_knn)
summary(regressor_knn)
#-------------------------RF-------------------------
set.seed(825)
regressor_RF <- train(RT ~ ., 
                      data = training_set, 
                      method = "ranger", 
                      trControl = fitControl)

print(regressor_RF)
summary(regressor_RF)
saveRDS(regressor_RF,  "regressorRF for mixed mode.rds")
#-------------------------Validation--------------------------------
#Predicting for the test set
pred_icr_test = predict(regressor_icr, newdata = test_set)
pred_LRBS_test = predict(regressor_LRBS, newdata = test_set)
pred_pls_test = predict(regressor_pls, newdata = test_set)
pred_MARS_test = predict(regressor_MARS, newdata = test_set) 
pred_knn_test = predict(regressor_knn, newdata = test_set) 
pred_RF_test = predict(regressor_RF, newdata = test_set) 

#Predicting tR for the training set 
pred_icr_training = predict(regressor_icr, newdata = training_set)
pred_LRBS_training = predict(regressor_LRBS, newdata = training_set)
pred_pls_training = predict(regressor_pls, newdata = training_set)
pred_MARS_training= predict(regressor_MARS, newdata = training_set) 
pred_knn_training= predict(regressor_knn, newdata = training_set) 
pred_RF_training = predict(regressor_RF, newdata = training_set) 

test_set<- test_set %>%
  mutate(pred_icr_test, pred_LRBS_test, pred_pls_test,  pred_MARS_test,  pred_knn_test, pred_RF_test)

training_set<- training_set %>%
  mutate(pred_icr_training, pred_LRBS_training, pred_pls_training,  pred_MARS_training,  pred_knn_training, pred_RF_training)


#Calculation the relative error 
#Predicting tR for the training set 
ICR <- c()
LRBS <- c()
pls <- c()
MARS <- c()
knn <- c()
RF <- c()
for( i in 1:113) {
  ICR[i] = test_set$RT[i]-test_set$pred_icr_test[i]
  LRBS[i] = test_set$RT[i]-test_set$pred_LRBS_test[i]
  pls[i] = test_set$RT[i]-test_set$pred_pls_test[i]
  MARS[i] = test_set$RT[i]-test_set$pred_MARS_test[i]
  knn[i] = test_set$RT[i]-test_set$pred_knn_test[i]
  RF[i] =  test_set$RT[i]-test_set$pred_RF_test[i]
}
test_set<- test_set %>%
  mutate(ICR , LRBS, pls ,MARS ,knn, RF )

for( i in 1:389) {
  ICR[i] = training_set$RT[i]-training_set$pred_icr_training[i]
  LRBS[i] = training_set$RT[i]-training_set$pred_LRBS_training[i]
  pls[i] =  training_set$RT[i]-training_set$pred_pls_training[i]
  MARS[i] = training_set$RT[i]-training_set$pred_MARS_training[i]
  knn[i] = training_set$RT[i]-training_set$pred_knn_training[i]
  RF[i] = training_set$RT[i]-training_set$pred_RF_training[i]
}
training_set<- training_set %>%
  mutate(ICR , LRBS, pls ,MARS ,knn, RF )


#calculating the RMSE for the test set 
RSME_icr_testset <- RMSE(test_set$pred_icr_test, test_set$RT)
RSME_LRBS_testset <- RMSE(test_set$pred_LRBS_test, test_set$RT)
RSME_pls_testset <- RMSE(test_set$pred_pls_test, test_set$RT)
RSME_MARS_testset <- RMSE(test_set$pred_MARS_test, test_set$RT)
RSME_knn_testset <- RMSE(test_set$pred_knn_test, test_set$RT)
RSME_RF_testset <- RMSE(test_set$pred_RF_test, test_set$RT)

RSME_testset <- c(RSME_icr_testset, RSME_LRBS_testset, RSME_pls_testset, RSME_MARS_testset,  RSME_knn_testset, RSME_RF_testset)

#calculating the RMSE for the training set 
RSME_icr_trainingset <- RMSE(training_set$pred_icr_training, training_set$RT)
RSME_LRBS_trainingset <- RMSE(training_set$pred_LRBS_training, training_set$RT)
RSME_pls_trainingset<- RMSE(training_set$pred_pls_training, training_set$RT)
RSME_MARS_trainingset <- RMSE(training_set$pred_MARS_training, training_set$RT)
RSME_knn_trainingset <- RMSE(training_set$pred_knn_training, training_set$RT)
RSME_RF_trainingset <- RMSE(training_set$pred_RF_training, training_set$RT)

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
                  legend.position = c(1,1), 
                  aspect.ratio = 1, text = element_text(family = font, size=fontsize),
                  plot.title = element_text(color=basecolor, size=14, face="bold"),
                  legend.text=element_text(family = font, size=fontsize, color=basecolor),
                  axis.text=element_text(family = font, size=fontsize, color=basecolor), 
                  panel.grid.major =  element_blank(),
                  panel.grid.minor=element_blank(), 
                  panel.border = element_blank())
p <- ggplot(test_set) + 
  geom_boxplot(aes(x="ICR", y=ICR))+ 
  geom_boxplot(aes(x="LRBS", y=LRBS))+
  geom_boxplot(aes(x="PLS", y=pls))+
  geom_boxplot(aes(x="MARS", y=MARS))+
  geom_boxplot(aes(x="KNN", y=knn))+
  geom_boxplot(aes(x="RF", y=RF))+
  xlab("Models")+
  ylab("Error (min)")

p+my_theme
ggsave("box plot for test set - Mixed mode .svg",
       width=10,
       height=10,
       units="cm")
p <- ggplot(training_set) + 
  geom_boxplot(aes(x="ICR", y=ICR))+ 
  geom_boxplot(aes(x="LRBS", y=LRBS))+
  geom_boxplot(aes(x="PLS", y=pls))+
  geom_boxplot(aes(x="MARS", y=MARS))+
  geom_boxplot(aes(x="KNN", y=knn))+
  geom_boxplot(aes(x="RF", y=RF))+
  xlab("Models")+
  ylab("Error (min)")

p+my_theme
ggsave("box plot for training set - Mixed mode .svg",
       width=10,
       height=10,
       units="cm")
p <-ggplot() +
  geom_point(mapping = aes(x = Models, y =RSME_trainingset), color =  "#3d5a80", size = 2, alpha = 1) +
  geom_point(mapping = aes(x = Models, y =RSME_testset),color ="#ee6c4d",  size = 2, alpha = 1) +
  ylab("RMSE (min)")+
  theme_bw()

p+my_theme
ggsave("Mixed mode with calculated logP_sum.svg",
       width=10,
       height=10,
       units="cm")
p <- ggplot() +
  geom_point(mapping = aes(x =training_set$RT, y =training_set$pred_RF_training), color =  "#3d5a80", size = 2, alpha = 1) +
  geom_point(mapping = aes(x = test_set$RT, y =test_set$pred_RF_test), color ="#ee6c4d",  size = 2, alpha = 1) +
  xlab("Predicted tR (min)")+
  ylab("Actual tR (min)")+ 
  theme_bw()
p+my_theme
ggsave("RF Mixed mode with calculated logP_sum.svg", 
       width=10,
       height=10, 
       units="cm")
print(RSME_RF_testset)
#2,27
print(RSME_RF_trainingset)
#0,42

Compound_name <- dataset$Compound_name
Smiles_dataset <- dataset$SMILES
SPLIT <- dataset$SPLIT
dataset <- dataset %>% 
  dplyr::select(-c(SPLIT, SMILES,  Compound_name))
dataset <- dataset %>% 
  mutate(tR_pred = predict(regressor_RF, newdata = dataset))

dataset <- dataset %>%
  mutate(Compound_name, Smiles_dataset, SPLIT)
dataset <- dataset %>% left_join(dataset4)
dataset <- data.frame(dataset, Error = dataset$RT - dataset$tR_pred)
p <- ggplot(data=dataset)+
  geom_point(mapping=aes(x=RT, y=tR_pred, color = LogP_sum), size=3, alpha=3/4)+
  xlab("RT")+
  ylab("PRED_RT")
p+ my_theme
ggsave("RT VS PRED RT COLOURED BASED ON logP_sum.svg",
       width=10,
       height=10,
       units="cm")
p <- ggplot(data=dataset)+
  geom_point(mapping=aes(x=logP.AH, y=LogP_sum, color =Error), size=3, alpha=3/4)+
  xlab("RT")+
  ylab("PRED_RT")
p+ my_theme
# ggsave("LogP.AH Vs LogP_sum.svg",
#        width=10,
#        height=10,
#        units="cm")
RF_training <- training_set %>% select(-c(ICR,LRBS, pls, MARS, knn))
names(RF_training)[names(RF_training) == "RF"] <- "Error"
RF_training <- RF_training %>% mutate(Model="RF")
ICR <- training_set %>% select(-c(RF,LRBS, pls, MARS, knn))
names(ICR)[names(ICR) == "ICR"] <- "Error"
ICR <- ICR %>% mutate(Model="ICR")
LRBS <- training_set %>% select(-c(RF,ICR, pls, MARS, knn))
names(LRBS)[names(LRBS) == "LRBS"] <- "Error"
LRBS <- LRBS %>% mutate(Model="LRBS")
pls <- training_set %>% select(-c(RF,ICR,LRBS, MARS, knn))
names(pls)[names(pls) == "pls"] <- "Error"
pls <- pls %>% mutate(Model="pls")
MARS <- training_set %>% select(-c(RF,ICR,LRBS, pls, knn))
names(MARS)[names(MARS) == "MARS"] <- "Error"
MARS <- MARS %>% mutate(Model="MARS")
knn <- training_set %>% select(-c(RF,ICR,LRBS, pls,MARS))
names(knn)[names(knn) == "knn"] <- "Error"
knn <- knn %>% mutate(Model="knn")

RF_training <- RF_training %>%
  bind_rows(ICR) %>% bind_rows(LRBS) %>% 
  bind_rows(pls) %>% bind_rows(MARS) %>% bind_rows(knn)

p <- ggplot(data = RF_training,
            aes(x = Model, y = Error)) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               binwidth = 0.04,
               #fill = "grey",
               alpha = 0.5
  )
p+my_theme
ggsave("dotplot for training set - Mixed mode .svg",
       width=10,
       height=10,
       units="cm")

RF_test <- test_set %>% select(-c(ICR,LRBS, pls, MARS, knn))
names(RF_test)[names(RF_test) == "RF"] <- "Error"
RF_test <- RF_test %>% mutate(Model="RF")
ICR <- test_set %>% select(-c(RF,LRBS, pls, MARS, knn))
names(ICR)[names(ICR) == "ICR"] <- "Error"
ICR <- ICR %>% mutate(Model="ICR")
LRBS <- test_set %>% select(-c(RF,ICR, pls, MARS, knn))
names(LRBS)[names(LRBS) == "LRBS"] <- "Error"
LRBS <- LRBS %>% mutate(Model="LRBS")
pls <- test_set %>% select(-c(RF,ICR,LRBS, MARS, knn))
names(pls)[names(pls) == "pls"] <- "Error"
pls <- pls %>% mutate(Model="pls")
MARS <- test_set %>% select(-c(RF,ICR,LRBS, pls, knn))
names(MARS)[names(MARS) == "MARS"] <- "Error"
MARS <- MARS %>% mutate(Model="MARS")
knn <- test_set %>% select(-c(RF,ICR,LRBS, pls,MARS))
names(knn)[names(knn) == "knn"] <- "Error"
knn <- knn %>% mutate(Model="knn")

RF_test<- RF_test %>%
  bind_rows(ICR) %>% bind_rows(LRBS) %>% 
  bind_rows(pls) %>% bind_rows(MARS) %>% bind_rows(knn)

p <- ggplot(data = RF_test,
            aes(x = Model, y = Error)) +
  geom_dotplot(binaxis = "y",
               stackdir = "center",
               binwidth = 0.04, #give it different values 
               #fill = "grey",
               alpha = 0.5
  )+ ylim(-10,10)
p+my_theme
ggplot(RF_test) + 
  geom_boxplot(aes(x=Model, y=Error))+
  geom_dotplot(aes(x=Model, y=Error), binaxis='y', stackdir='center')
ggsave("dotplot for test set - Mixed mode.svg",
       width=10,
       height=10,
       units="cm")













#--------------- calculating the RMSE for each LC condition -----------------
training2ACN <- training_set %>% filter(pH == 2.7 & Organic_modifier==1)
test2ACN <- test_set %>% filter(pH == 2.7 & Organic_modifier==1)
RMSE_train2ACN <- RMSE(training2ACN$pred_RF_training, training2ACN$RT)
RMSE_test2ACN <- RMSE(test2ACN$pred_RF_test, test2ACN$RT)


training2MEOH <- training_set %>% filter(pH == 2.7 & Organic_modifier==2)
test2MEOH<- test_set %>% filter(pH == 2.7 & Organic_modifier==2)
RMSE_train2MEOH <- RMSE(training2MEOH$pred_RF_training, training2MEOH$RT)
RMSE_test2MEOH <- RMSE(test2MEOH$pred_RF_test, test2MEOH$RT)


training5ACN <- training_set %>% filter(pH == 5 & Organic_modifier==1)
test5ACN <- test_set %>% filter(pH == 5 & Organic_modifier==1)
RMSE_train5ACN <- RMSE(training5ACN$pred_RF_training, training5ACN$RT)
RMSE_test5ACN <- RMSE(test5ACN$pred_RF_test, test5ACN$RT)


training5MEOH <- training_set %>% filter(pH == 5 & Organic_modifier==2)
test5MEOH<- test_set %>% filter(pH == 5 & Organic_modifier==2)
RMSE_train5MEOH <- RMSE(training5MEOH$pred_RF_training, training5MEOH$RT)
RMSE_test5MEOH <- RMSE(test5MEOH$pred_RF_test, test5MEOH$RT)

training8ACN <- training_set %>% filter(pH == 8 & Organic_modifier==1)
test8ACN <- test_set %>% filter(pH == 8 & Organic_modifier==1)
RMSE_train8ACN <- RMSE(training8ACN$pred_RF_training, training8ACN$RT)
RMSE_test8ACN <- RMSE(test8ACN$pred_RF_test, test8ACN$RT)


training8MEOH <- training_set %>% filter(pH == 8 & Organic_modifier==2)
test8MEOH<- test_set %>% filter(pH == 8 & Organic_modifier==2)
RMSE_train8MEOH <- RMSE(training8MEOH$pred_RF_training, training8MEOH$RT)
RMSE_test8MEOH <- RMSE(test8MEOH$pred_RF_test, test8MEOH$RT)


training10ACN <- training_set %>% filter(pH == 10 & Organic_modifier==1)
test10ACN <- test_set %>% filter(pH == 10 & Organic_modifier==1)
RMSE_train10ACN <- RMSE(training10ACN$pred_RF_training, training10ACN$RT)
RMSE_test10ACN <- RMSE(test10ACN$pred_RF_test, test10ACN$RT)


training10MEOH <- training_set %>% filter(pH == 10 & Organic_modifier==2)
test10MEOH<- test_set %>% filter(pH == 10 & Organic_modifier==2)
RMSE_train10MEOH <- RMSE(training10MEOH$pred_RF_training, training10MEOH$RT)
RMSE_test10MEOH <- RMSE(test10MEOH$pred_RF_test, test10MEOH$RT)


RMSE <- c("2.7ACN", "2.7MEOH", "5ACN",  "5MEOH",  "8ACN", "8MEOH", "10ACN", "10MEOH")
RMSE_training <- c(RMSE_train2ACN, RMSE_train2MEOH, RMSE_train5ACN, RMSE_train5MEOH, RMSE_train8ACN, RMSE_train8MEOH, RMSE_train10ACN, RMSE_train10MEOH)
RMSE_test <- c(RMSE_test2ACN, RMSE_test2MEOH, RMSE_test5ACN, RMSE_test5MEOH, RMSE_test8ACN, RMSE_test8MEOH, RMSE_test10ACN, RMSE_test10MEOH)

p <-ggplot() +
  geom_point(mapping = aes(x = RMSE, y =RMSE_training), color = "blue", size = 2, alpha = 0.25) +
  geom_point(mapping = aes(x = RMSE, y =RMSE_test),color = "red",  size = 2, alpha = 0.25) +
  ylab("RMSE Values")+
  theme_bw()
p+ my_theme
ggsave("RMSE for each LC.svg", 
       width=10,
       height=10, 
       units="cm")
