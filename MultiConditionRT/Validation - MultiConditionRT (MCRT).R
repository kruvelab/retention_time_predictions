#-------------------------- Libraries -----------------------------------------------
library(caret) #machine learinign workflow
library(leaps) #library for stepwise regression
library(MASS) #contains some important linear regression tools
library(caTools) #sample split is from this package
library(tidyverse) #helps us to write concise code
library(rcdk)
library(Metrics)
library(readxl)
setwd("C:/Users/amso1873/OneDrive - Kruvelab/Amina/R code")
Eawag <- read_excel("Eawag casnr.xlsx")
Miklos <- read_excel("Miklos casnr.xlsx")

Eawag <- Eawag %>%
  group_by(name, SMILES,mz,  IC, MW, casrn, Sample_type) %>%
  summarise(RT= mean(RT)) 

dataset2 <-  read_delim('Filter_RP_2010.csv',
                        delim = ",",
                        col_names = TRUE,
                        trim_ws = TRUE) 
dataset2 <- dataset2 %>%
  group_by(Compound_name, Column, Organic_modifier, pH, Buffer, SMILES ) %>%
  summarise(slope = mean(slope),
            RT_Miklos= mean(ret_time))  
dataset2 <- dataset2 %>% filter(pH==2.7)
dataset2 <- dataset2 %>% filter(Organic_modifier=="Methanol")
dataset2 <- dataset2 %>% filter(Buffer=="Formic acid")
dataset2 <- dataset2 %>% dplyr::select(-SMILES)
dataset2 <- dataset2 %>% left_join(Miklos)
dataset2 <- dataset2 %>% dplyr::select(-SMILES)
dataset2 <- dataset2 %>% left_join(Eawag)
dataset2 <- dataset2 %>% drop_na()
SelectionCalibrationCompounds <- read_excel("SelectionCalibrationCompounds.xlsx", 
                                            sheet = "Tabelle1")
names(SelectionCalibrationCompounds)[names(SelectionCalibrationCompounds) == "CAS_Nr"] <- "casrn"
SelectionCalibrationCompounds <- SelectionCalibrationCompounds %>% dplyr::select(-SMILES)
dataset2 <- dataset2 %>% left_join(SelectionCalibrationCompounds)
names(dataset2)[names(dataset2) == "JC_logD_pH3"] <- "logD at pH=3"
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

p1<- ggplot(data=dataset2) +
  geom_point(mapping = aes(x = RT, y=RT_Miklos, color=`logD at pH=3`), size = 3, alpha = 1) +
  xlab(italic("t")["R"]~"(min) external dataset")+
  ylab(italic("t")["R"]~"(min) internal dataset")+
  xlim(1,26)+
  ylim(1,26)+
  labs(color= substitute("log"~italic("D")["at pH=3"]))+
  theme_bw() 
p1+my_theme

# 27 common compounds from casnr
library(mgcv)
regressor1 <- gam(RT~RT_Miklos, family = gaussian(), data=dataset2, method="GCV.Cp", optimizer = c("outer", "newton"))
#----------------Predicting for the dataset----------

dataset <-  read_delim('Results_standards_Eawag.csv',
                       delim = ",",
                       col_names = TRUE,
                       trim_ws = TRUE) 
smiles <- dataset$SMILES

descs <- dataset[29:1245]
descs <- descs %>% mutate(smiles)
descs <- unique(descs)
names(descs)[names(descs) == "smiles"] <- "SMILES"
dataset <- dataset %>%
  group_by(Compound, organic_modifier, pH.aq. , mz, SMILES) %>%
  summarise(RT= mean(RT))
dataset <- dataset %>% left_join(descs)
names(dataset)[names(dataset) == "pH.aq."] <- "pH"
names(dataset)[names(dataset) == "organic_modifier"] <- "Organic_modifier"
dataset$Organic_modifier=factor(dataset$Organic_modifier,
                                levels = c('Acetonitrile','MeOH'),
                                labels=c(1,2))
dataset <- data.frame(dataset, Buffer="Formic acid")
dataset$Buffer=factor(dataset$Buffer,
                      levels = c('Formic acid'),
                      labels=c(1))
SelectionCalibrationCompounds <- read_excel("SelectionCalibrationCompounds.xlsx", 
                                            sheet = "Tabelle1")
logD <- data.frame(SMILES=SelectionCalibrationCompounds$SMILES )
logD <- logD %>% mutate(logD=SelectionCalibrationCompounds$JC_logD_pH3)
dataset <- dataset %>% left_join(logD)
# dataset <- dataset %>% filter(logD<6)
# dataset <- dataset %>% filter(logD>-3)

#------Reading the RF model -------
regressor_RF <- read_rds("regressorRF for RP 16082021.rds")
names(dataset)[names(dataset) == "logD"] <- "log_P_sum"
dataset[is.na(dataset)] <- 0
pred_RF= predict(regressor_RF, newdata = dataset) 
dataset <- dataset %>% mutate(pred_RF)

p2<- ggplot(data=dataset) +
  geom_point(mapping = aes(x = RT, y =pred_RF, color=log_P_sum), size = 3, alpha = 1) +
  xlab("Measured"~t[R]~"(min)")+
  ylab("Predicted"~t[R]~"(min)")+
  xlim(0,25)+
  ylim(0,25)+
  labs(color= substitute("log"~italic("D")["at pH=3"]))+
  theme_bw()

p2+my_theme
 
RSME_RF<- RMSE(dataset$pred_RF, dataset$RT)

names(dataset)[names(dataset) == "pred_RF"] <- "RT_Miklos"
Trans_pred_RF= predict(regressor1, newdata = dataset) 
dataset <- dataset %>% mutate(Trans_pred_RF)

p3<- ggplot(data=dataset) +
  geom_point(mapping = aes(x = RT, y =Trans_pred_RF, color=log_P_sum), size = 3, alpha = 1) +
  xlab("Measured"~t[R]~"(min)")+
  ylab("Predicted"~t[R]~"(min)")+
  xlim(0,25)+
  ylim(0,25)+
  labs(color= substitute("log"~italic("D")["at pH=3"]))+
  theme_bw()

p3+my_theme
RSME_RF2<- RMSE(dataset$Trans_pred_RF, dataset$RT)
dataset <- data.frame(dataset, error=abs(dataset$Trans_pred_RF-dataset$RT))
error_mean <- mean(abs(dataset$error))
error_median <- median(abs(dataset$error))
relative_mean <-  mean(abs(dataset$error/dataset$RT))*100
# 14.74%
relative_median <- median(abs(dataset$error/dataset$RT))*100
#6.45

dataset <-data.frame(dataset, error=abs(dataset$RT-dataset$Trans_pred_RF))
dataset4 <- dataset %>% filter( log_P_sum>0)
RSME_RF3<- RMSE(dataset4$Trans_pred_RF, dataset4$RT)
dataset4 <- data.frame(dataset4, error=abs(dataset4$Trans_pred_RF-dataset4$RT))
error_mean1 <- mean(abs(dataset4$error))
error_median1 <- median(abs(dataset4$error))
relative_error_mean1 <- mean(abs(dataset4$error)/dataset4$RT)
relative_error_median1 <- median(abs(dataset4$error/dataset4$RT))

dataset <- dataset %>% filter(error<1)

print(regressor1)
summary(regressor1)
#----------------------------------------Emma predictions--------------
Emma_list <- read_excel("Emma list.xlsx", col_names = TRUE)
Emma_list2 <- read_excel("logP.xlsx", col_names = TRUE) 

Emma_list <- Emma_list%>%
  left_join(Emma_list2)

for (i in 1:528) {
  if (is.na(Emma_list$`pKa base`[i])){
    Emma_list$`pKa base`[i]<--5
  }
}
for (i in 1:528) {
  if (is.na(Emma_list$`pKa acid`[i])){
    Emma_list$`pKa acid`[i]<-20
  }
}

Emma_list <- Emma_list %>% dplyr::select(-log_P_sum)

names(Emma_list)[names(Emma_list) == "logP AH"] <- "log_P_sum"
Emma_list$Organic_modifier=factor(Emma_list$Organic_modifier,
                                  levels = c(1),
                                  labels=c(1))
Emma_list$Buffer=factor(Emma_list$Buffer,
                        levels = c(1,6,7),
                        labels=c(1,6,7))
Emma <- read_excel("Emma casnr.xlsx")
Miklos <- read_excel("Miklos casnr.xlsx")
Emma <-unique(Emma[1:2])
Miklos <- Miklos %>% dplyr::select(-SMILES)
dataset3 <- read_excel("Miklos exp.xlsx")
dataset3$Organic_modifier=factor(dataset3$Organic_modifier,
                                levels = c(1),
                                labels=c(1))
dataset3$Buffer=factor(dataset3$Buffer,
                      levels = c(1,6,7),
                      labels=c(1,6,7))

dataset3 <- dataset3 %>% left_join(Miklos)

Emma_list <- Emma_list %>% left_join(Emma)
dataset3 <- dataset3 %>% dplyr::select(-Compound_name)

dataset3 <- dataset3 %>% left_join(Emma_list)
dataset3 <- dataset3 %>% drop_na()

names(dataset3)[names(dataset3) == "log_P_sum"] <- "LogP.AH"
p4 <- ggplot(data=dataset3)+
  geom_point(mapping=aes(x=RT, y=RT_Miklos, color=LogP.AH), size=3, alpha=1)+
  xlab(t[R]~"(min) external dataset")+
  ylab(t[R]~"(min) internal dataset")+
  xlim(0,25)+
  ylim(0,25)+
  labs(color= substitute("log"~italic("P")))+
  theme_bw()
p4+ my_theme 

dataset3 <- dataset3 %>% drop_na()
#173 common data points(58 common compound)
# ggsave("RT_EMMA VS RT_Miklos FOR EMMA colored based on logP_sum.svg",
#        width=10,
#        height=10,
#        units="cm")
library(mgcv)
regressor <- gam(RT~RT_Miklos, family = gaussian(), data=dataset3, method="GCV.Cp", optimizer = c("outer", "newton") )
descs <- read_delim('Emma_smiles_PaDEL_descs.csv',
                    delim = ",",
                    col_names = TRUE,
                    trim_ws = TRUE) 
names(descs)[names(descs) == "smiles"] <- "SMILES"
Emma_list <- Emma_list %>% left_join(descs)
setwd("C:/Users/amso1873/OneDrive - Kruvelab/Amina/R code")
regressor_RF <- read_rds("regressorRF for RP 16082021.rds")
names(Emma_list)[names(Emma_list) == "ASP-0"] <- "ASP.0"
names(Emma_list)[names(Emma_list) == "AVP-6"] <- "AVP.6"
# names(Emma_list)[names(Emma_list) == "ASP.2"] <- "ASP-2"
# names(Emma_list)[names(Emma_list) == "ASP.3"] <- "ASP-3"
# names(Emma_list)[names(Emma_list) == "ASP.4"] <- "ASP-4"
# names(Emma_list)[names(Emma_list) == "ASP.5"] <- "ASP-5"
names(Emma_list)[names(Emma_list) == "ASP-6"] <- "ASP.6"
names(Emma_list)[names(Emma_list) == "ASP-7"] <- "ASP.7"
names(Emma_list)[names(Emma_list) == "AVP-7"] <- "AVP.7"
names(Emma_list)[names(Emma_list) == "AVP-2"] <- "AVP.2"
names(Emma_list)[names(Emma_list) == "VCH-5"] <- "VCH.5"
names(Emma_list)[names(Emma_list) == "BCUTc-1h"] <- "BCUTc.1h"
names(Emma_list)[names(Emma_list) == "BCUTc-1l"] <- "BCUTc.1l"
names(Emma_list)[names(Emma_list) == "BCUTw-1h"] <- "BCUTw.1h"
names(Emma_list)[names(Emma_list) == "BCUTp-1l"] <- "BCUTp.1l"
names(Emma_list)[names(Emma_list) == "BCUTp-1h"] <- "BCUTp.1h"
names(Emma_list)[names(Emma_list) == "BCUTw-1l"] <- "BCUTw.1l"
names(Emma_list)[names(Emma_list) == "SCH-5"] <- "SCH.5"
names(Emma_list)[names(Emma_list) == "WTPT-4"] <- "WTPT.4"
names(Emma_list)[names(Emma_list) == "WTPT-5"] <- "WTPT.5"
names(Emma_list)[names(Emma_list) == "logP.AH"] <- "LogP_sum"
# names(Emma_list)[names(Emma_list) == "Organic_modifier"] <- "organic_modifier"
Emma_list <- Emma_list %>%
  drop_na()
RT_Miklos= predict(regressor_RF, newdata = Emma_list) 

Emma_list <- Emma_list %>% mutate(RT_Miklos)


Tran_pred_RF= predict(regressor, newdata = Emma_list) 
Emma_list <- Emma_list %>% mutate(Tran_pred_RF)
names(Emma_list)[names(Emma_list) == "log_P_sum"] <- "LogP.AH"
p5 <- ggplot(data=Emma_list)+
  geom_point(mapping=aes(x=RT, y=RT_Miklos, color=LogP.AH), size=3, alpha=3/4)+
  xlab("Actual "~t[R]~"(min)")+
  ylab("Predicted "~t[R]~"(min)")+
  xlim(0,25)+
  ylim(0,25)+
  theme_bw()+
  labs(color= substitute("log"~italic("P")))
p5+ my_theme
# Emma_list <- Emma_list %>% mutate(Tran_pred_RF)
Emma_list3 <- Emma_list %>% filter(LogP.AH>-3)
Emma_list3 <- Emma_list3 %>% filter(LogP.AH<5)

p6 <- ggplot(data=Emma_list)+
  geom_point(mapping=aes(x=RT, y=Tran_pred_RF, color=LogP.AH), size=3, alpha=3/4)+
  xlab("Measured"~t[R]~"(min)")+
  ylab("Predicted"~t[R]~"(min)")+
  xlim(0,25)+
  ylim(0,25)+
  labs(color= substitute("log"~italic("P")))+
  theme_bw()
p6+ my_theme
#NO EFFECT 
p7 <- ggplot(data=Emma_list3)+
  geom_point(mapping=aes(x=RT, y=RT_Miklos, color=LogP.AH), size=3, alpha=1)+
  xlab("Actual "~t[R]~"(min)")+
  ylab("Predicted "~t[R]~"(min)")+
  xlim(0,25)+
  ylim(0,25)+
  theme_bw()
p7+ my_theme
# write_delim(Emma_list,
#             'Emma predictions.csv',
#             delim = ",")
RMSE <- RMSE(Emma_list$RT, Emma_list$RT_Miklos)

RMSE_2 <- RMSE(Emma_list$RT, Emma_list$Tran_pred_RF)

#164 compounds 

#-------Isomers part--------
setwd("C:/Users/amso1873/OneDrive - Kruvelab/Amina/R code")
desc <-  read_delim('descs.csv',
                    delim = ",",
                    col_names = TRUE,
                    trim_ws = TRUE) 

Isomers <- read_excel("Isomers.xlsx", sheet = "SMILES ")
Isomers <- Isomers %>% left_join(desc)
names(Isomers)[names(Isomers) == "logKow"] <- "log_P_sum"
regressor_RF <- read_rds("regressorRF for RP16062021.rds")
Isomers <- data.frame(Isomers, pH =2.7)
Isomers <- data.frame(Isomers, Buffer=1)
Isomers$Buffer <- as.factor(Isomers$Buffer)
Isomers <- data.frame(Isomers, Organic_modifier=1)
Isomers$Organic_modifier <- as.factor(Isomers$Organic_modifier)
pred_RF= predict(regressor_RF, newdata = Isomers)
Isomers <- Isomers %>% mutate(pred_RF)
Isomers <- data.frame(Isomers, max=pred_RF+1.24)
Isomers <- data.frame(Isomers, min=pred_RF-1.24)
library(extrafont)
font_import()
font <- choose_font("Arial")
fontsize <- 12
basecolor <- "black"
my_theme2 <- theme(plot.background = element_blank(), 
                   panel.background = element_blank(), axis.line=element_line(  size=0.5), 
                   legend.background = element_blank(),
                   legend.title = element_text(),
                   legend.position = c(0.9,0.3), 
                   axis.text.x = element_text(angle =0, vjust = 0.5, hjust=0.5),
                   aspect.ratio = 1, text = element_text(family = font, size=fontsize),
                   plot.title = element_text(color=basecolor, size=14, face="bold"),
                   legend.text=element_text(family = font, size=fontsize, color=basecolor),
                   axis.text=element_text(family = font, size=fontsize, color=basecolor), 
                   panel.grid.major =  element_blank(),
                   panel.grid.minor=element_blank(), 
                   panel.border = element_blank())
p8 <- ggplot(data=Isomers) +
  geom_point(mapping = aes(x=Compound.name, y=pred_RF)) +
  geom_errorbar(aes(x=Compound.name, ymin=min,ymax=max))+
  xlab("Compound name")+
  ylab("Predicted" ~t[R]~" (min)")+coord_flip()
p8+my_theme2
ggsave("Isomers predictions .svg",
       width=20,
       height=10,
       units="cm")



library(cowplot)
plot_grid(p4+my_theme,  p6+my_theme , p1+my_theme,  p2+my_theme,
          labels =c("a) ", "b)  ","c)  ",
                    "d)  ","e)"),
          align = "h",
          hjust =  0.1, vjust = 1.5 ,
          label_x = 0.3,
          label_y = 0.999,
          nrow = 3,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("Figure 4 in manuscript.svg", height = 9, width=9)
Emma_list$Compound_name


plot_grid(p8+my_theme2,  p9+my_theme2 ,
          labels =c("a) ", "b)  "),
          align = "h",
          hjust =  0.1, vjust = 1.5 ,
          label_x = 0.3,
          label_y = 0.999,
          nrow = 1,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("Figure 6 in manuscript.svg", height = 5, width=13)

Isomers$Isomer.number <- as.factor(Isomers$Isomer.number)

Isomers1 <- Isomers %>% filter(set==1)
p8 <- ggplot(data=Isomers1) +
  geom_point(mapping = aes(x=Isomer.number, y=pred_RF)) +
  geom_errorbar(aes(x=Isomer.number, ymin=min,ymax=max))+
  xlab("Isomers ")+
  ylab("Predicted" ~italic("t")["R"]~"(min)")+coord_flip()
p8+my_theme2
Isomers3 <- Isomers %>% filter(set==2)
p9 <- ggplot(data=Isomers3) +
  geom_point(mapping = aes(x=Isomer.number, y=pred_RF)) +
  geom_errorbar(aes(x=Isomer.number, ymin=min,ymax=max))+
  xlab("Compound name")+
  ylab("Predicted" ~italic("t")["R"]~"(min)")+coord_flip()
p9+my_theme2
Isomers4 <- Isomers %>% filter(set==3)
p10 <- ggplot(data=Isomers4) +
  geom_point(mapping = aes(x=Isomer.number, y=pred_RF)) +
  geom_errorbar(aes(x=Isomer.number, ymin=min,ymax=max))+
  xlab("Compound name")+
  ylab("Predicted" ~italic("t")["R"]~"(min)")+coord_flip()
p10+my_theme2
Isomers5 <- Isomers %>% filter(set==4)
p11 <- ggplot(data=Isomers5) +
  geom_point(mapping = aes(x=Isomer.number, y=pred_RF)) +
  geom_errorbar(aes(x=Isomer.number, ymin=min,ymax=max))+
  xlab("Compound name")+
  ylab("Predicted" ~italic("t")["R"]~" (min)")+coord_flip()
p11+my_theme2
plot_grid(p4+my_theme,  p6+my_theme , p1+my_theme,  p2+my_theme,p8+my_theme2,  p9+my_theme2 , p10+my_theme2, p11+my_theme2,
          labels =c("a) ", "b)  ", "c)", "d)", "e) C"~["10"]~"H"~["21"]~"N"~["1"]~"O"~["1"], "f)C16H18N2O2 ", "g)C9H12O4", "h)C8H10O4"),
          align = "h",
          hjust =  0.1, vjust = 1.5 ,
          label_x = 0.3,
          label_y = 0.999,
          nrow = 2,
          label_size = 12,
          label_fontfamily = font,
          label_colour = basecolor,
          label_fontface = "plain")
ggsave("Figure 6 in manuscript 4 set  .svg", height = 8, width=15)
