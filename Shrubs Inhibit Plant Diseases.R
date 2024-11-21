setwd("D:\\谭晓丹\\谭晓丹的实验\\MQ\\2023MQ\\统计分析")


library(tidyverse)
library(lme4)
library(lmerTest)
library(readxl)
library(permute)
library(lattice)
library(vegan)
library(ggpubr)
library(sciplot)
library(MuMIn)
library("eoffice")
library(ggplot2)
library(plyr)
library(multcomp)

mytheme <- theme_test()+theme(axis.title.x = element_text(size = 20, face = "bold", colour = "black"),
                              axis.title.y = element_text(size = 20, face = "bold", colour = "black"),
                              axis.text.x = element_text(size = 20, colour = "black",hjust = 0.5),
                              axis.text.y = element_text(size = 20, colour = "black"))+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 0.5, linetype="solid"))+
  theme(panel.grid = element_blank())+
  theme(plot.margin = unit(x = c(0.2,0.4,0.2,0),units = "cm"))+
  theme(plot.title = element_text(size = 25,          #字体大小
                                  hjust = 0.5,          #字体左右的位置
                                  angle = 0))

col1 <- c("#7e0021","#214463","#43CD80","#fca311")#自定义图形颜色

na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}          #把NA变成零

##  Read data   ####

Community <- read.csv("Community.csv")

Population0 <- read.table("Population.txt",sep = '\t',header = TRUE) %>% mutate(Logitdiseaseseverity=log((DiseaseSeverity+0.001)/(1-DiseaseSeverity-0.001)))


## Linear mixed-effects mdoel  #### 

## Effects of different treatment on community level pathogen load  ####
model1 <- lmer(PL~Treatment+(1|Block),data=Community) 
summary(model1)
anova(model1)

model2 <- lmer(BioticPL~Treatment+(1|Block),data=Community) 
shapiro.test(resid(model2))
summary(model2)
anova(model2)

model3 <- lmer(NecrotroPL~Treatment+(1|Block),data=Community) 
shapiro.test(resid(model3))
summary(model3)
anova(model3)

## Effects of different treatment on population level pathogen load   ####

library(glmmTMB)
library(tidyverse)

data <- Population0 %>% dplyr::select(Treatment,DiseaseSeverity,Block,HistoryPathogen) %>%
  na.omit()
model <- glmmTMB(DiseaseSeverity~Treatment+(1|Block/Plot)+(1|Species), ziformula = ~Treatment+(1|Block/Plot)+(1|Species),
                 data=data %>% 
                   unite(Plot,Treatment,Block,remove = F),family = beta_family())
summary(model)

model <- glmmTMB(DiseaseSeverity~Treatment+(1|Block/Plot)+(1|Species), ziformula = ~Treatment+(1|Block/Plot)+(1|Species),
                 data=data %>% 
                   unite(Plot,Treatment,Block,remove = F) %>% 
                   filter(HistoryPathogen != 'Biotrophic') ,family = beta_family())
summary(model)

model <- glmmTMB(DiseaseSeverity~Treatment+(1|Block/Plot)+(1|Species), ziformula = ~Treatment+(1|Block/Plot),
                 data=data %>% 
                   unite(Plot,Treatment,Block,remove = F) %>% 
                   filter(HistoryPathogen != 'Necrotrophic') %>% 
                   dplyr::select(Block,Plot,Treatment,DiseaseSeverity,HistoryPathogen,Species) %>% na.omit() ,family = beta_family())
summary(model)


library(emmeans)

emmeans(model, specs = pairwise ~ Treatment, type = 'response', adjust = 'none')

## Effect of different treatment on abiotic and biotic factors ####
model1 <- lmer(log(Penetration)~Treatment+(1|Block),data=Community) 
summary(model1)
anova(model1)

model2 <- lmer(Temperature_m~Treatment+(1|Block),data=Community) 
summary(model2)
anova(model2)

model3 <- lmer(Moisture_m~Treatment+(1|Block),data=Community) 
summary(model3)
anova(model3)

model4 <- lmer(log(Litterbiomass)~Treatment+(1|Block),data=Community) 
summary(model4)
anova(model4)

model5 <- lmer(SR~Treatment+(1|Block),data=Community) 
summary(model5)
anova(model5)

model6 <- lmer(PlotAbovegroundBiomass~Treatment+(1|Block),data=Community) 
summary(model6)
anova(model6)

model1 <- lmer(SoilWater~Treatment+(1|Block),data=Community) 
summary(model1)
anova(model1)

model2 <- lmer(log(SoilC)~Treatment+(1|Block),data=Community) 
summary(model2)
anova(model2)

model3 <- lmer(log(SoilN)~Treatment+(1|Block),data=Community) 
summary(model3)
anova(model3)

model4 <- lmer(SoilAP~Treatment+(1|Block),data=Community) 
summary(model4)
anova(model4)


## Effects of abiotic and biotic factors on pathogen load  ####

model1 <- lmer(PL~Penetration+(1|Treatment)+(1|Block),Community)
summary(model1)
anova(model1)

model2 <- lmer(log(PL)~Temperature_m+(1|Treatment)+(1|Block),Community)
summary(model2)
anova(model2)

model3 <- lmer(log(PL)~Moisture_m +(1|Treatment)+(1|Block),Community)
summary(model3)
anova(model3)

model4 <- lmer(log(PL)~Litterbiomass +(1|Treatment)+(1|Block),Community)
summary(model4)
anova(model4)

model5 <- lmer(log(PL)~SR +(1|Treatment)+(1|Block),Community)
summary(model5)
anova(model5)

model6 <- lmer(PL~PlotAbovegroundBiomass +(1|Treatment)+(1|Block),Community)
summary(model6)
anova(model6)


## Model selection #### 
model = lmer(PL~PlotAbovegroundBiomass +SoilWater+SoilpH+SoilN+SoilAP +SR+Penetration+Litterbiomass+
               Temperature_m +Moisture_m +(1|Block),
             data = Community,na.action = "na.fail") 
summary(model)
shapiro.test(resid(model))
car::vif(model)#看共线性
model <- arm::standardize(model,standardize.y = TRUE) #数据标准化
p_msb <- dredge(model,evaluate = T,extra = c("R^2", F = function(x)
  summary(x)$fstatistic[[1]])) 
summary(p_msb)

p_msb

ma2 <- model.avg(p_msb)
summary(ma2)  #report results: importance

## SEM  ####
model1 <- psem(lme(PL~PlotAbovegroundBiomass+Penetration+SoilN,random = ~1|Block,Community),
               lme(SoilN~FS+GS+RS,random = ~1|Block,Community),
               lme(Penetration~FS+GS+RS,random = ~1|Block,Community),
               lme(PlotAbovegroundBiomass~FS+GS+RS,random = ~1|Block,Community))
summary(model1)


## Partial regression 

SEMPL <- list(lm(PL~PlotAbovegroundBiomass+SoilN+Penetration,random =~1|Block,data = SEM),
              lm(PlotAbovegroundBiomass~FS+GS+RS,random =~1|Block,data = SEM),
              lm(SoilN~FS+GS+RS,random =~1|Block,data = SEM),
              lm(Penetration~FS+GS+RS,random =~1|Block,data = SEM))

summary(SEMPL)

resids1 <- partialResid(PL ~ PlotAbovegroundBiomass ,SEMPL)
resids1

resids2 <- partialResid(PL ~ SoilN ,SEMPL)
resids2

resids3 <- partialResid(PL ~ Penetration ,SEMPL)
resids3





















