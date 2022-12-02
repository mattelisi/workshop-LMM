## ----------------------------------------------------------

rm(list=ls())
setwd("/mnt/D/Dropbox/sync/RHUL/stats_role/WORKSHOPS/workshop_LMM")
setwd("~/Dropbox/sync/RHUL/stats_role/WORKSHOPS/workshop_LMM")
library(ggplot2)
library(tidyverse)
library(lme4)
library(lmerTest)
library(lubridate)

## ----------------------------------------------------------
## 1 random intercepts and random slopes



## ----------------------------------------------------------
cd4 <- read_csv("./datasets/cd4/cd4.csv")
str(cd4)

cd4$time <- mdy(cd4$VDATE)
cd4$time[10] - min(cd4$time,na.rm=T)

cd4 %>%
  select(treatmnt, baseage, CD4PCT, time, newpid, VISIT, CD4CNT) %>%
  group_by(newpid) %>%
  filter(n()>1) %>%
  mutate(time = time - min(time,na.rm=T)) %>% write_csv(file="./datasets/cd4/cd4_mod.csv")

cd4 <- read.csv("./datasets/cd4/cd4_mod.csv")
str(cd4)

cd4 %>%
  mutate(time = time/365) %>% 
  write_csv(file="./datasets/cd4/cd4_mod.csv")

cd4 <- read.csv("./datasets/cd4/cd4_mod.csv")
str(cd4)

with(cd4, plot(VISIT, time))

cd4 %>%
  ggplot(aes(x=VISIT, y=CD4PCT))+
  geom_jitter()+
  geom_smooth(method="lm") +
  facet_wrap(~treatmnt, ncol=3)

cd4$treatmnt <- factor(cd4$treatmnt)
contrasts(cd4$treatmnt)


model2 <- lmer(CD4PCT ~ baseage + time + time:treatmnt + (time |newpid), cd4)
summary(model2)


# pred_df <- expand.grid(time = range(cd4$time, na.rm=T),
#                        treatmnt = unique(cd4$treatmnt),
#                        baseage = mean(cd4$baseage, na.rm=T))
# pred_df$CD4PCT <- predict(model2, re.form=NA, newdata=pred_df)
# pred_df
# 
# complete_index <- complete.cases(cd4[,c("time","treatmnt","baseage","newpid","CD4PCT")])
# cd4$CD4PCT_pred[row_index] <- predict(model2)
# 
# ggplot(cd4, aes(x=time))+ 
#   facet_grid(.~treatmnt)+ 
#   geom_line(aes(y=CD4PCT_pred,color=factor(newpid), group=newpid))+
#   scale_color_discrete(guide=F)+
#   geom_line(data=pred_df,aes(y=CD4PCT),size=2)
#             
            

ggplot(cd4, aes(x=time, y=CD4PCT))+ 
  facet_grid(.~treatmnt)+ 
  geom_point()+ 
  geom_line(data=pred_df,size=2,color="blue")

pred_df <- expand.grid(time = seq(0, max(cd4$time,na.rm=T), length.out=100),
                       treatmnt = unique(cd4$treatmnt),
                       baseage = mean(cd4$baseage, na.rm=T))

bootFUN <- function(model){
   pred_cd4 <- predict(model, re.form=NA, newdata=pred_df)
   return(pred_cd4)
}

boot_res <- bootMer(model2, bootFUN, nsim=250)
dim(boot_res$t)

alpha <- 0.05 # significance level
pred_df$lower <- apply(boot_res$t, 2, function(x){quantile(x, probs = alpha/2)})
pred_df$upper <- apply(boot_res$t, 2, function(x){quantile(x, probs = 1- alpha/2)})
pred_df$CD4PCT <- predict(model2, re.form=NA, newdata=pred_df)

ggplot(cd4, aes(x=time, y=CD4PCT))+ 
  facet_grid(.~treatmnt)+ 
  geom_point()+ 
  geom_ribbon(data=pred_df,aes(ymin=lower, ymax=upper, fill=treatmnt), alpha=0.5)+
  geom_line(data=pred_df,size=2,aes(color=treatmnt))




## ----------------------------------------------------------
d <- read_delim("./datasets/BSA_abortion/SOCATT.DAT", col_names=F)
colnames(d) <- c("empty","district","id","year","resp","party","sasc","gender","age","religion")
d %>%
  select(-empty) %>%
  mutate(year = year+1982,
         party=case_when(
           party == 1 ~ "conservative",
           party == 2 ~ "labour",
           party == 3 ~ "Liberal/SDP/Alliance",
           party == 4 ~ "other",
           party == 5 ~ "none"
         ),
         sasc = case_when(
           sasc==1 ~ "middle",
           sasc==2 ~ "upper working",
           sasc==3 ~ "lower working"
         ),
         gender = ifelse(gender==1,"male","female"),
         religion = case_when(
           religion==1 ~ "catholic",
           religion==2 ~ "protestant",
           religion==3 ~ "other",
           religion==4 ~ "none"
         )) -> d
  
write.csv(d,file="./datasets/BSA_abortion/BSA_abortion.csv", row.names=F)
d <- read.csv("./datasets/BSA_abortion/BSA_abortion.csv")
str(d)

## ----------------------------------------------------------
data(cake)
str(cake)

fm1 <- lmer(angle ~ temp + (1|recipe), cake)
summary(fm1)

str(cake)
## 'temp' is continuous, 'temperature' an ordered factor with 6 levels
(fm1 <- lmer(angle ~ recipe * temperature + (1|recipe:replicate), cake, REML= FALSE))
(fm2 <- lmer(angle ~ recipe + temperature + (1|recipe:replicate), cake, REML= FALSE))
(fm3 <- lmer(angle ~ recipe + temp + (1|recipe:replicate), cake, REML= FALSE))
## and now "choose" :
anova(fm3, fm2, fm1)


## ----------------------------------------------------------
library(MEMSS)
data(MathAchieve)
data(MathAchSchool)
str(MathAchieve)

data("MathAchieve", package="MEMSS")

schooldat <- left_join(MathAchieve, MathAchSchool[,colnames(MathAchSchool)!="MEANSES"], by="School")
str(schooldat)
colnames(schooldat) <- tolower(colnames(schooldat))
schooldat <- schooldat[,1:8]

write_csv(schooldat, file="./datasets/mathachieve.csv")

mathachieve <- read.csv("./datasets/mathachieve.csv")
str(mathachieve)

mathachieve$size <- scale(mathachieve$size) # standardize the size variable

math.1 <- lmer(mathach ~ meanses + sector + ses + size+ (ses | school), data=mathachieve)
summary(math.1)


#  with two additional interaction terms meanses:cses and cses:sector as fixed effects


math.2 <- lmer(mathach ~ meanses + sector + ses + size+  meanses:ses +(ses | school), data=mathachieve)
summary(math.2)

math.3 <- lmer(mathach ~ meanses + sector + ses + size+ ses:sector +(ses | school), data=mathachieve)
summary(math.3)


math.3 <- lmer(mathach ~ meanses + sector + ses + size+  sex +(ses | school), data=mathachieve)
summary(math.3)




ggplot(mathachieve, aes(x=sector, y=ses))+
  geom_boxplot()

ggplot(mathachieve, aes(x=sector, y=meanses))+
  geom_boxplot()

mathachieve %>%
  group_by(school, sector) %>%
  summarise(ses = mean(ses)) %>%
  ggplot(aes(x=sector, y=ses))+
  geom_violin() 

