rm(list=ls())
library(ggplot2)
library(tidyverse)

# problems with complete pooling 

set.seed(1)
n <- 12
n_obs <- 30
Z <- MASS::mvrnorm(n,mu=c(0,0), Sigma=matrix(c(1,-0.75,-0.75,1),ncol=2))
beta_1 <- rnorm(n, mean=0.7, sd=0.2)
sigma_y <- 0.25
sigma_x <- 0.5
simd <- expand.grid(obs=1:n_obs, id=1:n)
simd$x <- NA
simd$y <- NA
for(i in unique(simd$id)){
  simd$x[simd$id==i] <- Z[i,1] + rnorm(n_obs,0,sigma_x)
  simd$y[simd$id==i] <- Z[i,2] + simd$x[simd$id==i]*beta_1[i] - Z[i,1] + rnorm(n_obs,0,sigma_y)
}
simd$id <- paste("sj",simd$id,sep="")

ggplot(simd, aes(x=x, y=y, color=id))+
  geom_point()+
  scale_color_manual(values=rep("black",n))+
  geom_smooth(method="lm",se=F,aes(group=1))


ggplot(simd, aes(x=x, y=y, color=id))+
  geom_point()+
  scale_color_viridis_d()+
  geom_smooth(method="lm",se=F) + 
  geom_smooth(method="lm",se=F,aes(group=1))



# pooling in sleepstudy
library(lme4)

m.0 <- lm(Reaction ~ Days, data=sleepstudy)

d0 <- sleepstudy
d0$prediction <- predict(m.0)
d0$pooling <- "complete"


# estimate individual regressions
m.list <- lmList(Reaction ~ Days | Subject, data=sleepstudy)
d1 <- sleepstudy
d1$prediction <- predict(m.list)
d1$pooling <- "no"

d <- rbind(d0,d1)


ggplot(d, aes(x=Days,y=Reaction))+
  facet_wrap(~Subject,ncol=9)+
  geom_line(aes(y=prediction, color=pooling),size=0.8)+
  geom_point()+
  scale_x_continuous(breaks=seq(0,10,2))+
  labs(x="Days of sleep deprivation",y="Mean reaction times [ms]")


m.lmm <- lmer(Reaction ~ Days + (Days|Subject), data=sleepstudy)
d2 <- sleepstudy
d2$prediction <- predict(m.lmm)
d2$pooling <- "partial (LMM)"

d <- rbind(d0,d1, d2)

ggplot(d, aes(x=Days,y=Reaction))+
  facet_wrap(~Subject,ncol=9)+
  geom_line(aes(y=prediction, color=pooling),size=0.8)+
  geom_point()+
  scale_x_continuous(breaks=seq(0,10,2))+
  labs(x="Days of sleep deprivation",y="Mean reaction times [ms]")


## ----------------------------------------------------------
# extreme data loss

set.seed(1124)
  
  n_id <- length(unique(sleepstudy$Subject))
  n_obs <- rpois(n_id, lambda=0.9) + 2

  
  slice_consecutive <- function(x, n) {
    i_start <- sample(1:(nrow(x)-n), 1)
    i_end <- (i_start+n-1)
    return(x[i_start:i_end,])
  }
  
# sleepstudy %>%
#   group_split(Subject) %>%
#   map2_dfr(n_obs, ~ slice_sample(.x, n = .y)) -> d_miss

sleepstudy %>%
  group_split(Subject) %>%
  map2_dfr(n_obs, ~ slice_consecutive(.x, n = .y)) -> d_miss


d %>%
  filter(pooling=="no") %>%
ggplot(aes(x=Days,y=Reaction))+
  facet_wrap(~Subject,ncol=9)+
  geom_line(aes(y=prediction, color=pooling),size=0.4, lty=2)+
  geom_point(pch=21, color="dark grey")+
  geom_point(data=d_miss) +
  geom_smooth(data=d_miss,method="lm",fullrange=T,se=F) +
  scale_x_continuous(breaks=seq(0,10,2))+
  labs(x="Days of sleep deprivation",y="Mean reaction times [ms]")


m.lmm.ms <- lmer(Reaction ~ Days + (Days|Subject), data=d_miss)
summary(m.lmm.ms)

m.list.ms <- lmList(Reaction ~ Days | Subject, data=d_miss)
mean(coefficients(m.list.ms)$Days, na.rm=T)
mlisi::se(coefficients(m.list.ms)$Days)

d$prediction_miss <- NA
d$prediction_miss[d$pooling=="no"] <- predict(m.list.ms, newdata = d[d$pooling=="no",])
d$prediction_miss[d$pooling=="partial (LMM)"] <- predict(m.lmm.ms, newdata = d[d$pooling=="partial (LMM)",])



d %>% 
  filter(pooling!="complete") %>%
  ggplot(aes(x=Days,y=Reaction))+
  facet_wrap(~Subject,ncol=9)+
  geom_point(pch=21, color="dark grey")+
  geom_point(data=d_miss) +
  geom_smooth(method="lm",se=F,size=0.4, lty=2,color="dark grey") +
  geom_line(aes(y=prediction_miss, color=pooling),size=0.8)+
  scale_x_continuous(breaks=seq(0,10,2))+
  labs(x="Days of sleep deprivation",y="Mean reaction times [ms]")

# compare predictive r-squared
d %>% 
  filter(pooling!="complete") %>%
  anti_join(d_miss, by=c("Reaction", "Subject","Days")) %>% 
  mutate(sq_err = (prediction_miss - Reaction)^2) %>%
  group_by(pooling) %>%
  summarise(RMSE = sqrt(mean(sq_err)))


## VISUALIZE SHRINKGAGE

# analyze shrinking
m.single <- coef(m.list.ms)
par.mixed <- as.matrix(ranef(m.lmm.ms)$Subject) + mlisi::repmat(t(as.matrix(fixef(m.lmm.ms))),18,1)

# plot parameters from mixed model and individual fit
plot(m.single[,1], m.single[,2], xlab="Intercept",ylab="Slope",pch=19,cex.lab=1.2,col="dark grey",
     xlim=c(100,450),ylim=c(-60,60))

# draw ellipse illustrating covariance of random effects
vcov_m.1 <- matrix(as.vector(VarCorr(m.lmm.ms)$Subject),ncol=2)
mixtools::ellipse(c(mean(par.mixed[,1]), mean(par.mixed[,2])), sigma=vcov_m.1,alpha=0.05, col="grey", lty=2)
mixtools::ellipse(c(mean(par.mixed[,1]), mean(par.mixed[,2])), sigma=vcov_m.1,alpha=0.001, col="grey", lty=2)
mixtools::ellipse(c(mean(par.mixed[,1]), mean(par.mixed[,2])), sigma=vcov_m.1,alpha=0.5, col="grey", lty=2)

points(mean(m.single[,1]), mean(m.single[,2]),pch=19,col="dark grey",cex=2)
points(mean(par.mixed[,1]), mean(par.mixed[,2]),pch=21,col="black",cex=2,lwd=2)
text(m.single[,1], m.single[,2],labels=rownames(m.single),pos=1,cex=0.6)
points(par.mixed[,1], par.mixed[,2])
arrows(m.single[,1], m.single[,2],par.mixed[,1], par.mixed[,2], length=0.1)
legend("bottomright",c("single-subject fits","multilevel model"),pch=c(19,21),col=c("dark grey", "black"),bty="n",cex=0.6)



### multivariate normals

n_obs <- 250
r00 <- MASS::mvrnorm(n_obs,c(0,0), Sigma=matrix(data=c(1,0,0,1),nrow=2))
r03 <- MASS::mvrnorm(n_obs,c(0,0), Sigma=matrix(data=c(1,0.5,0.5,1),nrow=2))
r07 <- MASS::mvrnorm(n_obs,c(0,0), Sigma=matrix(data=c(1,0.85,0.85,1),nrow=2))
dat <- data.frame(rbind(r00, r03, r07))
dat$correlation <- c(rep("rho==0",n_obs),
                     rep("rho==0.5",n_obs),
                     rep("rho==0.85",n_obs))

ggplot(dat,aes(x=X1,y=X2))+
  geom_point()+
  facet_grid(.~correlation,labeller = label_parsed)


