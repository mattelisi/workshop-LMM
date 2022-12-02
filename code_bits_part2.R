
rm(list=ls())
setwd("/mnt/D/Dropbox/sync/RHUL/stats_role/WORKSHOPS/workshop_LMM")

library(lme4)
data("sleepstudy")
m.ML <- lmer(Reaction ~ Days + (Days|Subject), data=sleepstudy, REML=F)
m.REML <- lmer(Reaction ~ Days + (Days|Subject), data=sleepstudy, REML=T)

VarCorr(m.ML)
VarCorr(m.REML)


### BOOTSTRAP

# simulate(m.REML, nsim=100, re.form = ~0)
# predict(m.REML, newdata=data.frame(Days=0:9), re.form=NA)

nd <- data.frame(Days=0:9, Subject=1)
simulate(m.REML, newdata=nd, nsim=100, re.form = ~0)


summary(m.REML)

# function that calculate how many days of sleep deprivation
# are needed to cause a 30% slowing in response time
my_function <- function(model){
  th <- (fixef(model)[1]*0.3)/fixef(model)[2]
  return(unname(th))
}
my_function(m.REML)

boot_res <- bootMer(m.REML, FUN=my_function, nsim=1000, re.form = NULL, type="parametric")
saveRDS(boot_res,"./misc/sleep_boot_res.RDS")
hist(boot_res$t, breaks=30)

str(boot_res[1:3])

library(boot)
boot.ci(boot_res)
boot.ci(boot_res, type=c("perc"))


## 

library(MEMSS)
library(lmerTest)
data(Machines)
str(Machines)

contrasts(Machines$Machine) <- contr.treatment(levels(Machines$Machine), base=2)
contrasts(Machines$Machine)

machine.mod <- lmer(score ~ Machine + (Machine|Worker), Machines)
summary(machine.mod)

N <- 10
N_rep <- 3


# getME(machine.mod, c("theta","fixef","sigma"))
fit_par <- getME(machine.mod, c("theta","sigma"))
fit_par$beta <- fixef(machine.mod)

# random-effects parameter estimates: these are parameterized as the relative Cholesky factors of each random effect term

Omega <- VarCorr(machine.mod)
cohen_d <- fixef(machine.mod)["MachineC"] / sqrt(fit_par$sigma^2/(N_rep-1)
                                           + Omega$Worker["(Intercept)","(Intercept)"] 
                                           + 0.5*Omega$Worker["MachineC","MachineC"]
                                           + Omega$Worker["(Intercept)","MachineC"])

attr(Omega$Worker, "stddev")
attr(Omega$Worker, "correlation")

# https://stats.stackexchange.com/questions/155424/how-does-lme4-optimize-the-variance-parameters-of-random-effects-theta-vector
W <- Omega$Worker
attr(W, "stddev") <- NULL
attr(W, "correlation") <- NULL

VarCov2Theta <- function(X){
  X <- t(chol(X)) 
  X_unpacked <- X[lower.tri(X, diag=T)]
  return(X_unpacked)
}



getME(machine.mod, "ST")$Worker %*%  t(getME(machine.mod, "ST")$Worker)

VarCov2Theta(W)
getME(machine.mod, "theta")

rr <- ranef(machine.mod,condVar=TRUE)
pv <- attr(rr[[1]],"postVar")
str(pv)

# ------------------------ 

sim_d <- expand.grid(Worker = factor(1:N),
                  Machine = unique(Machines$Machine),
                  rep = 1:3)

contrasts(sim_d$Machine) <- contr.treatment(levels(sim_d$Machine), base=2)

sim_d$score <- simulate(~ Machine + (Machine|Worker),
                        nsim=1,
                        family=gaussian,
                        newdata = sim_d,
                        newparams=fit_par,
                        use.u = FALSE)

simulate_data <- function(N, fit_par){
  
  sim_d <- expand.grid(Worker = factor(1:N),
                       Machine = unique(Machines$Machine),
                       rep = 1:3)
  
  contrasts(sim_d$Machine) <- contr.treatment(levels(sim_d$Machine), base=2)
  
  sim_d$score <- simulate(~ Machine + (Machine|Worker),
                          nsim=1,
                          family=gaussian,
                          newdata = sim_d,
                          newparams=fit_par,
                          use.u = FALSE)$sim_1
  return(sim_d)
}


test_significance <- function(sim_d){
  mod <- lmer(score ~ Machine + (Machine|Worker), sim_d)
  p_value <- summary(mod)$coefficients["MachineC","Pr(>|t|)"]
  return(p_value)
}


N_sim <- 200
N_workers <- c(6, 7, 8, 9, 10, 12, 15)
sim_res <- {}
for(w in N_workers){
  for(i in 1:N_sim){
    p_val <- test_significance(simulate_data(w, fit_par))
    sim_res <- rbind(sim_res, data.frame(N=w, p=p_val))
  }
}
saveRDS(sim_res, "./misc/machine_sim01.RDS")

str(sim_res)
hist(sim_res$p, breaks=100)
abline(v=0.05)

library(ggplot2)
library(tidyverse)

sim_res$significant <- ifelse(sim_res$p<0.05,1,0)

sim_res %>%
  group_by(N) %>%
  summarise(SE = sqrt((mean(significant) * (1 - mean(significant)))/length(significant)),
            significant = mean(significant)) %>%
  ggplot(aes(x=N, y=significant))+
  geom_line(color="blue")+
  geom_errorbar(aes(ymin=significant-SE, ymax=significant+SE),width=0,color="blue")+
  geom_point(color="blue",size=2)+
  geom_hline(yintercept = 0.8,lty=2)+
  labs(y="power")
  

sim_res %>%
  group_by(N) %>%
  summarise(SE = sqrt((mean(significant) * (1 - mean(significant)))/length(significant)),
            significant = mean(significant)) %>%
  ggplot(aes(x=N, y=significant))+
  geom_smooth(data=sim_res,method="glm", method.args=list(family=binomial(logit)))+
  geom_errorbar(aes(ymin=significant-SE, ymax=significant+SE),width=0, color="blue")+
  geom_point(color="blue")+
  geom_hline(yintercept = 0.8,lty=2)+
  labs(y="power")



#### CHANGE EFFECT SIZE
cohen_sigma <- sqrt(fit_par$sigma^2/(3-1)
                    + Omega$Worker["(Intercept)","(Intercept)"] 
                    + 0.5*Omega$Worker["MachineC","MachineC"]
                    + Omega$Worker["(Intercept)","MachineC"])

new_beta_MachineC <- 0.5 * cohen_sigma
fit_par$beta["MachineC"] <- new_beta_MachineC

sim_res$cohen_d <- cohen_d
saveRDS(sim_res, "./misc/machine_sim01_d.RDS")

N_sim <- 200
N_workers <- c(6, 7, 8, 9, 10, 12, 15, 20, 25, 30, 40, 50)
sim_res <- {}
for(w in N_workers){
  for(i in 1:N_sim){
    p_val <- test_significance(simulate_data(w, fit_par))
    sim_res <- rbind(sim_res, data.frame(N=w, p=p_val))
  }
}
sim_res$cohen_d <- 0.5
sim_res$significant <- ifelse(sim_res$p<0.05,1,0)

sim_res <- sim_res[sim_res$N<=50,]
saveRDS(sim_res, "./misc/machine_sim02_d.RDS")

sim_res_all <- rbind(readRDS("./misc/machine_sim01_d.RDS"), readRDS("./misc/machine_sim02_d.RDS"))

library(pwr)

eqttest <- data.frame(N=5:100)
eqttest$power <- pwr.t.test(n=eqttest$N, d=0.5, type="paired")$power

pwr.t.test(n=30, d=0.5, type="paired")

sim_res_all %>%
  group_by(N,cohen_d) %>%
  summarise(SE = sqrt((mean(significant) * (1 - mean(significant)))/length(significant)),
            significant = mean(significant)) %>%
  ggplot(aes(x=N, y=significant, group=cohen_d, color=cohen_d))+
  geom_line(data=eqttest,aes(x=N,y=power),color="dark grey", lty=1)+
  geom_line()+
  geom_errorbar(aes(ymin=significant-SE, ymax=significant+SE),width=0)+
  geom_point(size=2)+
  coord_cartesian(xlim=c(5,50))+
  geom_hline(yintercept = 0.8,lty=2)+
  labs(y="power")


############## sensitivity

VarCov2Theta <- function(X, sigma){
  X <- t(chol(X/sigma^2)) 
  X_unpacked <- X[lower.tri(X, diag=T)]
  return(X_unpacked)
}

#
fit_par <- getME(machine.mod, c("theta","sigma"))
fit_par$beta <- fixef(machine.mod)

# random-effects parameter estimates: these are parameterized as the relative Cholesky factors of each random effect term

Omega <- VarCorr(machine.mod)$Worker
theta_names <- names(fit_par$theta)

N_sim <- 200
N_workers <- c(6, 7, 8, 9, 10, 12, 15)
variance_scaling <- c(1.5, 2)
sim_res3 <- {}
for(vs in variance_scaling){
for(w in N_workers){
  for(i in 1:N_sim){
    Omega_i <- (vs*diag(attr(Omega, "stddev"))) %*% attr(Omega, "correlation") %*% (vs*diag(attr(Omega, "stddev")))
    fit_par$theta <- VarCov2Theta(Omega_i, fit_par$sigma)
    names(fit_par$theta) <- theta_names
    p_val <- test_significance(simulate_data(w, fit_par))
    sim_res3 <- rbind(sim_res3, data.frame(N=w, p=p_val, vs=vs))
  }
}
}

sim01 <- readRDS("./misc/machine_sim01_d.RDS")
sim01$vs <- 1
sim_res3$cohen_d <- NA
sim_res3$significant <- ifelse(sim_res3$p<0.05,1,0)
sim03 <- rbind(sim01, sim_res3)
str(sim03)
saveRDS(sim03, "./misc/machine_sim03.RDS")

sim03 <- readRDS("./misc/machine_sim03.RDS")
sim03 %>%
  group_by(N,vs) %>%
  summarise(SE = sqrt((mean(significant) * (1 - mean(significant)))/length(significant)),
            significant = mean(significant)) %>%
  ggplot(aes(x=N, y=significant, group=vs, color=vs))+
  geom_line()+
  scale_color_continuous(name="variance\nscaling")+
  geom_errorbar(aes(ymin=significant-SE, ymax=significant+SE),width=0)+
  geom_point(size=2)+
  geom_hline(yintercept = 0.8,lty=2)+
  labs(y="power")

