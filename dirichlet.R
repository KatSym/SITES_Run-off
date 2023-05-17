data <- dat %>% 
  select(Incubation,Mes_ID,Treatment, aHF, aPF, aMFc) %>% 
  mutate(tot = aHF+aPF+aMFc,
         # calculate proportions
         h = aHF/tot,
         p = aPF/tot,
         m = aMFc/tot) %>% 
  as.data.frame()

library(DirichletReg)

d <- data
d$Y <- DR_data(data[,8:10])
plot(d$Y)

dir1 <- DirichReg(Y ~ Treatment * Incubation, d)
summary(dir1)
coef(dir1)



dbio <- bdat %>% 
  select(-Replicate) %>% 
  mutate(tot = HF+PF+MF,
         # calculate proportions
         h = HF/tot,
         p = PF/tot,
         m = MF/tot,
         
         # first log transform then proportion, I tried directly transforming 
         # the proportion but the Y variable must sum to 1
         lpf = log1p(PF),
         lhf = log1p(HF),
         lmf = log1p(MF),
         ltot = lpf+ lhf+ lmf,
         lh = lhf/ltot,
         lp = lpf/ltot,
         lm = lmf/ltot) 

dbio$Y <- cbind(p = dbio$p, h = dbio$h, m = dbio$m)
dbio$lY <- cbind(lp = dbio$lp, lh = dbio$lh, lm = dbio$lm)





# BAYESIAN


# https://numerilab.io/en/analysis_projects/DirichletElections 


bd2 <- brm(bf(Y ~ Treatment 
              * Incubation
              + (1|Mes_ID)),
           family = dirichlet(),
           chains = 4,
           iter = 2000,
           cores = 4,
           backend = "cmdstanr", 
           data = dbio)

# pp_check(bd2, ndraws = 100)
summary(bd2, prob = .9)
plot(conditional_effects(bd2, categorical = T, prob = .9), ask = FALSE)
plot(bd2, ask = F)

log.mdir <- brm(bf(lY ~ Treatment 
              * Incubation
              + (1|Mes_ID)),
           family = dirichlet(),
           chains = 4,
           iter = 2000,
           cores = 4,
           backend = "cmdstanr", 
           data = dbio)
summary(log.mdir, prob = .9)
plot(conditional_effects(log.mdir, categorical = T, prob = .9), ask = FALSE)
plot(log.mdir, ask = F)
