library(tidyverse)
library(readr)
library(stringr)

partial_data <- read.csv("partial_data.csv",
                         row.names = 1)

this = str_split(partial_data$Sample, 
                 "_", 
                 simplify = T)
Treatment = substr(this[,2], 1, 1)
 Mesocosm = substr(this[,2], 1, 2)
Replicate = substr(this[,2], 2, nchar(this[,2])-1)

x = (.1*.1*10) / (pi*(25/2)^2)

d = partial_data %>% 
  add_column(.after = "Sample",
             Time = this[,1],
             Treatment = Treatment,
              Mesocosm = Mesocosm,
             Replicate = Replicate,
             volume = .$Fields_counted*x*1e3) %>% 
  mutate(Treatment = fct_relevel(Treatment, 
                                 c("C","D","I","E")),
         pHF = HF_total/(HF_total+MF+PF),
         pMF = MF/(HF_total+MF+PF),
         pPF = PF/(HF_total+MF+PF)) %>% 
  filter(Time == "Ex2")

d$Y <- cbind( 
             pMF = d$pMF, pHF = d$pHF,
             pPF = d$pPF)

library(ggeffects)
m = glmmTMB(HF_total ~ 1
            + Treatment 
            + (1|Mesocosm),
            offset = log(Fields_counted),
            family = nbinom2(),
            data = d);summary(m)
m = glmmTMB(HF_total/volume ~ 0
            + Treatment 
            + (1|Mesocosm),
            #family = nbinom2(),
            data = d);summary(m)
check_model(m)
plot(ggpredict(m, terms = c("Treatment")))
library(DHARMa)
simulateResiduals(fittedModel = m, plot = T)

library(brms)


mb = brm(bf(HF_total ~ Treatment 
            + (1|Mesocosm)
            + offset(log(Fields_counted))),
         family = negbinomial(link = "log"),
         chains = 4,
         iter = 2000,
         cores = 4,
         backend = "cmdstanr", 
         data = d)

mb = brm(bf(HF_total/volume ~ 0
            + Treatment 
            + (1|Mesocosm)),
         chains = 4,
         iter = 20000,
         cores = 4,
         adapt_delta = 0.99,
         backend = "cmdstanr", 
         data = d)

mb = brm(bf(PF/volume ~ 1
            + Treatment 
            + (1|Mesocosm)),
         chains = 4,
         iter = 20000,
         cores = 4,
         adapt_delta = 0.99,
         backend = "cmdstanr", 
         data = d)

mb = brm(bf(HF_total|trials(MF+PF+HF_total) ~ 1
            + Treatment 
            + (1|Mesocosm)),
         family = binomial(),
         chains = 4,
         iter = 20000,
         cores = 4,
         adapt_delta = 0.95,
         backend = "cmdstanr", 
         data = d)

pp_check(mb, ndraws = 1e2)
plot(mb)
summary(mb)
plot(conditional_effects(mb), ask = FALSE)

mb = brm(bf(Y ~ 1
            + Treatment 
            + (1|Mesocosm)),
         family = dirichlet(),
         chains = 4,
         iter = 2000,
         cores = 4,
         adapt_delta = 0.99,
         backend = "cmdstanr", 
         data = d)

pp_check(mb, ndraws = 1e2)
plot(mb)
summary(mb)
plot(conditional_effects(mb, 
                         categorical = T), 
     ask = FALSE)


## predicted responses
pp <- predict(mb)
head(pp)

# DOES NOT WORK!!!
library(marginaleffects)
library(ggdist)
mfx = marginaleffects(mb,
                      re_formula = NA,
                      variables = "Treatment") %>% 
  posteriordraws()

ggplot(mfx, aes(x = draw,
                fill = contrast)) +
  stat_halfeye(slab_alpha = .5)  +
  labs(x = "difference from Control",
       y = "")
