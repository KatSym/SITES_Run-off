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




# BAYESIAN


# https://numerilab.io/en/analysis_projects/DirichletElections 


library(tidyverse)
library(rgdal)
library(knitr)
library(kableExtra)

data$Y <- cbind(h = data$h, m = data$m, 
               p = data$p)

library(brms)
bd1 <- bf(Y ~ Treatment 
          + Incubation
          + (1|Mes_ID),
           family = dirichlet)

get_prior(bd1, data = data) %>% 
  kable() %>% 
  kable_styling() %>%
  scroll_box(width = "100%", height = "300px")

curve(dnorm(x, 0, 1), from = -10, to = 10)
curve(dnorm(x, 0, 2), from = -10, to = 10, add = T, col = "blue")
curve(dnorm(x, 0, 10), from = -10, to = 10, add = T, col = "red")

## Prior chosen by brms
curve(dgamma(x, 0.01, rate = 0.01), to = 30, ylim = c(0, 0.040))
## Prior we will define
curve(dstudent_t(x, 3, 0, 10), to = 30, add = T, col = "red")

# I HAVE 0s SO i CANNOT RUN IT LIKE THAT

## Set prior with sd = 1
priors_1 <- c(prior(normal(0,1), class = b),
              prior(normal(0,1), class = Intercept),
              prior(student_t(3, 0, 10), class = phi))
## sample from the priors
ppp_1 <- brm(formula = bd1, prior = priors_1, 
             data = data, sample_prior = "only",
             chains = 1, cores = 1)
## Set prior with sd = 2
priors_2 <- c(prior(normal(0,2), class = b),
              prior(normal(0,1), class = Intercept),
              prior(student_t(3, 0, 10), class = phi))
## sample from the priors
ppp_2 <- brm(formula = bd1, prior = priors_2, 
             data = data, sample_prior = "only",
             chains = 1, cores = 1)
## Set prior with sd = 2
priors_10 <- c(prior(normal(0,10), class = b),
               prior(normal(0,1), class = Intercept),
               prior(student_t(3, 0, 10), class = phi))
## sample from the priors
ppp_10 <- brm(formula = bd1, prior = priors_10, 
              data = data, sample_prior = "only",
              chains = 1, cores = 1)







bd2 <- brm(bf(HF_dataal ~ Treatment 
              + Incubation
              + (1|Mesocosm),
           family = dirichlet(),
           chains = 4,
           iter = 2000,
           cores = 4,
           backend = "cmdstanr", 
           data = data)