library(brms)

# abundance

pm = brm(bf(aPF ~ Treatment 
            + Incubation
            + (1|Mesocosm)),
         family = lognormal(link = "identity", link_sigma = "log"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta=0.95),
         seed=543,
         backend = "cmdstanr", 
         data = dat)
pp_check(pm, ndraws = 100)
summary(pm, prob = .9)
plot(conditional_effects(pm, categorical = F, prob = .9), ask = FALSE)
plot(pm)


hm = brm(bf(aHF ~ Treatment 
            + Incubation
            + (1|Mesocosm)),
         family = lognormal(link = "identity", link_sigma = "log"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta=0.95),
         seed=543,
         backend = "cmdstanr", 
         data = dat)
pp_check(hm, ndraws = 100)
summary(hm, prob = .9)
plot(conditional_effects(hm, categorical = F, prob = .9), ask = FALSE)
plot(hm)


mm = brm(bf(aMFc ~ Treatment 
            + Incubation
            + (1|Mesocosm)),
         family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta=0.95),
         backend = "cmdstanr", 
         data = dat)
pp_check(mm, ndraws = 100)
summary(mm, prob = .9)
plot(conditional_effects(mm, categorical = F, prob = .9), ask = FALSE)
plot(mm)

# biomass
bh = brm(bf(HF ~ Treatment 
            + Incubation
            + (1|Mes_ID)),
         family = lognormal(link = "identity", link_sigma = "log"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta=0.95),
         backend = "cmdstanr", 
         data = bdat)
pp_check(bh, ndraws = 100)
summary(bh, prob = .9)
plot(conditional_effects(bh, categorical = F, prob = .9), ask = FALSE)
plot(bh)


bp = brm(bf(PF ~ Treatment 
            + Incubation
            + (1|Mes_ID)),
         family = lognormal(link = "identity", link_sigma = "log"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta=0.95),
         backend = "cmdstanr", 
         data = bdat)
pp_check(bp, ndraws = 100)
summary(bp, prob = .9)
plot(conditional_effects(bp, categorical = F, prob = .9), ask = FALSE)
plot(bp)

bm = brm(bf(MF ~ Treatment 
            + Incubation
            + (1|Mes_ID)),
         family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta=0.95),
         backend = "cmdstanr", 
         data = bdat)
pp_check(bm, ndraws = 100)
summary(bm, prob = .9)
plot(conditional_effects(bm, categorical = F, prob = .9), ask = FALSE)
plot(bm)
         
# ingestion rates 

mir = brm(bf(MFir ~ Treatment 
             + Incubation
             + (1|Mesocosm),
             hu ~ Treatment 
             + Incubation
             + (1|Mesocosm)),
          family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
          chains = 4,
          iter = 2000,
          cores = 4,
          control = list(adapt_delta=0.95),
          seed = 543,
          backend = "cmdstanr", 
          data = IR1)

pp_check(mir, ndraws = 100)
summary(mir, prob = .9)
## we want Rhat to be 1, Bulk_ESS & Tail_ESS to make sence withthe total post-warmup draws
## something is significant when CI have the same sign (+-)
plot(conditional_effects(mir, categorical = F, prob = .9), ask = FALSE)
plot(mir)

hir = brm(bf(HFir ~ Treatment 
             + Incubation
             + (1|Mesocosm),
             hu ~ Treatment 
             + Incubation
             + (1|Mesocosm)),
          family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
          chains = 4,
          iter = 2000,
          cores = 4,
          control = list(adapt_delta=0.95),
          backend = "cmdstanr", 
          data = IR1)

pp_check(hir, ndraws = 100)
summary(hir, prob = .9)
plot(conditional_effects(hir, categorical = F, prob = .9), ask = FALSE)
plot(hir)
