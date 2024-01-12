
library(tidyverse)
library(brms)
library(tidybayes)
library(modelr)
library(emmeans)

library(ggpubr)
library(grid)


# all models in Bayes

# data -----
load("all_data.RData")

mytr <- function(x){
  log10(x+1)
}

# all days 
# no transformation
env <- envir %>% 
  select(-c(TOC, DN, DOsat, DOconc)) %>% 
  relocate(c(ExpDay, Treatment), .before = Mes_ID) %>% 
  mutate(Nutr = TN*1000 + TP,
         # across(DOC:Nutr, scale),
         # scale all predictors
         # across(DOC:Nutr, mysc),
                # ~(.-mean(.)) / sd(.)
                ExpDay = factor(ExpDay))
  
dat <- data %>% 
  relocate(c(MF_Ir, HF_Ir, MF_Gr, HF_Gr), .after = biovol_MF) %>% 
  # log(x + 1) transormed
    mutate(across(4:12, mytr),
           ExpDay = factor(ExpDay)) 

# take out day 5
dat1 <- dat %>% filter(ExpDay !="5") %>% droplevels() %>% 
  mutate(Treatment = as.factor(Treatment))

zoop <-   rot %>% 
  select(-rotif.sp) %>% 
  filter(Mes_ID %in% c(1, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16)) %>% 
  group_by(ExpDay,Treatment, Mes_ID) %>% 
  summarise_all(mean) %>% 
  ungroup() 


# models -----

## background data -----
abs.m = brm(bf(a420 ~ 
                 + Treatment*ExpDay
               + (1|Mes_ID)),
            family = gaussian(link = "identity"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta = 0.99),
            seed = 543,
            backend = "cmdstanr", 
            data = env %>% filter(!is.na(a420)),
            file = "models/231029_abs",
            file_refit = "on_change"
)

pp_check(abs.m)
summary(abs.m)
conditional_effects(abs.m, effects = "Treatment:ExpDay")
conditional_effects(abs.m, effects = "ExpDay:Treatment")

emTrDay = emmeans(abs.m, ~ Treatment|ExpDay)
summary(emTrDay, point = "mean")
pTrDay = pairs(emTrDay)
summary(pTrDay, point = "mean")

emDayTr = emmeans(abs.m, ~ ExpDay|Treatment)
summary(emDayTr, point = "mean")
pDayTr = pairs(emDayTr)
summary(pDayTr, point = "mean")




par.m = brm(bf(PAR ~ 
                + Treatment*ExpDay
              + (1|Mes_ID)),
           family = gaussian(link = "identity"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta = 0.99),
           seed = 543,
           backend = "cmdstanr", 
           data = env %>% filter(!is.na(PAR)),
           file = "models/231025_par",
           file_refit = "on_change"
)

pp_check(par.m)
summary(par.m)
conditional_effects(par.m, effects = "Treatment:ExpDay")
conditional_effects(par.m, effects = "ExpDay:Treatment")

emTrDay = emmeans(par.m, ~ Treatment|ExpDay)
summary(emTrDay, point = "mean")
pTrDay = pairs(emTrDay)
summary(pTrDay, point = "mean")

emDayTr = emmeans(par.m, ~ ExpDay|Treatment)
summary(emDayTr, point = "mean")
pDayTr = pairs(emDayTr)
summary(pDayTr, point = "mean")

# 
# doc.m = brm(bf(DOC ~ 
#                 + Treatment*ExpDay
#               + (1|Mes_ID)),
#            family = gaussian(link = "identity"),
#            # family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
#            chains = 4,
#            iter = 2000,
#            # prior = dpr,
#            cores = 4,
#            control = list(adapt_delta = 0.99, 
#                           max_treedepth = 15),
#            seed = 543,
#            backend = "cmdstanr", 
#            data = env%>% filter(!is.na(DOC)),
#            file = "models/231026_doc",
#            file_refit = "on_change"
# )
# pp_check(doc.m)
# summary(doc.m)
# conditional_effects(doc.m, effects = "Treatment:ExpDay")
# conditional_effects(doc.m, effects = "ExpDay:Treatment")
# 
# emTrDay = emmeans(doc.m, ~ Treatment|ExpDay)
# summary(emTrDay, point = "mean")
# pTrDay = pairs(emTrDay)
# summary(pTrDay, point = "mean")
# 
# emDayTr = emmeans(doc.m, ~ ExpDay|Treatment)
# summary(emDayTr, point = "mean")
# pDayTr = pairs(emDayTr)
# summary(pDayTr, point = "mean")


chla.m <- brm(bf(Chla ~ 
                   + Treatment*ExpDay
                 + (1|Mes_ID)),
              family = gaussian(link = "identity"),
              chains = 4,
              iter = 2000,
              cores = 4,
              control = list(adapt_delta = 0.99,
                             max_treedepth = 12,
                             step_size = 0.2),
              seed = 543,
              backend = "cmdstanr", 
              data = env %>% filter(!is.na(Chla)),
              file = "models/231204_chla",
              # file = "models/231212_chla",
              file_refit = "on_change")

pp_check(chla.m)
summary(chla.m)
conditional_effects(chla.m, effects = "Treatment:ExpDay")
conditional_effects(chla.m, effects = "ExpDay:Treatment")

emTrDay = emmeans(chla.m, ~ Treatment|ExpDay)
summary(emTrDay, point = "mean")
pTrDay = pairs(emTrDay)

emDayTr = emmeans(chla.m, ~ ExpDay|Treatment)
summary(emDayTr, point = "mean")
pDayTr = pairs(emDayTr)

summary(pTrDay, point = "mean")
summary(pDayTr, point = "mean")


nut.m = brm(bf(Nutr ~ 
                 + Treatment*ExpDay
               + (1|Mes_ID)),
            family = gaussian(link = "identity"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 12,
                           step_size = 0.2),
            seed = 543,
            backend = "cmdstanr", 
            data = env%>% filter(!is.na(Nutr)),
            file = "models/231026_nutr",
            file_refit = "on_change")

pp_check(nut.m)
summary(nut.m)
conditional_effects(nut.m, effects = "Treatment:ExpDay")
conditional_effects(nut.m, effects = "ExpDay:Treatment")

emTrDay = emmeans(nut.m, ~ Treatment|ExpDay)
summary(emTrDay, point = "mean")
pTrDay = pairs(emTrDay)

emDayTr = emmeans(nut.m, ~ ExpDay|Treatment)
summary(emDayTr, point = "mean")
pDayTr = pairs(emDayTr)

summary(pTrDay, point = "mean")
summary(pDayTr, point = "mean")


tn.m = brm(bf(TN ~ 
                 + Treatment*ExpDay
               + (1|Mes_ID)),
            family = gaussian(link = "identity"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 12,
                           step_size = 0.2),
            seed = 543,
            backend = "cmdstanr", 
            data = env %>% filter(!is.na(TN)),
            file = "models/231120_TN",
            file_refit = "on_change"
)
pp_check(tn.m)
summary(tn.m)
conditional_effects(tn.m, effects = "Treatment:ExpDay")
conditional_effects(tn.m, effects = "ExpDay:Treatment")

emTrDay = emmeans(tn.m, ~ Treatment|ExpDay)
summary(emTrDay, point = "mean")
pTrDay = pairs(emTrDay)

emDayTr = emmeans(tn.m, ~ ExpDay|Treatment)
summary(emDayTr, point = "mean")
pDayTr = pairs(emDayTr)

summary(pTrDay, point = "mean")
summary(pDayTr, point = "mean")



tp.m = brm(bf(TP ~ 
                + Treatment*ExpDay
              + (1|Mes_ID)),
           family = gaussian(link = "identity"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta = 0.99,
                          max_treedepth = 12,
                          step_size = 0.2),
           seed = 543,
           backend = "cmdstanr", 
           data = env%>% filter(!is.na(TP)),
           file = "models/231120_TP",
           file_refit = "on_change"
)
pp_check(tp.m)
summary(tp.m)
conditional_effects(tp.m, effects = "Treatment:ExpDay")
conditional_effects(tp.m, effects = "ExpDay:Treatment")

emTrDay = emmeans(tp.m, ~ Treatment|ExpDay)
summary(emTrDay, point = "mean")
pTrDay = pairs(emTrDay)

emDayTr = emmeans(tp.m, ~ ExpDay|Treatment)
summary(emDayTr, point = "mean")
pDayTr = pairs(emDayTr)

summary(pTrDay, point = "mean")
summary(pDayTr, point = "mean")


## rotifers ----
rot.m = brm(bf(log1p(Rotif) ~ 
                 + Treatment*ExpDay
               + (1|Mes_ID)),
            family = gaussian(link = "identity"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 12,
                           step_size = 0.2),
            seed = 543,
            backend = "cmdstanr", 
            data = zoop %>% filter(!is.na(Rotif)),
            file = "models/231026_rotif",
            file_refit = "on_change"
)
pp_check(rot.m, ndraws = 100)
summary(rot.m)
conditional_effects(rot.m, effects = "Treatment:ExpDay")
conditional_effects(rot.m, effects = "ExpDay:Treatment")

emTrDay = emmeans(rot.m, ~ Treatment|ExpDay)
summary(emTrDay, point = "mean")
pTrDay = pairs(emTrDay)

emDayTr = emmeans(rot.m, ~ ExpDay|Treatment)
summary(emDayTr, point = "mean")
pDayTr = pairs(emDayTr)

summary(pTrDay, point = "mean")
summary(pDayTr, point = "mean")

## abundances ---- 

# heter. bact
bact.m = brm(bf(HB_abund ~ 
               + Treatment*ExpDay
             + (ExpDay|Mes_ID)),
          family = gaussian(link = "identity"),
          chains = 4,
          iter = 2000,
          cores = 4,
          control = list(adapt_delta = 0.99),
          seed = 543,
          backend = "cmdstanr", 
          data = dat,
          file = "models/231210_bacy.m"
)
pp_check(bact.m)
summary(bact.m)
conditional_effects(bact.m, effects = "Treatment:ExpDay")
conditional_effects(bact.m, effects = "ExpDay:Treatment")

dddd = emmeans(bact.m, ~ Treatment|ExpDay)
summary(dddd, point = "mean")
pppp = pairs(dddd)
summary(pppp, point = "mean")

dddd1 = emmeans(bact.m, ~ ExpDay|Treatment)
summary(dddd1, point = "mean")
pppp1 = pairs(dddd1)
summary(pppp1, point = "mean")

# cyanobacteria
cy.m = brm(bf(CY_abund ~ 
               + Treatment*ExpDay
             + (ExpDay|Mes_ID)),
          family = gaussian(link = "identity"),
          chains = 4,
          iter = 2000,
          cores = 4,
          control = list(adapt_delta = 0.99),
          seed = 543,
          backend = "cmdstanr", 
          data = dat,
          file = "models/231210_cy.m"
)
pp_check(cy.m)
summary(cy.m)
conditional_effects(cy.m, effects = "Treatment:ExpDay")
conditional_effects(cy.m, effects = "ExpDay:Treatment")

dddd = emmeans(cy.m, ~ Treatment|ExpDay)
summary(dddd, point = "mean")
pppp = pairs(dddd)
summary(pppp, point = "mean")

dddd1 = emmeans(cy.m, ~ ExpDay|Treatment)
summary(dddd1, point = "mean")
pppp1 = pairs(dddd1)
summary(pppp1, point = "mean")

# phototrophs
p.m = brm(bf(PF_abund ~ 
              + Treatment*ExpDay
              + (ExpDay|Mes_ID)),
           family = gaussian(link = "identity"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta = 0.99),
           seed = 543,
           backend = "cmdstanr", 
           data = dat,
           file = "models/231210_phot.m"
)
pp_check(p.m)
summary(p.m)
conditional_effects(p.m, effects = "Treatment:ExpDay")
conditional_effects(p.m, effects = "ExpDay:Treatment")

dddd = emmeans(p.m, ~ Treatment|ExpDay)
summary(dddd, point = "mean")
pppp = pairs(dddd)
summary(pppp, point = "mean")

dddd1 = emmeans(p.m, ~ ExpDay|Treatment)
summary(dddd, point = "mean")
pppp1 = pairs(dddd)
summary(pppp, point = "mean")

# heterotrophs
h.m = brm(bf(HF_abund ~ 
              + Treatment*ExpDay
              + (ExpDay|Mes_ID)),
           family = gaussian(link = "identity"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta = 0.99),
           seed = 543,
           backend = "cmdstanr", 
           data = dat,
           file = "models/231210_het.m"
)
pp_check(h.m)
summary(h.m)
conditional_effects(h.m, effects = "Treatment:ExpDay")
conditional_effects(h.m, effects = "ExpDay:Treatment")

dddd = emmeans(h.m, ~ Treatment|ExpDay)
summary(dddd, point = "mean")
pppp = pairs(dddd)
summary(pppp, point = "mean")

dddd1 = emmeans(h.m, ~ ExpDay|Treatment)
summary(dddd1, point = "mean")
pppp1 = pairs(dddd1)
summary(pppp1, point = "mean")

# mixotrophs
m.m1 = brm(bf(MF_abund ~ 
               + Treatment*ExpDay
             + (ExpDay|Mes_ID)),
          # family = gaussian(link = "identity"),
          family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
          chains = 4,
          iter = 2000,
          cores = 4,
          control = list(adapt_delta = 0.99),
          seed = 543,
          backend = "cmdstanr", 
          data = dat1,
          file = "models/231210_mix1.m"
)


pp_check(m.m1)

summary(m.m1)
conditional_effects(m.m1, effects = "Treatment:ExpDay")
conditional_effects(m.m1, effects = "ExpDay:Treatment")

dddd = emmeans(m.m1, ~ Treatment|ExpDay)
summary(dddd, point = "mean")
pppp = pairs(dddd)
summary(pppp, point = "mean")

dddd1 = emmeans(m.m1, ~ ExpDay|Treatment)
summary(dddd1, point = "mean")
pppp1 = pairs(dddd1)
summary(pppp1, point = "mean")


## Rates ----
#without transformation

# mix ingestion
mir.m1 = brm(bf(MF_Ir ~ 
               + Treatment*ExpDay
             + (ExpDay|Mes_ID)),
          # family = gaussian(link = "identity"),
          family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
          chains = 4,
          iter = 2000,
          cores = 4,
          control = list(adapt_delta = 0.99),
          seed = 543,
          backend = "cmdstanr", 
          data = dat1,
          file = "models/231025_mIR.m1",
          file_refit = "on_change"
)

mix.m <- brm(file = "models/231025_mIR.m")
pp_check(mir.m)
summary(mir.m1)

loo(mir.m, mir.m1, compare = T)

conditional_effects(mir.m1, effects = "Treatment:ExpDay")
conditional_effects(mir.m1, effects = "ExpDay:Treatment")

dddd = emmeans(mir.m1, ~ Treatment|ExpDay)
summary(dddd, point = "mean")
pppp = pairs(dddd)
summary(pppp, point = "mean")

dddd1 = emmeans(mir.m1, ~ ExpDay|Treatment)
summary(dddd1, point = "mean")
pppp1 = pairs(dddd1)
summary(pppp1, point = "mean")

# het ingestion
hir.m1 = brm(bf(HF_Ir ~ 
                 + Treatment*ExpDay
               + (ExpDay|Mes_ID)),
            # family = gaussian(link = "identity"),
            family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta = 0.99),
            seed = 543,
            backend = "cmdstanr", 
            data = dat1,
            file = "models/231025_hIR.m1",
            file_refit = "on_change"
)
pp_check(hir.m1)
summary(hir.m1)

conditional_effects(hir.m1, effects = "Treatment:ExpDay")
conditional_effects(hir.m1, effects = "ExpDay:Treatment")

dddd = emmeans(hir.m1, ~ Treatment|ExpDay)
summary(dddd, point = "mean")
pppp = pairs(dddd)
summary(pppp, point = "mean")

dddd1 = emmeans(hir.m1, ~ ExpDay|Treatment)
summary(dddd1, point = "mean")
pppp1 = pairs(dddd1)
summary(pppp1, point = "mean")

# mix grazing
mgr.m1 = brm(bf(MF_Gr ~ 
                 + Treatment*ExpDay
               + (ExpDay|Mes_ID)),
            # family = gaussian(link = "identity"),
            family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta = 0.99),
            seed = 543,
            backend = "cmdstanr", 
            data = dat1,
            file = "models/231025_mGR.m1",
            # file_refit = "on_change"
)
conditional_effects(mgr.m1, effects = "Treatment:ExpDay")
conditional_effects(mgr.m1, effects = "ExpDay:Treatment")

dddd = emmeans(mgr.m1, ~ Treatment|ExpDay)
summary(dddd, point = "mean")
pppp = pairs(dddd)
summary(pppp, point = "mean")

dddd1 = emmeans(mgr.m1, ~ ExpDay|Treatment)
summary(dddd1, point = "mean")
pppp1 = pairs(dddd1)
summary(pppp1, point = "mean")

# het ingestion
hgr.m1 = brm(bf(HF_Gr ~ 
                 + Treatment*ExpDay
               + (ExpDay|Mes_ID)),
            # family = gaussian(link = "identity"),
            family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta = 0.99),
            seed = 123,
            backend = "cmdstanr", 
            data = dat1,
            file = "models/231025_hGR.m1",
            # file_refit = "on_change"
)
pp_check(hgr.m1)
summary(hgr.m1)

conditional_effects(hgr.m1, effects = "Treatment:ExpDay")
conditional_effects(hgr.m1, effects = "ExpDay:Treatment")

dddd = emmeans(hgr.m1, ~ Treatment|ExpDay)
summary(dddd, point = "mean")
pppp = pairs(dddd)
summary(pppp, point = "mean")

dddd1 = emmeans(hgr.m1, ~ ExpDay|Treatment)
summary(dddd1, point = "mean")
pppp1 = pairs(dddd1)
summary(pppp1, point = "mean")

## Sizes -----

ps.m = brm(bf(biovol_PF ~ 
               + Treatment*ExpDay
             + (ExpDay|Mes_ID)),
          family = gaussian(link = "identity"),
          chains = 4,
          iter = 2000,
          cores = 4,
          control = list(adapt_delta = 0.99),
          seed = 543,
          backend = "cmdstanr", 
          data = dat,
          file = "models/231211_phot-biovol.m"
)
pp_check(ps.m)
summary(ps.m)
conditional_effects(ps.m, effects = "Treatment:ExpDay")
conditional_effects(ps.m, effects = "ExpDay:Treatment")

dddd = emmeans(ps.m, ~ Treatment|ExpDay)
summary(dddd, point = "mean")
pppp = pairs(dddd)
summary(pppp, point = "mean")

dddd1 = emmeans(ps.m, ~ ExpDay|Treatment)
summary(dddd1, point = "mean")
pppp1 = pairs(dddd1)
summary(pppp1, point = "mean")


hs.m = brm(bf(biovol_HF ~ 
               + Treatment*ExpDay
             + (ExpDay|Mes_ID)),
          family = gaussian(link = "identity"),
          chains = 4,
          iter = 2000,
          cores = 4,
          control = list(adapt_delta = 0.99),
          seed = 543,
          backend = "cmdstanr", 
          data = dat,
          file = "models/231211_het-biovol.m"
)
pp_check(hs.m)
summary(hs.m)
conditional_effects(hs.m, effects = "Treatment:ExpDay")
conditional_effects(hs.m, effects = "ExpDay:Treatment")

dddd = emmeans(hs.m, ~ Treatment|ExpDay)
summary(dddd, point = "mean")
pppp = pairs(dddd)
summary(pppp, point = "mean")

dddd1 = emmeans(hs.m, ~ ExpDay|Treatment)
summary(dddd1, point = "mean")
pppp1 = pairs(dddd1)
summary(pppp1, point = "mean")



ms.m1 = brm(bf(biovol_MF ~ 
               + Treatment*ExpDay
             + (ExpDay|Mes_ID)),
          family = gaussian(link = "identity"),
          chains = 4,
          iter = 2000,
          cores = 4,
          control = list(adapt_delta = 0.99),
          seed = 543,
          backend = "cmdstanr", 
          data = dat1,
          file = "models/231211_mix1-biovol.m"
)
pp_check(ms.m1)
summary(ms.m1)
conditional_effects(ms.m1, effects = "Treatment:ExpDay")
conditional_effects(ms.m1, effects = "ExpDay:Treatment")

dddd = emmeans(ms.m1, ~ Treatment|ExpDay)
summary(dddd, point = "mean")
pppp = pairs(dddd)
summary(pppp, point = "mean")

dddd1 = emmeans(ms.m1, ~ ExpDay|Treatment)
summary(dddd1, point = "mean")
pppp1 = pairs(dddd1)
summary(pppp1, point = "mean")

