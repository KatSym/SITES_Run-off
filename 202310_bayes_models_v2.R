
library(tidyverse)
library(brms)
library(tidybayes)
library(modelr)
library(emmeans)
library(loo)
library(easystats)
library(rempsyc)
library(ggpubr)


# all models in Bayes

# data
load("all_data.RData")

mysc <- function(x){
  (x-mean(x))/sd(x)
}
# 
# trt.cols2 <-  c(`C`= "#1a1a1a", #black - C
#                 `D`= "#20854e", #g - D
#                 `I`= "#0072b5", #blu - I
#                 `E`= "#bc3c29") #r - E
# all days 
# no transformation
env <- envir %>% 
  select(-c(TOC, DN, Chla, DOsat, DOconc)) %>% 
  relocate(c(ExpDay, Treatment), .before = Mes_ID) %>% 
  mutate(Nutr = TN*1000 + TP,
         # across(DOC:Nutr, scale),
         # scale all predisctors
         # across(DOC:Nutr, mysc),
                # ~(.-mean(.)) / sd(.)
                ExpDay = factor(ExpDay)
         )
  
dat <- data %>% 
  # log(x + 1) transormed
  relocate(c(M.Ir, H.Ir, M.Gr, H.Gr), .after = biovol_MF) %>% 
    mutate(across(4:16, log1p),
           ExpDay = factor(ExpDay
                           # , levels = c("5", "13", "21")
           )
) 
# take out day 5
dat1 <- dat %>% filter(ExpDay !="5") %>% droplevels() %>% 
  mutate(Treatment = as.factor(Treatment))

zoop <- full_join(rot1, zoop1, by = c("ExpDay", "Treatment", "Mes_ID")) %>% 
  ungroup() %>% 
  mutate(ExpDay = as.factor(ExpDay))

# dat2 <- dat %>% select(ExpDay, Treatment, PF_abund, MF_abund, HF_abund) %>% 
#   group_by(ExpDay,Treatment) %>% 
#   summarize_all(mean)
# 
# dat2$Ror = !is.na(env$Rot)

# models:

# # and then
# ap.m1 = brm(bf(PF_abund ~ 
#                 # 0 + 
#                 Nutr*PAR*ExpDay
#               + (1|Mes_ID)),
#            family = gaussian(link = "identity"),
#            chains = 4,
#            iter = 2000,
#            cores = 4,
#            control = list(adapt_delta = 0.95),
#            seed = 543,
#            backend = "cmdstanr", 
#            data = dat
# )
# summary(ap.m1)
# 
abs.m = brm(bf(a420 ~ 
                 + Treatment*ExpDay
               + (1|Mes_ID)),
            family = gaussian(link = "identity"),
            # family = mix,
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



pr <- get_prior(PAR ~ 
                  + Treatment*ExpDay
                + (1|Mes_ID), data = env, family = gaussian())

par.m = brm(bf(PAR ~ 
                + Treatment*ExpDay
              + (1|Mes_ID)),
           family = gaussian(link = "identity"),
           # family = mix,
           # prior = pr,
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

par <- report(par.m)
as.data.frame(par)

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

par.m %>% 
  spread_draws(condition_mean[condition]) %>%
  head(10)

sjPlot::tab_model(par.m, abs.m, nut.m)

# NOT GOOD

dpr <- get_prior(DOC ~ 
                   + Treatment*ExpDay
                 + (1|Mes_ID),
                 family = gaussian(link = "identity"),data = env)


doc.m = brm(bf(DOC ~ 
                + Treatment*ExpDay
              + (1|Mes_ID)),
           family = gaussian(link = "identity"),
           # family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
           chains = 4,
           iter = 2000,
           # prior = dpr,
           cores = 4,
           control = list(adapt_delta = 0.99, 
                          max_treedepth = 15),
           seed = 543,
           backend = "cmdstanr", 
           data = env%>% filter(!is.na(DOC)),
           file = "models/231026_doc",
           file_refit = "on_change"
)
pp_check(doc.m)
summary(doc.m)
conditional_effects(doc.m, effects = "Treatment:ExpDay")
conditional_effects(doc.m, effects = "ExpDay:Treatment")

emTrDay = emmeans(doc.m, ~ Treatment|ExpDay)
summary(emTrDay, point = "mean")
pTrDay = pairs(emTrDay)
summary(pTrDay, point = "mean")

emDayTr = emmeans(doc.m, ~ ExpDay|Treatment)
summary(emDayTr, point = "mean")
pDayTr = pairs(emDayTr)
summary(pDayTr, point = "mean")


nut.m = brm(bf(Nutr ~ 
                 + Treatment*ExpDay
               + (1|Mes_ID)),
            family = gaussian(link = "identity"),
            # family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
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
            file_refit = "on_change"
)
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
            # family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta = 0.99,
                           max_treedepth = 12,
                           step_size = 0.2),
            seed = 543,
            backend = "cmdstanr", 
            data = env%>% filter(!is.na(TN)),
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
           # family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
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



rot.m = brm(bf(log1p(Rotif) ~ 
                 + Treatment*ExpDay
               + (1|Mes_ID)),
            family = gaussian(link = "identity"),
            # family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
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




zoo.m = brm(bf(log1p(Zoopl) ~
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
            data = zoop %>% filter(!is.na(Zoopl)),
            file = "models/231127_zoopl",
            file_refit = "on_change"
)

pp_check(zoo.m, ndraws = 10)
summary(zoo.m)
conditional_effects(zoo.m, effects = "Treatment:ExpDay")
conditional_effects(zoo.m, effects = "ExpDay:Treatment")

emTrDay = emmeans(zoo.m, ~ Treatment|ExpDay)
summary(emTrDay, point = "mean")
pTrDay = pairs(emTrDay)

emDayTr = emmeans(zoo.m, ~ ExpDay|Treatment)
summary(emDayTr, point = "mean")
pDayTr = pairs(emDayTr)

summary(pTrDay, point = "mean")
summary(pDayTr, point = "mean")





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
          file = "models/231014_bacy.m"
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
          file = "models/231030_cy.m"
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
           file = "models/231014_phot.m"
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
           file = "models/231014_het.m"
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


m.m = brm(bf(MF_abund ~ 
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
          data = dat,
          file = "models/231014_mix.m"
)

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
          file = "models/231014_mix1.m"
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


## Ingestion rates

#without transformation

mir.m1 = brm(bf(M.Ir ~ 
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


hir.m1 = brm(bf(H.Ir ~ 
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


mgr.m1 = brm(bf(M.Gr ~ 
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

hgr.m1 = brm(bf(H.Gr ~ 
                 + Treatment*ExpDay
               + (ExpDay|Mes_ID)),
            # family = gaussian(link = "identity"),
            family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
            # family = hurdle_gamma(link = "log"),
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

## Sizes =====

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
          file = "models/231102_phot-biovol.m"
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
          file = "models/231102_het-biovol.m"
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

ms.m = brm(bf(biovol_MF ~ 
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
          file = "models/231102_mix-biovol.m"
)
pp_check(ms.m)
summary(ms.m)
conditional_effects(ms.m, effects = "Treatment:ExpDay")
conditional_effects(ms.m, effects = "ExpDay:Treatment")

dddd = emmeans(ms.m, ~ Treatment|ExpDay)
summary(dddd, point = "mean")
pppp = pairs(dddd)
summary(pppp, point = "mean")

dddd1 = emmeans(ms.m, ~ ExpDay|Treatment)
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
          file = "models/231102_mix1-biovol.m"
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



# Plots ======

matheme <- theme(
  legend.position = "none",
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black", linewidth = .3),
  axis.text = element_text(size = 11, colour = "black"),
  axis.title = element_text(size = 13),
  # axis.text.y = element_text(size = 12),
  # strip.text.y = element_text(size = 12)
) 

### facet plot------
bio_lst = list(p.m, h.m, m.m)
bio <- lapply(bio_lst, function(a) dat %>%
                group_by(Mes_ID,ExpDay, Treatment) %>% 
                add_epred_draws(a,
                                re_formula = NA,
                                # ndraws = 100
                ))
ab.df <- map_dfr(bio, ~ as.data.frame(.x), .id = "id") %>% 
  mutate(group = case_when(id == 1 ~ "Phototroph",
                           id == 2 ~ "Heterotroph",
                           id == 3 ~ "Mixotroph"))


abdatt <- dat %>% 
  select(ExpDay, Treatment, Mes_ID, PF_abund, HF_abund, MF_abund) %>% 
  pivot_longer(cols = c(PF_abund, HF_abund, MF_abund), names_to = "group", values_to = "abund") %>% 
  mutate(group = case_when(group == "PF_abund" ~ "Phototroph",
                           group == "HF_abund" ~ "Heterotroph",
                           group == "MF_abund" ~ "Mixotroph")) %>% 
  mutate(group = fct_relevel(group, c("Phototroph", "Heterotroph", "Mixotroph")))

ab.df %>% 
  ggplot(., aes(x = ExpDay,
                y = biomass,
                colour = Treatment)) +
  facet_grid(cols = vars(group),
             scales = "free_y") +
  stat_pointinterval(aes(y = (.epred)),
                     .width = c(0.95),
                     position = position_dodge(.5),
                     fatten_point = 3)+
  geom_point(data = abdatt, 
             aes(x = ExpDay, 
                 y = abund, 
                 colour = Treatment), 
             position = position_jitterdodge(dodge.width = .5),
             alpha = .5) +
  # scale_x_continuous(breaks = c(1, 2, 3), 
  #                    labels = c(5, 13, 21)) +
  # scale_y_continuous(labels = scales::label_number(scale = 1/1000))+
  scale_color_manual(values = trt.cols)+
  # scale_fill_manual(values = c("#E5E5E5", "#D6E6E5","#EBF0F9",  "#FCECEE")) +
matheme +
  labs(y = "Log(x+1) cells/mL",
       x = "Experimental day")

##----

# mod.plot <- function(mod){
#   dat %>% 
#     group_by(Mes_ID,ExpDay, Treatment) %>% 
#     # data_grid(ExpDay = seq_range(ExpDay, n = 101)) %>% 
#     add_epred_draws(mod,
#                     re_formula = NA,
#                     # ndraws = 100
#     ) %>% 
#     ggplot(., aes(x = ExpDay,
#                   y = mod$data[,1],
#                   colour = Treatment,
#                   # fill = Treatment
#     )) +
#     # geom_boxplot(aes(y = (.epred)), width= .7)+
#     # geom_violinhalf(aes(y = (.epred)))+
#     stat_pointinterval(aes(y = (.epred)),
#                        .width = c(0.95),
#                        position = position_dodge(.5),
#                        fatten_point = 3)+
#     geom_point(
#       data = dat,
#       aes(x = ExpDay,
#           y = PF_abund,
#           colour = Treatment),
#       inherit.aes = FALSE,
#       position = position_jitterdodge(dodge.width = .5),
#       alpha = .5) +
#     # scale_x_continuous(breaks = c(1, 2, 3), 
#     #                    labels = c(5, 13, 21)) +
#     scale_color_manual(values = trt.cols)+
#     scale_fill_manual(values = trt.cols)+
#     #                    aesthetics = c("colour")) +
#     # scale_fill_manual(values = c("#E5E5E5", "#D6E6E5","#EBF0F9",  "#FCECEE")) +
#     # theme_modern()
#     theme(panel.grid.minor = element_blank(),
#           panel.grid.major = element_blank(),
#           # panel.background = element_rect(fill = "grey98"),
#           panel.background = element_blank(),
#           axis.line = element_line(colour = "black", linewidth = .3),
#           axis.text = element_text(size = 11, colour = "black"),
#           axis.title = element_text(size = 11),
#           # axis.text.y = element_text(size = 12),
#           strip.text.y = element_text(size = 12)) +
#     labs(y = "Log(x+1) Abundance cells/mL",
#          x = "Experimental day")
# }
# 
# 
# plot_list <- map(bio_lst, mod.plot)
# 
# plt_df <- tibble(variable = bio_lst, plots = plot_list)

## backgr====

(par.pl <- env %>% 
   filter(!is.na(PAR)) %>% droplevels(ExpDay) %>% 
   group_by(Mes_ID,ExpDay, Treatment) %>% 
   add_epred_draws(par.m,
                   re_formula = NA,
   ) %>% 
   ggplot(., aes(x = ExpDay,
                 y = PAR,
                 colour = Treatment,
   )) +
   # geom_vline(aes(xintercept = 1.5), linetype = "dashed", color = "#e84855", alpha = 0.4)+
   stat_pointinterval(aes(y = (.epred)),
                      .width = c(0.95),
                      position = position_dodge(.5),
                      # fatten_point = 3,
                      linewidth = 2, 
                      show.legend = FALSE
   )+
   geom_point(
     data = env,
     aes(x = ExpDay,
         y = PAR,
         colour = Treatment,
         # shape = Treatment
     ),
     inherit.aes = FALSE,
     position = position_jitterdodge(dodge.width = .5),
     alpha = .35
   ) +
   scale_color_manual(values = trt.cols)+
   scale_fill_manual(values = trt.cols)+
   # theme_modern()+
   theme(
     # legend.position = "none",
     legend.title = element_blank(),
     legend.background	= element_blank(),
     legend.key	= element_blank(),
     legend.direction = "horizontal",
     panel.grid.minor = element_blank(),
     panel.grid.major = element_blank(),
     panel.background = element_blank(),
     axis.line = element_line(colour = "black", linewidth = .3),
     axis.text = element_text(size = 12, colour = "black"),
     axis.title = element_text(size = 14),
     # axis.text.y = element_text(size = 12),
     # strip.text.y = element_text(size = 12)
   ) +
   labs(y = "Log(x+1) cells/mL",
        x = "Experimental day")+
   guides(colour = guide_legend(override.aes = list(size=4,
                                                    color = trt.cols)))
)




## abundances-----
(phot.leg <- dat %>% 
    group_by(Mes_ID,ExpDay, Treatment) %>% 
    add_epred_draws(p.m,
                    re_formula = NA,
    ) %>% 
    ggplot(., aes(x = ExpDay,
                  y = PF_abund,
                  colour = Treatment,
    )) +
   # geom_vline(aes(xintercept = 1.5), linetype = "dashed", color = "#e84855", alpha = 0.4)+
    stat_pointinterval(aes(y = (.epred)),
                       .width = c(0.95),
                       position = position_dodge(.5),
                       # fatten_point = 3,
                       linewidth = 2, 
                       show.legend = FALSE
    )+
    geom_point(
      data = dat,
      aes(x = ExpDay,
          y = PF_abund,
          colour = Treatment,
          # shape = Treatment
          ),
      inherit.aes = FALSE,
      position = position_jitterdodge(dodge.width = .5),
      alpha = .35
      ) +
    scale_color_manual(values = trt.cols)+
    scale_fill_manual(values = trt.cols)+
    # theme_modern()+
    theme(
      # legend.position = "none",
      legend.title = element_blank(),
      legend.background	= element_blank(),
      legend.key	= element_blank(),
      legend.direction = "horizontal",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", linewidth = .3),
      axis.text = element_text(size = 12, colour = "black"),
      axis.title = element_text(size = 14),
      # axis.text.y = element_text(size = 12),
      # strip.text.y = element_text(size = 12)
    ) +
    labs(y = "Log(x+1) cells/mL",
         x = "Experimental day")+
    guides(colour = guide_legend(override.aes = list(size=4,
                                                     color = trt.cols)))
)

leg <- get_legend(phot.leg)
# ggsave("Plots/legend2.png", leg, dpi = 300)
phot <- phot.leg + theme(legend.position = "none")

(het <- dat %>% 
    group_by(Mes_ID,ExpDay, Treatment) %>% 
    add_epred_draws(h.m,
                    re_formula = NA,
    ) %>% 
    ggplot(., aes(x = ExpDay,
                  y = HF_abund,
                  colour = Treatment,
    )) +
    stat_pointinterval(aes(y = (.epred)),
                       .width = c(0.95),
                       position = position_dodge(.5),
                       # fatten_point = 3,
                       linewidth = 2
    )+
    geom_point(
      data = dat,
      aes(x = ExpDay,
          y = HF_abund,
          colour = Treatment),
      inherit.aes = FALSE,
      position = position_jitterdodge(dodge.width = .5),
      alpha = .35) +
    scale_color_manual(values = trt.cols)+
    scale_fill_manual(values = trt.cols)+
    # theme_modern()+
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", linewidth = .3),
      axis.text = element_text(size = 12, colour = "black"),
      axis.title = element_text(size = 14),
      # axis.text.y = element_text(size = 12),
      # strip.text.y = element_text(size = 12)
    ) +
    labs(y = NULL,
         x = "Experimental day")
)


(mix1 <- dat1 %>% 
  group_by(Mes_ID,ExpDay, Treatment) %>% 
  add_epred_draws(m.m1,
                  re_formula = NA,
  ) %>% 
  ggplot(., aes(x = ExpDay,
                y = MF_abund,
                colour = Treatment,
  )) +
  stat_pointinterval(aes(y = (.epred)),
                     .width = c(0.95),
                     position = position_dodge(.5),
                     # fatten_point = 3,
                     linewidth = 2
  )+
  geom_point(
    data = dat1,
    aes(x = ExpDay,
        y = MF_abund,
        colour = Treatment),
    inherit.aes = FALSE,
    position = position_jitterdodge(dodge.width = .5),
    alpha = .35) +
  scale_color_manual(values = trt.cols)+
  scale_x_discrete(breaks = c("5", "13", "21"),
                   limits = c("5", "13", "21"))+
  # theme_modern()+
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = .3),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14),
    # axis.text.y = element_text(size = 12),
    # strip.text.y = element_text(size = 12)
  ) +
  labs(y = NULL,
       x = "Experimental day")
)




p1 <- ggarrange(phot, mix1, het, ncol = 3
                # ,common.legend = T, legend.grob = leg, legend = "bottom"
                )+ bgcolor("white") 

ggsave("Plots/Nov2023/abund_20231103.tiff", p1, dpi=300) 




(bact <- dat %>% 
    group_by(Mes_ID,ExpDay, Treatment) %>% 
    add_epred_draws(bact.m,
                    re_formula = NA,
    ) %>% 
    ggplot(., aes(x = ExpDay,
                  y = HB_abund,
                  colour = Treatment,
    )) +
    stat_pointinterval(aes(y = (.epred)),
                       .width = c(0.95),
                       position = position_dodge(.5),
                       # fatten_point = 3,
                       linewidth = 2
    )+
    geom_point(
      data = dat,
      aes(x = ExpDay,
          y = HB_abund,
          colour = Treatment),
      inherit.aes = FALSE,
      position = position_jitterdodge(dodge.width = .5),
      alpha = .35) +
    scale_color_manual(values = trt.cols)+
    # theme_modern()+
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", linewidth = .3),
      axis.text = element_text(size = 12, colour = "black"),
      axis.title = element_text(size = 14),
      # axis.text.y = element_text(size = 12),
      # strip.text.y = element_text(size = 12)
    ) +
    labs(y = "Log(x+1) cells mL\u207b\u00b9",
         x = "Experimental day")
)
  (cyan <- dat %>% 
      group_by(Mes_ID,ExpDay, Treatment) %>% 
      add_epred_draws(cy.m,
                      re_formula = NA,
      ) %>% 
      ggplot(., aes(x = ExpDay,
                    y = CY_abund,
                    colour = Treatment,
      )) +
      stat_pointinterval(aes(y = (.epred)),
                         .width = c(0.95),
                         position = position_dodge(.5),
                         # fatten_point = 3,
                         linewidth = 2
      )+
      geom_point(
        data = dat,
        aes(x = ExpDay,
            y = CY_abund,
            colour = Treatment),
        inherit.aes = FALSE,
        position = position_jitterdodge(dodge.width = .5),
        alpha = .4) +
      scale_color_manual(values = trt.cols)+

      # theme_modern()+
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", linewidth = .3),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14),
        # axis.text.y = element_text(size = 12),
        # strip.text.y = element_text(size = 12)
      ) +
      labs(y = NULL,
           x = "Experimental day")
)

p2 <- ggarrange(bact, cyan, align = "h") + bgcolor("white")
ggsave("Plots/Nov2023/bact_20231030.tiff", p2, dpi=300)



env %>% 
  group_by(Mes_ID,ExpDay, Treatment) %>% 
  add_epred_draws(Rot,
                  re_formula = NA,
  ) %>% 
  ggplot(., aes(x = ExpDay,
                y = rot.m,
                colour = Treatment,
  )) +
  stat_pointinterval(aes(y = (.epred)),
                     .width = c(0.95),
                     position = position_dodge(.5),
                     # fatten_point = 3,
                     linewidth = 2
  )+
  geom_point(
    data = env,
    aes(x = ExpDay,
        y = log1p(Rot),
        colour = Treatment),
    inherit.aes = FALSE,
    position = position_jitterdodge(dodge.width = .5),
    alpha = .5) +
  scale_color_manual(values = trt.cols)+
  scale_fill_manual(values = trt.cols)+
  # theme_modern()+
  matheme+
  theme(    axis.title.y = ggtext::element_markdown())+
  #   # strip.text.y = element_text(size = 12)
  # ) +
  labs(y = NULL,
       x = "Experimental day")


## rates ----

(mir <- dat1 %>% 
    group_by(Mes_ID,ExpDay, Treatment) %>% 
    add_epred_draws(mir.m1,
                    re_formula = NA,
    ) %>% 
    ggplot(., aes(x = ExpDay,
                  y = M.Ir,
                  colour = Treatment,
    )) +
    stat_pointinterval(aes(y = (.epred)),
                       .width = c(0.95),
                       position = position_dodge(.5),
                       # fatten_point = 3,
                       linewidth = 2
    )+
    geom_point(
      data = dat1,
      aes(x = ExpDay,
          y = M.Ir,
          colour = Treatment),
      inherit.aes = FALSE,
      position = position_jitterdodge(dodge.width = .5),
      alpha = .5) +
    scale_color_manual(values = trt.cols)+
    scale_fill_manual(values = trt.cols)+
    # theme_modern()+
matheme+
    labs(y = "Ingestion rate (FLB cell\u207b\u00b9 h\u207b\u00b9)", 
      # y =paste0("<span style='font-size: 13pt'>Ingestion rate </span><span style='font-size: 11pt'>(FLB cell\u207b\u00b9 h\u207b\u00b9)</span>"),
         x = NULL))+
    theme(
      axis.title.y = ggtext::element_markdown(),
    #   # strip.text.y = element_text(size = 12)
    ) 

(mgr <- dat1 %>% 
    group_by(Mes_ID,ExpDay, Treatment) %>% 
    add_epred_draws(mgr.m1,
                    re_formula = NA,
    ) %>% 
    ggplot(., aes(x = ExpDay,
                  y = M.Gr,
                  colour = Treatment,
    )) +
    stat_pointinterval(aes(y = (.epred)),
                       .width = c(0.95),
                       position = position_dodge(.5),
                       # fatten_point = 3,
                       linewidth = 2
    )+
    geom_point(
      data = dat1,
      aes(x = ExpDay,
          y = M.Gr,
          colour = Treatment),
      inherit.aes = FALSE,
      position = position_jitterdodge(dodge.width = .5),
      alpha = .5) +
    scale_color_manual(values = trt.cols)+
    scale_fill_manual(values = trt.cols)+
    # theme_modern()+
matheme+
    labs(
      # y = paste0("<span style='font-size: 13pt'>Grazing rate </span><span style='font-size: 11pt'>(bacteria cell\u207b\u00b9 h\u207b\u00b9)</span>"),
      y = "Grazing rate (bacteria cell\u207b\u00b9 h\u207b\u00b9)",
         x = "Experimental day"))+
    theme(
      axis.title.y = ggtext::element_markdown(),
    #   # strip.text.y = element_text(size = 12)
    ) 

(hir <- dat1 %>% 
    group_by(Mes_ID,ExpDay, Treatment) %>% 
    add_epred_draws(hir.m1,
                    re_formula = NA,
    ) %>% 
    ggplot(., aes(x = ExpDay,
                  y = H.Ir,
                  colour = Treatment,
    )) +
    stat_pointinterval(aes(y = (.epred)),
                       .width = c(0.95),
                       position = position_dodge(.5),
                       # fatten_point = 3,
                       linewidth = 2
    )+
    geom_point(
      data = dat1,
      aes(x = ExpDay,
          y = H.Ir,
          colour = Treatment),
      inherit.aes = FALSE,
      position = position_jitterdodge(dodge.width = .5),
      alpha = .5) +
    scale_color_manual(values = trt.cols)+
    scale_fill_manual(values = trt.cols)+
    # theme_modern()+
matheme+
    theme(    axis.title.y = ggtext::element_markdown())+
    #   # strip.text.y = element_text(size = 12)
    # ) +
    labs(y = NULL,
         x = NULL))

(hgr <- dat1 %>% 
    group_by(Mes_ID,ExpDay, Treatment) %>% 
    add_epred_draws(hgr.m1,
                    re_formula = NA,
    ) %>% 
    ggplot(., aes(x = ExpDay,
                  y = H.Gr,
                  colour = Treatment,
    )) +
    stat_pointinterval(aes(y = (.epred)),
                       .width = c(0.95),
                       position = position_dodge(.5),
                       # fatten_point = 3,
                       linewidth = 2
    )+
    geom_point(
      data = dat1,
      aes(x = ExpDay,
          y = H.Gr,
          colour = Treatment),
      inherit.aes = FALSE,
      position = position_jitterdodge(dodge.width = .5),
      alpha = .5) +
    scale_color_manual(values = trt.cols)+
    scale_fill_manual(values = trt.cols)+
    # theme_modern()+
matheme+
    theme(    axis.title.y = ggtext::element_markdown())+
    #   # strip.text.y = element_text(size = 12)
    # ) +
    labs(y = NULL,
         x = "Experimental day"))

library(ggpubr)
p3 <- ggarrange(mir, hir, mgr, hgr, align = "hv") + bgcolor("white") 
  theme(    axis.title.y = ggtext::element_markdown())

ggsave("Plots/rates_2_20231101.png",p3, dpi = 300)
ggsave("Plots/Mir_20231101.png",mir, height = 3.63, width = 3.94, units = "in", dpi = 300)



plot_grid(mir, hir, mgr, hgr, align = "hv")
library(patchwork)
mir+ hir+ mgr+ hgr

library(ragg)
agg_tiff("Plots/rates_2_20231101.tiff",p3, width=7.88, height=7.26,units = "in", res = 72)

## biovolume ----

(phot.vol <- dat %>% 
   group_by(Mes_ID,ExpDay, Treatment) %>% 
   add_epred_draws(ps.m,
                   re_formula = NA,
   ) %>% 
   ggplot(., aes(x = ExpDay,
                 y = biovol_PF,
                 colour = Treatment,
   )) +
   stat_pointinterval(aes(y = (.epred)),
                      .width = c(0.95),
                      position = position_dodge(.5),
                      # fatten_point = 3,
                      linewidth = 2
   )+
   geom_point(
     data = dat,
     aes(x = ExpDay,
         y = biovol_PF,
         colour = Treatment),
     inherit.aes = FALSE,
     position = position_jitterdodge(dodge.width = .5),
     alpha = .4) +
   scale_color_manual(values = trt.cols)+
   scale_fill_manual(values = trt.cols)+
   # theme_modern()+
   theme(
     legend.position = "none",
     panel.grid.minor = element_blank(),
     panel.grid.major = element_blank(),
     panel.background = element_blank(),
     axis.line = element_line(colour = "black", linewidth = .3),
     axis.text = element_text(size = 11, colour = "black"),
     axis.title = element_text(size = 13),
     # axis.text.y = element_text(size = 12),
     # strip.text.y = element_text(size = 12)
   ) +
   labs(y = "Biovolume log(x+1) \u00b5m\u00b3 mL\u207b\u00b9",
        x = "Experimental day")
 )


(het.vol <- dat %>% 
   group_by(Mes_ID,ExpDay, Treatment) %>% 
   add_epred_draws(hs.m,
                   re_formula = NA,
   ) %>% 
   ggplot(., aes(x = ExpDay,
                 y = biovol_HF,
                 colour = Treatment,
   )) +
   stat_pointinterval(aes(y = (.epred)),
                      .width = c(0.95),
                      position = position_dodge(.5),
                      # fatten_point = 3,
                      linewidth = 2
   )+
   geom_point(
     data = dat,
     aes(x = ExpDay,
         y = biovol_HF,
         colour = Treatment),
     inherit.aes = FALSE,
     position = position_jitterdodge(dodge.width = .5),
     alpha = .4) +
   scale_color_manual(values = trt.cols)+
   scale_fill_manual(values = trt.cols)+
   # theme_modern()+
   theme(
     legend.position = "none",
     panel.grid.minor = element_blank(),
     panel.grid.major = element_blank(),
     panel.background = element_blank(),
     axis.line = element_line(colour = "black", linewidth = .3),
     axis.text = element_text(size = 11, colour = "black"),
     axis.title = element_text(size = 13),
     # axis.text.y = element_text(size = 12),
     # strip.text.y = element_text(size = 12)
   ) +
   labs(y = NULL,
        x = "Experimental day")
)

(mix.vol1 <- dat1 %>% 
   group_by(Mes_ID,ExpDay, Treatment) %>% 
   add_epred_draws(ms.m,
                   re_formula = NA,
   ) %>% 
   ggplot(., aes(x = ExpDay,
                 y = biovol_MF,
                 colour = Treatment,
   )) +
   stat_pointinterval(aes(y = (.epred)),
                      .width = c(0.95),
                      position = position_dodge(.5),
                      # fatten_point = 3,
                      linewidth = 2
   )+
   geom_point(
     data = dat1,
     aes(x = ExpDay,
         y = biovol_MF,
         colour = Treatment),
     inherit.aes = FALSE,
     position = position_jitterdodge(dodge.width = .5),
     alpha = .4) +
   scale_color_manual(values = trt.cols)+
   scale_x_discrete(breaks = c("5", "13", "21"),
                    limits = c("5", "13", "21"))+
   # theme_modern()+
   theme(
     legend.position = "none",
     panel.grid.minor = element_blank(),
     panel.grid.major = element_blank(),
     panel.background = element_blank(),
     axis.line = element_line(colour = "black", linewidth = .3),
     axis.text = element_text(size = 11, colour = "black"),
     axis.title = element_text(size = 13),
     # axis.text.y = element_text(size = 12),
     # strip.text.y = element_text(size = 12)
   ) +
   labs(y = NULL,
        x = "Experimental day")
)

p4 <- ggarrange(phot.vol, mix.vol1, het.vol, nrow = 1, align = "h")

ggsave("Plots/Nov2023/20231103_biovolumes.tiff", p4, dpi = 300)


