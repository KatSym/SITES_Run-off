
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
    mutate(across(4:16, mytr),
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

leg.v <- ggplot(dat, aes(x = ExpDay, y = PF_abund, color = Treatment))+
  geom_point()+
  lims(y = c(0,0))+
  scale_color_manual(values = trt.cols)+
  theme_void()+
  theme(legend.position = c(0.5,0.5),
        legend.text = element_text(size =  11),
        legend.background	= element_blank(),
        legend.key	= element_blank(),
        legend.title = element_blank()
        )+
  guides(colour = guide_legend(override.aes = list(size=4))) # or 4

leg.h <- ggplot(dat, aes(x = ExpDay, y = PF_abund, color = Treatment))+
  geom_point()+
  lims(y = c(0,0))+
  scale_color_manual(values = trt.cols)+
  theme_void()+
  theme(legend.position = c(0.5,0.5),
        legend.text = element_text(size =  11),
        legend.direction = "horizontal",
        legend.background	= element_blank(),
        legend.key	= element_blank(),
        legend.title = element_blank()
  )+
  guides(colour = guide_legend(override.aes = list(size=4))) # or 4
#--

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

nutr = envir %>% 
  select(-c(TOC, DN, DOsat, DOconc)) %>% 
  relocate(c(ExpDay, Treatment), .before = Mes_ID)


par.pl <- nutr %>% 
    filter(!is.na(PAR),
           ExpDay != 0) %>% 
   group_by(Mes_ID,ExpDay, Treatment) %>% 
   add_epred_draws(par.m,
                   re_formula = NA,
   ) %>% 
   ggplot(., aes(x = ExpDay,
                 y = PAR,
                 colour = Treatment,
   )) +
   # geom_vline(aes(xintercept = 7), linetype = "dashed", color = "#e84855", alpha = 0.4)+
  # geom_vline(aes(xintercept = 5), linetype = "dashed", color = "black", alpha = 0.4)+
  # geom_vline(aes(xintercept = 13), linetype = "dashed", color = "black", alpha = 0.4)+
  # geom_vline(aes(xintercept = 21), linetype = "dashed", color = "black", alpha = 0.4)+
     stat_lineribbon(aes(y = (.epred)),
                    .width = 0,
                    # point_interval = "mean_hdi",
                    position = position_dodge(.5),
                    size = .8
    )+
   stat_pointinterval(aes(y = (.epred)),
                      .width = c(0.95),
   # position = position_dodge(.5),
   fatten_point = 1.5,
                      linewidth = 1.5,
                      show.legend = FALSE
   )+
   scale_color_manual(values = trt.cols)+
    matheme+
   theme(
     axis.title.y = ggtext::element_markdown(),
   ) +
    labs(y = "<span style='font-size: 13pt'>PAR</span>
         <span style='font-size: 11pt'>(μmol m\u207b\u00b2 s\u207b\u00b9)</span>",
        x = "Experimental day")



(abs.pl <- nutr %>% 
    filter(!is.na(a420)) %>% 
    group_by(Mes_ID,ExpDay, Treatment) %>% 
    add_epred_draws(abs.m,
                    re_formula = NA,
    ) %>% 
    ggplot(., aes(x = ExpDay,
                  y = a420,
                  colour = Treatment,
    )) +
    # geom_vline(aes(xintercept = 7), linetype = "dashed", color = "#e84855", alpha = 0.4)+
    # geom_vline(aes(xintercept = 5), linetype = "dashed", color = "black", alpha = 0.4)+
    # geom_vline(aes(xintercept = 13), linetype = "dashed", color = "black", alpha = 0.4)+
    # geom_vline(aes(xintercept = 21), linetype = "dashed", color = "black", alpha = 0.4)+    # geom_vline(aes(xintercept = 21), linetype = "dashed", color = "black", alpha = 0.4)+
    stat_lineribbon(aes(y = (.epred)),
                    .width = 0,
                    # point_interval = "mean_hdi",
                    position = position_dodge(.5),
                    size = .8
    )+
    stat_pointinterval(aes(y = (.epred)),
                       .width = c(0.95),
                       # position = position_dodge(.5),
                       linewidth = 2,
                       show.legend = FALSE
    )+
    scale_color_manual(values = trt.cols)+
    matheme+
    theme(
      axis.title.y = ggtext::element_markdown(),
    ) +
    labs(y = "<span style='font-size: 13pt'>Abs<sub>420</sub></span>
         <span style='font-size: 11pt'>(m\u207b\u00b9)</span>",
         x = "Experimental day")
)


(tn.pl <- nutr %>% 
    filter(!is.na(TN)) %>% 
    group_by(Mes_ID,ExpDay, Treatment) %>% 
    add_epred_draws(tn.m,
                    re_formula = NA,
    ) %>% 
    ggplot(., aes(x = ExpDay,
                  y = TN,
                  colour = Treatment,
    )) +
    # geom_vline(aes(xintercept = 5), linetype = "dashed", color = "black", alpha = 0.4)+
    # geom_vline(aes(xintercept = 13), linetype = "dashed", color = "black", alpha = 0.4)+
    # geom_vline(aes(xintercept = 21), linetype = "dashed", color = "black", alpha = 0.4)+
    # geom_vline(aes(xintercept = 7), linetype = "dashed", color = "#e84855", alpha = 0.4)+
    stat_lineribbon(aes(y = (.epred)),
                    .width = 0,
                    # point_interval = "mean_hdi",
                    position = position_dodge(.5),
                    size = .8
    )+
    stat_pointinterval(aes(y = (.epred)),
                       .width = c(0.95),
                       # position = position_dodge(.5),
                       linewidth = 2,
                       show.legend = FALSE
    )+
    scale_color_manual(values = trt.cols)+
    matheme+
    theme(
      axis.title.y = ggtext::element_markdown(),
    ) +
    labs(y = "<span style='font-size: 13pt'>Total&nbsp;nitrogen</span>
         <span style='font-size: 11pt'>(mg L\u207b\u00b9)</span>",
         x = NULL)
)

(tp.pl <- nutr %>% 
    filter(!is.na(TP)) %>% 
    group_by(Mes_ID,ExpDay, Treatment) %>% 
    add_epred_draws(tp.m,
                    re_formula = NA,
    ) %>% 
    ggplot(., aes(x = ExpDay,
                  y = TP,
                  colour = Treatment,
    )) +
    # geom_vline(aes(xintercept = 5), linetype = "dashed", color = "black", alpha = 0.4)+
    # geom_vline(aes(xintercept = 13), linetype = "dashed", color = "black", alpha = 0.4)+
    # geom_vline(aes(xintercept = 21), linetype = "dashed", color = "black", alpha = 0.4)+
    # geom_vline(aes(xintercept = 7), linetype = "dashed", color = "#e84855", alpha = 0.4)+
    stat_lineribbon(aes(y = (.epred)),
                    .width = 0,
                    # point_interval = "mean_hdi",
                    position = position_dodge(.5),
                    size = .8
    )+
    stat_pointinterval(aes(y = (.epred)),
                       .width = c(0.95),
                       # position = position_dodge(.5),
                       linewidth = 2,
                       show.legend = FALSE
    )+
    scale_color_manual(values = trt.cols)+
    matheme+
    theme(
      axis.title.y = ggtext::element_markdown(),
    ) +
    labs(y = "<span style='font-size: 13pt'>Total  phosphorus</span>
         <span style='font-size: 11pt'>(\u00b5g L\u207b\u00b9)</span>",
         x = NULL)
)

(chla.pl <- nutr %>% 
    filter(!is.na(Chla)) %>% 
    group_by(Mes_ID,ExpDay, Treatment) %>% 
    add_epred_draws(chla.m,
                    re_formula = NA,
    ) %>% 
    ggplot(., aes(x = ExpDay,
                  y = Chla,
                  colour = Treatment,
    )) +
    # geom_vline(aes(xintercept = 5), linetype = "dashed", color = "black", alpha = 0.4)+
    # geom_vline(aes(xintercept = 13), linetype = "dashed", color = "black", alpha = 0.4)+
    # geom_vline(aes(xintercept = 21), linetype = "dashed", color = "black", alpha = 0.4)+
    # geom_vline(aes(xintercept = 7), linetype = "dashed", color = "#e84855", alpha = 0.4)+
    stat_lineribbon(aes(y = (.epred)),
                    .width = 0,
                    # point_interval = "mean_hdi",
                    position = position_dodge(.5),
                    size = .8
    )+
    stat_pointinterval(aes(y = (.epred)),
                       .width = c(0.95),
                       # position = position_dodge(.5),
                       linewidth = 2,
                       show.legend = FALSE
    )+
    scale_color_manual(values = trt.cols)+
    matheme+
    theme(
      axis.title.y = ggtext::element_markdown(),
    ) +
    labs(y = "<span style='font-size: 13pt'>Chl <i>a</i></span>
         <span style='font-size: 11pt'>(mg L\u207b\u00b9)</span>",
         x = NULL)
)

p.bckg <- ggarrange(treatm, par.pl + rremove("xlab"), abs.pl + rremove("xlab"),  
                    tn.pl, tp.pl, chla.pl, 
                    ncol = 3, nrow = 2, align = "hv",
                    labels = "auto"
                    # common.legend = T, legend.grob = leg.v, legend = "right"
                    ) 

p.bckg <- annotate_figure(p.bckg, 
                          bottom = textGrob("Experimental day", gp = gpar(fontsize = 13)))

p.bckg

ggsave("Plots/Dec2023/background_1812.png", dpi = 300, bg = "white")

## abundances-----
(phot <- dat %>% 
   group_by(Mes_ID,ExpDay, Treatment) %>% 
   add_epred_draws(p.m,
                   re_formula = NA,
   ) %>% 
   ggplot(., aes(x = ExpDay,
                 y = PF_abund,
                 colour = Treatment,
   )) +
   stat_pointinterval(aes(y = (.epred)),
                      .width = c(0.95),
                      position = position_dodge(.5),
                      # fatten_point = 3,
                      linewidth = 2, 
                      show.legend = FALSE
   ) +
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
   scale_color_manual(values = trt.cols) +
   scale_fill_manual(values = trt.cols) +
   labs(title = "Phototroph",
        y = "<span style='font-size: 13pt'>Abundance </span>
         <span style='font-size: 11pt'>log(x+1) cells mL\u207b\u00b9</span>",
        x = NULL)+
   matheme +
   theme(axis.title.y = ggtext::element_markdown(),
         axis.title.x = element_text(size = 13),
         plot.title = element_text(size=13, face="italic")
   ) 
)


het <- dat %>% 
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
    labs(title = "Heterotroph", 
         y = "<span style='font-size: 13pt'>Abundance </span>
         <span style='font-size: 11pt'>Log(x+1) cells mL\u207b\u00b9</span>",
         x = NULL)+
    matheme +
    theme(axis.title.y = ggtext::element_markdown(),
          axis.title.x = element_text(size = 13),
          plot.title = element_text(size=13, face="italic")
    ) 



mix1 <- dat1 %>% 
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
    labs(title = "Mixotroph",
         y = "<span style='font-size: 13pt'>Abundance </span>
         <span style='font-size: 11pt'>Log(x+1) cells mL\u207b\u00b9</span>",
         x = NULL) +
    matheme +
    theme(axis.title.y = ggtext::element_markdown(),
          axis.title.x = element_text(size = 13),
          plot.title = element_text(size=13, face="italic")
    ) 





p1 <- ggarrange(phot, 
                mix1 + rremove("ylab"), 
                het + rremove("ylab"), 
                ncol = 3, align = "h"
                # , common.legend = T, legend.grob = get_legend(leg.v), legend = "right"
                , labels = "auto"
                )
p1

p1.an <- annotate_figure(p1, 
                         bottom = textGrob("Experimental day", 
                                           gp = gpar(fontsize = 13)))
p1.h <- ggarrange(p1.an, leg.h, ncol = 1, align = "h", heights = c(3, 0.2))


ggsave("Plots/Dec2023/abund_20231206-hleg.tiff", p1.h, dpi=300, bg = "white") 




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
    labs(title = "Heterotrophic bacteria",
         y = "<span style='font-size: 13pt'>Abundance </span>
         <span style='font-size: 11pt'>log(x+1) cells mL\u207b\u00b9</span>",
         x = "Experimental Day") +
    matheme +
    theme(axis.title.y = ggtext::element_markdown(),
          axis.title.x = element_text(size = 13),
          plot.title = element_text(size=13, face="italic")
    ) 
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
    labs(title = "Cyanobacteria",
         y = "<span style='font-size: 13pt'>Abundance </span>
         <span style='font-size: 11pt'>log(x+1) cells mL\u207b\u00b9</span>",
         x = NULL) +
    matheme +
    theme( axis.title.y = ggtext::element_markdown(),
           axis.title.x = element_text(size = 13),
           plot.title = element_text(size=13, face="italic")
    ) 
)

p2 <- ggarrange(bact, cyan + rremove("ylab"), align = "h",
                labels = "auto"
                # , common.legend = T, legend.grob = get_legend(leg.v), legend = "right"
                )

p2.an <- annotate_figure(p2, 
                         bottom = textGrob("Experimental day", 
                                           gp = gpar(fontsize = 13)))
p2.h <- ggarrange(p2.an, leg.h, ncol = 1, align = "h", heights = c(3, 0.2))

p2.1 <- ggarrange(cyan,  bact,
                  ncol = 1, align = "v",
                labels = "auto"
                # , common.legend = T, legend.grob = get_legend(leg.h), legend = "bottom"
)



ggsave("Plots/Dec2023/bact-V_20231218-noleg.tiff", p2.1, dpi=300, bg = "white")



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
                 y = MF_Ir,
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
         y = MF_Ir,
         colour = Treatment),
     inherit.aes = FALSE,
     position = position_jitterdodge(dodge.width = .5),
     alpha = .35) +
   scale_color_manual(values = trt.cols)+
   scale_fill_manual(values = trt.cols)+
   # theme_modern()+
   labs(
     title = "Mixotroph",
     y ="<span style='font-size: 13pt'>Ingestion rate</span>
      <span style='font-size: 11pt'>(FLB cell\u207b\u00b9 h\u207b\u00b9)</span>",
     x = NULL)+
  matheme+
  theme(
    axis.title.y = ggtext::element_markdown(),
    axis.title.x = element_text(size = 13),
    plot.title = element_text(size=13, face="italic"))
  ) 

(mgr <- dat1 %>% 
    group_by(Mes_ID,ExpDay, Treatment) %>% 
    add_epred_draws(mgr.m1,
                    re_formula = NA,
    ) %>% 
    ggplot(., aes(x = ExpDay,
                  y = MF_Gr,
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
          y = MF_Gr,
          colour = Treatment),
      inherit.aes = FALSE,
      position = position_jitterdodge(dodge.width = .5),
      alpha = .35) +
    scale_color_manual(values = trt.cols)+
    scale_fill_manual(values = trt.cols)+
    # theme_modern()+
    labs(
      y ="<span style='font-size: 13pt'>Grazing rate</span>
         <span style='font-size: 11pt'>(bacteria cell\u207b\u00b9 h\u207b\u00b9)</span>",
      x = NULL)+
  matheme+
  theme(
    axis.title.y = ggtext::element_markdown(),
    axis.title.x = element_text(size = 13))
)

(hir <- dat1 %>% 
    group_by(Mes_ID,ExpDay, Treatment) %>% 
    add_epred_draws(hir.m1,
                    re_formula = NA,
    ) %>% 
    ggplot(., aes(x = ExpDay,
                  y = HF_Ir,
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
          y = HF_Ir,
          colour = Treatment),
      inherit.aes = FALSE,
      position = position_jitterdodge(dodge.width = .5),
      alpha = .35) +
    scale_color_manual(values = trt.cols)+
    scale_fill_manual(values = trt.cols)+
    # theme_modern()+
    labs(title = "Heterotroph",
         y = NULL,
         x = NULL)+
    matheme+
    theme(plot.title = element_text(size=13, face="italic"))
)

(hgr <- dat1 %>% 
    group_by(Mes_ID,ExpDay, Treatment) %>% 
    add_epred_draws(hgr.m1,
                    re_formula = NA,
    ) %>% 
    ggplot(., aes(x = ExpDay,
                  y = HF_Gr,
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
          y = HF_Gr,
          colour = Treatment),
      inherit.aes = FALSE,
      position = position_jitterdodge(dodge.width = .5),
      alpha = .35) +
    scale_color_manual(values = trt.cols)+
    scale_fill_manual(values = trt.cols)+
    # theme_modern()+
matheme+
    labs(y = NULL,
         x = NULL))

p3 <- ggarrange(mir, hir, mgr, hgr, 
                ncol = 2, nrow = 2,
                align = "hv",
                labels = "auto"
                # , common.legend = T, legend.grob = get_legend(leg.v), legend = "right"
                )
p3.an <- annotate_figure(p3, 
                          bottom = textGrob("Experimental day", 
                                            gp = gpar(fontsize = 13)))
p3.h <- ggarrange(p3.an, leg.h, ncol = 1, align = "h", heights = c(3, 0.2))


ggsave("Plots/Dec2023/rates_20231219-no.png",p3.an, dpi = 300, bg = "white")

# ragg::agg_tiff("Plots/rates_2_20231101.tiff",p3, width=7.88, height=7.26,units = "in", res = 72)

## biovolume ----

phot.vol <- dat %>% 
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
   labs(
     # title = "Phototroph", 
        y = "<span style='font-size: 13pt'>Biovolume </span>
         <span style='font-size: 11pt'>log(x+1) &nbsp;  \u00b5m\u00b3 mL\u207b\u00b9</span>",
        x = NULL)+
   matheme +
   theme(axis.title.y = ggtext::element_markdown(),
         axis.title.x = element_text(size = 13),
         plot.title = element_text(size=13, face="italic")
   ) 
 


het.vol <- dat %>% 
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
    labs(
      # title = "Heterotroph", 
         y = "<span style='font-size: 13pt'>Biovolume </span>
         <span style='font-size: 11pt'>log(x+1) &nbsp;  \u00b5m\u00b3 mL\u207b\u00b9</span>",
         x = NULL)+
    matheme +
    theme(axis.title.y = ggtext::element_markdown(),
          axis.title.x = element_text(size = 13),
          plot.title = element_text(size=13, face="italic")
    ) 


mix.vol1 <- dat1 %>% 
   group_by(Mes_ID,ExpDay, Treatment) %>% 
   add_epred_draws(ms.m1,
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
    labs(
      # title = "Mixotroph", 
         y = "<span style='font-size: 13pt'>Biovolume </span>
         <span style='font-size: 11pt'>log(x+1) &nbsp;  \u00b5m\u00b3 mL\u207b\u00b9</span>",
         x = NULL)+
    matheme +
    theme(axis.title.y = ggtext::element_markdown(),
          axis.title.x = element_text(size = 13),
          plot.title = element_text(size=13, face="italic")
    ) 


p4 <- ggarrange(phot.vol, 
                mix.vol1 + rremove("ylab"), 
                het.vol + rremove("ylab"), 
                ncol = 3, align = "h"
                # , common.legend = T, legend.grob = get_legend(leg.v), legend = "right"
                , labels = "auto"
                )

p4.an <- annotate_figure(p4, 
                         bottom = textGrob("Experimental day", 
                                           gp = gpar(fontsize = 13)))

p4.h <- ggarrange(p4.an, leg.h, ncol = 1, align = "h", heights = c(3, 0.2))


ggsave("Plots/Dec2023/biovol_20231206-hleg.tiff", p4.an, dpi=300, bg = "white") 

# combined abundance and biovolume
p5 <- ggarrange(phot, 
                mix1 + rremove("ylab"), 
                het + rremove("ylab"), 
                phot.vol, 
                mix.vol1 + rremove("ylab"), 
                het.vol + rremove("ylab"), 
                ncol = 3, nrow = 2, align = "hv"
                # , common.legend = T, legend.grob = get_legend(leg.v), legend = "right"
                , labels = "auto"
                )

p5.an <- annotate_figure(p5, 
                         bottom = textGrob("Experimental day", 
                                           gp = gpar(fontsize = 13)))

p5.h <- ggarrange(p5.an, leg.h, ncol = 1, align = "h", heights = c(3, 0.2))


ggsave("Plots/Dec2023/abund-biov_20231219-noleg.png", p5.an,  dpi=300, bg = "white") 

