library(brms)
library(tidybayes)
library(modelr)
library(emmeans)
# library(scales)
# library(loo)
# library(broom)
# library(broom.mixed) 

# abundance


ap.m = brm(bf(aPF ~ Treatment
            + Incubation
            + Treatment:Incubation
            + (1|Mesocosm)),
         family = lognormal(link = "identity", link_sigma = "log"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta = 0.95),
         seed = 543,
         backend = "cmdstanr", 
         data = dat,
         file = "models/ap.m"
         )
pp_check(ap.m, ndraws = 100)
summary(ap.m, prob = .9)
plot(conditional_effects(ap.m, categorical = F, prob = .9), ask = FALSE)
plot(ap.m, ask = F)

ap.m1 = brm(bf(aPF ~ Treatment
            + Incubation
            + Treatment:Incubation
            + (1|Mesocosm)
            + ar(p =1)),
         family = lognormal(link = "identity", link_sigma = "log"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta = 0.99),
         seed = 543,
         backend = "cmdstanr", 
         data = dat,
         # file = "models/ap.m1"
)
pp_check(ap.m1, ndraws = 100)
summary(ap.m1, prob = .9)
plot(conditional_effects(ap.m1, categorical = F, prob = .9), ask = FALSE)
plot(ap.m1, ask = F)


ap.m3 = brm(bf(aPF ~ Treatment
               + (1|Mesocosm)
               + ar(time = Incubation, p=1)),
            family = lognormal(link = "identity", link_sigma = "log"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta = 0.99),
            seed = 543,
            backend = "cmdstanr", 
            data = dat,
            # file = "models/ap.m1"
)





ap.m2 = brm(bf(aPF ~ Treatment
               + (1|Mesocosm)
               + (1|Incubation)),
            family = lognormal(link = "identity", link_sigma = "log"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta = 0.99),
            seed = 543,
            backend = "cmdstanr", 
            data = dat,
            # file = "models/ap.m1"
)
pp_check(ap.m2, ndraws = 100)
summary(ap.m2, prob = .9)
plot(conditional_effects(ap.m2, categorical = F, prob = .9), ask = FALSE)
plot(ap.m2, ask = F)


f1 <- loo(ap.m, is_method = "psis")
f2 <- loo(ap.m2, is_method = "psis")
loo_compare(f1, f2)

ap.cef <- conditional_effects(ap.m, categorical = F, prob = .9)
plot(ap.cef, plot = F)[[3]] +
  geom_line(size =1.2)+
  scale_color_manual(values = trt.cols, aesthetics = c("colour", "fill")) +
  geom_point(data = dat, aes(x = Incubation, y = aPF, colour = Treatment), 
             inherit.aes = FALSE, position = position_jitter(width=0.12),
             alpha = 0.5)+
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  labs(y = "Phototroph abundance cells/mL")

# slopes for each treatment
pm.emt= emtrends(ap.m, "Treatment", var = "Incubation")
summary(pm.emt, point.est = mean, level = .9)

# contrasts at intercept, or at specified points across the continuous predictor
pm.em = emmeans(ap.m, pairwise  ~ Treatment | Incubation)
pm.em = emmeans (ap.m, pairwise  ~ Treatment | Incubation,
                  at = list(Incubation = c(1, 2, 3)))
summary(pm.em, point.est = mean, level = .9)

# pairwise comparisons
pm.pairs = pairs(pm.emt, level = .9)
x <- summary(pm.pairs, point.est = mean, level = .9)

# post_apm <- posterior_samples(ap.m) %>% select(-lp__) %>% round(digits = 3)




ah.m = brm(bf(aHF ~ Treatment 
            * Incubation
            + (1|Mesocosm)),
         family = lognormal(link = "identity", link_sigma = "log"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta=0.95),
         seed=543,
         backend = "cmdstanr", 
         data = dat,
         file = "models/ah.m")
pp_check(ah.m, ndraws = 100)
summary(ah.m, prob = .9)
plot(conditional_effects(ah.m, categorical = F, prob = .9), ask = FALSE)
plot(ah.m)

ah.cef <- conditional_effects(ah.m, categorical = F, prob = .9)
plot(ah.cef, plot = F)[[3]] +
  geom_line(size =1.2)+
  scale_color_manual(values = trt.cols, aesthetics = c("colour", "fill")) +
  geom_point(data = dat, aes(x = Incubation, y = aHF, colour = Treatment), 
             inherit.aes = FALSE, position = position_jitter(width=0.12),
             alpha = 0.5)+
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  labs(y = "Heterotroph abundance cells/mL")

# slopes for each treatment
hm.emt= emtrends(ah.m, "Treatment", var = "Incubation")
summary(hm.emt, point.est = mean, level = .9)

# contrasts at intercept, or at specified points across the continuous predictor
hm.em = emmeans(hp.m, pairwise  ~ Treatment | Incubation)
hm.em = emmeans (hp.m, pairwise  ~ Treatment | Incubation,
                 at = list(Incubation = c(1, 2, 3)))
summary(hm.em, point.est = mean, level = .9)

# pairwise comparisons
hm.pairs = pairs(hm.emt, level = .9)
x <- summary(hm.pairs, point.est = mean, level = .9)



am.m = brm(bf(aMFc ~ Treatment 
            * Incubation
            + (1|Mesocosm)),
         family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta=0.95),
         seed=543,
         backend = "cmdstanr", 
         data = dat,
         file = "models/am.m")
pp_check(am.m, ndraws = 100)
summary(am.m, prob = .9)
plot(conditional_effects(am.m, categorical = F, prob = .9), ask = FALSE)
plot(am.m)


am.cef <- conditional_effects(am.m, categorical = F, prob = .9)
plot(am.cef, plot = F)[[3]] +
  geom_line(size =1.2)+
  scale_color_manual(values = trt.cols, aesthetics = c("colour", "fill")) +
  geom_point(data = dat, aes(x = Incubation, y = aMFc, colour = Treatment), 
             inherit.aes = FALSE, position = position_jitter(width=0.12),
             alpha = 0.5)+
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  labs(y = "Mixotroph abundance cells/mL")

# slopes for each treatment
mm.emt= emtrends(am.m, "Treatment", var = "Incubation")
summary(mm.emt, point.est = mean, level = .9)

# contrasts at intercept, or at specified points across the continuous predictor
mm.em = emmeans(am.m, pairwise  ~ Treatment | Incubation)
mm.em = emmeans (am.m, pairwise  ~ Treatment | Incubation,
                 at = list(Incubation = c(1, 2, 3)))
summary(mm.em, point.est = mean, level = .9)

# pairwise comparisons
mm.pairs = pairs(mm.emt, level = .9)
x <- summary(mm.pairs, point.est = mean, level = .9)




# biomass
bh.m = brm(bf(HF ~ Treatment 
            * Incubation
            + (1|Mes_ID)),
         family = lognormal(link = "identity", link_sigma = "log"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta=0.95),
         seed=543,
         backend = "cmdstanr", 
         data = bdat,
         file = "models/bh.m")
pp_check(bh.m, ndraws = 100) # not great
summary(bh.m, prob = .9)
plot(conditional_effects(bh.m, categorical = F, prob = .9), ask = FALSE)
plot(bh.m)

# plot
bh.cef <- conditional_effects(bh.m, categorical = F, prob = .9)
plot(bh.cef, plot = F)[[3]] +
  geom_line(size =1.2)+
  scale_color_manual(values = trt.cols, aesthetics = c("colour", "fill")) +
  geom_point(data = bdat, aes(x = Incubation, y = HF, colour = Treatment), 
             inherit.aes = FALSE, position = position_jitter(width=0.12),
             alpha = 0.5)+
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  labs(y = "Heterotroph biomass pg C/mL")
# slopes for each treatment
bh.emt= emtrends(bh.m, "Treatment", var = "Incubation")
summary(bh.emt, point.est = mean, level = .9)

# contrasts at intercept, or at specified points across the continuous predictor
bh.em = emmeans(bh.m, pairwise  ~ Treatment | Incubation)
bh.em = emmeans (bh.m, pairwise  ~ Treatment | Incubation,
                 at = list(Incubation = c(1, 2, 3)))
summary(bh.em, point.est = mean, level = .9)

# pairwise comparisons
bh.pairs = pairs(bh.emt, level = .9)
x <- summary(bh.pairs, point.est = mean, level = .9)





bp.m = brm(bf(PF ~ Treatment 
            * Incubation
            + (1|Mes_ID)),
         family = lognormal(link = "identity", link_sigma = "log"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta=0.99),
         seed=543,
         backend = "cmdstanr", 
         data = bdat, 
         file = "models/bp.m")
pp_check(bp.m, ndraws = 100)
summary(bp.m, prob = .9)
plot(conditional_effects(bp.m, categorical = F, prob = .9), ask = FALSE)
plot(bp.m)

# plot
bp.cef <- conditional_effects(bp.m, categorical = F, prob = .9)
plot(bp.cef, plot = F)[[3]] +
  geom_line(size =1.2)+
  scale_color_manual(values = trt.cols, aesthetics = c("colour", "fill")) +
  geom_point(data = bdat, aes(x = Incubation, y = PF, colour = Treatment), 
             inherit.aes = FALSE, position = position_jitter(width=0.12),
             alpha = 0.5)+
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  labs(y = "Photortoph biomass pg C/mL")

# slopes for each treatment
bp.emt= emtrends(bp.m, "Treatment", var = "Incubation")
summary(bp.emt, point.est = mean, level = .9)
# pairwise comparisons
bp.pairs = pairs(bp.emt, level = .9)
x <- summary(bp.pairs, point.est = mean, level = .9)


# check missing values
bm.m1 = brm(bf(MF ~ Treatment 
            * Incubation
            + (1|Mes_ID)),
        # family = lognormal(link = "identity", link_sigma = "log"),
         family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta=0.95),
         seed=543,
         backend = "cmdstanr", 
         data = bdat,
         file = "models/bm.m",
         # file_refit = "on_change"
         )
pp_check(bm.m1, ndraws = 100)
summary(bm.m1, prob = .9)
plot(conditional_effects(bm.m1, categorical = F, prob = .9), ask = FALSE)
plot(bm.m1)

# plot
bm.cef <- conditional_effects(bm.m1, categorical = F, prob = .9)
plot(bm.cef, plot = F)[[3]] +
  geom_line(size =1.2)+
  scale_color_manual(values = trt.cols, aesthetics = c("colour", "fill")) +
  geom_point(data = bdat, aes(x = Incubation, y = MF, colour = Treatment), 
             inherit.aes = FALSE, position = position_jitter(width=0.12),
             alpha = 0.5)+
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  labs(y = "Mixotroph biomass pg C/mL")


# slopes for each treatment
bm.emt= emtrends(bm.m1, "Treatment", var = "Incubation")
summary(bm.emt, point.est = mean, level = .9)
# pairwise comparisons
bm.pairs = pairs(bm.emt, level = .9)
x <- summary(bm.pairs, point.est = mean, level = .9)
         
# ingestion rates 

# if I remove incubation 1? -> NO
# ir23 <- IR1 %>% filter(Incubation!=1)

mir.m = brm(bf(MFir ~ Treatment 
             * Incubation
             + (1|Mes_ID),
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
          data = IR1, 
          file = "models/mir.m")

pp_check(mir.m, ndraws = 100)
summary(mir.m, prob = .9)
## we want Rhat to be 1, Bulk_ESS & Tail_ESS to make sence withthe total post-warmup draws
## something is significant when CI have the same sign (+-)
plot(conditional_effects(mir.m, categorical = F, prob = .9), ask = FALSE)
plot(mir.m)

# plot
mir.cef <- conditional_effects(mir.m, categorical = F, prob = .9)
plot(mir.cef, plot = F)[[3]] +
  geom_line(size =1.2)+
  scale_color_manual(values = trt.cols, aesthetics = c("colour", "fill")) +
  geom_point(data = IR1, aes(x = Incubation, y = MFir, colour = Treatment), 
             inherit.aes = FALSE, position = position_jitter(width=0.12),
             alpha = 0.5)+
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  labs(y = "MIxotroph ingestion rate FLB/cell Hr")

# slopes for each treatment
mir.emt= emtrends(mir.m, "Treatment", var = "Incubation")
summary(mir.emt, point.est = mean, level = .9)
# pairwise comparisons
mir.pairs = pairs(mir.emt, level = .9)
x <- summary(mir.pairs, point.est = mean, level = .9)





# n.ir1 <- IR1 %>% filter(HFir != 8)
hir.m = brm(bf(HFir ~ Treatment 
             * Incubation
             + (1|Mesocosm),
             hu ~ Treatment 
             + Incubation
             + (1|Mesocosm)),
          family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
          chains = 4,
          iter = 2000,
          cores = 4,
          control = list(adapt_delta=0.95),
          seed=543,
          backend = "cmdstanr", 
          data = IR1,
          file = "models/hir.m")

pp_check(hir.m, ndraws = 100)
summary(hir.m, prob = .9)
plot(conditional_effects(hir.m, categorical = F, prob = .9), ask = FALSE)
plot(hir.m)

# plot
hir.cef <- conditional_effects(hir.m, categorical = F, prob = .9)
plot(hir.cef, plot = F)[[3]] +
  geom_line(size =1.2)+
  scale_color_manual(values = trt.cols, aesthetics = c("colour", "fill")) +
  geom_point(data = n.ir1, aes(x = Incubation, y = HFir, colour = Treatment), 
             inherit.aes = FALSE, position = position_jitter(width=0.12),
             alpha = 0.5)+
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  labs(y = "Heterotroph ingestion rate FLB/cell Hr")

# slopes for each treatment
hir.emt= emtrends(hir.m, "Treatment", var = "Incubation")
summary(hir.emt, point.est = mean, level = .9)
# pairwise comparisons
hir.pairs = pairs(hir.emt, level = .9)
x <- summary(hir.pairs, point.est = mean, level = .9)



# grazing rates
mgr.m = brm(bf(Gm ~ Treatment 
               * Incubation
               + (1|Mes_ID),
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
            data = GR, 
            file = "models/mgr.m"
            )

pp_check(mgr.m, ndraws = 100)
summary(mgr.m, prob = .9)
plot(conditional_effects(mgr.m, categorical = F, prob = .9), ask = FALSE, main = "Mixotr grazing rate - flb/hb")
plot(mgr.m)

# plot
mgr.cef <- conditional_effects(mgr.m, categorical = F, prob = .9)
plot(mgr.cef, plot = F)[[3]] +
  geom_line(size =1.2)+
  scale_color_manual(values = trt.cols, aesthetics = c("colour", "fill")) +
  geom_point(data = GR, aes(x = Incubation, y = Gm, colour = Treatment), 
             inherit.aes = FALSE, position = position_jitter(width=0.12),
             alpha = 0.5)+
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  labs(y = "MIxotroph ingestion rate FLB/cell Hr")+
  ggtitle("Mixotr grazing rate - hb/flb")

# slopes for each treatment
mir.emt= emtrends(mir.m, "Treatment", var = "Incubation")
summary(mir.emt, point.est = mean, level = .9)
# pairwise comparisons
mir.pairs = pairs(mir.emt, level = .9)
x <- summary(mir.pairs, point.est = mean, level = .9)





hgr.m = brm(bf(Gh ~ Treatment 
               * Incubation
               + (1|Mesocosm),
               hu ~ Treatment 
               + Incubation
               + (1|Mesocosm)),
            family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta=0.95),
            seed=543,
            backend = "cmdstanr", 
            data = GR,
            file = "models/hir.m"
)
pp_check(hgr.m, ndraws = 100)
summary(hgr.m, prob = .9)
plot(conditional_effects(hgr.m, categorical = F, prob = .9), ask = FALSE)
plot(hgr.m)

##  plots =========
### abundance====
abund_lst = list(ap.m, ah.m, am.m)

ab <- lapply(abund_lst, function(a) dat %>%
                 group_by(Mesocosm, Treatment) %>% 
                 data_grid(Incubation = seq_range(Incubation, n = 101)) %>% 
                 add_epred_draws(a,
                                 re_formula = NA,
                                 # ndraws = 100
                 ))
ab.df <- map_dfr(ab, ~ as.data.frame(.x), .id = "id") %>% 
  mutate(group = case_when(id == 1 ~ "Phot",
                           id == 2 ~ "Het",
                           id == 3 ~ "Mix"))



datt <- dat %>% 
  select(Incubation, Treatment, Mesocosm, Mes_ID, aPF, aHF, aMFc) %>% 
  pivot_longer(cols = c(aHF, aPF, aMFc), names_to = "group", values_to = "abundance") %>% 
  mutate(group = case_when(group == "aHF" ~ "Het",
                           group == "aPF" ~ "Phot",
                           group == "aMFc" ~ "Mix"))


ab.df %>% 
  ggplot(., aes(x = Incubation,
                y = abundance,
                colour = Treatment,
                fill = Treatment)) +
  geom_point(data = datt, 
             aes(x = Incubation, 
                 y = abundance, 
                 colour = Treatment), 
             # inherit.aes = FALSE, 
             position = position_jitter(width = .02),
             alpha = .6) +
  facet_grid(rows = vars(group),
             scales = "free_y") +
  # geom_line(aes(y = .epred, 
  #               group = paste(Treatment, .draw)), alpha = .2) +
  stat_lineribbon(aes(y = (.epred)),
                  .width = .9,
                  point_interval = "mean_hdi",
                  size = 1,
                  alpha = .35,
                   # fill_ramp(from = trt.cols)
                  ) +
  scale_x_continuous(breaks = c(1, 2, 3), 
                     labels = c(5, 13, 21)) +
  scale_color_manual(values = trt.cols,
                     aesthetics = c("colour")) +
  # scale_color_manual(values =  c("#333333", "#0a8754", "#4472ca", "#e84855")) +
  scale_fill_manual(values = c("#b3b3b3", "#54ab87", "#7c9cda", "#f19199")) +
  # scale_fill_manual(values = c("#5e5e5e", "#54be86","#82a7ff",  "#ff7177")) +  
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11),
        # axis.text.y = element_text(size = 12),
        strip.text.y = element_text(size = 12)) +
  labs(y = "Abundance cells/mL",
       x = "Experimental day") 
ggsave("abund.png", dpi = 300)


### biomass=====
bio_lst = list(bp.m, bh.m, bm.m1)

bio <- lapply(bio_lst, function(a) bdat %>%
               group_by(Mes_ID, Treatment) %>% 
               data_grid(Incubation = seq_range(Incubation, n = 101)) %>% 
               add_epred_draws(a,
                               re_formula = NA,
                               # ndraws = 100
               ))
bio.df <- map_dfr(bio, ~ as.data.frame(.x), .id = "id") %>% 
  mutate(group = case_when(id == 1 ~ "Phototroph",
                           id == 2 ~ "Heterotroph",
                           id == 3 ~ "Mixotroph"))



bdatt <- bdat %>% 
  select(Incubation, Treatment, Mes_ID, PF, HF, MF) %>% 
  pivot_longer(cols = c(HF, PF, MF), names_to = "group", values_to = "biomass") %>% 
  mutate(group = case_when(group == "PF" ~ "Phototroph",
                           group == "HF" ~ "Heterotroph",
                           group == "MF" ~ "Mixotroph")) %>% 
  mutate(group = fct_relevel(group, c("Phototroph", "Heterotroph", "Mixotroph")))


bio.df %>% 
  ggplot(., aes(x = Incubation,
                y = biomass,
                colour = Treatment,
                fill = Treatment)) +
  geom_point(data = bdatt, 
             aes(x = Incubation, 
                 y = biomass, 
                 colour = Treatment), 
             # inherit.aes = FALSE, 
             position = position_jitter(width = .02),
             alpha = .6) +
  facet_grid(rows = vars(group),
             scales = "free_y") +
  # geom_line(aes(y = .epred, 
  #               group = paste(Treatment, .draw)), alpha = .2) +
  stat_lineribbon(aes(y = (.epred)),
                  .width = .9,
                  point_interval = "mean_hdi",
                  size = 1,
                  alpha = .35,
                  # fill_ramp(from = trt.cols)
  ) +
  scale_x_continuous(breaks = c(1, 2, 3), 
                     labels = c(5, 13, 21)) +
  # scale_color_manual(values = trt.cols) +
  scale_color_manual(values =  c("#000000", "#075f3b", "#30508d", "#8b2b33")) +
  # scale_color_manual(values =  c("#000000", "#05472c", "#243c6a", "#682026")) +
  scale_fill_manual(values = c("#b3b3b3", "#54ab87", "#7c9cda", "#f19199")) +
  # scale_fill_manual(values = c("#5e5e5e", "#54be86","#82a7ff",  "#ff7177")) +  
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11),
        # axis.text.y = element_text(size = 12),
        strip.text.y = element_text(size = 12)) +
  labs(y = "Biomass pg C/mL",
       x = "Experimental day") 
ggsave("abund.png", dpi = 300)
