library(brms)
library(tidybayes)
library(emmeans)
library(scales)
library(loo)
library(broom)
library(broom.mixed) 

# abundance

ap.m = brm(bf(aPF ~ Treatment
            + Incubation
            + Treatment:Incubation
            + (1|Mesocosm)),
         family = lognormal(link = "identity", link_sigma = "log"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta=0.95),
         seed=543,
         backend = "cmdstanr", 
         data = dat,
         file = "models/ap.m"
         )
pp_check(ap.m, ndraws = 100)
summary(ap.m, prob = .9)
plot(conditional_effects(ap.m, categorical = F, prob = .9), ask = FALSE)
plot(ap.m, ask = F)

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
         # file = "models/bm.m",
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
          file = "models/hir.m"
          )

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
            data = GR1, 
            # file = "models/mgr.m"
            )

pp_check(mgr.m, ndraws = 100)
summary(mgr.m, prob = .9)
plot(conditional_effects(mgr.m, categorical = F, prob = .9), ask = FALSE)
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
  labs(y = "MIxotroph ingestion rate FLB/cell Hr")

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
            data = GR1,
            # file = "models/hir.m"
)
pp_check(hgr.m, ndraws = 100)
summary(hgr.m, prob = .9)
plot(conditional_effects(hgr.m, categorical = F, prob = .9), ask = FALSE)
plot(hgr.m)
