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
summary(pm.pairs, point.est = mean, level = .9)

post_apm <- posterior_samples(ap.m) %>% select(-lp__) %>% round(digits = 3)




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


# check missing values
bm.m = brm(bf(MF ~ Treatment 
            * Incubation
            + (1|Mes_ID)),
         family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta=0.95),
         seed=543,
         backend = "cmdstanr", 
         data = bdat,
         file = "models/bm.m")
pp_check(bm.m, ndraws = 100)
summary(bm.m, prob = .9)
plot(conditional_effects(bm.m, categorical = F, prob = .9), ask = FALSE)
plot(bm.m)

# plot
bm.cef <- conditional_effects(bm.m, categorical = F, prob = .9)
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
         
# ingestion rates 

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

pp_check(hir, ndraws = 100)
summary(hir, prob = .9)
plot(conditional_effects(hir, categorical = F, prob = .9), ask = FALSE)
plot(hir)

# plot
hir.cef <- conditional_effects(hir.m, categorical = F, prob = .9)
plot(hir.cef, plot = F)[[3]] +
  geom_line(size =1.2)+
  scale_color_manual(values = trt.cols, aesthetics = c("colour", "fill")) +
  geom_point(data = IR1, aes(x = Incubation, y = HFir, colour = Treatment), 
             inherit.aes = FALSE, position = position_jitter(width=0.12),
             alpha = 0.5)+
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  labs(y = "Heterotroph ingestion rate FLB/cell Hr")
