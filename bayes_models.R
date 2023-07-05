library(brms)
library(tidybayes)
library(modelr)
library(emmeans)
# library(scales)
# library(loo)
# library(broom)
# library(broom.mixed) 


## Models ====

### abundance====


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

#### log model and plot ------
lap.m = brm(bf(log10(aPF) ~ Treatment
              + Incubation
              + Treatment:Incubation
              + (1|Mesocosm)),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta = 0.95),
           seed = 543,
           backend = "cmdstanr", 
           data = dat,
           file = "models/lap.m"
)
pp_check(lap.m, ndraws = 100)
summary(lap.m, prob = .9)
plot(conditional_effects(lap.m, categorical = F, prob = .9), ask = FALSE)
plot(lap.m, ask = F)


dat %>% 
  group_by(Mesocosm, Treatment) %>% 
  data_grid(Incubation = seq_range(Incubation, n = 101)) %>% 
  add_epred_draws(lap.m,
                  re_formula = NA,
                  # ndraws = 100
  ) %>% 
  ggplot(., aes(x = Incubation,
                y = log10(aPF),
                colour = Treatment,
                fill = Treatment)) +
  stat_lineribbon(aes(y = (.epred)),
                  .width = .9,
                  point_interval = "mean_hdi",
                  size = 1,
                  alpha = .7
  ) +
  geom_point(data = dat, 
             aes(x = Incubation, 
                 y = log10(aPF), 
                 colour = Treatment), 
             # inherit.aes = FALSE, 
             position = position_jitter(width = .02),
             alpha = .5) +
  scale_x_continuous(breaks = c(1, 2, 3), 
                     labels = c(5, 13, 21)) +
  scale_color_manual(values = trt.cols)+
  #                    aesthetics = c("colour")) +
  # scale_color_manual(values =  c("#000000", "#075f3b", "#30508d", "#8b2b33")) +
  # scale_fill_manual(values = c("#cccccc", "#cee7dd","#dae3f4",  "#fadadd")) +
  scale_fill_manual(values = c("#E5E5E5", "#D6E6E5","#EBF0F9",  "#FCECEE")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        # panel.background = element_rect(fill = "grey98"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = .3),
        axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11),
        # axis.text.y = element_text(size = 12),
        strip.text.y = element_text(size = 12)) +
  labs(y = "Log Abundance cells/mL",
       x = "Experimental day")+
ggsave("Plots/20230531_log-phot-abund.png", dpi = 300, width = 6.8, height = 2.2,units = "in")
#### tests----


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
# pm.em = emmeans(ap.m, pairwise  ~ Treatment | Incubation)
pm.em = emmeans (ap.m, pairwise  ~ Treatment | Incubation,
                  at = list(Incubation = c(1, 2, 3)))
summary(pm.em, point.est = mean, level = .9)

# pairwise comparisons
pm.pairs = pairs(pm.emt, level = .9)
summary(pm.pairs, point.est = mean, level = .9)

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
hm.em = emmeans(ah.m, pairwise  ~ Treatment | Incubation)
hm.em = emmeans (ah.m, pairwise  ~ Treatment | Incubation,
                 at = list(Incubation = c(1, 2, 3)))
summary(hm.em, point.est = mean, level = .9)

# pairwise comparisons
hm.pairs = pairs(hm.emt, level = .9)
summary(hm.pairs, point.est = mean, level = .9)



am.m = brm(bf(aMFc ~ Treatment 
            * Incubation
            + (1|Mesocosm)
            # + (1|Incubation)
            # + (1 + Incubation|Mesocosm)
            ),
         family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
         chains = 4,
         iter = 2000,
         cores = 4,
         control = list(adapt_delta=0.95),
         seed=543,
         backend = "cmdstanr", 
         data = dat,
         
         file = "models/am.m"
         )
pp_check(am.m, ndraws = 100)
summary(am.m, prob = .9)
plot(conditional_effects(am.m, categorical = F, prob = .9), ask = FALSE)
plot(am.m)

loo(am.m) # + (1 + Incubation|Mesocosm)
loo(am.m1)
loo_compare(loo(am.m), loo(am.m1))


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
summary(mm.pairs, point.est = mean, level = .9)

### biovol =====
bv = bvol %>% 
  # select(Incubation, Treatment, Mesocosm, GROUP, vol) %>% 
  pivot_wider(names_from = GROUP, values_from = vol) %>% 
  ungroup()

vol.p = brm(bf(PF ~ Treatment 
              * Incubation
              + (1|Mesocosm)),
           family = lognormal(link = "identity", link_sigma = "log"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta=0.95),
           seed=543,
           backend = "cmdstanr", 
           data = bv,
           file = "models/PF_biovol",
           file_refit = "on_change")
pp_check(vol.p, ndraws = 100)
summary(vol.p, prob = .9)
plot(conditional_effects(vol.p, categorical = F, prob = .9), ask = FALSE)
plot(vol.p)

# slopes for each treatment
vp.emt= emtrends(vol.p, "Treatment", var = "Incubation")
summary(vp.emt, point.est = mean, level = .9)
# pairwise comparisons
vp.pairs = pairs(vp.emt, level = .9)
summary(vp.pairs, point.est = mean, level = .9)
vp.em = emmeans (vol.p, pairwise  ~ Treatment | Incubation,
                  at = list(Incubation = c(1, 2, 3)))
summary(vp.em, point.est = mean, level = .9)


vol.h = brm(bf(HF ~ Treatment 
               * Incubation
               + (1|Mesocosm)),
            family = lognormal(link = "identity", link_sigma = "log"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta=0.95),
            seed=543,
            backend = "cmdstanr", 
            data = bv,
            file = "models/HF_biovol",
            file_refit = "on_change")
pp_check(vol.h, ndraws = 100)
summary(vol.h, prob = .9)
plot(conditional_effects(vol.h, categorical = F, prob = .9), ask = FALSE)
plot(vol.h)

# slopes for each treatment
vh.emt= emtrends(vol.h, "Treatment", var = "Incubation")
summary(vh.emt, point.est = mean, level = .9)
# pairwise comparisons
vh.pairs = pairs(vh.emt, level = .9)
summary(vh.pairs, point.est = mean, level = .9)
vh.em = emmeans (vol.h, pairwise  ~ Treatment | Incubation,
                 at = list(Incubation = c(1, 2, 3)))
summary(vh.em, point.est = mean, level = .9)



vol.m = brm(bf(MF ~ Treatment 
               * Incubation
               + (1|Mesocosm)),
            family = lognormal(link = "identity", link_sigma = "log"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta=0.95),
            seed=543,
            backend = "cmdstanr", 
            data = bv,
            file = "models/MF_biovol",
            file_refit = "on_change")
pp_check(vol.m, ndraws = 100)
summary(vol.m, prob = .9)
plot(conditional_effects(vol.m, categorical = F, prob = .9), ask = FALSE)
plot(vol.m)

# slopes for each treatment
vm.emt= emtrends(vol.m, "Treatment", var = "Incubation")
summary(vm.emt, point.est = mean, level = .9)
# pairwise comparisons
vm.pairs = pairs(vm.emt, level = .9)
summary(vm.pairs, point.est = mean, level = .9)
vm.em = emmeans (vol.m, pairwise  ~ Treatment | Incubation,
                 at = list(Incubation = c(1, 2, 3)))
summary(vm.em, point.est = mean, level = .9)


### biomass ====
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
         file = "models/bh.m",
         file_refit = "on_change")
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
summary(bh.pairs, point.est = mean, level = .9)





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
summary(bp.pairs, point.est = mean, level = .9)
bp.em = emmeans (bp.m, pairwise  ~ Treatment | Incubation,
                 at = list(Incubation = c(1, 2, 3)))
summary(bp.em, point.est = mean, level = .9)

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
summary(bm.pairs, point.est = mean, level = .9)
bm.em = emmeans (bm.m1, pairwise  ~ Treatment | Incubation,
                 at = list(Incubation = c(1, 2, 3)))
summary(bm.em, point.est = mean, level = .9)
         
### ingestion rates -----

# if I remove incubation 1? -> NO
# ir23 <- IR1 %>% filter(Incubation!=1)

mir.m = brm(bf(MFir ~ Treatment 
             * Incubation
             + (1|Mesocosm),
             hu ~ Treatment
             + Incubation
             + (1|Mesocosm)
             ),
          family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
          chains = 4,
          iter = 2000,
          cores = 4,
          control = list(adapt_delta=0.95),
          seed = 543,
          backend = "cmdstanr", 
          data = IR1
          # , 
          # file = "models/mir.m"
          )

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
summary(mir.pairs, point.est = mean, level = .9)
mir.em = emmeans (mir.m, pairwise  ~ Treatment | Incubation,
                 at = list(Incubation = c(1, 2, 3)))
summary(mir.em, point.est = mean, level = .9)





# n.ir1 <- IR1 %>% filter(HFir != 8)
hir.m = brm(bf(HFir ~ Treatment 
             * Incubation
             + (1|Mesocosm),
             hu ~ Treatment 
             + Incubation
             + (1|Mesocosm)),
            family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
          # family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
          chains = 4,
          iter = 2000,
          cores = 4,
          control = list(adapt_delta=0.95),
          seed=543,
          backend = "cmdstanr", 
          data = IR1
          # file = "models/hir.m",
          # file_refit = "on_change"
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
summary(hir.pairs, point.est = mean, level = .9)
hir.em = emmeans (hir.m, pairwise  ~ Treatment | Incubation,
                  at = list(Incubation = c(1, 2, 3)))
summary(hir.em, point.est = mean, level = .9)


### grazing rates -----
mgr.m = brm(bf(Gm ~ Treatment 
               * Incubation
               + (1|Mes_ID),
               hu ~ Treatment 
               + Incubation
               + (1|Mesocosm)),
            family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
            # family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
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
mgr.emt= emtrends(mgr.m, "Treatment", var = "Incubation")
summary(mgr.emt, point.est = mean, level = .9)
# pairwise comparisons
mgr.pairs = pairs(mgr.emt, level = .9)
summary(mgr.pairs, point.est = mean, level = .9)
mgr.em = emmeans (mgr.m, pairwise  ~ Treatment | Incubation,
                  at = list(Incubation = c(1, 2, 3)))
summary(mgr.em, point.est = mean, level = .9)


hgr.m = brm(bf(Gh ~ Treatment 
               * Incubation
               + (1|Mesocosm),
               hu ~ Treatment 
               + Incubation
               + (1|Mesocosm)),
            family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
            # family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta=0.95),
            seed=543,
            backend = "cmdstanr", 
            data = GR1,
            # file = "models/hir.m",
            # file_refit = "on_change"
)
pp_check(hgr.m, ndraws = 100)
summary(hgr.m, prob = .9)
plot(conditional_effects(hgr.m, categorical = F, prob = .9), ask = FALSE)
plot(hgr.m)

# plot
hgr.cef <- conditional_effects(hgr.m, categorical = F, prob = .9)
plot(hgr.cef, plot = F)[[3]] +
  geom_line(size =1.2)+
  scale_color_manual(values = trt.cols, aesthetics = c("colour", "fill")) +
  geom_point(data = GR, aes(x = Incubation, y = Gh, colour = Treatment), 
             inherit.aes = FALSE, position = position_jitter(width=0.12),
             alpha = 0.5)+
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  labs(y = "MIxotroph ingestion rate FLB/cell Hr")+
  ggtitle("Mixotr grazing rate - hb/flb")

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
  mutate(group = case_when(id == 1 ~ "Phototroph",
                           id == 2 ~ "Heterotroph",
                           id == 3 ~ "Mixotroph"))



datt <- dat %>% 
  select(Incubation, Treatment, Mesocosm, Mes_ID, aPF, aHF, aMFc) %>% 
  pivot_longer(cols = c(aHF, aPF, aMFc), names_to = "group", values_to = "abundance") %>% 
  mutate(group = case_when(group == "aHF" ~ "Heterotroph",
                           group == "aPF" ~ "Phototroph",
                           group == "aMFc" ~ "Mixotroph"))


ab.df %>% 
  ggplot(., aes(x = Incubation,
                y = abundance,
                colour = Treatment,
                fill = Treatment)) +
  facet_grid(rows = vars(group),
             scales = "free_y") +
  # geom_line(aes(y = .epred, 
  #               group = paste(Treatment, .draw)), alpha = .2) +
  stat_lineribbon(aes(y = (.epred)),
                  .width = .9,
                  point_interval = "mean_hdi",
                  size = 1,
                  alpha = .7,
                   # fill_ramp(from = trt.cols)
                  ) +
  geom_point(data = datt, 
             aes(x = Incubation, 
                 y = abundance, 
                 colour = Treatment), 
             # inherit.aes = FALSE, 
             position = position_jitter(width = .02),
             alpha = .5) +
  scale_y_continuous(breaks = waiver(),
                     # labels = `breaks`/1000
                     )+
  scale_x_continuous(breaks = c(1, 2, 3), 
                     labels = c(5, 13, 21)) +
  scale_color_manual(values = trt.cols)+
  #                    aesthetics = c("colour")) +
  # scale_color_manual(values =  c("#000000", "#075f3b", "#30508d", "#8b2b33")) +
  # scale_fill_manual(values = c("#cccccc", "#cee7dd","#dae3f4",  "#fadadd")) +
  scale_fill_manual(values = c("#E5E5E5", "#D6E6E5","#EBF0F9",  "#FCECEE")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        # panel.background = element_rect(fill = "grey98"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = .3),
        axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11),
        # axis.text.y = element_text(size = 12),
        strip.text.y = element_text(size = 12)) +
  labs(y = expression(paste("Abundance cells mL" ^{-1})),
       x = "Experimental day")
ggsave("Plots/20230531_abund.png", dpi = 300, width = 7.34, height = 5.93, units = "in")

### biovolume =====

bv_lst = list(vol.p, vol.h, vol.m)

bvl <- lapply(bv_lst, function(a) bv %>%
                group_by(Mesocosm, Treatment) %>% 
                data_grid(Incubation = seq_range(Incubation, n = 101)) %>% 
                add_epred_draws(a,
                                re_formula = NA,
                                # ndraws = 100
                ))
bv.df <- map_dfr(bvl, ~ as.data.frame(.x), .id = "id") %>% 
  mutate(group = case_when(id == 1 ~ "Phototroph",
                           id == 2 ~ "Heterotroph",
                           id == 3 ~ "Mixotroph")) 



bvdatt <- bv %>% 
  select(Incubation, Treatment, Mesocosm, PF, HF, MF) %>% 
  # filter(PF < 700) %>% 
  pivot_longer(cols = c(HF, PF, MF), names_to = "group", values_to = "biovol") %>% 
  mutate(group = case_when(group == "PF" ~ "Phototroph",
                           group == "HF" ~ "Heterotroph",
                           group == "MF" ~ "Mixotroph")) %>% 
  mutate(group = fct_relevel(group, c("Phototroph", "Heterotroph", "Mixotroph"))) %>% 
  drop_na()

bv.df %>% 
  ggplot(., aes(x = Incubation,
                y = biovol,
                colour = Treatment,
                fill = Treatment)) +
  facet_grid(rows = vars(group),
             scales = "free_y") +
  stat_lineribbon(aes(y = (.epred)),
                  .width = .9,
                  point_interval = "mean_hdi",
                  size = 1,
                  alpha = .7,
  ) +
  geom_point(data = bvdatt, 
             aes(x = Incubation, 
                 y = biovol, 
                 colour = Treatment), 
             position = position_jitter(width = .02),
             alpha = .5) +
  scale_x_continuous(breaks = c(1, 2, 3), 
                     labels = c(5, 13, 21)) +
  # change the scale so that PF difference is more visible
  # scale_y_continuous(limits = c(0, 500))+
  scale_color_manual(values = trt.cols)+
  scale_fill_manual(values = c("#E5E5E5", "#D6E6E5","#EBF0F9",  "#FCECEE")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = .3),
        axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11),
        strip.text.y = element_text(size = 12)) +
  labs(y = "Biovolume \u03bcm\u00b3",
       x = "Experimental day") 
ggsave("Plots/20230601_biovol.png", dpi = 300, width = 7.34, height = 5.93, units = "in")

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
  facet_grid(rows = vars(group),
             scales = "free_y") +
  # geom_line(aes(y = .epred, 
  #               group = paste(Treatment, .draw)), alpha = .2) +
  stat_lineribbon(aes(y = (.epred)),
                  .width = .9,
                  point_interval = "mean_hdi",
                  size = 1,
                  alpha = .7,
                  # fill_ramp(from = trt.cols)
  ) +
  geom_point(data = bdatt, 
             aes(x = Incubation, 
                 y = biomass, 
                 colour = Treatment), 
             # inherit.aes = FALSE, 
             position = position_jitter(width = .02),
             alpha = .5) +
  scale_x_continuous(breaks = c(1, 2, 3), 
                     labels = c(5, 13, 21)) +
  scale_y_continuous(labels = scales::label_number(scale = 1/1000))+
  # scale_color_manual(values =  c("#000000", "#075f3b", "#30508d", "#8b2b33")) +
  # scale_fill_manual(values = c("#cccccc", "#cee7dd","#dae3f4",  "#fadadd")) +
  scale_color_manual(values = trt.cols)+
  scale_fill_manual(values = c("#E5E5E5", "#D6E6E5","#EBF0F9",  "#FCECEE")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = .3),
        axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11),
        strip.text.y = element_text(size = 12)) +
  labs(y = expression(paste("Biomass ng C mL" ^{-1})),
       x = "Experimental day") +
ggsave("Plots/20230531_biomass.png", dpi = 300, width = 7.34, height = 5.93, units = "in")

### ingestion rates =====

ir_lst = list(mir.m, hir.m)

irl <- lapply(ir_lst, function(a) bdat %>%
                group_by(Mes_ID, Treatment) %>% 
                data_grid(Incubation = seq_range(Incubation, n = 101)) %>% 
                add_epred_draws(a,
                                re_formula = NA,
                                # ndraws = 100
                ))
ir.df <- map_dfr(irl, ~ as.data.frame(.x), .id = "id") %>% 
  mutate(group = case_when(id == 1 ~ "Mixotroph",
                           id == 2 ~ "Heterotroph"))



irdatt <- IR1 %>% 
  select(Incubation, Treatment, Mes_ID, HFir, MFir) %>% 
  pivot_longer(cols = c(HFir, MFir), names_to = "group", values_to = "ingrate") %>% 
  mutate(group = case_when(group == "HFir" ~ "Heterotroph",
                           group == "MFir" ~ "Mixotroph"))

ir.df %>% 
  ggplot(., aes(x = Incubation,
                y = ingrate,
                colour = Treatment,
                fill = Treatment)) +
  facet_grid(rows = vars(group),
             scales = "free_y") +
  stat_lineribbon(aes(y = (.epred)),
                  .width = .9,
                  point_interval = mean_qi,
                  size = 1,
                  alpha = .7,
  ) +
  geom_point(data = irdatt, 
             aes(x = Incubation, 
                 y = ingrate, 
                 colour = Treatment), 
             position = position_jitter(width = .02),
             alpha = .5) +
  scale_x_continuous(breaks = c(1, 2, 3), 
                     labels = c(5, 13, 21)) +
  # scale_color_manual(values =  c("#000000", "#075f3b", "#30508d", "#8b2b33")) +
  # scale_fill_manual(values = c("#cccccc", "#cee7dd","#dae3f4",  "#fadadd")) +
  scale_color_manual(values = trt.cols)+
  scale_fill_manual(values = c("#E5E5E5", "#D6E6E5","#EBF0F9",  "#FCECEE")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = .3),
        axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11),
        strip.text.y = element_text(size = 12)) +
  labs(y = expression(paste("Ingestion rate FLB cell" ^{-1}, "hr" ^{-1})),
       x = "Experimental day") 
ggsave("Plots/20230531_ingestion_rate.png", dpi = 300, width = 7.34, height = 5.93, units = "in")


### grazing rates =====

gr_lst = list(mgr.m, hgr.m)

grl <- lapply(gr_lst, function(a) bdat %>%
                group_by(Mes_ID, Treatment) %>% 
                data_grid(Incubation = seq_range(Incubation, n = 101)) %>% 
                add_epred_draws(a,
                                re_formula = NA,
                                # ndraws = 100
                ))
gr.df <- map_dfr(grl, ~ as.data.frame(.x), .id = "id") %>% 
  mutate(group = case_when(id == 1 ~ "Mixotroph",
                           id == 2 ~ "Heterotroph"))



grdatt <- GR1 %>% 
  select(Incubation, Treatment, Mes_ID, Gh, Gm) %>% 
  pivot_longer(cols = c(Gh, Gm), names_to = "group", values_to = "grrate") %>% 
  mutate(group = case_when(group == "Gh" ~ "Heterotroph",
                           group == "Gm" ~ "Mixotroph"))

gr.df %>% 
  ggplot(., aes(x = Incubation,
                y = grrate,
                colour = Treatment,
                fill = Treatment)) +
  facet_grid(rows = vars(group),
             scales = "free_y") +
  stat_lineribbon(aes(y = (.epred)),
                  .width = .9,
                  point_interval = mean_hdi,
                  size = 1,
                  alpha = .7,
  ) +
  geom_point(data = grdatt, 
             aes(x = Incubation, 
                 y = grrate, 
                 colour = Treatment), 
             position = position_jitter(width = .02),
             alpha = .5) +
  scale_x_continuous(breaks = c(1, 2, 3), 
                     labels = c(5, 13, 21)) +
  # scale_color_manual(values =  c("#000000", "#075f3b", "#30508d", "#8b2b33")) +
  # scale_fill_manual(values = c("#cccccc", "#cee7dd","#dae3f4",  "#fadadd")) +
  scale_color_manual(values = trt.cols)+
  scale_fill_manual(values = c("#E5E5E5", "#D6E6E5","#EBF0F9",  "#FCECEE")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = .3),
        axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11),
        strip.text.y = element_text(size = 12)) +
  labs(y = expression(paste("Grazing rate Bacteria hr" ^{-1})),
       x = "Experimental day") 
ggsave("Plots/20230531_grazing_rate-axis.png", dpi = 300, width = 7.34, height = 5.93, units = "in")
