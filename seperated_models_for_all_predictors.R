# remove 1st incubation

dat1 <- dat %>% filter(Incubation!=1)
bdat1 <- bdat %>% filter(Incubation!=1)
ir1 <- IR1 %>% filter(Incubation!=1)
gr1 <- GR %>% filter(Incubation!=1)

phot.ab = brm(bf(aPF ~ Treatment
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
           data = dat1,
           file = "models/phot.ab_inc23")
pp_check(phot.ab, ndraws = 100)
summary(phot.ab, prob = .9)
plot(phot.ab, ask = F)
plot(conditional_effects(phot.ab, categorical = F, prob = .9), ask = FALSE)


het.ab = brm(bf(aHF ~ Treatment 
              * Incubation
              + (1|Mesocosm)),
           family = lognormal(link = "identity", link_sigma = "log"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta=0.95),
           seed=543,
           backend = "cmdstanr", 
           data = dat1,
           file = "models/het.ab_inc23")
pp_check(het.ab, ndraws = 100)
summary(het.ab, prob = .9)
plot(het.ab, ask = F)
plot(conditional_effects(het.ab, categorical = F, prob = .9), ask = FALSE)


mix.ab = brm(bf(aMFc ~ Treatment 
              * Incubation
              + (1|Mesocosm)),
           family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta=0.95),
           seed=543,
           backend = "cmdstanr", 
           data = dat1,
           file = "models/mix.ab_inc23")
pp_check(mix.ab, ndraws = 100)
plot(mix.ab, ask = F)
summary(mix.ab, prob = .9)
plot(conditional_effects(mix.ab, categorical = F, prob = .9), ask = FALSE)


phot.bio = brm(bf(PF ~ Treatment 
              * Incubation
              + (1|Mes_ID)),
           family = lognormal(link = "identity", link_sigma = "log"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta=0.99),
           seed=543,
           backend = "cmdstanr", 
           data = bdat1, 
           file = "models/phot.bio_inc23")
pp_check(phot.bio, ndraws = 100)
plot(phot.bio, ask = F)
summary(phot.bio, prob = .9)
plot(conditional_effects(phot.bio, categorical = F, prob = .9), ask = FALSE)


het.bio = brm(bf(HF ~ Treatment 
              * Incubation
              + (1|Mes_ID)),
           family = lognormal(link = "identity", link_sigma = "log"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta=0.95),
           seed=543,
           backend = "cmdstanr", 
           data = bdat1,
           file = "models/het.bio_inc23")
pp_check(het.bio, ndraws = 100)
plot(het.bio, ask = F)
summary(het.bio, prob = .9)
plot(conditional_effects(het.bio, categorical = F, prob = .9), ask = FALSE)


mix.bio = brm(bf(MF ~ Treatment 
              * Incubation
              + (1|Mes_ID)),
           family = lognormal(link = "identity", link_sigma = "log"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta=0.99),
           seed=543,
           backend = "cmdstanr", 
           data = bdat1, 
           file = "models/mix.bio_inc23")
pp_check(mix.bio, ndraws = 100)
plot(mix.bio, ask = F)
summary(mix.bio, prob = .9)
plot(conditional_effects(mix.bio, categorical = F, prob = .9), ask = FALSE)


Mir= brm(bf(MFir ~ Treatment 
               * Incubation
               + (1|Mes_ID),
               hu ~ Treatment
               + Incubation
               + (1|Mes_ID)
            ),
            family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta=0.95),
            seed = 543,
            backend = "cmdstanr", 
            data = ir1, 
            file = "models/Mir_inc23",
         file_refit = "on_change"
         )
pp_check(Mir, ndraws = 100)
plot(Mir, ask = F)
summary(Mir, prob = .9)
plot(conditional_effects(Mir, categorical = F, prob = .9), ask = FALSE)

mir.emt= emtrends(Mir, "Treatment", var = "Incubation")
mir.pairs = pairs(mir.emt, level = .9)
summary(mir.pairs, point.est = mean, level = .9)
mir.em = emmeans (Mir, pairwise  ~ Treatment | Incubation,
                  at = list(Incubation = c(2, 3)))
summary(mir.em, point.est = mean, level = .9)



Hir = brm(bf(HFir ~ Treatment 
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
            data = ir1,
            file = "models/Hir_inc23",
          file_refit = "on_change"
          )
pp_check(Hir, ndraws = 100)
plot(Hir, ask = F)
summary(Hir, prob = .9)
plot(conditional_effects(Hir, categorical = F, prob = .9), ask = FALSE)

hir.emt= emtrends(Hir, "Treatment", var = "Incubation")
hir.pairs = pairs(hir.emt, level = .9)
summary(hir.pairs, point.est = mean, level = .9)
hir.em = emmeans (Hir, pairwise  ~ Treatment | Incubation,
                  at = list(Incubation = c(2, 3)))
summary(hir.em, point.est = mean, level = .9)


Mgr = brm(bf(Gm ~ Treatment 
               * Incubation
               + (1|Mes_ID)
             ,
               hu ~ Treatment
               + Incubation
               + (1|Mes_ID)
             ),
            family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
            chains = 4,
            iter = 2000,
            cores = 4,
            control = list(adapt_delta=0.95),
            seed = 543,
            backend = "cmdstanr", 
            data = gr1, 
            file = "models/Mgr_inc23"
          )
pp_check(Mgr, ndraws = 100)
plot(Mgr, ask = F)
summary(Mgr, prob = .9)
plot(conditional_effects(Mgr, categorical = F, prob = .9), ask = FALSE)

mgr.emt= emtrends(Mgr, "Treatment", var = "Incubation")
mgr.pairs = pairs(mgr.emt, level = .9)
summary(mgr.pairs, point.est = mean, level = .9)
mgr.em = emmeans (Mgr, pairwise  ~ Treatment | Incubation,
                  at = list(Incubation = c(2, 3)))
summary(mir.em, point.est = mean, level = .9)

Hgr= brm(bf(Gh ~ Treatment 
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
            data = gr1,
            file = "models/Hgr_inc23"
)
pp_check(Hgr, ndraws = 100)
plot(Hgr, ask = F)
summary(Hgr, prob = .9)
plot(conditional_effects(Hgr, categorical = F, prob = .9), ask = FALSE)

hgr.emt= emtrends(Hgr, "Treatment", var = "Incubation")
hgr.pairs = pairs(hgr.emt, level = .9)
summary(hgr.pairs, point.est = mean, level = .9)
hgr.em = emmeans (Hgr, pairwise  ~ Treatment | Incubation,
                  at = list(Incubation = c(2, 3)))
summary(hgr.em, point.est = mean, level = .9)



gr_lst = list(Mgr, Hgr)

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



grdatt <- gr1 %>% 
  select(Incubation, Treatment, Mes_ID, Gh, Gm) %>% 
  pivot_longer(cols = c(Gh, Gm), names_to = "group", values_to = "grrate") %>% 
  mutate(group = case_when(group == "Gh" ~ "Heterotroph",
                           group == "Gm" ~ "Mixotroph"))

gr.df %>% 
  ggplot(., aes(x = Incubation,
                y = grrate,
                colour = Treatment,
                fill = Treatment)) +
  geom_point(data = grdatt, 
             aes(x = Incubation, 
                 y = grrate, 
                 colour = Treatment), 
             # inherit.aes = FALSE, 
             position = position_jitter(width = .02),
             alpha = .5) +
  facet_grid(rows = vars(group),
             scales = "free_y") +
  # geom_line(aes(y = .epred, 
  #               group = paste(Treatment, .draw)), alpha = .2) +
  stat_lineribbon(aes(y = (.epred)),
                  .width = .9,
                  point_interval = mean_qi,
                  size = 1,
                  alpha = .35,
                  # fill_ramp(from = trt.cols)
  ) +
  scale_x_continuous(breaks = c( 2, 3), 
                     labels = c(13, 21),
                     limits = c(2, 3)) +
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
  labs(y = "Grazing rate Bacteria/hr",
       x = "Experimental day") 
ggsave("Plots/20230522_grazing_rate-1.png", dpi = 300)






## plot ====

abunddata = list(phot.ab, het.ab, mix.ab)

ab.m <- lapply(abunddata, function(a) dat %>%
                 group_by(Mesocosm, Treatment) %>% 
                 data_grid(Incubation = seq_range(Incubation, n = 101)) %>% 
                 add_epred_draws(a,
                                 re_formula = NA,
                                 # ndraws = 100
                 ))
ab.m <- map_dfr(ab.m, ~ as.data.frame(.x), .id = "id") %>% 
  mutate(group = case_when(id == 1 ~ "Phot",
                           id == 2 ~ "Het",
                           id == 3 ~ "Mix"))



datt <- dat1 %>% 
select(Incubation, Treatment, Mesocosm, Mes_ID, aPF, aHF, aMFc) %>% 
  pivot_longer(cols = c(aHF, aPF, aMFc), names_to = "group", values_to = "abundance") %>% 
  mutate(group = case_when(group == "aHF" ~ "Phot",
                           group == "aPF" ~ "Het",
                           group == "aMFc" ~ "Mix"))
  
ab.m %>% 
  ggplot(., aes(x = Incubation,
                y = abundance,
                colour = Treatment,
                fill = Treatment)) +
  geom_point(data = datt, 
             aes(x = Incubation, 
                 y = abundance, 
                 colour = Treatment), 
             # inherit.aes = FALSE, 
             position = position_jitter(width=0.02),
             alpha = 0.5)+
  facet_grid(rows = vars(group),
             scales = "free_y")+
  # geom_line(aes(y = .epred, 
  #               group = paste(Treatment, .draw)), alpha = .2) +
  stat_lineribbon(aes(y = (.epred)),
                  .width = .9,
                  point_interval = "mean_hdi",
                  alpha = .3) +
  # geom_lineribbon(aes(y = (.epred)),
  #                 # .width = .8,
  #                 point_interval = "mean_hdi") +
  scale_x_continuous(breaks = c(2,3), limits = c(1.5, 3), labels = c(13, 21))+
  scale_color_manual(values = trt.cols, 
                     aesthetics = c("colour", "fill")) +
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  labs(y = "Abundance cells/mL",
       x = "Experimental day")




# one by one---------------

## phot abund =====
phot.ab1 = brm(bf(aPF ~ Treatment
                 + (1|Mesocosm)),
              family = lognormal(link = "identity", link_sigma = "log"),
              chains = 4,
              iter = 2000,
              cores = 4,
              control = list(adapt_delta=0.99),
              seed=543,
              backend = "cmdstanr", 
              data = dat %>% filter(Incubation == 1),
              file = "models/phot.ab1",
              # file_refit = "on_change"
              )
pp_check(phot.ab1, ndraws = 100)
plot(phot.ab1, ask = F)
summary(phot.ab1, prob = .9)
plot(conditional_effects(phot.ab1, categorical = F, prob = .9), ask = FALSE)

phot.ab2 = brm(bf(aPF ~ Treatment
                  + (1|Mesocosm)),
               family = lognormal(link = "identity", link_sigma = "log"),
               chains = 4,
               iter = 2000,
               cores = 4,
               control = list(adapt_delta=0.95),
               seed=543,
               backend = "cmdstanr", 
               data = dat %>% filter(Incubation == 2),
               file = "models/phot.ab2",
               # file_refit = "on_change"
               )
pp_check(phot.ab2, ndraws = 100)
plot(phot.ab2, ask = F)
summary(phot.ab2, prob = .9)
plot(conditional_effects(phot.ab2, categorical = F, prob = .9), ask = FALSE)

phot.ab3 = brm(bf(aPF ~ Treatment
                  + (1|Mesocosm)),
               family = lognormal(link = "identity", link_sigma = "log"),
               chains = 4,
               iter = 2000,
               cores = 4,
               control = list(adapt_delta=0.95),
               seed=543,
               backend = "cmdstanr", 
               data = dat %>% filter(Incubation == 3),
               file = "models/phot.ab3"
)
pp_check(phot.ab3, ndraws = 100)
plot(phot.ab3, ask = F)
summary(phot.ab3, prob = .9)
plot(conditional_effects(phot.ab3, categorical = F, prob = .9), ask = FALSE)

## het abund
het.ab1 = brm(bf(aHF ~ Treatment
                  + (1|Mesocosm)),
               family = lognormal(link = "identity", link_sigma = "log"),
               chains = 4,
               iter = 2000,
               cores = 4,
               control = list(adapt_delta=0.95),
               seed=543,
               backend = "cmdstanr", 
               data = dat %>% filter(Incubation == 1),
               file = "models/het.ab1")
pp_check(het.ab1, ndraws = 100)
plot(het.ab1, ask = F)
summary(het.ab1, prob = .9)
plot(conditional_effects(het.ab1, categorical = F, prob = .9), ask = FALSE)

het.ab2 = brm(bf(aHF ~ Treatment
                 + (1|Mesocosm)),
              family = lognormal(link = "identity", link_sigma = "log"),
              chains = 4,
              iter = 2000,
              cores = 4,
              control = list(adapt_delta=0.95),
              seed=543,
              backend = "cmdstanr", 
              data = dat %>% filter(Incubation == 2),
              file = "models/het.ab2")
pp_check(het.ab2, ndraws = 100)
plot(het.ab2, ask = F)
summary(het.ab2, prob = .9)
plot(conditional_effects(het.ab2, categorical = F, prob = .9), ask = FALSE)

het.ab3 = brm(bf(aHF ~ Treatment
                 + (1|Mesocosm)),
              family = lognormal(link = "identity", link_sigma = "log"),
              chains = 4,
              iter = 2000,
              cores = 4,
              control = list(adapt_delta=0.98),
              seed=543,
              backend = "cmdstanr", 
              data = dat %>% filter(Incubation == 3),
              file = "models/het.ab3",
              # file_refit = "on_change"
              )
pp_check(het.ab3, ndraws = 100)
plot(het.ab3, ask = F)
summary(het.ab3, prob = .9)
plot(conditional_effects(het.ab3, categorical = F, prob = .9), ask = FALSE)

## mixotr abund
mix.ab1 = brm(bf(aMFc ~ Treatment
                 + (1|Mesocosm)),
              family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
              chains = 4,
              iter = 2000,
              cores = 4,
              control = list(adapt_delta=0.99),
              seed=543,
              backend = "cmdstanr", 
              data = dat %>% filter(Incubation == 1),
              file = "models/mix.ab1",
              # file_refit = "on_change"
              )
pp_check(mix.ab1, ndraws = 100)
plot(mix.ab1, ask = F)
summary(mix.ab1, prob = .9)
plot(conditional_effects(mix.ab1, categorical = F, prob = .9), ask = FALSE)

mix.ab2 = brm(bf(aMFc ~ Treatment
                 + (1|Mesocosm)),
              family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
              chains = 4,
              iter = 2000,
              cores = 4,
              control = list(adapt_delta=0.95),
              seed=543,
              backend = "cmdstanr", 
              data = dat %>% filter(Incubation == 2),
              file = "models/mix.ab2")
pp_check(mix.ab2, ndraws = 100)
plot(mix.ab2, ask = F)
summary(mix.ab2, prob = .9)
plot(conditional_effects(mix.ab2, categorical = F, prob = .9), ask = FALSE)

mix.ab3 = brm(bf(aMFc ~ Treatment
                 + (1|Mesocosm)),
              family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
              chains = 4,
              iter = 2000,
              cores = 4,
              control = list(adapt_delta=0.95),
              seed=543,
              backend = "cmdstanr", 
              data = dat %>% filter(Incubation == 3),
              file = "models/mix.ab3")
pp_check(mix.ab3, ndraws = 100)
plot(mix.ab3, ask = F)
summary(mix.ab3, prob = .9)
plot(conditional_effects(mix.ab3, categorical = F, prob = .9), ask = FALSE)

## abundance plot efforts ====

abunddata = list(phot.ab1, phot.ab2,phot.ab3,
                 het.ab1, het.ab2, het.ab3,
                 mix.ab1, mix.ab2, mix.ab3)

all.ab.m = lapply(abunddata, function(a) dat %>%
                    group_by(Mesocosm, Incubation) %>% 
                    data_grid(Treatment) %>%
                    add_epred_draws(a, dpar = TRUE, category = "tobgp"))
# all.ab.m <- map_dfr(all.ab.m, ~ as.data.frame(.x), .id = "id")



for (i in 1:length(all.ab.m)) {
  
  newd <- list()
  
  
<- all.ab.m[i] %>%
  as.data.frame() %>% 
  filter(Incubation == case_when(i <= 3 ~ i,
                                i > 3 & i <=6 ~ i-3,
                                i > 6 ~ i-6)) 
  
}
  

nd <- lapply(all.ab.m, function(q){
  as.data.frame() %>% 
    filter(Incubation == case_when(i <= 3 ~ i,
                                   i > 3 & i <=6 ~ i-3,
                                   i > 6 ~ i-6)) 
})
  
  
  
  ggplot(., aes(x = Incubation, y = .epred, color = Treatment)) +
  stat_pointinterval(position = position_dodge(width = .2)) +
  scale_size_continuous(guide = "none") +
  scale_color_manual(values = trt.cols) +
  theme(panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey98"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  labs(y = "Heterotroph abundance cells/mL")
  
}


if (x<=3){
  x
} else if (x<=6) {
  x-3
} else {
  x-6
}
  



## phototr biom
phot.bio1 = brm(bf(PF ~ Treatment
                  + (1|Mes_ID)),
               family = lognormal(link = "identity", link_sigma = "log"),
               chains = 4,
               iter = 2000,
               cores = 4,
               control = list(adapt_delta=0.99),
               seed=543,
               backend = "cmdstanr", 
               data = bdat %>% filter(Incubation == 1), 
               file = "models/phot.bio1")
pp_check(phot.bio1, ndraws = 100)
plot(phot.bio1, ask = F)
summary(phot.bio1, prob = .9)
plot(conditional_effects(phot.bio1, categorical = F, prob = .9), ask = FALSE)

phot.bio2 = brm(bf(PF ~ Treatment
                   + (1|Mes_ID)),
                family = lognormal(link = "identity", link_sigma = "log"),
                chains = 4,
                iter = 2000,
                cores = 4,
                control = list(adapt_delta=0.99),
                seed=543,
                backend = "cmdstanr", 
                data = bdat %>% filter(Incubation == 2), 
                file = "models/phot.bio2")
pp_check(phot.bio2, ndraws = 100)
plot(phot.bio2, ask = F)
summary(phot.bio2, prob = .9)
plot(conditional_effects(phot.bio2, categorical = F, prob = .9), ask = FALSE)

phot.bio3 = brm(bf(PF ~ Treatment
                   + (1|Mes_ID)),
                family = lognormal(link = "identity", link_sigma = "log"),
                chains = 4,
                iter = 2000,
                cores = 4,
                control = list(adapt_delta=0.95),
                seed=543,
                backend = "cmdstanr", 
                data = bdat %>% filter(Incubation == 3), 
                file = "models/phot.bio3")
pp_check(phot.bio3, ndraws = 100)
plot(phot.bio3, ask = F)
summary(phot.bio3, prob = .9)
plot(conditional_effects(phot.bio3, categorical = F, prob = .9), ask = FALSE)


## het biom
het.bio1 = brm(bf(HF ~ Treatment
                 + (1|Mes_ID)),
              family = lognormal(link = "identity", link_sigma = "log"),
              chains = 4,
              iter = 2000,
              cores = 4,
              control = list(adapt_delta=0.95),
              seed=543,
              backend = "cmdstanr", 
              data = bdat%>% filter(Incubation == 1),
              file = "models/het.bio1")
pp_check(het.bio1, ndraws = 100)
plot(het.bio1, ask = F)
summary(het.bio1, prob = .9)
plot(conditional_effects(het.bio1, categorical = F, prob = .9), ask = FALSE)

het.bio2 = brm(bf(HF ~ Treatment
                  + (1|Mes_ID)),
               family = lognormal(link = "identity", link_sigma = "log"),
               chains = 4,
               iter = 2000,
               cores = 4,
               control = list(adapt_delta=0.95),
               seed=543,
               backend = "cmdstanr", 
               data = bdat %>% filter(Incubation == 2),
               file = "models/het.bio2")
pp_check(het.bio2, ndraws = 100)
plot(het.bio2, ask = F)
summary(het.bio2, prob = .9)
plot(conditional_effects(het.bio2, categorical = F, prob = .9), ask = FALSE)

het.bio3 = brm(bf(HF ~ Treatment
                  + (1|Mes_ID)),
               family = lognormal(link = "identity", link_sigma = "log"),
               chains = 4,
               iter = 2000,
               cores = 4,
               control = list(adapt_delta=0.95),
               seed=543,
               backend = "cmdstanr", 
               data = bdat %>% filter(Incubation == 3),
               file = "models/het.bio3")
pp_check(het.bio3, ndraws = 100)
plot(het.bio3, ask = F)
summary(het.bio3, prob = .9)
plot(conditional_effects(het.bio3, categorical = F, prob = .9), ask = FALSE)

## mix biomass
mix.bio1 = brm(bf(MF ~ Treatment
                 + (1|Mes_ID)),
              family = lognormal(link = "identity", link_sigma = "log"),
              chains = 4,
              iter = 2000,
              cores = 4,
              control = list(adapt_delta=0.99),
              seed=543,
              backend = "cmdstanr", 
              data = bdat %>% filter(Incubation == 1), 
              file = "models/mix.bio1")
pp_check(mix.bio1, ndraws = 100)
plot(mix.bio1, ask = F)
summary(mix.bio1, prob = .9)
plot(conditional_effects(mix.bio1, categorical = F, prob = .9), ask = FALSE)

mix.bio2 = brm(bf(MF ~ Treatment
                  + (1|Mes_ID)),
               family = lognormal(link = "identity", link_sigma = "log"),
               chains = 4,
               iter = 2000,
               cores = 4,
               control = list(adapt_delta=0.95),
               seed=543,
               backend = "cmdstanr", 
               data = bdat %>% filter(Incubation == 2), 
               file = "models/mix.bio2")
pp_check(mix.bio2, ndraws = 100)
plot(mix.bio2, ask = F)
summary(mix.bio2, prob = .9)
plot(conditional_effects(mix.bio1, categorical = F, prob = .9), ask = FALSE)

mix.bio3 = brm(bf(MF ~ Treatment
                  + (1|Mes_ID)),
               family = lognormal(link = "identity", link_sigma = "log"),
               chains = 4,
               iter = 2000,
               cores = 4,
               control = list(adapt_delta=0.95),
               seed=543,
               backend = "cmdstanr", 
               data = bdat %>% filter(Incubation == 3), 
               file = "models/mix.bio3")
pp_check(mix.bio3, ndraws = 100)
plot(mix.bio3, ask = F)
summary(mix.bio3, prob = .9)
plot(conditional_effects(mix.bio3, categorical = F, prob = .9), ask = FALSE)

## mix ing rate 
Mir1 = brm(bf(MFir ~ Treatment
              + (1|Mes_ID),
              hu ~ Treatment
              + (1|Mesocosm)),
           family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta=0.95),
           seed = 543,
           backend = "cmdstanr", 
           data = IR1 %>% filter(Incubation == 1), 
           file = "models/Mir1"
)
pp_check(Mir1, ndraws = 100)
plot(Mir1, ask = F)
summary(Mir1, prob = .9)
plot(conditional_effects(Mir1, categorical = F, prob = .9), ask = FALSE)

Mir2 = brm(bf(MFir ~ Treatment
              + (1|Mes_ID),
              hu ~ Treatment
              + (1|Mesocosm)),
           family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta=0.99),
           seed = 543,
           backend = "cmdstanr", 
           data = IR1 %>% filter(Incubation == 2), 
           file = "models/Mir2",
           # file_refit = "on_change"
)
pp_check(Mir2, ndraws = 100)
plot(Mir2, ask = F)
summary(Mir2, prob = .9)
plot(conditional_effects(Mir2, categorical = F, prob = .9), ask = FALSE)

Mir3 = brm(bf(MFir ~ Treatment
              + (1|Mes_ID),
              hu ~ Treatment
              + (1|Mesocosm)),
           family = hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta=0.95),
           seed = 543,
           backend = "cmdstanr", 
           data = IR1 %>% filter(Incubation == 3), 
           file = "models/Mir3"
)
pp_check(Mir3, ndraws = 100)
plot(Mir3, ask = F)
summary(Mir3, prob = .9)
plot(conditional_effects(Mir3, categorical = F, prob = .9), ask = FALSE)

## het ing rate
Hir1 = brm(bf(HFir ~ Treatment
             + (1|Mesocosm),
             hu ~ Treatment
             + (1|Mesocosm)),
          family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
          chains = 4,
          iter = 2000,
          cores = 4,
          control = list(adapt_delta=0.95),
          seed=543,
          backend = "cmdstanr", 
          data =  IR1 %>% filter(Incubation == 1),
          file = "models/Hir1")
pp_check(Hir1, ndraws = 100)
plot(Hir1, ask = F)
summary(Hir1, prob = .9)
plot(conditional_effects(Hir1, categorical = F, prob = .9), ask = FALSE)

Hir2 = brm(bf(HFir ~ Treatment
              + (1|Mesocosm),
              hu ~ Treatment
              + (1|Mesocosm)),
           family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta=0.95),
           seed=543,
           backend = "cmdstanr", 
           data =  IR1 %>% filter(Incubation == 2),
           file = "models/Hir2")
pp_check(Hir2, ndraws = 100)
plot(Hir2, ask = F)
summary(Hir2, prob = .9)
plot(conditional_effects(Hir2, categorical = F, prob = .9), ask = FALSE)

Hir3 = brm(bf(HFir ~ Treatment
              + (1|Mesocosm),
              hu ~ Treatment
              + (1|Mesocosm)),
           family = hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit"),
           chains = 4,
           iter = 2000,
           cores = 4,
           control = list(adapt_delta=0.95),
           seed=543,
           backend = "cmdstanr", 
           data =  IR1 %>% filter(Incubation == 3),
           file = "models/Hir3")
pp_check(Hir3, ndraws = 100)
plot(Hir3, ask = F)
summary(Hir3, prob = .9)
plot(conditional_effects(Hir3, categorical = F, prob = .9), ask = FALSE)
