library(tidyverse)
library(brms)
library(tidybayes)
library(modelr) #?
library(emmeans)
library(ggpubr)
library(grid)

# data
load("all_data.RData")

mytr <- function(x){
  log10(x+1)
}

dat <- data %>% 
  relocate(c(M.Ir, H.Ir, M.Gr, H.Gr), .after = biovol_MF) %>% 
  # log(x + 1) transormed
  mutate(across(4:16, mytr),
         ExpDay = factor(ExpDay)) 

# take out day 5
dat1 <- dat %>% filter(ExpDay !="5") %>% droplevels() %>% 
  mutate(Treatment = as.factor(Treatment))

Treats <- read_xlsx("Data/SITES_experiment_planning.xlsx",
                    sheet = "Treatments",
                    range = "A1:D21", 
                    col_names = T) %>% 
  pivot_longer(cols = c("E", "I", "D"), names_to = "Treatment", values_to = "intensity")


# Same dataframe with the one that was used in the models, but here ExpDay is
# is numeric and not a factor
nutr = envir %>% 
  select(-c(TOC, DN, DOsat, DOconc)) %>% 
  relocate(c(ExpDay, Treatment), .before = Mes_ID)

# general theme
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

# vertical and horizontal legends
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


## Fig. 1 - background data ====

# treatment additions
treatm <- Treats %>%
  ggplot(., aes(x = Experimental_day, y = intensity, fill = Treatment))+
  geom_bar(stat = "identity", position=position_dodge(.85), width = .7,
           # alpha = .3
  ) +
  scale_fill_manual(values = trt.cols,
                    labels = c("Daily", "Intermediate", "Extreme")) +
  scale_x_continuous(limits = c(0, 22),
                     expand = c(0, 0),
                     breaks = c(0, 5, 10, 15, 20),
                     label = c(0, 5, 10, 15, 20)) +  
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 110)) +
  geom_vline(aes(xintercept = 5), linetype = "dashed", color = "black", alpha = 0.4)+
  geom_vline(aes(xintercept = 13), linetype = "dashed", color = "black", alpha = 0.4)+
  geom_vline(aes(xintercept = 21), linetype = "dashed", color = "black", alpha = 0.4)+
  labs(y = "Addition %",
       x = NULL) +
  matheme+
  theme(
    axis.title = element_text(size = 13),
    legend.position = "none",
    legend.title = element_blank())
ggsave("Plots/20231102_treatments.png", dpi = 300)


# models
# files are generated on the OTHER R SCRIPT
par.m <-  brm(file = "models/231025_par")
abs.m <-  brm(file = "models/231029_abs")
chla.m <- brm(file = "models/231204_chla")
tn.m <-  brm(file = "models/231120_TN")
tp.m <-  brm(file = "models/231120_TP")



par.pl <- nutr %>% 
  filter(!is.na(PAR),
         ExpDay != 0) %>% 
  group_by(Mes_ID,ExpDay, Treatment) %>% 
  add_epred_draws(par.m,
                  re_formula = NA
  ) %>% 
  ggplot(., aes(x = ExpDay,
                y = PAR,
                colour = Treatment
  )) +
  stat_lineribbon(aes(y = (.epred)),
                  .width = 0,
                  position = position_dodge(.5),
                  size = .8
  )+
  stat_pointinterval(aes(y = (.epred)),
                     .width = c(0.95),
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



abs.pl <- nutr %>% 
  filter(!is.na(a420)) %>% 
  group_by(Mes_ID,ExpDay, Treatment) %>% 
  add_epred_draws(abs.m,
                  re_formula = NA,
  ) %>% 
  ggplot(., aes(x = ExpDay,
                y = a420,
                colour = Treatment,
  )) +
  stat_lineribbon(aes(y = (.epred)),
                  .width = 0,
                  position = position_dodge(.5),
                  size = .8
  )+
  stat_pointinterval(aes(y = (.epred)),
                     .width = c(0.95),
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



tn.pl <- nutr %>% 
  filter(!is.na(TN)) %>% 
  group_by(Mes_ID,ExpDay, Treatment) %>% 
  add_epred_draws(tn.m,
                  re_formula = NA
  ) %>% 
  ggplot(., aes(x = ExpDay,
                y = TN,
                colour = Treatment
  )) +
  stat_lineribbon(aes(y = (.epred)),
                  .width = 0,
                  position = position_dodge(.5),
                  size = .8
  )+
  stat_pointinterval(aes(y = (.epred)),
                     .width = c(0.95),
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


tp.pl <- nutr %>% 
  filter(!is.na(TP)) %>% 
  group_by(Mes_ID,ExpDay, Treatment) %>% 
  add_epred_draws(tp.m,
                  re_formula = NA
  ) %>% 
  ggplot(., aes(x = ExpDay,
                y = TP,
                colour = Treatment
  )) +
  stat_lineribbon(aes(y = (.epred)),
                  .width = 0,
                  position = position_dodge(.5),
                  size = .8
  )+
  stat_pointinterval(aes(y = (.epred)),
                     .width = c(0.95),
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


chla.pl <- nutr %>% 
  filter(!is.na(Chla)) %>% 
  group_by(Mes_ID,ExpDay, Treatment) %>% 
  add_epred_draws(chla.m,
                  re_formula = NA
  ) %>% 
  ggplot(., aes(x = ExpDay,
                y = Chla,
                colour = Treatment
  )) +
  stat_lineribbon(aes(y = (.epred)),
                  .width = 0,
                  position = position_dodge(.5),
                  size = .8
  )+
  stat_pointinterval(aes(y = (.epred)),
                     .width = c(0.95),
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



# Fig. 2 - bacterial abundance ======


#models
bact.m = brm(file = "models/231210_bacy.m")
cy.m = brm(file = "models/231210_cy.m")

# Heterotrophic bacteria
bact <- dat %>% 
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

# cyanobacteria
cyan <- dat %>% 
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


## Fig. 3 - nanoflagellate abundance and bovolume ----

# models
p.m = brm(file = "models/231210_phot.m")
h.m = brm(file = "models/231210_het.m")
m.m1 = brm(file = "models/231210_mix1.m")
ps.m = brm(file = "models/231211_phot-biovol.m")
hs.m = brm(file = "models/231211_het-biovol.m")
ms.m1 = brm(file = "models/231211_mix1-biovol.m")

# abunndances

phot <- dat %>% 
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
                     position = position_dodge(.5)
                     linewidth = 2, 
                     show.legend = FALSE
  ) +
  geom_point(
    data = dat,
    aes(x = ExpDay,
        y = PF_abund,
        colour = Treatment
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
                     position = position_dodge(.5)
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

# biovolumes

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
  labs( y = "<span style='font-size: 13pt'>Biovolume </span>
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
  labs(
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
  labs(
    y = "<span style='font-size: 13pt'>Biovolume </span>
         <span style='font-size: 11pt'>log(x+1) &nbsp;  \u00b5m\u00b3 mL\u207b\u00b9</span>",
    x = NULL)+
  matheme +
  theme(axis.title.y = ggtext::element_markdown(),
        axis.title.x = element_text(size = 13),
        plot.title = element_text(size=13, face="italic")
  ) 

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

## Fig. 4 - ingestion and grazing rates =====

# models
mir.m1 = brm(file = "models/231025_mIR.m1")
hir.m1 = brm(file = "models/231025_hIR.m1")
mgr.m1 = brm(file = "models/231025_mGR.m1")
hgr.m1 = brm(file = "models/231025_hGR.m1")

mir <- dat1 %>% 
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
                      linewidth = 2
   )+
   geom_point(
     data = dat1,
     aes(x = ExpDay,
         y = M.Ir,
         colour = Treatment),
     inherit.aes = FALSE,
     position = position_jitterdodge(dodge.width = .5),
     alpha = .35) +
   scale_color_manual(values = trt.cols)+
   scale_fill_manual(values = trt.cols)+
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


mgr <- dat1 %>% 
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
                       linewidth = 2
    )+
    geom_point(
      data = dat1,
      aes(x = ExpDay,
          y = M.Gr,
          colour = Treatment),
      inherit.aes = FALSE,
      position = position_jitterdodge(dodge.width = .5),
      alpha = .35) +
    scale_color_manual(values = trt.cols)+
    scale_fill_manual(values = trt.cols)+
    labs(
      y ="<span style='font-size: 13pt'>Grazing rate</span>
         <span style='font-size: 11pt'>(bacteria cell\u207b\u00b9 h\u207b\u00b9)</span>",
      x = NULL)+
    matheme+
    theme(
      axis.title.y = ggtext::element_markdown(),
      axis.title.x = element_text(size = 13))


hir <- dat1 %>% 
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
                       linewidth = 2
    )+
    geom_point(
      data = dat1,
      aes(x = ExpDay,
          y = H.Ir,
          colour = Treatment),
      inherit.aes = FALSE,
      position = position_jitterdodge(dodge.width = .5),
      alpha = .35) +
    scale_color_manual(values = trt.cols)+
    scale_fill_manual(values = trt.cols)+
    labs(title = "Heterotroph",
         y = NULL,
         x = NULL)+
    matheme+
    theme(plot.title = element_text(size=13, face="italic"))


hgr <- dat1 %>% 
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
                       linewidth = 2
    )+
    geom_point(
      data = dat1,
      aes(x = ExpDay,
          y = H.Gr,
          colour = Treatment),
      inherit.aes = FALSE,
      position = position_jitterdodge(dodge.width = .5),
      alpha = .35) +
    scale_color_manual(values = trt.cols)+
    scale_fill_manual(values = trt.cols)+
    matheme+
    labs(y = NULL,
         x = NULL)

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


