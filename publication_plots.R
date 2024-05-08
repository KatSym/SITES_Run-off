library(tidyverse)
library(readxl)
library(brms)
library(tidybayes)
library(emmeans)
library(ggpubr)
library(grid)

# data
load("all_data.RData")

mytr <- function(x){
  log10(x+1)
}

dat <- data %>% 
  # log(x + 1) transormed
  mutate(across(4:11, mytr),
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
  # select(-c(TOC, DN, DOsat, DOconc)) %>% 
  relocate(c(ExpDay, Treatment), .before = Mes_ID)

# general theme
matheme <- theme(
  legend.position = "none",
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black", linewidth = .3),
  axis.text = element_text(size = 13, colour = "black"),
  axis.title = element_text(size = 15),
  # axis.text.y = element_text(size = 12),
  # strip.text.y = element_text(size = 12)
) 

# vertical and horizontal legends
leg.v <- ggplot(dat, aes(x = ExpDay, y = PF_abund, color = Treatment))+
  geom_point()+
  lims(y = c(0,0))+
  scale_color_manual(values = trt.cols,
                     labels = c("Control", "Daily", "Intermittent", "Extreme"))+
  theme_void()+
  theme(legend.position = c(0.5,0.5),
        legend.text = element_text(size =  13),
        legend.background	= element_blank(),
        legend.key	= element_blank(),
        legend.title = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=4))) # or 4

leg.h <- ggplot(dat, aes(x = ExpDay, y = PF_abund, color = Treatment))+
  geom_point()+
  lims(y = c(0,0))+
  scale_color_manual(values = trt.cols,
                     labels = c("Control", "Daily", "Intermittent", "Extreme"))+
  theme_void()+
  theme(legend.position = c(0.5,0.5),
        legend.text = element_text(size =  13,
          margin = margin(r = 30, unit = "pt")),
        legend.margin=margin(b = 0.6, unit='cm'),
        legend.direction = "horizontal",
        legend.background	= element_blank(),
        legend.key	= element_blank(),
        legend.title = element_blank())+
  guides(colour = guide_legend(override.aes = list(size=4)))# or 4

# make a two-line horizontal legend
leg.h2 <- leg.h +
  guides(colour = guide_legend(ncol=2, nrow=2, byrow=TRUE, override.aes = list(size=4)))



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
    axis.title = element_text(size = 15),
    legend.position = "none",
    legend.title = element_blank())
# ggsave("Plots/20231102_treatments.png", dpi = 300)


# models
# files are generated on the OTHER R SCRIPT
par.m <-  brm(file = "models/231025_par")
abs.m <-  brm(file = "models/231029_abs")
chla.m <- brm(file = "models/231204_chla")
tn.m <-  brm(file = "models/231120_TN")
tp.m <-  brm(file = "models/231120_TP")
no3.m <- brm(file = "models/240119_NO3")
no2.m <- brm(file = "models/240201_NO2")
nh.m <- brm(file = "models/240201_NH4")
po.m <- brm(file = "models/240201_PO4")


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
  # scale_x_continuous(expand = c(0,0))+
  matheme+
  theme(
    axis.title.y = ggtext::element_markdown(),
  ) +
  # geom_vline(aes(xintercept = 5), linetype = "dashed", color = "black", alpha = 0.4)+
  # geom_vline(aes(xintercept = 13), linetype = "dashed", color = "black", alpha = 0.4)+
  # geom_vline(aes(xintercept = 21), linetype = "dashed", color = "black", alpha = 0.4)+
  labs(y = "<span style='font-size: 15pt'>PAR</span>
         <span style='font-size: 12pt'>(μmol m\u207b\u00b2 s\u207b\u00b9)</span>",
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
  scale_x_continuous(expand = c(0,0))+
  matheme+
  theme(
    axis.title.y = ggtext::element_markdown(),
  ) +
  # geom_vline(aes(xintercept = 5), linetype = "dashed", color = "black", alpha = 0.4)+
  # geom_vline(aes(xintercept = 13), linetype = "dashed", color = "black", alpha = 0.4)+
  # geom_vline(aes(xintercept = 21), linetype = "dashed", color = "black", alpha = 0.4)+
  labs(y = "<span style='font-size: 15pt'>Abs<sub>420</sub></span>
         <span style='font-size: 12pt'>(m\u207b\u00b9)</span>",
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
  xlim(0, 21)+
  matheme+
  theme(
    axis.title.y = ggtext::element_markdown(),
  ) +
  # geom_vline(aes(xintercept = 5), linetype = "dashed", color = "black", alpha = 0.4)+
  # geom_vline(aes(xintercept = 13), linetype = "dashed", color = "black", alpha = 0.4)+
  # geom_vline(aes(xintercept = 21), linetype = "dashed", color = "black", alpha = 0.4)+
  labs(y = "<span style='font-size: 15pt'>Total&nbsp;nitrogen</span>
         <span style='font-size: 12pt'>(mg L\u207b\u00b9)</span>",
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
  xlim(0, 21)+
  matheme+
  theme(
    axis.title.y = ggtext::element_markdown(),
  ) +
  # geom_vline(aes(xintercept = 5), linetype = "dashed", color = "black", alpha = 0.4)+
  # geom_vline(aes(xintercept = 13), linetype = "dashed", color = "black", alpha = 0.4)+
  # geom_vline(aes(xintercept = 21), linetype = "dashed", color = "black", alpha = 0.4)+
  labs(y = "<span style='font-size: 15pt'>Total  phosphorus</span>
         <span style='font-size: 12pt'>(\u00b5g L\u207b\u00b9)</span>",
       x = NULL)

no3.pl <- nutr %>% 
  filter(!is.na(NO3)) %>% 
  group_by(Mes_ID,ExpDay, Treatment) %>% 
  add_epred_draws(no3.m,
                  re_formula = NA
  ) %>% 
  ggplot(., aes(x = ExpDay,
                y = NO3,
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
  scale_x_continuous(expand = c(0,0))+
  matheme+
  theme(
    axis.title.y = ggtext::element_markdown(),
  ) +
  geom_vline(aes(xintercept = 5), linetype = "dashed", color = "black", alpha = 0.4)+
  geom_vline(aes(xintercept = 13), linetype = "dashed", color = "black", alpha = 0.4)+
  geom_vline(aes(xintercept = 21), linetype = "dashed", color = "black", alpha = 0.4)+
  labs(y = "<span style='font-size: 15pt'>NO3</span>
         <span style='font-size: 12pt'>(\u00b5g L\u207b\u00b9)</span>",
       x = "Experimental day")


no2.pl <- nutr %>% 
  filter(!is.na(NO2)) %>% 
  group_by(Mes_ID,ExpDay, Treatment) %>% 
  add_epred_draws(no2.m,
                  re_formula = NA
  ) %>% 
  ggplot(., aes(x = ExpDay,
                y = NO2,
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
  scale_x_continuous(expand = c(0,0))+
  matheme+
  theme(
    axis.title.y = ggtext::element_markdown(),
  ) +
  
  labs(y = "<span style='font-size: 15pt'>NO<sup>-</sup><sub>2</sub></span>
         <span style='font-size: 12pt'>(\u00b5g L\u207b\u00b9)</span>",
       x = "Experimental day")

nh.pl <- nutr %>% 
  filter(!is.na(NH4)) %>% 
  group_by(Mes_ID,ExpDay, Treatment) %>% 
  add_epred_draws(nh.m,
                  re_formula = NA
  ) %>% 
  ggplot(., aes(x = ExpDay,
                y = NH4,
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
  scale_x_continuous(expand = c(0,0))+
  matheme+
  theme(
    axis.title.y = ggtext::element_markdown(),
  ) +
  geom_vline(aes(xintercept = 5), linetype = "dashed", color = "black", alpha = 0.4)+
  geom_vline(aes(xintercept = 13), linetype = "dashed", color = "black", alpha = 0.4)+
  geom_vline(aes(xintercept = 21), linetype = "dashed", color = "black", alpha = 0.4)+
  # NH<sup>+</sup><sub>4</sub></span>
  labs(y = "<span style='font-size: 15pt'>NH4</span>
         <span style='font-size: 12pt'>(\u00b5g L\u207b\u00b9)</span>",
       x = "Experimental day")



po.pl <- nutr %>% 
  filter(!is.na(PO4)) %>% 
  group_by(Mes_ID,ExpDay, Treatment) %>% 
  add_epred_draws(po.m,
                  re_formula = NA
  ) %>% 
  ggplot(., aes(x = ExpDay,
                y = PO4,
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
  scale_x_continuous(expand = c(0,0))+
  matheme+
  theme(
    axis.title.y = ggtext::element_markdown(),
  ) +
  geom_vline(aes(xintercept = 5), linetype = "dashed", color = "black", alpha = 0.4)+
  geom_vline(aes(xintercept = 13), linetype = "dashed", color = "black", alpha = 0.4)+
  geom_vline(aes(xintercept = 21), linetype = "dashed", color = "black", alpha = 0.4)+
  labs(y = "<span style='font-size: 15pt'>PO4</span>
         <span style='font-size: 12pt'>(\u00b5g L\u207b\u00b9)</span>",
       x = "Experimental day")



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
  labs(y = "<span style='font-size: 15pt'>Chl <i>a</i></span>
         <span style='font-size: 12pt'>(mg L\u207b\u00b9)</span>",
       x = NULL)


p.bckg <- ggarrange(treatm, par.pl + rremove("xlab"), abs.pl + rremove("xlab"),  
                    tn.pl, tp.pl, 
                    # chla.pl, 
                    ncol = 3, nrow = 2, align = "hv",
                    labels = "auto", font.label = list(size = 15),
                    common.legend = T, legend.grob = get_legend(leg.h), legend = "top"
) 

p.bckg <- annotate_figure(p.bckg, 
                          bottom = textGrob("Experimental day", gp = gpar(fontsize = 15)))

p.bckg

ggsave("Plots/Jan2024/background_1503-2.png", dpi = 300, bg = "white", 
       # width = 8, height = 5,
       scale = 1.18
       )


diss.nut <- ggarrange(no3.pl+ rremove("xlab"),
                      nh.pl+ rremove("xlab"), 
                      po.pl+ rremove("xlab"),  
                      ncol = 3, nrow = 1, align = "h",
                      labels = "auto", font.label = list(size = 15)
                      ,common.legend = T, legend.grob = get_legend(leg.h), legend = "top"
) 
diss.nut <- annotate_figure(diss.nut, 
                          bottom = textGrob("Experimental day", gp = gpar(fontsize = 15)))

ggsave("Plots/Jan2024/diss-nutr_1803-1.png", dpi = 300, bg = "white")

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
         y = "<span style='font-size: 12pt'>Abundance </span>
         <span style='font-size: 10pt'>log(x+1) cells mL\u207b\u00b9</span>",
         x = "Experimental Day") +
    matheme +
    theme(axis.title.y = ggtext::element_markdown(),
          axis.title.x = element_text(size = 12),
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
         y = "<span style='font-size: 12pt'>Abundance </span>
         <span style='font-size: 10pt'>log(x+1) cells mL\u207b\u00b9</span>",
         x = NULL) +
    matheme +
    theme( axis.title.y = ggtext::element_markdown(),
           axis.title.x = element_text(size = 12),
           plot.title = element_text(size=13, face="italic")
    ) 


p2.1 <- ggarrange(cyan,  bact,
                  ncol = 1, align = "v",
                  labels = "auto", font.label = list(size = 15),
                  common.legend = T, legend.grob = get_legend(leg.h2), legend = "top"
)

p2.2 <- ggarrange(cyan,  bact + rremove("xylab"),
                  nrow = 1, align = "h",
                  labels = "auto", font.label = list(size = 12),
                  common.legend = T, legend.grob = get_legend(leg.h), legend = "top"
)
p2.2an <- annotate_figure(p2.2, 
                          bottom = textGrob("Experimental day", gp = gpar(fontsize = 12)))

ggsave("Plots/Jan2024/bact_h.png", p2.2an, dpi=300, bg = "white")


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
                     position = position_dodge(.5),
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
       y = "<span style='font-size: 15pt'>Abundance </span>
         <span style='font-size: 13pt'>log(x+1) cells mL\u207b\u00b9</span>",
       x = NULL)+
  matheme +
  theme(axis.title.y = ggtext::element_markdown(),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(size=16, face="italic")
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
       y = "<span style='font-size: 15pt'>Abundance </span>
         <span style='font-size: 13pt'>Log(x+1) cells mL\u207b\u00b9</span>",
       x = NULL)+
  matheme +
  theme(axis.title.y = ggtext::element_markdown(),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(size=16, face="italic")
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
       y = "<span style='font-size: 15pt'>Abundance </span>
         <span style='font-size: 13pt'>Log(x+1) cells mL\u207b\u00b9</span>",
       x = NULL) +
  matheme +
  theme(axis.title.y = ggtext::element_markdown(),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(size=16, face="italic")
  ) 

# biovolumes

phot.vol <- dat %>% 
  group_by(Mes_ID,ExpDay, Treatment) %>% 
  add_epred_draws(ps.m,
                  re_formula = NA,
  ) %>% 
  ggplot(., aes(x = ExpDay,
                y = PF_biovol,
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
        y = PF_biovol,
        colour = Treatment),
    inherit.aes = FALSE,
    position = position_jitterdodge(dodge.width = .5),
    alpha = .4) +
  scale_color_manual(values = trt.cols)+
  scale_fill_manual(values = trt.cols)+
  labs( y = "<span style='font-size: 15pt'>Biovolume </span>
         <span style='font-size: 13pt'>log(x+1) &nbsp;  \u00b5m\u00b3 mL\u207b\u00b9</span>",
    x = NULL)+
  matheme +
  theme(axis.title.y = ggtext::element_markdown(),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(size=16, face="italic")
  ) 



het.vol <- dat %>% 
  group_by(Mes_ID,ExpDay, Treatment) %>% 
  add_epred_draws(hs.m,
                  re_formula = NA,
  ) %>% 
  ggplot(., aes(x = ExpDay,
                y = HF_biovol,
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
        y = HF_biovol,
        colour = Treatment),
    inherit.aes = FALSE,
    position = position_jitterdodge(dodge.width = .5),
    alpha = .4) +
  scale_color_manual(values = trt.cols)+
  scale_fill_manual(values = trt.cols)+
  labs(
    y = "<span style='font-size: 15pt'>Biovolume </span>
         <span style='font-size: 13pt'>log(x+1) &nbsp;  \u00b5m\u00b3 mL\u207b\u00b9</span>",
    x = NULL)+
  matheme +
  theme(axis.title.y = ggtext::element_markdown(),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(size=16, face="italic")
  ) 


mix.vol1 <- dat1 %>% 
  group_by(Mes_ID,ExpDay, Treatment) %>% 
  add_epred_draws(ms.m1,
                  re_formula = NA,
  ) %>% 
  ggplot(., aes(x = ExpDay,
                y = MF_biovol,
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
        y = MF_biovol,
        colour = Treatment),
    inherit.aes = FALSE,
    position = position_jitterdodge(dodge.width = .5),
    alpha = .4) +
  scale_color_manual(values = trt.cols)+
  scale_x_discrete(breaks = c("5", "13", "21"),
                   limits = c("5", "13", "21"))+
  labs(
    y = "<span style='font-size: 15pt'>Biovolume </span>
         <span style='font-size: 13pt'>log(x+1) &nbsp;  \u00b5m\u00b3 mL\u207b\u00b9</span>",
    x = NULL)+
  matheme +
  theme(axis.title.y = ggtext::element_markdown(),
        axis.title.x = element_text(size = 15),
        plot.title = element_text(size=16, face="italic")
  ) 

# combined abundance and biovolume
p5 <- ggarrange(phot, 
                mix1 + rremove("ylab"), 
                het + rremove("ylab"), 
                phot.vol, 
                mix.vol1 + rremove("ylab"), 
                het.vol + rremove("ylab"), 
                ncol = 3, nrow = 2, align = "hv"
                , labels = "auto", font.label = list(size = 15)
                , common.legend = T, legend.grob = get_legend(leg.h), legend = "top"
)

p5.an <- annotate_figure(p5, 
                         bottom = textGrob("Experimental day", 
                                           gp = gpar(fontsize = 15)))

# p5.h <- ggarrange(p5.an, leg.h, ncol = 1, align = "h", heights = c(3, 0.2))


ggsave("Plots/Jan2024/abund-biov_1503.png", p5.an,  dpi=300, bg = "white",
       scale = 1.18
       ) 

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
                 y = MF_Ir,
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
         y = MF_Ir,
         colour = Treatment),
     inherit.aes = FALSE,
     position = position_jitterdodge(dodge.width = .5),
     alpha = .35) +
   scale_color_manual(values = trt.cols)+
   scale_fill_manual(values = trt.cols)+
   labs(
     title = "Mixotroph",
     y ="<span style='font-size: 15pt'>Ingestion rate</span>
      <span style='font-size: 13pt'>(FLB cell\u207b\u00b9 h\u207b\u00b9)</span>",
     x = NULL)+
   matheme+
   theme(
     axis.title.y = ggtext::element_markdown(),
     axis.title.x = element_text(size = 15),
     plot.title = element_text(size=16, face="italic"))


mgr <- dat1 %>% 
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
    labs(
      y ="<span style='font-size: 15pt'>Grazing rate</span>
         <span style='font-size: 13pt'>(bacteria cell\u207b\u00b9 h\u207b\u00b9)</span>",
      x = NULL)+
    matheme+
    theme(
      axis.title.y = ggtext::element_markdown(),
      axis.title.x = element_text(size = 15))


hir <- dat1 %>% 
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
    labs(title = "Heterotroph",
         y = NULL,
         x = NULL)+
    matheme+
    theme(plot.title = element_text(size=16, face="italic"))


hgr <- dat1 %>% 
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
    matheme+
    labs(y = NULL,
         x = NULL)

p3 <- ggarrange(mir, hir, mgr, hgr, 
                ncol = 2, nrow = 2,
                align = "hv",
                labels = "auto", font.label = list(size = 15)
                , common.legend = T, legend.grob = get_legend(leg.h), legend = "top"
)
p3.an <- annotate_figure(p3, 
                         bottom = textGrob("Experimental day", 
                                           gp = gpar(fontsize = 13)))


ggsave("Plots/Jan2024/rates_20240112.png",p3.an, dpi = 300, bg = "white")


