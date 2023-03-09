
## set up ============
library(tidyverse)
library(readxl)

#needed for the plots 
trt.cols <- c(`C`= "#000000", #black - C
              `D`= "#0a8754", #g - D
              `I`= "#4472ca", #blu - I
              `E`= "#e84855") #r - E
mes.cols <- c(`1` = "#333333", `7` = "#000000", `16` = "#b3b3b3", #black - C
              `9` = "#e84855", `4`= "#8b2b33", `6` = "#f19199", #r - E
              `11` = "#0a8754",`8` = "#075f3b", `13` = "#54ab87", #g - D
              `12` = "#4472ca", `3`= "#30508d", `14` = "#7c9cda") #blu - I



group_names <- as_labeller(c(`aHF` = "Heterotrophs",
                             `aMF` = "Mixotrophs",
                             `aPF` = "Phototrophs"))
group_ir <- as_labeller(c(`HFir` = "Heterotrophs",
                          `MFir` = "Mixotrophs"))

## DATA ======
## and some data wrangling

## ABUNDANCE
master.dat <- read_xlsx("SITES_microscope_data_sheet.xlsx", 
                        sheet = "0.8_samples", 
                        range = "A1:T217", 
                        col_names = T)

# Area of square = 0.1 \* 0.1 (mm)
# Area of the filter = pi \* (21/2)^2^ (filter diameter is 25 mm)
# Volume filtered = 10 mL
# The volume of 1 field of view is (in mL):
Vp = 10 * (0.1 * 0.1) / (pi * (21/2)^2)

# for some reason all number columns are characters (possibly because of the NAs in most rows). Let's convert that 
master.dat[,14:20] <- lapply(master.dat[,14:20], as.numeric)
master.dat$Mesocosm <- as.factor(master.dat$Mesocosm)
master.dat$Mes_ID <- as.factor(master.dat$Mes_ID)

# some basic cleaning - wrangling
partial_dat <- master.dat %>% 
  select(-c(Date,Date_analysed,Box,Comments,`Photos/ImageJ`)) %>% 
  mutate(
    #maybe I need to do more factor stuff with other columns
    Treatment = fct_relevel(Treatment, c("C","D","I","E")),
    Tot_cells = HF_total + PF + MF,
    # calculate ingestion rates - NOW IT'S COUNTS
    # Mir = FLBinMF_total / (MF * 0.5),
    # Hir = FLBinHF_total / (HFwFLB * 0.5),
    # calculate abundances
    Vol = Fields_counted * Vp,
    aHF = HF_total/Vol,
    aMF = MF/Vol, 
    aPF = PF/Vol) %>% 
  
  # filter for Incubation 2 and T30 
  
  filter(Incubation != 1 
         # & Time_point == "Tend"
  )


## SIZE

size.dat <- read_xlsx("SITES_microscope_data_sheet.xlsx", 
                      sheet = "Sizes", 
                      range = "A1:U559", # CHANGE accordingly
                      col_names = T) %>% 
  
  select(Label, GROUP, Length) %>% 
  separate(Label, into = c("inc", "Bag", "time", "filter", "image"), sep = "_", remove = T) %>%
  select(-c(time, filter, image)) %>% 
  mutate(letter = rep(c("h", "d"), 279)) %>% 
  pivot_wider(names_from = letter, values_from = Length) %>% 
  unnest(cols = c(h, d)) %>% 
  mutate(vol = (pi/6)*d^2*h, # prolate spheroid in um^3
         Ccont = 0.216*vol^0.939, # Menden-Deuer & Lessard 2000, all non-diatom protists
         Incubation = as.factor(substring(inc, 3, 3)),
         Treatment = factor(substring(Bag, 1, 1), levels = c("C","D","I","E")),
         Mesocosm = as.factor(substring(Bag, 1, 2)),
         Replicate = substring(Bag, 3, 3)) %>% 
  group_by(Incubation, Treatment, GROUP) %>% 
  summarise(mean.cell.vol = mean(vol),
            sd.cell.vol = sd(vol))


biom = partial_dat %>% 
  filter(Time_point == "Tend") %>% 
  select(Treatment, Mesocosm, Mes_ID, Replicate, aHF, aMF, aPF) %>% 
  pivot_longer(cols = c(aHF, aPF, aMF), names_to = "group", values_to = "abundance") %>% 
  mutate(GROUP = substring(group, 2, 3))  %>% 
  inner_join(., size.dat, by = c("Treatment", "GROUP")) %>% 
  mutate(biomass = abundance * (0.216*mean.cell.vol^0.939)) # Menden-Deuer Lessard 2000, pgC/mL


## BACTERIA

bact.dat <- read_xlsx("SITES_microscope_data_sheet.xlsx", 
                      sheet = "0.2_samples",  
                      range = "A1:N37", # change accordingly
                      col_names = T) %>% 
  
  
  mutate(factor = (5 * Area)/(pi * (21/2)^2),
         HB_abund = HB/(Fields_counted*factor),
         FLB_abund = FLB/(Fields_counted*factor),
         CY_abund = CY/(Fields_counted*factor),
         FLB_perc = FLB_abund/HB_abund,
         
         Treatment = fct_relevel(Treatment, c("C","D","I","E")))

bact.dat$Mesocosm <- as.factor(bact.dat$Mesocosm)
bact.dat$Mes_ID <- as.factor(bact.dat$Mes_ID)


## INGESTION

ingest <- read_xlsx("SITES_microscope_data_sheet.xlsx", 
                    sheet = "Ingestion",  
                    range = "A1:L2233", # change accordingly
                    col_names = T) 

working <- ingest %>% 
  separate(Sample, into = c("inc", "Bag", "time", "filter"), sep = "_", remove = T) %>% 
  select(-time, -filter) %>% 
  group_by(inc, Treatment, Bag, Time_point, Grazer) %>% 
  # distinct(FLB_presence, Num_FLB, .keep_all = TRUE) 
  # mutate(nomnom = paste0(FLB_presence, "", Num_FLB))
  
  
  
  
  
  summarise(cell_num = sum(FLB_presence),
            flb_ingest = sum(Num_FLB, na.rm = T)) %>% # pay attention here! what's happening with the NAs?
  ungroup() %>% 
  pivot_wider(names_from = c(Time_point, Grazer), 
              values_from = c(cell_num, flb_ingest)) %>% 
  # correct ingestion rates
  mutate(MF_corr = cell_num_Tend_MF - cell_num_Tstart_MF,
         HF_corr = cell_num_Tend_HF - cell_num_Tstart_HF,
         HFir = flb_ingest_Tend_HF/(HF_corr*0.5),
         MFir = flb_ingest_Tend_MF/(MF_corr*0.5),
         Treatment = fct_relevel(Treatment, c("C","D","I","E"))) %>% 
  select(-contains(c("cell_num_T", "Tstart", "flb_ingest"))) 

working$Mes_ID <- as.factor(working$Mes_ID)

#if I replace with NA it gives NA in the summary
working[sapply(working, is.infinite)] <- 0

grazing =  left_join(working, bact.dat, by = c("Treatment", "Mes_ID", "Replicate"), copy = T) %>% 
  select(-c(9:20), -FLB_abund, -CY_abund) %>% 
  mutate(mCR = MFir/HB_abund, # in mL
         hCR = HFir/HB_abund)


## BASIC PLOTS =======

#abundance
partial_dat %>% 
  pivot_longer(cols = c(aHF, aPF, aMF), names_to = "group", values_to = "abundance") %>% 
  ggplot(., aes(x = Treatment, y = abundance)) +
  geom_boxplot(alpha = 0.8) +
  geom_point(aes(x = Treatment, y = abundance, color = Mes_ID),
             position = "jitter", 
             alpha = 0.6,
             # shape = 1, stroke = 0.7
  ) +
  ylab("Abundance cells/mL")+
  facet_grid(cols = vars(group), 
             labeller = group_names)+
  scale_color_manual(values = mes.cols)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey96"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))





partial_dat %>% 
  filter(Time_point=="Tend") %>% 
  pivot_longer(cols = c(aHF, aPF, aMF), names_to = "group", values_to = "abundance") %>%
  group_by(Incubation, Treatment, group) %>% 
  summarise(mean = mean(abundance), 
            sd = sd(abundance))%>% 
  mutate(Inc = as.factor(Incubation)
          # ,bag = as.factor(paste(Mesocosm, Replicate, sep = ""))
            ) %>% 
  
  ggplot(., aes(x = Inc, y = mean, colour = Treatment)) +
  geom_line(aes(group = Treatment))+
  # geom_point(aes(x = Incubation, y = mean, color = Treatment),
             # alpha = 0.6,
             # shape = 1, stroke = 0.7
  # ) +
  geom_pointrange(aes(ymin=mean - sd,
                      ymax=mean + sd),
                  size = .6, stat = "identity", show.legend = F, alpha = 0.7)+

  # geom_errorbar(aes(x = Incubation,
  #                   ymax = mean + sd,
  #                   ymin = mean - sd,
  #                   colour = Treatment), width = 0.05)+

  ylab("Abundance cells/mL")+
  facet_grid(rows = vars(group), 
             labeller = group_names,
             scales = "free_y")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey96"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  scale_color_manual(values = trt.cols)



partial_dat %>% 
  filter(Time_point=="Tend") %>%
  pivot_longer(cols = c(aHF, aPF, aMF), names_to = "group", values_to = "abundance") %>% 
  ggplot(., aes(x = Treatment, y = abundance)) +
  geom_boxplot(alpha = 0.8) +
  geom_point(aes(x = Treatment, y = abundance, color = Mes_ID),
             position = "jitter", 
             alpha = 0.6,
             # shape = 1, stroke = 0.7
  ) +
  ylab("Abundance cells/mL")+
  facet_grid(group~Incubation, 
             # labeller = group_names
             scales = "free_y")+
  scale_color_manual(values = mes.cols)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey96"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))




# mixotroph abundance is from the master dataframe, not corrected

# biomass
biom %>% 
  mutate(logB = log10(biomass)) %>% 
  ggplot(., aes(x = Treatment, y = c)) +
  geom_boxplot(alpha = 0.8) +
  geom_point(aes(x = Treatment, y = logB, color = Mes_ID),
             position = "jitter", 
             alpha = 0.6,
             # shape = 1, stroke = 0.7
  ) +
  ylab("Log pg C /mL")+
  facet_grid(cols = vars(group), 
             labeller = group_names)+
  scale_color_manual(values = mes.cols)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey96"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))

biom %>% 
  mutate(logB = log10(biomass),
         Inc = as.factor(Incubation)) %>% 
  group_by(Inc, Treatment, group) %>% 
  summarise(mean = mean(logB), 
            sd = sd(logB))%>% 
  ggplot(., aes(x = Inc, y = mean, colour = Treatment)) +
  geom_line(aes(group = Treatment))+
  # geom_point(aes(x = Incubation, y = mean, color = Treatment),
  # alpha = 0.6,
  # shape = 1, stroke = 0.7
  # ) +
  geom_pointrange(aes(ymin=mean - sd,
                      ymax=mean + sd),
                  size = .6, stat = "identity", show.legend = F, alpha = 0.7)+
  
  # geom_errorbar(aes(x = Incubation,
  #                   ymax = mean + sd,
  #                   ymin = mean - sd,
  #                   colour = Treatment), width = 0.05)+
  
  ylab("Log pg C /mL")+
  facet_grid(rows = vars(group), 
             labeller = group_names,
             scales = "free_y")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey96"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))+
  scale_color_manual(values = trt.cols)





# bacterial abundance
bact.dat %>% 
  pivot_longer(cols = c(HB_abund, FLB_abund, CY_abund), names_to = "group", values_to = "abundance") %>% 
  ggplot(., aes(x = Treatment, y = abundance)) +
  geom_boxplot(alpha = 0.8) +
  geom_point(aes(x = Treatment, y = abundance, color = Mes_ID),
             position = "jitter", 
             alpha = 0.6,
             # shape = 1, stroke = 0.7
  ) +
  ylab("Abundance cells/mL")+
  facet_grid(cols = vars(group))+
  scale_color_manual(values = mes.cols)+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey96"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))

# ingestion rate
working %>% 
  pivot_longer(cols = c(MFir, HFir), names_to = "group", values_to = "IR") %>% 
  ggplot(., aes(x = Treatment, y = IR)) +
  facet_grid(cols = vars(group), 
             labeller = group_ir)+
  # geom_boxplot(alpha = 0.8) +
  geom_point(aes(color = Mes_ID),
             position = "jitter", 
             alpha = 0.5,
             # shape = 1, stroke = 0.7
  ) +
  scale_color_manual(values = mes.cols)+
  geom_point(aes(color = Treatment), stat = "summary", fun = "mean",size = 4)+
  geom_pointrange(data =sumIR,
                  aes(x = sumIR$Treatment,
                      y = sumIR$IR.mean,
                      ymin=sumIR$IR.mean-sumIR$IR.sd,
                      ymax=sumIR$IR.mean+sumIR$IR.sd,
                      fill = sumIR$Treatment),
                  size = .6, stat = "identity", show.legend = F)+
  # geom_pointrange(stat = "summary", fun.y = mean, 
  #                 fun.ymin = mean-sd, fun.yman = mean+sd)
  scale_fill_manual(values = trt.cols)+
  ylab(expression("Ingestion rate"~"bacteria" ~ "cells"^{-1} ~"h"^{-1}))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey96"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))

## STATISTICS ========

library(glmmTMB)
library(ggeffects)
library(performance)
library(DHARMa)

# https://r.qcbs.ca/workshop07/book-en/choose-an-error-distribution.html

stat = partial_dat %>% 
  filter(Time_point=="Tend") %>% 
  left_join(working, bact.dat, by = c("Incubation", "Treatment", "Mes_ID", "Replicate"), copy = T)
  

hm = glmmTMB(aHF ~ 1
             + Treatment 
             + Incubation
             + (1|Mesocosm),
             # offset = Vol,
             family = poisson(),
             data = partial_dat %>% 
               filter(Time_point=="Tend"));summary(hm)

overdisp_fun(hm)

check_model(hm) # why error????

h.pred <- ggpredict(hm, terms = c("Treatment"))

plot(h.pred)
# simulateResiduals(fittedModel = hm, plot = T)
