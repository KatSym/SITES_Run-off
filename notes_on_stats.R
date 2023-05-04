
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
master.dat <- read_xlsx("SITES_microscope_data_working.xlsx", 
                        sheet = "0.8_samples", 
                        range = "A1:Q217", 
                        col_names = T)

# for some reason all number columns are characters (possibly because of the NAs in most rows). Let's convert that 
master.dat[,10:17] <- lapply(master.dat[,10:17], as.numeric)
master.dat$Mesocosm <- as.factor(master.dat$Mesocosm)
master.dat$Mes_ID <- as.factor(master.dat$Mes_ID)
# master.dat$Incubation <- as.factor(master.dat$Incubation)

# Area of square = 0.1 \* 0.1 (mm)
# Area of the filter = pi \* (21/2)^2^ (filter diameter is 25 mm)
# Volume filtered = 10 mL
# The volume of 1 field of view is (in mL):
Vp = 10 * (0.1 * 0.1) / (pi * (21/2)^2)

# some basic cleaning - wrangling
partial_dat <- master.dat %>% 
  select(-Date) %>% 
  mutate(
    Treatment = fct_relevel(Treatment, c("C","D","I","E")),
    # calculate abundances
    Vol = Fields_counted * Vp,
    aHF = HF/Vol,
    aMF = MF/Vol, 
    aPF = PF/Vol
    ) 
# replace NA values in the Cells_counted column for Inc 2 and Tstart with calculated values
rep.na <- master.dat %>% filter(Incubation == 2 & Time_point=="Tend")
celltofield = rep.na$Cells_counted/rep.na$Fields_counted

partial_dat <- within(partial_dat, Cells_counted <- ifelse(is.na(Cells_counted),
                                                           celltofield * Fields_counted,
                                                           Cells_counted))


## SIZE
        # THIS IS MESSY
# Load the size data. If the number of cells (rows) are more than 30 withing the Incubation - Treatment
# - GROUP grouping, select only 30. This is done in a new df. Then take the <30 rows group and merge them with 
# the new df.
size.dat <- read_xlsx("SITES_microscope_data_working.xlsx", 
                      sheet = "Sizes", 
                      range = "A1:C2043",
                      col_names = T) %>% 
  
  separate(Label, into = c("inc", "Bag", "time", "filter", "image"), sep = "_", remove = T) %>%
  select(-c(time, filter, image)) %>% 
  mutate(letter = rep(c("h", "d"), 1021)) %>% 
  pivot_wider(names_from = letter, values_from = Length) %>% 
  unnest(cols = c(h, d)) %>% 
  mutate(vol = (pi/6)*d^2*h, # prolate spheroid in um^3
         Ccont = 0.216*vol^0.939, # Menden-Deuer & Lessard 2000, all non-diatom protists
         Incubation = as.factor(substring(inc, 3, 3)),
         Treatment = factor(substring(Bag, 1, 1), levels = c("C","D","I","E")),
         Mesocosm = as.factor(substring(Bag, 1, 2)),
         Replicate = substring(Bag, 3, 3)) %>% 
  group_by(Incubation, Treatment, GROUP)

set.seed(2023)
sz.many <- size.dat %>% 
  filter(n()>30) %>% 
  slice_sample(n=30) 

sz.data <- size.dat %>% 
  filter(n()<=30) %>% 
  full_join(., sz.many) %>% 
  summarise(mean.cell.vol = mean(vol),
            sd.cell.vol = sd(vol))
# end of messiness

biom = partial_dat %>% 
  filter(Time_point == "Tend") %>% 
  select(Treatment, Mesocosm, Mes_ID, Replicate, aHF, aMF, aPF) %>% 
  pivot_longer(cols = c(aHF, aPF, aMF), names_to = "group", values_to = "abundance") %>% 
  mutate(GROUP = substring(group, 2, 3))  %>% 
  inner_join(., sz.data, by = c("Treatment", "GROUP")) %>% 
  mutate(biomass = abundance * (0.216*mean.cell.vol^0.939)) # Menden-Deuer Lessard 2000, pgC/mL


## BACTERIA

bact.dat <- read_xlsx("SITES_microscope_data_working.xlsx", 
                      sheet = "0.2_samples",  
                      range = "A1:M109", # change accordingly
                      col_names = T) %>% 
  
    mutate(factor = (5 * Area)/(pi * (21/2)^2),
         HB_abund = HB/(Fields_counted*factor),
         FLB_abund = FLB/(Fields_counted*factor),
         CY_abund = CY/(Fields_counted*factor),
         FLB_perc = FLB_abund/HB_abund,
         Treatment = fct_relevel(Treatment, c("C","D","I","E"))) %>% 
  select(-c(CY, FLB, HB, Fields_counted, factor, Area))

bact.dat$Mesocosm <- as.factor(bact.dat$Mesocosm)
bact.dat$Mes_ID <- as.factor(bact.dat$Mes_ID)


## INGESTION

ingest <- read_xlsx("SITES_microscope_data_working.xlsx", 
                    sheet = "Ingestion",  
                    range = "A1:K5361", # change accordingly
                    col_names = T) 

working <- ingest %>% 
  separate(Sample, into = c("inc", "Bag", "time", "filter"), sep = "_", remove = T) %>% 
  select(-time, -filter) %>% 
  mutate( Incubation = as.factor(substring(inc, 3, 3)),
          Treatment = factor(substring(Bag, 1, 1), levels = c("C","D","I","E")),
          Mesocosm = as.factor(substring(Bag, 1, 2)),
          Replicate = substring(Bag, 3, 3)) %>% 
  group_by(Incubation, Treatment, Mesocosm, Replicate, Time_point, Grazer) %>% 
  # distinct(FLB_presence, Num_FLB, .keep_all = TRUE) 
  # mutate(nomnom = paste0(FLB_presence, "", Num_FLB))
  
  summarise(cells_feeding = sum(FLB_presence),
            flb_ingest = sum(Num_FLB, na.rm = T)) %>% # pay attention here! what's happening with the NAs?
  ungroup() %>% 
  inner_join(., partial_dat %>%
                select(Incubation, Treatment, Mesocosm, Mes_ID, Replicate, Time_point, Cells_counted),
             by = c("Incubation", "Treatment", "Mesocosm", "Replicate", "Time_point")) %>%
  # mutate(flbc = 100*flb_ingest/Cells_counted) %>% 
  # for now
  select(-Cells_counted) %>% 
  #
  pivot_wider(names_from = c(Time_point, Grazer), 
              values_from = c(cells_feeding, flb_ingest
                              # , flbc
                              )) %>% 
  # correct ingestion rates
  mutate(
         # MF_corr1 = cell_num_Tend_MF - cell_num_Tstart_MF,
         # HF_corr1 = cell_num_Tend_HF - cell_num_Tstart_HF,
         MF_corr2 = flb_ingest_Tend_MF - flb_ingest_Tstart_MF,
         HF_corr2 = flb_ingest_Tend_HF - flb_ingest_Tstart_HF,
         MFir = MF_corr2/(cells_feeding_Tend_MF*0.5),
         HFir = HF_corr2/(cells_feeding_Tend_HF*0.5),
         Treatment = fct_relevel(Treatment, c("C","D","I","E"))) %>% 
  select(-contains(c("cell_num_T", "Tstart", "flb_ingest"))) 

working$Mes_ID <- as.factor(working$Mes_ID)

#if I replace with NA it gives NA in the summary
working[sapply(working, is.infinite)] <- 0
# what to do with negative values???

grazing =  left_join(working, bact.dat, by = c("Treatment", "Mes_ID", "Replicate"), copy = T) %>% 
  select(-c(9:20), -FLB_abund, -CY_abund) %>% 
  mutate(mCR = MFir/HB_abund, # in mL
         hCR = HFir/HB_abund)

# for plotting
sumIR <- working %>% 
  pivot_longer(cols = c(MFir, HFir), names_to = "group", values_to = "IR") %>% 
  group_by(Incubation, Treatment, group) %>% 
  summarise(IR.mean = mean(IR, na.rm =T), 
            IR.sd = sd(IR, na.rm =T))


## BASIC PLOTS =======

#abundance
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
  ggplot(., aes(x = Incubation, y = IR, colour = Treatment)) +
  facet_grid(rows = vars(group), 
             labeller = group_ir)+
  # geom_line(aes(group = Treatment))+
  # geom_boxplot(alpha = 0.8) +
  geom_point(aes(color = Mes_ID),
             position = "jitter", 
             alpha = 0.5,
             # shape = 1, stroke = 0.7
  ) +
  scale_color_manual(values = mes.cols)+
  geom_point(aes(color = Treatment), stat = "summary", fun = "mean",size = 4)+
  geom_pointrange(data =sumIR,
                  aes(x = sumIR$Incubation,
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
  filter(Time_point=="Tend")
  # left_join(IR1, bact.dat, by = c("Incubation", "Treatment", "Mes_ID", "Replicate"), copy = T)

ggplot(partial_dat, aes(x= HF))+
  geom_histogram()
ggplot(partial_dat, aes(x= aMF))+
  geom_histogram()
ggplot(partial_dat, aes(x= aPF))+
  geom_histogram()


hm = glmmTMB(aHF ~ 1
             + Treatment 
             + Incubation
             + (1|Mesocosm),
             # offset = Vol,
             family = nbinom2(),
             data = partial_dat %>% 
               filter(Time_point=="Tend"))
summary(hm)

check_model(hm) # why error????

simulationOutput <- simulateResiduals(hm, plot = F)
plot(simulationOutput, quantreg = T)
# NOT UNIFORM RESUDUALS WITH INCUBATION
plotResiduals(simulationOutput, form = stat$Incubation)
plotResiduals(simulationOutput, form = stat$Treatment)
testOutliers(simulationOutput)
testOverdispersion(simulationOutput)

plt <- ggpredict(hm, c( "Incubation", "Treatment"))
plot(plt, add.data = T)  
# ---

pm = glmmTMB(aPF ~ 1
             + Treatment 
             + Incubation
             + (1|Mesocosm),
             # offset = Vol,
             family = gaussian(link = "log"),
             data = partial_dat %>% 
               filter(Time_point=="Tend"))

# gaussian and nbimon don't have a huge difference in AIC etc. Both bad at deviation
summary(pm)

check_model(pm) # why error????

simulationOutput <- simulateResiduals(pm, plot = F)
plot(simulationOutput, quantreg = T)
# NOT UNIFORM RESUDUALS WITH INCUBATION
plotResiduals(simulationOutput, form = stat$Incubation)
plotResiduals(simulationOutput, form = stat$Treatment)
testOutliers(simulationOutput)
testOverdispersion(simulationOutput)

plt <- ggpredict(pm, c( "Incubation", "Treatment"))
plot(plt, add.data = T)  
# ---

# mixotrophs from master and not ingestion
mm = glmmTMB(aMF ~ 1
             + Treatment 
             + Incubation
             + (1|Mesocosm),
             # offset = Vol,
             family = gaussian(link = "log"),
             data = partial_dat %>% 
               filter(Time_point=="Tend"))
summary(mm)

check_model(mm) # why error????

simulationOutput <- simulateResiduals(mm, plot = F)
plot(simulationOutput, quantreg = T)
# NOT UNIFORM RESUDUALS WITH INCUBATION
plotResiduals(simulationOutput, form = stat$Incubation)
plotResiduals(simulationOutput, form = stat$Treatment)
testOutliers(simulationOutput)
testOverdispersion(simulationOutput)
testZeroInflation(simulationOutput)

plt <- ggpredict(mm, c( "Incubation", "Treatment"))
plot(plt, add.data = T) 

# mixotrophs from ingestion and not master
dat <- partial_dat %>% 
  filter(Time_point=="Tend") %>% 
  left_join(., IR1, by = c("Incubation", "Treatment", "Mesocosm", "Mes_ID", "Replicate")) %>% 
  mutate(aMFc = feeding_MF_corr/Vol) %>% 
  select(-c(Time_point, Fields_counted, Cells_counted, Vol, contains("FLB"), contains("feeding")))

dat$aMFc[dat$aMFc<0] <- 0

mm1 = glmmTMB(aMFc ~ 1
             + Treatment 
             + Incubation
             + (1|Mesocosm),
             # offset = Vol,
             family = gaussian(),
             data = dat)
summary(mm1)

check_model(mm1) # why error????

simulationOutput <- simulateResiduals(mm1, plot = F)
plot(simulationOutput, quantreg = T)
# NOT UNIFORM RESUDUALS WITH INCUBATION
plotResiduals(simulationOutput, form = stat$Incubation)
plotResiduals(simulationOutput, form = stat$Treatment)
testOutliers(simulationOutput)
testOverdispersion(simulationOutput)
testZeroInflation(simulationOutput)

plt <- ggpredict(mm1, c( "Incubation", "Treatment"))
plot(plt, add.data = T) 
