library(readxl)
library(tidyverse)


V = 10 * (0.1 * 0.1) / (pi * (21/2)^2)
master.dat <- read_xlsx("SITES_microscope_data_sheet.xlsx", 
                        sheet = 2, 
                        range = "A1:T217", 
                        col_names = T)

# for some reason all number columns are characters. Let's convert that 
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
    Mir = FLBinMF_total / (MF * 0.5),
    Hir = FLBinHF_total / (HFwFLB * 0.5),
    Vol = Fields_counted * V,
    # calculate abundances
    aHF = as.integer(as.character(HF_total/Vol)),
    aMF = MF/Vol, 
    aPF = PF/Vol) %>% 
  
  # filter for Incubation 2 and T30 
  
  filter(Incubation == 2 & Time_point == "Tend") %>% 
  drop_na()

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
group_ir <- as_labeller(c(`Hir` = "Heterotrophs",
                             `Mir` = "Mixotrophs"))


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


## Ingestion rate 

# ing.rate <- master.dat %>% 
#   select(-c(Date,Date_analysed,Box,Comments,`Photos/ImageJ`)) %>% 
#   filter(Incubation == 2 & Time_point == "Tend") %>% 
#   pivot_wider()
#   mutate(
#     #maybe I need to do more factor stuff with other columns
#     Treatment = fct_relevel(Treatment, c("C","D","I","E")),
#     Tot_cells = HF_total + PF + MF,
#     # calculate ingestion rates - NOW IT'S COUNTS
#     Mir = FLBinMF_total / (MF * 0.5),
#     Hir = FLBinHF_total / (HFwFLB * 0.5),
#     Vol = Fields_counted * V,
#     # calculate abundances
#     aHF = as.integer(as.character(HF_total/Vol)),
#     aMF = MF/Vol, 
#     aPF = PF/Vol) %>% 
#   
#   # filter for Incubation 2 and T30 
#   
#   drop_na()


trt.mean.ir <- partial_dat %>% 
  pivot_longer(cols = c(Mir, Hir), names_to = "group", values_to = "ingest_rates") %>% 
  group_by(Treatment, group) %>% 
  summarise(trt.mean = mean(ingest_rates), 
            trt.sd = sd(ingest_rates))


mes.mean.ir <- partial_dat %>% 
  pivot_longer(cols = c(Mir, Hir), names_to = "group", values_to = "ingest_rates") %>% 
  group_by(Treatment, Mes_ID, group) %>% 
  summarise(mean = mean(ingest_rates), 
            sd = sd(ingest_rates)) 


partial_dat %>% 
  pivot_longer(cols = c(Mir, Hir), names_to = "group", values_to = "ingest_rates") %>% 
  ggplot(., aes(x = Treatment, y = ingest_rates)) +
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
  geom_pointrange(data =trt.mean.ir,
             aes(x = trt.mean.ir$Treatment,
                 y = trt.mean.ir$trt.mean,
                 ymin=trt.mean.ir$trt.mean-trt.mean.ir$trt.sd,
                 ymax=trt.mean.ir$trt.mean+trt.mean.ir$trt.sd,
                fill = trt.mean.ir$Treatment),
             size = .6, stat = "identity", show.legend = F)+
  # geom_pointrange(stat = "summary", fun.y = mean, 
  #                 fun.ymin = mean-sd, fun.yman = mean+sd)
  scale_fill_manual(values = trt.cols)+
# geom_errorbar(data = trt.mean.ir,
#               aes(ymin=trt.mean.ir$trt.mean-trt.mean.ir$trt.sd,
#                   ymax=trt.mean.ir$trt.mean+trt.mean.ir$trt.sd),
#               width=.2,
#               position = position_dodge(0.05))+
ylab(expression("Ingestion rate"~"bacteria" ~ "cells"^{-1} ~"h"^{-1}))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey96"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))

