library(tidyverse)
library(readxl)

## Rotifers from Ben

rotifers <- read_xlsx("Data/rotifers_volvoxes.xlsx",
                      sheet = 1, range = "A1:K107", col_names = T) %>% 
  rename(ExpDay = `day of experiment`,
         Treatment = treatment,
         vol = `sample volume (ml)`) %>% 
  filter(mesocosm!="LE",
         ExpDay<21) %>% 
  mutate(
    # taking rotifer colonies as 1 individual, per L
    Rotif = 1000*(`small rotifers`+ asplanchna + kellicottia + `rotifer colony`)/vol,
    
    # small rotifers are probably Keratella
    small = 1000*`small rotifers`/vol,
    
    # rotifer colonies are probably Conochilus 
    coln = 1000*`rotifer colony`/vol,
    
    Kell = 1000*kellicottia/vol,
    Aspl = 1000*asplanchna/vol,
    Mes_ID = as.numeric(mesocosm)) %>%  
  pivot_longer(cols = c(small, coln, Kell, Aspl), 
               names_to = "rotif.sp", values_to = "abund") %>% 
  select(ExpDay, Treatment, Mes_ID, rotif.sp, abund, Rotif)

ggplot(rotifers, aes(x = as.factor(ExpDay), y = log(abund+1), fill = rotif.sp))+
  geom_bar(width = 0.8, stat = "identity", position = "dodge")+
  facet_wrap(~Treatment)

## Zooplankton from Benny

zoopl <- read_xlsx("Data/ZOOPLANKTON_FINAL.xlsx", 
                   sheet = 1, range = "A1:T981", col_names = T) %>% 
  rename(Mes_ID = Mesocosm) %>% 
  filter(Lake == "Erken",
         Treatment != "LE") %>% 
  mutate(
    ExpDay = case_when(Date == 20220706 ~ 0, # it's actually day -1
                       Date == 20220726 ~ 19,
                       Date == 20220810 ~ 34),
    # corrections
    Species = case_when(Species == "bosmian" ~ "bosmina",
                        Species == "small rotifer" ~ "S_rotifer",
                        Species == "large rotifer" ~ "L_rotifer",
                        .default = Species),
    # per L
    abund = 1000*(((Counted * Concentrated_Volume_ml)/Counted_Volume_ml) + Picked) / Volume_ml
  ) %>% 
  select(Mes_ID,Treatment, ExpDay, Sample, Group, Species, abund)

ggplot(zoopl, aes(x = as.factor(ExpDay), y = abund, fill = Group))+
  geom_bar(width = 0.8, stat = "identity", position = "dodge")+
  facet_wrap(~Treatment)

## combination to do other stuff
rot1 <- rotifers %>% 
  group_by(ExpDay,Treatment, Mes_ID) %>% 
  summarise_all(mean) %>% 
  select(-rotif.sp, -abund)

zoop1 <- zoopl %>% 
  # excluding rotifers - TO BE DISCUSSED
  filter(Group != "Rotifer",
         ExpDay < 30) %>% 
  group_by(ExpDay,Treatment, Mes_ID) %>% 
  summarise_all(mean) %>% 
  select(-Sample, -Species, -Group) %>% 
  rename(Zoopl = abund)

zoop <- full_join(rot1, zoop1, by = c("ExpDay", "Treatment", "Mes_ID")) %>% 
  filter(Mes_ID %in% c(1, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16)) %>% 
  ungroup() 
