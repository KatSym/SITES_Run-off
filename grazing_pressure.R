library(tidyverse)

load("all_data.RData")

rot1 = rot %>% 
  filter(Mes_ID %in% c(1, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16),
         ExpDay != 0, ExpDay < 21) %>% 
  mutate(ExpDay = case_when(ExpDay == 4 ~ 5,
                            ExpDay == 12 ~ 13,
                            ExpDay == 20 ~ 21)) %>% 
  select(ExpDay, Treatment, biovol) %>%   
  rename(rot.biovol = biovol) %>% 
  group_by(ExpDay, Treatment) %>% 
  summarise_all(mean)


flg.bio = data %>% 
  select(ExpDay, Treatment, contains("biovol")) %>% 
  group_by(ExpDay, Treatment) %>% 
  summarise_all(mean) %>% 
  left_join(., rot1, by = c("ExpDay", "Treatment")) %>% 
  mutate(grPF = rot.biovol/biovol_PF, 
         grMF = rot.biovol/biovol_MF,
         grHF = rot.biovol/biovol_HF)
