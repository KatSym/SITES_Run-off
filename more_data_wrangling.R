mydates <- as.Date(c("2022-07-11", "2022-07-19", "2022-07-27"))
mydays <- c(1, 3, 5) # sampling days

absorbance <- read_xlsx("Data/Erken DOM absorbence data.xlsx", 
                        sheet = "absorbance", 
                        range = "A2:G50",
                        col_names = T,
                        col_types = c("date", "guess", "guess", "numeric", "numeric", "numeric", "numeric")
                        ) %>% 
  # "filtered" the dates in Excel - it was horrible here
  filter(`Mesocosm ID` %in% dat$Mes_ID ) %>% 
  rename(MesID = `Mesocosm ID`,
         A420_5 = `Abs. 420 (5cm)`,
         A254_5 = `Abs. 254 (5cm)`,
         A420_1 = `Abs. 420 (1 cm)`,
         A254_1 = `Abs. 254 (1 cm)`,
         Mesocosm = `Treatment ID`) %>% 
  mutate(id = row_number())
# replace NAs with the date; it adds 26 to the origin date = 27
absorbance[sapply(absorbance, is.na)] <- as.Date(26, origin = "2022-07-01") 



Treats <- read_xlsx("Data/SITES_experiment_planning.xlsx",
                    sheet = "Treatments",
                    range = "A1:D21", 
                    col_names = T) %>% 
  pivot_longer(cols = c("E", "I", "D"), names_to = "Treatment", values_to = "intensity")



Treats %>% 
  ggplot(., aes(x = Experimental_day, y = intensity, colour = Treatment))+
  geom_bar()


TP_Chla <- read_xlsx("Data/TP Chl Erken and Bolmen.xlsx", 
                     sheet = 1, 
                     range = "A3:M103",
                     col_names = T,
                    ) %>% 
  select(Provtagningsdatum, Provplats...11, `Totalfosfor EV09 SFM (µg/l)`, `Klorofyll a (µg/l)`) %>% 
  rename(Date = Provtagningsdatum,
         mes = Provplats...11,
         TP = `Totalfosfor EV09 SFM (µg/l)`,
         Chla = `Klorofyll a (µg/l)`) %>% 
  mutate(Date = as.Date(Date)) %>% 
  separate(mes, into = c("text", "MesID")) %>% 
  select(-text) %>% 
  filter(MesID %in% dat$Mes_ID, Date %in% mydates) %>% 
  mutate(id = row_number(),
         MesID = as.numeric(MesID))

Dnut <- read_xlsx("Data/TOC_TN_DOC_DN_Erken_Bolmen_20230223.xlsx",
                  sheet = 2,
                  range = "A5:F131",
                  col_names = T,
                  col_types = c("numeric", "guess", "guess","guess","guess","guess")) %>% 
  filter(complete.cases(.), `Sample nr` %in% dat$Mes_ID, Day %in% mydays) %>% 
  rename(MesID = `Sample nr`) %>% 
  mutate(id = row_number())

fcBact <- read_xlsx("Data/Erken_bacterial_abundance_flow_cytometer.xlsx",
                    range = "A1:H161",
                    col_names = T) %>% 
  filter(Sampling %in% mydays) %>% 
  select(Sample_ID, Sampling, bacterial_cells_per_ml) %>% 
  separate(Sample_ID, into = c("smt", "code", "smte"), sep = "-", remove = T) %>% 
  mutate(Mesocosm = substring(code, 1, 2)) %>% 
  select(-c(smt, code, smte)) %>% 
  filter(Mesocosm %in% levels(dat$Mesocosm)) %>% 
  rename(Day = Sampling)


env <- full_join(Dnut, TP_Chla, absorbance, by = c("MesID", "id")) %>% 
  full_join(., absorbance, by = c("MesID", "id")) %>% 
  full_join(., fcBact, by = c("Day", "Mesocosm")) %>% 
  mutate(Treatment = substring(Mesocosm, 1, 1)) %>% 
  select(-id, -Date.y) %>% 
  rename(Bact = bacterial_cells_per_ml)
# TP & Chla in ug/L
# other nutr in mg/L
