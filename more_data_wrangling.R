mydates <- as.Date(c("2022-07-11", "2022-07-19", "2022-07-27"))
mydays <- c(1, 3, 5) # sampling days

absorbance <- read_xlsx("Data/Erken DOM absorbence data.xlsx", 
                        sheet = "absorbance", 
                        range = "A2:G50",
                        col_names = T,
                        col_types = c("date", "guess", "guess", "numeric", "numeric", "numeric", "numeric")
                        ) %>% 
  # "filtered" the dates in Excel - it was horrible here
  # filter(`Mesocosm ID` %in% dat$Mes_ID ) %>% 
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
  ggplot(., aes(x = Experimental_day, y = intensity, fill = Treatment))+
  geom_bar(stat = "identity", position=position_dodge(), width = .6) +
  scale_fill_manual(values = trt.cols) +
  scale_x_continuous(limits = c(4.5, 21),
                     # expand = c(0, 0),
                     breaks = c( 5, 13, 21)) +
  labs(y = "Addition %",
       x = "Experimental day") +
theme(panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      axis.text = element_text(size = 11, , colour = "black"),
      axis.title = element_text(size = 11),
      # axis.text.y = element_text(size = 12),
      strip.text.y = element_text(size = 12),
      axis.text.x = element_text(colour = c('black', 'black','black', 'darkred','black', 'darkred' ),
                                 face = c("plain", "plain", "plain", "bold", "plain", "bold"))) 
ggsave("Plots/20230523_treatments_PARTIAL.png", dpi = 300)






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
  # filter(MesID %in% dat$Mes_ID, Date %in% mydates) %>% 
  mutate(Day = case_when(Date == "2022-07-07" ~ 0,
                         Date == "2022-07-11" ~ 1,
                         Date == "2022-07-15" ~ 2,
                         Date == "2022-07-19" ~ 3,
                         Date == "2022-07-27" ~ 5)) %>% 
  mutate(id = row_number(),
         MesID = as.numeric(MesID))

Dnut <- read_xlsx("Data/TOC_TN_DOC_DN_Erken_Bolmen_20230223.xlsx",
                  sheet = 2,
                  range = "A5:F131",
                  col_names = T,
                  col_types = c("numeric", "numeric", "numeric","numeric","numeric","numeric")) %>% 
  filter(complete.cases(.), 
         # `Sample nr` %in% dat$Mes_ID, Day %in% mydays
         ) %>% 
  rename(MesID = `Sample nr`) %>% 
  mutate(id = row_number())

fcBact <- read_xlsx("Data/Erken_bacterial_abundance_flow_cytometer.xlsx",
                    range = "A1:H161",
                    col_names = T) %>% 
  # filter(Sampling %in% mydays) %>% 
  select(Sample_ID, Sampling, bacterial_cells_per_ml) %>% 
  separate(Sample_ID, into = c("smt", "code", "smte"), sep = "-", remove = T) %>% 
  mutate(Mesocosm = substring(code, 1, 2)) %>% 
  select(-c(smt, code, smte)) %>% 
  # filter(Mesocosm %in% levels(dat$Mesocosm)) %>% 
  rename(Day = Sampling)


env <- full_join(Dnut, TP_Chla, by = c("MesID", "id")) %>% 
  full_join(., absorbance, by = c("MesID", "id")) %>%
  full_join(., fcBact, by = c("Day", "Mesocosm")) %>%
  mutate(Treatment = substring(Mesocosm, 1, 1)) %>% 
  select(-id, -Date.y) %>% 
  rename(Bact = bacterial_cells_per_ml,
         Date = Date.x) %>% 
  as.data.frame()
# TP & Chla in ug/L
# other nutr in mg/L

# nutrient plot - comment the filter on the dataframes
nuts <- left_join(Dnut, TP_Chla, by = c("Day", "MesID")) %>% 
  filter(complete.cases(.)) %>% 
  mutate(Treatment = case_when(MesID == 1 | MesID == 7 | MesID == 10 | MesID == 16 ~ "C",
                               MesID == 2 | MesID == 8 | MesID == 11 | MesID == 13 ~ "D",
                               MesID == 3 | MesID == 5 | MesID == 12 | MesID == 14 ~ "I",
                               MesID == 4 | MesID == 6 | MesID == 9 | MesID == 15 ~ "E")) %>% 
  pivot_longer(cols = c("TOC", "TN", "TP"), names_to = "varbl", values_to = "val") %>% 
  mutate(val = log1p(val)) %>% 
  group_by(Day, Treatment, MesID, varbl) %>% 
  summarise_all(mean)

nuts %>% 
  ggplot(., aes(x = Day, y = val, colour = Treatment)) +
  geom_line(aes(group = varbl, linetype = varbl)) +
  facet_grid(Treatment~.)
  
  
