mydates <- as.Date(c("2022-07-11", "2022-07-19", "2022-07-27"))
mydays <- c(1, 3, 5) # sampling days

absorbance <- read_xlsx("Data/Erken DOM absorbence data.xlsx", 
                        # sheet = "absorbance",
                        sheet = 2,
                        # range = "A2:H50",
                        range = "A2:H184",
                        col_names = T,
                        col_types = c("date", "guess", "guess", "numeric", "numeric", "numeric", "numeric", "numeric")
                        ) %>% 
  # "filtered" the dates in Excel - it was horrible here
  # filter(`Mesocosm ID` %in% dat$Mes_ID ) %>% 
  rename(MesID = `Mesocosm ID`,
         A420_5 = `Abs. 420 (5cm)`,
         A254_5 = `Abs. 254 (5cm)`,
         A420_1 = `Abs. 420 (1 cm)`,
         A254_1 = `Abs. 254 (1 cm)`,
         Mesocosm = `Treatment ID`,
         ExpDay = Day) %>% 
  mutate(MesID = as.numeric(MesID))
# replace NAs with the date; it adds 26 to the origin date = 27
# absorbance[sapply(absorbance, is.na)] <- as.Date(26, origin = "2022-07-01") 



Treats <- read_xlsx("Data/SITES_experiment_planning.xlsx",
                    sheet = "Treatments",
                    range = "A1:D21", 
                    col_names = T) %>% 
  pivot_longer(cols = c("E", "I", "D"), names_to = "Treatment", values_to = "intensity")



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
  # change date to experimental day
  mutate(ExpDay = case_when(Date == "2022-07-07" ~ 0,
                            Date == "2022-07-11" ~ 4,
                            Date == "2022-07-15" ~ 8,
                            Date == "2022-07-19" ~ 12,
                            Date == "2022-07-27" ~ 20,
                            Date == "2022-08-12" ~ 36)) %>% 
  mutate(MesID = as.numeric(MesID))

Dnut <- read_xlsx("Data/TOC_TN_DOC_DN_Erken_Bolmen_20230223.xlsx",
                  sheet = 2,
                  range = "A5:F131",
                  col_names = T,
                  col_types = c("numeric", "numeric", "numeric","numeric","numeric","numeric")) %>% 
  filter(complete.cases(.), 
         # `Sample nr` %in% dat$Mes_ID, Day %in% mydays
         ) %>% 
  rename(MesID = `Sample nr`) %>% 
  # sampling day to experimental day
    mutate(ExpDay = case_when(Day == 0 ~ 0,
                            Day == 1 ~ 4,
                            Day == 2 ~ 8,
                            Day == 3 ~ 12,
                            Day == 5 ~ 20,
                            Day == 9 ~ 36))

fcBact <- read_xlsx("Data/Erken_bacterial_abundance_flow_cytometer.xlsx",
                    range = "A1:H161",
                    col_names = T) %>% 
  # filter(Sampling %in% mydays) %>% 
  select(Sample_ID, Sampling, bacterial_cells_per_ml) %>% 
  separate(Sample_ID, into = c("smt", "code", "smte"), sep = "-", remove = T) %>% 
  mutate(Mesocosm = substring(code, 1, 2),
         ExpDay = case_when(Sampling == 0 ~ 0,
                            Sampling == 1 ~ 4,
                            Sampling == 2 ~ 8,
                            Sampling == 3 ~ 12,
                            Sampling == 4 ~ 16,
                            Sampling == 5 ~ 20,
                            Sampling == 6 ~ 24,
                            Sampling == 7 ~ 28,
                            Sampling == 8 ~ 32,
                            Sampling == 9 ~ 36)) %>% 
  select(-c(smt, code, smte)) 
  # filter(Mesocosm %in% levels(dat$Mesocosm)) %>% 


# datetoday <- c("2022-07-07" ~ 0,"2022-07-08"~1,"2022-07-09"~2,"2022-07-10"~3,"2022-07-11"~4,
#                "2022-07-12"~5,"2022-07-13"~6,"2022-07-14"~7,"2022-07-15"~8,"2022-07-16"~9,
#                "2022-07-17"~10,"2022-07-18"~11,"2022-07-19"~12,"2022-07-20"~13,"2022-07-21"~14,
#                "2022-07-22"~15,"2022-07-23"~16,"2022-07-24"~17,"2022-07-25"~18,"2022-07-26"~19,
#                "2022-07-27"~20,"2022-07-28"~21,"2022-07-29"~22,"2022-07-30"~23,"2022-07-31"~24,
#                "2022-08-01"~25,"2022-08-02"~26,"2022-08-03"~27,"2022-08-04"~28,"2022-08-05"~29,
#                "2022-08-06"~29,"2022-08-07"~30,"2022-08-08"~31,"2022-08-09"~32,"2022-08-10"~33,
#                "2022-08-11"~34,"2022-08-12"~35)

backg <- read.csv("Data/Erken_Daily_avg_final_new_clean.csv",
                  header = T) %>% 
  pivot_longer(cols = 2:117, names_to = "varbl", values_to = "val") %>% 
  separate(varbl, into = c("vrbl", "MesID"), sep = "\\.") %>% 
  filter(TIMESTAMP %in% c("2022-07-11", "2022-07-12", "2022-07-19","2022-07-20", "2022-07-27", "2022-07-28")) %>% 
  pivot_wider(id_cols = c("TIMESTAMP", "MesID"), names_from = "vrbl", values_from = "val") %>% 
  select(-contains("lake")) %>% 
  filter(!is.na(MesID)) %>%   
  rename(PAR = PAR_Apogee,
         DOsat = DOsat_optode,
         Temp = Temp_optode,
         DOconc = DOconc_optode) %>% 
  mutate(ExpDay = case_when(TIMESTAMP == "2022-07-11" ~ 4,
                            TIMESTAMP == "2022-07-12" ~ 5,
                            TIMESTAMP == "2022-07-19" ~ 12,
                            TIMESTAMP == "2022-07-20" ~ 13,
                            TIMESTAMP == "2022-07-27" ~ 20,
                            TIMESTAMP == "2022-07-28" ~ 21),
         MesID = as.numeric(MesID),
         Treatment = case_when(MesID == 1 | MesID == 7 | MesID == 10 | MesID == 16 ~ "C",
                               MesID == 2 | MesID == 8 | MesID == 11 | MesID == 13 ~ "D",
                               MesID == 3 | MesID == 5 | MesID == 12 | MesID == 14 ~ "I",
                               MesID == 4 | MesID == 6 | MesID == 9 | MesID == 15 ~ "E")
         ) 

# arrange(TIMESTAMP) %>% 
  #   mutate(ExpDay = 0:n(TIMESTAMP))




env <- full_join(Dnut, TP_Chla, by = c("MesID", "ExpDay")) %>% 
  full_join(., absorbance, by = c("MesID", "ExpDay")) %>%
  full_join(., backg, by = c("MesID", "ExpDay")) %>%
  mutate(Treatment = substring(Mesocosm, 1, 1),
         Treatment = case_when(MesID == 1 | MesID == 7 | MesID == 10 | MesID == 16 ~ "C",
                               MesID == 2 | MesID == 8 | MesID == 11 | MesID == 13 ~ "D",
                               MesID == 3 | MesID == 5 | MesID == 12 | MesID == 14 ~ "I",
                               MesID == 4 | MesID == 6 | MesID == 9 | MesID == 15 ~ "E"),
         suva = (A254_1/DOC)*100) %>% # L/mg*m
  select(-contains("id", ignore.case = F), -contains("Date"), -Day, -TIMESTAMP, -Mesocosm) %>% 
  rename(Mes_ID = MesID) %>% 
  filter(!is.na(Mes_ID)) %>% 
  arrange(ExpDay) %>% 
  as.data.frame()
# TP & Chla in ug/L
# other nutr in mg/L



## plots ---------

# absorbance %>% 
#   mutate(Treatment = substring(Mesocosm, 1, 1)) %>% 
#   ggplot(., aes(x = Day, y = A254_1, colour = Treatment)) +
#   # geom_smooth(aes(group = Treatment, colour = Treatment, fill = Treatment)) +
#   geom_point(
#     aes(x = Day, 
#         y = A420_1, 
#         colour = Treatment), 
#     position = position_jitter(width = .02),
#     alpha = .5) +
#   # facet_grid(varbl~.,
#   #            scales = "free_y")+
#   # scale_x_continuous(breaks = c(0, 1, 2, 3, 5, 9),
#   #                    labels = c(0, 4, 8, 12, 20, 36)
#   # ) +
#   scale_color_manual(values = trt.cols)+
#   scale_fill_manual(values = c("#E5E5E5", "#D6E6E5", "#FCECEE","#EBF0F9")) +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.background = element_blank(),
#         legend.position = "none",
#         axis.line = element_line(colour = "black", size = .3),
#         axis.text = element_text(size = 11, colour = "black"),
#         axis.title = element_text(size = 11),
#         strip.text.y = element_text(size = 12)) +
#   labs(x = "Experimental day") 

# Treats %>%
#   ggplot(., aes(x = Experimental_day, y = intensity, fill = Treatment))+
#   geom_bar(stat = "identity", position=position_dodge(), width = .6,
#            # alpha = .3
#            ) +
#   scale_fill_manual(values = trt.cols, 
#                     labels = c("Daily", "Intermediate", "Extreme")) +
#   scale_x_continuous(limits = c(0, 21),
#                      expand = c(0, 0),
#                      breaks = c(0, 5, 10, 13, 20, 21)) +
#   labs(y = "Addition %",
#        x = "Experimental day") +
# theme(panel.grid.minor = element_blank(),
#       panel.grid.major = element_blank(),
#       panel.background = element_rect(fill='transparent'),
#       plot.background = element_rect(fill='transparent', color=NA),
#       axis.text = element_text(size = 11, , colour = "black"),
#       axis.title = element_text(size = 11),
#       legend.position = "bottom",
#       legend.title = element_blank(),
#       legend.background = element_rect(fill='transparent'),
#       # axis.text.y = element_text(size = 12),
#       strip.text.y = element_text(size = 12),
#       axis.text.x = element_text(colour = c('black', 'darkred','black', 'darkred','black', 'darkred' ),
#                                  face = c("plain", "bold", "plain", "bold", "plain", "bold"))) 
# ggsave("Plots/20230526_treatments.png", dpi = 300)
# 

# TP_Chla %>% 
#   mutate(Treatment = case_when(MesID == 1 | MesID == 7 | MesID == 10 | MesID == 16 ~ "C",
#                                MesID == 2 | MesID == 8 | MesID == 11 | MesID == 13 ~ "D",
#                                MesID == 3 | MesID == 5 | MesID == 12 | MesID == 14 ~ "I",
#                                MesID == 4 | MesID == 6 | MesID == 9 | MesID == 15 ~ "E")) %>% 
#   ggplot(., aes(x = Date, y = Chla, colour = Treatment)) +
#   geom_smooth(aes(group = Treatment, colour = Treatment, fill = Treatment)) +
#   geom_point(
#     aes(x = Date, 
#         y = Chla, 
#         colour = Treatment), 
#     position = position_jitter(width = .02),
#     alpha = .5) +
#   scale_x_continuous(
#     breaks = c("2022-07-07", "2022-07-11", "2022-07-15", "2022-07-19", ),
#                      labels = c(0, 4, 8, 12, 20, 36)
#   ) +
#   scale_color_manual(values = trt.cols)+
#   scale_fill_manual(values = c("#E5E5E5", "#D6E6E5", "#FCECEE","#EBF0F9")) +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.background = element_blank(),
#         legend.position = "none",
#         axis.line = element_line(colour = "black", size = .3),
#         axis.text = element_text(size = 11, colour = "black"),
#         axis.title = element_text(size = 11),
#         strip.text.y = element_text(size = 12)) +
#   labs(x = "Experimental day",
#        # y = "mg L\u207b\u2081"
#        y = "") 
# 
# # nutrient plot - comment the filter on the dataframes
# nuts <- left_join(Dnut, TP_Chla, by = c("Day", "MesID")) %>% 
#   # filter(complete.cases(.)) %>% 
#   mutate(Treatment = case_when(MesID == 1 | MesID == 7 | MesID == 10 | MesID == 16 ~ "C",
#                                MesID == 2 | MesID == 8 | MesID == 11 | MesID == 13 ~ "D",
#                                MesID == 3 | MesID == 5 | MesID == 12 | MesID == 14 ~ "I",
#                                MesID == 4 | MesID == 6 | MesID == 9 | MesID == 15 ~ "E"),
#          # ExDay =case_when(Day == 1 ~ 4,)
#          ) %>% 
#   pivot_longer(cols = c("DOC", "DN", "TP"), names_to = "varbl", values_to = "val") 
# # %>% 
# #   mutate(val = log1p(val)) %>% 
# #   group_by(Day, Treatment, MesID, varbl) %>% 
# #   summarise_all(mean)


# nuts %>% 
#   ggplot(., aes(x = Day, y = val, colour = Treatment)) +
#   geom_smooth(aes(group = Treatment, colour = Treatment, fill = Treatment)) +
#   geom_point(
#     aes(x = Day, 
#         y = val, 
#         colour = Treatment), 
#     position = position_jitter(width = .02),
#     alpha = .5) +
#   facet_grid(varbl~.,
#              scales = "free_y")+
#   scale_x_continuous(breaks = c(0, 1, 2, 3, 5, 9),
#     labels = c(0, 4, 8, 12, 20, 36)
#     ) +
#   scale_color_manual(values = trt.cols)+
#   scale_fill_manual(values = c("#E5E5E5", "#D6E6E5", "#FCECEE","#EBF0F9")) +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.background = element_blank(),
#         legend.position = "none",
#         axis.line = element_line(colour = "black", size = .3),
#         axis.text = element_text(size = 11, colour = "black"),
#         axis.title = element_text(size = 11),
#         strip.text.y = element_text(size = 12)) +
#   labs(x = "Experimental day",
#        # y = "mg L\u207b\u2081"
#        y = "") 
# ggsave("Plots/20230531_nutrients.png", dpi = 300)