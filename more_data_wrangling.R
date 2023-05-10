mydates <- as.Date(c("2022-07-11", "2022-07-19", "2022-07-27"))
mydays <- c(1, 3, 5) # sampling days

absorbance <- read_xlsx("Erken DOM absorbence data.xlsx", 
                        sheet = 1, 
                        range = "A2:G68",
                        col_names = T,
                        col_types = c("date", "guess", "guess", "numeric", "numeric", "numeric", "numeric")
                        ) %>% 
  filter(`Mesocosm ID` %in% dat$Mes_ID )
absorbance[sapply(absorbance, is.na)] <- 2022-07-27





TP_Chla <- read_xlsx("TP Chl Erken and Bolmen.xlsx", 
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
  filter(MesID %in% dat$Mes_ID, Date %in% mydates)


