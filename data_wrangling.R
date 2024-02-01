## set up ============
library(tidyverse)
library(readxl)

#needed for the plots 
trt.cols <- c(`C`= "#000000", #black - C
              `D`= "#0a8754", #g - D
              `I`= "#4472ca", #blu - I
              `E`= "#e84855") #r - E


## MICROSCOPE DATA ======

### ABUNDANCE ----
master.dat <- read_xlsx("Data/SITES_microscope_data_working.xlsx", 
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

### SIZE ----
# THIS IS MESSY
# Load the size data. If the number of cells (rows) are more than 30 withing the Incubation - Treatment
# - GROUP grouping, select only 30. This is done in a new df. Then take the <30 rows group and merge them with 
# the new df.
size.dat <- read_xlsx("Data/SITES_microscope_data_working.xlsx", 
                      sheet = "Sizes", 
                      range = "A1:C2053",
                      col_names = T) %>% 
  
  separate(Label, into = c("inc", "Bag", "time", "filter", "image"), sep = "_", remove = T) %>%
  select(-c(time, filter, image)) %>% 
  mutate(letter = rep(c("h", "d"), 1026)) %>% 
  pivot_wider(names_from = letter, values_from = Length) %>% 
  unnest(cols = c(h, d)) %>% 
  mutate(vol = (pi/6)*d^2*h, # prolate spheroid in µm³
         Ccont = 0.216*vol^0.939, # Menden-Deuer & Lessard 2000, all non-diatom protists
         Incubation = as.numeric(substring(inc, 3, 3)),
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


### INGESTION ----

ingest <- read_xlsx("Data/SITES_microscope_data_working.xlsx", 
                    sheet = "Ingestion",  
                    range = "A1:K5361", # change accordingly
                    col_names = T) %>% 
  
  separate(Sample, into = c("inc", "Bag", "time", "filter"), sep = "_", remove = T) %>% 
  select(-time, -filter) %>% 
  mutate( Incubation = as.numeric(substring(inc, 3, 3)),
          Treatment = factor(substring(Bag, 1, 1), levels = c("C","D","I","E")),
          Mesocosm = as.factor(substring(Bag, 1, 2)),
          Replicate = substring(Bag, 3, 3)) %>% 
  group_by(Incubation, Treatment, Mesocosm, Replicate, Time_point, Grazer) %>% 
  
  summarise(cells_feeding = sum(FLB_presence),
            flb_ingest = sum(Num_FLB, na.rm = T)) %>% # pay attention here! what's happening with the NAs?
  ungroup() %>% 
  inner_join(., partial_dat %>%
               select(Incubation, Treatment, Mesocosm, Mes_ID, Replicate, Time_point),
             by = c("Incubation", "Treatment", "Mesocosm", "Replicate", "Time_point")) %>%
  pivot_wider(names_from = c(Time_point, Grazer), 
              values_from = c(cells_feeding, flb_ingest
              )) 
ingest$Mes_ID <- as.factor(ingest$Mes_ID)

# correct ingestion rates
# subtract number of cells with FLB and number of FLB
IR1 <- ingest %>% 
  inner_join(., partial_dat %>%
               filter(Time_point == "Tend") %>% 
               select(Incubation, Treatment, Mesocosm, Mes_ID, Replicate, PF, HF), 
             by = c("Incubation", "Treatment", "Mesocosm", "Replicate", "Mes_ID")) %>% 
  mutate(
    pigm_cells_corr = PF + cells_feeding_Tend_MF - cells_feeding_Tstart_MF,
    het_cells_corr = HF + cells_feeding_Tend_HF - cells_feeding_Tstart_HF,
    feeding_MF_corr = cells_feeding_Tend_MF - cells_feeding_Tstart_MF,
    flb_MF_corr = flb_ingest_Tend_MF - flb_ingest_Tstart_MF,
    flb_HF_corr = flb_ingest_Tend_HF - flb_ingest_Tstart_HF, 
    MFir = flb_MF_corr/(pigm_cells_corr*0.5),
    HFir = flb_HF_corr/(het_cells_corr*0.5),
    Treatment = fct_relevel(Treatment, c("C","D","I","E"))) %>% 
  select(-contains(c("cells", "Tstart"))) 

#convert all infinite and negative values to 0
IR1[sapply(IR1, is.infinite)] <- 0
IR1[sapply(IR1, is.na)] <- 0
IR1$MFir[IR1$MFir<0] <- 0
IR1$HFir[IR1$HFir<0] <- 0

dat <- partial_dat %>% 
  filter(Time_point=="Tend") %>% 
  left_join(., IR1, by = c("Incubation", "Treatment", "Mesocosm", "Mes_ID", "Replicate")) %>% 
  mutate(aMFc = feeding_MF_corr/Vol) %>% 
  select(-c(Time_point, Fields_counted, Cells_counted, Vol, contains("FLB"), contains("feeding")))

dat$aMFc[dat$aMFc<0] <- 0

# get biovolume and biomass
bdat <- dat %>%
  select(Incubation, Treatment, Mes_ID, Replicate, aHF, aMFc, aPF) %>%
  pivot_longer(cols = c(aHF, aPF, aMFc), names_to = "agroup", values_to = "abundance") %>%
  mutate(GROUP = substring(agroup, 2, 3))  %>%
  inner_join(., sz.data, by = c("Incubation", "Treatment", "GROUP")) %>%
  mutate(biovol = abundance * mean.cell.vol, # µm³/mL
         biomass = abundance * (0.216*mean.cell.vol^0.939)) %>%  # Menden-Deuer Lessard 2000, pgC/mL
  select(Incubation, Treatment, Mes_ID, Replicate, GROUP, biomass, biovol) %>%
  pivot_wider(names_from = GROUP, values_from =c(biomass, biovol), values_fn = mean)


### BACTERIA ----

bact.dat <- read_xlsx("Data/SITES_microscope_data_working.xlsx", 
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

# grazing rate
GR1 <- IR1 %>% 
  left_join(., bact.dat, by = c("Incubation", "Treatment", "Mesocosm", "Mes_ID", "Replicate")) %>% 
  mutate(
    Gm = MFir*HB_abund/FLB_abund,
    Gh = HFir*HB_abund/FLB_abund) %>% 
  select(-contains(c("cells", "Tstart")), -Sample, -Time_point, -Date, -FLB_abund) 
GR1[sapply(GR1, is.infinite)] <- 0
GR1[sapply(GR1, is.nan)] <- 0
GR1$Gh[GR1$Gh<0] <- 0
GR1$Gm[GR1$Gm<0] <- 0

# remove(master.dat, partial_dat, size.dat, sz.many, sz.data, ingest, ING)

## BACKGROUND DATA =====

# absorbance
absorbance <- read_xlsx("Data/Erken DOM absorbence data.xlsx", 
                        sheet = 2,
                        range = "A2:H184",
                        col_names = T,
                        col_types = c("date", "guess", "guess", "numeric", 
                                      "numeric", "numeric", "numeric", "numeric")) %>% 
  # "filtered" the dates in Excel - it was horrible, different formats etc
  rename(Mes_ID = `Mesocosm ID`,
         A420_5 = `Abs. 420 (5cm)`,
         A254_5 = `Abs. 254 (5cm)`,
         A420_1 = `Abs. 420 (1 cm)`,
         A254_1 = `Abs. 254 (1 cm)`,
         ExpDay = Day) %>% 
  select(-Date, -`Treatment ID`) %>% 
  mutate(Mes_ID = as.numeric(Mes_ID),
         # use the 5cm cuvette and convert to m
         a420 = 100*A420_5/5)
# replace NAs with the date; it adds 26 to the origin date = 27
# absorbance[sapply(absorbance, is.na)] <- as.Date(26, origin = "2022-07-01") 

# Total phosphorus and Chla
TP_Chla <- read_xlsx("Data/TP Chl Erken and Bolmen.xlsx", 
                     sheet = 1, 
                     range = "A3:M103",
                     col_names = T) %>% 
  select(Provtagningsdatum, Provplats...11, `Totalfosfor EV09 SFM (µg/l)`, `Klorofyll a (µg/l)`) %>% 
  rename(Date = Provtagningsdatum,
         mes = Provplats...11,
         TP = `Totalfosfor EV09 SFM (µg/l)`,
         Chla = `Klorofyll a (µg/l)`) %>% 
  mutate(Date = as.Date(Date)) %>% 
  separate(mes, into = c("text", "Mes_ID")) %>% 
  select(-text) %>% 
  # filter(MesID %in% dat$Mes_ID, Date %in% mydates) %>% 
  # change date to experimental day
  mutate(ExpDay = case_when(Date == "2022-07-07" ~ 0,
                            Date == "2022-07-11" ~ 4,
                            Date == "2022-07-15" ~ 8,
                            Date == "2022-07-19" ~ 12,
                            Date == "2022-07-27" ~ 20,
                            Date == "2022-08-12" ~ 36)) %>% 
  mutate(Mes_ID = as.numeric(Mes_ID))

# Other nutrients
Dnut <- read_xlsx("Data/TOC_TN_DOC_DN_Erken_Bolmen_20230223.xlsx",
                  sheet = 2,
                  range = "A5:F131",
                  col_names = T,
                  col_types = c("numeric", "numeric", "numeric","numeric","numeric","numeric")) %>% 
  filter(complete.cases(.), 
         # `Sample nr` %in% dat$Mes_ID, Day %in% mydays
  ) %>% 
  rename(Mes_ID = `Sample nr`) %>% 
  # sampling day to experimental day
  mutate(ExpDay = case_when(Day == 0 ~ 0,
                            Day == 1 ~ 4,
                            Day == 2 ~ 8,
                            Day == 3 ~ 12,
                            Day == 5 ~ 20,
                            Day == 9 ~ 36))

# dissolved nutrients

dis.nut <- read_xlsx("Data/2022_nutrient_data_summary.xlsx",
                      sheet = 1,
                      range = "A1:Q103",
                      col_names = T) %>% 
  separate(`Mesocosm ID`, into = c("mes", "Mes_ID"), sep = " ") %>% 
  mutate(Mes_ID = as.integer(Mes_ID),
         ExpDay = case_when(`Sampling day` == 0 ~ 0,
                            `Sampling day` == 1 ~ 4,
                            `Sampling day` == 2 ~ 8,
                            `Sampling day` == 3 ~ 12,
                            `Sampling day` == 5 ~ 20,
                            `Sampling day` == 9 ~ 36)) %>% 
  filter(Site == "Erken", 
         Mes_ID != "Lake") %>% 
  rename(TP = `TP (ug/L)`,
         Chla = `Klorofyll a (µg/l)`,
         TOC = `TOC (mg/L)`,
         TN = `TN (mg/L)`,
         DOC = `DOC (mg/L)`,
         dTN = `diss TN (mg/L)`,
         PO4 = `PO4 (ug/L)`,
         NO3 = `NO3 (ug/L)`,
         NO2 = `NO2 (ug/L)`,
         NH4 = `NH4 (ug/L)`) %>% 
select(ExpDay, Mes_ID, TP, Chla, TOC, TN, DOC, dTN, PO4, NO3, NO2, NH4)

# assing something when values are below ditection limit
dis.nut[dis.nut == "<0,5"] <- "0.1"  
dis.nut[dis.nut == "<1"] <- "0.5"
dis.nut[dis.nut == "<3"] <- "2"
dis.nut[dis.nut == "NO RESULT"] <- NA 

# Sensor data
backg <- read.csv("Data/Erken_Daily_avg_final_new_clean.csv",
                  header = T) %>% 
  pivot_longer(cols = 2:117, names_to = "varbl", values_to = "val") %>% 
  separate(varbl, into = c("vrbl", "Mes_ID"), sep = "\\.") %>% 
  pivot_wider(id_cols = c("TIMESTAMP", "Mes_ID"), names_from = "vrbl", values_from = "val") %>% 
  select(-contains("lake")) %>% 
  filter(!is.na(Mes_ID)) %>%   
  rename(PAR = PAR_Apogee,
         DOsat = DOsat_optode,
         Temp = Temp_optode,
         DOconc = DOconc_optode) %>% 
  mutate(
    ExpDay = case_when(TIMESTAMP == "2022-07-07" ~ 0,
                       TIMESTAMP == "2022-07-08" ~ 1,
                       TIMESTAMP == "2022-07-09" ~ 2,
                       TIMESTAMP == "2022-07-10" ~ 3,
                       TIMESTAMP == "2022-07-11" ~ 4,
                       TIMESTAMP == "2022-07-12" ~ 5,
                       TIMESTAMP == "2022-07-13" ~ 6,
                       TIMESTAMP == "2022-07-14" ~ 7,
                       TIMESTAMP == "2022-07-15" ~ 8,
                       TIMESTAMP == "2022-07-16" ~ 9,
                       TIMESTAMP == "2022-07-17" ~ 10,
                       TIMESTAMP == "2022-07-18" ~ 11,
                       TIMESTAMP == "2022-07-19" ~ 12,
                       TIMESTAMP == "2022-07-20" ~ 13,
                       TIMESTAMP == "2022-07-21" ~ 14,
                       TIMESTAMP == "2022-07-22" ~ 15,
                       TIMESTAMP == "2022-07-23" ~ 16,
                       TIMESTAMP == "2022-07-24" ~ 17,
                       TIMESTAMP == "2022-07-25" ~ 18,
                       TIMESTAMP == "2022-07-26" ~ 19,
                       TIMESTAMP == "2022-07-27" ~ 20,
                       TIMESTAMP == "2022-07-28" ~ 21),
    Mes_ID = as.numeric(Mes_ID),
    Treatment = case_when(Mes_ID == 1 | Mes_ID == 7 | Mes_ID == 10 | Mes_ID == 16 ~ "C",
                          Mes_ID == 2 | Mes_ID == 8 | Mes_ID == 11 | Mes_ID == 13 ~ "D",
                          Mes_ID == 3 | Mes_ID == 5 | Mes_ID == 12 | Mes_ID == 14 ~ "I",
                          Mes_ID == 4 | Mes_ID == 6 | Mes_ID == 9 | Mes_ID == 15 ~ "E"),
    .keep = "unused"
  ) 

## ROTIFERS ====

rotifers <- read_xlsx("Data/rotifers_volvoxes.xlsx",
                      sheet = 1, range = "A1:K107", col_names = T) %>% 
  rename(ExpDay = `day of experiment`,
         Treatment = treatment,
         vol = `sample volume (ml)`) %>% 
  filter(mesocosm!="LE") %>% 
  mutate(
    # taking rotifer colonies as 1 individual, per L
    total = 1000*(`small rotifers`+ asplanchna + kellicottia + `rotifer colony`)/vol,
    # small rotifers are probably Keratella
    small = 1000*`small rotifers`/vol,
    # rotifer colonies are probably Conochilus 
    coln = 1000*`rotifer colony`/vol,
    Kell = 1000*kellicottia/vol,
    Aspl = 1000*asplanchna/vol,
    Mes_ID = as.numeric(mesocosm)) %>%  
  pivot_longer(cols = c(small, coln, Kell, Aspl), 
               names_to = "rotif.sp", values_to = "abund") %>% 
  select(ExpDay, Treatment, Mes_ID, rotif.sp, abund)


rot.meas <- read_xlsx("Data/Erken_zoo_measurements_2016-2019.xlsx", 
                      sheet = 1, range = "A1:O2023", col_names = T) %>% 
  select(-c(Site, Comment1, `Dyntaxa ID`, `Analysis method`, `Analysis laboratory`)) %>% 
  # keeping only rotifers in the summer months
  mutate(month = substr(as.character(Date), 3, 4)) %>% 
  filter(month %in% c("06", "07", "08"),
         Phylum == "Rotifera") %>% 
  # get the "volume" not biovolume
  mutate(indv.vol = `Biovol.(µm³/l)`/ `Density(n/l)`) %>% 
  # keeping only the genus
  separate(ScientificName, c("Genus", NA), remove = T) %>% 
  group_by(Genus) %>% 
  summarise(indv.vol.m = mean(indv.vol),
            indv.vol.sd = sd(indv.vol),
            biovol.m = mean(`Biovol.(µm³/l)`),
            biovol.sd = sd(`Biovol.(µm³/l)`),
            length.m = mean(`Mesured length (µm)`),
            length.sd = sd(`Mesured length (µm)`),
            FW.m = mean(`Freshweight(mg/l)`),
            FW.sd = sd(`Freshweight(mg/l)`)) %>% 
  as.data.frame()

rot <-  rotifers %>% 
  mutate(rot.vol = case_when(rotif.sp == "small" ~ 
                               rot.meas[rot.meas$Genus == "Keratella", "indv.vol.m"],
                             rotif.sp == "Aspl" ~ 
                               rot.meas[rot.meas$Genus == "Asplanchna", "indv.vol.m"],
                             rotif.sp == "coln" ~ 
                               rot.meas[rot.meas$Genus == "Conochilus", "indv.vol.m"],
                             rotif.sp == "Kell" ~ 
                               rot.meas[rot.meas$Genus == "Kellicottia", "indv.vol.m"]),
         # biovolume in µm³/l
         biovol = abund * rot.vol) 

## CREATE DATA FILES ====
# put it all together


envir <- 
   # # TP, Chla, NO3, NO2, PO4 in ug/L, other nutr in mg/L
  full_join(dis.nut, absorbance, by = c("Mes_ID", "ExpDay")) %>%
  full_join(., backg, by = c("Mes_ID", "ExpDay")) %>%
  mutate(Treatment = case_when(Mes_ID == 1 | Mes_ID == 7 | Mes_ID == 10 | Mes_ID == 16 ~ "C",
                               Mes_ID == 2 | Mes_ID == 8 | Mes_ID == 11 | Mes_ID == 13 ~ "D",
                               Mes_ID == 3 | Mes_ID == 5 | Mes_ID == 12 | Mes_ID == 14 ~ "I",
                               Mes_ID == 4 | Mes_ID == 6 | Mes_ID == 9 | Mes_ID == 15 ~ "E"),
         suva = (A254_1/DOC)*100) %>% # L/mg*m
  # select(-contains("id", ignore.case = F), -contains("Date")) %>% 
  # filter(!is.na(Mes_ID)) %>% 
  arrange(ExpDay) %>% 
  # keep only relevant mesocosms and time period
  filter(Mes_ID %in% c(1, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16),
         ExpDay<=21) %>% 
  select(-c(A254_5, A420_5, A254_1, A420_1, Trilux_Neph, Trilux_Chloro, Trilux_Phycocy)) %>% 
  as.data.frame()
#### NOTE! The env data don't match with my data on experimental day. They are a day before,
#### so change it to be the same depending on the analysis
## mutate(ExpDay = case_when(ExpDay == 4 ~ 5,
##                           ExpDay == 12 ~ 13,
##                           ExpDay == 20 ~ 21)) %>% 
## filter(!is.na(ExpDay))  

data = dat %>% 
  mutate(ExpDay = case_when(Incubation == 1 ~ 5,
                            Incubation == 2 ~ 13,
                            Incubation == 3 ~ 21)) %>% 
  left_join(., GR1, by = c("Incubation", "Treatment", "Mes_ID", "Replicate")) %>% 
  left_join(., bdat, by = c("Incubation", "Treatment", "Mes_ID", "Replicate")) %>% 
  select(-contains(".x"), -contains("Tend"), -contains("corr"), -id, -Sample) %>% 
  mutate(Mes_ID = as.numeric(as.character(Mes_ID))) %>%
  rename(MF_Ir = MFir.y,
         HF_Ir = HFir.y, 
         PF_abund = aPF,
         HF_abund = aHF,
         MF_abund = aMFc,
         MF_Gr = Gm,
         HF_Gr = Gh,
         Mesocosm = Mesocosm.y) %>% 
  as.data.frame() %>% 
  select(-contains(".y"), -c(MF, aMF, PF, HF, Replicate, Incubation, Mesocosm)) %>% 
  relocate(ExpDay, .before = Treatment)

# save the above data frame and the colours so you don't have to run all this every time
save(data, envir, rot, trt.cols, file = "all_data.RData")

# write csvs
# write.csv(data, "Data/microscope_data_topub.csv", row.names = F)
# write.csv(envir, "Data/envir_data_topub.csv", row.names = F)
# write.csv(rot, "Data/rotifer_data_topub.csv", row.names = F)
