## set up ============
library(tidyverse)
library(readxl)


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


## INGESTION

ingest <- read_xlsx("SITES_microscope_data_working.xlsx", 
                    sheet = "Ingestion",  
                    range = "A1:K5361", # change accordingly
                    col_names = T) 

ING <- ingest %>% 
  separate(Sample, into = c("inc", "Bag", "time", "filter"), sep = "_", remove = T) %>% 
  select(-time, -filter) %>% 
  mutate( Incubation = as.numeric(substring(inc, 3, 3)),
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
              )) 
ING$Mes_ID <- as.factor(ING$Mes_ID)

# correct ingestion rates
# subtract number of cells with FLB and number of FLB
IR1 <- ING %>% 
  mutate(
    feeding_MF_corr = cells_feeding_Tend_MF - cells_feeding_Tstart_MF,
    feeding_HF_corr = cells_feeding_Tend_HF - cells_feeding_Tstart_HF,
    flb_MF_corr = flb_ingest_Tend_MF - flb_ingest_Tstart_MF,
    flb_HF_corr = flb_ingest_Tend_HF - flb_ingest_Tstart_HF,
    MFir = flb_MF_corr/(feeding_MF_corr*0.5),
    HFir = flb_HF_corr/(feeding_HF_corr*0.5),
    Treatment = fct_relevel(Treatment, c("C","D","I","E"))) %>% 
  select(-contains(c("cells", "Tstart"))) 

#if I replace with NA it gives NA in the summary
IR1[sapply(IR1, is.infinite)] <- 0
IR1[sapply(IR1, is.na)] <- 0

dat <- partial_dat %>% 
  filter(Time_point=="Tend") %>% 
  left_join(., IR1, by = c("Incubation", "Treatment", "Mesocosm", "Mes_ID", "Replicate")) %>% 
  mutate(aMFc = feeding_MF_corr/Vol) %>% 
  select(-c(Time_point, Fields_counted, Cells_counted, Vol, contains("FLB"), contains("feeding")))

dat$aMFc[dat$aMFc<0] <- 0

bdat <- dat %>% 
  select(Treatment, Mes_ID, Replicate, aHF, aMFc, aPF) %>% 
  pivot_longer(cols = c(aHF, aPF, aMFc), names_to = "agroup", values_to = "abundance") %>% 
  mutate(GROUP = substring(agroup, 2, 3))  %>% 
  inner_join(., sz.data, by = c("Treatment", "GROUP")) %>% 
  mutate(biomass = abundance * (0.216*mean.cell.vol^0.939)) %>%  # Menden-Deuer Lessard 2000, pgC/mL
  select(Incubation, Treatment, Mes_ID, Replicate, GROUP, biomass) %>% 
  pivot_wider(names_from = GROUP, values_from = biomass, values_fn = mean)
dataa$MF[is.na(dataa$MF)] <- 0
