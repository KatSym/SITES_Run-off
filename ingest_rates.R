# ingestion rates
library(ggnewscale)

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
    cells_MF_corr = cells_feeding_Tend_MF - cells_feeding_Tstart_MF,
    cells_HF_corr = cells_feeding_Tend_HF - cells_feeding_Tstart_HF,
    flb_MF_corr = flb_ingest_Tend_MF - flb_ingest_Tstart_MF,
    flb_HF_corr = flb_ingest_Tend_HF - flb_ingest_Tstart_HF,
    MFir = flb_MF_corr/(cells_MF_corr*0.5),
    HFir = flb_HF_corr/(cells_HF_corr*0.5),
    Treatment = fct_relevel(Treatment, c("C","D","I","E"))) %>% 
  select(-contains(c("cells", "Tstart", "flb"))) 

#if I replace with NA it gives NA in the summary
IR1[sapply(IR1, is.infinite)] <- 0
IR1[sapply(IR1, is.na)] <- 0


# subtract number of FLB only
# IR2 <- ING %>% 
# mutate(
#   flb_MF_corr = flb_ingest_Tend_MF - flb_ingest_Tstart_MF,
#   flb_HF_corr = flb_ingest_Tend_HF - flb_ingest_Tstart_HF,
#   MFir = flb_MF_corr/(cells_feeding_Tend_MF*0.5),
#   HFir = flb_HF_corr/(cells_feeding_Tend_HF*0.5),
#   Treatment = fct_relevel(Treatment, c("C","D","I","E"))) %>% 
#   select(-contains(c("cells", "Tstart", "flb"))) 
# IR2[sapply(IR2, is.infinite)] <- 0
# IR2[sapply(IR2, is.na)] <- 0



IR1 %>% 
  pivot_longer(cols = c(MFir, HFir), names_to = "group", values_to = "IR") %>% 
  ggplot(., aes(x = Incubation, y = IR, colour = Treatment)) +
  facet_grid(rows = vars(group), 
             labeller = group_ir)+
  # geom_line(aes(group = Treatment))+
  # geom_boxplot(alpha = 0.8) +
  geom_point(aes(color = Mes_ID),
             position = "jitter", 
             alpha = 0.5,
             # shape = 1, stroke = 0.7,
             show.legend = F
  ) +
  scale_colour_manual(values = mes.cols, )+
  new_scale_colour()+
  geom_point(aes(color = Treatment), stat = "summary", fun = "mean", size = 4, alpha = .6)+
  scale_colour_manual(values = trt.cols)+
  ylab(expression("Ingestion rate"~"bacteria" ~ "cells"^{-1} ~"h"^{-1}))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "grey96"),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size= 12))


library(glmmTMB)
library(ggeffects)
library(performance)
library(DHARMa)

ggplot(IR1, aes(x= HFir))+
  geom_histogram()
# should I include interaction or not?
m.ir1 = glmmTMB(MFir ~ 1
             + Treatment 
             + Incubation
             + (1|Mes_ID)
             ,
             # offset = Vol,
             family = gaussian(),
             data = IR1)
summary(m.ir1)

# DHARMa workfow
simulationOutput <- simulateResiduals(m.ir1, plot = F)
plot(simulationOutput)
# NOT UNIFORM RESUDUALS WITH INCUBATION
plotResiduals(simulationOutput, form = IR1$Incubation)
plotResiduals(simulationOutput, form = IR1$Treatment)

testOverdispersion(simulationOutput)
testZeroInflation(simulationOutput) # I don't know what that means

check_model(m.ir1) 
plt <- ggpredict(m.ir1, c( "Incubation", "Treatment"))
plot(plt, add.data = T)             


IR1$HFir[IR1$HFir == 8] <- NA

h.ir1 = glmmTMB(HFir ~ 1
                + Treatment 
                + Incubation
                + (1|Mes_ID)
                ,
                # offset = Vol,
                family = gaussian(),
                data = IR1)
summary(h.ir1)
check_model(m.ir1) 

plt <- ggpredict(h.ir1, c( "Incubation", "Treatment"))
plot(plt, add.data = T)  
