library(tidyverse)

load("all_data.RData")

rot1 = rot %>% 
  filter(Mes_ID %in% c(1, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16),
         ExpDay != 0, ExpDay < 21,
         # exclude Asplanchna because it's too big
         rotif.sp != "Aspl"
         ) %>% 
  # change the days so comparisons are possible 
  mutate(ExpDay = case_when(ExpDay == 4 ~ 5,
                            ExpDay == 12 ~ 13,
                            ExpDay == 20 ~ 21),
         biomass = (62.57 * rot.vol^1.455) * abund ) %>% # Telesh 1998 table 3 and µg C/L
  select(ExpDay, Treatment, rot.vol, biovol, biomass) %>%   #
  rename(rot.biovol = biovol,
         rot.biomass = biomass)  


flg.bio <- data %>% 
  rowwise() %>% 
  mutate(fl.biov = sum(HF_biovol, MF_biovol, PF_biovol) * 1e3,
         fl.biom = sum(HF_biovol, MF_biovol, PF_biovol) * 1e-3) %>%  # pg C/mL to µg C/L
  select(ExpDay, Treatment, contains("fl.")) %>% 
  full_join(., rot1, by = c("ExpDay", "Treatment")) %>% 
  group_by(ExpDay, Treatment) %>%
  summarise_all(mean)


# Hansen et al. 1997 Table 9 ---
# calculate a from their numbers and then Imax and Cmax with our volume
LIng.a = -0.929 + 0.23*5.75
Imax = function(x) {10^(LIng.a - 0.23*log10(x))}

LC.a = -0.330 + 0.23*5.75
Cmax = 10^(LC.a - 0.23*log10(rot1$rot.vol))

# v is the volume filtered by the rotifers in an hour
# the 10e5 * the complicated power of 10 is the Cmax
# the 1e-12 is the conversion of µm³ to mL 
# the 24 to make it per day
r1 = flg.bio %>% mutate(
  # v = rot.biovol * 1e5 * 10^(LC.a - 0.23*log10(rot.vol)) * 1e-12 * 24,  # per L
                    i.bm = rot.biomass * Imax(rot.vol) * 24 * 1e-3,
                    perc.grz.bm = i.bm/fl.biom, # pg C/mL to µg C/L
                    # .keep = "unused"
                    i.bv = rot.biovol * Imax(rot.vol) * 24,
                    perc.grz.bv = 100*i.bv/fl.biov) # pg C/mL to µg C/L) %>% 
  group_by(ExpDay, Treatment) %>% 
  summarise_all(mean)

  

### NOT USED
# growth =  envir %>% 
#     mutate(ExpDay = case_when(ExpDay == 4 ~ 5,
#                             ExpDay == 12 ~ 13,
#                             ExpDay == 20 ~ 21,
#                             .default = ExpDay),
#            mu = 0.399*(PO4/(0.607+PO4))) %>% # Istvanovics et al 1994, Monod model and table III
#   filter(ExpDay %in% c(5, 13, 21),
#          !is.na(PO4)) %>% 
#   group_by(ExpDay, Treatment) %>% 
#   summarise_all(mean) %>%
#   full_join(.,r1, by = c("ExpDay", "Treatment")) %>% 
#   select(ExpDay, Treatment, flg.biom, biomass, mu, v, i) %>% 
#          mutate(grw = mu * flg.biom)  
# # need to match experimental days
# 
# 
# flg.bio = data %>% 
#   select(ExpDay, Treatment, contains("biomass")) %>% 
#   rowwise() %>% 
#   mutate(tot.fl = sum(HF_biovol, MF_biovol, PF_biovol), .keep = "unused") %>% 
#   select(ExpDay, Treatment, tot.fl) %>% 
#   group_by(ExpDay, Treatment) %>% 
#   summarise_all(mean) %>% 
#   pivot_wider(id_cols = "Treatment", names_from = c("ExpDay"), values_from = "tot.fl") %>% 
#   mutate(mu1 = (`13`-`5`)/(8*`13`),
#          mu2 = (`21`-`13`)/(8*`21`),
#          mu = (`21`-`5`)/(16*`21`))



ggplot(r1, aes(x = ExpDay, y = perc.grz.bv, colour = Treatment))+
  geom_point()+
  scale_color_manual(values = trt.cols)+
  labs(y = "<span style='font-size: 15pt'>Volume filtered by rotifers</span>
         <span style='font-size: 12pt'>(mL)</span>",
       x = "Experimental day")+
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = .3),
    axis.text = element_text(size = 13, colour = "black"),
    axis.title = element_text(size = 15),
    axis.title.y = ggtext::element_markdown()
    # axis.text.y = element_text(size = 12),
    # strip.text.y = element_text(size = 12)
  ) 


ggplot(rot %>% 
         filter(ExpDay < 22)
       , aes(x = ExpDay, y = log10(abund+1), fill = rotif.sp, group = rotif.sp))+
  geom_bar(stat="summary", width = 2.2)+
  # stat summary appliying a mean function but desplaying correct sum????
  facet_wrap(~Treatment)+
  # geom_bar(stat="identity", position = "dodge")+
  scale_fill_viridis_d(labels = c("Asplanchna", "Rot. colonies", "Kellicottia", "small rotifers")
                    )+
  labs(y = "<span style='font-size: 15pt'>Rotifer abundance</span>
         <span style='font-size: 12pt'>log(x+1) ind/L</span>",
       x = "Experimental day")+
  theme(
    # legend.position = "none",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = .3),
    axis.text = element_text(size = 13, colour = "black"),
    axis.title = element_text(size = 15),
    axis.title.y = ggtext::element_markdown()
    # axis.text.y = element_text(size = 12),
    # strip.text.y = element_text(size = 12)
  ) 

ggsave("Plots/Jan2024/rotif-comm-comp3.png", dpi=300, bg = "white")
