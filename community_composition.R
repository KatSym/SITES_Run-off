## Gabriela's data

comm <- read_xlsx("Data/Community_comp_400_countings.xlsx") %>% 
  rename(MesID = Mesocosm, 
         Mesocosm = Treatment) %>% 
  mutate(samplDay = as.numeric(substring(Sample, 3, 3)),
         Treatment = substring(Mesocosm, 1, 1),
         Ciliates = Mesodinium + Rimostrombidium + Ciliates,
         Bacillariophyta = Cyclotella + Stephanodiscus + Acanthoceras,
         Cyanobacteria = Snowella + Microcystis + Woronichinia + Dolichospermum + Dolichospermum2 +
           Merismopedia + Gloeocapsa + Snowella + Chroococcus,
         # Greens = Chlamydomonas + Eudorina+Desmodesmus+Closterium+Cosmarium+Tetraedron+Oocystis+
         #   Ankyra+Koliella+Coelastrum+Chlorophyceae+Ankistrodesmus+Paulschulzsia+Dictyosphaerium+
         #   Kirchneriella+Comasiella+Stauridium+Quadrigula,
        
         # these are classes
          Trebouxiophyceae = Oocystis+Koliella+Dictyosphaerium,
         Chlorophyceae = Chlamydomonas + Eudorina+Desmodesmus+Closterium+Cosmarium+Tetraedron+
           Ankyra+Coelastrum+Chlorophyceae+Ankistrodesmus+Paulschulzsia+
           Kirchneriella+Comasiella+Stauridium+Quadrigula,
         Chrysophyceae = Mallomonas + Ochromonas + Phaester,
         Cryptophyceae = Cryptomonas + Rhodomonas,
         .keep = "unused") %>% 
  rename(Klebsormidiophyceae = Elakatothrix,
         Euglenophyceae = Trachelomonas) %>% 
  # remove Total column becaute the calculation is wrong
  mutate(ExpDay = case_when(samplDay == 0 ~ 0,
                            samplDay == 3 ~ 12,
                            samplDay == 5 ~ 20,
                            samplDay == 9 ~ 36)) %>% 
  select(-No_species, -Total, -samplDay) %>% 
  relocate(c(ExpDay, Treatment), .after = MesID) 






comm %>% 
  # filter(ExpDay %in% c(4, 5, 12, 13, 20, 21)) %>% 
  pivot_longer(cols = Ciliates:Cryptophyceae, names_to = "Taxon", values_to = "count") %>%
  ggplot(., aes(x = as.factor(ExpDay), y = count, fill = Treatment)) +
  # geom_violindot(position = "dodge", stat = "identity")+
  geom_boxplot(position = "dodge2") +
  facet_wrap(vars(Taxon),
             # rows = 2, cols =5, 
             scales = "free"
             )+
  scale_fill_manual(values = trt.cols)+
  # scale_fill_manual(values = c("#E5E5E5", "#D6E6E5", "#FCECEE","#EBF0F9")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black", size = .3),
        # axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11),
        strip.text.y = element_text(size = 12)) +
  labs(x = "Experimental day",
       y = "counts")


comm %>% 
   group_by(ExpDay, Treatment) %>% 
   summarise_all(mean) %>% 
  ungroup() %>% 
   mutate(total = rowSums(.[,4:13])) %>% 
   relocate(total, .after = Treatment) %>% 
  mutate(across(Ciliates:Cryptophyceae, ~./total)) %>% 
   # rowwise() %>%
   # mutate(check = sum(across(Ciliates:Cryptophyceae)))
  
         
  pivot_longer(cols = 5:14, names_to = "Taxon", values_to = "percent") %>%
  select(-total) %>% 
  ggplot(., aes(x = as.factor(ExpDay), y = percent, fill = Taxon)) +
  geom_bar(stat="identity") +
   # geom_col(position = "fill", colour = "black")+
   facet_wrap(~Treatment, ncol = 1, strip.position = "left")+
  scale_fill_brewer(palette = "Set3")+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        # legend.position = "none",
        axis.line = element_line(colour = "black", size = .3),
        # axis.text = element_text(size = 11, colour = "black"),
        axis.title = element_text(size = 11),
        strip.text.y = element_text(size = 12))+ 
  labs(x = "Experimental day",
       y = "relative counts")
