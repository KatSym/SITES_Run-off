


example <- data.frame(
                      Date = rep(c("20220712", "20220720", "20220728"), each = 108),
                      Experiment = rep(c(1, 2, 3), each = 108),
                      Mesocosm = rep(rep(c("C1", "I1","E1", "E2","C2","D2", "E3", "D3", "I3",  "D4","I4", "C4"), 
                                     each=3), 3),
                      Replicate = rep(c("a", "b", "c"), 108),
                      Filter = rep(rep(c(0.2, 0.8), times = c(36, 72))),
                      Time_point = rep(c("T0", "T30"), times = c(72, 36))
   )

write.csv(example, "Erken_FLB_incubation_master_202207.csv", row.names = T)
getwd()






LUDAT <- data.frame(
                      Date = rep("20220722", 33),
                      Experiment = rep("LU", 33),
                      Mesocosm = rep("LU", 33),
                      Replicate = rep(c("a", "b", "c"), 33),
                      Time_point = rep(c("Incubation_start", "T0","T02", "T05", "T10", "T15", "T20", "T30",
                                         "T40", "T60", "T100"), each = 3)
                    )
write.csv(LUDAT, "LU.csv", row.names = T)
