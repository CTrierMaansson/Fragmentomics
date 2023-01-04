mer_count_NAC.1_cfChIP <- mer_count(NAC.1_cfChIP,3)
mer_count_NAC.2_cfChIP <- mer_count(NAC.2_cfChIP,3)
mer_count_NAC.3_cfChIP <- mer_count(NAC.3_cfChIP,3)
mer_count_NAC.4_cfChIP <- mer_count(NAC.4_cfChIP,3)
mer_count_NSC.1_cfChIP <- mer_count(NSC.1_cfChIP,3)
mer_count_NSC.2_cfChIP <- mer_count(NSC.2_cfChIP,3)
mer_count_NSC.3_cfChIP <- mer_count(NSC.3_cfChIP,3)
mer_count_NSC.4_cfChIP <- mer_count(NSC.4_cfChIP,3)
mer_count_SSC.1_cfChIP <- mer_count(SSC.1_cfChIP,3)
mer_count_SSC.2_cfChIP <- mer_count(SSC.2_cfChIP,3)
mer_count_SSC.3_cfChIP <- mer_count(SSC.3_cfChIP,3)
mer_count_SSC.4_cfChIP <- mer_count(SSC.4_cfChIP,3)
mer_count_HC.1_cfChIP <- mer_count(HC.1_cfChIP, 3)
mer_count_HC.2_cfChIP <- mer_count(HC.2_cfChIP, 3)
mer_count_HC.3_cfChIP <- mer_count(HC.3_cfChIP, 3)
mer_count_HC.4_cfChIP <- mer_count(HC.4_cfChIP, 3)
gc()

mer_count_NAC.1_input <- mer_count(NAC.1_input,3)
mer_count_NAC.2_input <- mer_count(NAC.2_input,3)
mer_count_NAC.3_input <- mer_count(NAC.3_input,3)
mer_count_NAC.4_input <- mer_count(NAC.4_input,3)
mer_count_NSC.1_input <- mer_count(NSC.1_input,3)
mer_count_NSC.2_input <- mer_count(NSC.2_input,3)
mer_count_NSC.3_input <- mer_count(NSC.3_input,3)
mer_count_NSC.4_input <- mer_count(NSC.4_input,3)
mer_count_SSC.1_input <- mer_count(SSC.1_input,3)
mer_count_SSC.2_input <- mer_count(SSC.2_input,3)
mer_count_SSC.3_input <- mer_count(SSC.3_input,3)
mer_count_SSC.4_input <- mer_count(SSC.4_input,3)
mer_count_HC.1_input <- mer_count(HC.1_input, 3)
mer_count_HC.2_input <- mer_count(HC.2_input, 3)
mer_count_HC.3_input <- mer_count(HC.3_input, 3)
mer_count_HC.4_input <- mer_count(HC.4_input, 3)
gc()

prop_mer_NAC.1_cfChIP <- prop_df(mer_count_NAC.1_cfChIP)
prop_mer_NAC.2_cfChIP <- prop_df(mer_count_NAC.2_cfChIP)
prop_mer_NAC.3_cfChIP <- prop_df(mer_count_NAC.3_cfChIP)
prop_mer_NAC.4_cfChIP <- prop_df(mer_count_NAC.4_cfChIP)
prop_mer_NSC.1_cfChIP <- prop_df(mer_count_NSC.1_cfChIP)
prop_mer_NSC.2_cfChIP <- prop_df(mer_count_NSC.2_cfChIP)
prop_mer_NSC.3_cfChIP <- prop_df(mer_count_NSC.3_cfChIP)
prop_mer_NSC.4_cfChIP <- prop_df(mer_count_NSC.4_cfChIP)
prop_mer_SSC.1_cfChIP <- prop_df(mer_count_SSC.1_cfChIP)
prop_mer_SSC.2_cfChIP <- prop_df(mer_count_SSC.2_cfChIP)
prop_mer_SSC.3_cfChIP <- prop_df(mer_count_SSC.3_cfChIP)
prop_mer_SSC.4_cfChIP <- prop_df(mer_count_SSC.4_cfChIP)
prop_mer_HC.1_cfChIP <- prop_df(mer_count_HC.1_cfChIP)
prop_mer_HC.2_cfChIP <- prop_df(mer_count_HC.2_cfChIP)
prop_mer_HC.3_cfChIP <- prop_df(mer_count_HC.3_cfChIP)
prop_mer_HC.4_cfChIP <- prop_df(mer_count_HC.4_cfChIP)

prop_mer_NAC.1_input <- prop_df(mer_count_NAC.1_input)
prop_mer_NAC.2_input <- prop_df(mer_count_NAC.2_input)
prop_mer_NAC.3_input <- prop_df(mer_count_NAC.3_input)
prop_mer_NAC.4_input <- prop_df(mer_count_NAC.4_input)
prop_mer_NSC.1_input <- prop_df(mer_count_NSC.1_input)
prop_mer_NSC.2_input <- prop_df(mer_count_NSC.2_input)
prop_mer_NSC.3_input <- prop_df(mer_count_NSC.3_input)
prop_mer_NSC.4_input <- prop_df(mer_count_NSC.4_input)
prop_mer_SSC.1_input <- prop_df(mer_count_SSC.1_input)
prop_mer_SSC.2_input <- prop_df(mer_count_SSC.2_input)
prop_mer_SSC.3_input <- prop_df(mer_count_SSC.3_input)
prop_mer_SSC.4_input <- prop_df(mer_count_SSC.4_input)
prop_mer_HC.1_input <- prop_df(mer_count_HC.1_input)
prop_mer_HC.2_input <- prop_df(mer_count_HC.2_input)
prop_mer_HC.3_input <- prop_df(mer_count_HC.3_input)
prop_mer_HC.4_input <- prop_df(mer_count_HC.4_input)
library(dplyr)
prop_mer_cfChIP <- prop_mer_NAC.1_cfChIP %>% 
    left_join(prop_mer_NAC.2_cfChIP, by = "res") %>% 
    left_join(prop_mer_NAC.3_cfChIP, by = "res") %>% 
    left_join(prop_mer_NAC.4_cfChIP, by = "res") %>% 
    left_join(prop_mer_NSC.1_cfChIP, by = "res") %>% 
    left_join(prop_mer_NSC.2_cfChIP, by = "res") %>% 
    left_join(prop_mer_NSC.3_cfChIP, by = "res") %>% 
    left_join(prop_mer_NSC.4_cfChIP, by = "res") %>% 
    left_join(prop_mer_SSC.1_cfChIP, by = "res") %>% 
    left_join(prop_mer_SSC.2_cfChIP, by = "res") %>% 
    left_join(prop_mer_SSC.3_cfChIP, by = "res") %>% 
    left_join(prop_mer_SSC.4_cfChIP, by = "res") %>% 
    left_join(prop_mer_HC.1_cfChIP, by = "res") %>%
    left_join(prop_mer_HC.2_cfChIP, by = "res") %>%
    left_join(prop_mer_HC.3_cfChIP, by = "res") %>%
    left_join(prop_mer_HC.4_cfChIP, by = "res")
prop_mer_input <- prop_mer_NAC.1_input %>% 
    left_join(prop_mer_NAC.2_input, by = "res") %>% 
    left_join(prop_mer_NAC.3_input, by = "res") %>% 
    left_join(prop_mer_NAC.4_input, by = "res") %>% 
    left_join(prop_mer_NSC.1_input, by = "res") %>% 
    left_join(prop_mer_NSC.2_input, by = "res") %>% 
    left_join(prop_mer_NSC.3_input, by = "res") %>% 
    left_join(prop_mer_NSC.4_input, by = "res") %>% 
    left_join(prop_mer_SSC.1_input, by = "res") %>% 
    left_join(prop_mer_SSC.2_input, by = "res") %>% 
    left_join(prop_mer_SSC.3_input, by = "res") %>% 
    left_join(prop_mer_SSC.4_input, by = "res") %>% 
    left_join(prop_mer_HC.1_input, by = "res") %>%
    left_join(prop_mer_HC.2_input, by = "res") %>%
    left_join(prop_mer_HC.3_input, by = "res") %>%
    left_join(prop_mer_HC.4_input, by = "res")

colnames(prop_mer_cfChIP) <- c("motif", 
                               "NAC.1_cfChIP",
                               "NAC.2_cfChIP",
                               "NAC.3_cfChIP",
                               "NAC.4_cfChIP",
                               "NSC.1_cfChIP",
                               "NSC.2_cfChIP",
                               "NSC.3_cfChIP",
                               "NSC.4_cfChIP",
                               "SSC.1_cfChIP",
                               "SSC.2_cfChIP",
                               "SSC.3_cfChIP",
                               "SSC.4_cfChIP",
                               "HC.1_cfChIP",
                               "HC.2_cfChIP",
                               "HC.3_cfChIP",
                               "HC.4_cfChIP")
colnames(prop_mer_input) <- c("motif", 
                              "NAC.1_input",
                              "NAC.2_input",
                              "NAC.3_input",
                              "NAC.4_input",
                              "NSC.1_input",
                              "NSC.2_input",
                              "NSC.3_input",
                              "NSC.4_input",
                              "SSC.1_input",
                              "SSC.2_input",
                              "SSC.3_input",
                              "SSC.4_input",
                              "HC.1_input",
                              "HC.2_input",
                              "HC.3_input",
                              "HC.4_input")
prop_mer_input
setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver")

write.table(prop_mer_cfChIP, file = "fragment end motif proportions cfChIP.txt", 
            sep = "\t",col.names = T)
write.table(prop_mer_input, file = "fragment end motif proportions input.txt", 
            sep = "\t",col.names = T)

count_mer_NAC.1_cfChIP <- count_df(mer_count_NAC.1_cfChIP)
count_mer_NAC.2_cfChIP <- count_df(mer_count_NAC.2_cfChIP)
count_mer_NAC.3_cfChIP <- count_df(mer_count_NAC.3_cfChIP)
count_mer_NAC.4_cfChIP <- count_df(mer_count_NAC.4_cfChIP)
count_mer_NSC.1_cfChIP <- count_df(mer_count_NSC.1_cfChIP)
count_mer_NSC.2_cfChIP <- count_df(mer_count_NSC.2_cfChIP)
count_mer_NSC.3_cfChIP <- count_df(mer_count_NSC.3_cfChIP)
count_mer_NSC.4_cfChIP <- count_df(mer_count_NSC.4_cfChIP)
count_mer_SSC.1_cfChIP <- count_df(mer_count_SSC.1_cfChIP)
count_mer_SSC.2_cfChIP <- count_df(mer_count_SSC.2_cfChIP)
count_mer_SSC.3_cfChIP <- count_df(mer_count_SSC.3_cfChIP)
count_mer_SSC.4_cfChIP <- count_df(mer_count_SSC.4_cfChIP)
count_mer_HC.1_cfChIP <- count_df(mer_count_HC.1_cfChIP)
count_mer_HC.2_cfChIP <- count_df(mer_count_HC.2_cfChIP)
count_mer_HC.3_cfChIP <- count_df(mer_count_HC.3_cfChIP)
count_mer_HC.4_cfChIP <- count_df(mer_count_HC.4_cfChIP)

count_mer_NAC.1_input <- count_df(mer_count_NAC.1_input)
count_mer_NAC.2_input <- count_df(mer_count_NAC.2_input)
count_mer_NAC.3_input <- count_df(mer_count_NAC.3_input)
count_mer_NAC.4_input <- count_df(mer_count_NAC.4_input)
count_mer_NSC.1_input <- count_df(mer_count_NSC.1_input)
count_mer_NSC.2_input <- count_df(mer_count_NSC.2_input)
count_mer_NSC.3_input <- count_df(mer_count_NSC.3_input)
count_mer_NSC.4_input <- count_df(mer_count_NSC.4_input)
count_mer_SSC.1_input <- count_df(mer_count_SSC.1_input)
count_mer_SSC.2_input <- count_df(mer_count_SSC.2_input)
count_mer_SSC.3_input <- count_df(mer_count_SSC.3_input)
count_mer_SSC.4_input <- count_df(mer_count_SSC.4_input)
count_mer_HC.1_input <- count_df(mer_count_HC.1_input)
count_mer_HC.2_input <- count_df(mer_count_HC.2_input)
count_mer_HC.3_input <- count_df(mer_count_HC.3_input)
count_mer_HC.4_input <- count_df(mer_count_HC.4_input)

count_mer_cfChIP <- count_mer_NAC.1_cfChIP %>% 
    left_join(count_mer_NAC.2_cfChIP, by = "res") %>% 
    left_join(count_mer_NAC.3_cfChIP, by = "res") %>% 
    left_join(count_mer_NAC.4_cfChIP, by = "res") %>% 
    left_join(count_mer_NSC.1_cfChIP, by = "res") %>% 
    left_join(count_mer_NSC.2_cfChIP, by = "res") %>% 
    left_join(count_mer_NSC.3_cfChIP, by = "res") %>% 
    left_join(count_mer_NSC.4_cfChIP, by = "res") %>% 
    left_join(count_mer_SSC.1_cfChIP, by = "res") %>% 
    left_join(count_mer_SSC.2_cfChIP, by = "res") %>% 
    left_join(count_mer_SSC.3_cfChIP, by = "res") %>% 
    left_join(count_mer_SSC.4_cfChIP, by = "res") %>% 
    left_join(count_mer_HC.1_cfChIP, by = "res") %>%
    left_join(count_mer_HC.2_cfChIP, by = "res") %>%
    left_join(count_mer_HC.3_cfChIP, by = "res") %>%
    left_join(count_mer_HC.4_cfChIP, by = "res")
count_mer_input <- count_mer_NAC.1_input %>% 
    left_join(count_mer_NAC.2_input, by = "res") %>% 
    left_join(count_mer_NAC.3_input, by = "res") %>% 
    left_join(count_mer_NAC.4_input, by = "res") %>% 
    left_join(count_mer_NSC.1_input, by = "res") %>% 
    left_join(count_mer_NSC.2_input, by = "res") %>% 
    left_join(count_mer_NSC.3_input, by = "res") %>% 
    left_join(count_mer_NSC.4_input, by = "res") %>% 
    left_join(count_mer_SSC.1_input, by = "res") %>% 
    left_join(count_mer_SSC.2_input, by = "res") %>% 
    left_join(count_mer_SSC.3_input, by = "res") %>% 
    left_join(count_mer_SSC.4_input, by = "res") %>% 
    left_join(count_mer_HC.1_input, by = "res") %>%
    left_join(count_mer_HC.2_input, by = "res") %>%
    left_join(count_mer_HC.3_input, by = "res") %>%
    left_join(count_mer_HC.4_input, by = "res")
colnames(count_mer_cfChIP) <- c("motif", 
                                "NAC.1_cfChIP",
                                "NAC.2_cfChIP",
                                "NAC.3_cfChIP",
                                "NAC.4_cfChIP",
                                "NSC.1_cfChIP",
                                "NSC.2_cfChIP",
                                "NSC.3_cfChIP",
                                "NSC.4_cfChIP",
                                "SSC.1_cfChIP",
                                "SSC.2_cfChIP",
                                "SSC.3_cfChIP",
                                "SSC.4_cfChIP",
                                "HC.1_cfChIP",
                                "HC.2_cfChIP",
                                "HC.3_cfChIP",
                                "HC.4_cfChIP")
colnames(count_mer_input) <- c("motif", 
                               "NAC.1_input",
                               "NAC.2_input",
                               "NAC.3_input",
                               "NAC.4_input",
                               "NSC.1_input",
                               "NSC.2_input",
                               "NSC.3_input",
                               "NSC.4_input",
                               "SSC.1_input",
                               "SSC.2_input",
                               "SSC.3_input",
                               "SSC.4_input",
                               "HC.1_input",
                               "HC.2_input",
                               "HC.3_input",
                               "HC.4_input")
write.table(count_mer_cfChIP, file = "fragment end motif count cfChIP.txt", 
            sep = "\t",col.names = T)
write.table(count_mer_input, file = "fragment end motif count input.txt", 
            sep = "\t",col.names = T)

NAC.1_active_inactive_df <- end_motif_active_inactive(enrichment_df,
                                                      grs,
                                                      "NAC.1",
                                                      "D:/Lung cancer input/PosDeduped/PosDeduped-A-1279-input.bam",
                                                      3)
NAC.2_active_inactive_df <- end_motif_active_inactive(enrichment_df,
                                                      grs,
                                                      "NAC.2",
                                                      "D:/Lung cancer input/PosDeduped/PosDeduped-B-1288-input.bam",
                                                      3)
NAC.3_active_inactive_df <- end_motif_active_inactive(enrichment_df,
                                                      grs,
                                                      "NAC.3",
                                                      "D:/Lung cancer input/PosDeduped/PosDeduped-C-1475-input.bam",
                                                      3)

NAC.4_active_inactive_df <- end_motif_active_inactive(enrichment_df,
                                                      grs,
                                                      "NAC.4",
                                                      "D:/Lung cancer input/PosDeduped/PosDeduped-D-1578-input.bam",
                                                      3)

NSC.1_active_inactive_df <- end_motif_active_inactive(enrichment_df,
                                                      grs,
                                                      "NSC.1",
                                                      "D:/Lung cancer input/PosDeduped/PosDeduped-E-439-input.bam",
                                                      3)

NSC.2_active_inactive_df <- end_motif_active_inactive(enrichment_df,
                                                      grs,
                                                      "NSC.2",
                                                      "D:/Lung cancer input/PosDeduped/PosDeduped-F-1449-input.bam",
                                                      3)

NSC.3_active_inactive_df <- end_motif_active_inactive(enrichment_df,
                                                      grs,
                                                      "NSC.3",
                                                      "D:/Lung cancer input/PosDeduped/PosDeduped-I-645-input.bam",
                                                      3)

NSC.4_active_inactive_df <- end_motif_active_inactive(enrichment_df,
                                                      grs,
                                                      "NSC.4",
                                                      "D:/Lung cancer input/PosDeduped/PosDeduped-J-1663-input.bam",
                                                      3)

SSC.1_active_inactive_df <- end_motif_active_inactive(enrichment_df,
                                                      grs,
                                                      "SSC.1",
                                                      "D:/Lung cancer input/PosDeduped/PosDeduped-G-514-input.bam",
                                                      3)

SSC.2_active_inactive_df <- end_motif_active_inactive(enrichment_df,
                                                      grs,
                                                      "SSC.2",
                                                      "D:/Lung cancer input/PosDeduped/PosDeduped-H-1169-input.bam",
                                                      3)

SSC.3_active_inactive_df <- end_motif_active_inactive(enrichment_df,
                                                      grs,
                                                      "SSC.3",
                                                      "D:/Lung cancer input/PosDeduped/PosDeduped-K-440-input.bam",
                                                      3)

SSC.4_active_inactive_df <- end_motif_active_inactive(enrichment_df,
                                                      grs,
                                                      "SSC.4",
                                                      "D:/Lung cancer input/PosDeduped/PosDeduped-L-1100-input.bam",
                                                      3)
HC.1_active_inactive_df <- end_motif_active_inactive(x = enrichment_df,
                                                     y = grs,
                                                     z = "HC.1",
                                                     i = "D:/Healthy input/PosDeduped/PosDeduped-Rask_kontrol_1_input.bam",
                                                     m = 3)
HC.2_active_inactive_df <- end_motif_active_inactive(x = enrichment_df,
                                                     y = grs,
                                                     z = "HC.2",
                                                     i = "D:/Healthy input/PosDeduped/PosDeduped-Rask_kontrol_2_input.bam",
                                                     m = 3)
HC.3_active_inactive_df <- end_motif_active_inactive(x = enrichment_df,
                                                     y = grs,
                                                     z = "HC.3",
                                                     i = "D:/Healthy input/PosDeduped/PosDeduped-Rask_kontrol_3_input.bam",
                                                     m = 3)
HC.4_active_inactive_df <- end_motif_active_inactive(x = enrichment_df,
                                                     y = grs,
                                                     z = "HC.4",
                                                     i = "D:/Healthy input/PosDeduped/PosDeduped-Rask_kontrol_4_input.bam",
                                                     m = 3)

collected_active_inactive_df <- NAC.1_active_inactive_df %>% 
    left_join(NAC.2_active_inactive_df, by = "motif") %>%
    left_join(NAC.3_active_inactive_df, by = "motif") %>%
    left_join(NAC.4_active_inactive_df, by = "motif") %>%
    left_join(NSC.1_active_inactive_df, by = "motif") %>%
    left_join(NSC.2_active_inactive_df, by = "motif") %>%
    left_join(NSC.3_active_inactive_df, by = "motif") %>%
    left_join(NSC.4_active_inactive_df, by = "motif") %>%
    left_join(SSC.1_active_inactive_df, by = "motif") %>%
    left_join(SSC.2_active_inactive_df, by = "motif") %>%
    left_join(SSC.3_active_inactive_df, by = "motif") %>%
    left_join(SSC.4_active_inactive_df, by = "motif") %>% 
    left_join(HC.1_active_inactive_df, by = "motif") %>%
    left_join(HC.2_active_inactive_df, by = "motif") %>%
    left_join(HC.3_active_inactive_df, by = "motif") %>%
    left_join(HC.4_active_inactive_df, by = "motif")

colnames(collected_active_inactive_df) <- c("motif",
                                            "NAC.1_active", "NAC.1_inactive",
                                            "NAC.2_active", "NAC.2_inactive",
                                            "NAC.3_active", "NAC.3_inactive",
                                            "NAC.4_active", "NAC.4_inactive",
                                            "NSC.1_active", "NSC.1_inactive",
                                            "NSC.2_active", "NSC.2_inactive",
                                            "NSC.3_active", "NSC.3_inactive",
                                            "NSC.4_active", "NSC.4_inactive",
                                            "SSC.1_active", "SSC.1_inactive",
                                            "SSC.2_active", "SSC.2_inactive",
                                            "SSC.3_active", "SSC.3_inactive",
                                            "SSC.4_active", "SSC.4_inactive",
                                            "HC.1_active", "HC.1_inactive",
                                            "HC.2_active", "HC.2_inactive",
                                            "HC.3_active", "HC.3_inactive",
                                            "HC.4_active", "HC.4_inactive")

collected_active_inactive_df

NAC.1_active_fraction <- motif_fraction(NAC.1_input, grs, active_motifs)
NAC.2_active_fraction <- motif_fraction(NAC.2_input, grs, active_motifs)
NAC.3_active_fraction <- motif_fraction(NAC.3_input, grs, active_motifs)
NAC.4_active_fraction <- motif_fraction(NAC.4_input, grs, active_motifs)
gc()
NSC.1_active_fraction <- motif_fraction(NSC.1_input, grs, active_motifs)
NSC.2_active_fraction <- motif_fraction(NSC.2_input, grs, active_motifs)
NSC.3_active_fraction <- motif_fraction(NSC.3_input, grs, active_motifs)
NSC.4_active_fraction <- motif_fraction(NSC.4_input, grs, active_motifs)
gc()
SSC.1_active_fraction <- motif_fraction(SSC.1_input, grs, active_motifs)
SSC.2_active_fraction <- motif_fraction(SSC.2_input, grs, active_motifs)
SSC.3_active_fraction <- motif_fraction(SSC.3_input, grs, active_motifs)
SSC.4_active_fraction <- motif_fraction(SSC.4_input, grs, active_motifs)
setwd("D:/Lung cancer cfChIP/Deduped")
Adeno1 <- gene_count("A-1279-cfChIP.bam",grs)
Adeno2 <- gene_count("B-1288-cfChIP.bam", grs)
Adeno3 <- gene_count("C-1475-cfChIP.bam", grs)
Adeno4 <- gene_count("D-1578-cfChIP.bam", grs)
Plano1 <- gene_count("E-439-cfChIP.bam", grs)
Plano2 <- gene_count("F-1449-cfChIP.bam", grs)
Plano3 <- gene_count("I-645-cfChIP.bam",grs)
Plano4 <- gene_count("J-1663-cfChIP.bam",grs)
SCLC1 <- gene_count("G-514-cfChIP.bam", grs)
SCLC2 <- gene_count("H-1169-cfChIP.bam", grs)
SCLC3 <- gene_count("K-440-cfChIP.bam", grs)
SCLC4 <- gene_count("L-1100-cfChIP.bam", grs)
setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
NAC.1_enrichment <- e.score(Adeno1, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
NAC.2_enrichment <- e.score(Adeno2, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
NAC.3_enrichment <- e.score(Adeno3, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
NAC.4_enrichment <- e.score(Adeno4, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
NSC.1_enrichment <- e.score(Plano1, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
NSC.2_enrichment <- e.score(Plano2, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
NSC.3_enrichment <- e.score(Plano3, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
NSC.4_enrichment <- e.score(Plano4, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
SSC.1_enrichment <- e.score(SCLC1, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
SSC.2_enrichment <- e.score(SCLC2, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
SSC.3_enrichment <- e.score(SCLC3, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
SSC.4_enrichment <- e.score(SCLC4, "Coverage of AVENIO genes.txt",gr("AVENIO_genes.txt"))
enrichment_df <- NAC.1_enrichment %>%
    left_join(NAC.2_enrichment, by = "genes") %>%
    left_join(NAC.3_enrichment, by = "genes") %>%
    left_join(NAC.4_enrichment, by = "genes") %>%
    left_join(NSC.1_enrichment, by = "genes") %>%
    left_join(NSC.2_enrichment, by = "genes") %>%
    left_join(NSC.3_enrichment, by = "genes") %>%
    left_join(NSC.4_enrichment, by = "genes") %>%
    left_join(SSC.1_enrichment, by = "genes") %>%
    left_join(SSC.2_enrichment, by = "genes") %>%
    left_join(SSC.3_enrichment, by = "genes") %>%
    left_join(SSC.4_enrichment, by = "genes")
colnames(enrichment_df) <- c("genes","NAC.1",
                             "NAC.2",
                             "NAC.3",
                             "NAC.4",
                             "NSC.1",
                             "NSC.2",
                             "NSC.3",
                             "NSC.4",
                             "SSC.1",
                             "SSC.2",
                             "SSC.3",
                             "SSC.4")

enrichment_df
active_motif_fraction_df <- NAC.1_active_fraction %>%
    left_join(NAC.2_active_fraction, by = "genes") %>%
    left_join(NAC.3_active_fraction, by = "genes") %>%
    left_join(NAC.4_active_fraction, by = "genes") %>%
    left_join(NSC.1_active_fraction, by = "genes") %>%
    left_join(NSC.2_active_fraction, by = "genes") %>%
    left_join(NSC.3_active_fraction, by = "genes") %>%
    left_join(NSC.4_active_fraction, by = "genes") %>%
    left_join(SSC.1_active_fraction, by = "genes") %>%
    left_join(SSC.2_active_fraction, by = "genes") %>%
    left_join(SSC.3_active_fraction, by = "genes") %>%
    left_join(SSC.4_active_fraction, by = "genes")
colnames(active_motif_fraction_df) <- c("genes",
                                        "NAC.1",
                                        "NAC.2",
                                        "NAC.3",
                                        "NAC.4",
                                        "NSC.1",
                                        "NSC.2",
                                        "NSC.3",
                                        "NSC.4",
                                        "SSC.1",
                                        "SSC.2",
                                        "SSC.3",
                                        "SSC.4")
setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver")

write.table(active_motif_fraction_df, file = "active motif fraction.txt", 
            sep = "\t",col.names = T)

NAC.1_short_long_fractions <- short_long_fragments_motifs(NAC.1_input, y = 150, m = 3)
gc()
NAC.2_short_long_fractions <- short_long_fragments_motifs(NAC.2_input, y = 150, m = 3)
gc()
NAC.3_short_long_fractions <- short_long_fragments_motifs(NAC.3_input, y = 150, m = 3)
gc()
NAC.4_short_long_fractions <- short_long_fragments_motifs(NAC.4_input, y = 150, m = 3)
gc()
NSC.1_short_long_fractions <- short_long_fragments_motifs(NSC.1_input, y = 150, m = 3)
gc()
NSC.2_short_long_fractions <- short_long_fragments_motifs(NSC.2_input, y = 150, m = 3)
gc()
NSC.3_short_long_fractions <- short_long_fragments_motifs(NSC.3_input, y = 150, m = 3)
gc()
NSC.4_short_long_fractions <- short_long_fragments_motifs(NSC.4_input, y = 150, m = 3)
gc()
SSC.1_short_long_fractions <- short_long_fragments_motifs(SSC.1_input, y = 150, m = 3)
gc()
SSC.2_short_long_fractions <- short_long_fragments_motifs(SSC.2_input, y = 150, m = 3)
gc()
SSC.3_short_long_fractions <- short_long_fragments_motifs(SSC.3_input, y = 150, m = 3)
gc()
SSC.4_short_long_fractions <- short_long_fragments_motifs(SSC.4_input, y = 150, m = 3)
gc()

short_long_fraction_df <- NAC.1_short_long_fractions %>% 
    left_join(NAC.2_short_long_fractions, by = "motif") %>%
    left_join(NAC.3_short_long_fractions, by = "motif") %>%
    left_join(NAC.4_short_long_fractions, by = "motif") %>%
    left_join(NSC.1_short_long_fractions, by = "motif") %>% 
    left_join(NSC.2_short_long_fractions, by = "motif") %>%
    left_join(NSC.3_short_long_fractions, by = "motif") %>%
    left_join(NSC.4_short_long_fractions, by = "motif") %>%
    left_join(SSC.1_short_long_fractions, by = "motif") %>% 
    left_join(SSC.2_short_long_fractions, by = "motif") %>%
    left_join(SSC.3_short_long_fractions, by = "motif") %>%
    left_join(SSC.4_short_long_fractions, by = "motif")
colnames(short_long_fraction_df) <- c("motif",
                                      "NAC.1_long",
                                      "NAC.1_short",
                                      "NAC.2_long",
                                      "NAC.2_short",
                                      "NAC.3_long",
                                      "NAC.3_short",
                                      "NAC.4_long",
                                      "NAC.4_short",
                                      "NSC.1_long",
                                      "NSC.1_short",
                                      "NSC.2_long",
                                      "NSC.2_short",
                                      "NSC.3_long",
                                      "NSC.3_short",
                                      "NSC.4_long",
                                      "NSC.4_short",
                                      "SSC.1_long",
                                      "SSC.1_short",
                                      "SSC.2_long",
                                      "SSC.2_short",
                                      "SSC.3_long",
                                      "SSC.3_short",
                                      "SSC.4_long",
                                      "SSC.4_short")
short_long_fraction_df
getwd()
write.table(short_long_fraction_df, file = "Fragment end motif fractions of short and long fragments.txt", 
            sep = "\t",col.names = T)

NAC.1_active_short_long <- 
    active_short_long_fragments_motifs(enrichment_df, "NAC.1",
                                       NAC.1_input, grs, y = 150, m = 3)
NAC.2_active_short_long <- 
    active_short_long_fragments_motifs(enrichment_df, "NAC.2",
                                       NAC.2_input, grs, y = 150, m = 3)
NAC.3_active_short_long <- 
    active_short_long_fragments_motifs(enrichment_df, "NAC.3",
                                       NAC.3_input, grs, y = 150, m = 3)
NAC.4_active_short_long <- 
    active_short_long_fragments_motifs(enrichment_df, "NAC.4",
                                       NAC.4_input, grs, y = 150, m = 3)
NSC.1_active_short_long <- 
    active_short_long_fragments_motifs(enrichment_df, "NSC.1",
                                       NSC.1_input, grs, y = 150, m = 3)
NSC.2_active_short_long <- 
    active_short_long_fragments_motifs(enrichment_df, "NSC.2",
                                       NSC.2_input, grs, y = 150, m = 3)
NSC.3_active_short_long <- 
    active_short_long_fragments_motifs(enrichment_df, "NSC.3",
                                       NSC.3_input, grs, y = 150, m = 3)
NSC.4_active_short_long <- 
    active_short_long_fragments_motifs(enrichment_df, "NSC.4",
                                       NSC.4_input, grs, y = 150, m = 3)
SSC.1_active_short_long <- 
    active_short_long_fragments_motifs(enrichment_df, "SSC.1",
                                       SSC.1_input, grs, y = 150, m = 3)
SSC.2_active_short_long <- 
    active_short_long_fragments_motifs(enrichment_df, "SSC.2",
                                       SSC.2_input, grs, y = 150, m = 3)
SSC.3_active_short_long <- 
    active_short_long_fragments_motifs(enrichment_df, "SSC.3",
                                       SSC.3_input, grs, y = 150, m = 3)
SSC.4_active_short_long <- 
    active_short_long_fragments_motifs(enrichment_df, "SSC.4",
                                       SSC.4_input, grs, y = 150, m = 3)

active_long_short_df <- NAC.1_active_short_long %>% 
    left_join(NAC.2_active_short_long, by = "motif") %>% 
    left_join(NAC.3_active_short_long, by = "motif") %>% 
    left_join(NAC.4_active_short_long, by = "motif") %>% 
    left_join(NSC.1_active_short_long, by = "motif") %>% 
    left_join(NSC.2_active_short_long, by = "motif") %>% 
    left_join(NSC.3_active_short_long, by = "motif") %>% 
    left_join(NSC.4_active_short_long, by = "motif") %>% 
    left_join(SSC.1_active_short_long, by = "motif") %>% 
    left_join(SSC.2_active_short_long, by = "motif") %>% 
    left_join(SSC.3_active_short_long, by = "motif") %>% 
    left_join(SSC.4_active_short_long, by = "motif")

colnames(active_long_short_df) <- c("motif",
                                    "NAC.1_bottom_long","NAC.1_bottom_short","NAC.1_top_long","NAC.1_top_short",
                                    "NAC.2_bottom_long","NAC.2_bottom_short","NAC.2_top_long","NAC.2_top_short",
                                    "NAC.3_bottom_long","NAC.3_bottom_short","NAC.3_top_long","NAC.3_top_short",
                                    "NAC.4_bottom_long","NAC.4_bottom_short","NAC.4_top_long","NAC.4_top_short",
                                    "NSC.1_bottom_long","NSC.1_bottom_short","NSC.1_top_long","NSC.1_top_short",
                                    "NSC.2_bottom_long","NSC.2_bottom_short","NSC.2_top_long","NSC.2_top_short",
                                    "NSC.3_bottom_long","NSC.3_bottom_short","NSC.3_top_long","NSC.3_top_short",
                                    "NSC.4_bottom_long","NSC.4_bottom_short","NSC.4_top_long","NSC.4_top_short",
                                    "SSC.1_bottom_long","SSC.1_bottom_short","SSC.1_top_long","SSC.1_top_short",
                                    "SSC.2_bottom_long","SSC.2_bottom_short","SSC.2_top_long","SSC.2_top_short",
                                    "SSC.3_bottom_long","SSC.3_bottom_short","SSC.3_top_long","SSC.3_top_short",
                                    "SSC.4_bottom_long","SSC.4_bottom_short","SSC.4_top_long","SSC.4_top_short")
write.table(active_long_short_df, file = "Fragment end motif fractions of short and long fragments in active and inactive genes.txt", 
            sep = "\t",col.names = T) 


NAC.1_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                            grs,
                                                                            "NAC.1",
                                                                            "D:/Lung cancer input/PosDeduped/PosDeduped-A-1279-input.bam")
NAC.2_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                            grs,
                                                                            "NAC.2",
                                                                            "D:/Lung cancer input/PosDeduped/PosDeduped-B-1288-input.bam")
NAC.3_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                            grs,
                                                                            "NAC.3",
                                                                            "D:/Lung cancer input/PosDeduped/PosDeduped-C-1475-input.bam")
NAC.4_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                            grs,
                                                                            "NAC.4",
                                                                            "D:/Lung cancer input/PosDeduped/PosDeduped-D-1578-input.bam")
NSC.1_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                            grs,
                                                                            "NSC.1",
                                                                            "D:/Lung cancer input/PosDeduped/PosDeduped-E-439-input.bam")
NSC.2_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                            grs,
                                                                            "NSC.2",
                                                                            "D:/Lung cancer input/PosDeduped/PosDeduped-F-1449-input.bam")
NSC.3_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                            grs,
                                                                            "NSC.3",
                                                                            "D:/Lung cancer input/PosDeduped/PosDeduped-I-645-input.bam")
NSC.4_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                            grs,
                                                                            "NSC.4",
                                                                            "D:/Lung cancer input/PosDeduped/PosDeduped-J-1663-input.bam")
SSC.1_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                            grs,
                                                                            "SSC.1",
                                                                            "D:/Lung cancer input/PosDeduped/PosDeduped-G-514-input.bam")
SSC.2_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                            grs,
                                                                            "SSC.2",
                                                                            "D:/Lung cancer input/PosDeduped/PosDeduped-H-1169-input.bam")
SSC.3_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                            grs,
                                                                            "SSC.3",
                                                                            "D:/Lung cancer input/PosDeduped/PosDeduped-K-440-input.bam")
SSC.4_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                            grs,
                                                                            "SSC.4",
                                                                            "D:/Lung cancer input/PosDeduped/PosDeduped-L-1100-input.bam")
HC.1_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                            grs,
                                                                            "HC.1",
                                                                            "D:/Healthy input/PosDeduped/PosDeduped-Rask_kontrol_1_input.bam")

HC.2_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                           grs,
                                                                           "HC.2",
                                                                           "D:/Healthy input/PosDeduped/PosDeduped-Rask_kontrol_2_input.bam")

HC.3_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                           grs,
                                                                           "HC.3",
                                                                           "D:/Healthy input/PosDeduped/PosDeduped-Rask_kontrol_3_input.bam")

HC.4_fragment_length_active_inactive <- fragment_length_active_inactive_df(enrichment_df,
                                                                           grs,
                                                                           "HC.4",
                                                                           "D:/Healthy input/PosDeduped/PosDeduped-Rask_kontrol_4_input.bam")

collected_fragment_length_active_inactive_healthy <- rbind(HC.1_fragment_length_active_inactive,
                                                           HC.2_fragment_length_active_inactive,
                                                           HC.3_fragment_length_active_inactive,
                                                           HC.4_fragment_length_active_inactive)

collected_fragment_length_active_inactive <- rbind(NAC.1_fragment_length_active_inactive,
      NAC.2_fragment_length_active_inactive,
      NAC.3_fragment_length_active_inactive,
      NAC.4_fragment_length_active_inactive,
      NSC.1_fragment_length_active_inactive,
      NSC.2_fragment_length_active_inactive,
      NSC.3_fragment_length_active_inactive,
      NSC.4_fragment_length_active_inactive,
      SSC.1_fragment_length_active_inactive,
      SSC.2_fragment_length_active_inactive,
      SSC.3_fragment_length_active_inactive,
      SSC.4_fragment_length_active_inactive)

collected_fragment_length_active_inactive_healthy <- rbind(HC.1_fragment_length_active_inactive,
                                                           HC.2_fragment_length_active_inactive,
                                                           HC.3_fragment_length_active_inactive,
                                                           HC.4_fragment_length_active_inactive)

write.table(collected_fragment_length_active_inactive_healthy, file = "Fragment lengths active inactive healthy.txt", 
            sep = "\t",col.names = T)

NAC.1_active_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.1","D:/Lung cancer input/PosDeduped/PosDeduped-A-1279-input.bam",3,active_motifs_cancer)
NAC.2_active_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.2","D:/Lung cancer input/PosDeduped/PosDeduped-B-1288-input.bam",3,active_motifs_cancer)
NAC.3_active_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.3","D:/Lung cancer input/PosDeduped/PosDeduped-C-1475-input.bam",3,active_motifs_cancer)
NAC.4_active_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.4","D:/Lung cancer input/PosDeduped/PosDeduped-D-1578-input.bam",3,active_motifs_cancer)
NSC.1_active_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.1","D:/Lung cancer input/PosDeduped/PosDeduped-E-439-input.bam",3,active_motifs_cancer)
NSC.2_active_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.2","D:/Lung cancer input/PosDeduped/PosDeduped-F-1449-input.bam",3,active_motifs_cancer)
NSC.3_active_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.3","D:/Lung cancer input/PosDeduped/PosDeduped-I-645-input.bam",3,active_motifs_cancer)
NSC.4_active_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.4","D:/Lung cancer input/PosDeduped/PosDeduped-J-1663-input.bam",3,active_motifs_cancer)
SSC.1_active_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.1","D:/Lung cancer input/PosDeduped/PosDeduped-G-514-input.bam",3,active_motifs_cancer)
SSC.2_active_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.2","D:/Lung cancer input/PosDeduped/PosDeduped-H-1169-input.bam",3,active_motifs_cancer)
SSC.3_active_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.3","D:/Lung cancer input/PosDeduped/PosDeduped-K-440-input.bam",3,active_motifs_cancer)
SSC.4_active_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.4","D:/Lung cancer input/PosDeduped/PosDeduped-L-1100-input.bam",3,active_motifs_cancer)

NAC.1_inactive_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.1","D:/Lung cancer input/PosDeduped/PosDeduped-A-1279-input.bam",3,inactive_motifs_cancer)
NAC.2_inactive_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.2","D:/Lung cancer input/PosDeduped/PosDeduped-B-1288-input.bam",3,inactive_motifs_cancer)
NAC.3_inactive_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.3","D:/Lung cancer input/PosDeduped/PosDeduped-C-1475-input.bam",3,inactive_motifs_cancer)
NAC.4_inactive_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.4","D:/Lung cancer input/PosDeduped/PosDeduped-D-1578-input.bam",3,inactive_motifs_cancer)
NSC.1_inactive_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.1","D:/Lung cancer input/PosDeduped/PosDeduped-E-439-input.bam",3,inactive_motifs_cancer)
NSC.2_inactive_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.2","D:/Lung cancer input/PosDeduped/PosDeduped-F-1449-input.bam",3,inactive_motifs_cancer)
NSC.3_inactive_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.3","D:/Lung cancer input/PosDeduped/PosDeduped-I-645-input.bam",3,inactive_motifs_cancer)
NSC.4_inactive_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.4","D:/Lung cancer input/PosDeduped/PosDeduped-J-1663-input.bam",3,inactive_motifs_cancer)
SSC.1_inactive_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.1","D:/Lung cancer input/PosDeduped/PosDeduped-G-514-input.bam",3,inactive_motifs_cancer)
SSC.2_inactive_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.2","D:/Lung cancer input/PosDeduped/PosDeduped-H-1169-input.bam",3,inactive_motifs_cancer)
SSC.3_inactive_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.3","D:/Lung cancer input/PosDeduped/PosDeduped-K-440-input.bam",3,inactive_motifs_cancer)
SSC.4_inactive_motif_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.4","D:/Lung cancer input/PosDeduped/PosDeduped-L-1100-input.bam",3,inactive_motifs_cancer)

df <- rbind(NAC.1_inactive_motif_quartiles,
            NAC.2_inactive_motif_quartiles,
            NAC.3_inactive_motif_quartiles,
            NAC.4_inactive_motif_quartiles,
            NSC.1_inactive_motif_quartiles,
            NSC.2_inactive_motif_quartiles,
            NSC.3_inactive_motif_quartiles,
            NSC.4_inactive_motif_quartiles,
            SSC.1_inactive_motif_quartiles,
            SSC.2_inactive_motif_quartiles,
            SSC.3_inactive_motif_quartiles,
            SSC.4_inactive_motif_quartiles)

write.table(df, file = "Inactive motif cancer fraction in cfChIP quatiles.txt", 
            sep = "\t",col.names = T)

df <- rbind(NAC.1_active_motif_quartiles,
            NAC.2_active_motif_quartiles,
            NAC.3_active_motif_quartiles,
            NAC.4_active_motif_quartiles,
            NSC.1_active_motif_quartiles,
            NSC.2_active_motif_quartiles,
            NSC.3_active_motif_quartiles,
            NSC.4_active_motif_quartiles,
            SSC.1_active_motif_quartiles,
            SSC.2_active_motif_quartiles,
            SSC.3_active_motif_quartiles,
            SSC.4_active_motif_quartiles)

write.table(df, file = "Active motif cancer fraction in cfChIP quatiles.txt", 
            sep = "\t",col.names = T)

NAC.1_active_active_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.1","D:/Lung cancer input/PosDeduped/PosDeduped-A-1279-input.bam",3,active_active)
NAC.2_active_active_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.2","D:/Lung cancer input/PosDeduped/PosDeduped-B-1288-input.bam",3,active_active)
NAC.3_active_active_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.3","D:/Lung cancer input/PosDeduped/PosDeduped-C-1475-input.bam",3,active_active)
NAC.4_active_active_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.4","D:/Lung cancer input/PosDeduped/PosDeduped-D-1578-input.bam",3,active_active)
NSC.1_active_active_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.1","D:/Lung cancer input/PosDeduped/PosDeduped-E-439-input.bam",3,active_active)
NSC.2_active_active_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.2","D:/Lung cancer input/PosDeduped/PosDeduped-F-1449-input.bam",3,active_active)
NSC.3_active_active_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.3","D:/Lung cancer input/PosDeduped/PosDeduped-I-645-input.bam",3,active_active)
NSC.4_active_active_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.4","D:/Lung cancer input/PosDeduped/PosDeduped-J-1663-input.bam",3,active_active)
SSC.1_active_active_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.1","D:/Lung cancer input/PosDeduped/PosDeduped-G-514-input.bam",3,active_active)
SSC.2_active_active_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.2","D:/Lung cancer input/PosDeduped/PosDeduped-H-1169-input.bam",3,active_active)
SSC.3_active_active_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.3","D:/Lung cancer input/PosDeduped/PosDeduped-K-440-input.bam",3,active_active)
SSC.4_active_active_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.4","D:/Lung cancer input/PosDeduped/PosDeduped-L-1100-input.bam",3,active_active)

df <- rbind(NAC.1_active_active_quartiles,
            NAC.2_active_active_quartiles,
            NAC.3_active_active_quartiles,
            NAC.4_active_active_quartiles,
            NSC.1_active_active_quartiles,
            NSC.2_active_active_quartiles,
            NSC.3_active_active_quartiles,
            NSC.4_active_active_quartiles,
            SSC.1_active_active_quartiles,
            SSC.2_active_active_quartiles,
            SSC.3_active_active_quartiles,
            SSC.4_active_active_quartiles)

write.table(df, file = "Active active fraction in cfChIP quatiles.txt", 
            sep = "\t",col.names = T)

NAC.1_inactive_inactive_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.1","D:/Lung cancer input/PosDeduped/PosDeduped-A-1279-input.bam",3,inactive_inactive)
NAC.2_inactive_inactive_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.2","D:/Lung cancer input/PosDeduped/PosDeduped-B-1288-input.bam",3,inactive_inactive)
NAC.3_inactive_inactive_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.3","D:/Lung cancer input/PosDeduped/PosDeduped-C-1475-input.bam",3,inactive_inactive)
NAC.4_inactive_inactive_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NAC.4","D:/Lung cancer input/PosDeduped/PosDeduped-D-1578-input.bam",3,inactive_inactive)
NSC.1_inactive_inactive_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.1","D:/Lung cancer input/PosDeduped/PosDeduped-E-439-input.bam",3,inactive_inactive)
NSC.2_inactive_inactive_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.2","D:/Lung cancer input/PosDeduped/PosDeduped-F-1449-input.bam",3,inactive_inactive)
NSC.3_inactive_inactive_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.3","D:/Lung cancer input/PosDeduped/PosDeduped-I-645-input.bam",3,inactive_inactive)
NSC.4_inactive_inactive_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"NSC.4","D:/Lung cancer input/PosDeduped/PosDeduped-J-1663-input.bam",3,inactive_inactive)
SSC.1_inactive_inactive_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.1","D:/Lung cancer input/PosDeduped/PosDeduped-G-514-input.bam",3,inactive_inactive)
SSC.2_inactive_inactive_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.2","D:/Lung cancer input/PosDeduped/PosDeduped-H-1169-input.bam",3,inactive_inactive)
SSC.3_inactive_inactive_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.3","D:/Lung cancer input/PosDeduped/PosDeduped-K-440-input.bam",3,inactive_inactive)
SSC.4_inactive_inactive_quartiles <- active_motif_quartiles_df(enrichment_df,grs,"SSC.4","D:/Lung cancer input/PosDeduped/PosDeduped-L-1100-input.bam",3,inactive_inactive)

df <- rbind(NAC.1_inactive_inactive_quartiles,
            NAC.2_inactive_inactive_quartiles,
            NAC.3_inactive_inactive_quartiles,
            NAC.4_inactive_inactive_quartiles,
            NSC.1_inactive_inactive_quartiles,
            NSC.2_inactive_inactive_quartiles,
            NSC.3_inactive_inactive_quartiles,
            NSC.4_inactive_inactive_quartiles,
            SSC.1_inactive_inactive_quartiles,
            SSC.2_inactive_inactive_quartiles,
            SSC.3_inactive_inactive_quartiles,
            SSC.4_inactive_inactive_quartiles)

write.table(df, file = "Inactive inactive fraction in cfChIP quatiles.txt", 
            sep = "\t",col.names = T)

NAC.1_sub150_active_motif_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NAC.1","D:/Lung cancer input/PosDeduped/PosDeduped-A-1279-input.bam",150,3,active_motifs_cancer)
NAC.2_sub150_active_motif_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NAC.2","D:/Lung cancer input/PosDeduped/PosDeduped-B-1288-input.bam",150,3,active_motifs_cancer)
NAC.3_sub150_active_motif_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NAC.3","D:/Lung cancer input/PosDeduped/PosDeduped-C-1475-input.bam",150,3,active_motifs_cancer)
NAC.4_sub150_active_motif_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NAC.4","D:/Lung cancer input/PosDeduped/PosDeduped-D-1578-input.bam",150,3,active_motifs_cancer)
NSC.1_sub150_active_motif_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NSC.1","D:/Lung cancer input/PosDeduped/PosDeduped-E-439-input.bam",150,3,active_motifs_cancer)
NSC.2_sub150_active_motif_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NSC.2","D:/Lung cancer input/PosDeduped/PosDeduped-F-1449-input.bam",150,3,active_motifs_cancer)
NSC.3_sub150_active_motif_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NSC.3","D:/Lung cancer input/PosDeduped/PosDeduped-I-645-input.bam",150,3,active_motifs_cancer)
NSC.4_sub150_active_motif_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NSC.4","D:/Lung cancer input/PosDeduped/PosDeduped-J-1663-input.bam",150,3,active_motifs_cancer)
SSC.1_sub150_active_motif_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"SSC.1","D:/Lung cancer input/PosDeduped/PosDeduped-G-514-input.bam",150,3,active_motifs_cancer)
SSC.2_sub150_active_motif_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"SSC.2","D:/Lung cancer input/PosDeduped/PosDeduped-H-1169-input.bam",150,3,active_motifs_cancer)
SSC.3_sub150_active_motif_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"SSC.3","D:/Lung cancer input/PosDeduped/PosDeduped-K-440-input.bam",150,3,active_motifs_cancer)
SSC.4_sub150_active_motif_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"SSC.4","D:/Lung cancer input/PosDeduped/PosDeduped-L-1100-input.bam",150,3,active_motifs_cancer)

df <- rbind(NAC.1_sub150_active_motif_quartiles,
            NAC.2_sub150_active_motif_quartiles,
            NAC.3_sub150_active_motif_quartiles,
            NAC.4_sub150_active_motif_quartiles,
            NSC.1_sub150_active_motif_quartiles,
            NSC.2_sub150_active_motif_quartiles,
            NSC.3_sub150_active_motif_quartiles,
            NSC.4_sub150_active_motif_quartiles,
            SSC.1_sub150_active_motif_quartiles,
            SSC.2_sub150_active_motif_quartiles,
            SSC.3_sub150_active_motif_quartiles,
            SSC.4_sub150_active_motif_quartiles)

write.table(df, file = "Sub 150 with active motif cancer fraction in cfChIP quatiles.txt", 
            sep = "\t",col.names = T)

NAC.1_sub150_active_active_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NAC.1","D:/Lung cancer input/PosDeduped/PosDeduped-A-1279-input.bam",150,3,active_active)
NAC.2_sub150_active_active_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NAC.2","D:/Lung cancer input/PosDeduped/PosDeduped-B-1288-input.bam",150,3,active_active)
NAC.3_sub150_active_active_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NAC.3","D:/Lung cancer input/PosDeduped/PosDeduped-C-1475-input.bam",150,3,active_active)
NAC.4_sub150_active_active_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NAC.4","D:/Lung cancer input/PosDeduped/PosDeduped-D-1578-input.bam",150,3,active_active)
NSC.1_sub150_active_active_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NSC.1","D:/Lung cancer input/PosDeduped/PosDeduped-E-439-input.bam",150,3,active_active)
NSC.2_sub150_active_active_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NSC.2","D:/Lung cancer input/PosDeduped/PosDeduped-F-1449-input.bam",150,3,active_active)
NSC.3_sub150_active_active_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NSC.3","D:/Lung cancer input/PosDeduped/PosDeduped-I-645-input.bam",150,3,active_active)
NSC.4_sub150_active_active_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"NSC.4","D:/Lung cancer input/PosDeduped/PosDeduped-J-1663-input.bam",150,3,active_active)
SSC.1_sub150_active_active_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"SSC.1","D:/Lung cancer input/PosDeduped/PosDeduped-G-514-input.bam",150,3,active_active)
SSC.2_sub150_active_active_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"SSC.2","D:/Lung cancer input/PosDeduped/PosDeduped-H-1169-input.bam",150,3,active_active)
SSC.3_sub150_active_active_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"SSC.3","D:/Lung cancer input/PosDeduped/PosDeduped-K-440-input.bam",150,3,active_active)
SSC.4_sub150_active_active_quartiles <- sub150_active_motif_fraction_quartiles_df(enrichment_df,grs,"SSC.4","D:/Lung cancer input/PosDeduped/PosDeduped-L-1100-input.bam",150,3,active_active)

df <- rbind(NAC.1_sub150_active_active_quartiles,
            NAC.2_sub150_active_active_quartiles,
            NAC.3_sub150_active_active_quartiles,
            NAC.4_sub150_active_active_quartiles,
            NSC.1_sub150_active_active_quartiles,
            NSC.2_sub150_active_active_quartiles,
            NSC.3_sub150_active_active_quartiles,
            NSC.4_sub150_active_active_quartiles,
            SSC.1_sub150_active_active_quartiles,
            SSC.2_sub150_active_active_quartiles,
            SSC.3_sub150_active_active_quartiles,
            SSC.4_sub150_active_active_quartiles)

write.table(df, file = "Sub 150 with active active fraction in cfChIP quatiles.txt", 
            sep = "\t",col.names = T)
