library(ggplot2)
th <- theme(
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 10, vjust = 1, hjust = 0.37),
    axis.text.x = element_text(angle = 0, size = 12),
    axis.text.y = element_text(angle = 0, size = 12),
    axis.title = element_text(size = 14, face = 'bold'),
    axis.title.x = element_text(size = 14, face = 'bold'),
    axis.title.y = element_text(size = 14, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    legend.text = element_text(size = 12), 
    title = element_text(size = 14, face = "bold"))
setwd("C:/Users/Christoffer/OneDrive/1PhD/R files/Important excel and text documents")
grs <- gr("AVENIO_genes.txt")
bad <- badgene("Complete table of all targets.txt",25)
multiple_region_genes <- multiple_region("Complete table of all targets.txt")
setwd("D:/Lung cancer cfChIP/PosDeduped")
NAC.1_cfChIP <- bamfile("PosDeduped-A-1279-cfChIP.bam",
                        "PosDeduped-A-1279-cfChIP.bam.bai")
NAC.2_cfChIP <- bamfile("PosDeduped-B-1288-cfChIP.bam",
                        "PosDeduped-B-1288-cfChIP.bam.bai")
NAC.3_cfChIP <- bamfile("PosDeduped-C-1475-cfChIP.bam",
                        "PosDeduped-C-1475-cfChIP.bam.bai")
NAC.4_cfChIP <- bamfile("PosDeduped-D-1578-cfChIP.bam",
                        "PosDeduped-D-1578-cfChIP.bam.bai")
NSC.1_cfChIP <- bamfile("PosDeduped-E-439-cfChIP.bam",
                        "PosDeduped-E-439-cfChIP.bam.bai")
NSC.2_cfChIP <- bamfile("PosDeduped-F-1449-cfChIP.bam",
                        "PosDeduped-F-1449-cfChIP.bam.bai")
NSC.3_cfChIP <- bamfile("PosDeduped-I-645-cfChIP.bam",
                        "PosDeduped-I-645-cfChIP.bam.bai")
NSC.4_cfChIP <- bamfile("PosDeduped-J-1663-cfChIP.bam",
                        "PosDeduped-J-1663-cfChIP.bam.bai")
SSC.1_cfChIP <- bamfile("PosDeduped-G-514-cfChIP.bam",
                        "PosDeduped-G-514-cfChIP.bam.bai")
SSC.2_cfChIP <- bamfile("PosDeduped-H-1169-cfChIP.bam",
                        "PosDeduped-H-1169-cfChIP.bam.bai")
SSC.3_cfChIP <- bamfile("PosDeduped-K-440-cfChIP.bam",
                        "PosDeduped-K-440-cfChIP.bam.bai")
SSC.4_cfChIP <- bamfile("PosDeduped-L-1100-cfChIP.bam",
                        "PosDeduped-L-1100-cfChIP.bam.bai")
setwd("D:/Healthy cfChIP/PosDeduped")
HC.1_cfChIP <- bamfile("PosDeduped-Rask_kontrol_1_cfChIP.bam",
                       "PosDeduped-Rask_kontrol_1_cfChIP.bam.bai")
HC.2_cfChIP <- bamfile("PosDeduped-Rask_kontrol_2_cfChIP.bam",
                       "PosDeduped-Rask_kontrol_2_cfChIP.bam.bai")
HC.3_cfChIP <- bamfile("PosDeduped-Rask_kontrol_3_cfChIP.bam",
                       "PosDeduped-Rask_kontrol_3_cfChIP.bam.bai")
HC.4_cfChIP <- bamfile("PosDeduped-Rask_kontrol_4_cfChIP.bam",
                       "PosDeduped-Rask_kontrol_4_cfChIP.bam.bai")
setwd("D:/Lung cancer input/PosDeduped")
NAC.1_input <- bamfile("PosDeduped-A-1279-input.bam",
                       "PosDeduped-A-1279-input.bam.bai")
NAC.2_input <- bamfile("PosDeduped-B-1288-input.bam",
                       "PosDeduped-B-1288-input.bam.bai")
NAC.3_input <- bamfile("PosDeduped-C-1475-input.bam",
                       "PosDeduped-C-1475-input.bam.bai")
NAC.4_input <- bamfile("PosDeduped-D-1578-input.bam",
                       "PosDeduped-D-1578-input.bam.bai")
NSC.1_input <- bamfile("PosDeduped-E-439-input.bam",
                       "PosDeduped-E-439-input.bam.bai")
NSC.2_input <- bamfile("PosDeduped-F-1449-input.bam",
                       "PosDeduped-F-1449-input.bam.bai")
NSC.3_input <- bamfile("PosDeduped-I-645-input.bam",
                       "PosDeduped-I-645-input.bam.bai")
NSC.4_input <- bamfile("PosDeduped-J-1663-input.bam",
                       "PosDeduped-J-1663-input.bam.bai")
SSC.1_input <- bamfile("PosDeduped-G-514-input.bam",
                       "PosDeduped-G-514-input.bam.bai")
SSC.2_input <- bamfile("PosDeduped-H-1169-input.bam",
                       "PosDeduped-H-1169-input.bam.bai")
SSC.3_input <- bamfile("PosDeduped-K-440-input.bam",
                       "PosDeduped-K-440-input.bam.bai")
SSC.4_input <- bamfile("PosDeduped-L-1100-input.bam",
                       "PosDeduped-L-1100-input.bam.bai")

setwd("D:/Healthy input/PosDeduped")
HC.1_input <- bamfile("PosDeduped-Rask_kontrol_1_input.bam",
                      "PosDeduped-Rask_kontrol_1_input.bam.bai")
HC.2_input <- bamfile("PosDeduped-Rask_kontrol_2_input.bam",
                      "PosDeduped-Rask_kontrol_2_input.bam.bai")
HC.3_input <- bamfile("PosDeduped-Rask_kontrol_3_input.bam",
                      "PosDeduped-Rask_kontrol_3_input.bam.bai")
HC.4_input <- bamfile("PosDeduped-Rask_kontrol_4_input.bam",
                      "PosDeduped-Rask_kontrol_4_input.bam.bai")

library(dplyr)
setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver")
prop_mer_input <- read.table("fragment end motif proportions input.txt", header = T) %>% select(-X)
prop_mer_cfChIP <- read.table("fragment end motif proportions cfChIP.txt", header = T) %>% select(-X)
count_mer_input <- read.table("fragment end motif counts input.txt", header = T) %>% select(-X)
count_mer_cfChIP <- read.table("fragment end motif counts cfChIP.txt", header = T) %>% select(-X)

prop_mer_input_cancer <- prop_mer_input[!grepl("HC.",colnames(prop_mer_input))]
prop_mer_cfChIP_cancer <- prop_mer_cfChIP[!grepl("HC.",colnames(prop_mer_cfChIP))]
prop_mer_input_healthy <- prop_mer_input[grepl("HC.",colnames(prop_mer_input))]
prop_mer_cfChIP_healthy <- prop_mer_cfChIP[grepl("HC.",colnames(prop_mer_cfChIP))]
prop_mer_input_healthy$motif <- prop_mer_input$motif
prop_mer_cfChIP_healthy$motif <- prop_mer_cfChIP$motif

umap_moitf(prop_mer_input, prop_mer_cfChIP, "input", "cfChIP", 15, 15)

setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver")
collected_active_inactive_df <- read.table("fragment end motif proportion of active and inactive fragments.txt", header = T)

volcano_motif("Input cancer", "cfChIP cancer", dif_motif_prop(prop_mer_input_cancer, prop_mer_cfChIP_cancer))
volcano_motif("Input healthy", "cfChIP healthy", dif_motif_prop(prop_mer_input_healthy, prop_mer_cfChIP_healthy))
volcano_motif("Input", "cfChIP", dif_motif_prop(prop_mer_input, prop_mer_cfChIP))
volcano_motif("Inactive genes", "Active genes", diff_active_inactive_end_motif(collected_active_inactive_df))

sequence_list <- list(sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "positive"),
                      sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "negative"),
                      sequence_end_motif(diff_active_inactive_end_motif(collected_active_inactive_df), 0.05, "positive"),
                      sequence_end_motif(diff_active_inactive_end_motif(collected_active_inactive_df), 0.05, "negative"))


sequence_end_motif_plot(sequence_list, c("cfChIP", "Input",
                                         "Active genes", "Inactive genes"))
active_motifs <- sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "positive")
active_motifs
inactive_motifs <- sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "negative")
setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver")
enrichment_df <- read.table("enrichment cancer.txt", header = T)
active_motif_fraction_df <- read.table("active motif fraction.txt", header = T)
active_motif_fraction_df


motif_fraction_vs_enrichment(NAC.1_active_fraction,
                             NAC.1_enrichment, "NAC.1", bad, multiple_region_genes)
motif_fraction_vs_enrichment(NAC.1_active_fraction,
                             NAC.1_enrichment, "NAC.1", bad)
motif_fraction_vs_enrichment(NAC.2_active_fraction,
                             NAC.2_enrichment, "NAC.2",bad, multiple_region_genes)
motif_fraction_vs_enrichment(NAC.2_active_fraction,
                             NAC.2_enrichment, "NAC.2",bad)
motif_fraction_vs_enrichment(NAC.3_active_fraction,
                             NAC.3_enrichment, "NAC.3",bad, multiple_region_genes)
motif_fraction_vs_enrichment(NAC.3_active_fraction,
                             NAC.3_enrichment, "NAC.3",bad)
motif_fraction_vs_enrichment(NAC.4_active_fraction,
                             NAC.4_enrichment, "NAC.4",bad, multiple_region_genes)
motif_fraction_vs_enrichment(NAC.4_active_fraction,
                             NAC.4_enrichment, "NAC.4",bad)
motif_fraction_vs_enrichment(NSC.1_active_fraction,
                             NSC.1_enrichment, "NSC.1", bad, multiple_region_genes)
motif_fraction_vs_enrichment(NSC.1_active_fraction,
                             NSC.1_enrichment, "NSC.1", bad)
motif_fraction_vs_enrichment(NSC.2_active_fraction,
                             NSC.2_enrichment, "NSC.2",bad, multiple_region_genes)
motif_fraction_vs_enrichment(NSC.2_active_fraction,
                             NSC.2_enrichment, "NSC.2",bad)
motif_fraction_vs_enrichment(NSC.3_active_fraction,
                             NSC.3_enrichment, "NSC.3",bad, multiple_region_genes)
motif_fraction_vs_enrichment(NSC.3_active_fraction,
                             NSC.3_enrichment, "NSC.3",bad)
motif_fraction_vs_enrichment(NSC.4_active_fraction,
                             NSC.4_enrichment, "NSC.4",bad, multiple_region_genes)
motif_fraction_vs_enrichment(NSC.4_active_fraction,
                             NSC.4_enrichment, "NSC.4",bad)
motif_fraction_vs_enrichment(SSC.1_active_fraction,
                             SSC.1_enrichment, "SSC.1", bad, multiple_region_genes)
motif_fraction_vs_enrichment(SSC.1_active_fraction,
                             SSC.1_enrichment, "SSC.1", bad)
motif_fraction_vs_enrichment(SSC.2_active_fraction,
                             SSC.2_enrichment, "SSC.2",bad, multiple_region_genes)
motif_fraction_vs_enrichment(SSC.2_active_fraction,
                             SSC.2_enrichment, "SSC.2",bad)
motif_fraction_vs_enrichment(SSC.3_active_fraction,
                             SSC.3_enrichment, "SSC.3",bad, multiple_region_genes)
motif_fraction_vs_enrichment(SSC.3_active_fraction,
                             SSC.3_enrichment, "SSC.3",bad)
motif_fraction_vs_enrichment(SSC.4_active_fraction,
                             SSC.4_enrichment, "SSC.4",bad, multiple_region_genes)
motif_fraction_vs_enrichment(SSC.4_active_fraction,
                             SSC.4_enrichment, "SSC.4",bad)

write.table(enrichment_df, file = "enrichment cancer.txt", 
            sep = "\t",col.names = T)    
getwd()
multiple_regions_correlations <- data.frame(
    sample = colnames(enrichment_df)[2:length(colnames(enrichment_df))],
    rho = c(0.51,0.49,0.31,0.22,0.0,0.05,0.25,0.05,0.21,0.58,0.21,0.18),
    p.value = c(0.004,0.007,0.095,0.241,0.972,0.785,0.178,0.7997,0.28,0.0011,0.256,0.347),
    alpha = c(57.65,49.51,38.26,17.31,19.22,20.77,20.03,-0.13,29.95,54.78,60.62,15.55)
)
multiple_regions_correlations$sample <- factor(multiple_regions_correlations$sample , levels=unique(multiple_regions_correlations$sample))
all_genes_correlations <- data.frame(
    sample = colnames(enrichment_df)[2:length(colnames(enrichment_df))],
    rho = c(0.43,0.58,0.44,0.32,0.2,0.32,0.28,0.15,0.31,0.50,0.37,0.31),
    p.value = c(0.0001,0.0001,0.0001,0.0001,0.011,0.0001,0.0003,0.0627,0.0001,0.0001,0.0001,0.0001),
    alpha = c(51.52,70.82,51.4,38.77,32.4,38.76,31.9,17.42,35.29,54.32,69.13,31.01))
all_genes_correlations
