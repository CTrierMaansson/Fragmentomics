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
prop_mer_input <- read.table("fragment end motif proportions input.txt", header = T)
prop_mer_cfChIP <- read.table("fragment end motif proportions cfChIP.txt", header = T)
count_mer_input <- read.table("fragment end motif counts input.txt", header = T)
count_mer_cfChIP <- read.table("fragment end motif counts cfChIP.txt", header = T)

prop_mer_input_cancer <- prop_mer_input[!grepl("HC.",colnames(prop_mer_input))]
prop_mer_cfChIP_cancer <- prop_mer_cfChIP[!grepl("HC.",colnames(prop_mer_cfChIP))]
prop_mer_input_healthy <- prop_mer_input[grepl("HC.",colnames(prop_mer_input))]
prop_mer_cfChIP_healthy <- prop_mer_cfChIP[grepl("HC.",colnames(prop_mer_cfChIP))]
prop_mer_input_healthy$motif <- prop_mer_input$motif
prop_mer_cfChIP_healthy$motif <- prop_mer_cfChIP$motif

umap_moitf(prop_mer_input, prop_mer_cfChIP, "input", "cfChIP", 15, 15)
umap_moitf(prop_mer_input_cancer, prop_mer_cfChIP_cancer, "Cancer input", "Cancer cfChIP", 15, 15, t = "cancer")
umap_moitf(prop_mer_input_healthy, prop_mer_cfChIP_healthy, "Healthy input", "Healthy cfChIP", 15, 5, t = "healthy")

collected_active_inactive_df <- read.table("fragment end motif proportion of active and inactive fragments.txt", header = T)
cancer_active_inactive <- collected_active_inactive_df[!grepl("HC", colnames(collected_active_inactive_df))]
healthy_active_inactive <- collected_active_inactive_df[grepl("HC", colnames(collected_active_inactive_df))]
healthy_active_inactive$motif <- collected_active_inactive_df$motif
setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver")
volcano_motif("Input healthy", "Input cancer", diff_gene_motiv_analysis(prop_mer_input_healthy, prop_mer_input_cancer))
volcano_motif("cfChIP healthy", "cfChIP cancer", diff_gene_motiv_analysis(prop_mer_cfChIP_healthy, prop_mer_cfChIP_cancer))
volcano_motif("Input cancer", "cfChIP cancer", dif_motif_prop(prop_mer_input_cancer, prop_mer_cfChIP_cancer))
volcano_motif("Input healthy", "cfChIP healthy", dif_motif_prop(prop_mer_input_healthy, prop_mer_cfChIP_healthy))
volcano_motif("Input", "cfChIP", dif_motif_prop(prop_mer_input, prop_mer_cfChIP))
volcano_motif("Cancer inactive genes", "Cancer active genes", diff_active_inactive_end_motif(cancer_active_inactive))
volcano_motif("Healthy inactive genes", "Healthy active genes", diff_active_inactive_end_motif(healthy_active_inactive))
volcano_motif("Inactive genes", "Active genes", diff_active_inactive_end_motif(collected_active_inactive_df))



sequence_list_all_samples <- list(sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "positive"),
                      sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "negative"),
                      sequence_end_motif(diff_active_inactive_end_motif(collected_active_inactive_df), 0.05, "positive"),
                      sequence_end_motif(diff_active_inactive_end_motif(collected_active_inactive_df), 0.05, "negative"))
sequence_list_cancer <- list(sequence_end_motif(dif_motif_prop(prop_mer_input_cancer, prop_mer_cfChIP_cancer), 0.05, "positive"),
                                  sequence_end_motif(dif_motif_prop(prop_mer_input_cancer, prop_mer_cfChIP_cancer), 0.05, "negative"),
                                  sequence_end_motif(diff_active_inactive_end_motif(cancer_active_inactive), 0.05, "positive"),
                                  sequence_end_motif(diff_active_inactive_end_motif(cancer_active_inactive), 0.05, "negative"))
sequence_list_healthy <- list(sequence_end_motif(dif_motif_prop(prop_mer_input_healthy, prop_mer_cfChIP_healthy), 0.05, "positive"),
                             sequence_end_motif(dif_motif_prop(prop_mer_input_healthy, prop_mer_cfChIP_healthy), 0.05, "negative"),
                             sequence_end_motif(diff_active_inactive_end_motif(healthy_active_inactive), 0.05, "positive"),
                             sequence_end_motif(diff_active_inactive_end_motif(healthy_active_inactive), 0.05, "negative"))


sequence_end_motif_plot(sequence_list_all_samples, c("cfChIP", "Input",
                                         "Active genes", "Inactive genes"))
sequence_end_motif_plot(sequence_list_cancer, c("Cancer cfChIP", "Cancer Input",
                                                     "Cancer Active genes", "Cancer Inactive genes"))
sequence_end_motif_plot(sequence_list_healthy, c("Healthy cfChIP", "Healthy Input",
                                                "Healthy Active genes", "Healthy Inactive genes"))
motif_venn(sequence_list_all_samples,
           c("cfChIP", "Input",
             "Active", "Inactive"), 
           "All samples")

motif_venn(sequence_list_cancer,
           c("cfChIP", "Input",
             "Active", "Inactive"), 
           "Cancer samples")
motif_venn(sequence_list_healthy,
           c("cfChIP", "Input",
             "Active", "Inactive"), 
           "Healthy samples")
active_motifs <- sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "positive")
active_motifs
inactive_motifs <- sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "negative")
setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver")
enrichment_df <- read.table("enrichment all.txt", header = T)
active_motif_fraction_df <- read.table("active motif fraction.txt", header = T)
active_motif_fraction_df


motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "SSC.4")],
                             enrichment_df[c("genes", "SSC.4")], "SSC.4", bad, multiple_region_genes)
motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "SSC.4")],
                             enrichment_df[c("genes", "SSC.4")], "SSC.4", bad)

NAC.1_input_sub150_reads <- short_long_fragments(NAC.1_input)
NAC.1_input_above150_reads <- short_long_fragments(NAC.1_input, z = FALSE)
NAC.2_input_sub150_reads <- short_long_fragments(NAC.2_input)
NAC.2_input_above150_reads <- short_long_fragments(NAC.2_input, z = FALSE)
NAC.3_input_sub150_reads <- short_long_fragments(NAC.3_input)
NAC.3_input_above150_reads <- short_long_fragments(NAC.3_input, z = FALSE)
NAC.4_input_sub150_reads <- short_long_fragments(NAC.4_input)
NAC.4_input_above150_reads <- short_long_fragments(NAC.4_input, z = FALSE)
NSC.1_input_sub150_reads <- short_long_fragments(NSC.1_input)
NSC.1_input_above150_reads <- short_long_fragments(NSC.1_input, z = FALSE)
NSC.2_input_sub150_reads <- short_long_fragments(NSC.2_input)
NSC.2_input_above150_reads <- short_long_fragments(NSC.2_input, z = FALSE)
NSC.3_input_sub150_reads <- short_long_fragments(NSC.3_input)
NSC.3_input_above150_reads <- short_long_fragments(NSC.3_input, z = FALSE)
NSC.4_input_sub150_reads <- short_long_fragments(NSC.4_input)
NSC.4_input_above150_reads <- short_long_fragments(NSC.4_input, z = FALSE)
SSC.1_input_sub150_reads <- short_long_fragments(SSC.1_input)
SSC.1_input_above150_reads <- short_long_fragments(SSC.1_input, z = FALSE)
SSC.2_input_sub150_reads <- short_long_fragments(SSC.2_input)
SSC.2_input_above150_reads <- short_long_fragments(SSC.2_input, z = FALSE)
SSC.3_input_sub150_reads <- short_long_fragments(SSC.3_input)
SSC.3_input_above150_reads <- short_long_fragments(SSC.3_input, z = FALSE)
SSC.4_input_sub150_reads <- short_long_fragments(SSC.4_input)
SSC.4_input_above150_reads <- short_long_fragments(SSC.4_input, z = FALSE)

getwd()
multiple_regions_correlations <- data.frame(
    sample = colnames(enrichment_df)[2:(length(colnames(enrichment_df))-4)],
    rho = c(0.51,0.49,0.31,0.22,0.0,0.05,0.25,0.05,0.21,0.58,0.21,0.18),
    p.value = c(0.004,0.007,0.095,0.241,0.972,0.785,0.178,0.7997,0.28,0.0011,0.256,0.347),
    alpha = c(57.65,49.51,38.26,17.31,19.22,20.77,20.03,-0.13,29.95,54.78,60.62,15.55)
)
multiple_regions_correlations$sample <- factor(multiple_regions_correlations$sample , levels=unique(multiple_regions_correlations$sample))
all_genes_correlations <- data.frame(
    sample = colnames(enrichment_df)[2:(length(colnames(enrichment_df))-4)],
    rho = c(0.43,0.58,0.44,0.32,0.2,0.32,0.28,0.15,0.31,0.50,0.37,0.31),
    p.value = c(0.0001,0.0001,0.0001,0.0001,0.011,0.0001,0.0003,0.0627,0.0001,0.0001,0.0001,0.0001),
    alpha = c(51.52,70.82,51.4,38.77,32.4,38.76,31.9,17.42,35.29,54.32,69.13,31.01))
all_genes_correlations

matrix_plot(multiple_regions_correlations,"Correlation statistics for multiple region genes")
matrix_plot(all_genes_correlations,"Correlation statistics for all genes")
setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver")
short_long_fraction_df <- read.table("Fragment end motif fractions of short and long fragments.txt", header = T)
prop_mer_short <- short_long_fraction_df[!grepl("_long",colnames(short_long_fraction_df))]
prop_mer_long <- short_long_fraction_df[!grepl("_short",colnames(short_long_fraction_df))]
prop_mer_short
prop_mer_long
volcano_motif("Cancer <150 bp fragments", "Cancer >150 bp fragments", dif_motif_prop(prop_mer_short,prop_mer_long))

sequence_list <- list("cfChIP" = sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "positive"),
                      "Input" = sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "negative"),
                      "Active_genes" = sequence_end_motif(diff_active_inactive_end_motif(collected_active_inactive_df), 0.05, "positive"),
                      "Inactive_genes" = sequence_end_motif(diff_active_inactive_end_motif(collected_active_inactive_df), 0.05, "negative"),
                      "Long_fragments" = sequence_end_motif(dif_motif_prop(prop_mer_short,prop_mer_long), 0.05, "positive"),
                      "Short_fragments" = sequence_end_motif(dif_motif_prop(prop_mer_short,prop_mer_long), 0.05, "negative"))
sequence_list_reduced <- list("cfChIP" = sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "positive"),
                      "Input" = sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "negative"),
                      "Active_genes" = sequence_end_motif(diff_active_inactive_end_motif(collected_active_inactive_df), 0.05, "positive"),
                      "Inactive_genes" = sequence_end_motif(diff_active_inactive_end_motif(collected_active_inactive_df), 0.05, "negative"))
sequence_end_motif_plot(sequence_list_reduced, c("cfChIP", "Input",
                                         "Active genes", "Inactive genes"))

`%ni%` <- Negate(`%in%`)
sequence_list$cfChIP[sequence_list$cfChIP %in% sequence_list$Short_fragments & 
                         sequence_list$cfChIP %ni% sequence_list$Input &
                         sequence_list$cfChIP %ni% sequence_list$Inactive_genes &
                         sequence_list$cfChIP %ni% sequence_list$Active_genes] 
sequence_list$Inactive_genes[sequence_list$Inactive_genes %in% sequence_list$Long_fragments]
sequence_list$Active_genes[sequence_list$Active_genes %in% sequence_list$cfChIP &
                               sequence_list$Active_genes %ni% sequence_list$Short_fragments &
                               sequence_list$Active_genes %ni% sequence_list$Long_fragments]

enrichment_df
setwd("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver")
active_long_short_df <- read.table("Fragment end motif fractions of short and long fragments in active and inactive genes.txt",
           header = TRUE)
active_long_short_df
df_top_motif <- active_long_short_df[!grepl("bottom",colnames(active_long_short_df))]
df_bottom_motif <- active_long_short_df[!grepl("top",colnames(active_long_short_df))]
df_top_motif_short <- df_top_motif[!grepl("long",colnames(df_top_motif))]
df_top_motif_long <- df_top_motif[!grepl("short",colnames(df_top_motif))]
df_bottom_motif_short <- df_bottom_motif[!grepl("long",colnames(df_bottom_motif))]
df_bottom_motif_long <- df_bottom_motif[!grepl("short",colnames(df_bottom_motif))]

gg1 <- volcano_motif("Active long", "Active short", dif_motif_prop(df_top_motif_long,df_top_motif_short))
gg2 <- volcano_motif("Inactive long", "Inactive short", dif_motif_prop(df_bottom_motif_long,df_bottom_motif_short))
gg3 <- volcano_motif("Inactive long", "Active short", dif_motif_prop(df_bottom_motif_long,df_top_motif_short))
gg4 <- volcano_motif("Inactive short", "Active short", dif_motif_prop(df_bottom_motif_short, df_top_motif_short))
gg5 <- volcano_motif("Inactive long", "Active long", dif_motif_prop(df_bottom_motif_long, df_top_motif_long))
gg6 <- volcano_motif("Inactive short", "Active long", dif_motif_prop(df_bottom_motif_short, df_top_motif_long))
f_plot <- cowplot::plot_grid(plotlist =  list(gg1,gg2,gg3,gg4,gg5,gg6), ncol = 2)
ggsave(filename = "combining fragment lengths and gene activity.png",
       width = 10000, height = 7500, units = "px",
       plot = f_plot,
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots",
       dpi = 350,
       device = "png")

enrichment_cancer <- enrichment_df[!grepl("HC", colnames(enrichment_df))]
enrichment_healthy <- enrichment_df[c(1,14,15,16,17)]

conf(enrichment_cancer, enrichment_healthy,
     v = "cancer", w = "healthy")
healthy_upregulated <- conf(enrichment_cancer, enrichment_healthy,
                            v = "cancer", w = "healthy") %>% filter(Log2FC < 0)
cancer_upregulated <- conf(enrichment_cancer, enrichment_healthy,
                            v = "cancer", w = "healthy") %>% filter(Log2FC > 0)

