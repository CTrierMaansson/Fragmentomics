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
setwd("C:/Users/chris/OneDrive/1PhD/R files/Important excel and text documents")
grs <- gr("AVENIO_genes.txt")
CH01_nuc_range_target <- coverages_nuc("HCC827.bedgraph", z = "GSE71378_CH01_avenio.bed")
target_df <- read.table("Complete table of all targets.txt", header = T)
target_grang <- GRanges(seqnames = target_df$seqnames,
                        ranges = IRanges(start = (target_df$start-1000),
                                         end = (target_df$end + 1000)),
                        strand = target_df$strand)

mcols(target_grang)$SYMBOL <- target_df$SYMBOL
bad <- badgene("Complete table of all targets.txt",25)
multiple_region_genes <- multiple_region("Complete table of all targets.txt")
setwd("C:/Users/chris/OneDrive/1PhD/Fragmentering/endemotiver")
enrichment_df <- read.table("enrichment all.txt", header = T)
pippin_enrichment_df <- read.table("pippin enrichment all.txt", header = T)
coding_exons <- exon_extract(grs,list("BRCA1" = "ENST00000357654.9",
                                      "BRCA2" = "ENST00000380152.8",
                                      "EGFR" = "ENST00000275493.7",
                                      "ERBB2" = "ENST00000584601.5",
                                      "KRAS" = "ENST00000256078.10",
                                      "MET" = "ENST00000318493.11"))
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
gc()
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
gc()

setwd("D:/Healthy input/PosDeduped")
HC.1_input <- bamfile("PosDeduped-Rask_kontrol_1_input.bam",
                      "PosDeduped-Rask_kontrol_1_input.bam.bai")
HC.2_input <- bamfile("PosDeduped-Rask_kontrol_2_input.bam",
                      "PosDeduped-Rask_kontrol_2_input.bam.bai")
HC.3_input <- bamfile("PosDeduped-Rask_kontrol_3_input.bam",
                      "PosDeduped-Rask_kontrol_3_input.bam.bai")
HC.4_input <- bamfile("PosDeduped-Rask_kontrol_4_input.bam",
                      "PosDeduped-Rask_kontrol_4_input.bam.bai")
HC.5_input <- bamfile("PosDeduped-Rask_kontrol_5_input.bam",
                      "PosDeduped-Rask_kontrol_5_input.bam.bai")
HC.6_input <- bamfile("PosDeduped-Rask_kontrol_6_input.bam",
                      "PosDeduped-Rask_kontrol_6_input.bam.bai")
HC.7_input <- bamfile("PosDeduped-Rask_kontrol_7_input.bam",
                      "PosDeduped-Rask_kontrol_7_input.bam.bai")
gc()

setwd("D:/Lung cancer Pippin/PosDeduped")
NAC.1_pippin <- bamfile("PosDeduped-pippin_NAC.1.bam",
                        "PosDeduped-pippin_NAC.1.bam.bai")
NAC.2_pippin <- bamfile("PosDeduped-pippin_NAC.2.bam",
                        "PosDeduped-pippin_NAC.2.bam.bai")
NAC.3_pippin <- bamfile("PosDeduped-pippin_NAC.3.bam",
                        "PosDeduped-pippin_NAC.3.bam.bai")
NAC.4_pippin <- bamfile("PosDeduped-pippin_NAC.4.bam",
                        "PosDeduped-pippin_NAC.4.bam.bai")
NSC.1_pippin <- bamfile("PosDeduped-pippin_NSC.1.bam",
                        "PosDeduped-pippin_NSC.1.bam.bai")
NSC.2_pippin <- bamfile("PosDeduped-pippin_NSC.2.bam",
                        "PosDeduped-pippin_NSC.2.bam.bai")
NSC.3_pippin <- bamfile("PosDeduped-pippin_NSC.3.bam",
                        "PosDeduped-pippin_NSC.3.bam.bai")
SSC.1_pippin <- bamfile("PosDeduped-pippin_SSC.1.bam",
                        "PosDeduped-pippin_SSC.1.bam.bai")
SSC.2_pippin <- bamfile("PosDeduped-pippin_SSC.2.bam",
                        "PosDeduped-pippin_SSC.2.bam.bai")
SSC.3_pippin <- bamfile("PosDeduped-pippin_SSC.3.bam",
                        "PosDeduped-pippin_SSC.3.bam.bai")
SSC.4_pippin <- bamfile("PosDeduped-pippin_SSC.4.bam",
                        "PosDeduped-pippin_SSC.4.bam.bai")
gc()
library(dplyr)
setwd("C:/Users/chris/OneDrive/1PhD/Fragmentering/endemotiver")
prop_mer_input <- read.table("fragment end motif proportions input.txt", header = T)
prop_mer_cfChIP <- read.table("fragment end motif proportions cfChIP.txt", header = T)
prop_mer_pippin <- read.table("fragment end motif proportions pippin.txt", header = T)
count_mer_input <- read.table("fragment end motif counts input.txt", header = T)
count_mer_cfChIP <- read.table("fragment end motif counts cfChIP.txt", header = T)

prop_mer_input_cancer <- prop_mer_input[!grepl("HC.",colnames(prop_mer_input))]
prop_mer_cfChIP_cancer <- prop_mer_cfChIP[!grepl("HC.",colnames(prop_mer_cfChIP))]
prop_mer_input_healthy <- prop_mer_input[grepl("HC.",colnames(prop_mer_input))]
prop_mer_cfChIP_healthy <- prop_mer_cfChIP[grepl("HC.",colnames(prop_mer_cfChIP))]
prop_mer_input_healthy$motif <- prop_mer_input$motif
prop_mer_cfChIP_healthy$motif <- prop_mer_cfChIP$motif

umap_moitf(prop_mer_input, prop_mer_cfChIP, "input", "cfChIP", 8, 15)
umap_moitf(prop_mer_input_cancer, prop_mer_cfChIP_cancer, "Cancer input", "Cancer cfChIP", 9, 10, t = "cancer")
umap_moitf(prop_mer_input_healthy, prop_mer_cfChIP_healthy, "Healthy input", "Healthy cfChIP", 3, 3, t = "healthy")

collected_active_inactive_df <- read.table("fragment end motif proportion of active and inactive fragments.txt", header = T)
cancer_active_inactive <- collected_active_inactive_df[!grepl("HC", colnames(collected_active_inactive_df))]
healthy_active_inactive <- collected_active_inactive_df[grepl("HC", colnames(collected_active_inactive_df))]
healthy_active_inactive$motif <- collected_active_inactive_df$motif
prop_mer_active <- collected_active_inactive_df[!grepl("inactive", colnames(collected_active_inactive_df))]
prop_mer_inactive <- collected_active_inactive_df[grepl("inactive", colnames(collected_active_inactive_df))]
prop_mer_inactive$motif <- collected_active_inactive_df$motif
umap_moitf(prop_mer_active, prop_mer_inactive, "active", "inactive", 3, 32, "Activity")
prop_mer_active_cancer <- prop_mer_active[!grepl("HC", colnames(prop_mer_active))]
prop_mer_inactive_cancer <- prop_mer_inactive[!grepl("HC", colnames(prop_mer_inactive))]
prop_mer_active_healthy <- prop_mer_active[grepl("HC", colnames(prop_mer_active))]
prop_mer_active_healthy$motif <- prop_mer_active$motif
prop_mer_inactive_healthy <- prop_mer_inactive[grepl("HC", colnames(prop_mer_inactive))]
prop_mer_inactive_healthy$motif <- prop_mer_inactive$motif

setwd("C:/Users/chris/OneDrive/1PhD/Fragmentering/endemotiver")
volcano_motif("Input healthy", "Input cancer", diff_gene_motiv_analysis(prop_mer_input_healthy, prop_mer_input_cancer))
volcano_motif("cfChIP healthy", "cfChIP cancer", diff_gene_motiv_analysis(prop_mer_cfChIP_healthy, prop_mer_cfChIP_cancer))
volcano_motif("Input cancer", "cfChIP cancer", dif_motif_prop(prop_mer_input_cancer, prop_mer_cfChIP_cancer))
volcano_motif("Input healthy", "cfChIP healthy", dif_motif_prop(prop_mer_input_healthy, prop_mer_cfChIP_healthy))
volcano_motif("Input", "cfChIP", dif_motif_prop(prop_mer_input, prop_mer_cfChIP))
volcano_motif("Cancer inactive genes", "Cancer active genes", diff_active_inactive_end_motif(cancer_active_inactive))
volcano_motif("Healthy inactive genes", "Healthy active genes", diff_active_inactive_end_motif(healthy_active_inactive))
volcano_motif("Inactive genes", "Active genes", diff_active_inactive_end_motif(collected_active_inactive_df))

volcano_motif("Input", "Size selected", diff_gene_motiv_analysis(prop_mer_input_cancer %>% select(-NSC.4_input), prop_mer_pippin))


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
sequence_list_correlation <- list(c("CCT", "GGA", "GGG", "GCA", "CGT", "CCA",
                                    "GGT", "GCC", "CGG", "CGA", "CCC", "GCT",
                                    "CCG", "GGC", "GCG", "CGC", "TGT", "ACA"),
                                  c("TCA", "AGT", "TGC", "ACG", "ACT", "TGG",
                                    "AGG", "TGA", "ACC", "TAT", "ATA", "TCG",
                                    "AGC", "TCT", "TCC", "AGA", "CAT", "GTA",
                                    "CTA", "GAT", "CAC", "GTG", "CTT", "GTT",
                                    "CAA", "GAA"),
                                  c("GAG", "CTC", "TTT", "AAA", "TAC", "ATG",
                                    "TAA", "ATT", "TTG", "GTC", "AAC", "TAG",
                                    "TTA", "ATC", "AAG", "TTC", "AAT", "GAC",
                                    "CTG", "CAG"))
sequence_end_motif_plot(sequence_list_correlation, c("Left", "Middle",
                                                 "Right"))

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


png("plots/Correlation matrix active cfChIP.png",
    width = 10000, height = 7500,
    units = "px",
    res = 1200
    )
# Code
clusters(prop_mer_cfChIP, prop_mer_active, c("cfChIP", "Active genes"))

# Close device
dev.off()

fragment_length_motif_plot(list(NAC.1_input, NAC.2_input,
                                NAC.3_input, NAC.4_input,
                                NSC.1_input, NSC.2_input,
                                NSC.3_input, NSC.4_input,
                                SSC.1_input, SSC.2_input,
                                SSC.3_input, SSC.4_input,
                                HC.1_input, HC.2_input,
                                HC.3_input, HC.4_input), y = list(sequence_list_correlation[[1]],
                                                                     sequence_list_correlation[[3]]),
                           z = c("GC group", "AT group"),
                           t = "GC group compared to AT group input")
ggsave(filename = "Fragment length based on end motif input.png",
       width = 10000, height = 7500, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots",
       dpi = 1200,
       device = "png")
gc()
fragment_length_motif_plot(list(NAC.1_cfChIP, NAC.2_cfChIP,
                                NAC.3_cfChIP, NAC.4_cfChIP,
                                NSC.1_cfChIP, NSC.2_cfChIP,
                                NSC.3_cfChIP, NSC.4_cfChIP,
                                SSC.1_cfChIP, SSC.2_cfChIP,
                                SSC.3_cfChIP, SSC.4_cfChIP,
                                HC.1_cfChIP, HC.2_cfChIP,
                                HC.3_cfChIP, HC.4_cfChIP), y = list(sequence_list_correlation[[1]],
                                                                  sequence_list_correlation[[3]]),
                           z = c("GC group", "AT group"),
                           t = "GC group compared to AT group cfChIP")
active_motifs <- sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "positive")
active_motifs
inactive_motifs <- sequence_end_motif(dif_motif_prop(prop_mer_input, prop_mer_cfChIP), 0.05, "negative")
inactive_motifs
setwd("C:/Users/chris/OneDrive/1PhD/Fragmentering/endemotiver")
active_motif_fraction_df <- read.table("active motif fraction.txt", header = T)
active_motif_fraction_df

motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NAC.1")],
                             enrichment_df[c("genes", "NAC.1")], "NAC.1", bad, multiple_region_genes)
motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NAC.1")],
                             enrichment_df[c("genes", "NAC.1")], "NAC.1", bad)

motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NAC.2")],
                             enrichment_df[c("genes", "NAC.2")], "NAC.2", bad, multiple_region_genes)
motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NAC.2")],
                             enrichment_df[c("genes", "NAC.2")], "NAC.2", bad)

motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NAC.3")],
                             enrichment_df[c("genes", "NAC.3")], "NAC.3", bad, multiple_region_genes)
motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NAC.3")],
                             enrichment_df[c("genes", "NAC.3")], "NAC.3", bad)

motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NAC.4")],
                             enrichment_df[c("genes", "NAC.4")], "NAC.4", bad, multiple_region_genes)
motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NAC.4")],
                             enrichment_df[c("genes", "NAC.4")], "NAC.4", bad)

motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NSC.1")],
                             enrichment_df[c("genes", "NSC.1")], "NSC.1", bad, multiple_region_genes)
motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NSC.1")],
                             enrichment_df[c("genes", "NSC.1")], "NSC.1", bad)


motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NSC.2")],
                             enrichment_df[c("genes", "NSC.2")], "NSC.2", bad, multiple_region_genes)
motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NSC.2")],
                             enrichment_df[c("genes", "NSC.2")], "NSC.2", bad)

motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NSC.3")],
                             enrichment_df[c("genes", "NSC.3")], "NSC.3", bad, multiple_region_genes)
motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NSC.3")],
                             enrichment_df[c("genes", "NSC.3")], "NSC.3", bad)

motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NSC.4")],
                             enrichment_df[c("genes", "NSC.4")], "NSC.4", bad, multiple_region_genes)
motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "NSC.4")],
                             enrichment_df[c("genes", "NSC.4")], "NSC.4", bad)

motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "SSC.1")],
                             enrichment_df[c("genes", "SSC.1")], "SSC.1", bad, multiple_region_genes)
motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "SSC.1")],
                             enrichment_df[c("genes", "SSC.1")], "SSC.1", bad)

motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "SSC.2")],
                             enrichment_df[c("genes", "SSC.2")], "SSC.2", bad, multiple_region_genes)
motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "SSC.2")],
                             enrichment_df[c("genes", "SSC.2")], "SSC.2", bad)

motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "SSC.3")],
                             enrichment_df[c("genes", "SSC.3")], "SSC.3", bad, multiple_region_genes)
motif_fraction_vs_enrichment(active_motif_fraction_df[c("genes", "SSC.3")],
                             enrichment_df[c("genes", "SSC.3")], "SSC.3", bad)

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
rho_multiple <- c(0.48, 0.44, 0.30, 0.21, -0.034, 0.024, 0.23, 0.061, 0.193, 0.55, 0.146, 0.13)
p_multiple <- c(0.0078, 0.016, 0.1025, 0.2505, 0.858, 0.899, 0.2137, 0.747, 0.3058, 0.0021, 0.4407, 0.478)
alpha_multiple <- c(60.79, 53.63, 39.94, 18.47, 17.98, 19, 21.08, 2.37, 31.79, 58.32, 60.28, 15.12)
rho_all <- c(0.42, 0.56, 0.41, 0.39, 0.17, 0.28, 0.25, 0.17, 0.298, 0.517, 0.39, 0.28)
p_all <- c(0.0001, 0.0001, 0.0001, 0.0001, 0.0366, 0.0001, 0.0014, 0.0318, 0.0001, 0.0001, 0.0001, 0.0001)
alpha_all <- c(47.27, 64.48, 44.55, 31.27, 24.43, 31.59, 25.1, 15.51, 30.61, 53.25, 66.22, 25.01)
multiple_regions_correlations <- data.frame(
    sample = colnames(enrichment_df)[2:(length(colnames(enrichment_df))-4)],
    rho = rho_multiple,
    p.value = p_multiple,
    alpha = alpha_multiple)
multiple_regions_correlations
all_genes_correlations <- data.frame(
    sample = colnames(enrichment_df)[2:(length(colnames(enrichment_df))-4)],
    rho = rho_all,
    p.value = p_all,
    alpha = alpha_all)
all_genes_correlations


matrix_plot(multiple_regions_correlations,"Correlation statistics for multiple region genes")
matrix_plot(all_genes_correlations,"Correlation statistics for all genes")
setwd("C:/Users/chris/OneDrive/1PhD/Fragmentering/endemotiver")
short_long_fraction_df <- read.table("Fragment end motif fractions of short and long fragments.txt", header = T)
prop_mer_short <- short_long_fraction_df[!grepl("_long",colnames(short_long_fraction_df))]
prop_mer_long <- short_long_fraction_df[!grepl("_short",colnames(short_long_fraction_df))]
prop_mer_short
prop_mer_long
volcano_motif("Cancer short", "Cancer long", dif_motif_prop(prop_mer_short,prop_mer_long))

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
setwd("C:/Users/chris/OneDrive/1PhD/Fragmentering/endemotiver")
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
       width = 7500, height = 10000, units = "px",
       plot = f_plot,
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots",
       dpi = 350,
       device = "png")
ggsave(filename = "Cancer active long vs. active short.png",
       width = 10000, height = 7500, units = "px",
       plot = gg1,
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots",
       dpi = 1200,
       device = "png")
ggsave(filename = "Cancer inactive long vs. inactive short.png",
       width = 10000, height = 7500, units = "px",
       plot = gg2,
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots",
       dpi = 1200,
       device = "png")
ggsave(filename = "Cancer inactive long vs. active short.png",
       width = 10000, height = 7500, units = "px",
       plot = gg3,
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots",
       dpi = 1200,
       device = "png")
ggsave(filename = "Cancer inactive short vs. active short.png",
       width = 10000, height = 7500, units = "px",
       plot = gg4,
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots",
       dpi = 1200,
       device = "png")
ggsave(filename = "Cancer inactive long vs. active long.png",
       width = 10000, height = 7500, units = "px",
       plot = gg5,
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots",
       dpi = 1200,
       device = "png")
ggsave(filename = "Cancer inactive short vs. active long.png",
       width = 10000, height = 7500, units = "px",
       plot = gg6,
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots",
       dpi = 1200,
       device = "png")



enrichment_cancer <- enrichment_df[!grepl("HC", colnames(enrichment_df))]
enrichment_healthy <- enrichment_df[c(1,14,15,16,17)]

conf(enrichment_cancer, enrichment_healthy,
     v = "cancer", w = "healthy")
healthy_upregulated <- conf(enrichment_cancer, enrichment_healthy,
                            v = "cancer", w = "healthy") %>% filter(Log2FC < 0)
cancer_upregulated <- conf(enrichment_cancer, enrichment_healthy,
                            v = "cancer", w = "healthy") %>% filter(Log2FC > 0)

cancer_input <- list(NAC.1_input, NAC.2_input,
                     NAC.3_input, NAC.4_input,
                     NSC.1_input, NSC.2_input,
                     NSC.3_input, NSC.4_input,
                     SSC.1_input, SSC.2_input,
                     SSC.3_input, SSC.4_input)
healthy_input <- list(HC.1_input, HC.2_input,
                      HC.3_input, HC.4_input)

healthy_upregulated_cancer_files <- diff_gene_epigenetic_analysis(healthy_upregulated, 
                                                                  grs, 
                                                                  cancer_input,
                                                                  150,
                                                                  3)
healthy_upregulated_healthy_files <- diff_gene_epigenetic_analysis(healthy_upregulated, 
                                                                   grs, 
                                                                   healthy_input,
                                                                   150,
                                                                   3)
cancer_upregulated_cancer_files <- diff_gene_epigenetic_analysis(cancer_upregulated, 
                                                                 grs, 
                                                                 cancer_input,
                                                                 150,
                                                                 3)
cancer_upregulated_healthy_files <- diff_gene_epigenetic_analysis(cancer_upregulated, 
                                                                  grs, 
                                                                  healthy_input,
                                                                  150,
                                                                  3)
colnames(cancer_upregulated_cancer_files[[2]])[1] <- "motif"
colnames(cancer_upregulated_healthy_files[[2]])[1] <- "motif"
colnames(healthy_upregulated_cancer_files[[2]])[1] <- "motif"
colnames(healthy_upregulated_healthy_files[[2]])[1] <- "motif"

volcano_motif("Healthy","Cancer",
              diff_gene_motiv_analysis(cancer_upregulated_healthy_files[[2]],cancer_upregulated_cancer_files[[2]]))
volcano_motif("Healthy","Cancer",
              diff_gene_motiv_analysis(healthy_upregulated_healthy_files[[2]],healthy_upregulated_cancer_files[[2]]))


shannon_entropy_plot(prop_mer_active,
                     prop_mer_inactive)
ggsave(filename = "Shannon entropy active inactive.png",
       width = 7500, height = 10000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots",
       dpi = 1200,
       device = "png")
shannon_entropy_plot(prop_mer_short,
                     prop_mer_long)
ggsave(filename = "Shannon entropy short long.png",
       width = 7500, height = 10000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots",
       dpi = 1200,
       device = "png")

transposed_shannon_entropy_plot(prop_mer_active,
                                prop_mer_inactive)
transposed_shannon_entropy_plot(prop_mer_short,
                                prop_mer_long)

name_list <- c("Input", "cfChIP",
               "Inactive", "Active",
               "Long", "Short")
stack_bar(sequence_list_correlation,
          cancer_dfs,
          name_list, 
          c("L", "M", "R"))
filter_motifs <- prop_mer_active$motif[prop_mer_active$motif %ni% c(active_motifs,
                                                                    inactive_motifs)]
stack_bar(list(active_motifs,
               filter_motifs,
               inactive_motifs),
          cancer_dfs,
          name_list, 
          c("Active", "Filtered", "Inacitve"))


cancer_dfs <- list(prop_mer_input_cancer,
                   prop_mer_cfChIP_cancer,
                   prop_mer_inactive_cancer,
                   prop_mer_active_cancer,
                   prop_mer_long,
                   prop_mer_short)

fragment_length(list(NAC.1_input, NAC.2_input,
                     NAC.3_input, NAC.4_input,
                     NSC.1_input, NSC.2_input,
                     NSC.3_input, NSC.4_input,
                     SSC.1_input, SSC.2_input,
                     SSC.3_input, SSC.4_input),
                list(HC.1_input, HC.2_input,
                     HC.3_input, HC.4_input),
                c("Cancer", "Healthy"), "Input")

gc()
fragment_length(list(NAC.1_input,
                     NAC.3_input, NAC.4_input,
                     NSC.1_input, NSC.2_input,
                     SSC.1_input, SSC.2_input,
                     SSC.3_input, SSC.4_input),
                list(NAC.2_input,NSC.3_input, NSC.4_input),
                c("ctDNA positive", "ctDNA negative"), "Input")

fragment_length(list(NAC.1_cfChIP, NAC.2_cfChIP,
                     NAC.3_cfChIP, NAC.4_cfChIP,
                     NSC.1_cfChIP, NSC.2_cfChIP,
                     NSC.3_cfChIP, NSC.4_cfChIP,
                     SSC.1_cfChIP, SSC.2_cfChIP,
                     SSC.3_cfChIP, SSC.4_cfChIP),
                list(NAC.1_input, NAC.2_input,
                     NAC.3_input, NAC.4_input,
                     NSC.1_input, NSC.2_input,
                     NSC.3_input, NSC.4_input,
                     SSC.1_input, SSC.2_input,
                     SSC.3_input, SSC.4_input),
                c("cfChIP", "Input"),"Cancer")
gc()
fragment_length(list(HC.1_cfChIP, HC.2_cfChIP,
                     HC.3_cfChIP, HC.4_cfChIP),
                list(HC.1_input, HC.2_input,
                     HC.3_input, HC.4_input),
                c("cfChIP", "Input"), "Healthy")
fragment_length(list(NAC.1_cfChIP, NAC.2_cfChIP,
                     NAC.3_cfChIP, NAC.4_cfChIP,
                     NSC.1_cfChIP, NSC.2_cfChIP,
                     NSC.3_cfChIP, NSC.4_cfChIP,
                     SSC.1_cfChIP, SSC.2_cfChIP,
                     SSC.3_cfChIP, SSC.4_cfChIP),
                list(HC.1_cfChIP, HC.2_cfChIP,
                     HC.3_cfChIP, HC.4_cfChIP),
                c("Cancer", "Healthy"),"cfChIP")

proportion_sub150_cancer_input <- proportion_sub_df(list(NAC.1_input, NAC.2_input,
                                                         NAC.3_input, NAC.4_input,
                                                         NSC.1_input, NSC.2_input,
                                                         NSC.3_input, NSC.4_input,
                                                         SSC.1_input, SSC.2_input,
                                                         SSC.3_input, SSC.4_input),
                                                    150,
                                                    c("NAC.1_input", "NAC.2_input",
                                                      "NAC.3_input", "NAC.4_input",
                                                      "NSC.1_input", "NSC.2_input",
                                                      "NSC.3_input", "NSC.4_input",
                                                      "SSC.1_input", "SSC.2_input",
                                                      "SSC.3_input", "SSC.4_input"))

proportion_sub150_healthy_input <- proportion_sub_df(list(HC.1_input, HC.2_input,
                                                          HC.3_input, HC.4_input,
                                                          HC.5_input, HC.6_input,
                                                          HC.7_input),
                                                    150,
                                                    c("HC.1_input", "HC.2_input",
                                                      "HC.3_input", "HC.4_input",
                                                      "HC.5_input", "HC.6_input",
                                                      "HC.7_input"))

proportion_sub150_ct_pos_input <- proportion_sub_df(list(NAC.1_input,
                                                         NAC.3_input, NAC.4_input,
                                                         NSC.1_input, NSC.2_input,
                                                         SSC.1_input, SSC.2_input,
                                                         SSC.3_input, SSC.4_input),
                                                     150,
                                                     c("NAC.1_input",
                                                       "NAC.3_input", "NAC.4_input",
                                                       "NSC.1_input", "NSC.2_input",
                                                       "SSC.1_input", "SSC.2_input",
                                                       "SSC.3_input", "SSC.4_input"))

proportion_sub150_ct_neg_input <- proportion_sub_df(list(NAC.2_input,NSC.3_input, 
                                                         NSC.4_input),
                                                    150,
                                                    c("NAC.2_input", "NSC.3_input", 
                                                      "NSC.4_input"))
proportion_sub_150_cancer_cfChIP <- proportion_sub_df(list(NAC.1_cfChIP, NAC.2_cfChIP,
                                                           NAC.3_cfChIP, NAC.4_cfChIP,
                                                           NSC.1_cfChIP, NSC.2_cfChIP,
                                                           NSC.3_cfChIP, NSC.4_cfChIP,
                                                           SSC.1_cfChIP, SSC.2_cfChIP,
                                                           SSC.3_cfChIP, SSC.4_cfChIP),
                                                      150,
                                                      c("NAC.1_cfChIP", "NAC.2_cfChIP",
                                                        "NAC.3_cfChIP", "NAC.4_cfChIP",
                                                        "NSC.1_cfChIP", "NSC.2_cfChIP",
                                                        "NSC.3_cfChIP", "NSC.4_cfChIP",
                                                        "SSC.1_cfChIP", "SSC.2_cfChIP",
                                                        "SSC.3_cfChIP", "SSC.4_cfChIP"))

proportion_sub150_healthy_cfChIP <- proportion_sub_df(list(HC.1_cfChIP, HC.2_cfChIP,
                                                          HC.3_cfChIP, HC.4_cfChIP),
                                                     150,
                                                     c("HC.1_cfChIP", "HC.2_cfChIP",
                                                       "HC.3_cfChIP", "HC.4_cfChIP"))
proportion_sub150_healthy_input_w_ChIP <- proportion_sub_df(list(HC.1_input, HC.2_input,
                                                                 HC.3_input, HC.4_input),
                                                            150,
                                                            c("HC.1_input", "HC.2_input",
                                                              "HC.3_input", "HC.4_input"))

proportion_sub_boxplot(list(proportion_sub150_healthy_input,
                            proportion_sub150_cancer_input),
                       c("Healthy", "Cancer"),
                       "Unpaired")
proportion_sub_boxplot(list(proportion_sub150_healthy_input,
                            proportion_sub150_ct_neg_input,
                            proportion_sub150_ct_pos_input),
                       c("Healthy",
                         "ctDNA negative",
                         "ctDNA positive"),
                       "Unpaired")
proportion_sub_boxplot(list(proportion_sub150_cancer_input,
                            proportion_sub_150_cancer_cfChIP),
                       c("Cancer input", "Cancer cfChIP"),
                       "Paired")

proportion_sub_boxplot(list(proportion_sub150_healthy_input_w_ChIP,
                            proportion_sub150_healthy_cfChIP),
                       c("healthy input", "Healthy cfChIP"),
                       "Paired")
proportion_sub_boxplot(list(proportion_sub150_healthy_cfChIP,
                            proportion_sub_150_cancer_cfChIP),
                       c("healthy cfChIP", "Cancer cfChIP"),
                       "Unpaired")


named_cancer_input_list <- list(NAC.1_input = NAC.1_input, 
                                NAC.2_input = NAC.2_input,
                                NAC.3_input = NAC.3_input, 
                                NAC.4_input = NAC.4_input,
                                NSC.1_input = NSC.1_input, 
                                NSC.2_input = NSC.2_input,
                                NSC.3_input = NSC.3_input, 
                                NSC.4_input = NSC.4_input,
                                SSC.1_input = SSC.1_input, 
                                SSC.2_input = SSC.2_input,
                                SSC.3_input = SSC.3_input, 
                                SSC.4_input = SSC.4_input)
named_cancer_cfChIP_list <- list(NAC.1_cfChIP = NAC.1_cfChIP, 
                                 NAC.2_cfChIP = NAC.2_cfChIP,
                                 NAC.3_cfChIP = NAC.3_cfChIP, 
                                 NAC.4_cfChIP = NAC.4_cfChIP,
                                 NSC.1_cfChIP = NSC.1_cfChIP, 
                                 NSC.2_cfChIP = NSC.2_cfChIP,
                                 NSC.3_cfChIP = NSC.3_cfChIP, 
                                 NSC.4_cfChIP = NSC.4_cfChIP,
                                 SSC.1_cfChIP = SSC.1_cfChIP, 
                                 SSC.2_cfChIP = SSC.2_cfChIP,
                                 SSC.3_cfChIP = SSC.3_cfChIP, 
                                 SSC.4_cfChIP = SSC.4_cfChIP)
setwd("C:/Users/chris/OneDrive/1PhD/Fragmentering/endemotiver")

fragment_length_wt_ctdna_input <- fragment_length_wt_ctdna_df("BL mutations input files.txt",
                                                              named_cancer_input_list)
fragment_length_wt_ctdna_input_vessies <- fragment_length_wt_ctdna_vessies_df("BL mutations input files.txt",
                                                              named_cancer_input_list)
fragment_length_wt_ctdna_plot(fragment_length_wt_ctdna_input)

MAF_in_bins(fragment_length_wt_ctdna_input)

NAC.1_input_ct_wt <-  fragment_length_wt_ctdna_patient("BL mutations input files.txt",
                                 list(NAC.1_input = NAC.1_input))
NAC.2_input_ct_wt <-  fragment_length_wt_ctdna_patient("BL mutations input files.txt",
                                                       list(NAC.2_input = NAC.2_input))
NAC.3_input_ct_wt <-  fragment_length_wt_ctdna_patient("BL mutations input files.txt",
                                                       list(NAC.3_input = NAC.3_input))
NAC.4_input_ct_wt <-  fragment_length_wt_ctdna_patient("BL mutations input files.txt",
                                                       list(NAC.4_input = NAC.4_input))
NSC.1_input_ct_wt <-  fragment_length_wt_ctdna_patient("BL mutations input files.txt",
                                                       list(NSC.1_input = NSC.1_input))
NSC.2_input_ct_wt <-  fragment_length_wt_ctdna_patient("BL mutations input files.txt",
                                                       list(NSC.2_input = NSC.2_input))
NSC.3_input_ct_wt <-  fragment_length_wt_ctdna_patient("BL mutations input files.txt",
                                                       list(NSC.3_input = NSC.3_input))
NSC.4_input_ct_wt <-  fragment_length_wt_ctdna_patient("BL mutations input files.txt",
                                                       list(NSC.4_input = NSC.4_input))
SSC.1_input_ct_wt <-  fragment_length_wt_ctdna_patient("BL mutations input files.txt",
                                                       list(SSC.1_input = SSC.1_input))
SSC.2_input_ct_wt <-  fragment_length_wt_ctdna_patient("BL mutations input files.txt",
                                                       list(SSC.2_input = SSC.2_input))
SSC.3_input_ct_wt <-  fragment_length_wt_ctdna_patient("BL mutations input files.txt",
                                                       list(SSC.3_input = SSC.3_input))
SSC.4_input_ct_wt <-  fragment_length_wt_ctdna_patient("BL mutations input files.txt",
                                                       list(SSC.4_input = SSC.4_input))
ct_wt_input_list <- list(NAC.1_input_ct_wt,
     NAC.2_input_ct_wt,
     NAC.3_input_ct_wt,
     NAC.4_input_ct_wt,
     NSC.1_input_ct_wt,
     NSC.2_input_ct_wt,
     NSC.3_input_ct_wt,
     NSC.4_input_ct_wt,
     SSC.1_input_ct_wt,
     SSC.2_input_ct_wt,
     SSC.3_input_ct_wt,
     SSC.4_input_ct_wt)

fragment_length_wt_ctdna_boxplot(ct_wt_input_list)
gc()
collected_fragment_length_active_inactive <- read.table("Fragment lengths active inactive.txt", header = T)
collected_fragment_length_active_inactive
collected_fragment_length_active_inactive_healthy <- read.table("Fragment lengths active inactive healthy.txt", header = T)

fragment_length_active_inactive_plot(collected_fragment_length_active_inactive,
                                     "Cancer")
fragment_length_active_inactive_plot(collected_fragment_length_active_inactive_healthy,
                                     "Healthy")



length_list_active_inactive_cancer <- list(fragment_length_active_inactive_df(enrichment_df,grs,"NAC.1","D:/Lung cancer input/PosDeduped/PosDeduped-A-1279-input.bam"),
     fragment_length_active_inactive_df(enrichment_df,grs,"NAC.2","D:/Lung cancer input/PosDeduped/PosDeduped-B-1288-input.bam"),
     fragment_length_active_inactive_df(enrichment_df,grs,"NAC.3","D:/Lung cancer input/PosDeduped/PosDeduped-C-1475-input.bam"),
     fragment_length_active_inactive_df(enrichment_df,grs,"NAC.4","D:/Lung cancer input/PosDeduped/PosDeduped-D-1578-input.bam"),
     fragment_length_active_inactive_df(enrichment_df,grs,"NSC.1","D:/Lung cancer input/PosDeduped/PosDeduped-E-439-input.bam"),
     fragment_length_active_inactive_df(enrichment_df,grs,"NSC.2","D:/Lung cancer input/PosDeduped/PosDeduped-F-1449-input.bam"),
     fragment_length_active_inactive_df(enrichment_df,grs,"NSC.3","D:/Lung cancer input/PosDeduped/PosDeduped-I-645-input.bam"),
     fragment_length_active_inactive_df(enrichment_df,grs,"NSC.4","D:/Lung cancer input/PosDeduped/PosDeduped-J-1663-input.bam"),
     fragment_length_active_inactive_df(enrichment_df,grs,"SSC.1","D:/Lung cancer input/PosDeduped/PosDeduped-G-514-input.bam"),
     fragment_length_active_inactive_df(enrichment_df,grs,"SSC.2","D:/Lung cancer input/PosDeduped/PosDeduped-H-1169-input.bam"),
     fragment_length_active_inactive_df(enrichment_df,grs,"SSC.3","D:/Lung cancer input/PosDeduped/PosDeduped-K-440-input.bam"),
     fragment_length_active_inactive_df(enrichment_df,grs,"SSC.4","D:/Lung cancer input/PosDeduped/PosDeduped-L-1100-input.bam"))

length_list_active_inactive_healthy <- list(fragment_length_active_inactive_df(enrichment_df,grs,"HC.1","D:/Healthy input/PosDeduped/PosDeduped-Rask_kontrol_1_input.bam"),
                                           fragment_length_active_inactive_df(enrichment_df,grs,"HC.2","D:/Healthy input/PosDeduped/PosDeduped-Rask_kontrol_2_input.bam"),
                                           fragment_length_active_inactive_df(enrichment_df,grs,"HC.3","D:/Healthy input/PosDeduped/PosDeduped-Rask_kontrol_3_input.bam"),
                                           fragment_length_active_inactive_df(enrichment_df,grs,"HC.4","D:/Healthy input/PosDeduped/PosDeduped-Rask_kontrol_4_input.bam"))
                                           

fragment_length_inactive_cancer <- fragment_length_active_inactive_sub_df(length_list_active_inactive_cancer, 150, "Inactive",
                                                                   c("NAC.1","NAC.2","NAC.3","NAC.4",
                                                                     "NSC.1","NSC.2","NSC.3","NSC.4",
                                                                     "SSC.1","SSC.2","SSC.3","SSC.4"))
fragment_length_active_cancer <- fragment_length_active_inactive_sub_df(length_list_active_inactive_cancer, 150, "Active",
                                                                 c("NAC.1","NAC.2","NAC.3","NAC.4",
                                                                   "NSC.1","NSC.2","NSC.3","NSC.4",
                                                                   "SSC.1","SSC.2","SSC.3","SSC.4"))

proportion_sub_boxplot(list(fragment_length_inactive_cancer,fragment_length_active_cancer),
                       c("Cancer Inactive", "Cancer Active"),
                       "Paired")

fragment_length_inactive_healthy <- fragment_length_active_inactive_sub_df(length_list_active_inactive_healthy, 150, "Inactive",
                                                                   c("HC.1","HC.2","HC.3","HC.4"))
fragment_length_active_healthy <- fragment_length_active_inactive_sub_df(length_list_active_inactive_healthy, 150, "Active",
                                                                           c("HC.1","HC.2","HC.3","HC.4"))

proportion_sub_boxplot(list(fragment_length_inactive_healthy,fragment_length_active_healthy),
                       c("Healthy Inactive", "healthy Active"),
                       "Paired")

collected_sub150_fraction_quartiles <- read.table("Sub150 fraction in cfChIP quatiles.txt", header = T)

sub150_fraction_quartiles_plot(collected_sub150_fraction_quartiles)


####input vs. cfChIP####
setwd("C:/Users/chris/OneDrive/1PhD/Fragmentering/endemotiver")
active_motifs_cancer <- sequence_end_motif(dif_motif_prop(prop_mer_input_cancer, prop_mer_cfChIP_cancer), 0.05, "positive")
active_motifs_cancer
inactive_motifs_cancer <- sequence_end_motif(dif_motif_prop(prop_mer_input_cancer, prop_mer_cfChIP_cancer), 0.05, "negative")
inactive_motifs_cancer

collected_active_motif_fraction_quartiles <- read.table("Active motif cancer fraction in cfChIP quatiles.txt", header = T)

active_motif_fraction_quartiles_plot(collected_active_motif_fraction_quartiles, "Active")

collected_inactive_motif_fraction_quartiles <- read.table("Inactive motif cancer fraction in cfChIP quatiles.txt", header = T)

active_motif_fraction_quartiles_plot(collected_inactive_motif_fraction_quartiles, "Inactive")

####Active vs. inactive####
active_active <- sequence_end_motif(dif_motif_prop(prop_mer_inactive, prop_mer_active), 0.05, "positive")
active_active
inactive_inactive <- sequence_end_motif(dif_motif_prop(prop_mer_inactive, prop_mer_active), 0.05, "negative")
inactive_inactive

collected_active_active_fraction_quartiles <- read.table("Active active fraction in cfChIP quatiles.txt", header = T)

active_motif_fraction_quartiles_plot(collected_active_active_fraction_quartiles, "Active")

collected_inactive_inactive_fraction_quartiles <- read.table("Inactive inactive fraction in cfChIP quatiles.txt", header = T)

active_motif_fraction_quartiles_plot(collected_inactive_inactive_fraction_quartiles, "Inactive")

active_active_cancer <- sequence_end_motif(dif_motif_prop(prop_mer_inactive_cancer, prop_mer_active_cancer), 0.05, "positive")
active_active_cancer
inactive_inactive_cancer <- sequence_end_motif(dif_motif_prop(prop_mer_inactive_cancer, prop_mer_active_cancer), 0.05, "negative")
inactive_inactive_cancer

collected_active_active_cancer_fraction_quartiles <- read.table("Active active cancer fraction in cfChIP quatiles.txt", header = T)

active_motif_fraction_quartiles_plot(collected_active_active_cancer_fraction_quartiles, "Active")

collected_inactive_inactive_cancer_fraction_quartiles <- read.table("Inactive inactive cancer fraction in cfChIP quatiles.txt", header = T)

active_motif_fraction_quartiles_plot(collected_inactive_inactive_cancer_fraction_quartiles, "Inactive")

####Combining length and motif####

#input vs. cfChIP
collected_sub150_active_motif_fraction_quartiles <- read.table("Sub 150 with active motif cancer fraction in cfChIP quatiles.txt", header = T)

sub150_active_motif_fraction_quartiles_plot(collected_sub150_active_motif_fraction_quartiles, "Active")

#Active vs. inactive
collected_sub150_active_active_fraction_quartiles <- read.table("Sub 150 with active active fraction in cfChIP quatiles.txt", header = T)

sub150_active_motif_fraction_quartiles_plot(collected_sub150_active_active_fraction_quartiles, "Active")


ggsave(filename = "Sub150 active motif fraction Cancer cfChIP quartiles.png",
       width = 10000, height = 7500, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots",
       dpi = 1200,
       device = "png")


NAC.1_nuc_dist_df_bed <- bed_based_cleavage_nuc_dist_data(top_bottom_genes(enrichment_df,"NAC.1"),
                                                          NAC.1_input, grs, CH01_nuc_range_target, target_grang)
gc()
NAC.2_nuc_dist_df_bed <- bed_based_cleavage_nuc_dist_data(top_bottom_genes(enrichment_df,"NAC.2"),
                                                          NAC.2_input, grs, CH01_nuc_range_target, target_grang)
gc()
NAC.3_nuc_dist_df_bed <- bed_based_cleavage_nuc_dist_data(top_bottom_genes(enrichment_df,"NAC.3"),
                                                          NAC.3_input, grs, CH01_nuc_range_target, target_grang)
gc()
NAC.4_nuc_dist_df_bed <- bed_based_cleavage_nuc_dist_data(top_bottom_genes(enrichment_df,"NAC.4"),
                                                          NAC.4_input, grs, CH01_nuc_range_target, target_grang)
gc()
NSC.1_nuc_dist_df_bed <- bed_based_cleavage_nuc_dist_data(top_bottom_genes(enrichment_df,"NSC.1"),
                                                          NSC.1_input, grs, CH01_nuc_range_target, target_grang)
gc()
NSC.2_nuc_dist_df_bed <- bed_based_cleavage_nuc_dist_data(top_bottom_genes(enrichment_df,"NSC.2"),
                                                          NSC.2_input, grs, CH01_nuc_range_target, target_grang)
gc()
NSC.3_nuc_dist_df_bed <- bed_based_cleavage_nuc_dist_data(top_bottom_genes(enrichment_df,"NSC.3"),
                                                          NSC.3_input, grs, CH01_nuc_range_target, target_grang)
gc()
NSC.4_nuc_dist_df_bed <- bed_based_cleavage_nuc_dist_data(top_bottom_genes(enrichment_df,"NSC.4"),
                                                          NSC.4_input, grs, CH01_nuc_range_target, target_grang)
gc()
SSC.1_nuc_dist_df_bed <- bed_based_cleavage_nuc_dist_data(top_bottom_genes(enrichment_df,"SSC.1"),
                                                          SSC.1_input, grs, CH01_nuc_range_target, target_grang)
gc()
SSC.2_nuc_dist_df_bed <- bed_based_cleavage_nuc_dist_data(top_bottom_genes(enrichment_df,"SSC.2"),
                                                          SSC.2_input, grs, CH01_nuc_range_target, target_grang)
gc()
SSC.3_nuc_dist_df_bed <- bed_based_cleavage_nuc_dist_data(top_bottom_genes(enrichment_df,"SSC.3"),
                                                          SSC.3_input, grs, CH01_nuc_range_target, target_grang)
gc()
SSC.4_nuc_dist_df_bed <- bed_based_cleavage_nuc_dist_data(top_bottom_genes(enrichment_df,"SSC.4"),
                                                          SSC.4_input, grs, CH01_nuc_range_target, target_grang)

collected_nuc_dist_df_bed <- rbind(NAC.1_nuc_dist_df_bed,
                                   NAC.2_nuc_dist_df_bed,
                                   NAC.3_nuc_dist_df_bed,
                                   NAC.4_nuc_dist_df_bed,
                                   NSC.1_nuc_dist_df_bed,
                                   NSC.2_nuc_dist_df_bed,
                                   NSC.3_nuc_dist_df_bed,
                                   NSC.4_nuc_dist_df_bed,
                                   SSC.1_nuc_dist_df_bed,
                                   SSC.2_nuc_dist_df_bed,
                                   SSC.3_nuc_dist_df_bed,
                                   SSC.4_nuc_dist_df_bed) %>% 
    mutate(sample = c(rep("NAC.1", nrow(NAC.1_nuc_dist_df_bed)),
                      rep("NAC.2", nrow(NAC.2_nuc_dist_df_bed)),
                      rep("NAC.3", nrow(NAC.3_nuc_dist_df_bed)),
                      rep("NAC.4", nrow(NAC.4_nuc_dist_df_bed)),
                      rep("NSC.1", nrow(NSC.1_nuc_dist_df_bed)),
                      rep("NSC.2", nrow(NSC.2_nuc_dist_df_bed)),
                      rep("NSC.3", nrow(NSC.3_nuc_dist_df_bed)),
                      rep("NSC.4", nrow(NSC.4_nuc_dist_df_bed)),
                      rep("SSC.1", nrow(SSC.1_nuc_dist_df_bed)),
                      rep("SSC.2", nrow(SSC.2_nuc_dist_df_bed)),
                      rep("SSC.3", nrow(SSC.3_nuc_dist_df_bed)),
                      rep("SSC.4", nrow(SSC.4_nuc_dist_df_bed))))

distance_nucleosome_plot(collected_nuc_dist_df_bed, TRUE)

distance_nucleosome_plot(collected_nuc_dist_df_bed, FALSE)

distance_nucleosome_boxplot(collected_nuc_dist_df_bed)

active_motifs_cancer <- sequence_end_motif(dif_motif_prop(prop_mer_input_cancer, prop_mer_cfChIP_cancer), 0.05, "positive")
inactive_motifs_cancer <- sequence_end_motif(dif_motif_prop(prop_mer_input_cancer, prop_mer_cfChIP_cancer), 0.05, "negative")
gc()
cfChIP_motif_lengths_NAC.1 <- length_fragment_with_motif_df(NAC.1_input,active_motifs_cancer,inactive_motifs_cancer)
cfChIP_motif_lengths_NAC.2 <- length_fragment_with_motif_df(NAC.2_input,active_motifs_cancer,inactive_motifs_cancer)
cfChIP_motif_lengths_NAC.3 <- length_fragment_with_motif_df(NAC.3_input,active_motifs_cancer,inactive_motifs_cancer)
cfChIP_motif_lengths_NAC.4 <- length_fragment_with_motif_df(NAC.4_input,active_motifs_cancer,inactive_motifs_cancer)
gc()
cfChIP_motif_lengths_NSC.1 <- length_fragment_with_motif_df(NSC.1_input,active_motifs_cancer,inactive_motifs_cancer)
cfChIP_motif_lengths_NSC.2 <- length_fragment_with_motif_df(NSC.2_input,active_motifs_cancer,inactive_motifs_cancer)
cfChIP_motif_lengths_NSC.3 <- length_fragment_with_motif_df(NSC.3_input,active_motifs_cancer,inactive_motifs_cancer)
cfChIP_motif_lengths_NSC.4 <- length_fragment_with_motif_df(NSC.4_input,active_motifs_cancer,inactive_motifs_cancer)
gc()
cfChIP_motif_lengths_SSC.1 <- length_fragment_with_motif_df(SSC.1_input,active_motifs_cancer,inactive_motifs_cancer)
cfChIP_motif_lengths_SSC.2 <- length_fragment_with_motif_df(SSC.2_input,active_motifs_cancer,inactive_motifs_cancer)
cfChIP_motif_lengths_SSC.3 <- length_fragment_with_motif_df(SSC.3_input,active_motifs_cancer,inactive_motifs_cancer)
cfChIP_motif_lengths_SSC.4 <- length_fragment_with_motif_df(SSC.4_input,active_motifs_cancer,inactive_motifs_cancer)
gc()
active_motif_lengths_NAC.1 <- motif_lengths_NAC.1
active_motif_lengths_NAC.2 <- motif_lengths_NAC.2
active_motif_lengths_NAC.3 <- motif_lengths_NAC.3
active_motif_lengths_NAC.4 <- motif_lengths_NAC.4
active_motif_lengths_NSC.1 <- motif_lengths_NSC.1
active_motif_lengths_NSC.2 <- motif_lengths_NSC.2
active_motif_lengths_NSC.3 <- motif_lengths_NSC.3
active_motif_lengths_NSC.4 <- motif_lengths_NSC.4
active_motif_lengths_SSC.1 <- motif_lengths_SSC.1
active_motif_lengths_SSC.2 <- motif_lengths_SSC.2
active_motif_lengths_SSC.3 <- motif_lengths_SSC.3
active_motif_lengths_SSC.4 <- motif_lengths_SSC.4


collected_cfChIP_motif_lengths <- rbind(cfChIP_motif_lengths_NAC.1,
                                        cfChIP_motif_lengths_NAC.2,
                                        cfChIP_motif_lengths_NAC.3,
                                        cfChIP_motif_lengths_NAC.4,
                                        cfChIP_motif_lengths_NSC.1,
                                        cfChIP_motif_lengths_NSC.2,
                                        cfChIP_motif_lengths_NSC.3,
                                        cfChIP_motif_lengths_NSC.4,
                                        cfChIP_motif_lengths_SSC.1,
                                        cfChIP_motif_lengths_SSC.2,
                                        cfChIP_motif_lengths_SSC.3,
                                        cfChIP_motif_lengths_SSC.4) %>% 
    mutate(patient_ID = c(rep("NAC.1", nrow(cfChIP_motif_lengths_NAC.1)),
                      rep("NAC.2", nrow(cfChIP_motif_lengths_NAC.2)),
                      rep("NAC.3", nrow(cfChIP_motif_lengths_NAC.3)),
                      rep("NAC.4", nrow(cfChIP_motif_lengths_NAC.4)),
                      rep("NSC.1", nrow(cfChIP_motif_lengths_NSC.1)),
                      rep("NSC.2", nrow(cfChIP_motif_lengths_NSC.2)),
                      rep("NSC.3", nrow(cfChIP_motif_lengths_NSC.3)),
                      rep("NSC.4", nrow(cfChIP_motif_lengths_NSC.4)),
                      rep("SSC.1", nrow(cfChIP_motif_lengths_SSC.1)),
                      rep("SSC.2", nrow(cfChIP_motif_lengths_SSC.2)),
                      rep("SSC.3", nrow(cfChIP_motif_lengths_SSC.3)),
                      rep("SSC.4", nrow(cfChIP_motif_lengths_SSC.4))))

length_fragment_with_motif_plot(collected_cfChIP_motif_lengths)

gc()
length_fragment_with_motif_boxplot(collected_cfChIP_motif_lengths)
setwd("C:/Users/chris/OneDrive/1PhD/Fragmentering/endemotiver")
fragment_length_input_cfChIP_mutated_genes <- length_of_mutated_genes_df("BL mutations input files.txt",
                                                                         named_cancer_input_list,
                                                                         named_cancer_cfChIP_list,grs)
length_of_mutated_genes_plot(fragment_length_input_cfChIP_mutated_genes,T)
length_of_mutated_genes_plot(fragment_length_input_cfChIP_mutated_genes,"individual")
length_of_mutated_genes_boxplot(fragment_length_input_cfChIP_mutated_genes)

cancer_healthy_di_nuc_fraction <- di_nucleosome_fraction_df(list(NAC.1_input=NAC.1_input,
                                                                 NAC.2_input=NAC.2_input,
                                                                 NAC.3_input=NAC.3_input,
                                                                 NAC.4_input=NAC.4_input,
                                                                 NSC.1_input=NSC.1_input,
                                                                 NSC.2_input=NSC.2_input,
                                                                 NSC.3_input=NSC.3_input,
                                                                 NSC.4_input=NSC.4_input,
                                                                 SSC.1_input=SSC.1_input,
                                                                 SSC.2_input=SSC.2_input,
                                                                 SSC.3_input=SSC.3_input,
                                                                 SSC.4_input=SSC.4_input),
                                                            list(HC.1_input=HC.1_input,
                                                                 HC.2_input=HC.2_input,
                                                                 HC.3_input=HC.3_input,
                                                                 HC.4_input=HC.4_input,
                                                                 HC.5_input=HC.5_input,
                                                                 HC.6_input=HC.6_input,
                                                                 HC.7_input=HC.7_input),
                                                            z = c("Cancer", "Healthy"))
di_nucleosome_fraction_boxplot(cancer_healthy_di_nuc_fraction,c("Healthy","Cancer"))

ct_pos_neg_healthy_di_nuc_fraction <- di_nucleosome_fraction_df(list(NAC.1_input = NAC.1_input,
                                                                     NAC.3_input = NAC.3_input, 
                                                                     NAC.4_input = NAC.4_input,
                                                                     NSC.1_input = NSC.1_input,
                                                                     NSC.2_input = NSC.2_input,
                                                                     SSC.1_input = SSC.1_input,
                                                                     SSC.2_input = SSC.2_input,
                                                                     SSC.3_input = SSC.3_input,
                                                                     SSC.4_input = SSC.4_input),
                                                                list(NAC.2_input = NAC.2_input,
                                                                     NSC.3_input = NSC.3_input, 
                                                                     NSC.4_input = NSC.4_input),
                                                                list(HC.1_input = HC.1_input,
                                                                     HC.2_input = HC.2_input,
                                                                     HC.3_input = HC.3_input,
                                                                     HC.4_input = HC.4_input,
                                                                     HC.5_input = HC.5_input,
                                                                     HC.6_input = HC.6_input,
                                                                     HC.7_input = HC.7_input),
                                                            z = c("ctDNA positive", "ctDNA negative","Healthy"))

di_nucleosome_fraction_boxplot(ct_pos_neg_healthy_di_nuc_fraction,c("Healthy", "ctDNA negative","ctDNA positive"))

input_cfChIP_di_nuc_fraction <- di_nucleosome_fraction_df(list(NAC.1_input=NAC.1_input,
                                                                 NAC.2_input=NAC.2_input,
                                                                 NAC.3_input=NAC.3_input,
                                                                 NAC.4_input=NAC.4_input,
                                                                 NSC.1_input=NSC.1_input,
                                                                 NSC.2_input=NSC.2_input,
                                                                 NSC.3_input=NSC.3_input,
                                                                 NSC.4_input=NSC.4_input,
                                                                 SSC.1_input=SSC.1_input,
                                                                 SSC.2_input=SSC.2_input,
                                                                 SSC.3_input=SSC.3_input,
                                                                 SSC.4_input=SSC.4_input),
                                                          list(NAC.1_cfChIP=NAC.1_cfChIP,
                                                               NAC.2_cfChIP=NAC.2_cfChIP,
                                                               NAC.3_cfChIP=NAC.3_cfChIP,
                                                               NAC.4_cfChIP=NAC.4_cfChIP,
                                                               NSC.1_cfChIP=NSC.1_cfChIP,
                                                               NSC.2_cfChIP=NSC.2_cfChIP,
                                                               NSC.3_cfChIP=NSC.3_cfChIP,
                                                               NSC.4_cfChIP=NSC.4_cfChIP,
                                                               SSC.1_cfChIP=SSC.1_cfChIP,
                                                               SSC.2_cfChIP=SSC.2_cfChIP,
                                                               SSC.3_cfChIP=SSC.3_cfChIP,
                                                               SSC.4_cfChIP=SSC.4_cfChIP),
                                                            z = c("Cancer input", "Cancer cfChIP"))
di_nucleosome_fraction_boxplot(input_cfChIP_di_nuc_fraction,c("Cancer input", "Cancer cfChIP"),T)

collected_fragment_length_active_inactive <- read.table("Fragment lengths active inactive.txt", header = T)

high_low_di_nuc_fraction <- di_nucleosome_high_low_fraction_df(collected_fragment_length_active_inactive)

di_nucleosome_fraction_boxplot(high_low_di_nuc_fraction,c("Cancer Low expression", "Cancer High expression"),T)

pippin_vs_cfChIP_df <- pippin_vs_cfChIP(enrichment_df,pippin_enrichment_df)

pippin_vs_cfChIP_df

pippin_vs_cfChIP_plot(pippin_vs_cfChIP_df)

pippin_vs_cfChIP_correlation(enrichment_df,pippin_enrichment_df)

cfChIP_pippin_correlations <- data.frame(sample = c("NAC.1","NAC.2","NAC.3","NAC.4",
                                                    "NSC.1","NSC.2","NSC.3",
                                                    "SSC.1","SSC.2","SSC.3","SSC.4"),
                                         rho = c(0.788,0.553,0.769,0.840,
                                                 0.678,0.756,0.843,
                                                 0.800,0.882,0.499,0.644))

pippin_vs_cfChIP_correlation_matrix(cfChIP_pippin_correlations)

fragment_length(list(NAC.1_pippin = NAC.1_pippin, 
                     NAC.2_pippin = NAC.2_pippin,
                     NAC.3_pippin = NAC.3_pippin, 
                     NAC.4_pippin = NAC.4_pippin,
                     NSC.1_pippin = NSC.1_pippin, 
                     NSC.2_pippin = NSC.2_pippin,
                     NSC.3_pippin = NSC.3_pippin, 
                     SSC.1_pippin = SSC.1_pippin, 
                     SSC.2_pippin = SSC.2_pippin,
                     SSC.3_pippin = SSC.3_pippin, 
                     SSC.4_pippin = SSC.4_pippin),
                list(NAC.1_input = NAC.1_input, 
                     NAC.2_input = NAC.2_input,
                     NAC.3_input = NAC.3_input, 
                     NAC.4_input = NAC.4_input,
                     NSC.1_input = NSC.1_input, 
                     NSC.2_input = NSC.2_input,
                     NSC.3_input = NSC.3_input, 
                     SSC.1_input = SSC.1_input, 
                     SSC.2_input = SSC.2_input,
                     SSC.3_input = SSC.3_input, 
                     SSC.4_input = SSC.4_input),
                c("Size selected", "Input"), "",F)

proportion_sub150_cancer_pippin <- proportion_sub_df(list(NAC.1_pippin, NAC.2_pippin,
                                                         NAC.3_pippin, NAC.4_pippin,
                                                         NSC.1_pippin, NSC.2_pippin,
                                                         NSC.3_pippin,
                                                         SSC.1_pippin, SSC.2_pippin,
                                                         SSC.3_pippin, SSC.4_pippin),
                                                    150,
                                                    c("NAC.1_pippin", "NAC.2_pippin",
                                                      "NAC.3_pippin", "NAC.4_pippin",
                                                      "NSC.1_pippin", "NSC.2_pippin",
                                                      "NSC.3_pippin",
                                                      "SSC.1_pippin", "SSC.2_pippin",
                                                      "SSC.3_pippin", "SSC.4_pippin"))

proportion_sub_boxplot(list(proportion_sub150_cancer_input %>% filter(sample != "NSC.4_input"),
                            proportion_sub150_cancer_pippin),
                       c("Input", "Size selected"),"Paired")

pippin_effect_MAF("pippin og input mutation.txt")

volcano_motif("Input", "Size selected", diff_gene_motiv_analysis(prop_mer_input_cancer %>% select(-NSC.4_input), prop_mer_pippin))


sequence_list_cancer_ish <- list(sequence_end_motif(dif_motif_prop(prop_mer_input_cancer, prop_mer_cfChIP_cancer), 0.05, "positive"),
                                 sequence_end_motif(dif_motif_prop(prop_mer_input_cancer, prop_mer_cfChIP_cancer), 0.05, "negative"),
                                 sequence_end_motif(diff_gene_motiv_analysis(prop_mer_input_cancer %>% dplyr::select(-NSC.4_input), prop_mer_pippin), 0.05, "positive"),
                                 sequence_end_motif(diff_gene_motiv_analysis(prop_mer_input_cancer %>% dplyr::select(-NSC.4_input), prop_mer_pippin), 0.05, "negative"))

motif_venn(sequence_list_cancer_ish,
           c("cfChIP", "Input",
             "Size selected", "Unselected"), 
           " ")

setwd("C:/Users/chris/OneDrive/1PhD/Fragmentering/endemotiver")
collected_in_silico_enrichment <- read.table("collected in silico enrichment.txt", header =T)
collected_in_silico_enrichment
collected_in_silico_fraction <- read.table("collected in silico fraction.txt", header =T)

in_silico_correlation_plot(NAC.1_in_silico_enrichment,enrichment_df,grs)

size_selected_cfChIP <- enrich_pvalue(64,12,3,12)
unselected_input <- enrich_pvalue(64,3,1,6)
size_selected_input <- enrich_pvalue(64,15,9,0)
unselected_cfChIP <- enrich_pvalue(64,24,7,0)
high_expression_cfChIP_cancer <- enrich_pvalue(64,10,7,14)
high_expression_cfChIP_healthy <- enrich_pvalue(64,4,13,2)
high_expression_cfChIP_all <- enrich_pvalue(64,10,10,14)

setwd("C:/Users/chris/OneDrive/1PhD/Fragmentering/endemotiver")
exon_Sub150_collected <- read.table("collected sub 150 in exons.txt", header =T)
sub150_exon_plot(exon_Sub150_collected,enrichment_df)

collected_fragment_endpoint_wt_ct <- fragment_endpoint_wt_ctdna_df("BL mutations input files.txt",
                                                                   named_cancer_input_list)
setwd("C:/Users/chris/OneDrive/1PhD/Fragmentering/endemotiver")

fragment_endpoint_wt_ctdna_plot(collected_fragment_endpoint_wt_ct,T)
fragment_endpoint_wt_ctdna_plot(collected_fragment_endpoint_wt_ct,F)
test_1 <- collected_fragment_endpoint_wt_ct %>% filter(group == "SSC.4_input") %>% filter(Sample == "WT") %>% 
    filter(-200 < relative & relative < 200)
hist(test_1$relative,breaks = 50)