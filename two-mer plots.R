####Fig 1 ####
fig1a <- fragment_length(list(NAC.1_input, NAC.2_input,
                              NAC.3_input, NAC.4_input,
                              NSC.1_input, NSC.2_input,
                              NSC.3_input, NSC.4_input,
                              SSC.1_input, SSC.2_input,
                              SSC.3_input, SSC.4_input),
                         list(HC.1_input, HC.2_input,
                              HC.3_input, HC.4_input),
                         c("Cancer", "Healthy"), "")

saveRDS(fig1a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1A.rds",compress = F)
fig1a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1A.rds")
fig1a


fig1b <- fragment_length(list(NAC.1_input,
                              NAC.3_input, NAC.4_input,
                              NSC.1_input, NSC.2_input,
                              SSC.1_input, SSC.2_input,
                              SSC.3_input, SSC.4_input),
                         list(NAC.2_input,NSC.3_input, NSC.4_input),
                         c("ctDNA positive", "ctDNA negative"), "")

saveRDS(fig1b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1B.rds",compress = F)
fig1b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1B.rds")
fig1b


fig1c <- proportion_sub_boxplot(list(proportion_sub150_healthy_input,
                            proportion_sub150_ct_neg_input,
                            proportion_sub150_ct_pos_input),
                       c("Healthy",
                         "ctDNA negative",
                         "ctDNA positive"),
                       "Unpaired")

saveRDS(fig1c, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1C.rds",compress = F)
fig1c <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1C.rds")
fig1c


fig1d <- fragment_length_wt_ctdna_plot(fragment_length_wt_ctdna_input)

saveRDS(fig1d, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1D.rds",compress = F)
fig1d <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1D.rds")
fig1d

fig1e <- MAF_in_bins(fragment_length_wt_ctdna_input)

saveRDS(fig1e, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1E.rds",compress = F)
fig1e <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1E.rds")
fig1e

fig1f <- fragment_length_wt_ctdna_boxplot(ct_wt_input_list)

saveRDS(fig1f, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1F.rds",compress = F)
fig1f <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1F.rds")
fig1f

fig1_list <- list(fig1a,fig1b,fig1c,fig1d,fig1e,fig1f)
gc()
ggarrange(plotlist = fig1_list,
          nrow = 2,
          ncol = ceiling(length(fig1_list)/2),
          labels = "AUTO")
?ggarrange
ggsave(filename = "fig1.png",
       width = 27000, height = 13000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1",
       dpi = 1200,
       device = "png")
####Fig 2 ####

fig2a <- fragment_length(list(NAC.1_cfChIP, NAC.2_cfChIP,
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
                         c("Cancer cfChIP", "Cancer Input"),"")
gc()
saveRDS(fig2a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2A.rds",compress = F)
fig2a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2A.rds")
fig2a
    

fig2b <- proportion_sub_boxplot(list(proportion_sub150_cancer_input,
                                     proportion_sub_150_cancer_cfChIP),
                                c("Cancer Input", "Cancer cfChIP"),
                                "Paired")

saveRDS(fig2b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2B.rds",compress = F)
fig2b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2B.rds")
fig2b

fig2c <- fragment_length_active_inactive_plot(collected_fragment_length_active_inactive,
                                              "Cancer")
saveRDS(fig2c, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2C.rds",compress = F)
fig2c <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2C.rds")
fig2c


fig2d <- proportion_sub_boxplot(list(fragment_length_inactive_cancer,fragment_length_active_cancer),
                       c("Cancer Inactive", "Cancer Active"),
                       "Paired", "activity")

saveRDS(fig2d, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2D.rds",compress = F)
fig2d <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2D.rds")
fig2d


fig2e <- sub150_fraction_quartiles_plot(collected_sub150_fraction_quartiles)

saveRDS(fig2e, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2E.rds",compress = F)
fig2e <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2E.rds")
fig2e

#fig2f <- 

fig2_list <- list(fig2a,fig2b,fig2c,fig2d,fig2e)
gc()
ggarrange(plotlist = fig2_list,
          nrow = 3,
          ncol = ceiling(length(fig2_list)/3),
          labels = "AUTO") + bgcolor("White")

ggsave(filename = "fig2.png",
       width = 13000, height = 27000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2",
       dpi = 1200,
       device = "png")

####Fig 3 ####

fig3a <- volcano_motif("Healthy", "Cancer", diff_gene_motiv_analysis(prop_mer_input_healthy, prop_mer_input_cancer))

saveRDS(fig3a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 3/fig3A.rds",compress = F)
fig3a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 3/fig3A.rds")
fig3a


fig3b <- volcano_motif("Input", "cfChIP", dif_motif_prop(prop_mer_input_cancer, prop_mer_cfChIP_cancer))

saveRDS(fig3b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 3/fig3B.rds",compress = F)
fig3b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 3/fig3B.rds")
fig3b


fig3c <- umap_moitf(prop_mer_input, prop_mer_cfChIP, "input", "cfChIP", 8, 15)

saveRDS(fig3c, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 3/fig3C.rds",compress = F)
fig3c <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 3/fig3C.rds")
fig3c

fig3d <- volcano_motif("Inactive genes", "Active genes", diff_active_inactive_end_motif(cancer_active_inactive))

saveRDS(fig3d, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 3/fig3D.rds",compress = F)
fig3d <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 3/fig3D.rds")
fig3d

fig3e <- umap_moitf(prop_mer_active, prop_mer_inactive, "Active", "Inactive", 3, 32, "Activity")

saveRDS(fig3e, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 3/fig3E.rds",compress = F)
fig3e <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 3/fig3E.rds")
fig3e


fig3_list <- list(fig3a,fig3b,fig3c,fig3d,fig3e)
gc()
ggarrange(plotlist = fig3_list,
          nrow = 2,
          ncol = ceiling(length(fig3_list)/2),
          labels = "AUTO") + bgcolor("White")  
ggsave(filename = "fig3.png",
       width = 27000, height = 13000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 3",
       dpi = 1200,
       device = "png")
####Fig 4 ####

fig4a <- active_motif_fraction_quartiles_plot(collected_active_motif_fraction_quartiles, "Active")

saveRDS(fig4a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 4/fig4A.rds",compress = F)
fig4a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 4/fig4A.rds")
fig4a


fig4b <- active_motif_fraction_quartiles_plot(collected_inactive_motif_fraction_quartiles, "Inactive")

saveRDS(fig4b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 4/fig4B.rds",compress = F)
fig4b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 4/fig4B.rds")
fig4b

fig4c <- sub150_active_motif_fraction_quartiles_plot(collected_sub150_active_motif_fraction_quartiles, "Active")

saveRDS(fig4c, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 4/fig4C.rds",compress = F)
fig4c <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 4/fig4C.rds")
fig4c

fig4_list <- list(fig4a,fig4b,fig4c)
gc()
ggarrange(plotlist = fig4_list,
          nrow = 2,
          ncol = ceiling(length(fig4_list)/2),
          labels = "AUTO") + bgcolor("White")  
ggsave(filename = "fig4.png",
       width = 20000, height = 13000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 4",
       dpi = 1200,
       device = "png")
#### Supplementary figures ####

####Supplementary figure 1####

s1a <- fragment_length(list(NAC.1_cfChIP, NAC.2_cfChIP,
                     NAC.3_cfChIP, NAC.4_cfChIP,
                     NSC.1_cfChIP, NSC.2_cfChIP,
                     NSC.3_cfChIP, NSC.4_cfChIP,
                     SSC.1_cfChIP, SSC.2_cfChIP,
                     SSC.3_cfChIP, SSC.4_cfChIP),
                list(HC.1_cfChIP, HC.2_cfChIP,
                     HC.3_cfChIP, HC.4_cfChIP),
                c("Cancer cfChIP", "Healthy cfChIP"),"")

saveRDS(s1a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S1/S1A.rds",compress = F)
s1a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S1/S1A.rds")
s1a


s1b <- proportion_sub_boxplot(list(proportion_sub150_healthy_cfChIP,
                            proportion_sub_150_cancer_cfChIP),
                       c("Healthy cfChIP", "Cancer cfChIP"),
                       "Unpaired")

saveRDS(s1b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S1/S1B.rds",compress = F)
s1b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S1/S1B.rds")
s1b

s1_list <- list(s1a,s1b)
gc()
ggarrange(plotlist = s1_list,
          nrow = 1,
          ncol = 2,
          labels = "AUTO") + bgcolor("White")  
ggsave(filename = "S1.png",
       width = 20000, height = 8000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S1",
       dpi = 1200,
       device = "png")

####Supplementary figure 2####

s2a <- fragment_length(list(HC.1_cfChIP, HC.2_cfChIP,
                     HC.3_cfChIP, HC.4_cfChIP),
                list(HC.1_input, HC.2_input,
                     HC.3_input, HC.4_input),
                c("Healthy cfChIP", "Healthy Input"), "")

saveRDS(s2a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S2/S2A.rds",compress = F)
s2a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S2/S2A.rds")
s2a


s2b <- proportion_sub_boxplot(list(proportion_sub150_healthy_input_w_ChIP,
                            proportion_sub150_healthy_cfChIP),
                       c("Healthy Input", "Healthy cfChIP"),
                       "Paired")

saveRDS(s2b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S2/S2B.rds",compress = F)
s2b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S2/S2B.rds")
s2b


s2c <- fragment_length_active_inactive_plot(collected_fragment_length_active_inactive_healthy,
                                            "Healthy")

saveRDS(s2c, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S2/S2C.rds",compress = F)
s2c <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S2/S2C.rds")
s2c

s2d <- proportion_sub_boxplot(list(fragment_length_inactive_healthy,fragment_length_active_healthy),
                              c("Healthy Inactive", "Healthy Active"),
                              "Paired", "activity")

saveRDS(s2d, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S2/S2D.rds",compress = F)
s2d <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S2/S2D.rds")
s2d

S2_list <- list(s2a,s2b,s2c,s2d)
gc()
ggarrange(plotlist = S2_list,
          nrow = 2,
          ncol = ceiling(length(S2_list)/2),
          labels = "AUTO") + bgcolor("White")  
ggsave(filename = "S2.png",
       width = 20000, height = 13000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S2",
       dpi = 1200,
       device = "png")


