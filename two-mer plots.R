library(ggpubr)
####Fig 1 ####
fig1a <- fragment_length(list(NAC.1_input = NAC.1_input, 
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
                              SSC.4_input = SSC.4_input),
                         list(HC.1_input = HC.1_input,
                              HC.2_input = HC.2_input,
                              HC.3_input = HC.3_input,
                              HC.4_input = HC.4_input,
                              HC.5_input = HC.5_input,
                              HC.6_input = HC.6_input,
                              HC.7_input = HC.7_input),
                         c("Cancer", "Healthy"), "",F)

saveRDS(fig1a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1A.rds",compress = F)
fig1a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1A.rds")
fig1a
gc()
fig1b <- .Last.value
fig1b <- fragment_length(list(NAC.1_input = NAC.1_input,
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
                         c("ctDNA positive", "ctDNA negative"), "",F)

saveRDS(fig1b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1B.rds",compress = F)
fig1b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1B.rds")
fig1b
gc()

fig1c <- proportion_sub_boxplot(list(proportion_sub150_healthy_input,
                            proportion_sub150_ct_neg_input,
                            proportion_sub150_ct_pos_input),
                       c("Healthy",
                         "ctDNA negative",
                         "ctDNA positive"),
                       "Unpaired") 

saveRDS(fig1c, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1C.rds",compress = F)
fig1c <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1/fig1C.rds")
fig1c <- fig1c + theme(legend.position = "none")
fig1c

fig1d <- fragment_length_wt_ctdna_plot(fragment_length_wt_ctdna_input,T)

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
fig1f <- fig1f+ theme(legend.position = "none")


fig1_list <- list(fig1a,fig1b,fig1c,fig1d,fig1e,fig1f)
gc()
ggarrange(plotlist = fig1_list,
          nrow = 2,
          ncol = ceiling(length(fig1_list)/2),
          widths = c(2,2,1),
          labels = "AUTO")
ggsave(filename = "fig1.png",
       width = 27000, height = 13000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 1",
       dpi = 1200,
       device = "png")
####Fig 2 ####

fig2a <- fragment_length(list(NAC.1_cfChIP = NAC.1_cfChIP, 
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
                              SSC.4_cfChIP = SSC.4_cfChIP),
                         list(NAC.1_input = NAC.1_input, 
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
                              SSC.4_input = SSC.4_input),
                         c("Cancer cfChIP", "Cancer Input"),"",T)
gc()
saveRDS(fig2a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2A.rds",compress = F)
fig2a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2A.rds")
fig2a
    

fig2b <- proportion_sub_boxplot(list(proportion_sub150_cancer_input,
                                     proportion_sub_150_cancer_cfChIP),
                                c("Cancer Input", "Cancer cfChIP"),"Paired")

saveRDS(fig2b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2B.rds",compress = F)
fig2b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2B.rds")
fig2b <- fig2b + theme(legend.position = "none")

fig2c <- fragment_length_active_inactive_plot(collected_fragment_length_active_inactive,
                                              "Cancer",T)
saveRDS(fig2c, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2C.rds",compress = F)
fig2c <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2C.rds")
fig2c


fig2d <- proportion_sub_boxplot(list(fragment_length_inactive_cancer,fragment_length_active_cancer),
                       c("Cancer Low epxressed", "Cancer High expressed"),
                       "Paired", "activity")

saveRDS(fig2d, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2D.rds",compress = F)
fig2d <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2D.rds")
fig2d <- fig2d + theme(legend.position = "none")


fig2e <- sub150_fraction_quartiles_plot(collected_sub150_fraction_quartiles)

saveRDS(fig2e, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2E.rds",compress = F)
fig2e <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2E.rds")
fig2e

fig2f <- ggplot(data = data.frame(x=seq(1:20),y = seq(21:40)))+
    geom_point(aes(x= x,y=y))+
    theme_bw()+
    labs(title = "Placeholder plot")+
    th

fig2g <- distance_nucleosome_plot(collected_nuc_dist_df_bed, TRUE)

saveRDS(fig2g, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2G.rds",compress = F)
fig2g <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2G.rds")
fig2g

fig2h <- distance_nucleosome_boxplot(collected_nuc_dist_df_bed)

saveRDS(fig2h, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2H.rds",compress = F)
fig2h <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 2/fig2H.rds")
fig2h <- fig2h + theme(legend.position = "none")

fig2_list <- list(fig2a,fig2b,fig2c,fig2d,fig2e,fig2f,fig2g,fig2h)
gc()
ggarrange(plotlist = fig2_list,
          nrow = 4,
          ncol = ceiling(length(fig2_list)/4),
          labels = "AUTO",
          widths = c(2,1)) + bgcolor("White")

ggsave(filename = "fig2.png",
       width = 20000, height = 25000, units = "px",
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

fig3d <- volcano_motif("Low expressed", "High expressed", diff_active_inactive_end_motif(cancer_active_inactive))

saveRDS(fig3d, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 3/fig3D.rds",compress = F)
fig3d <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 3/fig3D.rds")
fig3d

fig3e <- umap_moitf(prop_mer_active, prop_mer_inactive, "High expressed", "Low expressed", 3, 32, "Activity")

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

fig4c <- length_fragment_with_motif_plot(collected_cfChIP_motif_lengths)

saveRDS(fig4c, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 4/fig4C.rds",compress = F)
fig4c <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 4/fig4C.rds")
fig4c

fig4d <- length_fragment_with_motif_boxplot(collected_cfChIP_motif_lengths)

saveRDS(fig4d, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 4/fig4D.rds",compress = F)
fig4d <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 4/fig4D.rds")
fig4d <- fig4d + theme(legend.position = "none")

fig4e <- sub150_active_motif_fraction_quartiles_plot(collected_sub150_active_motif_fraction_quartiles, "Active")

saveRDS(fig4e, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 4/fig4E.rds",compress = F)
fig4e <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 4/fig4E.rds")
fig4e

fig4_list <- list(fig4a, fig4b, fig4c, fig4d, fig4e)
fig4_list1 <- list(fig4a,fig4b)
fig4_list2 <- list(fig4c,fig4d)
fig4_list3 <- list(fig4e, NULL)
gc()
ggarrange(plotlist = list(ggarrange(plotlist = fig4_list1, 
                                    ncol = ceiling(length(fig4_list1)/1),
                                    labels = c("A","B")),
          ggarrange(plotlist = fig4_list2,widths = c(2,1),
                    ncol = ceiling(length(fig4_list2)/1),
                    labels = c("C","D")),
          ggarrange(fig4e,NULL,ncol = ceiling(length(fig4_list3)/1),
                    labels = c("E",""))),
          nrow = 3,
          ncol = ceiling(length(fig4_list)/6) + bgcolor("White"))
ggarrange(fig4_list,
          nrow = 3,
          ncol = ceiling(length(fig4_list)/3),
          labels = c("AUTO") + bgcolor("White"))
ggsave(filename = "fig4.png",
       width = 15000, height = 20000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Fig. 4",
       dpi = 1200,
       device = "png")
#### Supplementary figures ####

####Supplementary figure 1####

s1a <- fragment_length(list(NAC.1_cfChIP = NAC.1_cfChIP, 
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
                            SSC.4_cfChIP = SSC.4_cfChIP),
                       list(HC.1_cfChIP = HC.1_cfChIP, 
                            HC.2_cfChIP = HC.2_cfChIP,
                            HC.3_cfChIP = HC.3_cfChIP, 
                            HC.4_cfChIP = HC.4_cfChIP),
                c("Cancer cfChIP", "Healthy cfChIP"),"",F)
s1a
saveRDS(s1a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S1/S1A.rds",compress = F)
s1a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S1/S1A.rds")
s1a


s1b <- proportion_sub_boxplot(list(proportion_sub150_healthy_cfChIP,
                            proportion_sub_150_cancer_cfChIP),
                       c("Healthy cfChIP", "Cancer cfChIP"),
                       "Unpaired")

saveRDS(s1b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S1/S1B.rds",compress = F)
s1b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S1/S1B.rds")
s1b <- s1b + theme(legend.position = "none")

s1_list <- list(s1a,s1b)
gc()
ggarrange(plotlist = s1_list,
          nrow = 1,
          ncol = 2,widths = c(2,1),
          labels = "AUTO") + bgcolor("White")  
ggsave(filename = "S1.png",
       width = 20000, height = 8000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S1",
       dpi = 1200,
       device = "png")

####Supplementary figure 2####
s2a <- fragment_length_wt_ctdna_plot(fragment_length_wt_ctdna_input,F)

saveRDS(s2a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S2/S2A.rds",compress = F)
s2a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S2/S2A.rds")
s2a

s2_list <- list(s2a)

ggarrange(plotlist = s2_list,
          nrow = 1,
          ncol = 1,
          labels = "AUTO") + bgcolor("White") 

ggsave(filename = "S2.png",
       width = 15000, height = 8000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S2",
       dpi = 1200,
       device = "png")

####Supplementary figure 3####

s3a <- fragment_length(list(NAC.1_cfChIP = NAC.1_cfChIP, 
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
                                     SSC.4_cfChIP = SSC.4_cfChIP),
                                list(NAC.1_input = NAC.1_input, 
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
                                     SSC.4_input = SSC.4_input),
                                c("Cancer cfChIP", "Cancer Input"),"","individual")

saveRDS(s3a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S3/S3A.rds",compress = F)
s3a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S3/S3A.rds")
s3a

s3b <- fragment_length_active_inactive_plot(collected_fragment_length_active_inactive,
                                            "Cancer","individual")

saveRDS(s3b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S3/S3B.rds",compress = F)
s3b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S3/S3B.rds")
s3b

gc()

s3_list <- list(s3a,s3b)

ggarrange(plotlist = s3_list,
          nrow = 2,
          ncol = ceiling(length(s3_list)/2),
          labels = "AUTO") + bgcolor("White") 

ggsave(filename = "S3.png",
       width = 25000, height = 20000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S3",
       dpi = 1200,
       device = "png")

####Supplementary figure 4####

s4a <- fragment_length(list(HC.1_cfChIP = HC.1_cfChIP,
                            HC.2_cfChIP = HC.2_cfChIP,
                            HC.3_cfChIP = HC.3_cfChIP,
                            HC.4_cfChIP = HC.4_cfChIP),
                list(HC.1_input = HC.1_input,
                     HC.2_input = HC.2_input,
                     HC.3_input = HC.3_input,
                     HC.4_input = HC.4_input),
                c("Healthy cfChIP", "Healthy Input"), "",T)

saveRDS(s4a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S4/S4A.rds",compress = F)
s4a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S4/S4A.rds")
s4a

s4b <- fragment_length(list(HC.1_cfChIP = HC.1_cfChIP,
                            HC.2_cfChIP = HC.2_cfChIP,
                            HC.3_cfChIP = HC.3_cfChIP,
                            HC.4_cfChIP = HC.4_cfChIP),
                       list(HC.1_input = HC.1_input,
                            HC.2_input = HC.2_input,
                            HC.3_input = HC.3_input,
                            HC.4_input = HC.4_input),
                       c("Healthy cfChIP", "Healthy Input"), "","individual")

s4b
saveRDS(s4b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S4/S4B.rds",compress = F)
s4b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S4/S4B.rds")
s4b

s4c <- proportion_sub_boxplot(list(proportion_sub150_healthy_input_w_ChIP,
                            proportion_sub150_healthy_cfChIP),
                       c("Healthy Input", "Healthy cfChIP"),
                       "Paired")

saveRDS(s4c, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S4/S4C.rds",compress = F)
s4c <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S4/S4C.rds")
s4c <- s4c + theme(legend.position = "none")
s4c

s4d <- fragment_length_active_inactive_plot(collected_fragment_length_active_inactive_healthy,
                                            "Healthy",T)

saveRDS(s4d, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S4/S4D.rds",compress = F)
s4d <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S4/S4D.rds")
s4d

s4e <- fragment_length_active_inactive_plot(collected_fragment_length_active_inactive_healthy,
                                            "Healthy","individual")

saveRDS(s4e, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S4/S4E.rds",compress = F)
s4e <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S4/S4E.rds")
s4e



s4f <- proportion_sub_boxplot(list(fragment_length_inactive_healthy,fragment_length_active_healthy),
                              c("Healthy Inactive", "Healthy Active"),
                              "Paired", "activity")

saveRDS(s4f, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S4/S4F.rds",compress = F)
s4f <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S4/S4F.rds")
s4f

s4_list <- list(s4a,s4b,s4c,s4d,s4e,s4f)
gc()
ggarrange(plotlist = s4_list,
          nrow = 2,
          ncol = ceiling(length(s4_list)/2),
          widths = c(2,2,1),
          labels = "AUTO") + bgcolor("White")  
ggsave(filename = "S4.png",
       width = 32000, height = 13000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S4",
       dpi = 1200,
       device = "png")

####Supplementary figure 5 ####

s5a <- di_nucleosome_fraction_boxplot(ct_pos_neg_healthy_di_nuc_fraction,c("Healthy", "ctDNA negative","ctDNA positive"))

saveRDS(s5a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S5/S5A.rds",compress = F)
s5a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S5/S5A.rds")
s5a

s5b <- di_nucleosome_fraction_boxplot(input_cfChIP_di_nuc_fraction,c("Cancer input", "Cancer cfChIP"),T)

saveRDS(s5b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S5/S5B.rds",compress = F)
s5b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S5/S5B.rds")
s5b

s5c <- di_nucleosome_fraction_boxplot(high_low_di_nuc_fraction,c("Cancer Low expression", "Cancer High expression"),T)

saveRDS(s5c, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S5/S5C.rds",compress = F)
s5c <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S5/S5C.rds")
s5c

s5_list <- list(s5a,s5b,s5c)

ggarrange(plotlist = s5_list,
          nrow = 1,
          ncol = ceiling(length(s5_list)/1),
          labels = "AUTO") + bgcolor("White") 

ggsave(filename = "S5.png",
       width = 17000, height = 7000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S5",
       dpi = 1200,
       device = "png")

####Supplementary figure 6 ####

s6a <- volcano_motif("Input", "cfChIP", dif_motif_prop(prop_mer_input_healthy, prop_mer_cfChIP_healthy))

saveRDS(s6a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S6/S6A.rds",compress = F)
s6a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S6/S6A.rds")
s6a


s6b <- volcano_motif("Input", "cfChIP", dif_motif_prop(prop_mer_input, prop_mer_cfChIP))

saveRDS(s6b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S6/S6B.rds",compress = F)
s6b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S6/S6B.rds")
s6b


s6c <- volcano_motif("Low expressed", "High expressed", diff_active_inactive_end_motif(healthy_active_inactive))

saveRDS(s6c, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S6/S6C.rds",compress = F)
s6c <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S6/S6C.rds")
s6c


s6d <- volcano_motif("Low expressed", "High expressed", diff_active_inactive_end_motif(collected_active_inactive_df))

saveRDS(s6d, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S6/S6D.rds",compress = F)
s6d <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S6/S6D.rds")
s6d

s6_list <- list(s6a,s6b,s6c,s6d)

gc()
ggarrange(plotlist = s6_list,
          nrow = 2,
          ncol = ceiling(length(s6_list)/2),
          labels = "AUTO") + bgcolor("White")

ggsave(filename = "S6.png",
       width = 20000, height = 13000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S6",
       dpi = 1200,
       device = "png")


####Supplementary figure 7 ####

s7a <- sequence_end_motif_plot(sequence_list_cancer, 
                               c("Cancer cfChIP", "Cancer Input",
                                "Cancer High expressed", "Cancer Low expressed"))

saveRDS(s7a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/8/S7A.rds",compress = F)
s7a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S7/S7A.rds")
s7a


s7b <- motif_venn(sequence_list_cancer,
                  c("cfChIP", "Input",
                    "High", "Low"), 
                  "Cancer samples")

saveRDS(s7b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S7/S7B.rds",compress = F)
s7b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S7/S7B.rds")
s7b


s7c <- sequence_end_motif_plot(sequence_list_healthy, 
                               c("Healthy cfChIP", "Healthy Input",
                                 "Healthy High expressed", "Healthy Low expressed"))

saveRDS(s7c, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S7/S7C.rds",compress = F)
s7c <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S7/S7C.rds")
s7c


s7d <- motif_venn(sequence_list_healthy,
                  c("cfChIP", "Input",
                    "High", "Low"), 
                  "Healthy samples")

saveRDS(s7d, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S7/S7D.rds",compress = F)
s7d <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S7/S7D.rds")
s7d


s7e <- sequence_end_motif_plot(sequence_list_all_samples, 
                               c("cfChIP", "Input",
                                 "High expressed", "Low expressed"))

saveRDS(s7e, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S7/S7E.rds",compress = F)
s7e <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S7/S7E.rds")
s7e

s7f <- motif_venn(sequence_list_all_samples,
                  c("cfChIP", "Input",
                    "High", "Low"), 
                  "All samples")

saveRDS(s7f, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/s7/s7F.rds",compress = F)
s7f <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/s7/s7F.rds")
s7f

s7_list <- list(s7a,s7b,s7c,s7d,s7e,s7f)

ggarrange(plotlist = s7_list,
          nrow = 3,
          ncol = ceiling(length(s7_list)/3),
          labels = "AUTO") + bgcolor("White")

ggsave(filename = "S7.png",
       width = 20000, height = 25000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S7",
       dpi = 1200,
       device = "png")

####Supplementary figure 8####
s8a <- high_motif_fraction_quartiles_plot(collected_active_active_cancer_fraction_quartiles, "Active")

saveRDS(s8a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S8/S8A.rds",compress = F)
s8a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S8/S8A.rds")
s8a

s8b <- high_motif_fraction_quartiles_plot(collected_inactive_inactive_cancer_fraction_quartiles, "Inactive")

saveRDS(s8b, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S8/S8B.rds",compress = F)
s8b <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S8/S8B.rds")
s8b

s8_list <- list(s8a,s8b)

ggarrange(plotlist = s8_list,
          nrow = 2,
          ncol = ceiling(length(s8_list)/2),
          labels = "AUTO") + bgcolor("White")

ggsave(filename = "S8.png",
       width = 10000, height = 13000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S8",
       dpi = 1200,
       device = "png")


####Supplementary figure 9####

s9a <- length_fragment_with_motif_plot_zoom(collected_cfChIP_motif_lengths)

saveRDS(s9a, file = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S9/S9A.rds",compress = F)
s9a <- readRDS("C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S9/S9A.rds")
s9a

s9_list <- list(s9a)

ggarrange(plotlist = s9_list,
          nrow = 1,
          ncol = ceiling(length(s9_list)/1),
          labels = "AUTO") + bgcolor("White")

ggsave(filename = "S9.png",
       width = 15000, height = 10000, units = "px",
       path = "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/endemotiver/plots/Supplementary/S9",
       dpi = 1200,
       device = "png")
