
library(cfDNAPro)
library(scales)
library(ggpubr)
library(ggplot2)
library(dplyr)

?callSize
my_data_path <- "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/picard.jar files/textfiles for cfDNAPro"

####Healthy input####
Healthy_input_plot <- cfDNAPro::callSize(path = my_data_path) %>%
  dplyr::filter(group == as.character("Healthy input")) %>%
  cfDNAPro::plotSingleGroup(vline = c(150,167,2*167))
Healthy_input_plot$prop_plot
Healthy_input_plot$cdf_plot
Healthy_input_plot$one_minus_cdf_plot


####Cancer input#####
Cancer_input_plot <- cfDNAPro::callSize(path = my_data_path) %>%
  dplyr::filter(group == as.character("Cancer input")) %>%
  cfDNAPro::plotSingleGroup(vline = c(150,167,2*167))
Cancer_input_plot$prop_plot
Cancer_input_plot$cdf_plot
Cancer_input_plot$one_minus_cdf_plot


####Healthy cfChIP####
Healthy_cfChIP_plot <- cfDNAPro::callSize(path = my_data_path) %>%
  dplyr::filter(group == as.character("Healthy cfChIP")) %>%
  cfDNAPro::plotSingleGroup(vline = c(150,167,2*167))
Healthy_cfChIP_plot$prop_plot
Healthy_cfChIP_plot$cdf_plot
Healthy_cfChIP_plot$one_minus_cdf_plot

?plotAllToOne

####Cancer cfChIP####
Cancer_cfChIP_plot <- cfDNAPro::callSize(path = my_data_path) %>%
  dplyr::filter(group == as.character("Cancer cfChIP")) %>%
  cfDNAPro::plotSingleGroup(vline = c(150,167,2*167))
Cancer_cfChIP_plot$prop_plot
Cancer_cfChIP_plot$cdf_plot
Cancer_cfChIP_plot$one_minus_cdf_plot

grp_list<-list("Healthy_input"="Healthy input",
               "Cancer_input"="Cancer input",
               "Healthy_cfChIP"="Healthy cfChIP",
               "Cancer_cfChIP"="Cancer cfChIP")
res <- cfDNAPro::callSize(path = my_data_path) %>% cfDNAPro::plotAllToOne(vline = c(150,167,2*167))

res$prop_plot

result<-sapply(grp_list, function(x){
  result <-callSize(path = my_data_path) %>%
    dplyr::filter(group==as.character(x)) %>%
    plotSingleGroup(vline = c(150,167,2*167))
}, simplify = FALSE)

####Input: Healthy vs. Cancer####
suppressWarnings(
  multiplex <-
    ggarrange(result$Healthy_input$prop_plot +
                theme(axis.title.x = element_blank()),
              result$Cancer_input$prop_plot +
                theme(axis.title = element_blank()),
              result$Healthy_input$cdf_plot,
              result$Cancer_input$cdf_plot +
                theme(axis.title.y = element_blank()),
              labels = c("Healthy input (n=7)", "Cancer input (n=12)"),
              label.x = 0.2,
              ncol = 2,
              nrow = 2))

multiplex

####cfChIP: Healthy vs. Cancer####
suppressWarnings(
  multiplex <-
    ggarrange(result$Healthy_cfChIP$prop_plot +
                theme(axis.title.x = element_blank()),
              result$Cancer_cfChIP$prop_plot +
                theme(axis.title = element_blank()),
              result$Healthy_cfChIP$cdf_plot,
              result$Cancer_cfChIP$cdf_plot +
                theme(axis.title.y = element_blank()),
              labels = c("Healthy cfChIP (n=4)", "Cancer cfChIP (n=12)"),
              label.x = 0.2,
              ncol = 2,
              nrow = 2))

multiplex


####Healthy: Input vs. cfChIP####
suppressWarnings(
  multiplex <-
    ggarrange(result$Healthy_input$prop_plot +
                theme(axis.title.x = element_blank()),
              result$Healthy_cfChIP$prop_plot +
                theme(axis.title = element_blank()),
              result$Healthy_input$cdf_plot,
              result$Healthy_cfChIP$cdf_plot +
                theme(axis.title.y = element_blank()),
              labels = c("Healthy input (n=7)", "Healthy cfChIP (n=4)"),
              label.x = 0.2,
              ncol = 2,
              nrow = 2))

multiplex


####Cancer: Input vs. cfChIP####
suppressWarnings(
  multiplex <-
    ggarrange(result$Cancer_input$prop_plot +
                theme(axis.title.x = element_blank()),
              result$Cancer_cfChIP$prop_plot +
                theme(axis.title = element_blank()),
              result$Cancer_input$cdf_plot,
              result$Cancer_cfChIP$cdf_plot +
                theme(axis.title.y = element_blank()),
              labels = c("Cancer input (n=12)", "Cancer cfChIP (n=12)"),
              label.x = 0.2,
              ncol = 2,
              nrow = 2))

multiplex



####Combining samples in groups, and plotting the median proportion of reads####
order <- c("Healthy input",
           "Cancer input",
           "Healthy cfChIP",
           "Cancer cfChIP")

compare_grps<-callMetrics(my_data_path) %>% plotMetrics(order=order, vline = c(150,167,167*2))


p1<-compare_grps$median_prop_plot +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12,face="bold")) +
  theme(legend.position = c(0.7, 0.5),
        legend.text = element_text( size = 11),
        legend.title = element_blank())
p1

p2<-compare_grps$median_cdf_plot +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  theme(axis.title=element_text(size=12,face="bold")) +
  theme(legend.position = c(0.7, 0.5),
        legend.text = element_text( size = 11),
        legend.title = element_blank())
p2

suppressWarnings(
  median_grps<-ggpubr::ggarrange(p1,
                                 p2,
                                 label.x = 0.3,
                                 ncol = 1,
                                 nrow = 2
  ))


median_grps


####Modal fragment size####
mode_bin <- callMode(my_data_path) %>% plotMode(order=order,hline = c(167,150,133))
suppressWarnings(print(mode_bin))

#Mode as stacked bar chart

mode_stacked <-
  callMode(my_data_path) %>%
  plotModeSummary(order=order,
                  mode_partition = list(c(165,166,167)))
mode_stacked <- mode_stacked + theme(legend.position = "top")
suppressWarnings(print(mode_stacked))
#Synes ikke dette plot virker særlig godt. Det viser egentlig antallet af samples
#som har den mode man skriver i mode_partition, eller som har en anden mode



####Inter-peak/ valley Distance####

inter_peak_dist<-callPeakDistance(path = my_data_path,  limit = c(50, 135)) %>%
  plotPeakDistance(order = order) +
  labs(y="Fraction") +
  theme(axis.title =  element_text(size=12,face="bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.5),
        legend.text = element_text(size = 11))

suppressWarnings(print(inter_peak_dist))


inter_valley_dist<-callValleyDistance(path = my_data_path,
                                     limit = c(50, 135)) %>%
  plotValleyDistance(order = order) +
  labs(y="Fraction") +
  theme(axis.title =  element_text(size=12,face="bold"),
        legend.title = element_blank(),
        legend.position = c(0.91, 0.5),
        legend.text = element_text(size = 11))

suppressWarnings(print(inter_valley_dist))



####Random####
peaks <- callPeakDistance(path = one_by_one)
valleys <- callValleyDistance(path = one_by_one)
plot_all <- callSize(path=one_by_one) %>% plotSingleGroup(vline = NULL)
file_name <- valleys$file_name[1]
str <- strsplit(file_name, split = "/")[1]
fil <- strsplit(str[[1]][10],split = ".", fixed = T)
tit <- fil[[1]][1]
plot_all$prop +
  coord_cartesian(xlim = c(90,135),ylim = c(0,0.0185)) +
  geom_vline(xintercept=peaks$insert_size, colour="red",linetype="dashed") +
  geom_vline(xintercept = valleys$insert_size,colour="blue")+
  ggtitle(paste(tit))
one_by_one <- "C:/Users/Christoffer/OneDrive/1PhD/Fragmentering/picard.jar files/one by one"

plot_all$prop_plot
sessionInfo()

