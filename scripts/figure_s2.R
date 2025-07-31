#---------------------------------------------------------------------
# Data import and workspace preparation
#---------------------------------------------------------------------

#Import the packages
library(ggplot2)
library(ggpubr)
library(mgcv)

#Import the data
source("scripts/data_import.R")



#Subset to only data from patients with paired blood and ff samples
paired_data <- data[data$ID %in% total_mps_merged$ID,]


#Generate the figure

FIG <-
  ggplot(data = paired_data,
         aes(x = sample,
             y = concentration,
             fill = sample)) +
  geom_boxplot(size = 0.1, outlier.size = 0.2, outlier.shape = 16, outlier.alpha = 0, alpha = 0.7) +
    facet_wrap(~category, scales = "free_y") +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2, color="black", fill="black") +
  geom_jitter(aes(col = sample), size = 1, shape = 16, position=position_jitter(height=0, width=0)) +
    geom_line(aes(group=patient), linewidth = 0.3, col = "grey20", alpha = 0.6) +
  scale_fill_manual(values = c("#edae49", "#00798c")) +
  scale_colour_manual(values = c("#edae49", "#00798c")) +
  ylab(expression(bold(MP~Concentration~(Particles~mL^-1)))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=14, family = "sans", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans", face = "bold", color = "black"),
        plot.title = element_text(hjust = -0.05, size = 14, family = "sans", face = "bold"),
        strip.text.x = element_text(size=10, family = "sans", face = "bold", color = "black"),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_x_discrete(breaks = c("BL", "FF"), labels = c("Blood","Follicular fluid"))


ggsave(FIG,
       file="figures/figure_s2.png",
       width = 6.86*1.6,
       height=6*1.6,
       units = "in",
       dpi = 600)
  