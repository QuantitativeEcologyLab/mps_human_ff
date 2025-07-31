#---------------------------------------------------------------------
# Data import and workspace preparation
#---------------------------------------------------------------------

#Import the packages
library(ggplot2)
library(ggpubr)
library(mgcv)

#Import the data
source("scripts/data_import.R")


#Remove two patients who opted to freeze their oocytes
total_mps_reduced <- total_mps[!total_mps$patient %in% c(55,59),]

#Build a "pregnancy boolean
total_mps_reduced$pregnancy <- 0
total_mps_reduced$pregnancy[total_mps_reduced$cycle_output == "pregnant"] <- 1

#---------------------------------------------------------------------
# Figure 1a - Stacked barplot of polymers in each patient
#---------------------------------------------------------------------

#subset to only the follicular fluid data
data_ff <- data[which(data$sample == "FF"),]


# Define the colors for each of the polymers
polymer_colors <- c(
  "EVA" = "#008d4d",
  "Cellulose" = "#fe2c8a",
  "PET" = "#41d75d",
  "Polyamide" = "#852c86",
  "Polyester" = "#939e00",
  "Polyethylene" = "#0182f3",
  "Polypropylene" = "#d5a400",
  "Polystyrene" = "#015695",
  "PVC" = "#ff856e",
  "Rubber" = "#00aea0",
  "Silicon" = "#a41132",
  "Teflon" = "#b6cf83",
  "Other" = "#ff85d2",
  "PEG" = "#6a4e05",
  "PVP" = "#dba4ff",
  "Polyvinylformal" = "#efbd86"
)

a <- 
  ggplot(data = data_ff, 
         aes(x = ID,
             y = concentration,
             fill = category)) +
  ggtitle("a") +
  geom_bar(stat = "identity", alpha = 0.9) +
  scale_fill_manual(labels = names(polymer_colors),
                    breaks = names(polymer_colors),
                    values = polymer_colors) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=8, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 12, family = "sans", face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=8, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.1, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_continuous(limits = c(0,580), expand = c(0,6)) +
  xlab(expression(bold(Patient))) +
  ylab(expression(bold(Total~MP~Concentration~(Particles~mL^-1)))) +
  guides(fill = guide_legend(ncol = 8, nrow = 2, byrow = TRUE))
  

#---------------------------------------------------------------------
# Figure 1b - Boxplots of FF MP concentrations and IVF outcomes
#---------------------------------------------------------------------

#Generate the figure
b <-
ggplot(data = total_mps_reduced,
       aes(x = cycle_output,
           y = concentration)) +
  ggtitle("b") +
  geom_boxplot(size = 0.1, outlier.size = 0.2, outlier.shape = 16, outlier.alpha = 0,
               fill = c("#fca311", "#14213d"), alpha = 0.7) +
  stat_summary(fun.y=mean, geom="point", shape=18, size=2, color="black", fill="black") +
  geom_jitter(aes(col = cycle_output), size = 0.5, shape = 16, position=position_jitter(height=0, width=0.1)) +
  scale_fill_manual(values = c("#fca311", "#14213d"), labels = c("Female", "Male")) +
  scale_colour_manual(values = c("#fca311", "#14213d"), labels = c("Female", "Male")) +
  ylab(expression(bold(Total~MP~Concentration~(Particles~mL^-1)))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans", face = "bold", color = "black"),
        plot.title = element_text(hjust = -0.05, size = 14, family = "sans", face = "bold"),
        #strip.text.x = element_text(size=6, family = "sans", face = "bold", color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_y_continuous(limits = c(0,580), expand = c(0,6)) +
  scale_x_discrete(breaks = c("not pregnant", "pregnant"), labels = c("Not pregnant","Pregnant"))





#---------------------------------------------------------------------
# Figure 1c - GLM of pregnancy vs. MP concentration
#---------------------------------------------------------------------


#Generate the figure
c <-
ggplot(data = total_mps_reduced,
       aes(x = concentration,
           y = pregnancy)) +
  ggtitle("c") +
  geom_smooth(method = "gam", formula = y ~ x, method.args = list(family = "binomial"), col = "grey30", alpha = 0.2) +
  geom_point(size = 1, shape = 16) +
  ylab(expression(bold(Probability~of~pregnancy))) + 
  xlab(expression(bold(Total~MP~Concentration~(Particles~mL^-1)))) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size=10, family = "sans", face = "bold"),
        axis.title.x = element_text(size=10, family = "sans", face = "bold"),
        axis.text.y = element_text(size=8, family = "sans"),
        axis.text.x  = element_text(size=10, family = "sans"),
        plot.title = element_text(hjust = -0.05, size = 14, family = "sans", face = "bold"),
        #strip.text.x = element_text(size=6, family = "sans", face = "bold", color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size=5, family = "sans", face = "bold"),
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.3, 'cm'),
        legend.spacing.y = unit(0.2, 'cm'),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = unit(c(0.2,0.1,0.2,0.2), "cm")) +
  scale_x_continuous(limits = c(0,520), expand = c(0,6))



#---------------------------------------------------------------------
# Compile and save
#---------------------------------------------------------------------

TOP <- ggarrange(a,b,
                 ncol=2, nrow=1,
                 widths = c(0.7, 0.3))

FIG <- ggarrange(TOP,c,
                 ncol=1, nrow=2)


#Save the figures
ggsave(FIG,
       width = 6.86*1.5, height = 6*1.5, units = "in",
       dpi = 600,
       bg = "transparent",
       file="figures/figure_1.png")

